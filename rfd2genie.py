#!/usr/bin/env python3

import os
import re
import argparse
import json
import csv
import glob
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor


def parse_pdb_structure(pdb_content):
    """
    Parse PDB content to extract residue information for validation
    
    Returns:
        dict: Dictionary mapping chains to sets of residue numbers
    """
    residues_by_chain = {}
    
    for line in pdb_content.splitlines():
        # Include both ATOM and HETATM records
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain_id not in residues_by_chain:
                    residues_by_chain[chain_id] = set()
                
                residues_by_chain[chain_id].add(res_num)
            except (ValueError, IndexError):
                # Skip lines that can't be parsed properly
                continue
    
    return residues_by_chain


def validate_motifs(motifs, residues_by_chain, strict=False):
    """
    Validate that all motif residues exist in the PDB structure
    
    Args:
        motifs (list): List of motif dictionaries
        residues_by_chain (dict): Dictionary of available residues by chain
        strict (bool): If True, fails on any missing residue; if False, warns but continues
        
    Returns:
        bool: True if validation passes, False otherwise
        list: List of validation error messages
    """
    all_valid = True
    errors = []
    
    for i, motif in enumerate(motifs):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        # Check if chain exists
        if chain not in residues_by_chain:
            error_msg = f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Chain not found in PDB"
            errors.append(error_msg)
            all_valid = False
            continue
            
        # Check if all residues in range exist
        missing_residues = []
        for res in range(start_res, end_res + 1):
            if res not in residues_by_chain[chain]:
                missing_residues.append(res)
        
        if missing_residues:
            if len(missing_residues) <= 5:
                missing_str = ", ".join(str(r) for r in missing_residues)
            else:
                missing_str = f"{len(missing_residues)} residues including {missing_residues[0]}, {missing_residues[1]}, ..."
                
            error_msg = f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Missing {missing_str}"
            errors.append(error_msg)
            
            # Only set all_valid to False if in strict mode
            if strict:
                all_valid = False
    
    return all_valid, errors


class RFDiffusionToGenie2Converter:
    def __init__(self, pdb_name="input_pdb", add_terminal_scaffolds=False, 
                 min_scaffold_length=5, max_scaffold_length=20,
                 min_total_length_factor=1.0, max_total_length_factor=1.5):
        self.pdb_name = pdb_name
        self.add_terminal_scaffolds = add_terminal_scaffolds
        self.min_scaffold_length = min_scaffold_length
        self.max_scaffold_length = max_scaffold_length
        self.min_total_length_factor = min_total_length_factor
        self.max_total_length_factor = max_total_length_factor
    
    def parse_rfdiffusion_format(self, rf_format):
        result = {
            "motifs": [],
            "linkers": [],
            "total_min_length": 0,
            "total_max_length": 0
        }
        
        parts = rf_format.split("/")
        
        for i, part in enumerate(parts):
            if not part.strip():
                continue  # Skip empty parts (e.g., trailing slash)
                
            if re.match(r'^\d+$', part):
                linker_length = int(part)
                min_length = max(5, linker_length - 5)
                max_length = linker_length + 5
                
                result["linkers"].append({
                    "min_length": min_length,
                    "max_length": max_length
                })
                
                result["total_min_length"] += min_length
                result["total_max_length"] += max_length
                continue
            
            motif_tag_match = re.search(r'\[(M\d+)\]', part)
            motif_tag = motif_tag_match.group(1) if motif_tag_match else f"M{len(result['motifs']) + 1}"
            
            clean_part = re.sub(r'\[M\d+\]', '', part)
            chain_match = re.match(r'^([A-Za-z])', clean_part)
            chain = chain_match.group(1) if chain_match else "A"
            
            range_match = re.search(r'(\d+)-(\d+)', clean_part)
            if not range_match:
                raise ValueError(f"Invalid motif format: {part}")
            
            start_res = int(range_match.group(1))
            end_res = int(range_match.group(2))
            
            # Ensure start_res <= end_res
            if start_res > end_res:
                start_res, end_res = end_res, start_res
                print(f"Warning: Swapped start and end residues for motif in '{part}'")
            
            group = "A" if len(result["motifs"]) < 2 else "B"
            
            result["motifs"].append({
                "chain": chain,
                "start_res": start_res,
                "end_res": end_res,
                "tag": motif_tag,
                "group": group
            })
            
            result["total_min_length"] += (end_res - start_res + 1)
            result["total_max_length"] += (end_res - start_res + 1)
        
        return result
    
    def convert_to_genie2_format(self, parsed_data):
        output = []
        
        name = self.pdb_name.split(".")[0] if "." in self.pdb_name else self.pdb_name
        output.append(f"REMARK 999 NAME   {name}_motifs")
        output.append(f"REMARK 999 PDB    {name}")
        
        if self.add_terminal_scaffolds:
            output.append(f"REMARK 999 INPUT     {self.min_scaffold_length:2d}  {self.max_scaffold_length:2d}")
        
        for i, motif in enumerate(parsed_data["motifs"]):
            output.append(
                f"REMARK 999 INPUT  {motif['chain']} {motif['start_res']:3d} {motif['end_res']:3d} {motif['group']}"
            )
            
            if i < len(parsed_data["linkers"]):
                linker = parsed_data["linkers"][i]
                output.append(
                    f"REMARK 999 INPUT     {linker['min_length']:2d}  {linker['max_length']:2d}"
                )
        
        if self.add_terminal_scaffolds:
            output.append(f"REMARK 999 INPUT     {self.min_scaffold_length:2d}  {self.max_scaffold_length:2d}")
        
        min_total_length = max(80, int(parsed_data["total_min_length"] * self.min_total_length_factor))
        max_total_length = max(200, int(parsed_data["total_max_length"] * self.max_total_length_factor))
        
        output.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total_length}")
        output.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total_length}")
        
        return "\n".join(output)
    
    def convert(self, rf_format):
        parsed_data = self.parse_rfdiffusion_format(rf_format)
        return self.convert_to_genie2_format(parsed_data), parsed_data


def read_specs_from_csv(csv_file):
    specs = {}
    try:
        with open(csv_file, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'pdb_file' in row and 'rf_format' in row:
                    specs[row['pdb_file']] = row['rf_format']
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return {}
    
    return specs


def read_specs_from_json(json_file):
    try:
        with open(json_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return {}


def process_single_pdb(args):
    pdb_path, rf_format, output_dir, options = args
    pdb_file = os.path.basename(pdb_path)
    pdb_name = os.path.splitext(pdb_file)[0]
    
    try:
        # Read PDB content first for validation
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        # Parse PDB structure for validation
        residues_by_chain = parse_pdb_structure(pdb_content)
        
        converter = RFDiffusionToGenie2Converter(
            pdb_name,
            add_terminal_scaffolds=options.get('add_terminal_scaffolds', False),
            min_scaffold_length=options.get('min_scaffold_length', 5),
            max_scaffold_length=options.get('max_scaffold_length', 20),
            min_total_length_factor=options.get('min_total_length_factor', 1.0),
            max_total_length_factor=options.get('max_total_length_factor', 1.5)
        )
        
        # Convert and get both the header and parsed data
        genie2_header, parsed_data = converter.convert(rf_format)
        
        # Always perform strict validation by default
        # Only use non-strict if explicitly requested with --lenient flag
        strict_validation = not options.get('lenient_validation', False)
        valid, errors = validate_motifs(parsed_data["motifs"], residues_by_chain, strict=strict_validation)
        
        # Print validation errors
        for error in errors:
            if strict_validation:
                print(f"Error - {pdb_file}: {error}")
            else:
                print(f"Warning - {pdb_file}: {error}")
        
        if not valid:
            if strict_validation:
                print(f"Error: Validation failed for {pdb_file}. Residues specified in motif definition do not exist in the PDB file. Skipping file.")
                return (pdb_path, False)
            else:
                print(f"Warning: Proceeding despite validation issues because --lenient flag was used.")
        
        # Create output file path
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        
        # Write combined content to output file
        with open(output_path, 'w') as f:
            f.write(genie2_header + "\n")
            f.write(pdb_content)
        
        return (pdb_path, True)
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return (pdb_path, False)


def process_pdb_files(pdb_sources, specs, output_dir, default_rf_format=None, 
                      parallel=True, max_workers=None, options=None):
    if options is None:
        options = {}
    
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_files = []
    if isinstance(pdb_sources, str):
        if os.path.isdir(pdb_sources):
            pdb_files = glob.glob(os.path.join(pdb_sources, "*.pdb"))
        elif os.path.isfile(pdb_sources) and pdb_sources.endswith('.pdb'):
            pdb_files = [pdb_sources]
    elif isinstance(pdb_sources, list):
        pdb_files = [f for f in pdb_sources if os.path.isfile(f) and f.endswith('.pdb')]
    
    if not pdb_files:
        print(f"No PDB files found")
        return (0, 0, [])
    
    tasks = []
    for pdb_path in pdb_files:
        pdb_file = os.path.basename(pdb_path)
        rf_format = specs.get(pdb_file, default_rf_format)
        
        if rf_format:
            tasks.append((pdb_path, rf_format, output_dir, options))
        else:
            print(f"Skipping {pdb_file} - no RF format specified and no default provided")
    
    results = []
    if parallel and len(tasks) > 1:
        print(f"Processing {len(tasks)} files in parallel...")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(process_single_pdb, tasks))
    else:
        print(f"Processing {len(tasks)} files sequentially...")
        results = [process_single_pdb(task) for task in tasks]
    
    success_count = sum(1 for _, success in results if success)
    failure_count = len(results) - success_count
    failed_files = [path for path, success in results if not success]
    
    return (success_count, failure_count, failed_files)


def main():
    parser = argparse.ArgumentParser(description='Convert RFDiffusion format to Genie2 format')
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--pdb_file', '-f', help='Path to a single PDB file')
    input_group.add_argument('--pdb_dir', '-d', help='Path to a directory containing PDB files')
    
    spec_group = parser.add_mutually_exclusive_group(required=True)
    spec_group.add_argument('--input', '-i', help='RFDiffusion format string (e.g., "A1-80[M1]/30/[M2]B81-100")')
    spec_group.add_argument('--csv', '-c', help='CSV file mapping PDB files to RFDiffusion formats')
    spec_group.add_argument('--json', '-j', help='JSON file mapping PDB files to RFDiffusion formats')
    
    parser.add_argument('--output', '-o', required=True, help='Output directory for Genie2 files')
    parser.add_argument('--sequential', '-s', action='store_true', help='Process files sequentially')
    parser.add_argument('--workers', '-w', type=int, default=None, help='Maximum number of parallel workers')
    parser.add_argument('--add_terminal_scaffolds', action='store_true', help='Add scaffolds at the beginning and end')
    parser.add_argument('--min_scaffold', type=int, default=5, help='Minimum scaffold length')
    parser.add_argument('--max_scaffold', type=int, default=20, help='Maximum scaffold length')
    parser.add_argument('--min_factor', type=float, default=1.0, help='Minimum total length factor')
    parser.add_argument('--max_factor', type=float, default=1.5, help='Maximum total length factor')
    parser.add_argument('--lenient', action='store_true', help='Enable lenient validation (only warn on missing residues)')
    
    args = parser.parse_args()
    
    specs = {}
    if args.csv:
        specs = read_specs_from_csv(args.csv)
    elif args.json:
        specs = read_specs_from_json(args.json)
    
    options = {
        'add_terminal_scaffolds': args.add_terminal_scaffolds,
        'min_scaffold_length': args.min_scaffold,
        'max_scaffold_length': args.max_scaffold,
        'min_total_length_factor': args.min_factor,
        'max_total_length_factor': args.max_factor,
        'lenient_validation': args.lenient
    }
    
    pdb_source = args.pdb_file if args.pdb_file else args.pdb_dir
    
    start_time = time.time()
    success, failure, failed_files = process_pdb_files(
        pdb_source, 
        specs, 
        args.output, 
        default_rf_format=args.input,
        parallel=not args.sequential,
        max_workers=args.workers,
        options=options
    )
    elapsed_time = time.time() - start_time
    
    print(f"\nProcessing complete in {elapsed_time:.2f} seconds")
    print(f"Successfully processed: {success} files")
    print(f"Failed to process: {failure} files")
    
    if failed_files:
        print("\nFailed files:")
        for file in failed_files:
            print(f"- {file}")


if __name__ == "__main__":
    main()
