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
        dict: Dictionary mapping chains to (min_res, max_res) tuples
    """
    residues_by_chain = {}
    chain_boundaries = {}
    
    for line in pdb_content.splitlines():
        # Include both ATOM and HETATM records
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain_id not in residues_by_chain:
                    residues_by_chain[chain_id] = set()
                    chain_boundaries[chain_id] = (float('inf'), -float('inf'))
                
                residues_by_chain[chain_id].add(res_num)
                
                # Update min and max residue numbers for this chain
                min_res, max_res = chain_boundaries[chain_id]
                chain_boundaries[chain_id] = (min(min_res, res_num), max(max_res, res_num))
            except (ValueError, IndexError):
                # Skip lines that can't be parsed properly
                continue
    
    return residues_by_chain, chain_boundaries


def correct_motif_numbering(motifs, chain_boundaries):
    """
    Attempt to correct motif numbering based on actual chain boundaries,
    but ONLY when the motif is completely outside boundaries or there's
    a clear shift pattern. Valid subsets of existing chains are not modified.
    
    Args:
        motifs (list): List of motif dictionaries
        chain_boundaries (dict): Dictionary of (min_res, max_res) tuples by chain
        
    Returns:
        list: Corrected list of motif dictionaries
        bool: Whether corrections were made
        str: Message describing corrections
    """
    corrected = False
    messages = []
    corrected_motifs = []
    
    for i, motif in enumerate(motifs):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        if chain in chain_boundaries:
            min_res, max_res = chain_boundaries[chain]
            
            # Only correct if the motif is completely outside the chain boundaries
            # or if there's a clear systematic shift
            new_start = start_res
            new_end = end_res
            
            # Case 1: Motif is completely before the chain starts
            if end_res < min_res:
                offset = min_res - end_res
                new_start = start_res + offset
                new_end = end_res + offset
                corrected = True
                messages.append(f"Adjusted motif {i+1} (completely outside chain): {start_res}-{end_res} to {new_start}-{new_end} (chain {chain} range: {min_res}-{max_res})")
            
            # Case 2: Motif is completely after the chain ends
            elif start_res > max_res:
                offset = max_res - start_res
                new_start = start_res + offset
                new_end = end_res + offset
                corrected = True
                messages.append(f"Adjusted motif {i+1} (completely outside chain): {start_res}-{end_res} to {new_start}-{new_end} (chain {chain} range: {min_res}-{max_res})")
            
            # Case 3: Clear systematic shift (both start and end are shifted by similar amount)
            # Only apply when BOTH ends are outside bounds and by similar amounts
            elif start_res < min_res and end_res > max_res:
                start_offset = min_res - start_res
                end_offset = end_res - max_res
                
                # If offsets are similar, might be a consistent numbering issue
                if abs(start_offset - end_offset) <= 2 and min(start_offset, end_offset) > 0:
                    avg_offset = (start_offset + end_offset) // 2
                    new_start = start_res + avg_offset
                    new_end = end_res - avg_offset
                    
                    corrected = True
                    messages.append(f"Adjusted motif {i+1} by offset {avg_offset}: {start_res}-{end_res} to {new_start}-{new_end} (chain {chain} range: {min_res}-{max_res})")
            
            corrected_motifs.append({
                "chain": chain,
                "start_res": new_start,
                "end_res": new_end,
                "tag": motif["tag"],
                "group": motif["group"]
            })
        else:
            # Chain not found, keep original
            corrected_motifs.append(motif.copy())
    
    return corrected_motifs, corrected, "\n".join(messages)


def validate_motifs(motifs, residues_by_chain, strict=False, auto_correct=False, chain_boundaries=None):
    """
    Validate that all motif residues exist in the PDB structure
    
    Args:
        motifs (list): List of motif dictionaries
        residues_by_chain (dict): Dictionary of available residues by chain
        strict (bool): If True, fails on any missing residue; if False, warns but continues
        auto_correct (bool): If True, attempts to fix numbering issues only when clearly necessary
        chain_boundaries (dict): Dictionary of (min_res, max_res) tuples by chain
        
    Returns:
        bool: True if validation passes, False otherwise
        list: List of validation error messages
        list: Potentially corrected motifs if auto_correct is True
    """
    all_valid = True
    errors = []
    corrected_motifs = motifs.copy()
    
    # First check for completely invalid motifs that might need correction
    completely_invalid_motifs = []
    for i, motif in enumerate(motifs):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        # Check if chain exists
        if chain not in residues_by_chain:
            completely_invalid_motifs.append(i)
            continue
            
        # Check if ANY residues in range exist
        any_valid = False
        for res in range(start_res, end_res + 1):
            if res in residues_by_chain[chain]:
                any_valid = True
                break
                
        if not any_valid:
            completely_invalid_motifs.append(i)
    
    # Only auto-correct when motifs are completely invalid
    correction_message = ""
    if auto_correct and chain_boundaries and completely_invalid_motifs:
        corrected_motifs, was_corrected, correction_message = correct_motif_numbering(motifs, chain_boundaries)
        if was_corrected:
            print("\nAuto-corrected motif numbering:")
            print(correction_message)
    
    # Validate all motifs (even after correction)
    for i, motif in enumerate(corrected_motifs):
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
    
    return all_valid, errors, corrected_motifs


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
            
            # For heterodimers: assign group based on chain identifier
            # First motif always gets group A
            # Motifs from different chains get different groups (A or B)
            if len(result["motifs"]) == 0:
                # First motif always gets group A
                group = "A"
            else:
                # Check if this motif is from a different chain than the first motif
                first_chain = result["motifs"][0]["chain"]
                group = "A" if chain == first_chain else "B"
            
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


def fix_pdb_for_salad(pdb_content, fix_residue_numbering=True, verbose=False):
    """
    Fix a PDB content specifically for SALAD to avoid the 
    'boolean index did not match shape of indexed array' error.
    
    This function uses two strategies:
    1. Renumbers residues to be continuous
    2. Ensures motif definitions cover full chains rather than just parts
    
    Args:
        pdb_content (str): The PDB file content as a string
        fix_residue_numbering (bool): Whether to renumber residues for SALAD compatibility
        verbose (bool): Whether to print detailed information
        
    Returns:
        str: The fixed PDB content
    """
    lines = pdb_content.splitlines()
    
    # Extract headers (REMARK lines) and non-REMARK, non-ATOM lines
    headers = []
    non_remark_lines = []
    
    for line in lines:
        if line.startswith('REMARK'):
            headers.append(line)
        elif not line.startswith(('ATOM', 'HETATM')):
            non_remark_lines.append(line)
    
    # Analyze chains and residues
    chains = {}
    atom_lines = []
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain_id = line[21]
            if chain_id not in chains:
                chains[chain_id] = {'residues': set(), 'mapping': {}}
            
            res_num = int(line[22:26])
            chains[chain_id]['residues'].add(res_num)
            atom_lines.append(line)
    
    # Extract motif definitions from REMARK lines
    motif_lines = []
    linker_lines = []
    other_remark_lines = []
    
    motif_defs = []
    
    for line in headers:
        if line.startswith('REMARK 999 INPUT') and len(line.split()) >= 6:
            parts = line.split()
            if len(parts[3]) == 1:  # This is a chain definition line
                try:
                    chain_id = parts[3]
                    start_res = int(parts[4])
                    end_res = int(parts[5])
                    group = parts[6] if len(parts) > 6 else "A"
                    
                    motif_defs.append({
                        'chain': chain_id,
                        'start': start_res,
                        'end': end_res,
                        'group': group
                    })
                    motif_lines.append(line)
                except (ValueError, IndexError):
                    other_remark_lines.append(line)
            else:
                # This might be a linker line
                linker_lines.append(line)
        else:
            other_remark_lines.append(line)
    
    if verbose:
        print("Original motif definitions:")
        for motif in motif_defs:
            print(f"Chain {motif['chain']}: {motif['start']}-{motif['end']} (Group {motif['group']})")
    
    # Create residue mappings to ensure continuous numbering
    for chain_id, chain_data in chains.items():
        sorted_residues = sorted(chain_data['residues'])
        for i, res_num in enumerate(sorted_residues, 1):
            chain_data['mapping'][res_num] = i
        
        if verbose:
            print(f"Chain {chain_id}: {len(sorted_residues)} residues")
            print(f"  Original range: {min(sorted_residues)}-{max(sorted_residues)}")
            print(f"  Remapped to: 1-{len(sorted_residues)}")
    
    # CRITICAL FIX FOR SALAD: Modify motifs to cover full chains
    # This is key to fixing the array size mismatch issue
    new_motif_lines = []
    
    for chain_id, chain_data in chains.items():
        min_res = min(chain_data['residues'])
        max_res = max(chain_data['residues'])
        
        # Find if this chain already has a motif
        existing_motif = None
        for motif in motif_defs:
            if motif['chain'] == chain_id:
                existing_motif = motif
                break
        
        # If no motif found for this chain, we need to add one
        if existing_motif is None:
            # Determine which group to assign - check other motifs
            group = "B" if len(motif_defs) > 0 and motif_defs[0]['group'] == "A" else "A"
            
            # Create a dummy motif covering the full chain
            if fix_residue_numbering:
                new_line = f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chain_data['residues']):3d} {group}"
            else:
                new_line = f"REMARK 999 INPUT  {chain_id} {min_res:3d} {max_res:3d} {group}"
            
            new_motif_lines.append(new_line)
            
            if verbose:
                print(f"Added new full-chain motif: Chain {chain_id}, {min_res}-{max_res} -> {1}-{len(chain_data['residues'])}, Group {group}")
        else:
            # Update existing motif to cover full chain
            if fix_residue_numbering:
                new_line = f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chain_data['residues']):3d} {existing_motif['group']}"
            else:
                new_line = f"REMARK 999 INPUT  {chain_id} {min_res:3d} {max_res:3d} {existing_motif['group']}"
            
            new_motif_lines.append(new_line)
            
            if verbose:
                print(f"Updated motif to full chain: Chain {chain_id}, {existing_motif['start']}-{existing_motif['end']} -> {min_res}-{max_res} -> {1}-{len(chain_data['residues'])}, Group {existing_motif['group']}")
    
    # Create new ATOM lines with renumbered residues if requested
    new_atoms = []
    if fix_residue_numbering:
        for line in atom_lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                old_res_num = int(line[22:26])
                new_res_num = chains[chain_id]['mapping'][old_res_num]
                
                # Replace the residue number in the line
                new_line = line[:22] + f"{new_res_num:4d}" + line[26:]
                new_atoms.append(new_line)
            else:
                new_atoms.append(line)
    else:
        new_atoms = atom_lines
    
    # Rebuild the file with modified headers
    new_remark_lines = other_remark_lines[:]
    
    # Add modified motif lines at the beginning of the REMARK blocks
    # Get the first two lines (NAME and PDB lines)
    name_line = None
    pdb_line = None
    
    for line in other_remark_lines:
        if "REMARK 999 NAME" in line:
            name_line = line
        elif "REMARK 999 PDB" in line:
            pdb_line = line
    
    # Rebuild the remarks in the correct order
    ordered_remarks = []
    if name_line:
        ordered_remarks.append(name_line)
    if pdb_line:
        ordered_remarks.append(pdb_line)
    
    # Add our modified motif lines
    ordered_remarks.extend(new_motif_lines)
    
    # Add the linker lines
    ordered_remarks.extend(linker_lines)
    
    # Add all other remarks except NAME and PDB lines
    for line in other_remark_lines:
        if line != name_line and line != pdb_line:
            ordered_remarks.append(line)
    
    # Combine everything into the final PDB content
    fixed_content = '\n'.join(ordered_remarks + new_atoms + non_remark_lines)
    
    return fixed_content


def process_single_pdb(args):
    pdb_path, rf_format, output_dir, options = args
    pdb_file = os.path.basename(pdb_path)
    pdb_name = os.path.splitext(pdb_file)[0]
    
    try:
        # Read PDB content first for validation
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        # Parse PDB structure for validation
        residues_by_chain, chain_boundaries = parse_pdb_structure(pdb_content)
        
        # Debug: print chain boundaries
        print(f"\nChain boundaries in {pdb_file}:")
        for chain, (min_res, max_res) in chain_boundaries.items():
            print(f"Chain {chain}: {min_res} - {max_res} ({len(residues_by_chain[chain])} residues)")
        
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
        
        # Print original motif definitions
        print(f"\nOriginal motif definitions:")
        for i, motif in enumerate(parsed_data["motifs"]):
            print(f"Motif {i+1}: Chain {motif['chain']}, Residues {motif['start_res']}-{motif['end_res']}")
        
        # Determine if validation should be strict or lenient
        strict_validation = not options.get('lenient_validation', False)
        
        # Determine if auto-correction should be used
        use_auto_correct = options.get('auto_correct', False)
        
        # Validate motifs (and potentially auto-correct if enabled)
        valid, errors, corrected_motifs = validate_motifs(
            parsed_data["motifs"], 
            residues_by_chain, 
            strict=strict_validation,
            auto_correct=use_auto_correct,
            chain_boundaries=chain_boundaries
        )
        
        # Update parsed data with corrected motifs if needed
        if use_auto_correct and corrected_motifs != parsed_data["motifs"]:
            parsed_data["motifs"] = corrected_motifs
            # Regenerate the header with corrected motifs
            genie2_header = converter.convert_to_genie2_format(parsed_data)
            
            # Print corrected motif definitions
            print(f"\nCorrected motif definitions:")
            for i, motif in enumerate(corrected_motifs):
                print(f"Motif {i+1}: Chain {motif['chain']}, Residues {motif['start_res']}-{motif['end_res']}")
        
        # Print validation errors/warnings
        for error in errors:
            if strict_validation:
                print(f"Error - {pdb_file}: {error}")
            else:
                print(f"Warning - {pdb_file}: {error}")
        
        # Handle validation failure
        if not valid and strict_validation:
            print(f"Error: Validation failed for {pdb_file}. Residues specified in motif definition do not exist in the PDB file.")
            print("You have several options:")
            print("1. Use --lenient flag to proceed despite missing residues")
            print("2. Use --auto_correct to attempt automatic fixing of numbering issues")
            print("3. Adjust your motif definition to match available residues")
            print("\nAvailable residue ranges:")
            for chain, (min_res, max_res) in chain_boundaries.items():
                print(f"Chain {chain}: {min_res} - {max_res}")
            return (pdb_path, False)
        elif not valid:
            print(f"Warning: Proceeding despite validation issues because --lenient flag was used.")
            print("Some motif residues might be missing in the final structure.")
        
        # Combine the PDB and header
        combined_content = genie2_header + "\n" + pdb_content
        
        # Always apply SALAD fix
        verbose = options.get('verbose', False)
        print(f"\nApplying SALAD compatibility fixes...")
        combined_content = fix_pdb_for_salad(
            combined_content, 
            fix_residue_numbering=True,
            verbose=verbose
        )
        
        # Create output file path
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        
        # Write combined content to output file
        with open(output_path, 'w') as f:
            f.write(combined_content)
        
        print(f"Successfully created {output_path}")
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
    parser = argparse.ArgumentParser(description='Convert RFDiffusion format to Genie2 format with SALAD compatibility')
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--pdb_file', '-f', help='Path to a single PDB file')
    input_group.add_argument('--pdb_dir', '-d', help='Path to a directory containing PDB files')
    
    spec_group = parser.add_mutually_exclusive_group(required=True)
    spec_group.add_argument('--input', '-i', help='RFDiffusion format string (e.g., "A1-80[M1]/30/[M2]B81-100")')
    spec_group.add_argument('--csv', '-c', help='CSV file mapping PDB files to RFDiffusion formats')
    spec_group.add_argument('--json', '-j', help='JSON file mapping PDB files to RFDiffusion formats')
    
    parser.add_argument('--output', '-o', required=True, help='Output directory for Genie2 files')
    
    # Processing options
    parser.add_argument('--sequential', '-s', action='store_true', help='Process files sequentially')
    parser.add_argument('--workers', '-w', type=int, default=None, help='Maximum number of parallel workers')
    
    # Scaffold options
    parser.add_argument('--add_terminal_scaffolds', action='store_true', help='Add scaffolds at the beginning and end')
    parser.add_argument('--min_scaffold', type=int, default=5, help='Minimum scaffold length')
    parser.add_argument('--max_scaffold', type=int, default=20, help='Maximum scaffold length')
    parser.add_argument('--min_factor', type=float, default=1.0, help='Minimum total length factor')
    parser.add_argument('--max_factor', type=float, default=1.5, help='Maximum total length factor')
    
    # Validation options
    parser.add_argument('--lenient', action='store_true', help='Enable lenient validation (only warn on missing residues)')
    parser.add_argument('--auto_correct', action='store_true', help='Enable automatic correction of motif numbering')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print detailed information')
    
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
        'lenient_validation': args.lenient,
        'auto_correct': args.auto_correct,
        'verbose': args.verbose
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
