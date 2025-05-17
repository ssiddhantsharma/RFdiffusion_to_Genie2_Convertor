#!/usr/bin/env python3

import os
import re
import argparse
import csv
import glob
from collections import defaultdict


def parse_pdb_structure(pdb_content):
    """Extract residue information and chain boundaries from PDB content"""
    residues_by_chain = {}
    chain_boundaries = {}
    atom_lines_by_residue = defaultdict(list)
    
    for line in pdb_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain_id not in residues_by_chain:
                    residues_by_chain[chain_id] = set()
                    chain_boundaries[chain_id] = (float('inf'), -float('inf'))
                
                residues_by_chain[chain_id].add(res_num)
                atom_lines_by_residue[(chain_id, res_num)].append(line)
                
                min_res, max_res = chain_boundaries[chain_id]
                chain_boundaries[chain_id] = (min(min_res, res_num), max(max_res, res_num))
            except (ValueError, IndexError):
                continue
    
    return residues_by_chain, chain_boundaries, atom_lines_by_residue


def validate_motifs(motifs, residues_by_chain):
    """Check if all specified residues exist in the PDB"""
    all_valid = True
    errors = []
    
    for i, motif in enumerate(motifs):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        if chain not in residues_by_chain:
            errors.append(f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Chain not found in PDB")
            all_valid = False
            continue
            
        missing_residues = []
        for res in range(start_res, end_res + 1):
            if res not in residues_by_chain[chain]:
                missing_residues.append(res)
        
        if missing_residues:
            if len(missing_residues) <= 5:
                missing_str = ", ".join(str(r) for r in missing_residues)
            else:
                missing_str = f"{len(missing_residues)} residues including {missing_residues[0]}, {missing_residues[1]}, ..."
                
            errors.append(f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Missing {missing_str}")
            all_valid = False
    
    return all_valid, errors


def fix_motifs(motifs, residues_by_chain):
    """Fix motifs to only include residues that exist in the PDB"""
    fixed_motifs = []
    
    for motif in motifs:
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        # Skip if chain doesn't exist
        if chain not in residues_by_chain:
            print(f"Warning: Chain {chain} not found in PDB, skipping motif")
            continue
        
        # Get the available residues in this chain
        available_residues = sorted(residues_by_chain[chain])
        if not available_residues:
            print(f"Warning: Chain {chain} exists but has no residues, skipping motif")
            continue
        
        # Adjust start_res to be within available residues
        if start_res not in residues_by_chain[chain]:
            closest_start = min(available_residues, key=lambda x: abs(x - start_res))
            print(f"Warning: Start residue {start_res} not found in chain {chain}, using {closest_start} instead")
            start_res = closest_start
        
        # Adjust end_res to be within available residues
        if end_res not in residues_by_chain[chain]:
            closest_end = min(available_residues, key=lambda x: abs(x - end_res))
            print(f"Warning: End residue {end_res} not found in chain {chain}, using {closest_end} instead")
            end_res = closest_end
        
        # Make sure start < end
        if start_res > end_res:
            start_res, end_res = end_res, start_res
            print(f"Warning: Swapped start and end residues for chain {chain}")
        
        # Create fixed motif
        fixed_motif = motif.copy()
        fixed_motif["start_res"] = start_res
        fixed_motif["end_res"] = end_res
        fixed_motifs.append(fixed_motif)
        
        print(f"Fixed motif: Chain {chain}, Residues {start_res}-{end_res}")
    
    return fixed_motifs


class MotifConverter:
    """Convert between RFDiffusion and Genie2 formats"""
    
    def __init__(self, pdb_name="input_pdb", add_terminal_scaffolds=False, 
                 min_scaffold_length=5, max_scaffold_length=20):
        self.pdb_name = pdb_name
        self.add_terminal_scaffolds = add_terminal_scaffolds
        self.min_scaffold_length = min_scaffold_length
        self.max_scaffold_length = max_scaffold_length
    
    def parse_rfdiffusion_format(self, rf_format):
        """Parse RFDiffusion format string into motifs and linkers"""
        result = {
            "motifs": [],
            "linkers": [],
            "total_min_length": 0,
            "total_max_length": 0
        }
        
        parts = rf_format.split("/")
        
        for i, part in enumerate(parts):
            if not part.strip():
                continue
                
            # Check if this is a simple linker (just a number)
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
            
            # Check if this is a linker range (e.g., "30-40")
            linker_range_match = re.match(r'^(\d+)-(\d+)$', part)
            if linker_range_match and not re.match(r'^[A-Za-z]', part):
                min_length = int(linker_range_match.group(1))
                max_length = int(linker_range_match.group(2))
                
                result["linkers"].append({
                    "min_length": min_length,
                    "max_length": max_length
                })
                
                result["total_min_length"] += min_length
                result["total_max_length"] += max_length
                continue
            
            # If we get here, it should be a motif
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
            
            if start_res > end_res:
                start_res, end_res = end_res, start_res
                print(f"Warning: Swapped start and end residues for motif in '{part}'")
            
            # Determine motif group (A for first chain, B for second chain)
            if len(result["motifs"]) == 0:
                group = "A"
            else:
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
        """Generate Genie2 format according to SALAD specifications"""
        output = []
        
        name = self.pdb_name.split(".")[0] if "." in self.pdb_name else self.pdb_name
        output.append(f"REMARK 999 NAME   {name}_motifs")
        output.append(f"REMARK 999 PDB    {name}")
        
        # Add N-terminal scaffold if requested
        if self.add_terminal_scaffolds:
            output.append(f"REMARK 999 INPUT    {self.min_scaffold_length:2d} {self.max_scaffold_length:2d}")
        
        # Add motifs and linkers
        for i, motif in enumerate(parsed_data["motifs"]):
            # Proper spacing for motif segments
            output.append(
                f"REMARK 999 INPUT  {motif['chain']} {motif['start_res']:3d} {motif['end_res']:3d} {motif['group']}"
            )
            
            # Add linker after motif if available
            if i < len(parsed_data["linkers"]):
                linker = parsed_data["linkers"][i]
                # Proper spacing for scaffold segments
                output.append(
                    f"REMARK 999 INPUT    {linker['min_length']:2d} {linker['max_length']:2d}"
                )
        
        # Add C-terminal scaffold if requested
        if self.add_terminal_scaffolds:
            output.append(f"REMARK 999 INPUT    {self.min_scaffold_length:2d} {self.max_scaffold_length:2d}")
        
        # Calculate total length constraints
        min_total_length = max(80, int(parsed_data["total_min_length"]))
        max_total_length = max(200, int(parsed_data["total_max_length"] * 1.5))
        
        # Ensure length limits make sense
        if min_total_length > max_total_length:
            max_total_length = min_total_length + 100
        
        # Proper spacing for total length specifications
        output.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total_length}")
        output.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total_length}")
        
        return "\n".join(output)
    
    def convert(self, rf_format):
        """Convert RFDiffusion format to Genie2 format"""
        parsed_data = self.parse_rfdiffusion_format(rf_format)
        return self.convert_to_genie2_format(parsed_data), parsed_data


def parse_genie2_motifs(pdb_content):
    """Extract motif and linker definitions from Genie2 format PDB headers"""
    motifs = []
    linkers = []
    min_total = None
    max_total = None
    
    for line in pdb_content.splitlines():
        if line.startswith('REMARK 999 INPUT') and len(line.split()) >= 4:
            parts = line.split()
            # This is a chain/motif definition with at least 4 parts
            if len(parts) >= 6 and len(parts[3]) == 1 and parts[3].isalpha():
                try:
                    chain_id = parts[3]
                    start_res = int(parts[4])
                    end_res = int(parts[5])
                    group = parts[6] if len(parts) > 6 else "A"
                    
                    motif_tag = f"M{len(motifs) + 1}"
                    motifs.append({
                        "chain": chain_id,
                        "start_res": start_res,
                        "end_res": end_res,
                        "tag": motif_tag,
                        "group": group
                    })
                except (ValueError, IndexError):
                    pass
            # This is a linker definition
            elif len(parts) >= 5:
                try:
                    min_length = int(parts[3])
                    max_length = int(parts[4])
                    linkers.append({
                        "min_length": min_length,
                        "max_length": max_length
                    })
                except (ValueError, IndexError):
                    pass
        elif line.startswith('REMARK 999 MINIMUM TOTAL LENGTH'):
            try:
                min_total = int(line.split()[-1])
            except (ValueError, IndexError):
                pass
        elif line.startswith('REMARK 999 MAXIMUM TOTAL LENGTH'):
            try:
                max_total = int(line.split()[-1])
            except (ValueError, IndexError):
                pass
    
    result = {"motifs": motifs, "linkers": linkers}
    
    if min_total is not None:
        result["min_total"] = min_total
    if max_total is not None:
        result["max_total"] = max_total
        
    return result


def fix_pdb_for_salad(pdb_content, verbose=False):
    """
    Fix a PDB content specifically for SALAD to avoid the 
    'boolean index did not match shape of indexed array' error.
    
    This function ensures:
    1. Motif definitions cover full chains
    2. Residues are continuously numbered
    3. All chains in the PDB have a motif definition
    """
    lines = pdb_content.splitlines()
    
    # Extract headers, atoms, and other content
    headers = [line for line in lines if line.startswith('REMARK')]
    atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
    other_lines = [line for line in lines if not line.startswith(('REMARK', 'ATOM', 'HETATM'))]
    
    # Analyze chains and residues
    chains = {}
    for line in atom_lines:
        chain_id = line[21]
        if chain_id not in chains:
            chains[chain_id] = {'residues': set(), 'mapping': {}}
        
        res_num = int(line[22:26])
        chains[chain_id]['residues'].add(res_num)
    
    # Create residue mappings for continuous numbering
    for chain_id, chain_data in chains.items():
        sorted_residues = sorted(chain_data['residues'])
        for i, res_num in enumerate(sorted_residues, 1):
            chain_data['mapping'][res_num] = i
        
        if verbose:
            print(f"Chain {chain_id}: {len(sorted_residues)} residues")
            print(f"  Original range: {min(sorted_residues)}-{max(sorted_residues)}")
            print(f"  Remapped to: 1-{len(sorted_residues)}")
    
    # Extract existing motif definitions
    motif_defs = []
    for line in headers:
        if line.startswith('REMARK 999 INPUT') and len(line.split()) >= 6:
            parts = line.split()
            if len(parts[3]) == 1 and parts[3].isalpha():
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
                except (ValueError, IndexError):
                    pass
    
    # Extract linker and other REMARK lines
    name_line = None
    pdb_line = None
    linker_lines = []
    other_remark_lines = []
    
    for line in headers:
        if line.startswith('REMARK 999 NAME'):
            name_line = line
        elif line.startswith('REMARK 999 PDB'):
            pdb_line = line
        elif line.startswith('REMARK 999 INPUT') and len(line.split()) >= 5:
            parts = line.split()
            # Check if this is a linker line (not a motif line)
            if len(parts) == 5 and not (len(parts[3]) == 1 and parts[3].isalpha()):
                linker_lines.append(line)
        elif not line.startswith(('REMARK 999 INPUT', 'REMARK 999 MINIMUM TOTAL LENGTH', 'REMARK 999 MAXIMUM TOTAL LENGTH')):
            other_remark_lines.append(line)
    
    # Create new motif lines ensuring full chain coverage
    new_motif_lines = []
    processed_chains = set(motif['chain'] for motif in motif_defs)
    
    # Update existing motifs to cover full chains
    for motif in motif_defs:
        chain_id = motif['chain']
        if chain_id in chains:
            # Create updated motif covering full chain
            new_line = f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chains[chain_id]['residues']):3d} {motif['group']}"
            new_motif_lines.append(new_line)
    
    # Add motifs for chains without existing motifs
    for chain_id in chains:
        if chain_id not in processed_chains:
            # Determine group
            group = "B" if motif_defs and motif_defs[0]['group'] == "A" else "A"
            new_line = f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chains[chain_id]['residues']):3d} {group}"
            new_motif_lines.append(new_line)
    
    # Renumber residues in ATOM lines
    new_atoms = []
    for line in atom_lines:
        chain_id = line[21]
        old_res_num = int(line[22:26])
        
        if chain_id in chains and old_res_num in chains[chain_id]['mapping']:
            new_res_num = chains[chain_id]['mapping'][old_res_num]
            new_line = line[:22] + f"{new_res_num:4d}" + line[26:]
            new_atoms.append(new_line)
        else:
            new_atoms.append(line)
    
    # Build new PDB content
    new_lines = []
    
    # Add name and PDB lines first
    if name_line:
        new_lines.append(name_line)
    if pdb_line:
        new_lines.append(pdb_line)
    
    # Add motif lines
    new_lines.extend(new_motif_lines)
    
    # Add linker lines
    new_lines.extend(linker_lines)
    
    # Add other REMARK lines
    for line in other_remark_lines:
        if not line.startswith(("REMARK 220")):  # Filter out non-standard REMARK lines
            new_lines.append(line)
    
    # Calculate and add total length if missing
    if not any(line.startswith("REMARK 999 MINIMUM TOTAL LENGTH") for line in new_lines):
        total_residues = sum(len(chain_data['residues']) for chain_data in chains.values())
        min_total = max(80, total_residues)
        new_lines.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total}")
    
    if not any(line.startswith("REMARK 999 MAXIMUM TOTAL LENGTH") for line in new_lines):
        if 'min_total' in locals():
            max_total = max(200, int(min_total * 1.5))
        else:
            total_residues = sum(len(chain_data['residues']) for chain_data in chains.values())
            max_total = max(200, int(total_residues * 1.5))
        new_lines.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total}")
    
    # Add ATOM lines and other content
    new_lines.extend(new_atoms)
    new_lines.extend(other_lines)
    
    return "\n".join(new_lines)


def process_file(pdb_path, rf_format, output_dir, verbose=False):
    """Process a single PDB file with given format string"""
    pdb_file = os.path.basename(pdb_path)
    pdb_name = os.path.splitext(pdb_file)[0]
    
    try:
        # Read PDB content
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        # Parse PDB structure for validation
        residues_by_chain, chain_boundaries = parse_pdb_structure(pdb_content)[:2]
        
        print(f"\nChain boundaries in {pdb_file}:")
        for chain, (min_res, max_res) in chain_boundaries.items():
            print(f"Chain {chain}: {min_res} - {max_res} ({len(residues_by_chain[chain])} residues)")
        
        output_content = None
        
        # Check if this is already a Genie2 formatted file
        if rf_format.lower() == "genie2" or rf_format.lower() == "auto":
            # Extract motif definitions from the PDB file
            genie2_data = parse_genie2_motifs(pdb_content)
            
            if genie2_data["motifs"]:
                print(f"\nDetected Genie2 format with {len(genie2_data['motifs'])} motifs")
                
                # Apply SALAD fixes directly to the PDB
                output_content = fix_pdb_for_salad(pdb_content, verbose)
            elif rf_format.lower() == "auto":
                print("No Genie2 motif definitions found, trying RFDiffusion format...")
            else:
                print("Error: No Genie2 motif definitions found in the PDB file")
                return False
        
        # If not already processed as Genie2, convert from RFDiffusion format
        if output_content is None:
            # Initialize converter
            converter = MotifConverter(
                pdb_name,
                add_terminal_scaffolds=False
            )
            
            # Convert RFDiffusion format to Genie2
            genie2_header, parsed_data = converter.convert(rf_format)
            
            print(f"\nOriginal motif definitions:")
            for i, motif in enumerate(parsed_data["motifs"]):
                print(f"Motif {i+1}: Chain {motif['chain']}, Residues {motif['start_res']}-{motif['end_res']}")
            
            # Validate motifs
            valid, errors = validate_motifs(parsed_data["motifs"], residues_by_chain)
            
            if not valid:
                for error in errors:
                    print(f"Error - {pdb_file}: {error}")
                
                # Auto-fix motifs
                print("Fixing invalid motif definitions...")
                fixed_motifs = fix_motifs(parsed_data["motifs"], residues_by_chain)
                
                if not fixed_motifs:
                    print("Error: No valid motifs could be created after fixing.")
                    return False
                
                # Generate new header with fixed motifs
                parsed_data["motifs"] = fixed_motifs
                genie2_header = converter.convert_to_genie2_format(parsed_data)
            
            # Combine header and PDB content
            combined_content = genie2_header + "\n" + pdb_content
            
            # Apply SALAD fixes
            output_content = fix_pdb_for_salad(combined_content, verbose)
        
        # Write to output file
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        with open(output_path, 'w') as f:
            f.write(output_content)
        
        print(f"Successfully created {output_path}")
        return True
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return False


def main():
    """Main function to process PDB files"""
    parser = argparse.ArgumentParser(description='Convert RFDiffusion format to Genie2 format with SALAD compatibility')
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--pdb_file', '-f', help='Path to a single PDB file')
    input_group.add_argument('--pdb_dir', '-d', help='Path to a directory containing PDB files')
    
    spec_group = parser.add_mutually_exclusive_group(required=True)
    spec_group.add_argument('--input', '-i', help='RFDiffusion format string or "genie2" to process existing Genie2 PDB')
    spec_group.add_argument('--csv', '-c', help='CSV file mapping PDB files to RFDiffusion formats')
    
    parser.add_argument('--output', '-o', required=True, help='Output directory for processed PDB files')
    parser.add_argument('--add_terminal_scaffolds', action='store_true', help='Add scaffolds at the beginning and end')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print detailed information')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Read specifications from CSV if provided
    specs = {}
    if args.csv:
        try:
            with open(args.csv, 'r', newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'pdb_file' in row and 'rf_format' in row:
                        specs[row['pdb_file']] = row['rf_format']
        except Exception as e:
            print(f"Error reading CSV file: {e}")
    
    # Process files
    print(f"\n--- Starting SALAD Motif Converter ---")
    
    pdb_files = []
    if args.pdb_file:
        pdb_files = [args.pdb_file]
    elif args.pdb_dir:
        pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
    
    if not pdb_files:
        print("No PDB files found")
        return
    
    print(f"Processing {len(pdb_files)} files...")
    success_count = 0
    failure_count = 0
    failed_files = []
    
    for pdb_path in pdb_files:
        pdb_file = os.path.basename(pdb_path)
        rf_format = specs.get(pdb_file, args.input)
        
        if rf_format:
            success = process_file(pdb_path, rf_format, args.output, args.verbose)
            if success:
                success_count += 1
            else:
                failure_count += 1
                failed_files.append(pdb_file)
        else:
            print(f"Skipping {pdb_file} - no RF format specified")
            failed_files.append(pdb_file)
    
    # Print summary
    print(f"\nProcessing complete")
    print(f"Successfully processed: {success_count} files")
    print(f"Failed to process: {failure_count} files")
    
    if failed_files:
        print("\nFailed files:")
        for file in failed_files:
            print(f"- {file}")
    
    print("\nOutput files are ready for SALAD 🧬")


if __name__ == "__main__":
    main()
