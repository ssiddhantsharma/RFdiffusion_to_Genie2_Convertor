#!/usr/bin/env python3
import os
import re
import argparse
import glob
import csv
from collections import defaultdict, OrderedDict

def parse_pdb(pdb_content):
    residues_by_chain = {}
    chain_boundaries = {}
    residue_atoms = defaultdict(list)
    chain_order = []
    
    for line in pdb_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain_id not in residues_by_chain:
                    residues_by_chain[chain_id] = set()
                    chain_boundaries[chain_id] = (float('inf'), -float('inf'))
                    chain_order.append(chain_id)
                
                residues_by_chain[chain_id].add(res_num)
                min_res, max_res = chain_boundaries[chain_id]
                chain_boundaries[chain_id] = (min(min_res, res_num), max(max_res, res_num))
                residue_atoms[(chain_id, res_num)].append(line)
            except:
                continue
    
    chain_order = list(OrderedDict.fromkeys(chain_order))
    return residues_by_chain, chain_boundaries, residue_atoms, chain_order

def parse_rf_format(rf_format, chain_order=None):
    motifs = []
    inter_motif_linkers = []
    initial_scaffold = None
    final_scaffold = None
    
    parts = [p.strip() for p in rf_format.split("/") if p.strip()]
    
    # Handle initial scaffold (N-terminal extension)
    if parts and re.match(r'^\d+(-\d+)?$', parts[0]) and not any(c.isalpha() for c in parts[0]):
        part = parts.pop(0)
        if "-" in part:
            min_length, max_length = map(int, part.split("-"))
        else:
            length = int(part)
            min_length = max_length = length
        initial_scaffold = {"min_length": min_length, "max_length": max_length}
    
    # Handle final scaffold (C-terminal extension)
    if parts and re.match(r'^\d+(-\d+)?$', parts[-1]) and not any(c.isalpha() for c in parts[-1]):
        part = parts.pop(-1)
        if "-" in part:
            min_length, max_length = map(int, part.split("-"))
        else:
            length = int(part)
            min_length = max_length = length
        final_scaffold = {"min_length": min_length, "max_length": max_length}
    
    # Initialize chain-to-group mapping using the chain order from PDB
    chain_to_group = {}
    if chain_order:
        for idx, chain in enumerate(chain_order):
            group = chr(ord('A') + idx % 26)
            chain_to_group[chain] = group
    
    # Process motifs and linkers in alternating pattern
    i = 0
    while i < len(parts):
        part = parts[i]
        
        # This is a motif definition
        if re.match(r'^[A-Za-z]', part) or re.search(r'\[(M\d+)\]', part):
            tag_match = re.search(r'\[(M\d+)\]', part)
            tag = tag_match.group(1) if tag_match else f"M{len(motifs) + 1}"
            
            clean_part = re.sub(r'\[M\d+\]', '', part)
            chain_match = re.match(r'^([A-Za-z])', clean_part)
            chain = chain_match.group(1) if chain_match else "A"
            
            range_match = re.search(r'(\d+)-(\d+)', clean_part)
            if not range_match:
                i += 1
                continue
            
            start_res = int(range_match.group(1))
            end_res = int(range_match.group(2))
            
            # Assign group based on chain
            if chain not in chain_to_group:
                next_group = chr(ord('A') + len(chain_to_group) % 26)
                chain_to_group[chain] = next_group
            
            group = chain_to_group[chain]
            
            # Check for explicit group override at end of motif string
            group_match = re.search(r'\s([A-Za-z])$', clean_part)
            if group_match:
                group = group_match.group(1)
            
            motifs.append({
                "chain": chain,
                "start_res": start_res,
                "end_res": end_res,
                "tag": tag,
                "group": group,
                "part_index": i
            })
            
        # This is a linker definition
        elif re.match(r'^\d+(-\d+)?$', part):
            if motifs:
                prev_motif_idx = len(motifs) - 1
                
                if "-" in part:
                    min_length, max_length = map(int, part.split("-"))
                else:
                    length = int(part)
                    min_length = max_length = length
                
                inter_motif_linkers.append({
                    "min_length": min_length,
                    "max_length": max_length,
                    "after_motif_index": prev_motif_idx
                })
        
        i += 1
    
    # Clean up temporary indices
    for motif in motifs:
        if "part_index" in motif:
            del motif["part_index"]
    
    return {
        "motifs": motifs,
        "inter_motif_linkers": inter_motif_linkers,
        "initial_scaffold": initial_scaffold,
        "final_scaffold": final_scaffold,
        "chain_to_group": chain_to_group
    }

def validate_motifs(motifs, residues_by_chain, verbose=False):
    """Verifies that specified residues exist in the PDB structure"""
    valid = True
    validation_results = []
    
    for i, motif in enumerate(motifs):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        # Check if chain exists
        if chain not in residues_by_chain:
            validation_results.append(f"Error: Chain {chain} not found in PDB")
            valid = False
            continue
        
        # Check for missing residues
        missing = []
        for res in range(start_res, end_res + 1):
            if res not in residues_by_chain[chain]:
                missing.append(res)
        
        if missing:
            msg = f"Error: Motif {i+1} (Chain {chain}:{start_res}-{end_res}) missing residues: {missing[:5]}"
            if len(missing) > 5:
                msg += f"... and {len(missing)-5} more"
            validation_results.append(msg)
            
            if verbose:
                available = sorted(list(residues_by_chain[chain]))
                res_msg = f"Available residues in chain {chain}: {available[:10]}"
                if len(available) > 10:
                    res_msg += f"... and {len(available)-10} more"
                validation_results.append(res_msg)
            
            valid = False
    
    # Print validation results
    if not valid:
        for result in validation_results:
            print(result)
    
    return valid

def calculate_length_constraints(parsed_data, residues_by_chain):
    """Calculate minimum and maximum total sequence lengths based on motifs and scaffolds"""
    
    # Start with motif residues (these are fixed lengths)
    motif_length = 0
    for motif in parsed_data["motifs"]:
        motif_length += motif["end_res"] - motif["start_res"] + 1
    
    # Add scaffold/linker lengths (these can be variable)
    min_scaffolds = 0
    max_scaffolds = 0
    
    # Initial scaffold (N-terminal extension)
    if parsed_data["initial_scaffold"]:
        min_scaffolds += parsed_data["initial_scaffold"]["min_length"]
        max_scaffolds += parsed_data["initial_scaffold"]["max_length"]
    
    # Inter-motif linkers
    for linker in parsed_data["inter_motif_linkers"]:
        min_scaffolds += linker["min_length"]
        max_scaffolds += linker["max_length"]
    
    # Final scaffold (C-terminal extension)
    if parsed_data["final_scaffold"]:
        min_scaffolds += parsed_data["final_scaffold"]["min_length"]
        max_scaffolds += parsed_data["final_scaffold"]["max_length"]
    
    min_total = motif_length + min_scaffolds
    max_total = motif_length + max_scaffolds
    
    return min_total, max_total

def generate_genie2_header(pdb_name, parsed_data, min_total, max_total):
    """
    Generate Genie2 header with precise column formatting for SALAD compatibility
    Based on official Genie2 format from 2b5i.pdb
    """
    genie2_header = []
    
    # PDB Name lines - exact column formatting required for SALAD
    genie2_header.append(f"REMARK 999 NAME   {pdb_name}")
    genie2_header.append(f"REMARK 999 PDB    {pdb_name}")
    
    # Add initial scaffold (N-terminal extension) if present
    if parsed_data["initial_scaffold"]:
        min_len = parsed_data["initial_scaffold"]["min_length"]
        max_len = parsed_data["initial_scaffold"]["max_length"]
        # Format: "REMARK 999 INPUT      5  15" (6 spaces after INPUT, 2-digit numbers, 2 spaces between)
        genie2_header.append(f"REMARK 999 INPUT      {min_len:2d}  {max_len:2d}")
    
    # Process motifs and any inter-motif linkers
    for i, motif in enumerate(parsed_data["motifs"]):
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        group = motif["group"]
        
        # Format: "REMARK 999 INPUT  A  11  23 B" (2 spaces after INPUT, 1 space after chain, 2 spaces between numbers, 1 space before group)
        genie2_header.append(f"REMARK 999 INPUT  {chain}  {start_res:2d}  {end_res:2d} {group}")
        
        # Add linker after this motif if present
        for linker in parsed_data["inter_motif_linkers"]:
            if linker["after_motif_index"] == i:
                min_len = linker["min_length"]
                max_len = linker["max_length"]
                # Format: "REMARK 999 INPUT     10  20" (5 spaces after INPUT)
                genie2_header.append(f"REMARK 999 INPUT     {min_len:2d}  {max_len:2d}")
    
    # Add final scaffold (C-terminal extension) if present
    if parsed_data["final_scaffold"]:
        min_len = parsed_data["final_scaffold"]["min_length"]
        max_len = parsed_data["final_scaffold"]["max_length"]
        genie2_header.append(f"REMARK 999 INPUT      {min_len:2d}  {max_len:2d}")
    
    # Add length constraints - exact formatting per documentation
    genie2_header.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total}")
    genie2_header.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total}")
    
    return "\n".join(genie2_header)

def ensure_column_formatting(line):
    """Fix existing REMARK 999 lines to have proper column formatting"""
    if line.startswith("REMARK 999 INPUT  ") and len(line) > 16:
        parts = line.split()
        if len(parts) >= 5 and len(parts[3]) == 1 and parts[3].isalpha():
            chain = parts[3]
            try:
                start_res = int(parts[4])
                end_res = int(parts[5])
                group = parts[6] if len(parts) > 6 else chain
                
                # Format to match official Genie2: "REMARK 999 INPUT  A  11  23 B"
                new_line = f"REMARK 999 INPUT  {chain}  {start_res:2d}  {end_res:2d} {group}"
                return new_line
            except (ValueError, IndexError):
                pass
        
    elif line.startswith("REMARK 999 INPUT     ") or line.startswith("REMARK 999 INPUT      "):
        parts = line.split()
        if len(parts) >= 4 and parts[3].isdigit():
            try:
                min_len = int(parts[3])
                max_len = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else min_len
                
                # Different spacing for initial/final scaffolds vs inter-motif linkers
                # Initial/final scaffolds: 6 spaces after INPUT
                # Inter-motif linkers: 5 spaces after INPUT
                # Default to 5 spaces (inter-motif) unless we can determine it's initial/final
                new_line = f"REMARK 999 INPUT     {min_len:2d}  {max_len:2d}"
                return new_line
            except (ValueError, IndexError):
                pass
            
    elif "MINIMUM TOTAL LENGTH" in line or "MAXIMUM TOTAL LENGTH" in line:
        parts = line.split()
        length_type = "MINIMUM" if "MINIMUM" in line else "MAXIMUM"
        value = parts[-1]
        
        # EXACT column positioning per SALAD requirements
        return f"REMARK 999 {length_type} TOTAL LENGTH      {value}"
        
    return line

def reorder_pdb_residues(pdb_content, parsed_data, do_reordering=False):
    """
    Reorders residues in the PDB file to match the order specified in RFDiffusion format.
    If do_reordering=False, just cleans existing REMARK 999 lines without changing order.
    """
    residues_by_chain, _, residue_atoms, _ = parse_pdb(pdb_content)
    
    if not do_reordering:
        # Don't reorder, just clean existing REMARK 999 lines
        clean_lines = []
        for line in pdb_content.splitlines():
            if not line.startswith("REMARK 999"):
                clean_lines.append(line)
        return "\n".join(clean_lines)
    
    # Reorder atoms according to motif specifications
    ordered_lines = []
    header_lines = []
    footer_lines = []
    
    # Extract header and footer lines
    for line in pdb_content.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith(("END", "TER")):
                footer_lines.append(line)
            elif not line.startswith("REMARK 999"):  # Skip existing REMARK 999 lines
                header_lines.append(line)
    
    ordered_lines.extend(header_lines)
    
    # Sort motifs by group, then chain, then start residue
    sorted_motifs = sorted(parsed_data["motifs"], key=lambda m: (m.get("group", "Z"), m["chain"], m["start_res"]))
    
    # Add atoms from each motif in order
    for motif in sorted_motifs:
        chain = motif["chain"]
        start_res = motif["start_res"]
        end_res = motif["end_res"]
        
        for res in range(start_res, end_res + 1):
            if (chain, res) in residue_atoms:
                for atom_line in residue_atoms[(chain, res)]:
                    ordered_lines.append(atom_line)
    
    ordered_lines.extend(footer_lines)
    return "\n".join(ordered_lines)

def process_file(pdb_path, rf_format, output_dir, verbose=False, reorder_residues=False):
    """
    Process a single PDB file with the specified format
    """
    try:
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        pdb_file = os.path.basename(pdb_path)
        pdb_name = os.path.splitext(pdb_file)[0]
        
        # Parse PDB structure
        residues_by_chain, chain_boundaries, _, chain_order = parse_pdb(pdb_content)
        
        if verbose:
            print(f"\nPDB file: {pdb_file}")
            print(f"Chains detected: {', '.join(chain_order)}")
            for chain, (min_res, max_res) in chain_boundaries.items():
                print(f"Chain {chain}: residues {min_res}-{max_res} ({len(residues_by_chain[chain])} total)")
        
        # Special case for existing Genie2 files
        if rf_format.lower() == "genie2":
            fixed_lines = []
            for line in pdb_content.splitlines():
                if line.startswith("REMARK 999"):
                    fixed_lines.append(ensure_column_formatting(line))
                else:
                    fixed_lines.append(line)
            
            output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
            with open(output_path, 'w') as f:
                f.write("\n".join(fixed_lines))
            
            print(f"Successfully processed existing Genie2 file: {output_path}")
            return True
        
        # Parse RF format specification
        parsed_data = parse_rf_format(rf_format, chain_order)
        
        if verbose:
            print("\nParsed RFDiffusion format:")
            print(f"  Chain to Group mapping: {parsed_data['chain_to_group']}")
            print(f"  Initial scaffold: {parsed_data['initial_scaffold']}")
            print(f"  Motifs ({len(parsed_data['motifs'])}):")
            for i, motif in enumerate(parsed_data["motifs"]):
                print(f"    {i+1}: Chain {motif['chain']} residues {motif['start_res']}-{motif['end_res']} (Group {motif['group']})")
            print(f"  Linkers ({len(parsed_data['inter_motif_linkers'])}):")
            for i, linker in enumerate(parsed_data["inter_motif_linkers"]):
                print(f"    {i+1}: After motif {linker['after_motif_index']+1}, length {linker['min_length']}-{linker['max_length']}")
            print(f"  Final scaffold: {parsed_data['final_scaffold']}")
        
        # Validate motifs against PDB structure
        if not validate_motifs(parsed_data["motifs"], residues_by_chain, verbose):
            print(f"Warning: Invalid motifs found in {pdb_file}")
            if not verbose:
                print("Run with --verbose for more details")
        
        # Calculate length constraints
        min_total, max_total = calculate_length_constraints(parsed_data, residues_by_chain)
        
        if verbose:
            print(f"Length constraints: MIN={min_total}, MAX={max_total}")
        
        # Generate Genie2 header
        header = generate_genie2_header(pdb_name, parsed_data, min_total, max_total)
        
        # Process residues - reorder if requested
        clean_pdb = reorder_pdb_residues(pdb_content, parsed_data, reorder_residues)
        
        # Combine header and PDB content
        combined_content = header + "\n" + clean_pdb
        
        # Write output file
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        with open(output_path, 'w') as f:
            f.write(combined_content)
        
        print(f"Successfully created {output_path}")
        return True
        
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return False

def process_csv(csv_path, pdb_dir, output_dir, verbose=False, reorder_residues=False):
    """Process multiple specifications from a CSV file"""
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            specs = list(reader)
        
        print(f"Processing {len(specs)} specifications from {csv_path}")
        success = 0
        
        for i, spec in enumerate(specs):
            pdb_file = spec.get('pdb_file')
            rf_format = spec.get('rf_format')
            
            if not pdb_file or not rf_format:
                print(f"Error: Row {i+1} missing required fields")
                continue
            
            pdb_path = os.path.join(pdb_dir, pdb_file) if pdb_dir else pdb_file
            
            if not os.path.exists(pdb_path):
                print(f"Error: PDB file not found: {pdb_path}")
                continue
            
            print(f"\nProcessing {pdb_file} with format: {rf_format}")
            if process_file(pdb_path, rf_format, output_dir, verbose, reorder_residues):
                success += 1
        
        print(f"\nProcessed {success} of {len(specs)} specifications successfully")
        return success
    
    except Exception as e:
        print(f"Error processing CSV file {csv_path}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 0

def main():
    parser = argparse.ArgumentParser(
        description='Convert RFDiffusion format to Genie2 format for SALAD compatibility'
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--pdb_file', help='Single PDB file')
    input_group.add_argument('-d', '--pdb_dir', help='Directory with PDB files')
    
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument('-i', '--input', help='RFDiffusion format or "genie2"')
    format_group.add_argument('-c', '--csv', help='CSV file with pdb_file,rf_format columns')
    
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('-r', '--reorder', action='store_true', help='Reorder residues according to motif order')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    if args.csv:
        pdb_dir = args.pdb_dir if args.pdb_dir else os.path.dirname(args.csv)
        process_csv(args.csv, pdb_dir, args.output, args.verbose, args.reorder)
    else:
        pdb_files = []
        if args.pdb_file:
            pdb_files = [args.pdb_file]
        elif args.pdb_dir:
            pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
        
        if not pdb_files:
            print("No PDB files found")
            return
        
        print(f"Processing {len(pdb_files)} files...")
        success = 0
        
        for pdb_path in pdb_files:
            if process_file(pdb_path, args.input, args.output, args.verbose, args.reorder):
                success += 1
        
        print(f"\nProcessed {success} of {len(pdb_files)} files successfully")

if __name__ == "__main__":
    main()
