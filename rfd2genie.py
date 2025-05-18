#!/usr/bin/env python3
import os
import re
import argparse
import glob
import csv
from collections import defaultdict


def parse_pdb(pdb_content):
    """Extract residue information from PDB content"""
    residues_by_chain = {}
    chain_boundaries = {}
    
    for line in pdb_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                
                if chain_id not in residues_by_chain:
                    residues_by_chain[chain_id] = set()
                    chain_boundaries[chain_id] = (float('inf'), -float('inf'))
                
                residues_by_chain[chain_id].add(res_num)
                min_res, max_res = chain_boundaries[chain_id]
                chain_boundaries[chain_id] = (min(min_res, res_num), max(max_res, res_num))
            except:
                continue
    
    return residues_by_chain, chain_boundaries

def parse_rf_format(rf_format):
    # Split into parts
    original_parts = rf_format.split("/")
    parts = [p.strip() for p in original_parts if p.strip()]
    
    # Initialize structures
    motifs = []
    inter_motif_linkers = []
    initial_scaffold = None
    final_scaffold = None
    total_length = 0
    
    # Extract initial and final scaffolds
    if parts and re.match(r'^\d+(-\d+)?$', parts[0]):
        part = parts[0]
        if "-" in part:
            min_length, max_length = map(int, part.split("-"))
        else:
            length = int(part)
            min_length = max(5, length-5)
            max_length = length+5
        
        initial_scaffold = {"min_length": min_length, "max_length": max_length}
        total_length += min_length
        parts = parts[1:]
    
    if parts and re.match(r'^\d+(-\d+)?$', parts[-1]):
        part = parts[-1]
        if "-" in part:
            min_length, max_length = map(int, part.split("-"))
        else:
            length = int(part)
            min_length = max(5, length-5)
            max_length = length+5
        
        final_scaffold = {"min_length": min_length, "max_length": max_length}
        total_length += min_length
        parts = parts[:-1]
    
    # Now parse motifs and linkers in alternating fashion
    motif_indices = []
    linker_indices = []
    
    for i, part in enumerate(parts):
        if re.match(r'^[A-Za-z]', part) or re.search(r'\[(M\d+)\]', part):
            motif_indices.append(i)
        elif re.match(r'^\d+(-\d+)?$', part):
            linker_indices.append(i)
    
    # Process each motif
    for i in motif_indices:
        part = parts[i]
        
        # Extract tag
        tag_match = re.search(r'\[(M\d+)\]', part)
        tag = tag_match.group(1) if tag_match else f"M{len(motifs) + 1}"
        
        # Clean part
        clean_part = re.sub(r'\[M\d+\]', '', part)
        
        # Extract chain and residue range
        chain_match = re.match(r'^([A-Za-z])', clean_part)
        chain = chain_match.group(1) if chain_match else "A"
        
        range_match = re.search(r'(\d+)-(\d+)', clean_part)
        if not range_match:
            raise ValueError(f"Invalid motif format: {part}")
        
        start_res = int(range_match.group(1))
        end_res = int(range_match.group(2))
        
        # Assign group (same chain = same group)
        existing_chains = {m["chain"]: m["group"] for m in motifs}
        if chain in existing_chains:
            group = existing_chains[chain]
        else:
            # New chain - assign next available group
            used_groups = set(existing_chains.values())
            next_group = chr(ord('A') + len(used_groups))
            group = next_group
        
        motifs.append({
            "chain": chain,
            "start_res": start_res,
            "end_res": end_res,
            "tag": tag,
            "group": group,
            "index_in_parts": i
        })
        
        total_length += (end_res - start_res + 1)
    
    # Process each linker
    for i in linker_indices:
        part = parts[i]
        
        # Extract min/max length
        if "-" in part:
            min_length, max_length = map(int, part.split("-"))
        else:
            length = int(part)
            min_length = max(5, length-5)
            max_length = length+5
        
        # Find which motif this linker follows
        prev_motif_index = None
        for j, motif in enumerate(motifs):
            if motif["index_in_parts"] < i:
                if prev_motif_index is None or motifs[prev_motif_index]["index_in_parts"] < motif["index_in_parts"]:
                    prev_motif_index = j
        
        # Find which motif this linker precedes
        next_motif_index = None
        for j, motif in enumerate(motifs):
            if motif["index_in_parts"] > i:
                if next_motif_index is None or motifs[next_motif_index]["index_in_parts"] > motif["index_in_parts"]:
                    next_motif_index = j
        
        # Only add as inter-motif linker if it's between two motifs
        if prev_motif_index is not None and next_motif_index is not None:
            inter_motif_linkers.append({
                "min_length": min_length,
                "max_length": max_length,
                "after_motif_index": prev_motif_index,
                "before_motif_index": next_motif_index
            })
            total_length += min_length
    
    # Clean up the motifs by removing the temporary index
    for motif in motifs:
        if "index_in_parts" in motif:
            del motif["index_in_parts"]
    
    return {
        "motifs": motifs,
        "inter_motif_linkers": inter_motif_linkers,
        "initial_scaffold": initial_scaffold,
        "final_scaffold": final_scaffold,
        "total_length": total_length
    }

def generate_genie2_header(parsed_data, pdb_name):
    output = []
    
    # Name and PDB
    name = pdb_name.split(".")[0] if "." in pdb_name else pdb_name
    output.append(f"REMARK 999 NAME   {name}")
    output.append(f"REMARK 999 PDB    {name}")
    
    # Add initial scaffold if present
    if parsed_data["initial_scaffold"]:
        scaffold = parsed_data["initial_scaffold"]
        scaffold_line = "REMARK 999 INPUT    "
        
        # Right-justify min length in columns 20-23
        min_len_str = str(scaffold["min_length"])
        scaffold_line += " " * (4 - len(min_len_str)) + min_len_str
        
        # Right-justify max length in columns 24-27
        max_len_str = str(scaffold["max_length"])
        scaffold_line += " " * (4 - len(max_len_str)) + max_len_str
        
        output.append(scaffold_line)
    
    # Process motifs and any inter-motif linkers
    for i, motif in enumerate(parsed_data["motifs"]):
        # Add motif
        motif_line = "REMARK 999 INPUT  "
        motif_line += motif["chain"]
        
        # Right-justify start residue in columns 20-23
        start_str = str(motif["start_res"])
        motif_line += " " * (4 - len(start_str)) + start_str
        
        # Right-justify end residue in columns 24-27
        end_str = str(motif["end_res"])
        motif_line += " " * (4 - len(end_str)) + end_str
        
        # Add group at column 29
        motif_line += " " + motif["group"]
        
        output.append(motif_line)
        
        # Check if there's a linker after this motif
        for linker in parsed_data["inter_motif_linkers"]:
            if linker["after_motif_index"] == i:
                # This linker follows the current motif
                linker_line = "REMARK 999 INPUT    "
                
                # Right-justify min length in columns 20-23
                min_len_str = str(linker["min_length"])
                linker_line += " " * (4 - len(min_len_str)) + min_len_str
                
                # Right-justify max length in columns 24-27
                max_len_str = str(linker["max_length"])
                linker_line += " " * (4 - len(max_len_str)) + max_len_str
                
                output.append(linker_line)
    
    # Add final scaffold if present
    if parsed_data["final_scaffold"]:
        scaffold = parsed_data["final_scaffold"]
        scaffold_line = "REMARK 999 INPUT    "
        
        # Right-justify min length in columns 20-23
        min_len_str = str(scaffold["min_length"])
        scaffold_line += " " * (4 - len(min_len_str)) + min_len_str
        
        # Right-justify max length in columns 24-27
        max_len_str = str(scaffold["max_length"])
        scaffold_line += " " * (4 - len(max_len_str)) + max_len_str
        
        output.append(scaffold_line)
    
    # Add length constraints
    min_total = max(80, parsed_data["total_length"])
    max_total = max(200, int(min_total * 1.5))
    
    min_len_line = "REMARK 999 MINIMUM TOTAL LENGTH      " + str(min_total)
    max_len_line = "REMARK 999 MAXIMUM TOTAL LENGTH      " + str(max_total)
    
    output.append(min_len_line)
    output.append(max_len_line)
    
    return "\n".join(output)

def fix_pdb_for_salad(pdb_content, verbose=False):
    lines = pdb_content.splitlines()
    
    # Separate different section types
    remark_999_lines = []
    other_remark_lines = []
    atom_lines = []
    other_lines = []
    
    for line in lines:
        if line.startswith('REMARK 999'):
            remark_999_lines.append(line)
        elif line.startswith('REMARK'):
            other_remark_lines.append(line)
        elif line.startswith(('ATOM', 'HETATM')):
            atom_lines.append(line)
        else:
            other_lines.append(line)
    
    # Process REMARK 999 lines to ensure exact column formatting
    fixed_remark_999_lines = []
    
    for line in remark_999_lines:
        if "INPUT" in line and len(line) > 16:
            parts = line.split()
            
            if len(parts) >= 5 and parts[3].isalpha() and len(parts[3]) == 1:
                # This is a motif segment line
                chain = parts[3]
                
                try:
                    start_res = int(parts[4])
                    end_res = int(parts[5])
                    group = parts[6] if len(parts) > 6 else chain
                    
                    # Reconstruct with exact column positioning
                    new_line = "REMARK 999 INPUT  "  # Columns 1-16
                    new_line += chain                 # Column 19
                    
                    # Right-justify start residue in columns 20-23
                    start_res_str = str(start_res)
                    new_line += " " * (4 - len(start_res_str)) + start_res_str
                    
                    # Right-justify end residue in columns 24-27
                    end_res_str = str(end_res)
                    new_line += " " * (4 - len(end_res_str)) + end_res_str
                    
                    # Column 29: Group
                    new_line += " " + group
                    
                    fixed_remark_999_lines.append(new_line)
                except (ValueError, IndexError):
                    # If we can't parse it properly, keep it as is
                    fixed_remark_999_lines.append(line)
                    
            elif len(parts) >= 4 and parts[3].isdigit():
                # This is a scaffold segment line
                try:
                    min_len = int(parts[3])
                    max_len = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else min_len
                    
                    # Reconstruct with exact column positioning
                    new_line = "REMARK 999 INPUT    "  # Columns 1-16
                    
                    # Right-justify min length in columns 20-23
                    min_len_str = str(min_len)
                    new_line += " " * (4 - len(min_len_str)) + min_len_str
                    
                    # Right-justify max length in columns 24-27
                    max_len_str = str(max_len)
                    new_line += " " * (4 - len(max_len_str)) + max_len_str
                    
                    fixed_remark_999_lines.append(new_line)
                except (ValueError, IndexError):
                    # If we can't parse it properly, keep it as is
                    fixed_remark_999_lines.append(line)
                    
            else:
                fixed_remark_999_lines.append(line)
        
        elif "MINIMUM TOTAL LENGTH" in line or "MAXIMUM TOTAL LENGTH" in line:
            # Length constraint line
            try:
                parts = line.split()
                length_type = "MINIMUM" if "MINIMUM" in line else "MAXIMUM"
                value = int(parts[-1])
                
                # Reconstruct with exact column positioning
                new_line = f"REMARK 999 {length_type} TOTAL LENGTH      {value}"
                
                fixed_remark_999_lines.append(new_line)
            except (ValueError, IndexError):
                # If we can't parse it properly, keep it as is
                fixed_remark_999_lines.append(line)
                
        else:
            # Keep other REMARK 999 lines as is
            fixed_remark_999_lines.append(line)
    
    # Reassemble the file with fixed formatting
    result = []
    result.extend(fixed_remark_999_lines)
    result.extend(other_remark_lines)
    result.extend(atom_lines)
    result.extend(other_lines)
    
    if verbose:
        print(f"Fixed {len(fixed_remark_999_lines)} REMARK 999 lines for exact column positioning")
    
    return "\n".join(result)


def process_file(pdb_path, rf_format, output_dir, verbose=False):
    try:
        # Read PDB content
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        pdb_file = os.path.basename(pdb_path)
        pdb_name = os.path.splitext(pdb_file)[0]
        
        # Parse PDB structure
        residues_by_chain, chain_boundaries = parse_pdb(pdb_content)
        
        if verbose:
            print(f"\nChain boundaries in {pdb_file}:")
            for chain, (min_res, max_res) in chain_boundaries.items():
                print(f"Chain {chain}: {min_res} - {max_res} ({len(residues_by_chain[chain])} residues)")
        
        # Handle already Genie2 formatted files
        if rf_format.lower() == "genie2":
            fixed_content = fix_pdb_for_salad(pdb_content, verbose)
            
            # Write to output file
            output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
            with open(output_path, 'w') as f:
                f.write(fixed_content)
            
            print(f"Successfully processed existing Genie2 file: {output_path}")
            return True
        
        # Process RFDiffusion format
        parsed_data = parse_rf_format(rf_format)
        
        if verbose:
            print("\nParsed RFDiffusion format:")
            
            # Show motifs
            print(f"  Motifs ({len(parsed_data['motifs'])}):")
            for i, motif in enumerate(parsed_data["motifs"]):
                print(f"    Motif {i+1}: Chain {motif['chain']}, "
                      f"Residues {motif['start_res']}-{motif['end_res']}, "
                      f"Group {motif['group']}")
            
            # Show linkers
            print(f"  Linkers/Extensions ({len(parsed_data['linkers'])}):")
            for i, linker in enumerate(parsed_data["linkers"]):
                position = "N-terminal" if i == 0 and i < len(parsed_data["motifs"]) else \
                           "C-terminal" if i >= len(parsed_data["motifs"]) else \
                           f"Between motifs {i} and {i+1}"
                print(f"    Linker {i+1} ({position}): "
                      f"{linker['min_length']}-{linker['max_length']} residues")
        
        # Generate Genie2 header
        genie2_header = generate_genie2_header(parsed_data, pdb_name)
        
        # Validate motifs against PDB structure
        for i, motif in enumerate(parsed_data["motifs"]):
            chain = motif["chain"]
            start_res = motif["start_res"]
            end_res = motif["end_res"]
            
            if chain not in residues_by_chain:
                print(f"Warning: Chain {chain} not found in PDB")
                continue
                
            missing = []
            for res in range(start_res, end_res + 1):
                if res not in residues_by_chain[chain]:
                    missing.append(res)
            
            if missing:
                print(f"Warning: Motif {i+1} missing residues: {missing[:5]}..." +
                      ("" if len(missing) <= 5 else f" and {len(missing)-5} more"))
                if verbose:
                    print(f"  Available residues in chain {chain}: " +
                          f"{sorted(list(residues_by_chain[chain]))[:10]}..." +
                          ("" if len(residues_by_chain[chain]) <= 10 else " and more"))
        
        # Remove any existing REMARK 999 lines from the original PDB
        clean_pdb_lines = [line for line in pdb_content.splitlines() 
                          if not line.startswith('REMARK 999')]
        clean_pdb = "\n".join(clean_pdb_lines)
        
        # Combine header and cleaned content
        combined_content = genie2_header + "\n" + clean_pdb
        
        # Write output
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        with open(output_path, 'w') as f:
            f.write(combined_content)
        
        print(f"Successfully created {output_path}")
        
        # If verbose, show the generated header
        if verbose:
            print("\nGenerated Genie2 header:")
            for line in genie2_header.splitlines():
                print(f"  {line}")
        
        return True
    
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return False


def process_csv(csv_path, pdb_dir, output_dir, verbose=False):
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
                print(f"Error: Row {i+1} missing required fields (pdb_file, rf_format)")
                continue
            
            pdb_path = os.path.join(pdb_dir, pdb_file) if pdb_dir else pdb_file
            
            if not os.path.exists(pdb_path):
                print(f"Error: PDB file not found: {pdb_path}")
                continue
            
            print(f"\nProcessing {pdb_file} with format: {rf_format}")
            if process_file(pdb_path, rf_format, output_dir, verbose):
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
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Convert RFDiffusion format to Genie2 format for SALAD compatibility',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--pdb_file', help='Single PDB file')
    input_group.add_argument('-d', '--pdb_dir', help='Directory with PDB files')
    
    # Format options
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument('-i', '--input', help='RFDiffusion format (e.g., "A1-80[M1]/30/[M2]B81-100") or "genie2"')
    format_group.add_argument('-c', '--csv', help='CSV file with pdb_file,rf_format columns')
    
    # Output and other options
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Process files
    if args.csv:
        # Process multiple files with different formats from CSV
        pdb_dir = args.pdb_dir if args.pdb_dir else os.path.dirname(args.csv)
        process_csv(args.csv, pdb_dir, args.output, args.verbose)
    else:
        # Process single file or multiple files with same format
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
            if process_file(pdb_path, args.input, args.output, args.verbose):
                success += 1
        
        print(f"\nProcessed {success} of {len(pdb_files)} files successfully")
    
    print("Files are ready for SALAD 🧬")

if __name__ == "__main__":
    main()
