#!/usr/bin/env python3
import os
import re
import argparse
import glob
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
    """Parse RFDiffusion format into motifs and linkers"""
    motifs = []
    linkers = []
    total_length = 0
    
    parts = rf_format.split("/")
    
    for i, part in enumerate(parts):
        if not part.strip():
            continue
            
        # Check if this is a linker
        if re.match(r'^\d+$', part):
            # Single number format: 30
            length = int(part)
            linkers.append({"min_length": max(5, length-5), "max_length": length+5})
            total_length += length
            continue
        
        # Check if this is a linker range
        linker_match = re.match(r'^(\d+)-(\d+)$', part)
        if linker_match and not re.match(r'^[A-Za-z]', part):
            # Range format: 30-40
            min_length = int(linker_match.group(1))
            max_length = int(linker_match.group(2))
            linkers.append({"min_length": min_length, "max_length": max_length})
            total_length += min_length
            continue
        
        # This is a motif
        # Extract motif tag [M1], [M2], etc.
        tag_match = re.search(r'\[(M\d+)\]', part)
        tag = tag_match.group(1) if tag_match else f"M{len(motifs) + 1}"
        
        # Remove tag for further processing
        clean_part = re.sub(r'\[M\d+\]', '', part)
        
        # Extract chain and residue range
        chain_match = re.match(r'^([A-Za-z])', clean_part)
        chain = chain_match.group(1) if chain_match else "A"
        
        range_match = re.search(r'(\d+)-(\d+)', clean_part)
        if not range_match:
            raise ValueError(f"Invalid motif format: {part}")
        
        start_res = int(range_match.group(1))
        end_res = int(range_match.group(2))
        
        # Determine group (A for first chain, B for others)
        if len(motifs) == 0:
            group = "A"
        else:
            first_chain = motifs[0]["chain"]
            group = "A" if chain == first_chain else "B"
        
        motifs.append({
            "chain": chain,
            "start_res": start_res,
            "end_res": end_res,
            "tag": tag,
            "group": group
        })
        
        total_length += (end_res - start_res + 1)
    
    return {"motifs": motifs, "linkers": linkers, "total_length": total_length}


def generate_genie2_header(parsed_data, pdb_name):
    """Generate properly formatted Genie2 header with exact column spacing"""
    output = []
    
    # Name and PDB
    name = pdb_name.split(".")[0] if "." in pdb_name else pdb_name
    output.append(f"REMARK 999 NAME   {name}_motifs")
    output.append(f"REMARK 999 PDB    {name}")
    
    # Motifs and linkers
    for i, motif in enumerate(parsed_data["motifs"]):
        # Motif line with exact column spacing for SALAD compatibility
        output.append(
            f"REMARK 999 INPUT  {motif['chain']} {motif['start_res']:3d} {motif['end_res']:3d} {motif['group']}"
        )
        
        # Add linker after motif if available
        if i < len(parsed_data["linkers"]):
            linker = parsed_data["linkers"][i]
            # Scaffold segment with exact column spacing
            output.append(
                f"REMARK 999 INPUT    {linker['min_length']:2d} {linker['max_length']:2d}"
            )
    
    # Length constraints
    min_total = max(80, parsed_data["total_length"])
    max_total = max(200, int(min_total * 1.5))
    
    output.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total}")
    output.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total}")
    
    return "\n".join(output)


def fix_pdb_for_salad(pdb_content, verbose=False):
    """
    Fix PDB for SALAD compatibility:
    1. Ensure motifs cover full chains
    2. Renumber residues continuously
    3. Set proper formatting
    """
    lines = pdb_content.splitlines()
    
    # Separate headers, atoms, and other content
    headers = [line for line in lines if line.startswith('REMARK')]
    atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
    other_lines = [line for line in lines if not line.startswith(('REMARK', 'ATOM', 'HETATM'))]
    
    # Extract chain and residue info
    chains = {}
    for line in atom_lines:
        chain_id = line[21]
        if chain_id not in chains:
            chains[chain_id] = {'residues': set(), 'mapping': {}}
        
        res_num = int(line[22:26])
        chains[chain_id]['residues'].add(res_num)
    
    # Create continuous residue numbering
    for chain_id, chain_data in chains.items():
        sorted_residues = sorted(chain_data['residues'])
        for i, res_num in enumerate(sorted_residues, 1):
            chain_data['mapping'][res_num] = i
    
    # Extract existing motif and header info
    motif_defs = []
    name_line = None
    pdb_line = None
    linker_lines = []
    
    for line in headers:
        if line.startswith('REMARK 999 NAME'):
            name_line = line
        elif line.startswith('REMARK 999 PDB'):
            pdb_line = line
        elif line.startswith('REMARK 999 INPUT') and len(line.split()) >= 6:
            parts = line.split()
            if len(parts[3]) == 1 and parts[3].isalpha():
                try:
                    chain_id = parts[3]
                    group = parts[6] if len(parts) > 6 else "A"
                    motif_defs.append({'chain': chain_id, 'group': group})
                except:
                    pass
        elif line.startswith('REMARK 999 INPUT') and len(line.split()) >= 5:
            parts = line.split()
            if not (len(parts[3]) == 1 and parts[3].isalpha()):
                linker_lines.append(line)
    
    # Create new motifs for full chain coverage
    new_lines = []
    if name_line: new_lines.append(name_line)
    if pdb_line: new_lines.append(pdb_line)
    
    # Create/update motifs for all chains
    processed_chains = set(m['chain'] for m in motif_defs)
    
    # First handle existing motifs
    for motif in motif_defs:
        chain_id = motif['chain']
        if chain_id in chains:
            new_lines.append(f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chains[chain_id]['residues']):3d} {motif['group']}")
    
    # Then add any missing chains
    for chain_id in chains:
        if chain_id not in processed_chains:
            group = "B" if motif_defs and motif_defs[0]['group'] == "A" else "A"
            new_lines.append(f"REMARK 999 INPUT  {chain_id} {1:3d} {len(chains[chain_id]['residues']):3d} {group}")
    
    # Add linker lines
    new_lines.extend(linker_lines)
    
    # Add length constraints if missing
    if not any("MINIMUM TOTAL LENGTH" in line for line in new_lines):
        total_residues = sum(len(chain_data['residues']) for chain_data in chains.values())
        min_total = max(80, total_residues)
        max_total = max(200, int(min_total * 1.5))
        new_lines.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total}")
        new_lines.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total}")
    
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
    
    # Combine everything
    new_lines.extend(new_atoms)
    new_lines.extend(other_lines)
    
    return "\n".join(new_lines)


def process_file(pdb_path, rf_format, output_dir, verbose=False):
    """Process a single PDB file with the specified format"""
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
                print(f"Warning: Motif {i+1} missing residues: {missing[:5]}...")
        
        # Combine header and content
        combined_content = genie2_header + "\n" + pdb_content
        
        # Fix for SALAD compatibility
        fixed_content = fix_pdb_for_salad(combined_content, verbose)
        
        # Write output
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        with open(output_path, 'w') as f:
            f.write(fixed_content)
        
        print(f"Successfully created {output_path}")
        return True
    
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return False


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Convert RFDiffusion to Genie2 format for SALAD compatibility')
    
    # Required arguments
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--pdb_file', help='Single PDB file')
    input_group.add_argument('-d', '--pdb_dir', help='Directory with PDB files')
    
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument('-i', '--input', help='RFDiffusion format (e.g., "A1-80[M1]/30-40/B81-100[M2]") or "genie2"')
    
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Process files
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
