#!/usr/bin/env python3

import os
import re
import argparse
import csv
import glob
from pathlib import Path
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
            error_msg = f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Chain not found in PDB"
            errors.append(error_msg)
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
                
            error_msg = f"Motif {i+1} (chain {chain}, residues {start_res}-{end_res}): Missing {missing_str}"
            errors.append(error_msg)
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


class MotifScaffoldConverter:
    def __init__(self, pdb_name="input_pdb", add_terminal_scaffolds=False, 
                 min_scaffold_length=5, max_scaffold_length=20,
                 min_total_length_factor=1.0, max_total_length_factor=1.5,
                 auto_fix_motifs=True):
        self.pdb_name = pdb_name
        self.add_terminal_scaffolds = add_terminal_scaffolds
        self.min_scaffold_length = min_scaffold_length
        self.max_scaffold_length = max_scaffold_length
        self.min_total_length_factor = min_total_length_factor
        self.max_total_length_factor = max_total_length_factor
        self.auto_fix_motifs = auto_fix_motifs
    
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
            output.append(f"REMARK 999 INPUT     {self.min_scaffold_length:2d}  {self.max_scaffold_length:2d}")
        
        # Add motifs and linkers
        for i, motif in enumerate(parsed_data["motifs"]):
            # Format according to specifications: proper spacing and right-justified numbers
            output.append(
                f"REMARK 999 INPUT  {motif['chain']} {motif['start_res']:3d} {motif['end_res']:3d} {motif['group']}"
            )
            
            # Add linker after motif if available
            if i < len(parsed_data["linkers"]):
                linker = parsed_data["linkers"][i]
                output.append(
                    f"REMARK 999 INPUT     {linker['min_length']:2d}  {linker['max_length']:2d}"
                )
        
        # Add C-terminal scaffold if requested
        if self.add_terminal_scaffolds:
            output.append(f"REMARK 999 INPUT     {self.min_scaffold_length:2d}  {self.max_scaffold_length:2d}")
        
        # Calculate total length constraints
        min_total_length = max(80, int(parsed_data["total_min_length"] * self.min_total_length_factor))
        max_total_length = max(200, int(parsed_data["total_max_length"] * self.max_total_length_factor))
        
        # Ensure length limits make sense
        if min_total_length > max_total_length:
            print(f"Warning: Minimum length ({min_total_length}) exceeds maximum length ({max_total_length}). Adjusting maximum.")
            max_total_length = min_total_length + 100
        
        # Print a warning if total length is very large
        if max_total_length > 500:
            print(f"Warning: Maximum total length is {max_total_length}, which may result in long computation times.")
        
        # Format according to specifications: left-justified numbers
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
        if line.startswith('REMARK 999 INPUT') and len(line.split()) >= 4:
            parts = line.split()
            if len(parts) >= 6 and len(parts[3]) == 1 and parts[3].isalpha():
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


def check_and_fix_constraint_consistency(pdb_content):
    """Check and fix constraint consistency for SALAD"""
    lines = pdb_content.splitlines()
    
    # Parse residue information from PDB
    residues_by_chain, chain_boundaries, _ = parse_pdb_structure(pdb_content)
    
    # Parse motif definitions
    genie2_data = parse_genie2_motifs(pdb_content)
    
    # Check if motifs need fixing
    valid, errors = validate_motifs(genie2_data["motifs"], residues_by_chain)
    
    # If everything is valid, no need to fix
    if valid:
        print("All motifs are valid, no fixes needed")
        return pdb_content
    
    print("Found motif issues, applying fixes...")
    for error in errors:
        print(f"  {error}")
    
    # Fix motifs
    fixed_motifs = fix_motifs(genie2_data["motifs"], residues_by_chain)
    
    # If no valid motifs remain after fixing, return original content
    if not fixed_motifs:
        print("Error: No valid motifs could be created after fixing. Returning original content.")
        return pdb_content
    
    # Rebuild the header with fixed motifs
    header_lines = []
    atom_lines = []
    footer_lines = []
    name_line = None
    pdb_line = None
    atom_section_started = False
    atom_section_ended = False
    
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            atom_section_started = True
            atom_lines.append(line)
        elif atom_section_started and not atom_section_ended:
            if line.startswith(("TER", "END", "CONECT")):
                atom_lines.append(line)
            else:
                atom_section_ended = True
                footer_lines.append(line)
        elif atom_section_ended:
            footer_lines.append(line)
        elif line.startswith("REMARK 999 NAME"):
            name_line = line
        elif line.startswith("REMARK 999 PDB"):
            pdb_line = line
        elif not line.startswith(("REMARK 999 INPUT", "REMARK 999 MINIMUM TOTAL LENGTH", "REMARK 999 MAXIMUM TOTAL LENGTH")):
            header_lines.append(line)
    
    # Rebuild the file with fixed headers
    new_lines = []
    
    # Add name and PDB lines first
    if name_line:
        new_lines.append(name_line)
    if pdb_line:
        new_lines.append(pdb_line)
    
    # Add fixed motifs
    for motif in fixed_motifs:
        new_lines.append(f"REMARK 999 INPUT  {motif['chain']} {motif['start_res']:3d} {motif['end_res']:3d} {motif['group']}")
        
    # Add linkers
    for i, linker in enumerate(genie2_data["linkers"]):
        if i < len(fixed_motifs) - 1:  # Only add linkers between motifs
            new_lines.append(f"REMARK 999 INPUT     {linker['min_length']:2d}  {linker['max_length']:2d}")
    
    # Calculate and add new total length
    min_total = sum(motif['end_res'] - motif['start_res'] + 1 for motif in fixed_motifs)
    max_total = min_total + 100  # Add some buffer
    
    if genie2_data["linkers"]:
        linker_min = sum(linker['min_length'] for linker in genie2_data["linkers"][:len(fixed_motifs)-1])
        linker_max = sum(linker['max_length'] for linker in genie2_data["linkers"][:len(fixed_motifs)-1])
        min_total += linker_min
        max_total = min_total + linker_max
    
    # Ensure minimum viable protein length
    min_total = max(80, min_total)
    max_total = max(200, max_total)
    
    # Use original values if present
    if "min_total" in genie2_data and genie2_data["min_total"] >= min_total:
        min_total = genie2_data["min_total"]
    if "max_total" in genie2_data and genie2_data["max_total"] >= max_total:
        max_total = genie2_data["max_total"]
    
    # Ensure max is greater than min
    if max_total <= min_total:
        max_total = min_total + 100
    
    new_lines.append(f"REMARK 999 MINIMUM TOTAL LENGTH      {min_total}")
    new_lines.append(f"REMARK 999 MAXIMUM TOTAL LENGTH      {max_total}")
    
    # Add the rest of the content
    new_lines.extend(header_lines)
    new_lines.extend(atom_lines)
    new_lines.extend(footer_lines)
    
    return "\n".join(new_lines)


def read_specs_from_csv(csv_file):
    """Read mapping of PDB files to RF formats from CSV"""
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

def process_single_pdb(pdb_path, rf_format, output_dir, options):
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
        
        # Check if this is already a Genie2 formatted file
        if rf_format.lower() == "genie2" or rf_format.lower() == "auto":
            # Extract motif definitions from the PDB file
            genie2_data = parse_genie2_motifs(pdb_content)
            
            if genie2_data["motifs"]:
                print(f"\nDetected Genie2 format with {len(genie2_data['motifs'])} motifs")
                for i, motif in enumerate(genie2_data["motifs"]):
                    print(f"Motif {i+1}: Chain {motif['chain']}, Residues {motif['start_res']}-{motif['end_res']}")
                
                # Check and fix constraint consistency
                fixed_content = check_and_fix_constraint_consistency(pdb_content)
                
                # Apply SALAD fix
                verbose = options.get('verbose', False)
                print(f"\nApplying SALAD compatibility fixes...")
                fixed_content = fix_pdb_for_salad(
                    fixed_content, 
                    fix_residue_numbering=True,
                    verbose=verbose
                )
                
                # Write to output file
                output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
                with open(output_path, 'w') as f:
                    f.write(fixed_content)
                
                print(f"Successfully created {output_path}")
                return (pdb_path, True)
            elif rf_format.lower() == "auto":
                print("No Genie2 motif definitions found, trying RFDiffusion format...")
            else:
                print("Error: No Genie2 motif definitions found in the PDB file")
                return (pdb_path, False)
        
        # Initialize converter with options
        converter = MotifScaffoldConverter(
            pdb_name,
            add_terminal_scaffolds=options.get('add_terminal_scaffolds', False),
            min_scaffold_length=options.get('min_scaffold_length', 5),
            max_scaffold_length=options.get('max_scaffold_length', 20),
            min_total_length_factor=options.get('min_total_length_factor', 1.0),
            max_total_length_factor=options.get('max_total_length_factor', 1.5)
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
            
            # Auto-fix motifs if requested
            if options.get('auto_fix_motifs', True):
                print("Fixing invalid motif definitions...")
                parsed_data["motifs"] = fix_motifs(parsed_data["motifs"], residues_by_chain)
                
                # Regenerate the header with fixed motifs
                converter = MotifScaffoldConverter(
                    pdb_name,
                    add_terminal_scaffolds=options.get('add_terminal_scaffolds', False),
                    min_scaffold_length=options.get('min_scaffold_length', 5),
                    max_scaffold_length=options.get('max_scaffold_length', 20),
                    min_total_length_factor=options.get('min_total_length_factor', 1.0),
                    max_total_length_factor=options.get('max_total_length_factor', 1.5)
                )
                
                # Create a mock RF format string from the fixed motifs and linkers
                mock_rf_format = ""
                for i, motif in enumerate(parsed_data["motifs"]):
                    if i > 0:
                        # Add linker
                        if i-1 < len(parsed_data["linkers"]):
                            linker = parsed_data["linkers"][i-1]
                            mock_rf_format += f"/{linker['min_length']}-{linker['max_length']}/"
                        else:
                            mock_rf_format += "/10-20/"
                    
                    # Add motif
                    mock_rf_format += f"{motif['chain']}{motif['start_res']}-{motif['end_res']}[{motif['tag']}]"
                
                genie2_header, _ = converter.convert(mock_rf_format)
            else:
                print(f"Error: Validation failed for {pdb_file}. Use --auto_fix_motifs to automatically fix issues.")
                print("\nAvailable residue ranges:")
                for chain, (min_res, max_res) in chain_boundaries.items():
                    print(f"Chain {chain}: {min_res} - {max_res}")
                return (pdb_path, False)
        
        # Combine header and PDB content
        combined_content = genie2_header + "\n" + pdb_content
        
        # Apply SALAD fix
        verbose = options.get('verbose', False)
        print(f"\nApplying SALAD compatibility fixes...")
        fixed_content = fix_pdb_for_salad(
            combined_content, 
            fix_residue_numbering=True,
            verbose=verbose
        )
        
        # Write to output file
        output_path = os.path.join(output_dir, f"{pdb_name}_genie2.pdb")
        with open(output_path, 'w') as f:
            f.write(fixed_content)
        
        print(f"Successfully created {output_path}")
        return (pdb_path, True)
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return (pdb_path, False)
    
def process_pdb_files(pdb_sources, specs, output_dir, default_rf_format=None, options=None):
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
    
    print(f"Processing {len(pdb_files)} files...")
    success_count = 0
    failure_count = 0
    failed_files = []
    
    for pdb_path in pdb_files:
        pdb_file = os.path.basename(pdb_path)
        rf_format = specs.get(pdb_file, default_rf_format)
        
        if rf_format:
            result = process_single_pdb(pdb_path, rf_format, output_dir, options)
            
            if result[1]:
                success_count += 1
            else:
                failure_count += 1
                failed_files.append(pdb_path)
        else:
            print(f"Skipping {pdb_file} - no RF format specified")
            failed_files.append(pdb_path)
    
    return (success_count, failure_count, failed_files)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert RFDiffusion format to Genie2 format with SALAD compatibility')
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--pdb_file', '-f', help='Path to a single PDB file')
    input_group.add_argument('--pdb_dir', '-d', help='Path to a directory containing PDB files')
    
    spec_group = parser.add_mutually_exclusive_group(required=True)
    spec_group.add_argument('--input', '-i', help='RFDiffusion format string or "genie2" to process existing Genie2 PDB')
    spec_group.add_argument('--csv', '-c', help='CSV file mapping PDB files to RFDiffusion formats')
    
    parser.add_argument('--output', '-o', required=True, help='Output directory for processed PDB files')
    
    parser.add_argument('--add_terminal_scaffolds', action='store_true', help='Add scaffolds at the beginning and end')
    parser.add_argument('--min_scaffold', type=int, default=5, help='Minimum scaffold length')
    parser.add_argument('--max_scaffold', type=int, default=20, help='Maximum scaffold length')
    parser.add_argument('--min_factor', type=float, default=1.0, help='Minimum total length factor')
    parser.add_argument('--max_factor', type=float, default=1.5, help='Maximum total length factor')
    
    parser.add_argument('--auto_fix_motifs', action='store_true', default=True, help='Automatically fix motif definitions')
    parser.add_argument('--no_fix', action='store_true', help='Disable automatic fixing of motif definitions')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print detailed information')
    
    args = parser.parse_args()
    
    # Read specifications from CSV if provided
    specs = {}
    if args.csv:
        specs = read_specs_from_csv(args.csv)
    
    # Prepare options
    options = {
        'add_terminal_scaffolds': args.add_terminal_scaffolds,
        'min_scaffold_length': args.min_scaffold,
        'max_scaffold_length': args.max_scaffold,
        'min_total_length_factor': args.min_factor,
        'max_total_length_factor': args.max_factor,
        'auto_fix_motifs': args.auto_fix_motifs and not args.no_fix,
        'verbose': args.verbose
    }
    
    # Determine PDB source
    pdb_source = args.pdb_file if args.pdb_file else args.pdb_dir
    
    # Process files
    print(f"\n--- Starting SALAD Motif Converter ---")
    print(f"Processing PDB source: {pdb_source}")
    
    success, failure, failed_files = process_pdb_files(
        pdb_source, 
        specs, 
        args.output, 
        default_rf_format=args.input,
        options=options
    )
    
    # Print summary
    print(f"\nProcessing complete")
    print(f"Successfully processed: {success} files")
    print(f"Failed to process: {failure} files")
    
    if failed_files:
        print("\nFailed files:")
        for file in failed_files:
            print(f"- {os.path.basename(file)}")
    
    print("\nOutput files are ready for SALAD 🧬")

if __name__ == "__main__":
    main()
