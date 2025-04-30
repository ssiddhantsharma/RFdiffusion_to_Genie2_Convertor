# RFDiffusion to Genie2 Converter

A Python tool to convert RFDiffusion-style motif specifications to [Genie2](https://github.com/aqlaboratory/genie2) format for protein design. This tool is particularly useful for preparing input files for [SALAD's](https://github.com/mjendrusch/salad) multi-motif scaffolding and heterodimer design.

## Installation

```bash
git clone https://github.com/ssiddhantsharma/rfdiffusion-to-genie2.git
cd rfdiffusion-to-genie2
chmod +x rfd2genie.py
```

## Basic Usage

### Converting a Single PDB File

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir
```

### Converting Multiple PDB Files (Same Specification)

```bash
python rfd2genie.py --pdb_dir pdb_directory --input "A1-80[M1]/30/[M2]B81-100" --output output_dir
```

## RFDiffusion Format Syntax

The input format follows this pattern: `[Chain][Start]-[End][Motif Tag]/[Linker]/[Chain][Start]-[End][Motif Tag]`

### Examples

- **Multiple chains with linker:** `A1-80[M1]/30/[M2]B81-100`
- **Multiple motifs on same chain:** `A1-29[M1]/30/[M2]A39-50/40/[M3]A2-49`
- **Heterodimer (no linker):** `A1-30[M1]/B89-95[M2]`
- **Heterotrimer:** `A1-30[M1]/B10-40[M2]/C5-25[M3]`

### Components

- **Chain:** Single letter (e.g., A, B)
- **Residue Range:** Start-End (e.g., 1-80)
- **Motif Tag** (optional): [M1], [M2], etc.
- **Linker:** Numeric value representing amino acid length

## Group Assignment for Multi-Chain Complexes

The converter intelligently assigns groups based on chain identifiers:

- Motifs from the same chain are assigned to the same group
- Motifs from different chains are assigned to different groups
- First chain gets group "A", all other chains get group "B"

This is particularly useful for heterodimer, heterotrimer, or other multi-chain complex designs where you want to maintain separate entities.

## Advanced Options

### Terminal Scaffolds

Add flexible regions at the N and C termini of your design to improve stability and folding:

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --add_terminal_scaffolds
```

Terminal scaffolds are useful for:
- Preventing edge effects at protein termini
- Improving solubility of the designed protein
- Adding buffer regions that can be removed later if needed

### Scaffold Length

Control the minimum and maximum length of terminal scaffolds (default: 5-20 residues):

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --min_scaffold 10 --max_scaffold 30
```

Adjusting scaffold length helps:
- Create longer terminal regions for enhanced stability
- Set tighter constraints when working with size-limited designs
- Balance flexibility with structural integrity

### Total Length Factors

Adjust the minimum and maximum total sequence length calculation by multiplying the sum of motif and linker lengths:

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --min_factor 1.2 --max_factor 2.0
```

Length factors are crucial for:
- Ensuring sufficient sequence length for proper folding
- Allowing the design algorithm enough flexibility to create novel structures
- Constraining total protein size for expression or application requirements

### Processing Multiple Files with Different Specifications

Create a CSV file (e.g., `specs.csv`):
```
pdb_file,rf_format
protein1.pdb,A1-80[M1]/30/[M2]B81-100
protein2.pdb,A1-29[M1]/30/[M2]A39-50
```

Then run:
```bash
python rfd2genie.py --pdb_dir pdb_directory --csv specs.csv --output output_dir
```

Alternatively, use a JSON file:
```bash
python rfd2genie.py --pdb_dir pdb_directory --json specs.json --output output_dir
```

### Automatic Motif Numbering Correction

If your PDB file has residue numbering that doesn't match your motif definition, you can use the auto-correct feature:

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --auto_correct
```

Auto-correction only activates when:
- A motif is completely outside the available residue range in a chain
- There's a clear systematic shift in numbering (e.g., off by +1 or -1)

## SALAD Compatibility

The converter now automatically makes output files compatible with SALAD by:

1. Ensuring residue numbering is continuous (1,2,3,...) within each chain
2. Expanding motif definitions to cover full chains while preserving which parts should be fixed vs. designable
3. Adjusting array shapes to avoid SALAD's size mismatch errors

These compatibility fixes are applied automatically for all outputs. This eliminates common errors when using the generated files with SALAD's multi-motif scaffolding.

## Integration with SALAD

The converter is designed to work seamlessly with SALAD's multi-motif scaffolding:

1. Convert your PDB file to Genie2 format:
   ```bash
   python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output genie2_pdbs
   ```

2. Use the converted file with SALAD:
   ```bash
   python salad/training/eval_motif_benchmark.py \
       --config multimotif_vp \
       --params params/multimotif_vp-200k.jax \
       --out_path designed_proteins/ \
       --num_steps 500 --out_steps 400 --prev_threshold 0.8 \
       --num_designs 10 --timescale_pos "cosine(t)" \
       --template genie2_pdbs/your_protein_genie2.pdb
   ```

## Key Features

- **Built-in SALAD Compatibility**: All output files are automatically compatible with SALAD
- **Chain-Based Group Assignment**: Automatically assigns different groups for different chains (ideal for heterodimers)
- **Strict Residue Validation**: Ensures that specified residues exist in the PDB file
- **Full-Chain Motif Support**: Properly handles the chain boundaries for SALAD compatibility
- **Multi-Chain Support**: Handles motifs spanning multiple chains
- **Flexible Linkers**: Supports variable-length linkers between motifs
- **HETATM Support**: Handles non-standard residues in motif definitions
- **Parallel Processing**: Efficiently processes multiple PDB files simultaneously
