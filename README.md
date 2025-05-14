# RFDiffusion to Genie2 Converter

A Python tool to convert RFDiffusion-style motif specifications to [Genie2](https://github.com/aqlaboratory/genie2) format for protein design. This tool is particularly useful for preparing input files for both Genie2's motif scaffolding and [SALAD's](https://github.com/mjendrusch/salad) multi-motif scaffolding and heterodimer design.

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
- **Linker:** Numeric value representing amino acid length or a range (e.g., 30 or 30-40)

## Format Conversion: RFDiffusion to Genie2

The script converts RFDiffusion format into Genie2's precise column-specific format:

| Specification | Columns | Data | Justification | Data Type |
|---------------|---------|------|--------------|-----------|
| Motif segment | 1-16 | "REMARK 999 INPUT" | - | string |
| | 19 | Chain index of motif segment in the PDB file | - | string |
| | 20-23 | Starting residue index of motif segment | right | int |
| | 24-27 | Ending residue index of motif segment | right | int |
| | 29 | Motif group that the segment belongs to | - | string |
| Scaffold segment | 1-16 | "REMARK 999 INPUT" | - | string |
| | 20-23 | Minimum length of scaffold segment | right | int |
| | 24-27 | Maximum length of scaffold segment | right | int |
| Minimum sequence length | 1-31 | "REMARK 999 MINIMUM TOTAL LENGTH" | - | string |
| | 38-40 | Minimum sequence length | left | int |
| Maximum sequence length | 1-31 | "REMARK 999 MAXIMUM TOTAL LENGTH" | - | string |
| | 38-40 | Maximum sequence length | left | int |

For example, the RFDiffusion format `A1-80[M1]/30/[M2]B81-100` converts to:

```
REMARK 999 NAME   protein_name_motifs
REMARK 999 PDB    protein_name
REMARK 999 INPUT  A   1  80 A
REMARK 999 INPUT     30  35
REMARK 999 INPUT  B  81 100 B
REMARK 999 MINIMUM TOTAL LENGTH      116
REMARK 999 MAXIMUM TOTAL LENGTH      215
```

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

## Using Converted Files

### Running with Genie2

After converting your files, you can directly use them with Genie2 for motif scaffolding:

1. Place your converted PDB files in the appropriate Genie2 data directory (e.g., `data/multimotifs` for multi-motif scaffolding)

2. Run Genie2's scaffolding algorithm:
   ```bash
   python genie/sample_scaffold.py --name base --epoch 30 --scale 0.4 --outdir results/my_designs --datadir my_converted_pdbs --num_samples 100
   ```

   Key parameters:
   - `--name`: Model name (default: "base" for Genie2's base model)
   - `--epoch`: Model epoch (e.g., 30)
   - `--scale`: Sampling noise scale (between 0 and 1)
   - `--outdir`: Output directory for results
   - `--datadir`: Directory containing your converted PDB files
   - `--num_samples`: Number of designs to generate per motif (default: 100)

3. Your results will be stored in the specified output directory.

#### Multi-Motif Scaffolding Example

For multi-motif scaffolding with Genie2, run:

```bash
# First convert your PDB files
python rfd2genie.py --pdb_dir my_motifs --input "A1-80[M1]/30/[M2]B81-100" --output genie2_pdbs

# Then run Genie2 on the converted files
python genie/sample_scaffold.py --name base --epoch 30 --scale 0.4 --outdir results/my_designs --datadir genie2_pdbs --num_samples 1000
```

### Running with SALAD

The converter makes all output files automatically compatible with SALAD:

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

## SALAD Compatibility Features

The converter automatically makes output files compatible with SALAD by:

1. Ensuring residue numbering is continuous (1,2,3,...) within each chain
2. Expanding motif definitions to cover full chains while preserving which parts should be fixed vs. designable
3. Adjusting array shapes to avoid SALAD's size mismatch errors

These compatibility fixes are applied automatically for all outputs. This eliminates common errors when using the generated files with SALAD's multi-motif scaffolding.

## Key Features

- **Compatible with Both Genie2 and SALAD**: Generate files that work with both protein design platforms
- **Chain-Based Group Assignment**: Automatically assigns different groups for different chains (ideal for heterodimers)
- **Strict Residue Validation**: Ensures that specified residues exist in the PDB file
- **Full-Chain Motif Support**: Properly handles the chain boundaries for SALAD compatibility
- **Multi-Chain Support**: Handles motifs spanning multiple chains
- **Flexible Linkers**: Supports variable-length linkers between motifs
