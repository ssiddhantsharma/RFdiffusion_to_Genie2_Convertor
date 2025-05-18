# RFDiffusion to Genie2 Converter

A Python tool to convert RFDiffusion-style motif specifications to [Genie2](https://github.com/aqlaboratory/genie2) format for protein design. This tool is particularly useful for preparing input files for both Genie2's motif scaffolding and [SALAD's](https://github.com/mjendrusch/salad) multi-motif scaffolding.

## Installation

```bash
git clone https://github.com/ssiddhantsharma/rfdiffusion-to-genie2.git
cd rfdiffusion-to-genie2
chmod +x rfd2genie.py
```

## Basic Usage

### Converting a Single PDB File

```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --verbose
```

## RFDiffusion Format Syntax

The input format follows this pattern: `[N-term extension]/[Chain][Start]-[End][Motif Tag]/[Linker]/[Chain][Start]-[End][Motif Tag]/[C-term extension]`

### Examples

- **Multiple chains with linker:** `A1-80[M1]/30/[M2]B81-100`
- **Multiple motifs on same chain:** `A1-29[M1]/30/[M2]A39-50/40/[M3]A2-49`
- **Heterodimer (no linker):** `A1-30[M1]/B89-95[M2]`
- **Heterotrimer:** `A1-30[M1]/B10-40[M2]/C5-25[M3]`
- **With N-terminal extension:** `20/A1-80[M1]/30/[M2]B81-100`
- **With C-terminal extension:** `A1-80[M1]/30/[M2]B81-100/40`
- **With both N and C-terminal extensions:** `15/A1-80[M1]/30/[M2]B81-100/25`
- **Variable length extensions/linkers:** `10-20/A1-80[M1]/25-35/[M2]B81-100/15-25`

### Components

- **Chain:** Single letter (e.g., A, B)
- **Residue Range:** Start-End (e.g., 1-80)
- **Motif Tag** (optional): [M1], [M2], etc.
- **Linker/Extension:** Numeric value representing amino acid length or a range (e.g., 30 or 30-40)
- **N-terminal Extension:** Number at beginning of format string followed by "/"
- **C-terminal Extension:** Number at end of format string preceded by "/"

### N-terminal Extensions

Add a number at the beginning of your format string:
```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "20/A1-80[M1]/30/[M2]B81-100" --output output_dir --verbose
```
This will generate a 20-residue N-terminal extension before the first motif.

### C-terminal Extensions

Add a number at the end of your format string:
```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100/40" --output output_dir --verbose
```
This will generate a 40-residue C-terminal extension after the last motif.

### Combined Terminal Extensions

Use both N and C terminal extensions:
```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "15/A1-80[M1]/30/[M2]B81-100/25" --output output_dir --verbose
```
This creates a 15-residue N-terminal extension and a 25-residue C-terminal extension.

### Variable Length Extensions

For more flexible designs, specify ranges:
```bash
python rfd2genie.py --pdb_file your_protein.pdb --input "10-20/A1-80[M1]/30/[M2]B81-100/15-25" --output output_dir --verbose
```
This allows the N-terminal extension to be 10-20 residues and the C-terminal extension to be 15-25 residues.

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

With N and C terminal extensions (`15/A1-80[M1]/30/[M2]B81-100/25`):

```
REMARK 999 NAME   protein_name_motifs
REMARK 999 PDB    protein_name
REMARK 999 INPUT     15  20
REMARK 999 INPUT  A   1  80 A
REMARK 999 INPUT     30  35
REMARK 999 INPUT  B  81 100 B
REMARK 999 INPUT     25  30
REMARK 999 MINIMUM TOTAL LENGTH      156
REMARK 999 MAXIMUM TOTAL LENGTH      265
```

### Important Genie2 Format Considerations

The converter addresses several aspects of the Genie2 format:

1. **Column Precision**: Generates output following Genie2's column-specific formatting requirements with proper justification (right/left). This is absolutely critical for SALAD compatibility.

2. **Group Assignment for Multi-Motif Scaffolding**: Assigns a unique group letter (A, B, C, etc.) to each motif to control their relative positioning and orientation during scaffolding. This is essential for multi-motif scaffolding where each motif needs to move independently.

3. **PDB Name Integration**: Places the PDB name only in the `REMARK 999 PDB` line, not in motif segment definitions.

4. **Residue Validation**: Verifies that specified residues exist in the PDB structure.

5. **Residue Reordering**: Automatically reorders residues in the PDB file to match the order specified in the RFDiffusion format string, ensuring compatibility with Genie2's requirements.

6. **Auto-correction**: Fixes common motif definition issues.

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

### Processing Multiple Files with Different Specifications

Create a CSV file (e.g., `specs.csv`):
```
pdb_file,rf_format
protein1.pdb,A1-80[M1]/30/[M2]B81-100
protein2.pdb,A1-29[M1]/30/[M2]A39-50
protein3.pdb,15/A1-30[M1]/25/B10-40[M2]/20
```

Then run:
```bash
python rfd2genie.py --pdb_dir pdb_directory --csv specs.csv --output output_dir --verbose
```

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
python rfd2genie.py --pdb_dir my_motifs --input "A1-80[M1]/30/[M2]B81-100" --output genie2_pdbs --verbose

# Then run Genie2 on the converted files
python genie/sample_scaffold.py --name base --epoch 30 --scale 0.4 --outdir results/my_designs --datadir genie2_pdbs --num_samples 1000
```

### Running with SALAD

The converter ensures optimal compatibility with SALAD by implementing several critical fixes:

1. **Exact Formatting**: The script strictly adheres to the column-specific formatting requirements of Genie2, with precise spacing, justification, and positioning. This is crucial because SALAD uses exact column positions for parsing, and even a small formatting deviation can cause hours-long processing or failure.

2. **Full Chain Coverage**: Ensures motif definitions cover entire chains rather than just segments, which resolves the common "boolean index did not match shape of indexed array" error in SALAD. This fix prevents SALAD from running indefinitely or for extremely long periods.

3. **Consistent Residue Numbering**: Residues are renumbered to be continuous from 1 to N for each chain, which is what SALAD expects. This eliminates array size mismatches.

4. **Removal of Extraneous REMARK Lines**: Filters out non-standard REMARK lines that could confuse SALAD's parser.

5. **Chain and Group Consistency**: Ensures chain IDs in motif definitions match those in the PDB structure, and that motif groups are assigned correctly.

#### SALAD Multi-Motif Scaffolding Example:

1. Convert your PDB file to Genie2 format:
   ```bash
   python rfd2genie.py --pdb_file your_protein.pdb --input "15/A1-80[M1]/30/[M2]B81-100/25" --output genie2_pdbs --verbose
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

- **Precise Format Compliance**: Implements exact Genie2 format requirements with proper spacing and justification
- **Chain Completeness Validation**: Ensures all chains in the PDB have proper motif definitions
- **Efficient Processing**: Streamlined code that prevents SALAD from running indefinitely
- **Robust Error Handling**: Identifies and resolves common issues that cause SALAD failures
- **Full Multi-Chain Support**: Properly handles heterodimers, heterotrimers, and other multi-chain complexes
- **Flexible Terminal Extensions**: Supports both N-terminal and C-terminal extensions without script modifications
