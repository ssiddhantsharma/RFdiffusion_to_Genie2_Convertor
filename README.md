# RFDiffusion to Genie2 Converter

A Python tool to convert RFDiffusion-style motif specifications to [Genie2](https://github.com/aqlaboratory/genie2) format for protein design. This tool is particularly useful for preparing input files for [SALAD's](https://github.com/mjendrusch/salad) multi-motif scaffolding.

## Installation
```bash
git clone https://github.com/ssiddhantsharma/rfdiffusion-to-genie2.git
cd rfdiffusion-to-genie2
chmod +x genie2.py
```

## Basic Usage

### Converting a Single PDB File

```bash
python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir
```

### Converting Multiple PDB Files (Same Specification)

```bash
python genie2.py --pdb_dir pdb_directory --input "A1-80[M1]/30/[M2]B81-100" --output output_dir
```

## RFDiffusion Format Syntax

The input format follows this pattern: `[Chain][Start]-[End][Motif]/[Linker]/[Chain][Start]-[End][Motif]`

Examples:
- `A1-80[M1]/30/[M2]B81-100` - Two motifs on different chains with a 30-residue linker
- `A1-29[M1]/30/[M2]A39-50/40/[M3]A2-49` - Three motifs on the same chain with linkers

Components:
- **Chain**: Single letter (e.g., A, B)
- **Residue Range**: Start-End (e.g., 1-80)
- **Motif Tag** (optional): [M1], [M2], etc.
- **Linker**: Numeric value representing amino acid length

## Advanced Options

### Terminal Scaffolds

By default, no terminal scaffolds are added. To add them:

```bash
python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --add_terminal_scaffolds
```

### Scaffold Length

Control scaffold length range (default: 5-20):

```bash
python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --min_scaffold 10 --max_scaffold 30
```

### Total Length Factors

Adjust the minimum and maximum total sequence length calculation:

```bash
python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --min_factor 1.2 --max_factor 2.0
```

### Processing Multiple Files with Different Specifications

Create a CSV file (e.g., `specs.csv`):
```
pdb_file,rf_format
protein1.pdb,A1-80[M1]/30/[M2]B81-100
protein2.pdb,A1-29[M1]/30/[M2]A39-50
```

Then run:
```bash
python genie2.py --pdb_dir pdb_directory --csv specs.csv --output output_dir
```

Alternatively, use a JSON file:
```bash
python genie2.py --pdb_dir pdb_directory --json specs.json --output output_dir
```

### Sequential Processing

Force sequential (non-parallel) processing:

```bash
python genie2.py --pdb_dir pdb_directory --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --sequential
```

### Validation Options

Use lenient validation (continue with warnings for missing residues):

```bash
python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output output_dir --lenient
```

## Output Format

For each input PDB file, the script generates a corresponding file in the output directory with "_genie2.pdb" appended to the original filename. Each output file contains:

1. Genie2 header with motif and scaffold specifications
2. Original PDB structure data

Example header:
```
REMARK 999 NAME   protein_motifs
REMARK 999 PDB    protein
REMARK 999 INPUT  A   1  80 A
REMARK 999 INPUT     25  35
REMARK 999 INPUT  B  81 100 B
REMARK 999 MINIMUM TOTAL LENGTH      125
REMARK 999 MAXIMUM TOTAL LENGTH      250
```

## Key Features and Implementation Notes

- **Residue Order Matching**: The converter ensures that motif order in the Genie2 specification matches the order defined in the RFDiffusion format specification.

- **Multi-Chain Support**: The script correctly handles motifs spanning multiple chains (e.g., `A1-80[M1]/30/[M2]B81-100`).

- **PDB Name Handling**: The converter automatically extracts and uses the correct PDB name in the Genie2 format.

- **Overlapping Motifs**: The script properly handles overlapping motif specifications (e.g., `A1-29[M1]/30/[M2]A10-40`).

- **Group Assignment**: Motifs are assigned to groups (A or B) based on their position in the specification, following Genie2 conventions.

- **Residue Validation**: The script validates that the specified motif residues actually exist in the PDB file and provides warnings for missing residues.

- **HETATM Support**: Handles non-standard residues (HETATMs) in motif definitions.

## Validation Features

The script performs strict validation by default to ensure compatibility between your motif specifications and the actual PDB content:

- **Strict Validation (Default)**: The script verifies that all specified residues exist in the PDB file and will fail if any residues are missing. This helps prevent issues when using the output with SALAD.

- **Lenient Validation (Optional)**: If you need to proceed despite missing residues, you can use the `--lenient` flag to allow the conversion to complete with warnings only.

- **Automatic Corrections**: The script automatically handles issues like swapped start/end residue numbers and trailing slashes in the format specification.

## Integration with SALAD

This converter is designed to work seamlessly with SALAD's multi-motif scaffolding capabilities. SALAD requires motif-annotated PDB files in Genie2 format, which is exactly what this converter produces.

### Using with SALAD

1. Convert your PDB file(s) to Genie2 format:
   ```bash
   python genie2.py --pdb_file your_protein.pdb --input "A1-80[M1]/30/[M2]B81-100" --output genie2_pdbs
   ```

2. Use the converted file with SALAD:
   ```bash
   python salad/training/eval_motif_benchmark.py --template genie2_pdbs/your_protein_genie2.pdb --num_designs 10 --num_steps 20 --config your_config.yml --params your_params.pt
   ```

### SALAD Multi-motif Scaffolding

SALAD provides two scripts for multi-motif scaffolding:
- `salad/training/eval_motif_benchmark.py` - For models trained with multi-motif conditioning
- `salad/training/eval_motif_benchmark_nocond.py` - Uses structure-editing to scaffold motifs without explicit motif conditioning

Both scripts accept a `--template` argument that points to your Genie2-formatted PDB file:
```bash
python salad/training/eval_motif_benchmark.py --template your_protein_genie2.pdb [other_options]
```

The converter ensures that your PDB files are correctly formatted for SALAD, allowing you to seamlessly integrate RFDiffusion-style specifications into the SALAD workflow.
