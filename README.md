# SNP Array Primer Design in Brassica juncea

## Overview
A Python pipeline for designing gene-specific PCR primers in *Brassica juncea* (Varuna). Given a target gene ID, it extracts the gene's genomic sequence from a GFF3/FASTA pair, generates candidate primer pairs with Primer3, and validates each candidate's uniqueness against the genome using BLAST 

## How it works
1. **Extract gene sequence** — parses the GFF file for the target gene's coordinates and pulls the corresponding sequence from the genome FASTA.
2. **Design candidates** — generates up to 100 forward/reverse primer pairs with Primer3 (18–23 bp, 50–60% GC, 800–2000 bp product size).
3. **Filter for GC clamp** — keeps only candidates where both primers start with a G or C base.
4. **Validate uniqueness** — BLASTs each candidate against a local nucleotide database (`blastn-short`) and confirms all hits fall on the target chromosome.
5. **Report** — writes the first unique, GC-clamped primer pair found to a tab-separated output file.

## Requirements
- Python 3
- [`primer3-py`](https://pypi.org/project/primer3-py/)
- [`biopython`](https://pypi.org/project/biopython/)
- NCBI BLAST+ (`blastn` on your `PATH`) with a pre-built local BLAST database

Install the Python dependencies:
```bash
pip install primer3-py biopython
```

## Input
- **Genome FASTA** — `Brassica_juncea_Varuna.genome_pk.fasta`
- **GFF3 annotation** — `Brassica_juncea_Varuna_pk.gff`
- **Local BLAST database** — built from the same genome (`makeblastdb -in <fasta> -dbtype nucl -out Brassica_juncea_varuna`)

## Usage
Edit the file paths and target gene ID at the top of `Primerdesign_script.py`:
```python
gff_file = "Brassica_juncea_Varuna_pk.gff"
fasta_file = "Brassica_juncea_Varuna.genome_pk.fasta"
blast_db = "Brassica_juncea_varuna"
output_file = "unique_primers.txt"
target_gene_id = "BjuVaA01g004457"
```
Then run:
```bash
python Primerdesign_script.py
```

## Output
A tab-separated file (`unique_primers.txt`) containing the gene ID and its validated forward/reverse primer pair:
```
GeneID	Forward_Primer	Reverse_Primer
BjuVaA01g004457	...	...
```

## Notes
- The script processes one gene per run; to screen multiple genes, loop over gene IDs or extend the script accordingly.
- Only the first passing primer pair is reported, not all valid candidates.

## Author
**Srishti**
