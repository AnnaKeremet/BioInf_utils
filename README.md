# `BioInf Utils`🧬

A collection of utilities for processing bioinformatics files: .FASTQ, .FASTA, and BLAST output.

## Functions:

- **`filter_fastq()`** — Filters .FASTQ reads by GC-content, sequence length and mean Phred+33 quality. Saves results to the `filtered/` directory;
- **`convert_multiline_fasta_to_oneline()`** — Converts multiline .FASTA sequences into single-line format;
- **`parse_blast_output()`** — Extracts the best protein hit description for each BLAST query and saves the results alphabetically sorted.

## Usage:

Import supplementary functions from `main.py` script:

```python
from main import filter_fastq, convert_multiline_fasta_to_oneline, parse_blast_output

# Filter FASTQ
filter_fastq("input.fastq", "output.fastq",gc_bounds=(20, 80), length_bounds=(50, 200),quality_threshold=30)

# Convert FASTA
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")

# Parse BLAST results
parse_blast_output("blast_results.txt", "proteins.txt")
```

## Project Structure:

```text
main.py
helpers/
    fastq.py
    bio_files.py
README.md
```