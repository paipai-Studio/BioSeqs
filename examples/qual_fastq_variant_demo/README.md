# QUAL / FASTQ Variant Demo

Demonstrates QUAL format, FASTQ variant encodings, quality score conversion, and PairedFastaQualIterator:

- QUAL format parsing and writing
- FASTQ Illumina (ASCII offset 64) parsing and writing
- FASTQ Solexa (Solexa scores with ASCII offset 64) parsing and writing
- Quality score conversion (Phred ↔ Solexa)
- PairedFastaQualIterator (FASTA + QUAL matching)

## Run

```bash
moon run examples/qual_fastq_variant_demo/main.mbt
```