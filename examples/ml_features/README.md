# Machine Learning Feature Extraction

Demonstrates extracting features from DNA and protein sequences
for machine learning tasks.

## Features

- **DNA Features**:
  - Nucleotide composition (A/T/G/C)
  - k-mer frequency (di/tri/tetra nucleotides)
  - Canonical k-mer frequency
  - GC content
  - Comprehensive DNA feature vector

- **Protein Features**:
  - Amino acid composition
  - Dipeptide/tripeptide composition
  - Physical-chemical properties (hydrophobicity, charge, polarity)
  - Molecular weight
  - Secondary structure propensity (Chou-Fasman)
  - Position-specific features
  - Comprehensive protein feature vector

## Run

```bash
moon run examples/ml_features/main.mbt
```
