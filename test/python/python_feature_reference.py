#!/usr/bin/env python3
"""
Python reference implementation for machine learning feature extraction.
Used for comparison with MoonBit implementation.

Run: python test/python/python_feature_reference.py
"""

import sys
from collections import defaultdict

def kmer_frequency(seq, k, normalize=True):
    """Compute k-mer frequency for a sequence"""
    n = len(seq)
    if k <= 0 or k > n:
        return {}
    
    freq = defaultdict(float)
    total = n - k + 1
    
    for i in range(total):
        kmer = seq[i:i+k]
        freq[kmer] += 1.0
    
    if normalize and total > 0:
        for kmer in freq:
            freq[kmer] /= total
    
    return dict(freq)

def reverse_complement_str(s):
    """Reverse complement of DNA string"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement.get(c, c) for c in reversed(s))

def dna_canonical_kmer_frequency(seq, k, normalize=True):
    """Compute canonical k-mer frequency"""
    n = len(seq)
    if k <= 0 or k > n:
        return {}
    
    freq = defaultdict(float)
    total = n - k + 1
    
    for i in range(total):
        kmer = seq[i:i+k]
        rc = reverse_complement_str(kmer)
        canonical = min(kmer, rc)
        freq[canonical] += 1.0
    
    if normalize and total > 0:
        for kmer in freq:
            freq[kmer] /= total
    
    return dict(freq)

def nucleotide_composition(seq):
    """Compute nucleotide composition"""
    n = len(seq)
    if n == 0:
        return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0, 'other': 0.0}
    
    counts = {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0, 'other': 0.0}
    
    for c in seq:
        c_upper = c.upper()
        if c_upper in counts:
            counts[c_upper] += 1.0
        else:
            counts['other'] += 1.0
    
    for key in counts:
        counts[key] /= n
    
    return counts

def amino_acid_composition(seq):
    """Compute amino acid composition"""
    n = len(seq)
    if n == 0:
        return {}
    
    freq = defaultdict(float)
    for c in seq:
        freq[c.upper()] += 1.0
    
    for aa in freq:
        freq[aa] /= n
    
    return dict(freq)

def dipeptide_composition(seq):
    """Compute dipeptide composition"""
    n = len(seq)
    if n < 2:
        return {}
    
    freq = defaultdict(float)
    total = n - 1
    
    for i in range(total):
        dipeptide = seq[i:i+2].upper()
        freq[dipeptide] += 1.0
    
    for dp in freq:
        freq[dp] /= total
    
    return dict(freq)

def tripeptide_composition(seq):
    """Compute tripeptide composition"""
    n = len(seq)
    if n < 3:
        return {}
    
    freq = defaultdict(float)
    total = n - 2
    
    for i in range(total):
        tripeptide = seq[i:i+3].upper()
        freq[tripeptide] += 1.0
    
    for tp in freq:
        freq[tp] /= total
    
    return dict(freq)

# Kyte-Doolittle hydrophobicity scale
AA_HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

AA_CHARGE = {
    'A': 0.0, 'R': 1.0, 'N': 0.0, 'D': -1.0, 'C': 0.0,
    'E': -1.0, 'Q': 0.0, 'G': 0.0, 'H': 0.0, 'I': 0.0,
    'L': 0.0, 'K': 1.0, 'M': 0.0, 'F': 0.0, 'P': 0.0,
    'S': 0.0, 'T': 0.0, 'W': 0.0, 'Y': 0.0, 'V': 0.0
}

AA_POLARITY = {
    'A': 0.0, 'R': 1.0, 'N': 1.0, 'D': 1.0, 'C': 0.5,
    'E': 1.0, 'Q': 1.0, 'G': 0.0, 'H': 1.0, 'I': 0.0,
    'L': 0.0, 'K': 1.0, 'M': 0.0, 'F': 0.0, 'P': 0.0,
    'S': 1.0, 'T': 1.0, 'W': 0.5, 'Y': 0.5, 'V': 0.0
}

def avg_hydrophobicity(seq):
    """Compute average hydrophobicity"""
    n = len(seq)
    if n == 0:
        return 0.0
    total = sum(AA_HYDROPHOBICITY.get(c.upper(), 0.0) for c in seq)
    return total / n

def avg_charge(seq):
    """Compute average charge"""
    n = len(seq)
    if n == 0:
        return 0.0
    total = sum(AA_CHARGE.get(c.upper(), 0.0) for c in seq)
    return total / n

def avg_polarity(seq):
    """Compute average polarity"""
    n = len(seq)
    if n == 0:
        return 0.0
    total = sum(AA_POLARITY.get(c.upper(), 0.0) for c in seq)
    return total / n

def molecular_weight(seq):
    """Compute molecular weight"""
    n = len(seq)
    if n == 0:
        return 0.0
    
    weight = 18.0  # H2O
    
    AA_WEIGHTS = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
        'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
    }
    
    for c in seq:
        weight += AA_WEIGHTS.get(c.upper(), 0.0) - 18.0
    
    return weight

def secondary_structure_propensity(seq):
    """Compute secondary structure propensity"""
    n = len(seq)
    if n == 0:
        return (0.0, 0.0, 0.0)
    
    helix_p = {
        'A': 1.42, 'R': 1.05, 'N': 0.67, 'D': 1.01, 'C': 0.70,
        'E': 1.51, 'Q': 1.11, 'G': 0.57, 'H': 1.00, 'I': 1.08,
        'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
        'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 1.07, 'V': 1.06
    }
    
    sheet_p = {
        'A': 0.83, 'R': 0.74, 'N': 0.89, 'D': 0.54, 'C': 1.19,
        'E': 0.37, 'Q': 1.11, 'G': 0.75, 'H': 0.87, 'I': 1.60,
        'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
        'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
    }
    
    helix = sum(helix_p.get(c.upper(), 0.0) for c in seq)
    sheet = sum(sheet_p.get(c.upper(), 0.0) for c in seq)
    
    return (helix / n, sheet / n, (n - helix - sheet) / n)

def protein_feature_vector(seq):
    """Compute comprehensive protein feature vector"""
    features = []
    
    aa_order = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    aa_comp = amino_acid_composition(seq)
    for aa in aa_order:
        features.append(aa_comp.get(aa, 0.0))
    
    dipep = dipeptide_composition(seq)
    for k in list(dipep.keys())[:20]:
        features.append(dipep[k])
    
    features.append(avg_hydrophobicity(seq))
    features.append(avg_charge(seq))
    features.append(avg_polarity(seq))
    features.append(molecular_weight(seq))
    
    h, s, c = secondary_structure_propensity(seq)
    features.append(h)
    features.append(s)
    features.append(c)
    
    return features

def dna_feature_vector(seq):
    """Compute comprehensive DNA feature vector"""
    features = []
    
    comp = nucleotide_composition(seq)
    features.append(comp['A'])
    features.append(comp['T'])
    features.append(comp['G'])
    features.append(comp['C'])
    
    dinuc = kmer_frequency(seq, 2)
    dinuc_order = ["AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC",
                   "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC"]
    for d in dinuc_order:
        features.append(dinuc.get(d, 0.0))
    
    trinuc = kmer_frequency(seq, 3)
    for k in list(trinuc.keys())[:20]:
        features.append(trinuc[k])
    
    features.append(comp['G'] + comp['C'])
    features.append(len(seq))
    
    return features

def main():
    """Run all tests"""
    test_cases = []
    
    # k-mer frequency
    result = kmer_frequency("ATGCATGC", 2)
    test_cases.append(("kmer_frequency", str(result)))
    
    # nucleotide composition
    result = nucleotide_composition("ATGCATGC")
    test_cases.append(("nucleotide_composition", str(result)))
    
    # amino acid composition
    result = amino_acid_composition("MRKKA")
    test_cases.append(("amino_acid_composition", str(result)))
    
    # dipeptide composition
    result = dipeptide_composition("MRKK")
    test_cases.append(("dipeptide_composition", str(result)))
    
    # avg hydrophobicity
    result = avg_hydrophobicity("MLKK")
    test_cases.append(("avg_hydrophobicity", str(result)))
    
    # avg charge
    result = avg_charge("RK")
    test_cases.append(("avg_charge", str(result)))
    
    # molecular weight
    result = molecular_weight("A")
    test_cases.append(("molecular_weight", str(result)))
    
    # secondary structure
    result = secondary_structure_propensity("MLKK")
    test_cases.append(("secondary_structure", str(result)))
    
    # protein feature vector
    result = protein_feature_vector("MLKK")
    test_cases.append(("protein_feature_vector", str(len(result))))
    
    # DNA feature vector
    result = dna_feature_vector("ATGCATGC")
    test_cases.append(("dna_feature_vector", str(len(result))))
    
    # Print results
    for name, value in test_cases:
        print(f'"{name}" "{value}"')

if __name__ == "__main__":
    main()
