#!/usr/bin/env python3
"""
Python reference script for scikit-bio alignment tests.
Generates expected outputs for MoonBit implementation.
"""

import sys
sys.path.insert(0, '/tmp/skbio_env/lib/python3.12/site-packages')

from skbio import DNA, RNA, Protein
from skbio.alignment import (
    TabularMSA,
    global_pairwise_align_nucleotide,
    local_pairwise_align_nucleotide,
    global_pairwise_align_protein,
    local_pairwise_align_protein,
)

def test_dna_creation():
    dna = DNA("ACGT")
    print(f"DNA creation: '{dna}', length={len(dna)}")

def test_rna_creation():
    rna = RNA("ACGU")
    print(f"RNA creation: '{rna}', length={len(rna)}")

def test_protein_creation():
    prot = Protein("ACDEFGHIKLMNPQRSTVWY")
    print(f"Protein creation: '{prot}', length={len(prot)}")

def test_tabular_msa_basic():
    dna1 = DNA("ACGT")
    dna2 = DNA("ACCT")
    msa = TabularMSA([dna1, dna2])
    print(f"TabularMSA shape: {msa.shape}")
    print(f"Number of records: {len(msa)}")
    print(f"Alignment length: {msa.shape.position}")

def test_tabular_msa_column():
    dna1 = DNA("ACGT")
    dna2 = DNA("ACCT")
    msa = TabularMSA([dna1, dna2])
    print(f"Column 0: '{msa.iloc[:, 0]}'")
    print(f"Column 1: '{msa.iloc[:, 1]}'")
    print(f"Column 2: '{msa.iloc[:, 2]}'")
    print(f"Column 3: '{msa.iloc[:, 3]}'")

def test_consensus():
    dna1 = DNA("ACGT")
    dna2 = DNA("ACCT")
    dna3 = DNA("AGGT")
    msa = TabularMSA([dna1, dna2, dna3])
    print(f"Consensus: '{msa.consensus()}'")

def test_global_pairwise_align_nucleotide():
    seq1 = DNA("ACGT")
    seq2 = DNA("AGT")
    
    msa, score, pos = global_pairwise_align_nucleotide(seq1, seq2)
    print(f"Global alignment score: {score}")
    print(f"Positions: {pos}")
    print(f"Aligned seq1: '{msa[0]}'")
    print(f"Aligned seq2: '{msa[1]}'")

def test_local_pairwise_align_nucleotide():
    seq1 = DNA("ACGTACGT")
    seq2 = DNA("CGT")
    
    msa, score, pos = local_pairwise_align_nucleotide(seq1, seq2)
    print(f"Local alignment score: {score}")
    print(f"Positions: {pos}")
    print(f"Aligned seq1: '{msa[0]}'")
    print(f"Aligned seq2: '{msa[1]}'")

def test_global_pairwise_align_protein():
    seq1 = Protein("ACDE")
    seq2 = Protein("ACE")
    
    msa, score, pos = global_pairwise_align_protein(seq1, seq2)
    print(f"Global protein alignment score: {score}")
    print(f"Positions: {pos}")
    print(f"Aligned seq1: '{msa[0]}'")
    print(f"Aligned seq2: '{msa[1]}'")

def test_local_pairwise_align_protein():
    seq1 = Protein("ACDEFGHI")
    seq2 = Protein("EFG")
    
    msa, score, pos = local_pairwise_align_protein(seq1, seq2)
    print(f"Local protein alignment score: {score}")
    print(f"Positions: {pos}")
    print(f"Aligned seq1: '{msa[0]}'")
    print(f"Aligned seq2: '{msa[1]}'")

def test_conservation():
    dna1 = DNA("ACGT")
    dna2 = DNA("ACCT")
    dna3 = DNA("AGGT")
    msa = TabularMSA([dna1, dna2, dna3])
    conservation = msa.conservation()
    print(f"Conservation scores: {conservation}")

def test_alignment_with_gaps():
    seq1 = DNA("ATCGATCG")
    seq2 = DNA("ATATCG")
    
    msa, score, pos = global_pairwise_align_nucleotide(seq1, seq2)
    print(f"Global alignment with gaps score: {score}")
    print(f"Aligned seq1: '{msa[0]}'")
    print(f"Aligned seq2: '{msa[1]}'")

if __name__ == "__main__":
    print("=" * 60)
    print("scikit-bio alignment reference tests")
    print("=" * 60)
    
    print("\n--- DNA/RNA/Protein creation ---")
    test_dna_creation()
    test_rna_creation()
    test_protein_creation()
    
    print("\n--- TabularMSA basic operations ---")
    test_tabular_msa_basic()
    
    print("\n--- TabularMSA column access ---")
    test_tabular_msa_column()
    
    print("\n--- Consensus ---")
    test_consensus()
    
    print("\n--- Conservation ---")
    test_conservation()
    
    print("\n--- Global pairwise alignment (nucleotide) ---")
    test_global_pairwise_align_nucleotide()
    
    print("\n--- Local pairwise alignment (nucleotide) ---")
    test_local_pairwise_align_nucleotide()
    
    print("\n--- Global pairwise alignment (protein) ---")
    test_global_pairwise_align_protein()
    
    print("\n--- Local pairwise alignment (protein) ---")
    test_local_pairwise_align_protein()
    
    print("\n--- Alignment with gaps ---")
    test_alignment_with_gaps()
    
    print("\n" + "=" * 60)
    print("All tests completed!")
    print("=" * 60)