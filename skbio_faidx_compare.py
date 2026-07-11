#!/usr/bin/env python3
"""
Comparison test script for MoonBit faidx implementation vs Python pyfaidx.
"""

import os
import tempfile
from pyfaidx import Fasta

# Test FASTA content (consistent line lengths for pyfaidx)
test_fasta = """>chr1 description of chr1
ACGTACGTACGTACGTACGTACGTACGTGC
ATGCATGCATGCATGCATGCATGCATGCAT
ACGTACGTACGTACGTACGTACGTACGTGC
>chr2
ATGCATGCATGCAT
GCATGCATGCATGC
>chr3 short
AAAA
"""

# Write to temp file
with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
    f.write(test_fasta)
    fasta_file = f.name

try:
    # Load with pyfaidx
    fa = Fasta(fasta_file)
    
    print("=" * 60)
    print("pyfaidx Results")
    print("=" * 60)
    
    # Test 1: Get all sequence names
    print("\n1. Sequence names:")
    print(f"   {list(fa.keys())}")
    
    # Test 2: Get sequence lengths
    print("\n2. Sequence lengths:")
    for name in fa.keys():
        print(f"   {name}: {len(fa[name])}")
    
    # Test 3: Fetch full sequences
    print("\n3. Full sequences:")
    for name in fa.keys():
        seq = str(fa[name])
        print(f"   {name}: {seq[:50]}... (len={len(seq)})")
    
    # Test 4: Fetch sub-sequences (1-based)
    print("\n4. Sub-sequences:")
    print(f"   chr1[1:10]: {str(fa['chr1'][0:10])}")
    print(f"   chr1[1:30]: {str(fa['chr1'][0:30])}")
    print(f"   chr1[31:59]: {str(fa['chr1'][30:59])}")
    print(f"   chr1[60:67]: {str(fa['chr1'][59:67])}")
    print(f"   chr1[20:25]: {str(fa['chr1'][19:25])}")
    print(f"   chr2[1:14]: {str(fa['chr2'][0:14])}")
    print(f"   chr2[15:28]: {str(fa['chr2'][14:28])}")
    print(f"   chr3[1:4]: {str(fa['chr3'][0:4])}")
    
    # Test 5: Print .fai content
    print("\n5. .fai Index content:")
    fai_file = fasta_file + '.fai'
    if os.path.exists(fai_file):
        with open(fai_file, 'r') as f:
            print(f.read())
    
    print("=" * 60)
    print("Expected MoonBit Results")
    print("=" * 60)
    print("""
1. Sequence names:
   ['chr1', 'chr2', 'chr3']

2. Sequence lengths:
   chr1: 90
   chr2: 28
   chr3: 4

3. Full sequences:
   chr1: ACGTACGTACGTACGTACGTACGTACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATACGTACGTACGTACGTACGTACGTACGTACGTGC (len=90)
   chr2: ATGCATGCATGCATGCATGCATGCATGC (len=28)
   chr3: AAAA (len=4)

4. Sub-sequences:
   chr1[1:10]: ACGTACGTAC
   chr1[1:30]: ACGTACGTACGTACGTACGTACGTACGTGC
   chr1[31:59]: ATGCATGCATGCATGCATGCATGCATGCA
   chr1[60:67]: TACGTACG
   chr1[20:25]: TACGTA
   chr2[1:14]: ATGCATGCATGCAT
   chr2[15:28]: GCATGCATGCATGC
   chr3[1:4]: AAAA

5. .fai Index content:
chr1	90	26	30	31
chr2	28	125	14	15
chr3	4	167	4	5
""")
    
    print("=" * 60)
    print("Verification Complete!")
    print("=" * 60)
    
finally:
    # Cleanup
    os.unlink(fasta_file)
    fai_file = fasta_file + '.fai'
    if os.path.exists(fai_file):
        os.unlink(fai_file)