#!/usr/bin/env python3
"""Compare pyfaidx results with MoonBit FastaIndex implementation."""
import os
import tempfile

FASTA_CONTENT = """>chr1
ACGTACGTACGT
>chr2
TTTTAAAA
>chr3
GATTACAGATTACAGATTACA
"""

MULTILINE_FASTA = """>chr1
ACGT
ACGT
ACGT
>chr2
TT
TTAA
AA
"""


def test_pyfaidx():
    """Test with pyfaidx if available."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        print("pyfaidx not installed, skipping comparison")
        return False

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(FASTA_CONTENT)
        fasta_path = f.name

    try:
        fasta = Fasta(fasta_path)

        results = []

        # Test full sequence retrieval
        results.append(("chr1_full", str(fasta['chr1'])))
        results.append(("chr2_full", str(fasta['chr2'])))
        results.append(("chr3_full", str(fasta['chr3'])))

        # Test subsequence retrieval
        results.append(("chr1_0_4", str(fasta['chr1'][0:4])))
        results.append(("chr1_4_8", str(fasta['chr1'][4:8])))
        results.append(("chr1_0_12", str(fasta['chr1'][0:12])))

        # Test lengths
        results.append(("chr1_len", len(fasta['chr1'])))
        results.append(("chr2_len", len(fasta['chr2'])))
        results.append(("chr3_len", len(fasta['chr3'])))

        # Test names
        results.append(("names", list(fasta.keys())))

        for name, value in results:
            print(f"  {name}: {value}")

        print("\n=== FAI file content ===")
        with open(fasta_path + '.fai') as f:
            print(f.read())

        return True

    finally:
        os.unlink(fasta_path)
        if os.path.exists(fasta_path + '.fai'):
            os.unlink(fasta_path + '.fai')


def test_multiline():
    """Test multiline FASTA with pyfaidx."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        return False

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(MULTILINE_FASTA)
        fasta_path = f.name

    try:
        fasta = Fasta(fasta_path)

        print("\n=== Multiline FASTA tests ===")
        print(f"  chr1_full: {str(fasta['chr1'])}")
        print(f"  chr1_len: {len(fasta['chr1'])}")
        print(f"  chr2_full: {str(fasta['chr2'])}")
        print(f"  chr2_len: {len(fasta['chr2'])}")

        # Subsequence across lines
        print(f"  chr1_3_9: {str(fasta['chr1'][3:9])}")

        print("\n=== Multiline FAI file content ===")
        with open(fasta_path + '.fai') as f:
            print(f.read())

        return True

    finally:
        os.unlink(fasta_path)
        if os.path.exists(fasta_path + '.fai'):
            os.unlink(fasta_path + '.fai')


if __name__ == '__main__':
    print("=== pyfaidx Reference Tests ===")
    print("\n--- Basic FASTA tests ---")
    test_pyfaidx()
    print("\n--- Multiline FASTA tests ---")
    test_multiline()
