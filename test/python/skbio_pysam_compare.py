#!/usr/bin/env python3
"""
Comparison test script for MoonBit SAM/VCF implementation vs Python pysam.
"""

import sys
import os

# Test data for SAM
test_sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:2000
@RG\tID:test\tSM:sample1
@CO\tThis is a comment
read1\t0\tchr1\t100\t60\t80M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~!\tNM:i:0\tMD:Z:80\tAS:i:80
read2\t16\tchr1\t200\t40\t50M20D30M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcd\tNM:i:2\tMD:Z:48A20D28C2
read3\t73\tchr2\t50\t50\t10S60M10S\t=\t100\t150\tCGTACGTACGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGATGCATGCATG\t\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~!\"#$%&'()*+,-./012\tNM:i:1
read4\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;
"""

# Test data for VCF
test_vcf_content = """##fileformat=VCFv4.2
##fileDate=20240101
##source=TestVCF
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
chr1\t100\trs123\tA\tT\t60.5\tPASS\tAC=1;AF=0.333;DP=12\tGT:AD\t0/1:3,5\t1/1:0,8\t0/0:6,0
chr1\t200\t.\tC\tG,A\t45.0\t.\tAC=2,1;AF=0.5,0.25;DP=8\tGT:AD\t0/1:2,3\t1/2:1,2,1\t0/0:4,0
chr2\t500\trs456\tAT\tA\t30.0\tLowQual\tDP=5\tGT\t0/1\t0/0\t1/1
chr2\t600\t.\tG\t.\t.\tFAIL\t.\tGT\t./.\t0/0\t./.
"""

def run_sam_tests():
    """Run SAM comparison tests using pysam."""
    try:
        import pysam
    except ImportError:
        print("ERROR: pysam not installed. Please install with: pip install pysam")
        return False
    
    print("=" * 60)
    print("SAM Comparison Tests")
    print("=" * 60)
    
    # Write test SAM file
    with open("/tmp/test.sam", "w") as f:
        f.write(test_sam_content)
    
    # Open with pysam
    samfile = pysam.AlignmentFile("/tmp/test.sam", "r")
    
    # Test header
    print("\n1. Header Tests:")
    n_headers = sum(1 for _ in samfile.header)
    print(f"   Number of header lines: {n_headers}")
    
    # Get sequence names
    seq_names = list(samfile.references)
    print(f"   Sequence names: {seq_names}")
    
    # Get sequence lengths
    seq_lengths = {ref: samfile.lengths[i] for i, ref in enumerate(samfile.references)}
    print(f"   Sequence lengths: {seq_lengths}")
    
    # Test records
    print("\n2. Record Tests:")
    records = list(samfile)
    print(f"   Number of records: {len(records)}")
    
    for i, rec in enumerate(records):
        print(f"\n   Record {i+1} ({rec.query_name}):")
        print(f"     flag: {rec.flag}")
        print(f"     rname: {rec.reference_name}")
        print(f"     pos: {rec.pos + 1} (1-based)")
        print(f"     mapq: {rec.mapping_quality}")
        print(f"     cigar: {rec.cigarstring}")
        print(f"     is_paired: {rec.is_paired}")
        print(f"     is_proper_pair: {rec.is_proper_pair}")
        print(f"     is_unmapped: {rec.is_unmapped}")
        print(f"     is_reverse: {rec.is_reverse}")
        print(f"     is_read1: {rec.is_read1}")
        print(f"     is_read2: {rec.is_read2}")
        print(f"     alignment_length: {rec.query_alignment_length}")
        
        # Get tags
        nm_tag = rec.get_tag("NM") if rec.has_tag("NM") else None
        md_tag = rec.get_tag("MD") if rec.has_tag("MD") else None
        print(f"     NM tag: {nm_tag}")
        print(f"     MD tag: {md_tag}")
    
    samfile.close()
    print("\nSAM tests completed successfully!")
    return True

def run_vcf_tests():
    """Run VCF comparison tests using pysam."""
    try:
        import pysam
    except ImportError:
        print("ERROR: pysam not installed. Please install with: pip install pysam")
        return False
    
    print("\n" + "=" * 60)
    print("VCF Comparison Tests")
    print("=" * 60)
    
    # Write test VCF file
    with open("/tmp/test.vcf", "w") as f:
        f.write(test_vcf_content)
    
    # Open with pysam
    vcf_file = pysam.VariantFile("/tmp/test.vcf", "r")
    
    # Test header
    print("\n1. Header Tests:")
    print(f"   Version: {vcf_file.header.version}")
    
    # Get sample names
    samples = list(vcf_file.header.samples)
    print(f"   Samples: {samples}")
    
    # Test records
    print("\n2. Record Tests:")
    records = list(vcf_file)
    print(f"   Number of records: {len(records)}")
    
    for i, rec in enumerate(records):
        print(f"\n   Record {i+1}:")
        print(f"     chrom: {rec.chrom}")
        print(f"     pos: {rec.pos}")
        print(f"     id: {rec.id}")
        print(f"     ref: {rec.ref}")
        print(f"     alt: {[str(a) for a in rec.alts]}")
        print(f"     qual: {rec.qual}")
        print(f"     filter: {rec.filter.keys()}")
        
        # INFO fields
        print(f"     INFO fields:")
        for key in rec.info.keys():
            print(f"       {key}: {rec.info[key]}")
        
        # FORMAT fields
        if rec.format:
            print(f"     FORMAT: {list(rec.format.keys())}")
        
        # Samples
        for sample in samples:
            gt = rec.samples[sample].get("GT")
            ad = rec.samples[sample].get("AD")
            print(f"     {sample}: GT={gt}, AD={ad}")
    
    vcf_file.close()
    print("\nVCF tests completed successfully!")
    return True

def main():
    """Main entry point."""
    print("pysam Comparison Test Script")
    print("=" * 60)
    
    sam_ok = run_sam_tests()
    vcf_ok = run_vcf_tests()
    
    print("\n" + "=" * 60)
    if sam_ok and vcf_ok:
        print("All comparison tests completed successfully!")
        print("\nPlease compare the output above with MoonBit test results.")
        return 0
    else:
        print("Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
