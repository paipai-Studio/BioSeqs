///|
/// Test data: a simple VCF file for testing.
let test_vcf_content = "##fileformat=VCFv4.2\n##fileDate=20240101\n##source=TestVCF\n##reference=GRCh38\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\nchr1\t100\trs123\tA\tT\t60.5\tPASS\tAC=1;AF=0.333;DP=12\tGT:AD\t0/1:3,5\t1/1:0,8\t0/0:6,0\nchr1\t200\t.\tC\tG,A\t45.0\t.\tAC=2,1;AF=0.5,0.25;DP=8\tGT:AD\t0/1:2,3\t1/2:1,2,1\t0/0:4,0\nchr2\t500\trs456\tAT\tA\t30.0\tLowQual\tDP=5\tGT\t0/1\t0/0\t1/1\nchr2\t600\t.\tG\t.\t.\tFAIL\t.\tGT\t./.\t0/0\t./.\n"

///|
test "vcf_parse_basic" {
  let vcf = @src.parse_vcf(test_vcf_content)

  assert_eq(vcf.num_meta_headers(), 11)
  assert_eq(vcf.len(), 4)
  assert_eq(vcf.get_version(), Some("VCFv4.2"))

  let sample_names = vcf.get_sample_names()
  assert_eq(sample_names.length(), 3)
  assert_eq(sample_names[0], "sample1")
  assert_eq(sample_names[1], "sample2")
  assert_eq(sample_names[2], "sample3")
}

///|
test "vcf_parse_record" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let records = vcf.records

  let r1 = records[0]
  assert_eq(r1.chrom, "chr1")
  assert_eq(r1.pos, 100)
  assert_eq(r1.ref, "A")
  assert_eq(r1.alt.length(), 1)
  assert_eq(r1.alt[0], "T")
  assert_eq(r1.qual, Some(60.5))

  let r3 = records[2]
  assert_eq(r3.chrom, "chr2")
  assert_eq(r3.pos, 500)
  assert_eq(r3.ref, "AT")
  assert_eq(r3.alt[0], "A")

  let r4 = records[3]
  assert_eq(r4.qual, None)
}

///|
test "vcf_parse_info" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let records = vcf.records

  let r1 = records[0]
  assert_eq(r1.info.get("AC"), Some("1"))
  assert_eq(r1.info.get("AF"), Some("0.333"))
  assert_eq(r1.info.get("DP"), Some("12"))

  let r2 = records[1]
  assert_eq(r2.info.get("AC"), Some("2,1"))
  assert_eq(r2.info.get("AF"), Some("0.5,0.25"))

  let r4 = records[3]
  assert_true(r4.info.get("AC").is_none())
}

///|
test "vcf_parse_filter" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let records = vcf.records

  let r1 = records[0]
  assert_eq(r1.filter.length(), 1)
  assert_eq(r1.filter[0], "PASS")

  let r2 = records[1]
  assert_eq(r2.filter.length(), 0)

  let r3 = records[2]
  assert_eq(r3.filter.length(), 1)
  assert_eq(r3.filter[0], "LowQual")

  let r4 = records[3]
  assert_eq(r4.filter.length(), 1)
  assert_eq(r4.filter[0], "FAIL")
}

///|
test "vcf_parse_samples" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let records = vcf.records

  let r1 = records[0]
  let gt1 = r1.get_genotype("sample1")
  assert_true(gt1.is_some())
  assert_eq(gt1.unwrap().length(), 2)

  let gt2 = r1.get_genotype("sample2")
  assert_true(gt2.is_some())
}

///|
test "vcf_variant_types" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let records = vcf.records

  let r1 = records[0]
  assert_true(r1.is_snp())
  assert_true(!r1.is_insertion())
  assert_true(!r1.is_deletion())

  let r3 = records[2]
  assert_true(!r3.is_snp())
  assert_true(!r3.is_insertion())
  assert_true(r3.is_deletion())
}

///|
test "vcf_meta_headers" {
  let vcf = @src.parse_vcf(test_vcf_content)
  let headers = vcf.meta_headers

  let mut found_contig = false
  let mut found_format = false
  let mut i = 0
  while i < headers.length() {
    if headers[i].key == "contig" {
      found_contig = true
    }
    if headers[i].key == "FORMAT" {
      found_format = true
    }
    i = i + 1
  }

  assert_true(found_contig)
  assert_true(found_format)
}
