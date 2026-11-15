///|
/// BAM/VCF round-trip 序列化测试

// ---------- VCF round-trip ----------

let roundtrip_vcf_text : String =
  "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##contig=<ID=chr2,length=800>\n" +
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
  "##INFO=<ID=CSQ,Number=A,Type=String,Description=\"consequence\">\n" +
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\n" +
  "chr1\t100\trs1\tA\tT\t60\tPASS\tDP=10;CSQ=missense\tGT\t0/1\t0/0\n" +
  "chr1\t200\trs2\tC\tG,A\t50\tPASS\tDP=20\tGT\t1/1\t0/1\n" +
  "chr2\t50\t.\tG\tA\t.\tLowQual\t.\tGT\t0/0\t./.\n"

///|
test "vcf_roundtrip_idempotent_text" {
  let text1 = roundtrip_vcf_text
  let vcf = @src.parse_vcf(text1)
  let text2 = @src.vcf_to_text(vcf)
  let vcf2 = @src.parse_vcf(text2)
  // records 数应该一致
  assert_eq(vcf2.records.length(), 3)
  // 第一条基本字段
  let r1 = vcf2.records[0]
  assert_eq(r1.chrom, "chr1")
  assert_eq(r1.pos, 100)
  assert_eq(r1.ref_seq, "A")
  assert_eq(r1.alt.length(), 1)
  assert_eq(r1.alt[0], "T")
  // INFO 保留
  match r1.info_get("DP") {
    Some(v) => assert_eq(v, "10")
    None => abort("DP missing after round-trip")
  }
  match r1.info_get("CSQ") {
    Some(v) => assert_eq(v, "missense")
    None => abort("CSQ missing")
  }
  // FORMAT + samples
  match r1.format {
    Some(f) => { assert_eq(f[0], "GT") }
    None => abort("FORMAT missing")
  }
  match r1.samples.get("sampleA") {
    Some(arr) => assert_eq(arr[0], "0/1")
    None => abort("sampleA missing")
  }
  match r1.samples.get("sampleB") {
    Some(arr) => assert_eq(arr[0], "0/0")
    None => abort("sampleB missing")
  }
}

///|
test "vcf_roundtrip_multi_allelic" {
  let vcf = @src.parse_vcf(roundtrip_vcf_text)
  let text = @src.vcf_to_text(vcf)
  let vcf2 = @src.parse_vcf(text)
  let r2 = vcf2.records[1]  // rs2 chr1:200 REF=C ALT=G,A
  assert_eq(r2.alt.length(), 2)
  assert_eq(r2.alt[0], "G")
  assert_eq(r2.alt[1], "A")
}

///|
test "vcf_roundtrip_missing_qual_filter_id" {
  let vcf = @src.parse_vcf(roundtrip_vcf_text)
  let text = @src.vcf_to_text(vcf)
  let vcf2 = @src.parse_vcf(text)
  let r3 = vcf2.records[2]
  // QUAL="." → None
  match r3.qual {
    Some(_) => abort("expected None QUAL")
    None => ()
  }
  // FILTER=LowQual（不是 PASS）
  assert_eq(r3.filter.length(), 1)
  assert_eq(r3.filter[0], "LowQual")
  // ID 为空 → []
  assert_eq(r3.id.length(), 0)
}

///|
test "vcf_roundtrip_empty_vcf" {
  let vcf = @src.parse_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
  let text = @src.vcf_to_text(vcf)
  let vcf2 = @src.parse_vcf(text)
  assert_eq(vcf2.records.length(), 0)
}

///|
test "vcf_roundtrip_no_meta_headers_fallback" {
  // 用 parse_vcf 构造一份只带 fileformat 没有 contig meta header 的 VCF
  let vcf = @src.parse_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
  let text = @src.vcf_to_text(vcf)
  assert_true(text.contains("CHROM"))
  assert_true(text.contains("POS"))
}

// ---------- BAM → SAM round-trip ----------

let roundtrip_sam_text : String =
  "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:1000\n" +
  "@SQ\tSN:chr2\tLN:800\n" +
  "@RG\tID:lib1\tSM:sample1\n" +
  "read1\t0\tchr1\t101\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\tRG:Z:lib1\tNM:i:0\n" +
  "read2\t16\tchr1\t201\t50\t8M\t*\t0\t0\tTGCATGCA\t!!!!!!!!\n" +
  "pair_fwd\t67\tchr1\t301\t40\t8M\t=\t401\t150\tGGCCGGCC\t!!!!!!!!\tRG:Z:lib1\n" +
  "pair_rev\t131\tchr1\t401\t40\t8M\t=\t301\t-150\tCCGGCCGG\t!!!!!!!!\tRG:Z:lib1\n"

///|
test "bam_to_sam_basic_records" {
  let bam = @src.sam_to_bam(roundtrip_sam_text)
  let sam = @src.bam_to_sam_text(bam)
  assert_true(sam.contains("@HD"))
  assert_true(sam.contains("@SQ"))
  assert_true(sam.contains("@RG"))
  assert_true(sam.contains("read1"))
  assert_true(sam.contains("read2"))
  assert_true(sam.contains("pair_fwd"))
  assert_true(sam.contains("pair_rev"))
}

///|
test "bam_roundtrip_pos_conversion" {
  let bam = @src.sam_to_bam(roundtrip_sam_text)
  let sam = @src.bam_to_sam_text(bam)
  // SAM pos 是 1-based；BamRecord.pos 是 0-based；bam_to_sam 要 +1 回去
  assert_true(sam.contains("\t101\t"))  // chr1:100 (0-based) → SAM 101
  assert_true(sam.contains("\t201\t"))
  assert_true(sam.contains("\t301\t"))
  assert_true(sam.contains("\t401\t"))
}

///|
test "bam_to_sam_flags_roundtrip" {
  let bam = @src.sam_to_bam(roundtrip_sam_text)
  let sam = @src.bam_to_sam_text(bam)
  // flag 0 → "0"
  assert_true(sam.contains("\t0\tchr1\t"))
  // flag 16 → "16"
  assert_true(sam.contains("\t16\tchr1\t"))
  // flag 67 → "67"
  assert_true(sam.contains("\t67\tchr1\t"))
  // flag 131 → "131"
  assert_true(sam.contains("\t131\tchr1\t"))
}

///|
test "bam_to_sam_ref_names" {
  let bam = @src.sam_to_bam(roundtrip_sam_text)
  let sam = @src.bam_to_sam_text(bam)
  // chr1 作为 ref 出现在 record 行 + chr2 作为 @SQ header
  assert_true(sam.contains("\tchr1\t"))
  assert_true(sam.contains("@SQ\tSN:chr2"))
}

///|
test "bam_to_sam_tags_output" {
  let bam = @src.sam_to_bam(roundtrip_sam_text)
  let sam = @src.bam_to_sam_text(bam)
  // RG tag 应该输出 RG:Z:lib1
  assert_true(sam.contains("RG:Z:lib1"))
  // NM tag 应该输出 NM:i:0
  assert_true(sam.contains("NM:i:0"))
}

///|
test "bam_roundtrip_sam_to_bam_to_sam" {
  let bam1 = @src.sam_to_bam(roundtrip_sam_text)
  let sam1 = @src.bam_to_sam_text(bam1)
  let bam2 = @src.sam_to_bam(sam1)
  let sam2 = @src.bam_to_sam_text(bam2)
  // 记录数一致
  assert_eq(bam2.records.length(), 4)
  // qname 都在
  let mut all_qnames_ok = true
  let mut ri = 0
  while ri < bam2.records.length() {
    let qn = bam2.records[ri].qname
    if qn != "read1" && qn != "read2" && qn != "pair_fwd" && qn != "pair_rev" {
      all_qnames_ok = false
    }
    ri = ri + 1
  }
  assert_true(all_qnames_ok)
  // flag round-trip
  assert_eq(bam2.records[0].flag, 0)
  assert_eq(bam2.records[1].flag, 16)
  assert_eq(bam2.records[2].flag, 67)
  assert_eq(bam2.records[3].flag, 131)
  // seq 保留
  assert_eq(bam2.records[0].seq, "ACGTACGT")
}

///|
test "bam_roundtrip_empty_header_uses_defaults" {
  // 构造没有 header.text 的 BamFile — 序列化时应 fallback
  let bam = @src.sam_to_bam(
    "@HD\tVN:1.6\tSO:unsorted\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "r1\t0\tchr1\t101\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  )
  let sam = @src.bam_to_sam_text(bam)
  assert_true(sam.contains("@HD"))
  assert_true(sam.contains("@SQ\tSN:chr1"))
}

///|
test "bam_roundtrip_unmapped_read_star_fields" {
  let sam_text : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "unmap\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_text)
  let sam_out = @src.bam_to_sam_text(bam)
  assert_true(sam_out.contains("unmap"))
  // unmapped 的 RNAME 应该输出 "*"
  assert_true(sam_out.contains("\t*\t"))
}

///|
test "bam_roundtrip_rnext_eq_symbol" {
  let sam_text : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "pair\t67\tchr1\t101\t60\t8M\t=\t201\t100\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_text)
  let sam_out = @src.bam_to_sam_text(bam)
  // RNEXT 应该输出 "=" 表示 mate 在同一个 ref
  assert_true(sam_out.contains("\t=\t"))
}
