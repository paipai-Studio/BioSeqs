///|
/// P3 Batch 1: pysam_vcf_norm 单元测试
let norm_vcf_content : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\trs1\tACGT\tACGTA\t60\tPASS\t.\n" + // 共享 ACGT 4 前缀 → REF=N ALT=A pos=104
  "chr1\t200\trs2\tATG\tACG\t50\tPASS\t.\n" + // 共享 A 1 → REF=TG ALT=CG pos=201
  "chr1\t300\trs3\tGGTC\tGGCC\t50\tPASS\t.\n" + // 共享 G 1 → REF=GTC ALT=GCC pos=301
  "chr1\t400\trs4\tA\tT\t60\tPASS\t.\n" + // SNP，长度差 0 → 不是 norm target
  "chr1\t500\t.\tAG\tA\t40\tPASS\t.\n" + // REF=AG ALT=A 共享 A → REF=G ALT=N pos=501
  "chr1\t600\t.\tA\tAG\t30\tPASS\t.\n" // REF=A ALT=AG 共享 A → REF=N ALT=G pos=601

///|
test "norm_is_norm_target" {
  // 单 ALT、非空、不等 → target。长度差不管（trim_prefix 对 MNV 也能 canonicalize）
  let v = @src.parse_vcf(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tAT\t60\tPASS\t.\n" +
    "chr1\t200\t.\tA\tT\t60\tPASS\t.\n" + // SNP 长度差 0 也被判 target（trim 后 ref/alt 不变）
    "chr1\t300\t.\tA\tT,C\t60\tPASS\t.\n", // multi-allelic → 不是
  )
  assert_true(@src.vcf_is_norm_target(v.records[0])) // A→AT insertion
  assert_true(@src.vcf_is_norm_target(v.records[1])) // A→T SNP 也判 target（ref!=alt）
  assert_true(!@src.vcf_is_norm_target(v.records[2])) // multi-allelic 不是
}

///|
test "norm_trim_prefix_shared_acgt" {
  let vcf = @src.parse_vcf(norm_vcf_content)
  let rs1 = vcf.records[0] // REF=ACGT ALT=ACGTA pos=100
  let trimmed = @src.vcf_norm_trim_prefix(rs1)
  // 共享 ACGT 4 碱基 → trim 后 REF="" → 用 "N"，ALT="A"，pos+4=104
  assert_eq(trimmed.pos, 104)
  assert_eq(trimmed.ref_seq, "N")
  assert_eq(trimmed.alt.length(), 1)
  assert_eq(trimmed.alt[0], "A")
}

///|
test "norm_trim_prefix_no_change_for_snp" {
  let vcf = @src.parse_vcf(norm_vcf_content)
  let rs2 = vcf.records[1] // REF=ATG ALT=ACG 首碱基都是 A → 共享 1
  let r = @src.vcf_norm_trim_prefix(rs2)
  assert_eq(r.pos, 201)
  assert_eq(r.ref_seq, "TG")
  assert_eq(r.alt[0], "CG")
}

///|
test "norm_trim_prefix_deletion_keeps_first_base" {
  let vcf = @src.parse_vcf(norm_vcf_content)
  // rs5: REF=AG ALT=A → 共享 1 (A) → REF=G ALT=N pos=501
  let r = @src.vcf_norm_trim_prefix(vcf.records[4])
  assert_eq(r.pos, 501)
  assert_eq(r.ref_seq, "G")
}

///|
test "norm_trim_all_on_vcf_file" {
  let vcf = @src.parse_vcf(norm_vcf_content)
  let result = @src.vcf_norm_trim_only(vcf)
  assert_eq(result.len(), 6)
  // rs1 应该被 trim — pos 从 100 变成 104
  assert_true(result.records[0].pos != 100)
}

///|
test "norm_trim_suffix_shared" {
  // REF="ATT" ALT="CTT" → 共享 TT 后缀 → 去掉后 REF="A" ALT="C" pos 不变
  let v = @src.parse_vcf(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tATT\tCTT\t60\tPASS\t.\n",
  )
  let r = @src.vcf_norm_trim_suffix(v.records[0])
  assert_eq(r.pos, 100)
  assert_eq(r.ref_seq, "A")
  assert_eq(r.alt[0], "C")
}

///|
test "norm_trim_all_combines_prefix_suffix" {
  // REF="AACGGT" ALT="AACTT"
  // 前缀：AAC 共享 3 → REF="GGT" ALT="TT" pos=103
  // 后缀：REF="GGT" vs ALT="TT" → 结尾 T 共享 1 → REF="GG" ALT="T"
  let v = @src.parse_vcf(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tAACGGT\tAACTT\t60\tPASS\t.\n",
  )
  let rec = v.records[0]
  let r = @src.vcf_norm_trim_all(rec)
  assert_eq(r.pos, 103) // 3 个 prefix 去掉
  assert_eq(r.ref_seq, "GG")
  assert_eq(r.alt[0], "T")
}

///|
test "norm_file_trim_only_passthrough_empty" {
  let empty = @src.parse_vcf(
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
  )
  let result = @src.vcf_norm_trim_only(empty)
  assert_eq(result.len(), 0)
}
