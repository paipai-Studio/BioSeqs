///|
/// P3: pysam_vcf_annotate 单元测试
let annot_vcf_content : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##INFO=<ID=CSQ,Number=A,Type=String,Description=\"consequence\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\trs1\tA\tT\t60.0\tPASS\t.\n" +
  "chr1\t200\trs2\tC\tG\t50.0\tPASS\t.\n" +
  "chr1\t400\trs3\tG\tA\t40.0\tPASS\t.\n"

///|
let annot_bed_text : String = "chr1\t95\t105\texon1\t+1\n" +
  "chr1\t195\t210\texon2\t+2\n" +
  "chr1\t300\t350\tintergenic\t+3\n" +
  "chr2\t50\t60\tchrom2\t+4\n"

///|
test "annotate_parse_bed_annot" {
  let intervals = @src.parse_bed_annot(annot_bed_text)
  assert_eq(intervals.length(), 4)
  assert_eq(intervals[0].chrom, "chr1")
  assert_eq(intervals[0].start, 95)
  assert_eq(intervals[0].end, 105)
  assert_eq(intervals[0].value, "exon1") // 默认 value_col=4
}

///|
test "annotate_parse_bed_value_col_0" {
  let intervals = @src.parse_bed_annot(annot_bed_text, value_col=0)
  assert_eq(intervals[0].value, "chr1:95-105")
}

///|
test "annotate_bed_overlap_mode" {
  let vcf = @src.parse_vcf(annot_vcf_content)
  // rs1 pos=100 落在 exon1 [95,105) 内
  // rs2 pos=200 落在 exon2 [195,210) 内
  // rs3 pos=400 不在任何 bed 区间
  let result = @src.vcf_annotate_bed(
    vcf,
    annot_bed_text,
    "CSQ",
    mode=@src.annot_match_mode_overlap(),
  )
  assert_eq(result.len(), 3)
  assert_eq(result.records[0].info.get("CSQ"), Some("exon1"))
  assert_eq(result.records[1].info.get("CSQ"), Some("exon2"))
  assert_true(result.records[2].info.get("CSQ") is None)
}

///|
test "annotate_bed_exact_mode" {
  let vcf = @src.parse_vcf(annot_vcf_content)
  // Exact 模式要求 pos == bed.start。
  // exon1 start=95, rs1 pos=100 → 不匹配
  let result = @src.vcf_annotate_bed(
    vcf,
    annot_bed_text,
    "CSQ",
    mode=@src.annot_match_mode_exact(),
  )
  assert_eq(result.len(), 3)
  assert_true(result.records[0].info.get("CSQ") is None)
}

///|
test "annotate_from_vcf_copy_info" {
  let anno_text : String = "##fileformat=VCFv4.2\n" +
    "##contig=<ID=chr1,length=1000>\n" +
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"allele freq\">\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tT\t.\t.\tAF=0.45\n" +
    "chr1\t200\t.\tC\tG\t.\t.\tAF=0.33\n" +
    "chr1\t999\t.\tT\tC\t.\t.\tAF=0.10\n"
  let anno = @src.parse_vcf(anno_text)
  let vcf = @src.parse_vcf(annot_vcf_content)
  let result = @src.vcf_annotate_from_vcf(vcf, anno, ["AF"])
  assert_eq(result.len(), 3)
  assert_eq(result.records[0].info.get("AF"), Some("0.45"))
  assert_eq(result.records[1].info.get("AF"), Some("0.33"))
  assert_true(result.records[2].info.get("AF") is None) // 400 != 999
}

///|
test "annotate_source_info_entry" {
  let vcf = @src.parse_vcf(annot_vcf_content)
  let source = @src.vcf_annot_source_from_bed(
    @src.parse_bed_annot(annot_bed_text),
  )
  let result = @src.vcf_annotate_info_from_source(vcf, source, "GENE")
  assert_eq(result.len(), 3)
  assert_eq(result.records[0].info.get("GENE"), Some("exon1"))
}

///|
test "annotate_summary_and_rate" {
  let vcf = @src.parse_vcf(annot_vcf_content)
  let result = @src.vcf_annotate_bed(vcf, annot_bed_text, "CSQ")
  let (ann, unch) = @src.vcf_annotate_summary(vcf, result, "CSQ")
  assert_eq(ann, 2)
  assert_eq(unch, 1)
  assert_true(
    @moonbit.assert_close_e6(@src.vcf_annotate_rate(result, "CSQ"), 2.0 / 3.0),
  )
}

///|
test "annotate_id_from_source" {
  let vcf = @src.parse_vcf(annot_vcf_content)
  let source = @src.vcf_annot_source_from_bed(
    @src.parse_bed_annot(annot_bed_text),
  )
  let result = @src.vcf_annotate_id_from_source(vcf, source)
  assert_eq(result.records[0].id.length(), 1)
  assert_eq(result.records[0].id[0], "exon1")
  assert_eq(result.records[2].id.length(), 1)
  assert_eq(result.records[2].id[0], "rs3") // 没命中 bed，原 id 保留
}
