// ===========================================================================
// pedigree_qc_test.mbt — 家系 .fam 解析 / 孟德尔错误 / 一致性 / ROH 测试.
// ===========================================================================

///|
fn ped_trio_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tC\n" +
  "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\n" +
  "chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/1\n" +
  "chr1\t30\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/0\t0/0\n" +
  "chr1\t40\t.\tT\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n"
}

///|
fn ped_roh_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
  "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t300\t.\tG\tA\t.\t.\t.\tGT\t1/1\n" +
  "chr1\t400\t.\tT\tC\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/1\n" +
  "chr1\t600\t.\tC\tT\t.\t.\t.\tGT\t0/0\n"
}

///|
fn ped_fam_text() -> String {
  "F1 F 0 0 1 1\n" + "F1 M 0 0 2 1\n" + "F1 C F M 1 2\n"
}

// ===========================================================================
// .fam 解析 / trio 提取 / 校验
// ===========================================================================

///|
test "pedigree_parse_reads_individuals" {
  let pedigree = @src.pedigree_parse(ped_fam_text())
  assert_eq(pedigree.individuals.length(), 3)
  let child = pedigree.individuals[2]
  assert_eq(child.family_id, "F1")
  assert_eq(child.id, "C")
  assert_eq(child.father, Some("F"))
  assert_eq(child.mother, Some("M"))
  assert_eq(child.sex, 1)
  assert_eq(child.phenotype, "2")
  // "0" 父/母视为缺失.
  assert_true(pedigree.individuals[0].father is None)
}

///|
test "pedigree_trios_extracts_core_families" {
  let pedigree = @src.pedigree_parse(ped_fam_text())
  let trios = @src.pedigree_trios(pedigree)
  assert_eq(trios.length(), 1)
  assert_eq(trios[0].child, "C")
  assert_eq(trios[0].father, "F")
  assert_eq(trios[0].mother, "M")
}

///|
test "pedigree_validate_accepts_consistent_pedigree" {
  let pedigree = @src.pedigree_parse(ped_fam_text())
  let ok = try {
    @src.pedigree_validate(pedigree)
    true
  } catch {
    @src.PedigreeError(_) => false
  }
  assert_true(ok)
}

///|
test "pedigree_validate_rejects_unknown_parent" {
  let bad = @src.pedigree_parse("F1 C X Y 1 1\n")
  let raised = try {
    @src.pedigree_validate(bad)
    false
  } catch {
    @src.PedigreeError(_) => true
  }
  assert_true(raised)
}

// ===========================================================================
// 孟德尔错误
// ===========================================================================

///|
test "pedigree_mendelian_check_flags_inconsistent_site" {
  let vcf = @src.parse_vcf(ped_trio_vcf())
  let pedigree = @src.pedigree_parse(ped_fam_text())
  let trio = @src.pedigree_trios(pedigree)[0]
  let qc = @src.pedigree_mendelian_check(vcf, trio)
  assert_eq(qc.sites_checked, 4)
  assert_eq(qc.errors.length(), 1)
  let err = qc.errors[0]
  assert_eq(err.chrom, "chr1")
  assert_eq(err.pos, 40)
  assert_eq(err.child, "C")
  assert_eq(err.father_gt, "0/0")
  assert_eq(err.mother_gt, "0/0")
  assert_eq(err.child_gt, "0/1")
}

///|
test "pedigree_check_all_summarizes" {
  let vcf = @src.parse_vcf(ped_trio_vcf())
  let pedigree = @src.pedigree_parse(ped_fam_text())
  let all = @src.pedigree_check_all(vcf, pedigree)
  assert_eq(all.length(), 1)
  assert_true(
    @src.pedigree_qc_summary(vcf, pedigree).contains("mendelian_errors=1"),
  )
}

// ===========================================================================
// 一致性
// ===========================================================================

///|
test "pedigree_concordance_counts_matches" {
  let vcf = @src.parse_vcf(ped_trio_vcf())
  let result = @src.pedigree_concordance(vcf, "F", "M")
  assert_eq(result.sites_compared, 4)
  assert_eq(result.concordant, 2)
  assert_eq(result.discordant, 2)
  assert_true((result.concordance - 0.5).abs() < 1.0e-9)
}

///|
test "pedigree_concordance_self_is_perfect" {
  let vcf = @src.parse_vcf(ped_trio_vcf())
  let result = @src.pedigree_concordance(vcf, "F", "F")
  assert_eq(result.concordance, 1.0)
}

// ===========================================================================
// ROH
// ===========================================================================

///|
test "pedigree_roh_detects_homozygous_run" {
  let vcf = @src.parse_vcf(ped_roh_vcf())
  let result = @src.pedigree_roh(vcf, sample="S", min_markers=3, max_gap=1000)
  assert_eq(result.total_markers, 6)
  assert_eq(result.segments.length(), 1)
  let seg = result.segments[0]
  assert_eq(seg.chrom, "chr1")
  assert_eq(seg.start, 100)
  assert_eq(seg.end, 400)
  assert_eq(seg.n_markers, 4)
}

///|
test "pedigree_roh_respects_min_markers" {
  let vcf = @src.parse_vcf(ped_roh_vcf())
  let result = @src.pedigree_roh(vcf, sample="S", min_markers=5, max_gap=1000)
  assert_eq(result.segments.length(), 0)
}

// ===========================================================================
// JSON 汇总
// ===========================================================================

///|
test "pedigree_qc_report_json_contains_sections" {
  let vcf = @src.parse_vcf(ped_trio_vcf())
  let pedigree = @src.pedigree_parse(ped_fam_text())
  let json = @src.pedigree_qc_report_json(pedigree, vcf)
  assert_true(json.contains("\"trios\""))
  assert_true(json.contains("\"mendelian_errors\":1"))
  assert_true(json.contains("\"roh\""))
  assert_true(@src.parse_json(json) is Some(_))
}
