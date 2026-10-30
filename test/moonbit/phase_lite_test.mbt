// ===========================================================================
// phase_lite_test.mbt — 读段支持的 VCF 定相测试.
// ===========================================================================

///|
fn phase_lite_test_seq(len : Int32, edits : Array[(Int32, Char)]) -> String {
  let chars : Array[Char] = Array::make(len, 'A')
  for edit in edits {
    chars[edit.0] = edit.1
  }
  let sb = StringBuilder(size_hint=len)
  for c in chars {
    sb.write_char(c)
  }
  sb.to_string()
}

///|
fn phase_lite_test_vcf() -> @src.VcfFile {
  @src.parse_vcf(
    "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
    "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" +
    "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n" +
    "chr1\t100\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
  )
}

///|
fn phase_lite_test_reads() -> Array[@src.PhaseRead] {
  [
    @src.phase_read_new(
      "r1",
      "chr1",
      0,
      "40M",
      phase_lite_test_seq(40, [(9, 'G'), (19, 'C')]),
    ),
    @src.phase_read_new(
      "r2",
      "chr1",
      0,
      "40M",
      phase_lite_test_seq(40, [(9, 'A'), (19, 'T')]),
    ),
    @src.phase_read_new(
      "r3",
      "chr1",
      0,
      "40M",
      phase_lite_test_seq(40, [(9, 'G'), (19, 'T')]),
    ),
  ]
}

///|
test "phase_lite_builds_block_and_phases_snps" {
  let result = @src.phase_vcf_records(
    phase_lite_test_vcf(),
    phase_lite_test_reads(),
    "S1",
  )
  assert_eq(result.blocks, 1)
  assert_eq(result.calls.length(), 3)
  assert_eq(result.calls[0].chrom, "chr1")
  assert_eq(result.calls[0].pos, 10)
  assert_eq(result.calls[0].genotype, "1|0")
  assert_eq(result.calls[0].block, Some(0))
  assert_eq(result.calls[1].genotype, "0|1")
  assert_eq(result.calls[1].block, Some(0))
  // 无读段覆盖的位点保持未定相.
  assert_eq(result.calls[2].genotype, "0/1")
  assert_eq(result.calls[2].block, None)
}

///|
test "phase_lite_tracks_read_conflicts" {
  let result = @src.phase_vcf_records(
    phase_lite_test_vcf(),
    phase_lite_test_reads(),
    "S1",
  )
  // 3 条读段各产生 1 条约束,其中 r3 与最终相位冲突.
  assert_true((result.conflict_rate - 1.0 / 3.0).abs() < 1.0e-9)
}

///|
test "phase_lite_softclip_projection" {
  // 前 5bp soft clip 后从参考 0 起 30M,应仍覆盖两个位点.
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
    "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" +
    "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
  )
  let reads : Array[@src.PhaseRead] = [
    @src.phase_read_new(
      "clip",
      "chr1",
      0,
      "5S30M",
      phase_lite_test_seq(35, [(14, 'G'), (24, 'C')]),
    ),
  ]
  let result = @src.phase_vcf_records(vcf, reads, "S1")
  assert_eq(result.blocks, 1)
  assert_eq(result.calls[0].genotype, "1|0")
  assert_eq(result.calls[1].genotype, "0|1")
}

///|
test "phase_lite_skips_missing_sample" {
  let result = @src.phase_vcf_records(
    phase_lite_test_vcf(),
    phase_lite_test_reads(),
    "NotThere",
  )
  assert_eq(result.blocks, 0)
  for call in result.calls {
    assert_eq(call.block, None)
  }
}

///|
test "phase_lite_json_and_tsv" {
  let result = @src.phase_vcf_records(
    phase_lite_test_vcf(),
    phase_lite_test_reads(),
    "S1",
  )
  let json = @src.phase_result_json(result)
  assert_true(json.contains("\"conflict_rate\""))
  assert_true(json.contains("\"genotype\":\"1|0\""))
  assert_true(@src.parse_json(json) is Some(_))

  let tsv = @src.phase_result_tsv(result)
  let lines = tsv.split("\n").to_array()
  assert_eq(lines[0].to_owned(), "chrom\tpos\tref\talt\tgenotype\tblock")
  assert_eq(lines.length(), 5)
}
