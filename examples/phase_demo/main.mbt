///|
/// 读段支持的 VCF 定相演示:
/// 杂合 SNP + 跨位点读段 -> 相位约束图 -> phased GT / 单倍型块.
fn main {
  println("=== 读段支持的 VCF 定相 (read-backed phasing) ===")

  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
    "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" +
    "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n" +
    "chr1\t100\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
  )

  // 读段 1: 位点 10=G(alt), 20=C(ref) -> 两等位基因反相.
  // 读段 2: 位点 10=A(ref), 20=T(alt) -> 同上.
  // 读段 3: 位点 10=G(alt), 20=T(alt) -> 与最终相位冲突(模拟测序噪声).
  let reads : Array[@src.PhaseRead] = [
    @src.phase_read_new("r1", "chr1", 0, "40M", seq(40, [(9, 'G'), (19, 'C')])),
    @src.phase_read_new("r2", "chr1", 0, "40M", seq(40, [(9, 'A'), (19, 'T')])),
    @src.phase_read_new("r3", "chr1", 0, "40M", seq(40, [(9, 'G'), (19, 'T')])),
  ]

  let result = @src.phase_vcf_records(vcf, reads, "S1")
  println("\n单倍型块数量 = \{result.blocks}")
  println("读段冲突率   = \{result.conflict_rate}")
  println("\n定相结果:")
  println(@src.phase_result_tsv(result))
  println("JSON 输出:")
  println(@src.phase_result_json(result))
}

///|
fn seq(len : Int32, edits : Array[(Int32, Char)]) -> String {
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
