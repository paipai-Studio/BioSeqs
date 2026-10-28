///|
/// BEDPE 交互对格式演示: 解析/统计/区间判断.
fn main raise {
  let content = "chr1\t1000\t2000\tchr1\t5000\t6000\tloop1\t780\t+\t-\nchr1\t3000\t3500\tchr1\t8000\t9000\tloop2\t520\t.\t.\nchr1\t100\t200\tchr2\t500\t600\tsv_deletion\t900\t+\t-\n"
  println("=== 解析 BEDPE ===")
  let pairs = @src.parse_bedpe(content)
  println("  交互对数: " + pairs.length().to_string())
  for p in pairs {
    println(
      "  " +
      p.name +
      ": " +
      p.chrom1 +
      ":" +
      p.start1.to_string() +
      "-" +
      p.end1.to_string() +
      " <-> " +
      p.chrom2 +
      ":" +
      p.start2.to_string() +
      "-" +
      p.end2.to_string() +
      " score=" +
      p.score.to_string(),
    )
  }
  println("=== 同染色体过滤与 loop 长度 ===")
  let same = @src.bedpe_filter_same_chrom(pairs)
  for p in same {
    println("  " + p.name + " loop 长度: " + p.loop_length().to_string())
  }
  println("=== 染色体对计数 ===")
  let counts = @src.bedpe_count_chrom_pairs(pairs)
  for entry in counts {
    println("  " + entry.0 + ": " + entry.1.to_string())
  }
  println("=== 距离直方图(0-5000, 步长 1000) ===")
  let hist = @src.bedpe_distance_histogram(same, 0, 1000, 5)
  for i, n in hist {
    println(
      "  [" +
      (i * 1000).to_string() +
      ", " +
      ((i + 1) * 1000).to_string() +
      "): " +
      n.to_string(),
    )
  }
  println("=== 区间重叠判断 ===")
  let p0 = pairs[0]
  println(
    "  loop1 与 chr1:5500-6500 重叠: " +
    p0.overlaps_region("chr1", 5500, 6500).to_string(),
  )
  println(
    "  loop1 与 chr1:2500-3500 重叠: " +
    p0.overlaps_region("chr1", 2500, 3500).to_string(),
  )
  println("=== 写出往返 ===")
  let out = @src.write_bedpe(pairs)
  let reparsed = @src.parse_bedpe(out)
  println("  往返条数: " + reparsed.length().to_string())
  println("  loop1 score 保持: " + reparsed[0].score.to_string())
}
