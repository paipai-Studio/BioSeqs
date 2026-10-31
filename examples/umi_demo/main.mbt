///|
/// UMI 分子标签演示: header 提取 -> 5' 剪除 -> 去重 -> 纠错 -> 共识.
fn main raise {
  println(
    "=== UMI 提取 / 去重 / 纠错 (umi_tools 风格, 本地完成) ===",
  )

  // 模拟测序 read: UMI 通过下划线后缀携带, 5' 端含 UMI 碱基.
  let headers = ["@read1_ACGTACGT", "@read2_ACGTACGT", "@read3_ACGAACGT"]
  let raw_seqs = ["ACGTACGTAAAA", "ACGTACGTAATA", "ACGAACGTAAAT"]
  let raw_quals = ["IIIIIIIIIIII", "IIIIIIIIIIII", "IIIIIIIIIIII"]

  println("\n1. header 解析 + UMI 剪除")
  let reads : Array[@src.UmiRead] = Array(capacity=0)
  for i = 0; i < headers.length(); i = i + 1 {
    let parsed = @src.umi_parse_header(headers[i], @src.Underscore)
    let q = @src.umi_qualities_from_ascii(raw_quals[i])
    let (seq, qual) = @src.umi_clip_read(raw_seqs[i], q, parsed.umi.length())
    println("  \{parsed.name}: UMI=\{parsed.umi}, 剪除后=\{seq}")
    reads.push(@src.UmiRead::{
      id: parsed.name,
      sequence: seq,
      qualities: qual,
      umi: parsed.umi,
      chrom: "chr1",
      pos: 100,
      strand: 1,
    })
  }

  println("\n2. 按 (UMI, locus) 去重并生成共识")
  let families = @src.umi_dedup(reads)
  for f in families {
    println(
      "  UMI=\{f.umi} pos=\{f.pos} family_size=\{f.size} consensus=\{f.consensus}",
    )
  }

  println("\n3. directional UMI 纠错 (编辑距离 <= 1)")
  let corrected = @src.umi_directional_correct(families)
  for f in corrected {
    println("  代表 UMI=\{f.umi} 合并后 size=\{f.size}")
  }

  println("\n4. 汇总统计")
  let stats = @src.umi_stats(corrected)
  println(
    "  families=\{stats.n_families}, reads=\{stats.n_reads}, singletons=\{stats.singleton_families}, mean_quality=\{stats.mean_quality}",
  )

  println("\n5. JSON 输出")
  println(@src.umi_families_to_json(corrected))

  println("完成: UMI 链路全程本地处理, 数据不出本机.")
}
