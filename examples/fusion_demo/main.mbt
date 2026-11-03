///|
/// 融合基因 caller 演示: SAM SA 标签 -> split-read 断点簇 -> 融合候选.
fn main raise {
  println(
    "=== 融合基因 caller (split-read 断点簇, STAR-Fusion 风格) ===",
  )

  // 模拟 SAM: 6 条 split read
  // split1-3: chr1:1016 -- chr2:5016 (融合 A, 支持 3)
  // split4-5: chr1:2016 -- chr2:8016 (融合 B, 支持 2)
  // split6: 单读(低于支持阈值)
  let sam_header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:1000000\n"
  let mk = fn(
    qname : String,
    chrom : String,
    pos : Int32,
    sa : String,
  ) -> String {
    qname +
    "\t0\t" +
    chrom +
    "\t" +
    pos.to_string() +
    "\t60\t16M16S\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tSA:Z:" +
    sa
  }
  let reads : Array[String] = [
    mk("split1", "chr1", 1000, "chr2,5000,+,16S16M,60,0;"),
    mk("split2", "chr1", 1000, "chr2,5000,+,16S16M,60,0;"),
    mk("split3", "chr1", 1005, "chr2,5005,+,16S16M,60,0;"),
    mk("split4", "chr1", 2000, "chr2,8000,+,16S16M,60,0;"),
    mk("split5", "chr1", 2000, "chr2,8000,+,16S16M,60,0;"),
    mk("split6", "chr1", 5000, "chr2,9000,+,16S16M,60,0;"),
  ]
  let sb = StringBuilder()
  sb.write_string(sam_header)
  sb.write_char('\n')
  for r in reads {
    sb.write_string(r)
    sb.write_char('\n')
  }
  let sam = sb.to_string()
  println("SAM 输入: \{reads.length()} 条 split read")

  // 1. 断点对提取
  println("\n1. split-read 断点对提取(SA 标签)")
  let js = @src.fusion_junctions_from_sam(sam)
  println("  共 \{js.length()} 个 junction")
  for j in js {
    println(
      "  \{j.qname}: \{j.chrom_a}:\{j.break_a} -- \{j.chrom_b}:\{j.break_b}",
    )
  }

  // 2. 断点簇聚类与候选
  println("\n2. 融合候选(断点簇, min_support=2)")
  let candidates = @src.fusion_call(js, min_support=2)
  println("  " + @src.fusion_candidates_summary(candidates))
  for c in candidates {
    println("  " + c.summary())
  }

  // 3. read-through 过滤演示
  println("\n3. read-through 过滤(同链近距)")
  let rt_js : Array[@src.FusionJunction] = []
  for i = 0; i < 5; i = i + 1 {
    rt_js.push({
      qname: "rt" + i.to_string(),
      chrom_a: "chr1",
      break_a: 100000,
      strand_a: 0,
      chrom_b: "chr1",
      break_b: 150000,
      strand_b: 0,
    })
  }
  let rt_candidates = @src.fusion_call(
    rt_js,
    min_support=2,
    readthrough_max=200000,
  )
  for c in rt_candidates {
    println("  " + c.summary() + " (被过滤)")
  }

  // 4. TSV 输出
  println("\n4. STAR-Fusion 风格 TSV")
  let tsv = @src.fusion_write_tsv(candidates)
  let lines = tsv.split("\n").to_array()
  for l in lines {
    let s = l.to_owned()
    if s.length() > 0 {
      println("  " + s)
    }
  }

  // 5. BEDPE 输出(与 parse_bedpe 闭环)
  println("\n5. BEDPE 输出(与 bedpe.mbt 闭环)")
  let bedpe = @src.fusion_write_bedpe(candidates)
  let content_lines = bedpe.split("\n").to_array()
  let filtered : Array[String] = []
  for l in content_lines {
    let s = l.to_owned()
    if s.length() > 0 && !s.has_prefix("#") {
      filtered.push(s)
    }
  }
  let joined : String = filtered.join("\n") + "\n"
  let pairs = @src.parse_bedpe(joined)
  println("  BEDPE 解析: \{pairs.length()} 条记录闭环")

  println(
    "\n完成: SA 解析 -> 断点簇 -> 融合候选 -> TSV/BEDPE 全链路本地完成.",
  )
}
