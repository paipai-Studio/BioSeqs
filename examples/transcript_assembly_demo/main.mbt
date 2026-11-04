///|
/// transcript_assembly_demo — StringTie 风格转录本组装演示:
/// 合成 SAM(共享外显子的双 isoform + 未剪接读段) ->
/// junction 提取 -> 剪接图贪心流 -> isoform(GTF/FPKM) -> JSON.

///|
fn td_sam() -> String {
  // E1 [10,100), E2 [200,300), E3 [400,440), E4 [600,640)
  // junction (100,200) x12, (300,400) x4, (300,600) x2, E1 未剪接 x20
  let sb = StringBuilder()
  sb.write_string("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n")
  let long_seq : String = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" +
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" +
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" +
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
  for i = 0; i < 12; i = i + 1 {
    sb.write_string(
      "j1_\{i}\t0\tchr1\t81\t60\t20M100N40M\t*\t0\t0\t\{long_seq}\t*\n",
    )
  }
  for i = 0; i < 4; i = i + 1 {
    sb.write_string(
      "j2_\{i}\t0\tchr1\t281\t60\t20M100N40M\t*\t0\t0\t\{long_seq}\t*\n",
    )
  }
  for i = 0; i < 2; i = i + 1 {
    sb.write_string(
      "j3_\{i}\t0\tchr1\t281\t60\t20M300N40M\t*\t0\t0\t\{long_seq}\t*\n",
    )
  }
  for i = 0; i < 20; i = i + 1 {
    sb.write_string("u_\{i}\t0\tchr1\t11\t60\t90M\t*\t0\t0\t\{long_seq}\t*\n")
  }
  sb.to_string()
}

///|
fn main {
  let sam = td_sam()
  println(
    "1. 合成 SAM: 38 条读段(共享外显子双 isoform + 未剪接覆盖)",
  )
  // 2. junction 提取
  let junctions = @src.ta_junctions_from_sam("chr1", sam)
  println("2. junction 提取: \{junctions.length()} 条")
  for j in junctions {
    println("  (\{j.donor}, \{j.acceptor}) x \{j.count}")
  }
  // 3. 组装
  let ts = @src.ta_assemble_transcripts("chr1", sam, strand="+")
  println("3. 组装 \{ts.length()} 条 isoform:")
  for t in ts {
    let mut ex = ""
    for e in t.exons {
      ex = ex + "[\{e.0}, \{e.1}) "
    }
    println(
      "  \{t.id}: exons \{ex}| cov=\{t.coverage} | reads=\{t.read_count} | fpkm=\{t.fpkm}",
    )
  }
  // 4. 汇总与 GTF
  let stats = @src.ta_assembly_stats("chr1", sam, ts)
  println(
    "4. 汇总: reads=\{stats.n_reads}, junctions=\{stats.n_junctions}, transcripts=\{stats.n_transcripts}",
  )
  println(@src.ta_assembly_stats_line(stats))
  let gtf = @src.ta_transcripts_to_gtf(ts)
  println("5. GTF (前 6 行, 共 \{gtf.length()} 字节):")
  let mut shown = 0
  for line in gtf.split("\n") {
    if line != "" && shown < 6 {
      println("  " + line.to_owned())
      shown = shown + 1
    }
  }
  // 6. JSON
  let json = @src.ta_transcripts_json(ts)
  println("6. JSON 导出: \{json.length()} 字节")
  println(
    "演示完成: SAM -> junction -> 剪接图流 -> isoform GTF/FPKM 闭环",
  )
}
