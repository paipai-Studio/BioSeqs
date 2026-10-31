///|
/// SV 长读长检测演示.
///
/// 展示基于 minimap 比对的 SV 检出流程:
/// 1. 解析 CIGAR 字符串
/// 2. 从 CIGAR 提取 SV 候选(DEL/INS)
/// 3. 合并相同位置的 SV
/// 4. 输出检测结果

///|
/// 主函数: 演示 SV 检测流程.
fn main {
  println("=== SV 长读长检测演示 ===")

  // 示例 CIGAR 字符串: 10M50D20M 表示 10bp 匹配, 50bp 缺失, 20bp 匹配
  let cigar_str = "10M50D20M"
  println("\nCIGAR 字符串: " + cigar_str)
  let cigar = @src.parse_cigar(cigar_str)
  println("解析结果: " + cigar.length().to_string() + " 个 CIGAR 元素")

  // 检测 SV
  let params = @src.SvParams::default()
  println(
    "\n检测参数: min_support=" +
    params.min_support.to_string() +
    ", min_sv_len=" +
    params.min_sv_len.to_string(),
  )

  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  println("\n检测到的 SV 候选数: " + calls.length().to_string())
  for c in calls {
    println(
      "  - " +
      @src.sv_type_label(c.sv_type()) +
      " at chr1:" +
      c.start().to_string() +
      "-" +
      c.end().to_string() +
      " (support=" +
      c.support_reads().to_string() +
      ")",
    )
  }

  // 插入检测示例
  let ins_cigar = @src.parse_cigar("10M100I20M")
  println("\n插入检测示例: 10M100I20M")
  let ins_calls = @src.sv_candidates_from_sam("chr1", 100, ins_cigar, params)
  for c in ins_calls {
    println(
      "  - " +
      @src.sv_type_label(c.sv_type()) +
      " at chr1:" +
      c.start().to_string() +
      " (support=" +
      c.support_reads().to_string() +
      ")",
    )
  }

  // 合并示例
  println("\n合并相同位置 SV 示例:")
  let all_calls : Array[@src.SvCall] = []
  for c in calls {
    all_calls.push(c)
  }
  for c in calls {
    all_calls.push(c)
  }
  let merged = @src.merge_sv_candidates(all_calls)
  println("  合并前: " + all_calls.length().to_string() + " 个候选")
  println("  合并后: " + merged.length().to_string() + " 个候选")
  if merged.length() > 0 {
    println("  支持读数: " + merged[0].support_reads().to_string())
  }

  println("\n演示完成")
}
