///|
/// HLA 分型演示.
///
/// 展示基于 NGS 读取的 HLA 基因座高分辨率分型流程:
/// 1. 加载 HLA 参考序列数据库
/// 2. 对模拟读取进行分型
/// 3. 输出 GL String 格式结果

///|
/// 主函数: 演示 HLA 分型流程.
fn main {
  let db = @src.hla_reference_db()
  println(
    "HLA 参考数据库: " + db.length().to_string() + " 个等位基因",
  )
  for entry in db {
    println(
      "  - " +
      entry.allele() +
      " (" +
      entry.locus() +
      ", exon " +
      entry.exon().to_string() +
      ")",
    )
  }

  // 模拟读取: 来自 HLA-A*02:01 的片段
  let ref02 = db[0].seq
  let reads : Array[String] = [
    ref02[0:50].to_owned(),
    ref02[50:100].to_owned(),
    ref02[100:150].to_owned(),
    ref02[150:200].to_owned(),
  ]

  println("\n模拟读取数: " + reads.length().to_string())
  println("逐条读取分型结果:")
  for read in reads {
    let call = @src.hla_type_read(read)
    println(
      "  - " +
      call.allele() +
      " (score=" +
      call.score().to_string() +
      ", " +
      call.matches().to_string() +
      "/" +
      call.total().to_string() +
      " match)",
    )
  }

  let gl = @src.hla_type_batch(reads)
  println("\n批量分型结果(GL String): " + gl)

  // 杂合子示例
  let ref01 = db[1].seq
  let het_reads : Array[String] = [
    ref02[0:50].to_owned(),
    ref01[0:50].to_owned(),
    ref02[50:100].to_owned(),
  ]
  let het_gl = @src.hla_type_batch(het_reads)
  println("\n杂合子示例 GL String: " + het_gl)
}
