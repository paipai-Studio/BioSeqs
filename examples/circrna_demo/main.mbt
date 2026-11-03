///|
/// circrna_demo — find_circ 风格环状 RNA 检测演示:
/// 合成 SAM(注入 2 个 circRNA 的 BSJ split-read + 线性融合读段 +
/// 微小环) -> 反向剪接检测 -> 过滤 -> TSV/JSON.

///|
fn crd_read(qname : String, main_pos : Int32, sa_pos : Int32) -> String {
  let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
  "\{qname}\t0\tchr1\t\{main_pos}\t60\t60M\t*\t0\t0\t\{seq}\t*\tSA:Z:chr1,\{sa_pos},+,60M,60,0\n"
}

///|
fn main {
  println("1. find_circ 风格环状 RNA(反向剪接)检测演示")
  let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
  let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n" +
    // circRNA A [1000, 2000): 4 条 BSJ 读段(同断点簇)
    crd_read("a1", 1941, 1001) +
    crd_read("a2", 1941, 1001) +
    crd_read("a3", 1941, 1001) +
    crd_read("a4", 1941, 1001) +
    // circRNA B [3000, 3200): 2 条 BSJ 读段
    crd_read("b1", 3141, 3001) +
    crd_read("b2", 3141, 3001) +
    // 线性融合读段(SA 在下游, 不构成 BSJ) x3
    crd_read("f1", 501, 601) +
    crd_read("f2", 501, 601) +
    crd_read("f3", 501, 601) +
    // 微小环 [4000, 4040) 长 40 < 50 被过滤
    "t1\t0\tchr1\t4031\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tSA:Z:chr1,4001,+,10M,60,0\n" +
    "t2\t0\tchr1\t4031\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tSA:Z:chr1,4001,+,10M,60,0\n" +
    // 单读环 [5000, 5600) support 1 < 2 被过滤
    crd_read("s1", 5541, 5001) +
    // 无 SA 普通读段
    "n1\t0\tchr1\t100\t60\t60M\t*\t0\t0\t\{seq}\t*\n"
  println(
    "2. 合成 SAM: 2 个 circRNA BSJ 读段簇 + 线性融合/微小环/单读干扰",
  )
  // 3. 检测
  let records = @src.circrna_detect(sam)
  println("3. 检测结果:")
  for r in records {
    println("  \{r.summary()}")
  }
  // 4. 导出
  let tsv = @src.circrna_tsv(records)
  let json = @src.circrna_json(records)
  println("4. TSV \{tsv.length()} 字节, JSON \{json.length()} 字节")
  println(@src.circrna_summary_line(records))
  println(
    "演示完成: SAM split-read -> BSJ 判定 -> 聚合过滤 -> circRNA 记录闭环",
  )
}
