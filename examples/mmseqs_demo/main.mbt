///|
/// mmseqs2 风格本地蛋白质搜索演示: 合成激酶家族库 → k-mer 索引 →
/// prefilter(ungapped x-drop) → SW 精化 + E-value → linclust 预聚类.
fn main raise {
  let rng = @src.Rng::from_int(42)
  let aas = "ACDEFGHIKLMNPQRSTVWY"
  // 1. 合成库: 15 条激酶家族(Walker-A motif 嵌入随机背景) + 85 条随机
  let db = @src.MmseqsDB::new()
  for i in 0..<15 {
    let bg = StringBuilder()
    for _j in 0..<30 {
      bg.write_char(aas.unsafe_get(rng.next_int(20)).unsafe_to_char())
    }
    let ed = StringBuilder()
    for _j in 0..<30 {
      ed.write_char(aas.unsafe_get(rng.next_int(20)).unsafe_to_char())
    }
    ignore(
      db.add(
        "MK" + bg.to_string() + "GHVDHGKTT" + ed.to_string() + "LMNL",
        id="kinase" + i.to_string(),
      ),
    )
  }
  for i in 0..<85 {
    let rs = StringBuilder()
    for _j in 0..<60 {
      rs.write_char(aas.unsafe_get(rng.next_int(20)).unsafe_to_char())
    }
    ignore(db.add(rs.to_string(), id="rand" + i.to_string()))
  }
  let query = "MKTGHVDHGKTTLKMNLVEG"
  println("=== mmseqs2 风格蛋白质搜索 (库 100 条, k=6) ===")
  println("  库残基总数: " + db.total_residues().to_string())
  // 2. 建索引
  let idx = @src.MmseqsIndex::build(db, k=6)
  println(
    "  k-mer 索引: " +
    idx.num_kmers().to_string() +
    " 个不同 6-mer, k=" +
    idx.k().to_string(),
  )
  // 3. prefilter: ungapped x-drop 快速筛查
  let pre = @src.mmseqs_prefilter(idx, query)
  println("=== prefilter (ungapped x-drop, 分数降序 top5) ===")
  let top_n = if pre.length() < 5 { pre.length() } else { 5 }
  for i in 0..<top_n {
    let p = pre[i]
    let tid = db.seq_id(p.target()).unwrap()
    println(
      "  " +
      tid +
      "  score " +
      p.ungapped_score().to_string() +
      "  seeds " +
      p.seeds().to_string() +
      "  diag " +
      p.diagonal().to_string(),
    )
  }
  // 4. SW 精化 + E-value
  let hits = @src.mmseqs_search(idx, query)
  println("=== SW 精化命中 (E-value 升序 top5) ===")
  let hit_n = if hits.length() < 5 { hits.length() } else { 5 }
  for i in 0..<hit_n {
    let h = hits[i]
    println(
      "  " +
      h.target_id() +
      "  raw " +
      h.raw_score().to_string() +
      "  bits " +
      h.bitscore().to_string() +
      "  E " +
      h.evalue().to_string() +
      "  id " +
      (h.identity() * 100.0).to_string() +
      "%",
    )
  }
  if hits.length() > 0 {
    println("=== 最优命中比对 ===")
    println(hits[0].format())
  }
  // 5. 预聚类(独立演示库: 10 个"物种"基序列 × 各 10 条近重复拷贝,
  //    linclust 的本用途是去冗余而非远程同源聚类)
  let cdb = @src.MmseqsDB::new()
  for s in 0..<10 {
    let base = StringBuilder()
    for _j in 0..<50 {
      base.write_char(aas.unsafe_get(rng.next_int(20)).unsafe_to_char())
    }
    let bs = base.to_string()
    for c in 0..<10 {
      let arr : Array[Int32] = []
      for ch in bs {
        arr.push(ch.to_int())
      }
      let mpos = rng.next_int(arr.length())
      arr[mpos] = aas.unsafe_get(rng.next_int(20)).to_int()
      let sb = StringBuilder()
      for v in arr {
        sb.write_char(Char::from_int(v))
      }
      ignore(
        cdb.add(sb.to_string(), id="sp" + s.to_string() + "_" + c.to_string()),
      )
    }
  }
  // 代表自身含 1 突变 + 拷贝 1 突变 → 拷贝与代表差 2 位,
  // shared ≈ 33/45 ≈ 0.73, coverage 取 0.7
  let assign = @src.mmseqs_cluster(cdb, coverage=0.7)
  let mut nclusters = 0
  for a in assign {
    if a + 1 > nclusters {
      nclusters = a + 1
    }
  }
  let mut first_size = 0
  for a in assign {
    if a == assign[0] {
      first_size = first_size + 1
    }
  }
  println("=== linclust 风格预聚类 (10 物种 × 10 近重复拷贝) ===")
  println(
    "  100 条序列 → " +
    nclusters.to_string() +
    " 簇; 首簇 " +
    first_size.to_string() +
    " 条 (coverage=0.7)",
  )
  println("")
  println(
    "注: 本模块为 mmseqs2 风格筛查级实现(5bit 精确 k-mer/单点对角",
  )
  println("    延伸/ungapped 统计 E-value), 非其替身.")
}
