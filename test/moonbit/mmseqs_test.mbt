///|
test "mmseqs db construction" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTAYIAKQR", id="a"))
  ignore(db.add("MKTAYIAKQRQ", id="b"))
  assert_eq(db.size(), 2)
  assert_eq(db.seq(0), Some("MKTAYIAKQR"))
  assert_eq(db.seq_id(1), Some("b"))
  assert_true(db.seq(2) is None)
  assert_true(db.seq_id(-1) is None)
  assert_eq(db.total_residues(), 21)
}

///|
test "mmseqs db duplicate id raises" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MK", id="a"))
  let msg = try {
    ignore(db.add("MK", id="a"))
    "NO_ERROR"
  } catch {
    MmseqsError(m) => m
  }
  assert_true(msg.contains("duplicate sequence id"))
}

///|
test "mmseqs db default id" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKT"))
  ignore(db.add("MKT"))
  assert_eq(db.seq_id(0), Some("seq0"))
  assert_eq(db.seq_id(1), Some("seq1"))
}

///|
test "mmseqs index build k bounds" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTAYIAKQR"))
  let msg = try {
    ignore(@src.MmseqsIndex::build(db, k=2))
    "NO_ERROR"
  } catch {
    MmseqsError(m) => m
  }
  assert_true(msg.contains("k must be in 3..12"))
  let msg2 = try {
    ignore(@src.MmseqsIndex::build(db, k=13))
    "NO_ERROR"
  } catch {
    MmseqsError(m) => m
  }
  assert_true(msg2.contains("k must be in 3..12"))
}

///|
test "mmseqs index kmer exactness" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("ACDEFGHIKL", id="s1"))
  ignore(db.add("XXACDEFGHIKL", id="s2"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  assert_eq(idx.k(), 6)
  // ACDEFG 在两条序列各出现一次 → 2 个位置
  let key = @src.mmseqs_kmer_key_of("ACDEFG", 0, 6).unwrap()
  let hits = idx.table[key]
  assert_eq(hits.length(), 2)
  assert_eq(hits[0].0, 0)
  assert_eq(hits[1].0, 1)
  // 不同 k-mer 不同 key
  assert_true(
    @src.mmseqs_kmer_key_of("ACDEFG", 0, 6) !=
    @src.mmseqs_kmer_key_of("ACDEFH", 0, 6),
  )
  // 小写容错: 同一 key
  assert_eq(
    @src.mmseqs_kmer_key_of("acdefg", 0, 6),
    @src.mmseqs_kmer_key_of("ACDEFG", 0, 6),
  )
  // 非标字符返回 None
  assert_true(@src.mmseqs_kmer_key_of("ACDEFZ", 0, 6) is None)
}

///|
test "mmseqs index num kmers" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("AAAAAA"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  assert_eq(idx.num_kmers(), 1)
}

///|
test "mmseqs prefilter homolog found" {
  let db = @src.MmseqsDB::new()
  // 家族: Walker-A P-loop motif 嵌入不同背景
  ignore(db.add("MKTGHVDHGKTTLKMNL", id="kinase1"))
  ignore(db.add("VWQRTGHVDHGKTTVWQERY", id="kinase2"))
  ignore(db.add("PLMNTGHVDHGKTTPLMNEE", id="kinase3"))
  // 无关随机
  ignore(db.add("WFCYPMIQSVLWRF", id="rand1"))
  ignore(db.add("PQLSMDCVYHTPQL", id="rand2"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  let query = "MKTGHVDHGKTTLKM"
  let pre = @src.mmseqs_prefilter(idx, query, min_score=20.0)
  assert_true(pre.length() >= 3)
  // 最优命中应为家族成员
  let targets : Array[Int32] = []
  for p in pre {
    targets.push(p.target())
  }
  assert_true(targets.contains(0))
  assert_true(targets.contains(1))
  assert_true(targets.contains(2))
  // 排序: 分数降序
  let mut i = 1
  while i < pre.length() {
    assert_true(pre[i - 1].ungapped_score() >= pre[i].ungapped_score())
    i = i + 1
  }
}

///|
test "mmseqs prefilter no hit for unrelated" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTGHVDHGKTTLKMNL", id="kinase1"))
  ignore(db.add("VWQRTGHVDHGKTTVWQERY", id="kinase2"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  let query = "WFCYPMIQSVLWRFAA"
  let pre = @src.mmseqs_prefilter(idx, query, min_score=30.0)
  assert_eq(pre.length(), 0)
}

///|
test "mmseqs search ranks homolog first" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTGHVDHGKTTLKMNL", id="kinase1"))
  ignore(db.add("VWQRTGHVDHGKTTVWQERY", id="kinase2"))
  ignore(db.add("PLMNTGHVDHGKTTPLMNEE", id="kinase3"))
  ignore(db.add("WFCYPMIQSVLWRF", id="rand1"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  let hits = @src.mmseqs_search(idx, "MKTGHVDHGKTTLKM")
  assert_true(hits.length() >= 1)
  let top = hits[0]
  assert_eq(top.target_id(), "kinase1")
  // 完全包含关系 → identity 高
  assert_true(top.identity() > 0.9)
  assert_true(top.evalue() < 1.0e-3)
  assert_true(top.bitscore() > 20.0)
  // E-value 升序
  let mut i = 1
  while i < hits.length() {
    assert_true(hits[i - 1].evalue() <= hits[i].evalue())
    i = i + 1
  }
}

///|
test "mmseqs search no false positive" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("WFCYPMIQSVLWRF", id="rand1"))
  ignore(db.add("PQLSMDCVYHTPQL", id="rand2"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  let hits = @src.mmseqs_search(idx, "MKTAYIAKQRQISV")
  assert_eq(hits.length(), 0)
}

///|
test "mmseqs search evalue threshold filter" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTGHVDHGKTTLKMNL", id="kinase1"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  // 宽松阈值 1e10: 通过 prefilter 的都会保留
  let loose = @src.mmseqs_search(
    idx,
    "MKTGHVDHGKTTLKM",
    evalue_threshold=1.0e10,
  )
  // 严格阈值 1e-30: 全部滤掉
  let strict = @src.mmseqs_search(
    idx,
    "MKTGHVDHGKTTLKM",
    evalue_threshold=1.0e-30,
  )
  assert_true(loose.length() >= 1)
  assert_eq(strict.length(), 0)
}

///|
test "mmseqs search max targets caps hits" {
  let db = @src.MmseqsDB::new()
  let rng = @src.Rng::from_int(7)
  // 10 条同族序列(同一 motif 不同背景)
  for i in 0..<10 {
    let bg : Array[String] = []
    for _j in 0..<6 {
      bg.push("XQEP")
    }
    ignore(
      db.add(
        "MKT" + bg[i % 6] + "GHVDHGKTT" + i.to_string() + "LKMNL",
        id="k" + i.to_string(),
      ),
    )
    ignore(rng.next_int(10))
  }
  let idx = @src.MmseqsIndex::build(db, k=6)
  let hits = @src.mmseqs_search(
    idx,
    "MKTGHVDHGKTTLKMNL",
    evalue_threshold=1.0e10,
    max_targets=3,
  )
  assert_true(hits.length() <= 3)
}

///|
test "mmseqs hit format" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTGHVDHGKTTLKMNL", id="kinase1"))
  let idx = @src.MmseqsIndex::build(db, k=6)
  let hits = @src.mmseqs_search(idx, "MKTGHVDHGKTTLKM")
  assert_true(hits.length() >= 1)
  let f = hits[0].format()
  assert_true(f.contains("query"))
  assert_true(f.contains("target"))
  assert_true(f.contains("kinase1") || f.contains("GHVDHGKTT"))
}

///|
test "mmseqs cluster groups homologs" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKTGHVDHGKTTLKMNLAAAA", id="k1"))
  ignore(db.add("MKTGHVDHGKTTLKMNLBBBB", id="k2"))
  ignore(db.add("MKTGHVDHGKTTLKMNLCCCC", id="k3"))
  ignore(db.add("WFCYPMIQSVLWRFAAWFCYPM", id="r1"))
  ignore(db.add("PQLSMDCVYHTPQLPQLSMDC", id="r2"))
  let assign = @src.mmseqs_cluster(db, coverage=0.5)
  assert_eq(assign.length(), 5)
  // 三条同族同簇
  assert_eq(assign[0], assign[1])
  assert_eq(assign[1], assign[2])
  // 无关序列独立成簇(与家族不同)
  assert_true(assign[3] != assign[0])
  assert_true(assign[4] != assign[0])
  assert_true(assign[3] != assign[4])
}

///|
test "mmseqs cluster coverage validation" {
  let db = @src.MmseqsDB::new()
  ignore(db.add("MKT"))
  let msg = try {
    ignore(@src.mmseqs_cluster(db, coverage=0.0))
    "NO_ERROR"
  } catch {
    MmseqsError(m) => m
  }
  assert_true(msg.contains("coverage must be in"))
  let msg2 = try {
    ignore(@src.mmseqs_cluster(db, k=2))
    "NO_ERROR"
  } catch {
    MmseqsError(m) => m
  }
  assert_true(msg2.contains("k must be in 3..12"))
}

///|
test "mmseqs cluster empty db" {
  let db = @src.MmseqsDB::new()
  let assign = @src.mmseqs_cluster(db)
  assert_eq(assign.length(), 0)
}

///|
test "mmseqs search determinism" {
  let db = @src.MmseqsDB::new()
  let rng = @src.Rng::from_int(42)
  let aas = "ACDEFGHIKLMNPQRSTVWY"
  // 同族: motif + 随机背景
  for i in 0..<5 {
    let bg = StringBuilder()
    for _j in 0..<20 {
      bg.write_char(aas.unsafe_get(rng.next_int(20)).unsafe_to_char())
    }
    ignore(db.add("GHVDHGKTT" + bg.to_string(), id="k" + i.to_string()))
  }
  let idx = @src.MmseqsIndex::build(db, k=6)
  let h1 = @src.mmseqs_search(idx, "GHVDHGKTTMKTAY")
  let h2 = @src.mmseqs_search(idx, "GHVDHGKTTMKTAY")
  assert_eq(h1.length(), h2.length())
  for i in 0..<h1.length() {
    assert_eq(h1[i].target(), h2[i].target())
    assert_eq(h1[i].evalue(), h2[i].evalue())
    assert_eq(h1[i].raw_score(), h2[i].raw_score())
  }
}
