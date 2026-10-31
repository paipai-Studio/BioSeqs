// ===========================================================================
// kallisto_lite_test.mbt — 伪比对转录本定量测试
// ===========================================================================

///|
fn kallisto_test_txs() -> (Array[String], Array[String]) {
  // t1/t2 共享中段(SHARED), 两端独特; t3 全独特
  let names : Array[String] = ["t1", "t2", "t3"]
  let seqs : Array[String] = [
    "AAAAACCCCCGGGGGTTTTT" + "ATATGCGCGCTATAT" + "AACCTTGGAACCTTGGAA",
    "TTTTTGGGGGCCCCCAAAAA" + "ATATGCGCGCTATAT" + "GGTTCCAAGGTTCCAAGG",
    "ACACACGTGTGTATATATCGCGCGTATATCGCGCGTATAT",
  ]
  (names, seqs)
}

///|
test "kallisto build index" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  assert_eq(index.k, 9)
  assert_eq(index.tx_names.length(), 3)
  // 有效长度 = len - k + 1
  assert_true(index.eff_lens[0] > 44.0 && index.eff_lens[0] < 46.0)
}

///|
test "kallisto build index mismatch raises" {
  let names : Array[String] = ["t1", "t2"]
  let seqs : Array[String] = ["ACGTACGT"]
  let mut raised = false
  try {
    let idx = @src.kallisto_build_index(names, seqs, k=9)
    ignore(idx)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "kallisto pseudoalign unique read" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // t3 独特区 read -> 只兼容 t3
  let read = "ACACACGTGTGTATATATCGCGCG"
  let ids = @src.kallisto_pseudoalign_read(index, read)
  assert_eq(ids, [2])
  // t1 独特区 read -> 只兼容 t1
  let read1 = "AAAAACCCCCGGGGGTTTTTATAT"
  let ids1 = @src.kallisto_pseudoalign_read(index, read1)
  assert_eq(ids1, [0])
}

///|
test "kallisto pseudoalign shared read ambiguous" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // 共享区 read -> 兼容 {t1, t2}
  let read = "ATATGCGCGCTATAT"
  let ids = @src.kallisto_pseudoalign_read(index, read)
  assert_eq(ids, [0, 1])
}

///|
test "kallisto pseudoalign unmapped read" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // 索引中不存在的 k-mer
  let read = "TTTTAAAAACCCCRRRRGGGGTTTTAAAAACCCC"
  let ids = @src.kallisto_pseudoalign_read(index, read)
  assert_eq(ids.length(), 0)
}

///|
test "kallisto equivalence classes" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  let reads : Array[String] = [
    "AAAAACCCCCGGGGGTTTTTATAT", // t1 unique
     "AAAAACCCCCGGGGGTTTTTATAT", // t1 unique x2
     "ATATGCGCGCTATAT", // shared
     "ACACACGTGTGTATATATCGCGCG", // t3 unique
  ]
  let ecs = @src.kallisto_pseudoalign(index, reads)
  assert_eq(ecs.length(), 3)
  // 单转录本类在前
  assert_eq(ecs[0].tx_ids, [0])
  assert_eq(ecs[0].count, 2)
  assert_eq(ecs[1].tx_ids, [2])
  assert_eq(ecs[1].count, 1)
  // 多转录本类在后
  assert_eq(ecs[2].tx_ids, [0, 1])
  assert_eq(ecs[2].count, 1)
}

///|
test "kallisto em assigns unique reads exactly" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // 全部独特 read: EM 应精确分配
  let reads : Array[String] = []
  for i = 0; i < 10; i = i + 1 {
    reads.push("AAAAACCCCCGGGGGTTTTTATAT")
  }
  for i = 0; i < 5; i = i + 1 {
    reads.push("ACACACGTGTGTATATATCGCGCG")
  }
  let quant = @src.kallisto_quantify(index, reads)
  assert_eq(quant.n_pseudoaligned, 15)
  assert_eq(quant.n_ec, 2)
  assert_true(quant.est_counts[0] > 9.99 && quant.est_counts[0] < 10.01)
  assert_true(quant.est_counts[2] > 4.99 && quant.est_counts[2] < 5.01)
  assert_true(quant.est_counts[1] < 0.01)
}

///|
test "kallisto em splits ambiguous reads" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // 仅共享区 read x4: 等长 t1/t2 -> 各分一半
  let reads : Array[String] = []
  for i = 0; i < 4; i = i + 1 {
    reads.push("ATATGCGCGCTATAT")
  }
  let quant = @src.kallisto_quantify(index, reads)
  assert_eq(quant.n_pseudoaligned, 4)
  assert_true(quant.est_counts[0] > 1.99 && quant.est_counts[0] < 2.01)
  assert_true(quant.est_counts[1] > 1.99 && quant.est_counts[1] < 2.01)
  // t1/t2 等有效长度 -> TPM 相等
  let tpm_diff = quant.tpm[0] - quant.tpm[1]
  let ad = if tpm_diff >= 0.0 { tpm_diff } else { -tpm_diff }
  assert_true(ad < 1.0)
}

///|
test "kallisto em unique context biases assignment" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  // t1 独特 read x8 + 共享 read x4: 共享应多数流向 t1
  let reads : Array[String] = []
  for i = 0; i < 8; i = i + 1 {
    reads.push("AAAAACCCCCGGGGGTTTTTATAT")
  }
  for i = 0; i < 4; i = i + 1 {
    reads.push("ATATGCGCGCTATAT")
  }
  let quant = @src.kallisto_quantify(index, reads)
  assert_true(quant.est_counts[0] > 9.9)
  assert_true(quant.est_counts[1] < 2.1)
  // 总计数守恒
  let mut total = 0.0
  for c in quant.est_counts {
    total = total + c
  }
  assert_true(total > 11.99 && total < 12.01)
}

///|
test "kallisto tpm sums to one million" {
  let (names, seqs) = kallisto_test_txs()
  let index = @src.kallisto_build_index(names, seqs, k=9)
  let reads : Array[String] = [
    "AAAAACCCCCGGGGGTTTTTATAT", "ATATGCGCGCTATAT", "ACACACGTGTGTATATATCGCGCG", "ACACACGTGTGTATATATCGCGCG",
    "GGTTCCAAGGTTCCAAGGATATGC",
  ]
  let quant = @src.kallisto_quantify(index, reads)
  let mut sum = 0.0
  for t in quant.tpm {
    sum = sum + t
  }
  assert_true(sum > 999999.0 && sum < 1000001.0)
}

///|
test "kallisto run one step entry" {
  let (names, seqs) = kallisto_test_txs()
  let quant = @src.kallisto_run(names, seqs, ["AAAAACCCCCGGGGGTTTTTATAT"], k=9)
  assert_eq(quant.n_reads, 1)
  assert_eq(quant.n_pseudoaligned, 1)
  assert_true(quant.est_counts[0] > 0.99)
}

///|
test "kallisto quant sf format and tximport roundtrip" {
  let (names, seqs) = kallisto_test_txs()
  let reads : Array[String] = [
    "AAAAACCCCCGGGGGTTTTTATAT", "AAAAACCCCCGGGGGTTTTTATAT", "ATATGCGCGCTATAT", "ACACACGTGTGTATATATCGCGCG",
  ]
  let quant = @src.kallisto_run(names, seqs, reads, k=9)
  let sf = quant.to_quant_sf()
  let lines = sf.split("\n").to_array()
  assert_eq(lines[0].to_owned(), "Name\tLength\tEffectiveLength\tTPM\tNumReads")
  // 与 tximport.parse_quant_sf 闭环
  let (tx_names, abundance, counts, lengths) = @src.parse_quant_sf(sf)
  assert_eq(tx_names.length(), 3)
  assert_eq(tx_names[0], "t1")
  assert_eq(lengths[0], 53.0)
  assert_true(counts[0] > 2.0 && counts[0] < 3.0)
  assert_true(abundance[0] > 0.0)
  assert_true(abundance.length() == counts.length())
}

///|
test "kallisto summary string" {
  let (names, seqs) = kallisto_test_txs()
  let quant = @src.kallisto_run(names, seqs, ["ACACACGTGTGTATATATCGCGCG"], k=9)
  let s = quant.summary()
  assert_true(s.contains("kallisto-lite"))
  assert_true(s.contains("pseudoaligned"))
}
