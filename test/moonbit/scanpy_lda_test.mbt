///|
/// Tests for the LDA topic model (collapsed Gibbs).

///|
test "scanpy_lda_fit" {
  // two clear topics: docs 0/1 heavy on genes 0-1, docs 2/3 heavy on genes 2-3
  let x = [
    [20.0, 20.0, 1.0, 1.0],
    [20.0, 20.0, 1.0, 1.0],
    [1.0, 1.0, 20.0, 20.0],
    [1.0, 1.0, 20.0, 20.0],
  ]
  let res = @src.scanpy_lda_fit(x, n_topics=2, n_iter=50, seed=7)
  assert_eq(res.doc_topic.length(), 4)
  assert_eq(res.doc_topic[0].length(), 2)
  assert_eq(res.topic_term.length(), 2)
  assert_eq(res.topic_term[0].length(), 4)
  // rows are probability distributions
  for d in 0..<4 {
    let s = res.doc_topic[d][0] + res.doc_topic[d][1]
    assert_true(@moonbit.assert_num_close(s, 1.0, 1.0e-9))
  }
  for t in 0..<2 {
    let mut s = 0.0
    for g in 0..<4 {
      s = s + res.topic_term[t][g]
    }
    assert_true(@moonbit.assert_num_close(s, 1.0, 1.0e-9))
  }
  // docs 0 and 1 share one dominant topic, docs 2 and 3 the other
  let dom0 = if res.doc_topic[0][0] > res.doc_topic[0][1] { 0 } else { 1 }
  let dom2 = if res.doc_topic[2][0] > res.doc_topic[2][1] { 0 } else { 1 }
  assert_true(dom0 != dom2)
  assert_true(res.doc_topic[0][dom0] > 0.8)
  assert_true(res.doc_topic[2][dom2] > 0.8)
  // dominant topic of docs 0/1 is concentrated on genes 0-1
  let mass_01 = res.topic_term[dom0][0] + res.topic_term[dom0][1]
  let mass_23 = res.topic_term[dom2][2] + res.topic_term[dom2][3]
  assert_true(mass_01 > 0.8)
  assert_true(mass_23 > 0.8)
}

///|
test "scanpy_lda_topics_adata" {
  let x = [[5.0, 5.0, 0.0], [5.0, 5.0, 0.0], [0.0, 0.0, 5.0]]
  let adata = @src.AnnData::new(x, ["d0", "d1", "d2"], ["A", "B", "C"])
  let res = @src.scanpy_lda_topics(adata, n_topics=2, n_iter=30, seed=3)
  assert_eq(res.n_topics, 2)
  // results stored on the AnnData
  match adata.obsm.get("X_lda") {
    Some(m) => {
      assert_eq(m.length(), 3)
      assert_eq(m[0].length(), 2)
    }
    None => assert_true(false)
  }
  match adata.varm.get("lda_topics") {
    Some(m) => {
      assert_eq(m.length(), 3) // n_vars
      assert_eq(m[0].length(), 2) // n_topics
    }
    None => assert_true(false)
  }
}
