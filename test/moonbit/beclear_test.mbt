///|
/// Tests for BEclear — batch effect detection and correction.

///|
test "beclear_matrix construction" {
  let m = @src.beclear_sample_matrix()
  assert_eq(m.n_rows(), 3)
  assert_eq(m.n_cols(), 4)
  assert_eq(m.row_names[0], "gene1")
  assert_eq(m.col_names[3], "s4")
}

///|
test "beclear_detect_batch_effect finds effect" {
  let m = @src.beclear_sample_matrix()
  let summaries = @src.detect_batch_effect(m, threshold=0.05)
  assert_eq(summaries.length(), 3)
  let mut n_affected = 0
  for s in summaries {
    if s.affected {
      n_affected += 1
    }
  }
  assert_true(n_affected > 0)
}

///|
test "beclear_detect_batch_effect no effect with high threshold" {
  let m = @src.beclear_sample_matrix()
  let summaries = @src.detect_batch_effect(m, threshold=100.0)
  for s in summaries {
    assert_true(not(s.affected))
  }
}

///|
test "beclear_affected_blocks non-empty" {
  let m = @src.beclear_sample_matrix()
  let blocks = @src.beclear_affected_blocks(m, threshold=0.05)
  assert_true(blocks.length() > 0)
}

///|
test "beclear_corrects batch effect" {
  let m = @src.beclear_sample_matrix()
  let corrected = @src.beclear(m, threshold=0.05)
  assert_eq(corrected.length(), 3)
  assert_eq(corrected[0].length(), 4)
  let score_before = @src.beclear_score(m, m.data)
  let score_after = @src.beclear_score(m, corrected)
  assert_true(score_after <= score_before)
}

///|
test "beclear_score zero for no effect" {
  let row_names = ["g1"]
  let col_names = ["s1", "s2"]
  let data = [[1.0, 1.0]]
  let batch = {
    let b : Map[String, Int] = Map([])
    b.set("s1", 1)
    b.set("s2", 1)
    b
  }
  let m = @src.BEclearMatrix::new(
    row_names~,
    col_names~,
    data~,
    batch~,
  )
  let score = @src.beclear_score(m, data)
  assert_true(score < 1.0e-9)
}

///|
test "beclear_corrected preserves dimensions" {
  let m = @src.beclear_sample_matrix()
  let corrected = @src.beclear(m, threshold=0.05)
  assert_eq(corrected.length(), m.n_rows())
  for row in corrected {
    assert_eq(row.length(), m.n_cols())
  }
}

///|
test "beclear_detect_batch_effect batch medians" {
  let m = @src.beclear_sample_matrix()
  let summaries = @src.detect_batch_effect(m, threshold=0.0)
  let s0 = summaries[0]
  assert_true(s0.batch_medians.contains(1))
  assert_true(s0.batch_medians.contains(2))
}

///|
test "beclear_idempotent on clean matrix" {
  let row_names = ["g1", "g2"]
  let col_names = ["s1", "s2", "s3"]
  let data = [
    [5.0, 5.1, 4.9],
    [3.0, 3.1, 2.9],
  ]
  let batch = {
    let b : Map[String, Int] = Map([])
    b.set("s1", 1)
    b.set("s2", 1)
    b.set("s3", 1)
    b
  }
  let m = @src.BEclearMatrix::new(
    row_names~,
    col_names~,
    data~,
    batch~,
  )
  let corrected = @src.beclear(m, threshold=0.05)
  for i in 0..<2 {
    for j in 0..<3 {
      assert_true((corrected[i][j] - data[i][j]).abs() < 1.0e-9)
    }
  }
}