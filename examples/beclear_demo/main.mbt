///|
fn main {
  println("=== BEclear Demo ===")

  let m = @src.beclear_sample_matrix()
  println("\n1. Input matrix: " + m.n_rows().to_string() + " genes x " + m.n_cols().to_string() + " samples")
  for i in 0..<m.n_rows() {
    let mut row = ""
    for j in 0..<m.n_cols() {
      row = row + m.data[i][j].to_string() + "  "
    }
    println("   " + m.row_names[i] + ": " + row)
  }

  println("\n2. Batch assignment")
  for j in 0..<m.n_cols() {
    let b = m.batch.get(m.col_names[j])
    match b {
      Some(batch_id) => println("   " + m.col_names[j] + " -> batch " + batch_id.to_string())
      None => println("   " + m.col_names[j] + " -> unassigned")
    }
  }

  println("\n3. Batch effect detection (threshold=0.05)")
  let summaries = @src.detect_batch_effect(m, threshold=0.05)
  for s in summaries {
    println(
      "   " + s.row_name +
      ": max_dev=" + s.max_deviation.to_string() +
      ", affected=" + s.affected.to_string(),
    )
  }

  println("\n4. Affected blocks")
  let blocks = @src.beclear_affected_blocks(m, threshold=0.05)
  for b in blocks {
    println("   row " + b.0.to_string() + ", batch " + b.1.to_string())
  }

  println("\n5. Corrected matrix")
  let corrected = @src.beclear(m, threshold=0.05)
  for i in 0..<corrected.length() {
    let mut row = ""
    for j in 0..<corrected[i].length() {
      row = row + corrected[i][j].to_string() + "  "
    }
    println("   " + m.row_names[i] + ": " + row)
  }

  println("\n6. BEscore before vs after")
  let score_before = @src.beclear_score(m, m.data)
  let score_after = @src.beclear_score(m, corrected)
  println("   Before: " + score_before.to_string())
  println("   After:  " + score_after.to_string())
}