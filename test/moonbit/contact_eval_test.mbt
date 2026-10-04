// ===========================================================================
// contact_eval_test.mbt — Tests for DCA contact map vs 3D structure comparison
// ===========================================================================

///|
/// Build a simple test PDB with 10 residues in a compact arrangement.
/// Residues 0-4 form a helix-like arrangement, 5-9 fold back.
/// Contacts: (0,8), (1,7), (2,6), (3,5) are close in space.
fn make_test_pdb() -> String {
  "MODEL     1\n" +
  "ATOM      1  CA  ALA A   1       1.000   1.000   0.000  1.00 20.00           C\n" +
  "ATOM      2  CA  ALA A   2       3.000   1.000   0.000  1.00 20.00           C\n" +
  "ATOM      3  CA  ALA A   3       5.000   1.000   0.000  1.00 20.00           C\n" +
  "ATOM      4  CA  ALA A   4       7.000   1.000   0.000  1.00 20.00           C\n" +
  "ATOM      5  CA  ALA A   5       9.000   1.000   0.000  1.00 20.00           C\n" +
  "ATOM      6  CA  ALA A   6       9.000   4.000   0.000  1.00 20.00           C\n" +
  "ATOM      7  CA  ALA A   7       7.000   4.000   0.000  1.00 20.00           C\n" +
  "ATOM      8  CA  ALA A   8       5.000   4.000   0.000  1.00 20.00           C\n" +
  "ATOM      9  CA  ALA A   9       3.000   4.000   0.000  1.00 20.00           C\n" +
  "ATOM     10  CA  ALA A  10       1.000   4.000   0.000  1.00 20.00           C\n" +
  "ENDMDL\n" +
  "END\n"
}

///|
/// Get the test chain.
fn get_test_chain() -> @src.Chain {
  let pdb = make_test_pdb()
  let s = @src.parse_pdb(pdb)
  let chains = s.get_chains()
  chains[0]
}

test "extract_true_contacts basic" {
  let chain = get_test_chain()
  let contacts = @src.extract_true_contacts(chain, threshold=5.0, min_seq_dist=3)
  // Should have contacts between the two strands
  assert_true(contacts.length() > 0)
  // All pairs should have i < j
  for c in contacts {
    let (i, j) = c
    assert_true(i < j)
  }
}

test "extract_true_contacts threshold" {
  let chain = get_test_chain()
  // With very large threshold, many contacts
  let many = @src.extract_true_contacts(chain, threshold=100.0, min_seq_dist=1)
  // With small threshold, few contacts
  let few = @src.extract_true_contacts(chain, threshold=2.0, min_seq_dist=1)
  assert_true(many.length() > few.length())
}

test "extract_true_contacts min_seq_dist" {
  let chain = get_test_chain()
  let close = @src.extract_true_contacts(chain, threshold=100.0, min_seq_dist=1)
  let far = @src.extract_true_contacts(chain, threshold=100.0, min_seq_dist=5)
  // Far should have fewer contacts (excludes |i-j| < 5)
  assert_true(far.length() <= close.length())
}

test "evaluate_contact_prediction perfect" {
  let true_contacts = [(0, 5), (1, 6), (2, 7), (3, 8)]
  let predicted = [(0, 5), (1, 6), (2, 7), (3, 8)]
  let eval = @src.evaluate_contact_prediction(predicted, true_contacts)
  assert_true((eval.precision - 1.0).abs() < 1.0e-10)
  assert_true((eval.recall - 1.0).abs() < 1.0e-10)
  assert_true((eval.f1 - 1.0).abs() < 1.0e-10)
  assert_eq(eval.tp, 4)
  assert_eq(eval.fp, 0)
  assert_eq(eval.fn_count, 0)
}

test "evaluate_contact_prediction no overlap" {
  let true_contacts = [(0, 5), (1, 6), (2, 7), (3, 8)]
  let predicted = [(0, 1), (2, 3), (4, 5), (6, 7)]
  let eval = @src.evaluate_contact_prediction(predicted, true_contacts)
  assert_true(eval.precision.abs() < 1.0e-10)
  assert_true(eval.recall.abs() < 1.0e-10)
  assert_eq(eval.tp, 0)
  assert_eq(eval.fp, 4)
  assert_eq(eval.fn_count, 4)
}

test "evaluate_contact_prediction partial" {
  let true_contacts = [(0, 5), (1, 6), (2, 7), (3, 8)]
  let predicted = [(0, 5), (1, 6), (0, 1), (2, 3)]
  let eval = @src.evaluate_contact_prediction(predicted, true_contacts)
  // 2 TP, 2 FP, 2 FN
  assert_eq(eval.tp, 2)
  assert_eq(eval.fp, 2)
  assert_eq(eval.fn_count, 2)
  assert_true((eval.precision - 0.5).abs() < 1.0e-10)
  assert_true((eval.recall - 0.5).abs() < 1.0e-10)
}

test "evaluate_contact_prediction empty predicted" {
  let true_contacts = [(0, 5), (1, 6)]
  let predicted : Array[(Int, Int)] = []
  let eval = @src.evaluate_contact_prediction(predicted, true_contacts)
  assert_true(eval.precision.abs() < 1.0e-10)
  assert_true(eval.recall.abs() < 1.0e-10)
  assert_eq(eval.tp, 0)
}

test "evaluate_contact_prediction empty true" {
  let true_contacts : Array[(Int, Int)] = []
  let predicted = [(0, 5), (1, 6)]
  let eval = @src.evaluate_contact_prediction(predicted, true_contacts)
  assert_true(eval.precision.abs() < 1.0e-10)
  assert_eq(eval.fp, 2)
}

test "residue_ca_distance" {
  let chain = get_test_chain()
  // Residues 0 and 1: (1,1,0) to (3,1,0) = 2.0
  let d = @src.residue_ca_distance(chain, 0, 1)
  assert_true((d - 2.0).abs() < 1.0e-6)
  // Residues 0 and 8: (1,1,0) to (3,4,0) = sqrt(4+9) = sqrt(13) ≈ 3.606
  let d2 = @src.residue_ca_distance(chain, 0, 8)
  assert_true((d2 - 3.605551).abs() < 1.0e-4)
}

test "evaluate_with_distances" {
  let chain = get_test_chain()
  let true_contacts = @src.extract_true_contacts(chain, threshold=5.0, min_seq_dist=3)
  let predicted = true_contacts
  let eval = @src.evaluate_with_distances(predicted, true_contacts, chain)
  // Perfect prediction
  assert_true((eval.precision - 1.0).abs() < 1.0e-10)
  // Mean TP distance should be positive
  assert_true(eval.mean_tp_dist > 0.0)
  // No false positives
  assert_true(eval.mean_fp_dist.abs() < 1.0e-10)
}

test "dca_structure_comparison" {
  let chain = get_test_chain()
  // Create a simple score matrix that favors true contacts
  let n = 10
  let scores = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    scores.push(Array::make(n, 0.0))
  }
  // Give high scores to true contact pairs
  let true_pairs = [(0, 8), (1, 7), (2, 6), (3, 5)]
  for p in true_pairs {
    let (i, j) = p
    scores[i][j] = 10.0
    scores[j][i] = 10.0
  }
  // Give low scores to non-contacts
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j && scores[i][j] == 0.0 {
        scores[i][j] = 0.1
      }
    }
  }
  let eval = @src.dca_structure_comparison(scores, chain, 4, threshold=5.0, min_seq_dist=3)
  // Should have high precision (all 4 predicted are true contacts)
  assert_true(eval.precision > 0.5)
  assert_eq(eval.num_predicted, 4)
}

test "top_k_precision_curve" {
  let true_contacts = [(0, 5), (1, 6), (2, 7), (3, 8)]
  // Score matrix: true contacts have score 10, others 1
  let n = 10
  let scores = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    scores.push(Array::make(n, 0.0))
  }
  for p in true_contacts {
    let (i, j) = p
    scores[i][j] = 10.0
    scores[j][i] = 10.0
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j && scores[i][j] == 0.0 {
        scores[i][j] = 1.0
      }
    }
  }
  let curve = @src.top_k_precision_curve(scores, true_contacts, [1, 2, 4], 3)
  assert_eq(curve.length(), 3)
  // Top 4 should have precision 1.0 (all true contacts)
  let (_, prec4) = curve[2]
  assert_true((prec4 - 1.0).abs() < 1.0e-10)
}

test "contact_map_jaccard identical" {
  let m = [[true, false], [false, true]]
  let j = @src.contact_map_jaccard(m, m)
  assert_true((j - 1.0).abs() < 1.0e-10)
}

test "contact_map_jaccard disjoint" {
  let m1 = [[true, false], [false, false]]
  let m2 = [[false, false], [false, true]]
  let j = @src.contact_map_jaccard(m1, m2)
  assert_true(j.abs() < 1.0e-10)
}

test "contact_map_jaccard partial" {
  let m1 = [[true, true], [true, false]]
  let m2 = [[true, false], [false, true]]
  // Intersection: (0,0) = 1
  // Union: (0,0), (0,1), (1,0), (1,1) = 4
  let j = @src.contact_map_jaccard(m1, m2)
  assert_true((j - 0.25).abs() < 1.0e-10)
}

test "contact_map_overlap" {
  let m1 = [[true, true], [true, false]]
  let m2 = [[true, false], [false, true]]
  // Intersection: 1
  // min(|m1|, |m2|) = min(3, 2) = 2
  let o = @src.contact_map_overlap(m1, m2)
  assert_true((o - 0.5).abs() < 1.0e-10)
}

test "coupling_distance_correlation" {
  // Perfect anti-correlation: high score → short distance
  let scores = [[0.0, 10.0, 1.0], [10.0, 0.0, 5.0], [1.0, 5.0, 0.0]]
  let dist = [[0.0, 1.0, 5.0], [1.0, 0.0, 2.0], [5.0, 2.0, 0.0]]
  let (pearson, spearman) = @src.coupling_distance_correlation(scores, dist, 1)
  // Score and 1/distance should be positively correlated
  assert_true(pearson > 0.0)
}

test "stratified_evaluation" {
  let chain = get_test_chain()
  let n = 10
  let scores = Array::make(n, Array::make(n, 1.0))
  // Give some contacts high scores
  scores[0][8] = 10.0
  scores[8][0] = 10.0
  scores[1][7] = 10.0
  scores[7][1] = 10.0
  let result = @src.stratified_evaluation(scores, chain, [3, 5], 4, threshold=5.0)
  assert_eq(result.length(), 2)
  let (sep1, _) = result[0]
  assert_eq(sep1, 3)
  let (sep2, _) = result[1]
  assert_eq(sep2, 5)
}

test "contact_eval_report" {
  let chain = get_test_chain()
  let n = 10
  let scores = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    scores.push(Array::make(n, 0.0))
  }
  scores[0][8] = 10.0
  scores[8][0] = 10.0
  scores[1][7] = 10.0
  scores[7][1] = 10.0
  scores[2][6] = 10.0
  scores[6][2] = 10.0
  scores[3][5] = 10.0
  scores[5][3] = 10.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j && scores[i][j] == 0.0 {
        scores[i][j] = 0.1
      }
    }
  }
  let report = @src.contact_eval_report(scores, chain, 4, threshold=5.0, min_seq_dist=3)
  assert_true(report.overall.precision >= 0.0)
  assert_true(report.overall.precision <= 1.0)
  assert_eq(report.stratified.length(), 3)
  assert_true(report.pearson_r >= -1.0)
  assert_true(report.pearson_r <= 1.0)
  assert_true(report.spearman_r >= -1.0)
  assert_true(report.spearman_r <= 1.0)
  assert_true(report.jaccard >= 0.0)
  assert_true(report.jaccard <= 1.0)
}

test "rank_dca_pairs" {
  let scores = [[0.0, 5.0, 1.0], [5.0, 0.0, 3.0], [1.0, 3.0, 0.0]]
  let ranked = @src.rank_dca_pairs(scores, 2, 1)
  assert_eq(ranked.length(), 2)
  // Top pair should be (0, 1) with score 5.0
  let (key, score) = ranked[0]
  let (i, j) = (key / 100000, key % 100000)
  assert_eq(i, 0)
  assert_eq(j, 1)
  assert_true((score - 5.0).abs() < 1.0e-10)
}
