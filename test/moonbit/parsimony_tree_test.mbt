///|
/// Tests for maximum parsimony tree construction (parsimony_tree.mbt),
/// ported from Biopython ``Bio.Phylo.TreeConstruction``.  All expected
/// scores, neighbor orderings, branch lengths and final trees are frozen
/// from Biopython 1.88 (see /tmp/gen_parsimony_expected.py).

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn pt_leaf(name : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(bl))
}

fn pt_leaf0(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

fn pt_internal(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(branch_length=Some(bl), clades=children)
}

fn pt_internal0(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(clades=children)
}

fn pt_tree(root : @src.Clade, rooted : Bool) -> @src.Tree {
  @src.Tree::new(root, rooted=rooted)
}

fn pt_map(pairs : Array[(String, String)]) -> Map[String, String] {
  let m : Map[String, String] = Map([])
  for p in pairs {
    m[p.0] = p.1
  }
  m
}

///|
/// Render the numeric token in ``cs[lo..hi)`` with exactly 6 decimal digits,
/// rounding from the 7th. Both MoonBit's raw Double Show (``0.1``, ``2``)
/// and Biopython's shortest round-trip repr normalize to the same fixed form.
fn pt_emit_fixed6(
  cs : Array[Char],
  lo : Int,
  hi : Int,
  sb : StringBuilder,
) -> Unit {
  let mut p = lo
  if p < hi && cs[p] == '-' {
    sb.write_char('-')
    p = p + 1
  }
  let int_digits : Array[Char] = Array(capacity=0)
  let frac_digits : Array[Char] = Array(capacity=0)
  let mut in_frac = false
  while p < hi {
    let c = cs[p]
    if c == '.' {
      in_frac = true
    } else if c >= '0' && c <= '9' {
      if in_frac {
        frac_digits.push(c)
      } else {
        int_digits.push(c)
      }
    }
    p = p + 1
  }
  if int_digits.length() == 0 {
    int_digits.push('0')
  }
  let out_frac : Array[Char] = Array(capacity=0)
  let mut q = 0
  while q < 6 {
    if q < frac_digits.length() {
      out_frac.push(frac_digits[q])
    } else {
      out_frac.push('0')
    }
    q = q + 1
  }
  let mut carry = frac_digits.length() > 6 &&
    frac_digits[6].to_int() >= '5'.to_int()
  let mut k = 5
  while carry && k >= 0 {
    if out_frac[k].to_int() < '9'.to_int() {
      out_frac[k] = (out_frac[k].to_int() + 1).unsafe_to_char()
      carry = false
    } else {
      out_frac[k] = '0'
      k = k - 1
    }
  }
  if carry {
    let mut r = int_digits.length() - 1
    let mut int_carry = true
    while int_carry && r >= 0 {
      if int_digits[r].to_int() < '9'.to_int() {
        int_digits[r] = (int_digits[r].to_int() + 1).unsafe_to_char()
        int_carry = false
      } else {
        int_digits[r] = '0'
        r = r - 1
      }
    }
    if int_carry {
      sb.write_char('1')
    }
  }
  for c in int_digits {
    sb.write_char(c)
  }
  sb.write_char('.')
  for c in out_frac {
    sb.write_char(c)
  }
}

///|
/// Normalize every ``:number`` branch length in a Newick string to 6 fixed
/// decimal digits, so MoonBit raw-Double output can be compared against the
/// 6-decimal strings frozen from Biopython.
fn pt_norm6(s : String) -> String {
  let sb = StringBuilder::new()
  let cs = s.to_array()
  let mut i = 0
  while i < cs.length() {
    let c = cs[i]
    if c == ':' {
      sb.write_char(':')
      let mut j = i + 1
      while j < cs.length() {
        let t = cs[j]
        if (t >= '0' && t <= '9') || t == '.' || t == '-' {
          j = j + 1
        } else {
          break
        }
      }
      pt_emit_fixed6(cs, i + 1, j, sb)
      i = j
    } else {
      sb.write_char(c)
      i = i + 1
    }
  }
  sb.to_string()
}

///|
/// Replay the Biopython NNI hill climb, freezing every round: start score,
/// neighbor scores in order, the accepted neighbor index (-1 = none) and the
/// accepted tree Newick. Returns the final tree.
fn pt_assert_climb(
  start : @src.Tree,
  alignment : Map[String, String],
  scorer : @src.ParsimonyScorer,
  expected : Array[(Int, Array[Int], Int, String)],
) -> @src.Tree raise Error {
  let mut work = start
  let mut round = 0
  while round < expected.length() {
    let exp = expected[round]
    let mut best_score = scorer.get_score(work, alignment)
    assert_eq(best_score, exp.0)
    let neighbors = @src.nni_get_neighbors(work)
    let scores : Array[Int] = Array(capacity=0)
    let mut accepted = -1
    let mut i = 0
    while i < neighbors.length() {
      let s = scorer.get_score(neighbors[i], alignment)
      scores.push(s)
      if s < best_score {
        best_score = s
        work = neighbors[i]
        accepted = i
      }
      i = i + 1
    }
    assert_eq(scores, exp.1)
    assert_eq(accepted, exp.2)
    if accepted >= 0 {
      assert_eq(pt_norm6(work.root.to_newick()), exp.3)
    } else {
      break
    }
    round = round + 1
  }
  work
}

// ============================================================================
// Fixture A: rooted balanced 4-taxon tree, Fitch
// ============================================================================

///|
test "pt_balanced_four_taxon_scores_and_neighbors" {
  let t = pt_tree(
    pt_internal([
      pt_internal([pt_leaf("A", 0.1), pt_leaf("B", 0.2)], 0.3),
      pt_internal([pt_leaf("C", 0.4), pt_leaf("D", 0.5)], 0.6),
    ], 0.0),
    true,
  )
  let aln = pt_map([
    ("A", "AATA"),
    ("B", "AATC"),
    ("C", "GCGG"),
    ("D", "GCGT"),
  ])
  let scorer = @src.ParsimonyScorer::new()
  assert_eq(scorer.get_score(t, aln), 6)
  let neighbors = @src.nni_get_neighbors(t)
  assert_eq(neighbors.length(), 2)
  // neighbor 1 (left_left + right_right)
  assert_eq(scorer.get_score(neighbors[0], aln), 9)
  assert_eq(
    pt_norm6(neighbors[0].root.to_newick()),
    "((A:0.100000,D:0.500000):0.300000,(C:0.400000,B:0.200000):0.600000):0.000000;",
  )
  // neighbor 2 (left_left + right_left)
  assert_eq(scorer.get_score(neighbors[1], aln), 9)
  assert_eq(
    pt_norm6(neighbors[1].root.to_newick()),
    "((A:0.100000,C:0.400000):0.300000,(B:0.200000,D:0.500000):0.600000):0.000000;",
  )
  // no improving neighbor: search returns the starting topology unchanged
  let searcher = @src.NNITreeSearcher::new(scorer)
  let result = try { searcher.search(t, aln) } catch {
    @src.ParsimonyTreeError(message) => abort(message)
  }
  assert_eq(result.is_rooted(), true)
  assert_eq(scorer.get_score(result, aln), 6)
  assert_eq(
    result.root.to_newick(lengths=false),
    "((A,B),(C,D));",
  )
}

// ============================================================================
// Fixture A2: rooted pectinate 4-taxon tree
// ============================================================================

///|
test "pt_pectinate_four_taxon_neighbors" {
  let t = pt_tree(
    pt_internal([
      pt_internal([
        pt_internal([pt_leaf("A", 0.1), pt_leaf("B", 0.15)], 0.2),
        pt_leaf("C", 0.25),
      ], 0.3),
      pt_leaf("D", 0.4),
    ], 0.0),
    true,
  )
  let aln = pt_map([
    ("A", "AATA"),
    ("B", "AATC"),
    ("C", "GCGG"),
    ("D", "GCGT"),
  ])
  let scorer = @src.ParsimonyScorer::new()
  assert_eq(scorer.get_score(t, aln), 6)
  let neighbors = @src.nni_get_neighbors(t)
  assert_eq(neighbors.length(), 2)
  assert_eq(
    pt_norm6(neighbors[0].root.to_newick()),
    "(((A:0.100000,C:0.250000):0.200000,B:0.150000):0.300000,D:0.400000):0.000000;",
  )
  assert_eq(
    pt_norm6(neighbors[1].root.to_newick()),
    "(((C:0.250000,B:0.150000):0.200000,A:0.100000):0.300000,D:0.400000):0.000000;",
  )
  assert_eq(scorer.get_score(neighbors[0], aln), 9)
  assert_eq(scorer.get_score(neighbors[1], aln), 9)
}

// ============================================================================
// Fixture G: 5-taxon Fitch, two-round hill climb with distinct lengths
// ============================================================================

///|
test "pt_fitch_two_round_climb" {
  let t = pt_tree(
    pt_internal([
      pt_leaf("S1", 0.11),
      pt_internal([
        pt_leaf("S2", 0.12),
        pt_internal([
          pt_leaf("S3", 0.13),
          pt_internal([pt_leaf("S4", 0.14), pt_leaf("S5", 0.15)], 0.25),
        ], 0.23),
      ], 0.21),
    ], 0.02),
    true,
  )
  let aln = pt_map([
    ("S1", "AAAAA"),
    ("S2", "AAGAA"),
    ("S3", "GTATT"),
    ("S4", "GTAAT"),
    ("S5", "TTTTT"),
  ])
  let scorer = @src.ParsimonyScorer::new()
  assert_eq(scorer.get_score(t, aln), 8)
  let neighbors = @src.nni_get_neighbors(t)
  assert_eq(neighbors.length(), 4)
  assert_eq(
    pt_norm6(neighbors[0].root.to_newick()),
    "(S1:0.110000,((S4:0.140000,S5:0.150000):0.250000,(S3:0.130000,S2:0.120000):0.230000):0.210000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[1].root.to_newick()),
    "(S1:0.110000,(S3:0.130000,(S2:0.120000,(S4:0.140000,S5:0.150000):0.250000):0.230000):0.210000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[2].root.to_newick()),
    "(S1:0.110000,(S2:0.120000,(S5:0.150000,(S4:0.140000,S3:0.130000):0.250000):0.230000):0.210000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[3].root.to_newick()),
    "(S1:0.110000,(S2:0.120000,(S4:0.140000,(S3:0.130000,S5:0.150000):0.250000):0.230000):0.210000):0.020000;",
  )
  let expected : Array[(Int, Array[Int], Int, String)] = [
    (
      8,
      [11, 11, 8, 7],
      3,
      "(S1:0.110000,(S2:0.120000,(S4:0.140000,(S3:0.130000,S5:0.150000):0.250000):0.230000):0.210000):0.020000;",
    ),
    (7, [10, 10, 8, 8], -1, ""),
  ]
  let final_tree = pt_assert_climb(t, aln, scorer, expected)
  assert_eq(scorer.get_score(final_tree, aln), 7)
  assert_eq(
    pt_norm6(final_tree.root.to_newick()),
    "(S1:0.110000,(S2:0.120000,(S4:0.140000,(S3:0.130000,S5:0.150000):0.250000):0.230000):0.210000):0.020000;",
  )
  // end-state neighbor scores [10, 10, 8, 8]
  let final_neighbors = @src.nni_get_neighbors(final_tree)
  assert_eq(
    final_neighbors.map(fn(x) -> Int raise @src.ParsimonyTreeError {
      scorer.get_score(x, aln)
    }),
    [10, 10, 8, 8],
  )
}

// ============================================================================
// Fixture B: 5-taxon unrooted NJ start -> midpoint root -> Fitch NNI
// ============================================================================

fn pt_biopython_msa() -> Array[(String, String)] {
  [
    ("Alpha", "AACGTGGCCACAT"),
    ("Beta", "AAGGTCGCCACAC"),
    ("Gamma", "CAGTTCGCCACAA"),
    ("Delta", "GAGATTTCCGCCT"),
    ("Epsilon", "GAGATCTCCGCCC"),
  ]
}

///|
test "pt_nj_unrooted_midpoint_and_search" {
  let pairs = pt_biopython_msa()
  let aln = pt_map(pairs)
  let calculator = try {
    @src.DistanceCalculator::new(model="identity")
  } catch {
    @src.SeqError(message) => abort(message)
  }
  let dm = try {
    calculator.get_distance(pairs)
  } catch {
    @src.SeqError(message) => abort(message)
  }
  let names : Array[String] = pairs.map(fn(kv) { kv.0 })
  let builder = @src.TreeConstructor::new(tree_method="nj")
  let nj = pt_tree(builder.build_tree(dm, names), false)
  assert_eq(nj.is_rooted(), false)
  // scorer auto-roots unrooted trees at their midpoint (Biopython behaviour)
  let scorer = @src.ParsimonyScorer::new()
  assert_eq(scorer.get_score(nj, aln), 13)

  let mid = nj.root_at_midpoint()
  assert_eq(mid.is_rooted(), true)
  assert_eq(
    mid.root.to_newick(lengths=false),
    "(((Alpha,Beta),Gamma),(Delta,Epsilon));",
  )
  // frozen midpoint branch lengths (Biopython 1.88)
  let left = mid.root.clades[0]
  let right = mid.root.clades[1]
  assert_true((left.branch_length.unwrap() - 0.07532051282051289).abs() < 0.000000000001)
  assert_true((right.branch_length.unwrap() - 0.20352564102564097).abs() < 0.000000000001)
  assert_true((left.clades[0].branch_length.unwrap() - 0.04807692307692307).abs() < 0.000000000001)
  assert_true((left.clades[0].clades[0].branch_length.unwrap() - 0.18269230769230765).abs() < 0.000000000001)
  assert_true((left.clades[0].clades[1].branch_length.unwrap() - 0.04807692307692307).abs() < 0.000000000001)
  assert_true((left.clades[1].branch_length.unwrap() - 0.14423076923076922).abs() < 0.000000000001)
  assert_true((right.clades[0].branch_length.unwrap() - 0.10256410256410259).abs() < 0.000000000001)
  assert_true((right.clades[1].branch_length.unwrap() - 0.051282051282051266).abs() < 0.000000000001)

  // midpoint neighbor scores [16, 16, 13, 13], none strictly better
  let mid_scores =
    @src.nni_get_neighbors(mid).map(fn(x) -> Int raise @src.ParsimonyTreeError {
      scorer.get_score(x, aln)
    })
  assert_eq(mid_scores, [16, 16, 13, 13])
  let searcher = @src.NNITreeSearcher::new(scorer)
  let result = try { searcher.search(nj, aln) } catch {
    @src.ParsimonyTreeError(message) => abort(message)
  }
  assert_eq(result.is_rooted(), true)
  assert_eq(scorer.get_score(result, aln), 13)
  assert_eq(
    result.root.to_newick(lengths=false),
    "(((Alpha,Beta),Gamma),(Delta,Epsilon));",
  )
}

///|
test "pt_default_upgma_start_and_constructor" {
  let pairs = pt_biopython_msa()
  let aln = pt_map(pairs)

  // No starting tree: the constructor builds Biopython's identity-UPGMA
  // start internally (equal-weight row averaging, last-equal tie wins) and no
  // NNI neighbor is strictly better, so the result is that tree unchanged.
  let scorer = @src.ParsimonyScorer::new()
  let ctor = @src.ParsimonyTreeConstructor::new(
    @src.NNITreeSearcher::new(@src.ParsimonyScorer::new()),
  )
  let built = try {
    ctor.build_tree(aln)
  } catch {
    @src.ParsimonyTreeError(message) => abort(message)
  }
  assert_eq(built.is_rooted(), true)
  // frozen Biopython 1.88 UPGMA identity tree, full branch lengths
  assert_eq(
    pt_norm6(built.root.to_newick()),
    "((Epsilon:0.076923,Delta:0.076923):0.187500,((Gamma:0.115385,Beta:0.115385):0.038462,Alpha:0.153846):0.110577):0.000000;",
  )
  assert_eq(
    built.root.to_newick(lengths=false),
    "((Epsilon,Delta),((Gamma,Beta),Alpha));",
  )
  let eps_delta = built.root.clades[0]
  let rest = built.root.clades[1]
  assert_true((eps_delta.branch_length.unwrap() - 0.18749999999999994).abs() < 0.000000000001)
  assert_true((eps_delta.clades[0].branch_length.unwrap() - 0.07692307692307693).abs() < 0.000000000001)
  assert_true((eps_delta.clades[1].branch_length.unwrap() - 0.07692307692307693).abs() < 0.000000000001)
  assert_true((rest.branch_length.unwrap() - 0.11057692307692304).abs() < 0.000000000001)
  assert_true((rest.clades[0].branch_length.unwrap() - 0.038461538461538464).abs() < 0.000000000001)
  assert_true((rest.clades[0].clades[0].branch_length.unwrap() - 0.11538461538461536).abs() < 0.000000000001)
  assert_true((rest.clades[0].clades[1].branch_length.unwrap() - 0.11538461538461536).abs() < 0.000000000001)
  assert_true((rest.clades[1].branch_length.unwrap() - 0.15384615384615383).abs() < 0.000000000001)
  assert_eq(scorer.get_score(built, aln), 13)
}

// ============================================================================
// Fixture C: 6-taxon Sankoff matrix, five-round climb (52 -> 32)
// ============================================================================

///|
test "pt_sankoff_five_round_climb" {
  let start = pt_tree(
    pt_internal([
      pt_internal([
        pt_internal([
          pt_internal([pt_leaf("S1", 0.11), pt_leaf("S3", 0.13)], 0.31),
          pt_leaf("S4", 0.14),
        ], 0.32),
        pt_leaf("S5", 0.15),
      ], 0.33),
      pt_internal([pt_leaf("S6", 0.16), pt_leaf("S2", 0.12)], 0.34),
    ], 0.02),
    true,
  )
  let aln = pt_map([
    ("S1", "AAAAA"),
    ("S2", "AAGAA"),
    ("S3", "GTATT"),
    ("S4", "GTAAT"),
    ("S5", "TTTTT"),
    ("S6", "TTTGT"),
  ])
  // transitions A<->G, C<->T cost 1; transversions cost 5
  let matrix = @src.ParsimonyMatrix::new(
    ["A", "C", "G", "T"],
    [
      [0.0, 5.0, 1.0, 5.0],
      [5.0, 0.0, 5.0, 1.0],
      [1.0, 5.0, 0.0, 5.0],
      [5.0, 1.0, 5.0, 0.0],
    ],
  )
  let fitch = @src.ParsimonyScorer::new()
  let scorer = @src.ParsimonyScorer::new(matrix=Some(matrix))
  assert_eq(fitch.get_score(start, aln), 12)
  assert_eq(scorer.get_score(start, aln), 52)

  let neighbors = @src.nni_get_neighbors(start)
  assert_eq(neighbors.length(), 6)
  let nscores = neighbors.map(fn(x) -> Int raise @src.ParsimonyTreeError {
    scorer.get_score(x, aln)
  })
  assert_eq(nscores, [44, 52, 54, 54, 52, 52])
  assert_eq(
    pt_norm6(neighbors[0].root.to_newick()),
    "((((S1:0.110000,S3:0.130000):0.310000,S4:0.140000):0.320000,S2:0.120000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[1].root.to_newick()),
    "((((S1:0.110000,S3:0.130000):0.310000,S4:0.140000):0.320000,S6:0.160000):0.330000,(S5:0.150000,S2:0.120000):0.340000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[2].root.to_newick()),
    "((((S1:0.110000,S3:0.130000):0.310000,S5:0.150000):0.320000,S4:0.140000):0.330000,(S6:0.160000,S2:0.120000):0.340000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[3].root.to_newick()),
    "(((S5:0.150000,S4:0.140000):0.320000,(S1:0.110000,S3:0.130000):0.310000):0.330000,(S6:0.160000,S2:0.120000):0.340000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[4].root.to_newick()),
    "((((S1:0.110000,S4:0.140000):0.310000,S3:0.130000):0.320000,S5:0.150000):0.330000,(S6:0.160000,S2:0.120000):0.340000):0.020000;",
  )
  assert_eq(
    pt_norm6(neighbors[5].root.to_newick()),
    "((((S4:0.140000,S3:0.130000):0.310000,S1:0.110000):0.320000,S5:0.150000):0.330000,(S6:0.160000,S2:0.120000):0.340000):0.020000;",
  )

  let expected : Array[(Int, Array[Int], Int, String)] = [
    (
      52,
      [44, 52, 54, 54, 52, 52],
      0,
      "((((S1:0.110000,S3:0.130000):0.310000,S4:0.140000):0.320000,S2:0.120000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
    ),
    (
      44,
      [52, 52, 44, 44, 44, 43],
      5,
      "((((S4:0.140000,S3:0.130000):0.310000,S1:0.110000):0.320000,S2:0.120000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
    ),
    (
      43,
      [52, 52, 43, 33, 44, 44],
      3,
      "(((S2:0.120000,S1:0.110000):0.320000,(S4:0.140000,S3:0.130000):0.310000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
    ),
    (
      33,
      [42, 42, 43, 43, 32, 33],
      4,
      "((S3:0.130000,(S4:0.140000,(S2:0.120000,S1:0.110000):0.320000):0.310000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
    ),
    (32, [38, 42, 33, 33, 43, 43], -1, ""),
  ]
  let final_tree = pt_assert_climb(start, aln, scorer, expected)
  assert_eq(
    pt_norm6(final_tree.root.to_newick()),
    "((S3:0.130000,(S4:0.140000,(S2:0.120000,S1:0.110000):0.320000):0.310000):0.330000,(S6:0.160000,S5:0.150000):0.340000):0.020000;",
  )
  assert_eq(scorer.get_score(final_tree, aln), 32)
  assert_eq(fitch.get_score(final_tree, aln), 8)

  // the public searcher reaches the same optimum
  let result = try {
    @src.NNITreeSearcher::new(scorer).search(start, aln)
  } catch {
    @src.ParsimonyTreeError(message) => abort(message)
  }
  assert_eq(scorer.get_score(result, aln), 32)
}

// ============================================================================
// Fixture E: root_at_midpoint
// ============================================================================

///|
test "pt_root_at_midpoint_unrooted_trees" {
  // 3-tip trifurcating unrooted tree with unequal branches
  let star = pt_tree(
    pt_internal0([
      pt_leaf("A", 1.0),
      pt_leaf("B", 2.0),
      pt_leaf("C", 4.0),
    ]),
    false,
  )
  let m1 = star.root_at_midpoint()
  assert_eq(pt_norm6(m1.root.to_newick()), "((B:2.000000,A:1.000000):1.000000,C:3.000000);")

  // 4-tip case; diameter is tip D to any of A/B/C, midpoint on the D edge
  let t2 = pt_tree(
    pt_internal0([
      pt_internal([pt_leaf("A", 0.1), pt_leaf("B", 0.1)], 0.2),
      pt_leaf("C", 0.3),
      pt_leaf("D", 1.5),
    ]),
    false,
  )
  let m2 = t2.root_at_midpoint()
  assert_eq(
    pt_norm6(m2.root.to_newick()),
    "(((A:0.100000,B:0.100000):0.200000,C:0.300000):0.600000,D:0.900000);",
  )
  assert_eq(m2.is_rooted(), true)
  // midpoint rooting does not mutate the input tree
  assert_eq(star.is_rooted(), false)
  assert_eq(
    star.root.to_newick(lengths=false),
    "(A,B,C);",
  )
}

// ============================================================================
// Fixture F: error conditions
// ============================================================================

///|
test "pt_rejects_non_bifurcating_and_mismatched_names" {
  let scorer = @src.ParsimonyScorer::new()
  // internal trifurcation (root may be trifurcating, deeper nodes may not)
  let poly = pt_tree(
    pt_internal0([
      pt_internal0([pt_leaf0("A"), pt_leaf0("B")]),
      pt_internal0([
        pt_leaf0("C"),
        pt_internal0([pt_leaf0("D"), pt_leaf0("E"), pt_leaf0("F")]),
      ]),
    ]),
    true,
  )
  let aln1 = pt_map([
    ("A", "X"),
    ("B", "Y"),
    ("C", "Z"),
    ("D", "W"),
    ("E", "V"),
    ("F", "U"),
  ])
  let msg1 = try {
    ignore(scorer.get_score(poly, aln1))
    "NO_ERROR"
  } catch {
    @src.ParsimonyTreeError(message) => message
  }
  assert_eq(msg1, "The tree provided should be bifurcating.")

  // taxon name sets differ
  let t = pt_tree(
    pt_internal0([
      pt_internal0([pt_leaf0("A"), pt_leaf0("B")]),
      pt_internal0([pt_leaf0("C"), pt_leaf0("D")]),
    ]),
    true,
  )
  let aln2 = pt_map([
    ("A", "X"),
    ("B", "Y"),
    ("C", "Z"),
    ("E", "W"),
  ])
  let msg2 = try {
    ignore(scorer.get_score(t, aln2))
    "NO_ERROR"
  } catch {
    @src.ParsimonyTreeError(message) => message
  }
  assert_eq(
    msg2,
    "Taxon names of the input tree should be the same with the alignment.",
  )
}

///|
test "pt_scorer_does_not_mutate_input_tree" {
  let t = pt_tree(
    pt_internal([
      pt_internal([pt_leaf("A", 0.1), pt_leaf("B", 0.2)], 0.3),
      pt_internal([pt_leaf("C", 0.4), pt_leaf("D", 0.5)], 0.6),
    ], 0.0),
    true,
  )
  let before = t.root.to_newick()
  let aln = pt_map([
    ("A", "AATA"),
    ("B", "AATC"),
    ("C", "GCGG"),
    ("D", "GCGT"),
  ])
  let scorer = @src.ParsimonyScorer::new()
  assert_eq(scorer.get_score(t, aln), 6)
  assert_eq(scorer.get_score(t, aln), 6)
  assert_eq(t.root.to_newick(), before)
}
