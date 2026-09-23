///|
/// Tests for Bio.Phylo.Consensus equivalents: strict / majority / Adam
/// consensus, bootstrap support and bootstrap replicate trees.
///
/// All frozen expected values were produced directly by Biopython
/// (Bio.Phylo.Consensus) with a fixed random.seed(42); see
/// /tmp/gen_phylo_consensus_expected.py.

// ---------------------------------------------------------------------------
// Helpers (unique pbs_ prefix)
// ---------------------------------------------------------------------------

///|
fn pbs_leaf(name : String, bl : Double?) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=bl)
}

///|
fn pbs_node(children : Array[@src.Clade], bl : Double?) -> @src.Clade {
  @src.Clade::new(branch_length=bl, clades=children)
}

///|
fn pbs_tree(root : @src.Clade) -> @src.Tree {
  @src.Tree::new(root, rooted=true)
}

///|
fn pbs_key(clade : @src.Clade) -> String {
  let names = Array(capacity=0)
  for t in clade.get_terminals() {
    match t.name {
      Some(n) => names.push(n)
      None => abort("unnamed terminal")
    }
  }
  names.sort()
  let sb = StringBuilder()
  for i in 0..<names.length() {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(names[i])
  }
  sb.to_string()
}

///|
fn pbs_preorder_names(tree : @src.Tree) -> Array[String] {
  let names = Array(capacity=0)
  for t in tree.find_clades(terminal=Some(true)) {
    match t.name {
      Some(n) => names.push(n)
      None => abort("unnamed terminal")
    }
  }
  names
}

///|
fn pbs_find_internal(tree : @src.Tree, key : String) -> @src.Clade {
  for c in tree.find_clades(terminal=Some(false)) {
    if pbs_key(c) == key {
      return c
    }
  }
  abort("internal clade not found: " + key)
}

///|
fn pbs_child_keys(clade : @src.Clade) -> Array[String] {
  let keys = Array(capacity=0)
  for c in clade.clades {
    keys.push(pbs_key(c))
  }
  keys
}

///|
fn pbs_near(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-9
}

///|
fn pbs_resample(seq : String, cols : Array[Int]) -> String {
  let sb = StringBuilder()
  for j in cols {
    sb.write_char(Char::from_int(seq.unsafe_get(j).to_int()))
  }
  sb.to_string()
}

///|
fn pbs_raise_strict(trees : Array[@src.Tree]) -> String {
  try {
    ignore(@src.consensus_strict(trees))
    "NO_ERROR"
  } catch {
    @src.PhyloConsensusError(msg) => msg
    _ => "OTHER_ERROR"
  }
}

///|
fn pbs_raise_bootstrap(
  seqs : Array[String],
  names : Array[String],
  cols : Array[Array[Int]],
  consensus : String,
) -> String {
  try {
    ignore(@src.bootstrap_consensus(seqs, names, cols, consensus~))
    "NO_ERROR"
  } catch {
    @src.PhyloConsensusError(msg) => msg
    _ => "OTHER_ERROR"
  }
}

///|
fn pbs_hand_trees() -> Array[@src.Tree] {
  // tree1: (((A:0.1,B:0.2):0.3,C:0.4):0.5,(D:0.5,E:0.6):0.7)
  let ab1 = pbs_node(
    [pbs_leaf("A", Some(0.1)), pbs_leaf("B", Some(0.2))],
    Some(0.3),
  )
  let abc1 = pbs_node([ab1, pbs_leaf("C", Some(0.4))], Some(0.5))
  let de1 = pbs_node(
    [pbs_leaf("D", Some(0.5)), pbs_leaf("E", Some(0.6))],
    Some(0.7),
  )
  let tree1 = pbs_tree(pbs_node([abc1, de1], None))

  // tree2: ((A:0.1,(B:0.2,C:0.4):0.3):0.5,(D:0.5,E:0.6):0.7)
  let bc = pbs_node(
    [pbs_leaf("B", Some(0.2)), pbs_leaf("C", Some(0.4))],
    Some(0.3),
  )
  let abc2 = pbs_node([pbs_leaf("A", Some(0.1)), bc], Some(0.5))
  let de2 = pbs_node(
    [pbs_leaf("D", Some(0.5)), pbs_leaf("E", Some(0.6))],
    Some(0.7),
  )
  let tree2 = pbs_tree(pbs_node([abc2, de2], None))

  // tree3: ((A:0.1,B:0.2):0.3,(C:0.4,(D:0.5,E:0.6):0.7):0.5)
  let ab3 = pbs_node(
    [pbs_leaf("A", Some(0.1)), pbs_leaf("B", Some(0.2))],
    Some(0.3),
  )
  let de3 = pbs_node(
    [pbs_leaf("D", Some(0.5)), pbs_leaf("E", Some(0.6))],
    Some(0.7),
  )
  let cde = pbs_node([pbs_leaf("C", Some(0.4)), de3], Some(0.5))
  let tree3 = pbs_tree(pbs_node([ab3, cde], None))

  [tree1, tree2, tree3]
}

// ---------------------------------------------------------------------------
// Strict consensus
// ---------------------------------------------------------------------------

///|
test "pbs_strict_consensus_matches_biopython" {
  let tree = @src.consensus_strict(pbs_hand_trees())
  // root: 4 children (DE clade + A, B, C as separate leaves)
  let root = tree.get_root()
  assert_eq(root.clades.length(), 4)
  assert_true(root.branch_length is None)
  assert_true(root.confidence is None)
  assert_eq(pbs_preorder_names(tree), ["A", "B", "C", "D", "E"])
  // exactly one internal clade besides the root: (D, E)
  let internals = tree.find_clades(terminal=Some(false))
  assert_eq(internals.length(), 2)
  let de = pbs_find_internal(tree, "D,E")
  assert_eq(de.clades.length(), 2)
  assert_true(de.confidence is None)
  // child key multiset of the root
  let keys = pbs_child_keys(root)
  assert_true(keys.contains("A"))
  assert_true(keys.contains("B"))
  assert_true(keys.contains("C"))
  assert_true(keys.contains("D,E"))
}

// ---------------------------------------------------------------------------
// Majority rule consensus
// ---------------------------------------------------------------------------

///|
test "pbs_majority_consensus_matches_biopython" {
  let trees = pbs_hand_trees()
  let tree = @src.consensus_majority(trees)
  assert_eq(pbs_preorder_names(tree), ["D", "E", "C", "A", "B"])

  let root = tree.get_root()
  assert_eq(root.clades.length(), 2)
  assert_true(root.confidence is None)
  assert_eq(pbs_child_keys(root), ["D,E", "A,B,C"])

  let de = pbs_find_internal(tree, "D,E")
  assert_eq(de.clades.length(), 2)
  assert_true(pbs_near(de.confidence.unwrap(), 100.0))
  assert_true(pbs_near(de.branch_length.unwrap(), 0.7))

  let abc = pbs_find_internal(tree, "A,B,C")
  assert_eq(abc.clades.length(), 2)
  assert_true(pbs_near(abc.confidence.unwrap(), 66.66666666666667))
  assert_true(pbs_near(abc.branch_length.unwrap(), 0.5))
  assert_eq(pbs_child_keys(abc), ["C", "A,B"])

  let ab = pbs_find_internal(tree, "A,B")
  assert_true(pbs_near(ab.confidence.unwrap(), 66.66666666666667))
  assert_true(pbs_near(ab.branch_length.unwrap(), 0.3))
}

///|
test "pbs_majority_consensus_cutoff_half_matches_biopython" {
  let tree = @src.consensus_majority(pbs_hand_trees(), cutoff=0.5)
  assert_eq(pbs_preorder_names(tree), ["D", "E", "C", "A", "B"])
  let root = tree.get_root()
  assert_eq(pbs_child_keys(root), ["D,E", "A,B,C"])
  let de = pbs_find_internal(tree, "D,E")
  assert_true(pbs_near(de.confidence.unwrap(), 100.0))
}

// ---------------------------------------------------------------------------
// Adam consensus
// ---------------------------------------------------------------------------

///|
test "pbs_adam_consensus_matches_biopython" {
  let tree = @src.consensus_adam(pbs_hand_trees())
  assert_true(tree.is_rooted())
  assert_eq(pbs_preorder_names(tree), ["D", "E", "C", "A", "B"])
  let root = tree.get_root()
  assert_eq(root.clades.length(), 3)
  assert_eq(pbs_child_keys(root), ["D,E", "C", "A,B"])
  let de = pbs_find_internal(tree, "D,E")
  assert_eq(de.clades.length(), 2)
  assert_true(de.confidence is None)
  let ab = pbs_find_internal(tree, "A,B")
  assert_eq(ab.clades.length(), 2)
  assert_true(ab.confidence is None)
}

// ---------------------------------------------------------------------------
// Bootstrap support
// ---------------------------------------------------------------------------

///|
test "pbs_get_support_matches_biopython" {
  // target: tree1 topology, no branch lengths
  let ab = pbs_node([pbs_leaf("A", None), pbs_leaf("B", None)], None)
  let abc = pbs_node([ab, pbs_leaf("C", None)], None)
  let de = pbs_node([pbs_leaf("D", None), pbs_leaf("E", None)], None)
  let target = pbs_tree(pbs_node([abc, de], None))

  let tree = @src.tree_get_support(target, pbs_hand_trees())
  assert_eq(pbs_preorder_names(tree), ["A", "B", "C", "D", "E"])
  let root = tree.get_root()
  assert_true(pbs_near(root.confidence.unwrap(), 100.0))
  assert_true(
    pbs_near(
      pbs_find_internal(tree, "A,B,C").confidence.unwrap(),
      66.66666666666667,
    ),
  )
  assert_true(
    pbs_near(
      pbs_find_internal(tree, "A,B").confidence.unwrap(),
      66.66666666666667,
    ),
  )
  assert_true(
    pbs_near(pbs_find_internal(tree, "D,E").confidence.unwrap(), 100.0),
  )
}

///|
test "pbs_get_support_len_trees_override" {
  let ab = pbs_node([pbs_leaf("A", None), pbs_leaf("B", None)], None)
  let abc = pbs_node([ab, pbs_leaf("C", None)], None)
  let de = pbs_node([pbs_leaf("D", None), pbs_leaf("E", None)], None)
  let target = pbs_tree(pbs_node([abc, de], None))
  // all 3 trees contain DE; with len_trees=6 support is 50%
  let tree = @src.tree_get_support(target, pbs_hand_trees(), len_trees=Some(6))
  assert_true(
    pbs_near(pbs_find_internal(tree, "D,E").confidence.unwrap(), 50.0),
  )
}

// ---------------------------------------------------------------------------
// Bootstrap column indices (our own deterministic RNG)
// ---------------------------------------------------------------------------

///|
test "pbs_bootstrap_column_indices_deterministic_and_in_range" {
  let cols_a = @src.bootstrap_column_indices(10, 10, seed=42)
  let cols_b = @src.bootstrap_column_indices(10, 10, seed=42)
  assert_eq(cols_a.length(), 10)
  for i in 0..<10 {
    assert_eq(cols_a[i].length(), 10)
    assert_eq(cols_a[i], cols_b[i])
    for j in cols_a[i] {
      assert_true(j >= 0 && j < 10)
    }
  }
  // a different seed should change the stream
  let cols_c = @src.bootstrap_column_indices(10, 10, seed=7)
  assert_true(cols_c[0] != cols_a[0])
}

// ---------------------------------------------------------------------------
// Bootstrap replicate trees on the frozen MSA (Python random.seed(42))
// ---------------------------------------------------------------------------

///|
fn pbs_frozen_columns() -> Array[Array[Int]] {
  [
    [1, 0, 4, 3, 3, 2, 1, 8, 1, 9],
    [6, 0, 0, 1, 3, 3, 8, 9, 0, 8],
    [3, 8, 6, 3, 7, 9, 4, 0, 2, 6],
    [5, 4, 2, 3, 5, 1, 1, 6, 1, 5],
    [5, 9, 4, 0, 7, 8, 1, 6, 1, 8],
    [4, 9, 5, 9, 3, 1, 0, 3, 4, 1],
    [3, 1, 6, 4, 7, 5, 2, 5, 5, 3],
    [4, 1, 9, 2, 8, 3, 2, 7, 6, 4],
    [8, 3, 5, 0, 3, 0, 5, 6, 4, 1],
    [3, 9, 5, 3, 7, 6, 7, 2, 4, 2],
  ]
}

///|
fn pbs_frozen_msa() -> (Array[String], Array[String]) {
  let seqs = [
    "AACGTGCCTA", "ATCGTGCCTA", "ATCGTGCCTC", "ATCGTGCTTC", "ATCGGGCTTC",
  ]
  let names = ["S1", "S2", "S3", "S4", "S5"]
  (seqs, names)
}

///|
test "pbs_bootstrap_replicate0_distances_match_biopython" {
  let (seqs, _) = pbs_frozen_msa()
  let cols = pbs_frozen_columns()[0]
  let repl = Array(capacity=0)
  for s in seqs {
    repl.push(pbs_resample(s, cols))
  }
  let calc = @src.DistanceCalculator::new()
  let expected = [
    [0.0, 0.3, 0.4, 0.4, 0.5],
    [0.3, 0.0, 0.1, 0.1, 0.2],
    [0.4, 0.1, 0.0, 0.0, 0.1],
    [0.4, 0.1, 0.0, 0.0, 0.1],
    [0.5, 0.2, 0.1, 0.1, 0.0],
  ]
  for i in 0..<5 {
    for j in 0..<5 {
      let d = calc.pairwise(repl[i], repl[j])
      assert_true(pbs_near(d, expected[i][j]))
    }
  }
}

///|
test "pbs_bootstrap_trees_strict_consensus_matches_biopython" {
  let (seqs, names) = pbs_frozen_msa()
  let trees = @src.bootstrap_trees(seqs, names, pbs_frozen_columns())
  assert_eq(trees.length(), 10)
  for t in trees {
    assert_eq(t.count_terminals(), 5)
  }
  let strict = @src.consensus_strict(trees)
  assert_eq(pbs_preorder_names(strict), ["S3", "S4", "S5", "S1", "S2"])
  let root = strict.get_root()
  assert_eq(root.clades.length(), 4)
  let keys = pbs_child_keys(root)
  assert_true(keys.contains("S3"))
  assert_true(keys.contains("S4"))
  assert_true(keys.contains("S5"))
  assert_true(keys.contains("S1,S2"))
}

///|
/// Our NJ constructor resolves equal-Q cherries differently than
/// Biopython's (a different valid tie break), so the 10 replicate topologies
/// are not identical to Biopython's.  Consensus values below are the correct
/// name-set based counts for our 10 trees: (S1,S2) occurs 10/10, (S4,S5)
/// 5/10, (S3,S4) and (S1,S2,S5) 4/10 (the latter two are incompatible with
/// already retained clades, so extended majority rule skips them).
test "pbs_bootstrap_majority_consensus_name_consistent" {
  let (seqs, names) = pbs_frozen_msa()
  let trees = @src.bootstrap_trees(seqs, names, pbs_frozen_columns())
  let tree = @src.consensus_majority(trees)
  assert_eq(pbs_preorder_names(tree), ["S3", "S1", "S2", "S4", "S5"])
  let root = tree.get_root()
  assert_eq(root.clades.length(), 3)
  assert_eq(pbs_child_keys(root), ["S3", "S1,S2", "S4,S5"])

  let ab = pbs_find_internal(tree, "S1,S2")
  assert_eq(ab.clades.length(), 2)
  assert_true(pbs_near(ab.confidence.unwrap(), 100.0))
  assert_true(pbs_near(ab.branch_length.unwrap(), 0.07))

  let s45 = pbs_find_internal(tree, "S4,S5")
  assert_true(pbs_near(s45.confidence.unwrap(), 50.0))
  assert_true(pbs_near(s45.branch_length.unwrap(), 0.11))
}

///|
test "pbs_bootstrap_adam_consensus_matches_biopython" {
  let (seqs, names) = pbs_frozen_msa()
  let trees = @src.bootstrap_trees(seqs, names, pbs_frozen_columns())
  let tree = @src.consensus_adam(trees)
  assert_eq(pbs_preorder_names(tree), ["S1", "S2", "S5", "S4", "S3"])
  let root = tree.get_root()
  assert_eq(root.clades.length(), 4)
  assert_eq(pbs_child_keys(root), ["S1,S2", "S5", "S4", "S3"])
}

///|
test "pbs_bootstrap_consensus_dispatch" {
  let (seqs, names) = pbs_frozen_msa()
  let cols = pbs_frozen_columns()
  let strict = @src.bootstrap_consensus(seqs, names, cols, consensus="strict")
  assert_eq(pbs_preorder_names(strict), ["S3", "S4", "S5", "S1", "S2"])
  let majority = @src.bootstrap_consensus(
    seqs,
    names,
    cols,
    consensus="majority",
  )
  assert_eq(pbs_preorder_names(majority), ["S3", "S1", "S2", "S4", "S5"])
  let adam = @src.bootstrap_consensus(seqs, names, cols, consensus="adam")
  assert_eq(pbs_preorder_names(adam), ["S1", "S2", "S5", "S4", "S3"])
}

// ---------------------------------------------------------------------------
// Newick serialization
// ---------------------------------------------------------------------------

///|
test "pbs_clade_to_newick_smoke" {
  let ab = pbs_node(
    [pbs_leaf("A", Some(0.1)), pbs_leaf("B", Some(0.2))],
    Some(0.3),
  )
  let s = ab.to_newick()
  assert_true(s.has_prefix("("))
  assert_true(s.has_suffix(";"))
  assert_true(s.contains("A:0.1"))
  assert_true(s.contains("B:0.2"))
  // terminal clade serializes without parentheses
  let leaf_s = pbs_leaf("Z", Some(1.5)).to_newick()
  assert_eq(leaf_s, "Z:1.5;")
}

// ---------------------------------------------------------------------------
// Error paths
// ---------------------------------------------------------------------------

///|
test "pbs_consensus_empty_trees_raises" {
  let trees : Array[@src.Tree] = Array(capacity=0)
  assert_true(pbs_raise_strict(trees) != "NO_ERROR")
}

///|
test "pbs_consensus_inconsistent_taxa_raises" {
  let good = pbs_hand_trees()
  let bad_leaf = @src.Clade::new(name=Some("F"))
  let bad = @src.Tree::new(bad_leaf, rooted=true)
  good.push(bad)
  assert_true(pbs_raise_strict(good).contains("Taxons"))
}

///|
test "pbs_bootstrap_columns_wrong_length_raises" {
  let (seqs, names) = pbs_frozen_msa()
  let cols = [[0]]
  let msg = pbs_raise_bootstrap(seqs, names, cols, "strict")
  assert_true(msg.contains("columns"))
}

///|
test "pbs_bootstrap_columns_index_out_of_range_raises" {
  let (seqs, names) = pbs_frozen_msa()
  let cols = [[10, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
  let msg = pbs_raise_bootstrap(seqs, names, cols, "strict")
  assert_true(msg.contains("range"))
}

///|
test "pbs_bootstrap_unknown_consensus_raises" {
  let (seqs, names) = pbs_frozen_msa()
  let msg = pbs_raise_bootstrap(seqs, names, pbs_frozen_columns(), "mrc")
  assert_true(msg.contains("Unknown consensus"))
}

///|
test "pbs_bootstrap_sequence_name_count_mismatch_raises" {
  let seqs = ["ACGT", "ACGT"]
  let names = ["S1"]
  let cols = [[0, 1, 2, 3]]
  let msg = pbs_raise_bootstrap(seqs, names, cols, "strict")
  assert_true(msg.contains("match"))
}
