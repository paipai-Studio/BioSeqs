///|
/// Phylo consensus / bootstrap demo (Bio.Phylo.Consensus equivalents):
/// strict, majority-rule and Adam consensus trees, bootstrap support
/// annotation, and nonparametric bootstrap over a sequence alignment.
///
/// Run: moon run examples/phylo_bootstrap_demo
fn main raise {
  println("=== Bio.Phylo.Consensus Demo ===")
  println("")

  // Three rooted trees over the same five leaves (Biopython tutorial set)
  let trees = hand_trees()

  println("--- 1. Strict consensus ---")
  let strict = @bio.consensus_strict(trees)
  println("  " + strict.get_root().to_newick())
  println("")

  println("--- 2. Majority-rule consensus (percent support on labels) ---")
  let majority = @bio.consensus_majority(trees)
  println("  " + majority.get_root().to_newick(lengths=false))
  println("")

  println("--- 3. Adam (compatibility) consensus ---")
  let adam = @bio.consensus_adam(trees)
  println("  " + adam.get_root().to_newick(lengths=false))
  println("")

  println("--- 4. Bootstrap support on a target topology ---")
  let supported = @bio.tree_get_support(target_tree(), trees)
  println("  " + supported.get_root().to_newick(lengths=false))
  println("")

  println("--- 5. Alignment bootstrap (100 NJ replicates, identity model) ---")
  let names = ["S1", "S2", "S3", "S4", "S5"]
  let msa = [
    "AACGTGCCTA", "ATCGTGCCTA", "ATCGTGCCTC", "ATCGTGCTTC", "ATCGGGCTTC",
  ]
  for i in 0..<names.length() {
    println("  \{names[i]}  \{msa[i]}")
  }
  let columns = @bio.bootstrap_column_indices(msa[0].length(), 100, seed=42)
  let replicates = @bio.bootstrap_trees(msa, names, columns)
  println("  built \{replicates.length()} bootstrap replicate trees")
  let boot = @bio.bootstrap_consensus(msa, names, columns, consensus="majority")
  let boot_newick = boot.get_root().to_newick(lengths=false)
  println("  majority consensus: " + boot_newick)
}

///|
fn leaf(name : String, bl : Double?) -> @bio.Clade {
  @bio.Clade::new(name=Some(name), branch_length=bl)
}

///|
fn node(children : Array[@bio.Clade], bl : Double?) -> @bio.Clade {
  @bio.Clade::new(branch_length=bl, clades=children)
}

///|
fn rooted(root : @bio.Clade) -> @bio.Tree {
  @bio.Tree::new(root, rooted=true)
}

///|
fn hand_trees() -> Array[@bio.Tree] {
  // (((A:0.1,B:0.2):0.3,C:0.4):0.5,(D:0.5,E:0.6):0.7)
  let ab1 = node([leaf("A", Some(0.1)), leaf("B", Some(0.2))], Some(0.3))
  let abc1 = node([ab1, leaf("C", Some(0.4))], Some(0.5))
  let de1 = node([leaf("D", Some(0.5)), leaf("E", Some(0.6))], Some(0.7))
  let tree1 = rooted(node([abc1, de1], None))

  // ((A:0.1,(B:0.2,C:0.4):0.3):0.5,(D:0.5,E:0.6):0.7)
  let bc = node([leaf("B", Some(0.2)), leaf("C", Some(0.4))], Some(0.3))
  let abc2 = node([leaf("A", Some(0.1)), bc], Some(0.5))
  let de2 = node([leaf("D", Some(0.5)), leaf("E", Some(0.6))], Some(0.7))
  let tree2 = rooted(node([abc2, de2], None))

  // ((A:0.1,B:0.2):0.3,(C:0.4,(D:0.5,E:0.6):0.7):0.5)
  let ab3 = node([leaf("A", Some(0.1)), leaf("B", Some(0.2))], Some(0.3))
  let de3 = node([leaf("D", Some(0.5)), leaf("E", Some(0.6))], Some(0.7))
  let cde = node([leaf("C", Some(0.4)), de3], Some(0.5))
  let tree3 = rooted(node([ab3, cde], None))

  [tree1, tree2, tree3]
}

///|
fn target_tree() -> @bio.Tree {
  let ab = node([leaf("A", None), leaf("B", None)], None)
  let abc = node([ab, leaf("C", None)], None)
  let de = node([leaf("D", None), leaf("E", None)], None)
  rooted(node([abc, de], None))
}
