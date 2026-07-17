///|
/// Phylogenetic tree example: parse Newick, traverse, and analyze.
///
/// Run: moon run examples/phylo_demo/main.mbt
fn main {
  // === Parse Newick tree ===
  let newick = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
  let tree = @bio.parse_newick(newick)

  println("=== Newick Tree ===")
  println("  Input: \{newick}")

  // === Tree statistics ===
  let leaf_count = @bio.get_leaf_count(tree)
  println("\n=== Tree Statistics ===")
  println("  Terminal nodes (leaves): \{leaf_count}")

  let node_count = @bio.get_node_count(tree)
  let nonterminals = node_count - leaf_count
  println("  Non-terminal nodes: \{nonterminals}")

  // === Get terminal clades ===
  let leaves = @bio.get_leaves(tree)
  println("\n=== Terminal Clades ===")
  for i = 0; i < leaves.length(); i = i + 1 {
    println("  Leaf \{i}: \{leaves[i]}")
  }

  // === Find node by name ===
  let found = @bio.find_node(tree, "A")
  match found {
    Some(n) => println("\n=== Find Node 'A' ===\n  Found: \{n.name}")
    None => println("\n=== Find Node 'A' ===\n  Not found")
  }

  // === Calculate tree depth ===
  let depth = @bio.get_tree_depth(tree)
  println("\n=== Tree Depth ===")
  println("  Maximum depth: \{depth}")

  // === Total branch length ===
  let total_len = @bio.total_branch_length(tree)
  println("\n=== Total Branch Length ===")
  println("  \{total_len}")

  // === Write Newick ===
  let newick_out = @bio.write_newick(tree, true)
  println("\n=== Serialized Newick ===")
  println("  \{newick_out}")

  // === Prune tree ===
  let pruned = @bio.prune_tree(tree, ["A", "B"])
  match pruned {
    Some(p) => {
      let pruned_newick = @bio.write_newick(p, true)
      println("\n=== Pruned Tree (A, B) ===")
      println("  \{pruned_newick}")
    }
    None => println("\n=== Pruned Tree ===\n  Failed")
  }
}
