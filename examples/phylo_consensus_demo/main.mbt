///|
/// Phylo.Consensus example: Consensus tree building algorithms.
///
/// Run: moon run examples/phylo_consensus_demo/main.mbt
fn main {
  println("=== Bio.Phylo.Consensus Demo ===\n")

  println("1. Parse Newick format:")
  let newick1 = "((A,B),(C,D));"
  let tree1 = @bio.newick_to_tree(newick1)
  println("   Input: \{newick1}")
  match tree1.root {
    Some(root) => {
      println("   Root: internal node with \{root.children.length()} children")
      let leaves = root.get_leaves()
      let mut leaves_str = "   Leaves: "
      let mut i = 0
      while i < leaves.length() {
        if i > 0 {
          leaves_str = leaves_str + ", "
        }
        leaves_str = leaves_str + leaves[i]
        i = i + 1
      }
      println(leaves_str)
    }
    None => println("   Failed to parse")
  }

  println("\n2. Newick roundtrip:")
  let tree2 = @bio.newick_to_tree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
  let newick_out = tree2.to_newick()
  println("   Input: ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
  println("   Output: \{newick_out}")

  println("\n3. Create example trees:")
  let trees = @bio.create_example_trees()
  println("   Created \{trees.length()} trees for consensus analysis")
  let mut i = 0
  while i < trees.length() {
    println("   Tree \{i + 1}: \{trees[i].to_newick()}")
    i = i + 1
  }

  println("\n4. Extract splits from tree:")
  match trees[0].root {
    Some(root) => {
      let leaves = root.get_leaves()
      let splits = @bio.get_all_splits(root, leaves)
      println("   Splits in tree 1: \{splits.length()}")
      let mut j = 0
      while j < splits.length() {
        println("     Split \{j + 1}: \{splits[j].hash()}")
        j = j + 1
      }
    }
    None => ()
  }

  println("\n5. Build Majority Rule Consensus (threshold=0.5):")
  let majority = @bio.build_majority_consensus(trees)
  println("   Majority consensus tree:")
  println("     \{majority.to_newick()}")

  println("\n6. Build Strict Consensus (threshold=1.0):")
  let strict = @bio.build_strict_consensus(trees)
  println("   Strict consensus tree:")
  println("     \{strict.to_newick()}")

  println("\n7. Build custom threshold consensus:")
  let custom = @bio.build_consensus(trees, 0.6)
  println("   Consensus with threshold=0.6:")
  println("     \{custom.to_newick()}")

  println("\n8. Calculate support for specific split:")
  let split = @bio.Split::new(["A", "B"])
  let support = @bio.calculate_consensus_support(trees, split)
  println("   Support for split (A,B): \{support * 100.0}%")

  println("\n9. Get all unique splits from trees:")
  let all_splits = @bio.get_all_splits_from_trees(trees)
  println("   Unique splits across all trees: \{all_splits.length()}")
  let mut k = 0
  while k < all_splits.length() {
    let supp = @bio.calculate_consensus_support(trees, all_splits[k])
    println("     Split \{k + 1}: \{all_splits[k].hash()} (support: \{supp * 100.0}%)")
    k = k + 1
  }

  println("\n10. Tree structure visualization:")
  match majority.root {
    Some(root) => {
      println("   Majority consensus tree structure:")
      println(root.to_string())
    }
    None => ()
  }

  println("\n=== Demo Complete ===")
}