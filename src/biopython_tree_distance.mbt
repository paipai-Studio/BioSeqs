/// biopython_tree_distance.mbt — Phylogenetic tree distance metrics
/// Inspired by Biopython Bio.Phylo: Robinson-Foulds, tree bipartitions

///|
/// A simple tree representation for distance calculations.
/// Tree is represented as a set of bipartitions (splits).
/// Each bipartition is a sorted set of taxon indices on one side.
pub struct BiopythonTreeBipartitions {
  taxa : Array[String]
  splits : Array[Array[Bool]] // each split is a boolean array (true = side A)
}

///|
/// Parse a simple Newick-like tree to extract bipartitions.
/// This is a simplified parser that handles unrooted trees.
/// Input: comma-separated taxa list for star tree, or Newick string.
pub fn biopython_tree_bipartitions(
  taxa : Array[String],
  newick : String,
) -> BiopythonTreeBipartitions {
  // For simplicity, extract clades from Newick by finding parentheses groups
  // Each internal node defines a bipartition
  let n = taxa.length()
  let splits : Array[Array[Bool]] = Array(capacity=0)
  // Simple approach: scan for taxon names in each parenthetical group
  let mut depth = 0
  let mut group_start = -1
  for i = 0; i < newick.length(); i = i + 1 {
    let ch = newick[i].to_int()
    if ch == 40 { // '('
      depth = depth + 1
      if depth == 1 {
        group_start = i + 1
      }
    } else if ch == 41 { // ')'
      if depth == 1 && group_start >= 0 {
        // Extract taxa in this group
        let group = newick[group_start:i].to_owned()
        let split : Array[Bool] = Array::make(n, false)
        for t = 0; t < n; t = t + 1 {
          if group.contains(taxa[t]) {
            split[t] = true
          }
        }
        // Only add non-trivial splits
        let mut count = 0
        for t = 0; t < n; t = t + 1 {
          if split[t] {
            count = count + 1
          }
        }
        if count > 1 && count < n - 1 {
          splits.push(split)
        }
      }
      depth = depth - 1
    }
  }
  { taxa, splits, }
}

///|
/// Robinson-Foulds distance between two trees.
/// Count of bipartitions in one tree but not the other (symmetric difference).
pub fn biopython_robinson_foulds(
  tree1 : BiopythonTreeBipartitions,
  tree2 : BiopythonTreeBipartitions,
) -> Int32 {
  let n = tree1.taxa.length()
  let mut dist = 0
  // Check splits in tree1 not in tree2
  for i = 0; i < tree1.splits.length(); i = i + 1 {
    if !biopython_split_in_tree(tree1.splits[i], tree2, n) {
      dist = dist + 1
    }
  }
  // Check splits in tree2 not in tree1
  for i = 0; i < tree2.splits.length(); i = i + 1 {
    if !biopython_split_in_tree(tree2.splits[i], tree1, n) {
      dist = dist + 1
    }
  }
  dist
}

///|
/// Check if a split exists in a tree (considering both orientations).
fn biopython_split_in_tree(
  split : Array[Bool],
  tree : BiopythonTreeBipartitions,
  n : Int32,
) -> Bool {
  for i = 0; i < tree.splits.length(); i = i + 1 {
    let other = tree.splits[i]
    // Check if splits match (same or complement)
    let mut same = true
    let mut comp = true
    for j = 0; j < n; j = j + 1 {
      if split[j] != other[j] {
        same = false
      }
      if split[j] == other[j] {
        comp = false
      }
    }
    if same || comp {
      return true
    }
  }
  false
}

///|
/// Normalized Robinson-Foulds distance (0 to 1).
pub fn biopython_robinson_foulds_normalized(
  tree1 : BiopythonTreeBipartitions,
  tree2 : BiopythonTreeBipartitions,
) -> Double {
  let rf = biopython_robinson_foulds(tree1, tree2)
  let max_rf = tree1.splits.length() + tree2.splits.length()
  if max_rf == 0 {
    0.0
  } else {
    rf.to_double() / max_rf.to_double()
  }
}

///|
/// Branch score distance (sum of squared branch length differences).
/// Simplified: uses topology only (branch lengths = 1).
pub fn biopython_branch_score_distance(
  tree1 : BiopythonTreeBipartitions,
  tree2 : BiopythonTreeBipartitions,
) -> Double {
  let n = tree1.taxa.length()
  let mut dist = 0.0
  // For each split in tree1 not in tree2, add 1.0 (squared branch length)
  for i = 0; i < tree1.splits.length(); i = i + 1 {
    if !biopython_split_in_tree(tree1.splits[i], tree2, n) {
      dist = dist + 1.0
    }
  }
  for i = 0; i < tree2.splits.length(); i = i + 1 {
    if !biopython_split_in_tree(tree2.splits[i], tree1, n) {
      dist = dist + 1.0
    }
  }
  dist
}
