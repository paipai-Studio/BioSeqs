///|
/// Bio.Phylo.TreeMetrics - Advanced phylogenetic tree metrics.
/// Inspired by Bio.Phylo's tree statistics and comparison functions.

///|
/// Calculate the total branch length of a clade (sum of all branch lengths).
pub fn Clade::total_branch_length(self : Clade) -> Double {
  let mut total = self.branch_length.unwrap_or(0.0)
  for clade in self.clades {
    total = total + clade.total_branch_length()
  }
  total
}

///|
/// Calculate the maximum depth (height) of the tree from root to deepest leaf.
pub fn Clade::max_depth(self : Clade) -> Double {
  let depth = self.branch_length.unwrap_or(0.0)
  if self.clades.length() == 0 {
    return depth
  }
  let mut max_child = 0.0
  for clade in self.clades {
    let child_depth = clade.max_depth()
    if child_depth > max_child {
      max_child = child_depth
    }
  }
  depth + max_child
}

///|
/// Calculate the maximum depth of a Tree.
pub fn Tree::max_depth(self : Tree) -> Double {
  self.root.max_depth()
}

///|
/// Calculate the total branch length of a Tree.
pub fn Tree::total_branch_length(self : Tree) -> Double {
  self.root.total_branch_length()
}

///|
/// Get all bipartitions (splits) in the tree.
/// A bipartition is a set of two leaf sets created by removing an edge.
/// Returns an array of (left_leaves, right_leaves) tuples.
pub fn Clade::get_bipartitions(
  self : Clade,
) -> Array[(Array[String], Array[String])] {
  let result : Array[(Array[String], Array[String])] = Array::new()
  let all_leaves = self.get_leaf_names()
  self.collect_bipartitions(all_leaves, result)
  result
}

///|
/// Helper to collect bipartitions.
fn Clade::collect_bipartitions(
  self : Clade,
  all_leaves : Array[String],
  result : Array[(Array[String], Array[String])],
) -> Unit {
  for clade in self.clades {
    let clade_leaves = clade.get_leaf_names()
    if clade_leaves.length() > 0 && clade_leaves.length() < all_leaves.length() {
      // Build the complement set
      let complement : Array[String] = Array::new()
      let mut i = 0
      while i < all_leaves.length() {
        let mut found = false
        let mut j = 0
        while j < clade_leaves.length() {
          if all_leaves[i] == clade_leaves[j] {
            found = true
            break
          }
          j = j + 1
        }
        if not(found) {
          complement.push(all_leaves[i])
        }
        i = i + 1
      }
      result.push((clade_leaves, complement))
    }
    clade.collect_bipartitions(all_leaves, result)
  }
}

///|
/// Get all leaf names (terminal clade names) in this subtree.
pub fn Clade::get_leaf_names(self : Clade) -> Array[String] {
  let result : Array[String] = Array::new()
  self.collect_leaf_names(result)
  result
}

///|
/// Helper to collect leaf names.
fn Clade::collect_leaf_names(self : Clade, result : Array[String]) -> Unit {
  if self.is_terminal() {
    match self.name {
      Some(n) => result.push(n)
      None => ()
    }
  } else {
    for clade in self.clades {
      clade.collect_leaf_names(result)
    }
  }
}

///|
/// Calculate the Robinson-Foulds distance between two trees.
/// RF distance = number of bipartitions in tree1 but not tree2
///              + number of bipartitions in tree2 but not in tree1
pub fn robinson_foulds_distance(tree1 : Tree, tree2 : Tree) -> Int {
  let biparts1 = tree1.root.get_bipartitions()
  let biparts2 = tree2.root.get_bipartitions()

  // Normalize bipartitions (sort each side, then take the smaller set as "left")
  let normalized1 = normalize_bipartitions(biparts1)
  let normalized2 = normalize_bipartitions(biparts2)

  // Count shared bipartitions
  let mut shared = 0
  for bp1 in normalized1 {
    for bp2 in normalized2 {
      if bipartition_equals(bp1, bp2) {
        shared = shared + 1
        break
      }
    }
  }

  normalized1.length() - shared + (normalized2.length() - shared)
}

///|
/// Normalize a bipartition by ensuring the smaller set is always first.
fn normalize_bipartitions(
  biparts : Array[(Array[String], Array[String])],
) -> Array[(Array[String], Array[String])] {
  let result : Array[(Array[String], Array[String])] = Array::new()
  for bp in biparts {
    let (left, right) = bp
    if left.length() <= right.length() {
      let sorted_left = sort_strings(left)
      let sorted_right = sort_strings(right)
      result.push((sorted_left, sorted_right))
    } else {
      let sorted_left = sort_strings(right)
      let sorted_right = sort_strings(left)
      result.push((sorted_left, sorted_right))
    }
  }
  result
}

///|
/// Sort an array of strings alphabetically (simple bubble sort).
fn sort_strings(arr : Array[String]) -> Array[String] {
  let result = arr.copy()
  let n = result.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      if result[j] > result[j + 1] {
        let temp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Check if two bipartitions are equal (after normalization).
fn bipartition_equals(
  bp1 : (Array[String], Array[String]),
  bp2 : (Array[String], Array[String]),
) -> Bool {
  let (left1, right1) = bp1
  let (left2, right2) = bp2
  arrays_equal(left1, left2) && arrays_equal(right1, right2)
}

///|
/// Check if two string arrays are equal (same elements).
fn arrays_equal(a : Array[String], b : Array[String]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let mut i = 0
  while i < a.length() {
    if a[i] != b[i] {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Calculate the Colless index for tree balance.
/// Measures how balanced a tree is. Lower = more balanced.
pub fn Clade::colless_index(self : Clade) -> Int {
  if self.is_terminal() {
    return 0
  }
  let mut total = 0
  if self.clades.length() == 2 {
    let left = self.clades[0].count_terminals()
    let right = self.clades[1].count_terminals()
    let diff = if left > right { left - right } else { right - left }
    total = total + diff
  }
  for clade in self.clades {
    total = total + clade.colless_index()
  }
  total
}

///|
/// Calculate the Colless index for a Tree.
pub fn Tree::colless_index(self : Tree) -> Int {
  self.root.colless_index()
}

///|
/// Calculate patristic distance (sum of branch lengths along path) between two leaves.
pub fn Clade::patristic_distance(
  self : Clade,
  name1 : String,
  name2 : String,
) -> Double? {
  match (self.get_path(name1), self.get_path(name2)) {
    (Some(path1), Some(path2)) => {
      // Find common ancestor position by checking if paths diverge
      // Since both paths start from root, find the last common index
      let mut i = 0
      let min_len = if path1.length() < path2.length() {
        path1.length()
      } else {
        path2.length()
      }
      // Compare leaf names at the end of paths to detect divergence
      while i < min_len {
        // Check if the remaining paths lead to the same leaf
        let p1_name = path1[path1.length() - 1].name
        let p2_name = path2[path2.length() - 1].name
        if p1_name == p2_name {
          break
        }
        // Check if current nodes are the same by comparing their leaf sets
        let leaves1 = path1[i].get_leaf_names()
        let leaves2 = path2[i].get_leaf_names()
        let mut same = true
        if leaves1.length() == leaves2.length() {
          let mut k = 0
          while k < leaves1.length() {
            let mut found = false
            let mut m = 0
            while m < leaves2.length() {
              if leaves1[k] == leaves2[m] {
                found = true
                break
              }
              m = m + 1
            }
            if not(found) {
              same = false
              break
            }
            k = k + 1
          }
        } else {
          same = false
        }
        if not(same) {
          break
        }
        i = i + 1
      }
      // Sum branch lengths from common ancestor to each leaf
      let mut dist = 0.0
      for j in i..<path1.length() {
        dist = dist + path1[j].branch_length.unwrap_or(0.0)
      }
      for j in i..<path2.length() {
        dist = dist + path2[j].branch_length.unwrap_or(0.0)
      }
      Some(dist)
    }
    _ => None
  }
}

///|
/// Calculate patristic distance between two leaves in a Tree.
pub fn Tree::patristic_distance(
  self : Tree,
  name1 : String,
  name2 : String,
) -> Double? {
  self.root.patristic_distance(name1, name2)
}

///|
/// Compute the patristic distance matrix between all pairs of leaves.
pub fn Clade::patristic_distance_matrix(
  self : Clade,
) -> Array[(String, String, Double)] {
  let result : Array[(String, String, Double)] = Array::new()
  let leaves = self.get_leaf_names()
  let n = leaves.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      match self.patristic_distance(leaves[i], leaves[j]) {
        Some(d) => result.push((leaves[i], leaves[j], d))
        None => ()
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Compute the patristic distance matrix for a Tree.
pub fn Tree::patristic_distance_matrix(
  self : Tree,
) -> Array[(String, String, Double)] {
  self.root.patristic_distance_matrix()
}

///|
/// Count the number of internal nodes (non-terminal clades).
pub fn Clade::count_internal(self : Clade) -> Int {
  if self.is_terminal() {
    0
  } else {
    let mut count = 1
    for clade in self.clades {
      count = count + clade.count_internal()
    }
    count
  }
}

///|
/// Count the number of internal nodes in a Tree.
pub fn Tree::count_internal(self : Tree) -> Int {
  self.root.count_internal()
}

///|
/// Create a demo tree for testing.
/// Tree structure: ((A:1, B:1):0.5, (C:1, D:1):0.5)
pub fn create_demo_tree() -> Tree {
  let leaf_a = Clade::new(branch_length=Some(1.0), name=Some("A"))
  let leaf_b = Clade::new(branch_length=Some(1.0), name=Some("B"))
  let leaf_c = Clade::new(branch_length=Some(1.0), name=Some("C"))
  let leaf_d = Clade::new(branch_length=Some(1.0), name=Some("D"))

  let clade_ab = Clade::new(branch_length=Some(0.5), clades=[leaf_a, leaf_b])
  let clade_cd = Clade::new(branch_length=Some(0.5), clades=[leaf_c, leaf_d])
  let root = Clade::new(clades=[clade_ab, clade_cd])
  Tree::new(root, rooted=true, name=Some("demo"))
}

///|
/// Create a second demo tree for Robinson-Foulds testing.
/// Tree structure: ((A:1, C:1):0.5, (B:1, D:1):0.5)
pub fn create_demo_tree2() -> Tree {
  let leaf_a = Clade::new(branch_length=Some(1.0), name=Some("A"))
  let leaf_b = Clade::new(branch_length=Some(1.0), name=Some("B"))
  let leaf_c = Clade::new(branch_length=Some(1.0), name=Some("C"))
  let leaf_d = Clade::new(branch_length=Some(1.0), name=Some("D"))

  let clade_ac = Clade::new(branch_length=Some(0.5), clades=[leaf_a, leaf_c])
  let clade_bd = Clade::new(branch_length=Some(0.5), clades=[leaf_b, leaf_d])
  let root = Clade::new(clades=[clade_ac, clade_bd])
  Tree::new(root, rooted=true, name=Some("demo2"))
}
