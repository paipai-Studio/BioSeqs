///|
/// DendroPy `calculate` port — phylogenetic tree measures (balance and shape
/// statistics, tree height, and the unweighted Robinson–Foulds distance),
/// operating on BioSeqs' `phylo.Clade` trees.

///|
fn dendropy_count_leaves(clade : Clade) -> Int {
  if clade.is_terminal() {
    1
  } else {
    let mut c = 0
    for ch in clade.clades {
      c = c + dendropy_count_leaves(ch)
    }
    c
  }
}

///|
/// Colless' tree imbalance index (un-normalized): the sum over internal
/// bifurcating nodes of the absolute difference between the leaf counts of
/// their two subtrees.
pub fn dendropy_colless(clade : Clade) -> Int {
  dendropy_colless_rec(clade).1
}

///|
fn dendropy_colless_rec(clade : Clade) -> (Int, Int) {
  if clade.is_terminal() {
    return (1, 0)
  }
  let mut total = 0
  let mut csum = 0
  let child_leaves : Array[Int] = Array::new()
  for ch in clade.clades {
    let (leaves, cs) = dendropy_colless_rec(ch)
    total = total + leaves
    csum = csum + cs
    child_leaves.push(leaves)
  }
  if child_leaves.length() == 2 {
    csum = csum + (child_leaves[0] - child_leaves[1]).abs()
  }
  (total, csum)
}

///|
/// Sackin's index (un-normalized): the sum over terminal nodes of their depth
/// (number of ancestors, i.e., edges from the root).
pub fn dendropy_sackin(clade : Clade) -> Int {
  dendropy_sackin_rec(clade, 0)
}

///|
fn dendropy_sackin_rec(clade : Clade, depth : Int) -> Int {
  if clade.is_terminal() {
    depth
  } else {
    let mut s = 0
    for ch in clade.clades {
      s = s + dendropy_sackin_rec(ch, depth + 1)
    }
    s
  }
}

///|
/// The N̄ statistic: the average number of ancestors per terminal node
/// (Sackin's index normalized by the number of leaves).
pub fn dendropy_n_bar(clade : Clade) -> Double {
  let n = dendropy_count_leaves(clade)
  if n == 0 {
    0.0
  } else {
    dendropy_sackin(clade).to_double() / n.to_double()
  }
}

///|
/// The B1 statistic: the sum over internal nodes (excluding the root) of the
/// reciprocal of the maximum number of nodes between that node and a tip.
pub fn dendropy_b1(clade : Clade) -> Double {
  dendropy_b1_rec(clade, true).1
}

///|
fn dendropy_b1_rec(clade : Clade, is_root : Bool) -> (Double, Double) {
  if clade.is_terminal() {
    return (0.0, 0.0)
  }
  let mut max_child = 0.0
  let mut b1_sum = 0.0
  for ch in clade.clades {
    let (mi, bs) = dendropy_b1_rec(ch, false)
    if mi > max_child {
      max_child = mi
    }
    b1_sum = b1_sum + bs
  }
  let mi = max_child + 1.0
  if is_root {
    (mi, b1_sum)
  } else {
    (mi, b1_sum + 1.0 / mi)
  }
}

///|
/// Tree height: the maximum root-to-tip sum of branch lengths.
pub fn dendropy_tree_height(clade : Clade) -> Double {
  dendropy_height_rec(clade)
}

///|
fn dendropy_height_rec(clade : Clade) -> Double {
  if clade.is_terminal() {
    0.0
  } else {
    let mut max_h = 0.0
    for ch in clade.clades {
      let bl = match ch.branch_length {
        None => 0.0
        Some(b) => b
      }
      let h = dendropy_height_rec(ch) + bl
      if h > max_h {
        max_h = h
      }
    }
    max_h
  }
}

///|
fn dendropy_sort_strings(arr : Array[String]) -> Array[String] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] > key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn dendropy_join(names : Array[String]) -> String {
  let sb = StringBuilder::new()
  let mut first = true
  for n in names {
    if not(first) {
      sb.write_char(',')
    }
    sb.write_string(n)
    first = false
  }
  sb.to_string()
}

///|
/// Sorted leaf names of a subtree.
fn dendropy_leaf_names(clade : Clade) -> Array[String] {
  if clade.is_terminal() {
    match clade.name {
      None => []
      Some(n) => [n]
    }
  } else {
    let out : Array[String] = Array::new()
    for ch in clade.clades {
      for n in dendropy_leaf_names(ch) {
        out.push(n)
      }
    }
    dendropy_sort_strings(out)
  }
}

///|
/// Complement of `sub` within `all` (both sorted).
fn dendropy_complement(all : Array[String], sub : Array[String]) -> Array[String] {
  let out : Array[String] = Array::new()
  let mut si = 0
  for x in all {
    if si < sub.length() && sub[si] == x {
      si = si + 1
    } else {
      out.push(x)
    }
  }
  out
}

///|
/// Collect canonical split keys (non-trivial bipartitions) of a tree.
fn dendropy_split_rec(
  clade : Clade,
  all_names : Array[String],
  n : Int,
) -> (Array[String], Map[String, Bool]) {
  if clade.is_terminal() {
    match clade.name {
      None => ([], Map([], capacity=0))
      Some(x) =>
        ([x], Map([], capacity=0))
    }
  } else {
    let keys : Map[String, Bool] = Map([], capacity=0)
    let sub : Array[String] = Array::new()
    for ch in clade.clades {
      let (names, child_keys) = dendropy_split_rec(ch, all_names, n)
      for nm in names {
        sub.push(nm)
      }
      // merge child_keys into keys
      for (k, _) in child_keys {
        keys.set(k, true)
      }
    }
    let sorted_sub = dendropy_sort_strings(sub)
    let m = sorted_sub.length()
    if m > 1 && m < n {
      let comp = dendropy_complement(all_names, sorted_sub)
      let key_sub = dendropy_join(sorted_sub)
      let key_comp = dendropy_join(comp)
      let canonical = if m < comp.length() {
        key_sub
      } else if m > comp.length() {
        key_comp
      } else if key_sub < key_comp {
        key_sub
      } else {
        key_comp
      }
      keys.set(canonical, true)
    }
    (sorted_sub, keys)
  }
}

///|
/// Unweighted Robinson–Foulds distance (symmetric difference) between two
/// trees on the same leaf set (`dendropy.calculate.treecompare.symmetric_difference`).
pub fn dendropy_symmetric_difference(clade1 : Clade, clade2 : Clade) -> Int {
  let all_names = dendropy_leaf_names(clade1)
  let n = all_names.length()
  let (_, s1) = dendropy_split_rec(clade1, all_names, n)
  let (_, s2) = dendropy_split_rec(clade2, all_names, n)
  // count keys in s1 not in s2 + keys in s2 not in s1
  let mut count = 0
  for (k, _) in s1 {
    match s2.get(k) {
      None => count = count + 1
      Some(_) => ()
    }
  }
  for (k, _) in s2 {
    match s1.get(k) {
      None => count = count + 1
      Some(_) => ()
    }
  }
  count
}