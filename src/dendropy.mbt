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
  let child_leaves : Array[Int] = Array(capacity=0)
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
    let out : Array[String] = Array(capacity=0)
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
  let out : Array[String] = Array(capacity=0)
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
    let sub : Array[String] = Array(capacity=0)
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
}///|
/// Canonical split → branch-length map: for each non-trivial bipartition, the
/// length of the branch subtending the defining internal node (0.0 if absent).
fn dendropy_split_lengths_rec(
  clade : Clade,
  all_names : Array[String],
  n : Int,
) -> (Array[String], Map[String, Double]) {
  if clade.is_terminal() {
    match clade.name {
      None => ([], Map([], capacity=0))
      Some(x) => ([x], Map([], capacity=0))
    }
  } else {
    let keys : Map[String, Double] = Map([], capacity=0)
    let sub : Array[String] = Array(capacity=0)
    for ch in clade.clades {
      let (names, child_keys) = dendropy_split_lengths_rec(ch, all_names, n)
      for nm in names {
        sub.push(nm)
      }
      for (k, v) in child_keys {
        keys.set(k, v)
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
      let bl = match clade.branch_length {
        None => 0.0
        Some(b) => b
      }
      keys.set(canonical, bl)
    }
    (sorted_sub, keys)
  }
}

///|
/// Edge-weighted Robinson–Foulds distance: the sum, over all bipartitions in
/// either tree, of the absolute difference in branch lengths (a missing
/// bipartition contributes its own length).
pub fn dendropy_weighted_robinson_foulds(clade1 : Clade, clade2 : Clade) -> Double {
  let all = dendropy_leaf_names(clade1)
  let n = all.length()
  let (_, m1) = dendropy_split_lengths_rec(clade1, all, n)
  let (_, m2) = dendropy_split_lengths_rec(clade2, all, n)
  let mut d = 0.0
  let seen : Map[String, Bool] = Map([], capacity=0)
  for (k, v1) in m1 {
    let v2 = match m2.get(k) {
      None => 0.0
      Some(x) => x
    }
    d = d + (v1 - v2).abs()
    seen.set(k, true)
  }
  for (k, v2) in m2 {
    match seen.get(k) {
      Some(_) => ()
      None => d = d + v2.abs()
    }
  }
  d
}

///|
/// Find the lowest common ancestor of two leaf names.
fn dendropy_lca_find(clade : Clade, name1 : String, name2 : String) -> Clade? {
  if clade.is_terminal() {
    match clade.name {
      Some(nn) => if nn == name1 || nn == name2 { Some(clade) } else { None }
      None => None
    }
  } else {
    let found : Array[Clade] = Array(capacity=0)
    for ch in clade.clades {
      match dendropy_lca_find(ch, name1, name2) {
        None => ()
        Some(c) => found.push(c)
      }
    }
    if found.length() == 2 {
      Some(clade)
    } else if found.length() == 1 {
      Some(found[0])
    } else {
      None
    }
  }
}

///|
/// Distance from `clade` down to leaf `name`, accumulating branch lengths.
fn dendropy_down(clade : Clade, name : String, acc : Double) -> Double? {
  if clade.is_terminal() {
    match clade.name {
      Some(nn) => if nn == name { Some(acc) } else { None }
      None => None
    }
  } else {
    for ch in clade.clades {
      let bl = match ch.branch_length {
        None => 0.0
        Some(b) => b
      }
      match dendropy_down(ch, name, acc + bl) {
        None => ()
        Some(d) => return Some(d)
      }
    }
    None
  }
}

///|
/// Patristic distance between two taxa: the sum of branch lengths along the
/// path between them (`dendropy.calculate.treemeasure.patristic_distance`).
pub fn dendropy_patristic_distance(clade : Clade, name1 : String, name2 : String) -> Double {
  match dendropy_lca_find(clade, name1, name2) {
    None => 0.0
    Some(lca) =>
      match (dendropy_down(lca, name1, 0.0), dendropy_down(lca, name2, 0.0)) {
        (Some(d1), Some(d2)) => d1 + d2
        _ => 0.0
      }
  }
}

///|
/// Pairwise patristic-distance matrix, returning the sorted leaf names and the
/// dense distance matrix.
pub fn dendropy_patristic_distance_matrix(
  clade : Clade,
) -> (Array[String], Array[Array[Double]]) {
  let names = dendropy_leaf_names(clade)
  let nn = names.length()
  let matrix : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < nn; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < nn; j = j + 1 {
      row.push(dendropy_patristic_distance(clade, names[i], names[j]))
    }
    matrix.push(row)
  }
  (names, matrix)
}

///|
fn dendropy_sort_doubles_desc(arr : Array[Double]) -> Array[Double] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] < key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Node height (distance to the deepest leaf); for an ultrametric tree this is
/// uniform over the leaves below the node.
fn dendropy_node_height(clade : Clade) -> Double {
  if clade.is_terminal() {
    0.0
  } else {
    let mut max_h = 0.0
    for ch in clade.clades {
      let bl = match ch.branch_length {
        None => 0.0
        Some(b) => b
      }
      let h = dendropy_node_height(ch) + bl
      if h > max_h {
        max_h = h
      }
    }
    max_h
  }
}

///|
/// Collect the heights of bifurcating internal nodes (speciation times) and
/// return the number of leaves.
fn dendropy_gamma_collect(clade : Clade, ages : Array[Double]) -> Int {
  if clade.is_terminal() {
    1
  } else {
    let mut total = 0
    for ch in clade.clades {
      total = total + dendropy_gamma_collect(ch, ages)
    }
    if clade.clades.length() == 2 {
      ages.push(dendropy_node_height(clade))
    }
    total
  }
}

///|
/// Pybus & Harvey (2000) gamma statistic of diversification-rate constancy.
/// Assumes an ultrametric, bifurcating tree.
pub fn dendropy_pybus_harvey_gamma(clade : Clade) -> Double {
  let ages : Array[Double] = Array(capacity=0)
  let n = dendropy_gamma_collect(clade, ages)
  if n < 3 || ages.length() == 0 {
    return 0.0
  }
  let sorted = dendropy_sort_doubles_desc(ages)
  let g : Array[Double] = Array(capacity=0)
  let mut older = sorted[0]
  for i = 1; i < sorted.length(); i = i + 1 {
    g.push(older - sorted[i])
    older = sorted[i]
  }
  g.push(older)
  let mut t = 0.0
  let mut accum = 0.0
  for i = 2; i < n; i = i + 1 {
    let idx = i - 2
    t = t + i.to_double() * g[idx]
    accum = accum + t
  }
  t = t + n.to_double() * g[n - 2]
  let nmt = (n - 2).to_double()
  let numerator = accum / nmt - t / 2.0
  let c = t * (1.0 / (12.0 * nmt)).sqrt()
  if c == 0.0 {
    0.0
  } else {
    numerator / c
  }
}