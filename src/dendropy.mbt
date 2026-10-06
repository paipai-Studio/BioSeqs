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
  let sb = StringBuilder()
  let mut first = true
  for n in names {
    if !(first) {
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

// ============================================================================
// DendroPy treemeasure / treecompare 补齐：treeness / euclidean / fusco_cronk
// ============================================================================

///|
/// Canonical split key for a clade: sorted, comma-joined terminal names.
fn dendropy_split_key(clade : Clade) -> String {
  let names : Array[String] = Array(capacity=0)
  for t in clade.get_terminals() {
    match t.name {
      Some(n) => names.push(n)
      None => ()
    }
  }
  let sorted = dendropy_sort_strings(names)
  let sb = StringBuilder()
  for n in sorted {
    sb.write_string(n)
    sb.write_char(',')
  }
  sb.to_string()
}

///|
/// Treeness: proportion of total tree length that is internal branch length
/// (`dendropy.calculate.treemeasure.treeness`).
pub fn dendropy_treeness(clade : Clade) -> Double {
  let mut internal = 0.0
  let mut external = 0.0
  for leaf in clade.get_terminals() {
    match leaf.branch_length {
      Some(l) => external = external + l
      None => ()
    }
  }
  for nonterm in clade.get_nonterminals() {
    match nonterm.branch_length {
      Some(l) => internal = internal + l
      None => ()
    }
  }
  let total = internal + external
  if total == 0.0 {
    0.0
  } else {
    internal / total
  }
}

///|
fn dendropy_collect_splits(clade : Clade, out : Map[String, Double]) -> Unit {
  if clade.clades.length() == 0 {
    return
  }
  let key = dendropy_split_key(clade)
  let length = match clade.branch_length {
    Some(l) => l
    None => 0.0
  }
  out.set(key, length)
  for c in clade.clades {
    dendropy_collect_splits(c, out)
  }
}

///|
/// Euclidean branch-length distance between two trees (Felsenstein 2004 /
/// `dendropy.calculate.treecompare.euclidean_distance`): sqrt of the sum of
/// squared branch-length differences over the union of bipartitions (a split
/// present in one tree only is treated as length 0 in the other).
pub fn dendropy_euclidean_distance(clade1 : Clade, clade2 : Clade) -> Double {
  let s1 : Map[String, Double] = Map([], capacity=0)
  let s2 : Map[String, Double] = Map([], capacity=0)
  dendropy_collect_splits(clade1, s1)
  dendropy_collect_splits(clade2, s2)
  let mut sum_sq = 0.0
  for (key, l1) in s1 {
    let l2 = match s2.get(key) {
      Some(l) => l
      None => 0.0
    }
    let d = l1 - l2
    sum_sq = sum_sq + d * d
  }
  for (key, l2) in s2 {
    if !s1.contains(key) {
      sum_sq = sum_sq + l2 * l2
    }
  }
  sum_sq.sqrt()
}

///|
/// Fusco-Cronk balance statistic: mean over bifurcating interior nodes of
/// `Im = |m_L - m_R| / (m_L + m_R - 2)`, where m_L/m_R are the child clade's
/// terminal counts (Fusco & Cronk 1995 tree-shape imbalance).
pub fn dendropy_fusco_cronk(clade : Clade) -> Double {
  let mut sum_im = 0.0
  let mut count = 0
  for nonterm in clade.get_nonterminals() {
    if nonterm.clades.length() == 2 {
      let ml = nonterm.clades[0].count_terminals()
      let mr = nonterm.clades[1].count_terminals()
      let diff = if ml > mr { ml - mr } else { mr - ml }
      let denom = ml + mr - 2
      if denom > 0 {
        sum_im = sum_im + diff.to_double() / denom.to_double()
        count = count + 1
      }
    }
  }
  if count == 0 {
    0.0
  } else {
    sum_im / count.to_double()
  }
}


///|
// =========================================================================
// s-shape statistic (Bortolussi et al. 2005 / Blum & François 2006)
//   σ(T) = Σ_{u ∈ internal} log(n_u - 1)
//   where n_u = number of leaves in the pending subtree rooted at u.
// Source: apTreeshape/R/shape.statistic.R, treebalance::sShapeI
// =========================================================================

///|
/// Returns (total_leaves, sum_of_log(n_u - 1)) for the subtree.
fn dendropy_shape_stat_rec(clade : Clade) -> (Int, Double) {
  if clade.is_terminal() {
    return (1, 0.0)
  }
  let mut total_leaves = 0
  let mut sigma = 0.0
  for ch in clade.clades {
    let (leaves, s) = dendropy_shape_stat_rec(ch)
    total_leaves = total_leaves + leaves
    sigma = sigma + s
  }
  // Add log(n_u - 1) for this internal node (if n_u > 1)
  if total_leaves > 1 {
    sigma = sigma + @math.ln(total_leaves.to_double() - 1.0)
  }
  (total_leaves, sigma)
}

///|
/// s-shape statistic (raw): Σ log(n_u - 1) over all internal nodes u,
/// where n_u is the number of leaves in the pending subtree at u.
/// Uses natural logarithm. For a single leaf returns 0.
pub fn dendropy_shape_statistic(clade : Clade) -> Double {
  dendropy_shape_stat_rec(clade).1
}

///|
/// s-shape statistic with Yule-model normalization.
/// σ'_yule = (σ - 1.204*(n+1) + log(n) + 2) / sqrt(0.168*(n+1) - 0.71)
/// where n = number of internal nodes (= leaves - 1 for binary trees).
/// Under the Yule model this is approximately Standard Normal.
pub fn dendropy_shape_statistic_yule(clade : Clade) -> Double {
  let n_leaves = dendropy_count_leaves(clade)
  let n = n_leaves - 1 // internal nodes for binary tree
  if n < 1 {
    return 0.0
  }
  let sigma = dendropy_shape_statistic(clade)
  let numerator = sigma - 1.204 * (n.to_double() + 1.0) + @math.ln(n.to_double()) + 2.0
  let denom_sq = 0.168 * (n.to_double() + 1.0) - 0.71
  if denom_sq <= 0.0 {
    return 0.0
  }
  numerator / denom_sq.sqrt()
}

///|
/// s-shape statistic with PDA-model normalization.
/// σ'_pda = (σ - 2.03*(n+1) + 3.545*sqrt(n)) / sqrt(1.570*n*log(n) - 5.674*n + 3.602*sqrt(n) + 14.915)
/// where n = number of internal nodes (= leaves - 1 for binary trees).
/// Under the PDA (uniform) model this is approximately Standard Normal.
pub fn dendropy_shape_statistic_pda(clade : Clade) -> Double {
  let n_leaves = dendropy_count_leaves(clade)
  let n = n_leaves - 1
  if n < 1 {
    return 0.0
  }
  let nd = n.to_double()
  let sigma = dendropy_shape_statistic(clade)
  let sqrt_n = nd.sqrt()
  let numerator = sigma - 2.03 * (nd + 1.0) + 3.545 * sqrt_n
  let denom_sq = 1.570 * nd * @math.ln(nd) - 5.674 * nd + 3.602 * sqrt_n + 14.915
  if denom_sq <= 0.0 {
    return 0.0
  }
  numerator / denom_sq.sqrt()
}

///|
// =========================================================================
// Normalized Colless index (Yule / PDA)
//   Source: DendroPy colless_tree_imbalance, apTreeshape colless.R
// =========================================================================

/// Euler–Mascheroni constant γ ≈ 0.5772156649
let dendropy_euler_gamma : Double = 0.5772156649015329

///|
/// Colless index with Yule-model normalization.
/// Ic' = (Ic - n*log(n) - n*(γ - 1 - log(2))) / n
/// where n = number of leaves. Under Yule, Ic' ≈ Standard Normal.
pub fn dendropy_colless_yule(clade : Clade) -> Double {
  let n = dendropy_count_leaves(clade)
  if n < 2 {
    return 0.0
  }
  let nd = n.to_double()
  let ic = dendropy_colless(clade).to_double()
  (ic - nd * @math.ln(nd) - nd * (dendropy_euler_gamma - 1.0 - @math.ln(2.0))) / nd
}

///|
/// Colless index with PDA-model normalization.
/// Ic_pda = Ic / n^(3/2)
pub fn dendropy_colless_pda(clade : Clade) -> Double {
  let n = dendropy_count_leaves(clade)
  if n < 2 {
    return 0.0
  }
  let nd = n.to_double()
  let ic = dendropy_colless(clade).to_double()
  ic / nd.sqrt() * nd
}

///|
// =========================================================================
// Normalized Sackin index (Yule / PDA)
//   Source: DendroPy sackin_index, apTreeshape sackin.R
// =========================================================================

///|
/// Sackin index with Yule-model normalization.
/// Is' = (S - 2n * Σ_{j=2}^{n} 1/j) / n
/// where n = number of leaves. Under Yule, Is' ≈ Standard Normal.
pub fn dendropy_sackin_yule(clade : Clade) -> Double {
  let n = dendropy_count_leaves(clade)
  if n < 2 {
    return 0.0
  }
  let nd = n.to_double()
  let s = dendropy_sackin(clade).to_double()
  // Compute harmonic sum H_n - 1 = Σ_{j=2}^{n} 1/j
  let mut h = 0.0
  for j in 2..=n {
    h = h + 1.0 / j.to_double()
  }
  (s - 2.0 * nd * h) / nd
}

///|
/// Sackin index with PDA-model normalization.
/// Is_pda = S / n^(3/2)
pub fn dendropy_sackin_pda(clade : Clade) -> Double {
  let n = dendropy_count_leaves(clade)
  if n < 2 {
    return 0.0
  }
  let nd = n.to_double()
  let s = dendropy_sackin(clade).to_double()
  s / nd.sqrt() * nd
}


///|
// =========================================================================
// Colless-like index (Mir, Rotger & Rossello 2018)
//   Generalizes Colless to multifurcating trees.
//   CL(T) = Σ_{v ∈ internal} D(fs(T_{v1}), ..., fs(T_{vk}))
//   where fs(T) = Σ_{v ∈ V(T)} f(deg+(v)) is the f.size,
//   and D is a dissimilarity (mdm, var, or sd).
// Source: treebalance::collesslikeI, Mir et al. PLOSONE 2018
// =========================================================================

///|
/// Weight function f for f.size computation.
pub(all) enum CollessLikeFSize {
  FExp // f(x) = e^x
  FLn  // f(x) = ln(x + e)
}

///|
/// Dissimilarity D for balance value computation.
pub(all) enum CollessLikeDissim {
  DMdm // mean deviation from median: Σ|x_i - median| / n
  DVar // sample variance: Σ(x_i - mean)^2 / (n-1)
  DSd  // sample standard deviation: sqrt(var)
}

///|
/// Evaluate the weight function f at integer x.
fn cl_fsize_eval(f : CollessLikeFSize, x : Int) -> Double {
  match f {
    FExp => @math.exp(x.to_double())
    FLn => @math.ln(x.to_double() + @math.exp(1.0))
  }
}

///|
/// Evaluate the dissimilarity D on a vector of doubles.
fn cl_dissim_eval(d : CollessLikeDissim, xs : Array[Double]) -> Double {
  let n = xs.length()
  if n <= 1 {
    return 0.0
  }
  // Compute mean
  let mut sum = 0.0
  for x in xs {
    sum = sum + x
  }
  let mean = sum / n.to_double()
  match d {
    DMdm => {
      // Median: sort copy, pick middle (or average of two middles)
      let sorted = Array(capacity=n)
      for x in xs {
        sorted.push(x)
      }
      sorted.sort_by(fn(a, b) { compare(a, b) })
      let median = if n % 2 == 1 {
        sorted[n / 2]
      } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
      }
      let mut dev = 0.0
      for x in xs {
        dev = dev + (x - median).abs()
      }
      dev / n.to_double()
    }
    DVar => {
      let mut ss = 0.0
      for x in xs {
        ss = ss + (x - mean) * (x - mean)
      }
      ss / (n - 1).to_double()
    }
    DSd => {
      let mut ss = 0.0
      for x in xs {
        ss = ss + (x - mean) * (x - mean)
      }
      (ss / (n - 1).to_double()).sqrt()
    }
  }
}

///|
/// Compute f.size of a subtree: fs(T) = Σ_{v ∈ V(T)} f(deg+(v)).
/// Returns (f_size, ) for the subtree rooted at clade.
fn cl_compute_fsize(
  clade : Clade,
  f : CollessLikeFSize,
) -> Double {
  let outdeg = clade.clades.length()
  let mut fs = cl_fsize_eval(f, outdeg)
  for ch in clade.clades {
    fs = fs + cl_compute_fsize(ch, f)
  }
  fs
}

///|
/// Colless-like index CL(T) for a rooted (possibly multifurcating) tree.
///
/// CL(T) = Σ_{v ∈ internal} D(fs(T_{v1}), ..., fs(T_{vk}))
///
/// where fs is the f.size with weight function `f`, and D is the
/// dissimilarity. For binary trees with f=exp and D=mdm, this reduces
/// to a variant of the Colless index.
///
/// Reference: Mir, Rotger & Rossello. PLOSONE 13(9):e0203401, 2018.
pub fn dendropy_colless_like(
  clade : Clade,
  f : CollessLikeFSize,
  dissim : CollessLikeDissim,
) -> Double {
  cl_colless_like_rec(clade, f, dissim)
}

///|
fn cl_colless_like_rec(
  clade : Clade,
  f : CollessLikeFSize,
  dissim : CollessLikeDissim,
) -> Double {
  if clade.is_terminal() {
    return 0.0
  }
  // Compute balance value at this node
  let child_fsizes : Array[Double] = Array(capacity=0)
  for ch in clade.clades {
    child_fsizes.push(cl_compute_fsize(ch, f))
  }
  let bal = cl_dissim_eval(dissim, child_fsizes)
  // Recurse into children
  let mut total = bal
  for ch in clade.clades {
    total = total + cl_colless_like_rec(ch, f, dissim)
  }
  total
}


///|
// =========================================================================
// Cherry count
//   A cherry is a pair of sister leaf nodes (an internal node with
//   exactly two children, both leaves). This is a basic tree shape
//   statistic used in ete3 and DendroPy.
// =========================================================================

///|
/// Count the number of cherries in a rooted tree.
/// A cherry is an internal node with exactly two children that are both leaves.
pub fn dendropy_cherry_count(clade : Clade) -> Int {
  if clade.is_terminal() {
    return 0
  }
  let mut count = 0
  // Check if this node is a cherry
  if clade.clades.length() == 2 {
    if clade.clades[0].is_terminal() && clade.clades[1].is_terminal() {
      count = 1
    }
  }
  // Recurse into children
  for ch in clade.clades {
    count = count + dendropy_cherry_count(ch)
  }
  count
}


///|
// =========================================================================
// B2 index (Shao & Sokal 1990)
//   B2(T) = -Σ_{x ∈ leaves} p_x * log(p_x)
//   where p_x = Π_{v ∈ anc(x)} 1/|child(v)| is the probability of
//   reaching leaf x under equiprobable branching.
//   Default log base 2 (standard for binary trees).
// Reference: Shao & Sokal, Systematic Zoology 39(3):266, 1990.
// =========================================================================

///|
/// Compute B2 index with log base 2 (standard for binary trees).
pub fn dendropy_b2(clade : Clade) -> Double {
  dendropy_b2_base(clade, 2.0)
}

///|
/// Compute B2 index with a custom log base.
pub fn dendropy_b2_base(clade : Clade, log_base : Double) -> Double {
  // Collect (probability) for each leaf
  let probs : Array[Double] = Array(capacity=0)
  b2_collect_probs(clade, 1.0, probs)
  if probs.length() <= 1 {
    return 0.0
  }
  let log_b = @math.ln(log_base) // ln(base) for change of base
  let mut b2 = 0.0
  for p in probs {
    if p > 0.0 {
      b2 = b2 - p * (@math.ln(p) / log_b)
    }
  }
  b2
}

///|
/// Recursively collect leaf probabilities.
fn b2_collect_probs(
  clade : Clade,
  prob : Double,
  probs : Array[Double],
) -> Unit {
  if clade.is_terminal() {
    probs.push(prob)
    return
  }
  let n_children = clade.clades.length()
  let child_prob = prob / n_children.to_double()
  for ch in clade.clades {
    b2_collect_probs(ch, child_prob, probs)
  }
}

///|
/// Collect only non-trivial internal bipartitions (both sides ≥ 2 taxa).
/// Used by Kuhner-Felsenstein distance; skips the root and any split where
/// one side has a single taxon.
fn dendropy_collect_internal_splits(
  clade : Clade,
  total_taxa : Int,
  is_root : Bool,
  out : Map[String, Double]
) -> Unit {
  if clade.clades.length() == 0 {
    return
  }
  if not(is_root) {
    let k = clade.count_terminals()
    // Non-trivial: 2 ≤ k ≤ total-2
    if k >= 2 && k <= total_taxa - 2 {
      let key = dendropy_split_key(clade)
      let length = match clade.branch_length {
        Some(l) => l
        None => 0.0
      }
      out.set(key, length)
    }
  }
  for c in clade.clades {
    dendropy_collect_internal_splits(c, total_taxa, false, out)
  }
}

///|
/// Kuhner-Felsenstein branch-score distance (Kuhner & Felsenstein 1994):
/// sqrt of the sum of squared branch-length differences over the union of
/// **non-trivial internal bipartitions** only (both sides ≥ 2 taxa).
/// Unlike `dendropy_euclidean_distance`, pendant branches and the root are
/// excluded — only splits that affect topology contribute.
pub fn dendropy_kuhner_felsenstein_distance(
  clade1 : Clade,
  clade2 : Clade
) -> Double {
  let n1 = clade1.count_terminals()
  let n2 = clade2.count_terminals()
  // Trees must have the same number of taxa to be comparable
  if n1 != n2 {
    return -1.0
  }
  let s1 : Map[String, Double] = Map([], capacity=0)
  let s2 : Map[String, Double] = Map([], capacity=0)
  dendropy_collect_internal_splits(clade1, n1, true, s1)
  dendropy_collect_internal_splits(clade2, n2, true, s2)
  let mut sum_sq = 0.0
  for (key, l1) in s1 {
    let l2 = match s2.get(key) {
      Some(l) => l
      None => 0.0
    }
    let d = l1 - l2
    sum_sq = sum_sq + d * d
  }
  for (key, l2) in s2 {
    if not(s1.contains(key)) {
      sum_sq = sum_sq + l2 * l2
    }
  }
  sum_sq.sqrt()
}
