///|
/// Kingman coalescent simulation and tree visualization.
/// Port of DendroPy's coalescent simulation and tree drawing utilities.

// =========================================================================
// Kingman Coalescent Simulation
//   Standard n-coalescent: start with n lineages, each pair coalesces
//   at rate 1. With k lineages, total rate = k*(k-1)/2.
//   Waiting time ~ Exp(k*(k-1)/2).
//   For diploid with effective population size Ne: rate = k*(k-1)/(2*Ne).
// =========================================================================

///|
/// Result of a coalescent simulation.
pub(all) struct CoalescentResult {
  tree : Clade                // the simulated coalescent tree
  coalescent_times : Array[Double]  // times of each coalescent event
  total_depth : Double       // time from present to MRCA
  n_lineages : Int            // number of starting lineages
}

///|
/// Simulate a Kingman coalescent tree with n lineages.
/// Uses the standard coalescent (Ne=1, haploid).
/// Branch lengths are in coalescent time units.
pub fn phylo_simulate_coalescent(n : Int, rng : Rng) -> CoalescentResult {
  phylo_simulate_coalescent_ne(n, 1.0, rng)
}

///|
/// Simulate a Kingman coalescent tree with effective population size Ne.
/// For diploid: waiting time with k lineages ~ Exp(k*(k-1)/(2*Ne)).
pub fn phylo_simulate_coalescent_ne(
  n : Int,
  ne : Double,
  rng : Rng
) -> CoalescentResult {
  if n <= 0 {
    return CoalescentResult::{
      tree: Clade::new(name=Some("t1")),
      coalescent_times: [],
      total_depth: 0.0,
      n_lineages: 0,
    }
  }
  if ne <= 0.0 {
    return CoalescentResult::{
      tree: Clade::new(name=Some("t1")),
      coalescent_times: [],
      total_depth: 0.0,
      n_lineages: n,
    }
  }
  if n == 1 {
    return CoalescentResult::{
      tree: Clade::new(name=Some("t1")),
      coalescent_times: [],
      total_depth: 0.0,
      n_lineages: 1,
    }
  }
  // Active lineages: each is a Clade with a recorded time
  let active : Array[CoalescentLineage] = Array(capacity=n)
  for i in 0..<n {
    active.push(CoalescentLineage::{
      clade: Clade::new(name=Some("t" + (i + 1).to_string())),
      time: 0.0,
    })
  }
  let coalescent_times : Array[Double] = Array(capacity=n - 1)
  let mut current_time = 0.0
  let mut k = n
  while k > 1 {
    // Waiting time ~ Exp(k*(k-1)/(2*Ne))
    let rate = k.to_double() * (k - 1).to_double() / (2.0 * ne)
    let u = rng.next_double()
    // Avoid log(0)
    let u_safe = if u < 1.0e-10 { 1.0e-10 } else { u }
    let wait_time = -@math.ln(u_safe) / rate
    current_time = current_time + wait_time
    coalescent_times.push(current_time)
    // Pick two random distinct lineages
    let i1 = rng.randint(k)
    let mut i2 = rng.randint(k)
    if i2 == i1 {
      i2 = (i1 + 1 + rng.randint(k - 1)) % k
    }
    let idx1 = if i1 < i2 { i1 } else { i2 }
    let idx2 = if i1 < i2 { i2 } else { i1 }
    // Get the two lineages
    let lin1 = active[idx1]
    let lin2 = active[idx2]
    // Set branch lengths
    let bl1 = current_time - lin1.time
    let bl2 = current_time - lin2.time
    let child1 = set_branch_length(lin1.clade, bl1)
    let child2 = set_branch_length(lin2.clade, bl2)
    // Create new internal node
    let new_node = Clade::new(
      name=Some("c" + coalescent_times.length().to_string()),
      clades=[child1, child2],
    )
    // Replace idx1 with new node, remove idx2
    active[idx1] = CoalescentLineage::{ clade: new_node, time: current_time }
    // Remove idx2 by swapping with last and popping
    let last = k - 1
    if idx2 != last {
      active[idx2] = active[last]
    }
    ignore(active.pop())
    k = k - 1
  }
  let root = active[0].clade
  CoalescentResult::{
    tree: root,
    coalescent_times,
    total_depth: current_time,
    n_lineages: n,
  }
}

///|
/// Internal helper: a lineage in the coalescent process.
struct CoalescentLineage {
  clade : Clade
  time : Double
}

///|
/// Set branch length on a clade, returning a new clade.
fn set_branch_length(clade : Clade, bl : Double) -> Clade {
  Clade::new(
    branch_length=Some(bl),
    name=clade.name,
    confidence=clade.confidence,
    clades=clade.clades,
  )
}

///|
/// Simulate multiple coalescent replicates and return the mean tree depth.
pub fn phylo_coalescent_mean_depth(
  n : Int,
  ne : Double,
  replicates : Int,
  rng : Rng
) -> Double {
  if replicates <= 0 {
    return 0.0
  }
  let mut total = 0.0
  for _ in 0..<replicates {
    let result = phylo_simulate_coalescent_ne(n, ne, rng)
    total = total + result.total_depth
  }
  total / replicates.to_double()
}

///|
/// Expected total depth of Kingman coalescent with n lineages.
/// E[T_total] = 2 * (1 - 1/n) for haploid (Ne=1).
pub fn phylo_coalescent_expected_depth(n : Int) -> Double {
  if n <= 1 {
    return 0.0
  }
  2.0 * (1.0 - 1.0 / n.to_double())
}

///|
/// Expected number of segregating sites (Watterson's estimator).
/// E[S] = theta * H_{n-1}, where H_{n-1} = sum_{i=1}^{n-1} 1/i.
pub fn phylo_watterson_expected(n : Int, theta : Double) -> Double {
  if n <= 1 {
    return 0.0
  }
  let mut harmonic = 0.0
  for i in 1..<n {
    harmonic = harmonic + 1.0 / i.to_double()
  }
  theta * harmonic
}

// =========================================================================
// Tree Visualization
//   Enhanced ASCII tree drawing with branch lengths, support values,
//   rectangular and radial layouts.
// =========================================================================

///|
/// Draw a tree in rectangular layout with branch lengths.
/// Example output:
///   t1 ───────────┐
///                  ├──────────┐
///   t2 ───────────┘           │
///                              ├────────── root
///   t3 ───────────────────────┘
pub fn phylo_draw_rectangular(tree : Clade) -> String {
  let lines : Array[String] = Array(capacity=0)
  phylo_rect_helper(tree, "", 0, lines)
  let sb = StringBuilder()
  for i in 0..<lines.length() {
    sb.write_string(lines[i])
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn phylo_rect_helper(
  node : Clade,
  prefix : String,
  depth : Int,
  lines : Array[String]
) -> Unit {
  let bl_str = match node.branch_length {
    Some(bl) => {
      let bl_rounded = (bl * 100.0).to_int().to_double() / 100.0
      " " + bl_rounded.to_string()
    }
    None => ""
  }
  let conf_str = match node.confidence {
    Some(c) => {
      let c_int = (c + 0.5).to_int()
      "/" + c_int.to_string()
    }
    None => ""
  }
  if node.is_terminal() {
    let name = match node.name {
      Some(n) => n
      None => ""
    }
    lines.push(prefix + "── " + name + bl_str)
  } else {
    let name = match node.name {
      Some(n) => n
      None => ""
    }
    // Process children first
    let child_count = node.clades.length()
    for i in 0..<child_count {
      let child = node.clades[i]
      let is_last = i == child_count - 1
      let connector = if is_last { "└" } else { "├" }
      let child_prefix = if is_last { prefix + "   " } else { prefix + "│  " }
      if i == 0 {
        // First child: add branch info on the connector line
        if name != "" || conf_str != "" {
          lines.push(prefix + connector + bl_str + conf_str + " " + name)
        }
      }
      phylo_rect_helper(child, child_prefix, depth + 1, lines)
    }
  }
}

///|
/// Draw a tree in a compact ASCII format with branch lengths.
/// This is an enhanced version of Tree::draw_ascii that includes
/// branch lengths and confidence values.
pub fn phylo_draw_ascii_enhanced(tree : Clade) -> String {
  let sb = StringBuilder()
  phylo_draw_ascii_helper(tree, "", true, sb)
  sb.to_string()
}

///|
fn phylo_draw_ascii_helper(
  node : Clade,
  prefix : String,
  is_last : Bool,
  buf : StringBuilder
) -> Unit {
  let name = match node.name {
    Some(n) => n
    None => ""
  }
  let bl_str = match node.branch_length {
    Some(bl) => {
      let bl_rounded = (bl * 1000.0).to_int().to_double() / 1000.0
      ":" + bl_rounded.to_string()
    }
    None => ""
  }
  let conf_str = match node.confidence {
    Some(c) => {
      let c_int = (c + 0.5).to_int()
      " (" + c_int.to_string() + ")"
    }
    None => ""
  }
  buf.write_string(prefix)
  if is_last {
    buf.write_string("└── ")
  } else {
    buf.write_string("├── ")
  }
  buf.write_string(name)
  buf.write_string(bl_str)
  buf.write_string(conf_str)
  buf.write_char('\n')
  let child_count = node.clades.length()
  for i in 0..<child_count {
    let child = node.clades[i]
    let child_is_last = i == child_count - 1
    let child_prefix = if is_last { prefix + "    " } else { prefix + "│   " }
    phylo_draw_ascii_helper(child, child_prefix, child_is_last, buf)
  }
}

///|
/// Draw a tree as a cladogram (equal branch lengths, no proportional).
/// Useful for displaying tree topology without branch length distortion.
pub fn phylo_draw_cladogram(tree : Clade) -> String {
  let sb = StringBuilder()
  phylo_draw_cladogram_helper(tree, "", true, sb)
  sb.to_string()
}

///|
fn phylo_draw_cladogram_helper(
  node : Clade,
  prefix : String,
  is_last : Bool,
  buf : StringBuilder
) -> Unit {
  let name = match node.name {
    Some(n) => n
    None => ""
  }
  buf.write_string(prefix)
  if is_last {
    buf.write_string("└── ")
  } else {
    buf.write_string("├── ")
  }
  buf.write_string(name)
  buf.write_char('\n')
  let child_count = node.clades.length()
  for i in 0..<child_count {
    let child = node.clades[i]
    let child_is_last = i == child_count - 1
    let child_prefix = if is_last { prefix + "    " } else { prefix + "│   " }
    phylo_draw_cladogram_helper(child, child_prefix, child_is_last, buf)
  }
}

///|
/// Draw a tree with a horizontal bar chart showing branch lengths.
/// Each leaf gets a horizontal bar proportional to its total depth.
pub fn phylo_draw_barplot(tree : Clade, max_bar_width : Int) -> String {
  // Collect leaves with their depths
  let leaves : Array[(String, Double)] = Array(capacity=0)
  phylo_collect_leaf_depths(tree, 0.0, leaves)
  // Find max depth
  let mut max_depth = 1.0
  for (_, d) in leaves {
    if d > max_depth {
      max_depth = d
    }
  }
  let sb = StringBuilder()
  for (name, depth) in leaves {
    let bar_len = (depth / max_depth * max_bar_width.to_double()).to_int()
    sb.write_string(name)
    // Pad name to align bars
    let name_len = name.length()
    let pad = if name_len < 20 { 20 - name_len } else { 1 }
    for _ in 0..<pad {
      sb.write_char(' ')
    }
    for _ in 0..<bar_len {
      sb.write_char('█')
    }
    sb.write_string(" ")
    let d_rounded = (depth * 100.0).to_int().to_double() / 100.0
    sb.write_string(d_rounded.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn phylo_collect_leaf_depths(
  node : Clade,
  current_depth : Double,
  out : Array[(String, Double)]
) -> Unit {
  let bl = match node.branch_length {
    Some(b) => b
    None => 0.0
  }
  let new_depth = current_depth + bl
  if node.is_terminal() {
    let name = match node.name {
      Some(n) => n
      None => ""
    }
    out.push((name, new_depth))
  } else {
    for ch in node.clades {
      phylo_collect_leaf_depths(ch, new_depth, out)
    }
  }
}

///|
/// Format a coalescent result as a summary string.
pub fn phylo_coalescent_summary(result : CoalescentResult) -> String {
  let sb = StringBuilder()
  sb.write_string("Kingman Coalescent Simulation:\n")
  sb.write_string("  Lineages: ")
  sb.write_string(result.n_lineages.to_string())
  sb.write_string("\n")
  sb.write_string("  Coalescent events: ")
  sb.write_string(result.coalescent_times.length().to_string())
  sb.write_string("\n")
  sb.write_string("  Total depth (T_MRCA): ")
  let d = (result.total_depth * 10000.0).to_int().to_double() / 10000.0
  sb.write_string(d.to_string())
  sb.write_string("\n")
  sb.write_string("  Expected depth: ")
  let exp = phylo_coalescent_expected_depth(result.n_lineages)
  let exp_r = (exp * 10000.0).to_int().to_double() / 10000.0
  sb.write_string(exp_r.to_string())
  sb.write_string("\n")
  sb.write_string("  Event times:\n")
  for i in 0..<result.coalescent_times.length() {
    sb.write_string("    k=")
    sb.write_string((result.n_lineages - i).to_string())
    sb.write_string(" → k=")
    sb.write_string((result.n_lineages - i - 1).to_string())
    sb.write_string(": t=")
    let t = (result.coalescent_times[i] * 10000.0).to_int().to_double() /
      10000.0
    sb.write_string(t.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

// =========================================================================
// Site Frequency Spectrum (SFS)
//   Expected SFS under standard neutral model: E[ξ_i] = θ / i
//   Folded SFS: E[η_i] = θ/i + θ/(n-i)  for i < n/2
//   Plus classic test statistics: Watterson θ, Tajima D, Fu & Li D/F.
// =========================================================================

///|
/// Expected unfolded SFS: array of E[ξ_i] for i = 1..n-1.
/// Under the standard neutral coalescent, E[ξ_i] = θ / i.
pub fn phylo_expected_sfs(n : Int, theta : Double) -> Array[Double] {
  if n <= 1 {
    return []
  }
  let result : Array[Double] = Array(capacity=n - 1)
  for i in 1..<n {
    result.push(theta / i.to_double())
  }
  result
}

///|
/// Expected folded SFS: array of E[η_i] for i = 1..floor(n/2).
/// η_i = ξ_i + ξ_{n-i} for i < n/2, η_i = ξ_i for i = n/2.
pub fn phylo_expected_sfs_folded(n : Int, theta : Double) -> Array[Double] {
  if n <= 1 {
    return []
  }
  let half = n / 2
  let result : Array[Double] = Array(capacity=half)
  for i in 1..<(half + 1) {
    let val = if i == n - i {
      // i == n/2 (only when n is even and i = n/2)
      theta / i.to_double()
    } else {
      theta / i.to_double() + theta / (n - i).to_double()
    }
    result.push(val)
  }
  result
}

///|
/// Observe SFS from aligned sequences.
/// `sequences` is an array of equal-length strings (DNA).
/// Ancestral state is taken as the first sequence's base at each position.
/// Returns array of counts for derived allele frequency i = 1..n-1.
pub fn phylo_observed_sfs(sequences : Array[String]) -> Array[Int] {
  let n = sequences.length()
  if n <= 1 {
    return []
  }
  let seq_len = sequences[0].length()
  let sfs : Array[Int] = Array(capacity=n - 1)
  for _ in 1..<n {
    sfs.push(0)
  }
  for pos in 0..<seq_len {
    let ancestral = sequences[0].char_to_at(pos)
    let mut derived_count = 0
    for i in 1..<n {
      if sequences[i].char_to_at(pos) != ancestral {
        derived_count = derived_count + 1
      }
    }
    if derived_count > 0 && derived_count < n {
      sfs[derived_count - 1] = sfs[derived_count - 1] + 1
    }
  }
  sfs
}

///|
/// Observe folded SFS from aligned sequences.
/// For each segregating site, the minor allele count is min(derived, n-derived).
/// Returns array of counts for minor allele frequency i = 1..floor(n/2).
pub fn phylo_observed_sfs_folded(sequences : Array[String]) -> Array[Int] {
  let n = sequences.length()
  if n <= 1 {
    return []
  }
  let seq_len = sequences[0].length()
  let half = n / 2
  let sfs : Array[Int] = Array(capacity=half)
  for _ in 1..<(half + 1) {
    sfs.push(0)
  }
  for pos in 0..<seq_len {
    // Count allele frequencies
    let counts : Map[Int, Int] = Map([], capacity=4)
    for i in 0..<n {
      let c = sequences[i].char_to_at(pos)
      match counts.get(c) {
        Some(v) => counts.set(c, v + 1)
        None => counts.set(c, 1)
      }
    }
    if counts.size() == 2 {
      // Biallelic: get minor count
      let mut minor = n
      for (_, v) in counts {
        if v < minor {
          minor = v
        }
      }
      if minor >= 1 && minor <= half {
        sfs[minor - 1] = sfs[minor - 1] + 1
      }
    }
  }
  sfs
}

///|
/// Helper: get character code at position.
fn String::char_to_at(self : String, pos : Int) -> Int {
  self.unsafe_get(pos).to_int()
}

///|
/// Count segregating sites from aligned sequences.
pub fn phylo_segregating_sites(sequences : Array[String]) -> Int {
  let n = sequences.length()
  if n <= 1 {
    return 0
  }
  let seq_len = sequences[0].length()
  let mut count = 0
  for pos in 0..<seq_len {
    let first = sequences[0].char_to_at(pos)
    let mut seg = false
    for i in 1..<n {
      if sequences[i].char_to_at(pos) != first {
        seg = true
      }
    }
    if seg {
      count = count + 1
    }
  }
  count
}

///|
/// Count singletons (derived allele frequency = 1) from aligned sequences.
pub fn phylo_count_singletons(sequences : Array[String]) -> Int {
  let n = sequences.length()
  if n <= 1 {
    return 0
  }
  let seq_len = sequences[0].length()
  let mut count = 0
  for pos in 0..<seq_len {
    let ancestral = sequences[0].char_to_at(pos)
    let mut derived_count = 0
    for i in 1..<n {
      if sequences[i].char_to_at(pos) != ancestral {
        derived_count = derived_count + 1
      }
    }
    if derived_count == 1 {
      count = count + 1
    }
  }
  count
}

///|
/// Watterson's estimator: θ_W = S / H_{n-1}.
pub fn phylo_watterson_theta(s : Int, n : Int) -> Double {
  if n <= 1 || s == 0 {
    return 0.0
  }
  let mut harmonic = 0.0
  for i in 1..<n {
    harmonic = harmonic + 1.0 / i.to_double()
  }
  s.to_double() / harmonic
}

///|
/// Nucleotide diversity π: average pairwise differences.
/// π = (2 / (n*(n-1))) * Σ_{i<j} d(i,j)
pub fn phylo_pi(sequences : Array[String]) -> Double {
  let n = sequences.length()
  if n <= 1 {
    return 0.0
  }
  let seq_len = sequences[0].length()
  let mut total_diff = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      for pos in 0..<seq_len {
        if sequences[i].char_to_at(pos) != sequences[j].char_to_at(pos) {
          total_diff = total_diff + 1
        }
      }
    }
  }
  let num_pairs = n.to_double() * (n - 1).to_double() / 2.0
  total_diff.to_double() / num_pairs
}

///|
/// Tajima's D test statistic.
/// D = (π - S/a1) / sqrt(e1*S + e2*S*(S-1))
pub fn phylo_tajima_d(s : Int, n : Int, pi : Double) -> Double {
  if n <= 1 || s == 0 {
    return 0.0
  }
  let a1 = sfs_harmonic(n - 1)
  let a2 = sfs_harmonic_sq(n - 1)
  let b1 = (n.to_double() + 1.0) / (3.0 * (n.to_double() - 1.0))
  let b2 = 2.0 * (
    n.to_double() * n.to_double() + n.to_double() + 3.0
  ) / (9.0 * n.to_double() * (n.to_double() - 1.0))
  let c1 = b1 - 1.0 / a1
  let c2 = b2 - (n.to_double() + 2.0) / (a1 * n.to_double()) + a2 / (a1 * a1)
  let e1 = c1 / a1
  let e2 = c2 / (a1 * a1 + a2)
  let s_d = s.to_double()
  let numerator = pi - s_d / a1
  let denominator = sfs_sqrt(e1 * s_d + e2 * s_d * (s_d - 1.0))
  if denominator == 0.0 {
    return 0.0
  }
  numerator / denominator
}

///|
/// Fu and Li's D test statistic.
/// D_FL = (S - a_n*η) / sqrt(u_D*S + v_D*S²)
/// where η = number of singletons (external mutations).
pub fn phylo_fu_li_d(s : Int, n : Int, singletons : Int) -> Double {
  if n <= 2 || s == 0 {
    return 0.0
  }
  let a_n = sfs_harmonic(n - 1)
  let b_n = sfs_harmonic_sq(n - 1) // sum of 1/i^2
  let v_d = 1.0 + a_n * a_n / (b_n * (n.to_double() - 2.0)) -
    2.0 * a_n / ((n.to_double() - 1.0) * b_n)
  let u_d = a_n - 1.0 - v_d
  let s_d = s.to_double()
  let eta = singletons.to_double()
  let numerator = s_d - a_n * eta
  let denominator = sfs_sqrt(u_d * s_d + v_d * s_d * s_d)
  if denominator == 0.0 {
    return 0.0
  }
  numerator / denominator
}

///|
/// Fu and Li's F test statistic.
/// F_FL = (π - η) / sqrt(u_F*S + v_F*S²)
pub fn phylo_fu_li_f(
  s : Int,
  n : Int,
  pi : Double,
  singletons : Int
) -> Double {
  if n <= 2 || s == 0 {
    return 0.0
  }
  let a_n = sfs_harmonic(n - 1)
  let b_n = sfs_harmonic_sq(n - 1) // sum of 1/i^2
  let v_f = 1.0 + (n.to_double() * n.to_double() - 3.0 * n.to_double() + 2.0) /
    (3.0 * (n.to_double() - 1.0) * (n.to_double() - 2.0)) -
    2.0 * a_n / (n.to_double() - 1.0) +
    2.0 * a_n * a_n / ((n.to_double() - 1.0) * (n.to_double() - 1.0) * b_n)
  let u_f = 1.0 + (n.to_double() + 1.0) / (3.0 * (n.to_double() - 1.0)) -
    4.0 * (n.to_double() + 1.0) / (
      (n.to_double() - 1.0) * (n.to_double() - 1.0)
    ) * (a_n - 1.0) - v_f
  let s_d = s.to_double()
  let eta = singletons.to_double()
  let numerator = pi - eta
  let denominator = sfs_sqrt(u_f * s_d + v_f * s_d * s_d)
  if denominator == 0.0 {
    return 0.0
  }
  numerator / denominator
}

///|
/// Fay and Wu's H test statistic.
/// H = π - θ_H, where θ_H = (2 / (n(n-1))) * Σ i² * ξ_i
/// `sfs` is the observed unfolded SFS (counts for i=1..n-1).
pub fn phylo_fay_wu_h(sfs : Array[Int], n : Int, pi : Double) -> Double {
  if n <= 1 || sfs.length() == 0 {
    return 0.0
  }
  let mut theta_h = 0.0
  for i in 0..<sfs.length() {
    let freq = (i + 1).to_double()
    theta_h = theta_h + freq * freq * sfs[i].to_double()
  }
  theta_h = theta_h * 2.0 / (n.to_double() * (n.to_double() - 1.0))
  pi - theta_h
}

///|
/// Zeng's E test statistic (E = θ_π - θ_L).
/// θ_L = (1 / (n-1)) * Σ i * ξ_i
pub fn phylo_zeng_e(sfs : Array[Int], n : Int, pi : Double) -> Double {
  if n <= 1 || sfs.length() == 0 {
    return 0.0
  }
  let mut theta_l = 0.0
  for i in 0..<sfs.length() {
    let freq = (i + 1).to_double()
    theta_l = theta_l + freq * sfs[i].to_double()
  }
  theta_l = theta_l / (n.to_double() - 1.0)
  pi - theta_l
}

///|
/// Sturge's rule for optimal number of bins for SFS histogram.
/// k = 1 + log2(n)
pub fn phylo_sfs_num_bins(n : Int) -> Int {
  if n <= 1 {
    return 0
  }
  // log2(n) = ln(n) / ln(2)
  let log2n = @math.ln(n.to_double()) / @math.ln(2.0)
  1 + log2n.to_int()
}

///|
/// Format SFS as a table string.
pub fn phylo_sfs_to_string(sfs : Array[Double]) -> String {
  let sb = StringBuilder()
  sb.write_string("Site Frequency Spectrum:\n")
  sb.write_string("  i    E[ξ_i]\n")
  sb.write_string("  ───────────\n")
  for i in 0..<sfs.length() {
    sb.write_string("  ")
    let freq = (i + 1).to_string()
    sb.write_string(freq)
    // Pad
    let pad = if freq.length() < 4 { 4 - freq.length() } else { 1 }
    for _ in 0..<pad {
      sb.write_char(' ')
    }
    let val = (sfs[i] * 10000.0).to_int().to_double() / 10000.0
    sb.write_string(val.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

// --- helper functions ---

///|
/// Harmonic number H_k = sum_{i=1}^{k} 1/i.
fn sfs_harmonic(k : Int) -> Double {
  let mut h = 0.0
  for i in 1..<(k + 1) {
    h = h + 1.0 / i.to_double()
  }
  h
}

///|
/// Sum of 1/i^2 for i=1..k.
fn sfs_harmonic_sq(k : Int) -> Double {
  let mut h = 0.0
  for i in 1..<(k + 1) {
    let d = i.to_double()
    h = h + 1.0 / (d * d)
  }
  h
}

///|
/// Square root via Newton's method.
fn sfs_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  for _ in 0..<50 {
    guess = (guess + x / guess) / 2.0
  }
  guess
}

// =========================================================================
// SFS Stats Context (P1 optimization)
//   Inspired by Scanpy's "compute once, query many" pattern.
//   phylo_observed_sfs / phylo_segregating_sites / phylo_count_singletons /
//   phylo_pi each scan the entire alignment independently. When computing
//   multiple test statistics (Tajima D + Fu&Li D/F + Fay&Wu H + Zeng E),
//   the alignment is scanned 4 times redundantly.
//
//   SfsStats computes SFS, S, π, η in a single pass, then exposes all
//   test statistics as methods that reuse the cached values.
// =========================================================================

///|
/// Cached SFS statistics context: computes SFS, S, π, singletons in one
/// pass over the alignment, then provides all neutrality test statistics.
///
/// Usage:
///   let stats = SfsStats::from_sequences(seqs)
///   let d = stats.tajima_d()
///   let f = stats.fu_li_f()
///   let h = stats.fay_wu_h()
///   let e = stats.zeng_e()
pub struct SfsStats {
  n : Int                        // number of sequences
  sfs : Array[Int]               // observed SFS (ξ_1, ξ_2, ..., ξ_{n-1})
  s : Int                        // segregating sites = Σ ξ_i
  pi : Double                    // nucleotide diversity
  singletons : Int               // η = ξ_1
}

///|
/// Construct an SfsStats by scanning the alignment once.
/// Computes SFS, S, π, and singletons in a single pass.
pub fn SfsStats::from_sequences(sequences : Array[String]) -> SfsStats {
  let n = sequences.length()
  if n <= 1 {
    return SfsStats::{ n, sfs: [], s: 0, pi: 0.0, singletons: 0 }
  }
  let seq_len = sequences[0].length()
  // SFS: ξ_1 .. ξ_{n-1}
  let sfs : Array[Int] = Array(capacity=n - 1)
  for _ in 1..<n {
    sfs.push(0)
  }
  let mut s = 0
  let mut singletons = 0
  // For π: accumulate pairwise differences
  let mut total_pairwise_diff = 0
  // Single pass over all positions
  for pos in 0..<seq_len {
    let ancestral = sequences[0].char_to_at(pos)
    let mut derived_count = 0
    // For π: count pairwise differences at this position
    // d(i,j) = 1 if seq[i][pos] != seq[j][pos]
    // Instead of nested loop per position, collect allele counts
    let allele_counts : Map[Int, Int] = Map([], capacity=4)
    for i in 0..<n {
      let c = sequences[i].char_to_at(pos)
      match allele_counts.get(c) {
        Some(v) => allele_counts.set(c, v + 1)
        None => allele_counts.set(c, 1)
      }
    }
    // Pairwise differences at this position = total pairs - same-allele pairs
    // = C(n,2) - Σ C(count[c], 2)
    let total_pairs = n * (n - 1) / 2
    let mut same_pairs = 0
    for (_, cnt) in allele_counts {
      same_pairs = same_pairs + cnt * (cnt - 1) / 2
    }
    total_pairwise_diff = total_pairwise_diff + (total_pairs - same_pairs)
    // SFS and singletons (using first sequence as ancestral)
    for i in 1..<n {
      if sequences[i].char_to_at(pos) != ancestral {
        derived_count = derived_count + 1
      }
    }
    if derived_count > 0 && derived_count < n {
      sfs[derived_count - 1] = sfs[derived_count - 1] + 1
      s = s + 1
      if derived_count == 1 {
        singletons = singletons + 1
      }
    }
  }
  let num_pairs = n.to_double() * (n - 1).to_double() / 2.0
  let pi = if num_pairs > 0.0 {
    total_pairwise_diff.to_double() / num_pairs
  } else {
    0.0
  }
  SfsStats::{ n, sfs, s, pi, singletons }
}

///|
/// Number of sequences.
pub fn SfsStats::n_sequences(self : SfsStats) -> Int {
  self.n
}

///|
/// Segregating sites S = Σ ξ_i.
pub fn SfsStats::segregating_sites(self : SfsStats) -> Int {
  self.s
}

///|
/// Observed SFS array (ξ_1, ..., ξ_{n-1}).
pub fn SfsStats::sfs(self : SfsStats) -> Array[Int] {
  self.sfs
}

///|
/// Nucleotide diversity π.
pub fn SfsStats::nucleotide_diversity(self : SfsStats) -> Double {
  self.pi
}

///|
/// Singleton count η = ξ_1.
pub fn SfsStats::singleton_count(self : SfsStats) -> Int {
  self.singletons
}

///|
/// Watterson's θ_W = S / H_{n-1}.
pub fn SfsStats::watterson_theta(self : SfsStats) -> Double {
  phylo_watterson_theta(self.s, self.n)
}

///|
/// Tajima's D from cached S and π.
pub fn SfsStats::tajima_d(self : SfsStats) -> Double {
  phylo_tajima_d(self.s, self.n, self.pi)
}

///|
/// Fu & Li's D from cached S and singletons.
pub fn SfsStats::fu_li_d(self : SfsStats) -> Double {
  phylo_fu_li_d(self.s, self.n, self.singletons)
}

///|
/// Fu & Li's F from cached S, π, and singletons.
pub fn SfsStats::fu_li_f(self : SfsStats) -> Double {
  phylo_fu_li_f(self.s, self.n, self.pi, self.singletons)
}

///|
/// Fay & Wu's H from cached SFS and π.
pub fn SfsStats::fay_wu_h(self : SfsStats) -> Double {
  phylo_fay_wu_h(self.sfs, self.n, self.pi)
}

///|
/// Zeng's E from cached SFS and π.
pub fn SfsStats::zeng_e(self : SfsStats) -> Double {
  phylo_zeng_e(self.sfs, self.n, self.pi)
}
