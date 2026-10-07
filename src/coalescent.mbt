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
      i2 = (i2 + 1) % k
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
  let depths : Array[Int] = Array(capacity=0)
  phylo_rect_helper(tree, "", 0, lines, depths)
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
  lines : Array[String],
  depths : Array[Int]
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
    depths.push(depth)
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
          depths.push(depth)
        }
      }
      phylo_rect_helper(child, child_prefix, depth + 1, lines, depths)
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
