///|
/// Motif stacking and clustering visualization, ported from
/// Bioconductor's `motifStack` package.
///
/// The reference package computes pairwise similarities between
/// position count/weight matrices, clusters them into a tree, and
/// stacks aligned sequence logos in the order given by that tree.
/// This port provides Pearson-distance matrices, average-linkage
/// hierarchical clustering, leaf-order extraction and ASCII logo
/// stacking over DNA PWMs (4 x L, rows A C G T).

///|
/// A named DNA motif with a 4 x L PWM.
pub struct StackMotif {
  name : String
  pwm : Array[Array[Double]]
} derive(Eq, Show, Debug)

///|
pub fn StackMotif::new(
  name~ : String,
  pwm~ : Array[Array[Double]],
) -> StackMotif {
  StackMotif::{ name, pwm }
}

///|
/// Length of the motif in columns.
pub fn StackMotif::len(self : StackMotif) -> Int {
  if self.pwm.is_empty() {
    0
  } else {
    self.pwm[0].length()
  }
}

///|
/// Flatten a 4 x L PWM into a length-4L vector for distance
/// computation.
pub fn motifstack_flatten(motif : StackMotif) -> Array[Double] {
  let flat : Array[Double] = Array::new()
  for row in motif.pwm {
    for v in row {
      flat.push(v)
    }
  }
  flat
}

///|
/// Pearson correlation of two equal-length vectors.
fn motifstack_pearson(x : Array[Double], y : Array[Double]) -> Double {
  let n = if x.length() < y.length() { x.length() } else { y.length() }
  if n < 2 {
    return 0.0
  }
  let mut sx = 0.0
  let mut sy = 0.0
  for i in 0..<n {
    sx += x[i]
    sy += y[i]
  }
  let mx = sx / n.to_double()
  let my = sy / n.to_double()
  let mut sxy = 0.0
  let mut sxx = 0.0
  let mut syy = 0.0
  for i in 0..<n {
    let dx = x[i] - mx
    let dy = y[i] - my
    sxy += dx * dy
    sxx += dx * dx
    syy += dy * dy
  }
  let denom = motifstack_sqrt(sxx * syy)
  if denom <= 0.0 {
    0.0
  } else {
    sxy / denom
  }
}

///|
fn motifstack_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  for _ in 0..<50 {
    guess = (guess + x / guess) / 2.0
  }
  guess
}

///|
/// Pearson distance (1 - r) between two motifs.
pub fn motifstack_distance(
  m1 : StackMotif,
  m2 : StackMotif,
) -> Double {
  let r = motifstack_pearson(
    motifstack_flatten(m1),
    motifstack_flatten(m2),
  )
  1.0 - r
}

///|
/// Symmetric pairwise distance matrix over a motif collection.
pub fn motifstack_distance_matrix(
  motifs : Array[StackMotif],
) -> Array[Array[Double]] {
  let n = motifs.length()
  let matrix : Array[Array[Double]] = Array::new()
  for _i in 0..<n {
    matrix.push(Array::make(n, 0.0))
  }
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d = motifstack_distance(motifs[i], motifs[j])
      matrix[i][j] = d
      matrix[j][i] = d
    }
  }
  matrix
}

///|
/// A node in the motif clustering tree: either a leaf holding a motif
/// index or an internal merge with a linkage height.
pub(all) enum MotifStackNode {
  Leaf(Int)
  Merge(MotifStackNode, MotifStackNode, Double)
}

///|
/// Collect leaves of a subtree in order.
pub fn MotifStackNode::leaves(self : MotifStackNode) -> Array[Int] {
  let out : Array[Int] = Array::new()
  match self {
    MotifStackNode::Leaf(i) => out.push(i)
    MotifStackNode::Merge(l, r, _) => {
      let left = l.leaves()
      for v in left {
        out.push(v)
      }
      let right = r.leaves()
      for v in right {
        out.push(v)
      }
    }
  }
  out
}

///|
/// Linkage height of a node (0.0 for leaves).
pub fn MotifStackNode::height(self : MotifStackNode) -> Double {
  match self {
    MotifStackNode::Leaf(_) => 0.0
    MotifStackNode::Merge(_, _, h) => h
  }
}

///|
/// Average-linkage hierarchical clustering over a symmetric distance
/// matrix.  Returns the root node, or None for an empty collection.
pub fn motifstack_hclust(
  distances : Array[Array[Double]],
) -> MotifStackNode? {
  let n = distances.length()
  if n == 0 {
    return None
  }
  if n == 1 {
    return Some(MotifStackNode::Leaf(0))
  }
  let clusters : Array[MotifStackNode] = Array::new()
  let members : Array[Array[Int]] = Array::new()
  for i in 0..<n {
    clusters.push(MotifStackNode::Leaf(i))
    members.push([i])
  }
  while clusters.length() > 1 {
    let m = clusters.length()
    let mut best_i = 0
    let mut best_j = 1
    let mut best_d = 1.0e18
    for a in 0..<m {
      let hi = a + 1
      for b in hi..<m {
        let mut sum = 0.0
        let mut count = 0
        for ia in members[a] {
          let jb = b
          for ib in members[jb] {
            sum += distances[ia][ib]
            count += 1
          }
        }
        let avg = if count == 0 {
          0.0
        } else {
          sum / count.to_double()
        }
        if avg < best_d {
          best_d = avg
          best_i = a
          best_j = b
        }
      }
    }
    let node = MotifStackNode::Merge(
      clusters[best_i],
      clusters[best_j],
      best_d,
    )
    let merged_members = members[best_i] + members[best_j]
    let new_clusters : Array[MotifStackNode] = Array::new()
    let new_members : Array[Array[Int]] = Array::new()
    for a in 0..<m {
      if a != best_i && a != best_j {
        new_clusters.push(clusters[a])
        new_members.push(members[a])
      }
    }
    new_clusters.push(node)
    new_members.push(merged_members)
    clusters.clear()
    for c in new_clusters {
      clusters.push(c)
    }
    members.clear()
    for mm in new_members {
      members.push(mm)
    }
  }
  Some(clusters[0])
}

///|
/// Stack order for a motif collection: leaf order of its clustering
/// tree.  Falls back to identity order when the tree cannot be built.
pub fn motifstack_stack_order(motifs : Array[StackMotif]) -> Array[Int] {
  let distances = motifstack_distance_matrix(motifs)
  match motifstack_hclust(distances) {
    Some(root) => root.leaves()
    None => []
  }
}

///|
/// Per-column information content of a DNA motif:
/// IC_j = sum_i p * log2(p / 0.25).
pub fn motifstack_column_ic(motif : StackMotif) -> Array[Double] {
  let width = motif.len()
  let ics : Array[Double] = Array::new()
  for j in 0..<width {
    let mut ic = 0.0
    for i in 0..<4 {
      let p = motif.pwm[i][j]
      if p > 0.0 {
        ic += p * (motifstack_log2(p) + 2.0)
      }
    }
    ics.push(ic)
  }
  ics
}

///|
/// log2 via bisection.
fn motifstack_log2(x : Double) -> Double {
  if x <= 0.0 {
    return -20.0
  }
  let mut lo = -20.0
  let mut hi = 20.0
  for _ in 0..<60 {
    let mid = (lo + hi) / 2.0
    let pow = motifstack_pow2(mid)
    if pow < x {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
fn motifstack_pow2(x : Double) -> Double {
  let ln2 = 0.6931471805599453
  motifstack_exp(x * ln2)
}

///|
fn motifstack_exp(x : Double) -> Double {
  if x > 40.0 {
    return motifstack_exp(40.0) * motifstack_exp(x - 40.0)
  }
  if x < -40.0 {
    return 0.0
  }
  let mut sum = 1.0
  let mut term = 1.0
  for i in 1..<30 {
    term = term * x / i.to_double()
    sum += term
    let abs_term = if term < 0.0 { 0.0 - term } else { term }
    if abs_term < 1.0e-12 {
      break
    }
  }
  sum
}

///|
/// Render one motif as an ASCII logo block: `rows` lines where column
/// j shows the dominant base scaled by per-column IC, then a position
/// ruler line.
pub fn motifstack_ascii_logo(
  motif : StackMotif,
  rows~ : Int = 6,
) -> String {
  let bases = ["A", "C", "G", "T"]
  let width = motif.len()
  let ics = motifstack_column_ic(motif)
  let mut max_ic = 0.0
  for ic in ics {
    if ic > max_ic {
      max_ic = ic
    }
  }
  let mut out = ""
  if max_ic <= 0.0 {
    for _r in 0..<rows {
      out = out + "\n"
    }
    return out
  }
  let mut row = rows - 1
  while row >= 0 {
    let mut line = ""
    for j in 0..<width {
      let mut best_i = 0
      let mut best_p = -1.0
      for i in 0..<4 {
        if motif.pwm[i][j] > best_p {
          best_p = motif.pwm[i][j]
          best_i = i
        }
      }
      let scaled = ics[j] / max_ic
      let threshold = (row.to_double() + 0.5) / rows.to_double()
      if scaled >= threshold {
        line = line + bases[best_i]
      } else {
        line = line + " "
      }
    }
    out = out + line + "\n"
    row -= 1
  }
  let mut ruler = ""
  for j in 0..<width {
    ruler = ruler + (j + 1).to_string()
  }
  out = out + ruler + "\n"
  out
}

///|
/// Stack several motifs in clustering leaf order: for each motif a
/// name line followed by its ASCII logo block.
pub fn motifstack_stack_ascii(
  motifs : Array[StackMotif],
  rows~ : Int = 6,
) -> String {
  let order = motifstack_stack_order(motifs)
  let mut out = ""
  for idx in order {
    let motif = motifs[idx]
    out = out + "== " + motif.name + " ==\n"
    out = out + motifstack_ascii_logo(motif, rows~)
  }
  out
}

///|
/// Two related PWMs plus one divergent motif, for demos and tests.
pub fn motifstack_sample_motifs() -> Array[StackMotif] {
  let m1 = StackMotif::new(
    name="GATA1",
    pwm=[
      [0.1, 0.1, 0.8, 0.1, 0.1],
      [0.1, 0.1, 0.1, 0.1, 0.8],
      [0.7, 0.1, 0.1, 0.7, 0.1],
      [0.1, 0.7, 0.0, 0.1, 0.0],
    ],
  )
  let m2 = StackMotif::new(
    name="GATA2",
    pwm=[
      [0.1, 0.1, 0.7, 0.1, 0.1],
      [0.1, 0.1, 0.2, 0.1, 0.7],
      [0.7, 0.1, 0.1, 0.7, 0.1],
      [0.1, 0.7, 0.0, 0.1, 0.0],
    ],
  )
  let m3 = StackMotif::new(
    name="CTCF",
    pwm=[
      [0.6, 0.1, 0.1, 0.2],
      [0.1, 0.1, 0.1, 0.7],
      [0.2, 0.7, 0.1, 0.1],
      [0.1, 0.1, 0.7, 0.0],
    ],
  )
  [m1, m2, m3]
}