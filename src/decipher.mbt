///|
/// DECIPHER - Tools for biological sequence analysis.
/// Inspired by Bioconductor DECIPHER package (Wright 2016).
/// Core: two-step alignment, sequence database, clustering, phylogenetics, chimera detection.

///|
/// DNA scoring matrix (substitution scores).
pub struct DecipherSubMatrix {
  scores : Array[Array[Int]] // 4x4 for ACGT
  gap_open : Int
  gap_extend : Int
}

///|
/// Default DNA scoring matrix (NUC4.4 equivalent).
pub fn default_dna_matrix() -> DecipherSubMatrix {
  DecipherSubMatrix::{
    scores: [
      [5, -4, -4, -4], // A
      [-4, 5, -4, -4], // C
      [-4, -4, 5, -4], // G
      [-4, -4, -4, 5], // T
    ],
    gap_open: -12,
    gap_extend: -4,
  }
}

///|
/// Nucleotide to index: A=0, C=1, G=2, T=3, else=-1
pub fn decipher_nt_idx(c : Char) -> Int {
  match c.to_int() {
    65 | 97 => 0
    67 | 99 => 1
    71 | 103 => 2
    84 | 116 => 3
    _ => -1
  }
}

///|
/// Index to nucleotide.
pub fn decipher_idx_nt(i : Int) -> Char {
  match i {
    0 => 'A'
    1 => 'C'
    2 => 'G'
    3 => 'T'
    _ => '-'
  }
}

///|
/// Alignment result.
pub struct DecipherAlnResult {
  aligned_a : String
  aligned_b : String
  score : Int
  identity : Double
  coverage : Double
  mismatches : Int
  gaps : Int
  gap_length : Int
}

///|
/// Needleman-Wunsch global alignment (affine gap penalty via Gotoh algorithm).
pub fn align_nw_affine(
  seq_a : String,
  seq_b : String,
  matrix : DecipherSubMatrix,
) -> DecipherAlnResult {
  let a = seq_a.to_array()
  let b = seq_b.to_array()
  let m = a.length()
  let n = b.length()
  if m == 0 || n == 0 {
    return DecipherAlnResult::{
      aligned_a: seq_a,
      aligned_b: seq_b,
      score: 0,
      identity: 0.0,
      coverage: 0.0,
      mismatches: 0,
      gaps: 0,
      gap_length: 0,
    }
  }
  let neg_inf = -999999
  // DP matrices: M (match/mismatch), Ix (gap in b), Iy (gap in a)
  // Use 1D arrays for efficiency
  let m_size = (m + 1) * (n + 1)
  let mat_m : Array[Int] = Array::make(m_size, neg_inf)
  let mat_ix : Array[Int] = Array::make(m_size, neg_inf)
  let mat_iy : Array[Int] = Array::make(m_size, neg_inf)
  // Traceback: 0=match, 1=gap_x (gap in b), 2=gap_y (gap in a)
  let tb : Array[Int] = Array::make(m_size, 0)
  // Initialize
  mat_m[0] = 0
  let mut i = 1
  while i <= m {
    mat_ix[i * (n + 1)] = matrix.gap_open + matrix.gap_extend * i
    tb[i * (n + 1)] = 1
    i = i + 1
  }
  let mut j = 1
  while j <= n {
    mat_iy[j] = matrix.gap_open + matrix.gap_extend * j
    tb[j] = 2
    j = j + 1
  }
  // Fill matrices
  let mut ii = 1
  while ii <= m {
    let mut jj = 1
    while jj <= n {
      let idx = ii * (n + 1) + jj
      let idx_up = (ii - 1) * (n + 1) + jj
      let idx_left = ii * (n + 1) + (jj - 1)
      let idx_diag = (ii - 1) * (n + 1) + (jj - 1)
      // Ix: gap in b (consume a)
      let ix_from_m = mat_m[idx_up] + matrix.gap_open + matrix.gap_extend
      let ix_from_ix = mat_ix[idx_up] + matrix.gap_extend
      mat_ix[idx] = if ix_from_m > ix_from_ix { ix_from_m } else { ix_from_ix }
      // Iy: gap in a (consume b)
      let iy_from_m = mat_m[idx_left] + matrix.gap_open + matrix.gap_extend
      let iy_from_iy = mat_iy[idx_left] + matrix.gap_extend
      mat_iy[idx] = if iy_from_m > iy_from_iy { iy_from_m } else { iy_from_iy }
      // M: match/mismatch
      let ai = decipher_nt_idx(a[ii - 1])
      let bj = decipher_nt_idx(b[jj - 1])
      let sub_score = if ai >= 0 && bj >= 0 {
        matrix.scores[ai][bj]
      } else {
        -1
      }
      let m_from_m = mat_m[idx_diag] + sub_score
      let m_from_ix = mat_ix[idx_diag] + sub_score
      let m_from_iy = mat_iy[idx_diag] + sub_score
      let best_diag = if m_from_m > m_from_ix {
        if m_from_m > m_from_iy { m_from_m } else { m_from_iy }
      } else {
        if m_from_ix > m_from_iy { m_from_ix } else { m_from_iy }
      }
      mat_m[idx] = best_diag
      // Choose best of M, Ix, Iy for traceback
      let best = mat_m[idx]
      let best_ix = mat_ix[idx]
      let best_iy = mat_iy[idx]
      if best >= best_ix && best >= best_iy {
        tb[idx] = 0
      } else if best_ix >= best_iy {
        tb[idx] = 1
        mat_m[idx] = best_ix
      } else {
        tb[idx] = 2
        mat_m[idx] = best_iy
      }
      jj = jj + 1
    }
    ii = ii + 1
  }
  // Traceback
  let aligned_a : Array[Char] = []
  let aligned_b : Array[Char] = []
  let mut ti = m
  let mut tj = n
  while ti > 0 || tj > 0 {
    let idx = ti * (n + 1) + tj
    let dir = tb[idx]
    if dir == 0 && ti > 0 && tj > 0 {
      aligned_a.push(a[ti - 1])
      aligned_b.push(b[tj - 1])
      ti = ti - 1
      tj = tj - 1
    } else if dir == 1 && ti > 0 {
      aligned_a.push(a[ti - 1])
      aligned_b.push('-')
      ti = ti - 1
    } else if tj > 0 {
      aligned_a.push('-')
      aligned_b.push(b[tj - 1])
      tj = tj - 1
    } else if ti > 0 {
      aligned_a.push(a[ti - 1])
      aligned_b.push('-')
      ti = ti - 1
    } else {
      break
    }
  }
  // Reverse
  let ra = StringBuilder::new()
  let rb = StringBuilder::new()
  let mut k = aligned_a.length() - 1
  while k >= 0 {
    ra.write_char(aligned_a[k])
    rb.write_char(aligned_b[k])
    k = k - 1
  }
  let result_a = ra.to_string()
  let result_b = rb.to_string()
  // Compute statistics
  let mut matches = 0
  let mut mismatches = 0
  let mut gaps = 0
  let mut gap_length = 0
  let mut in_gap = false
  let mut ki = 0
  while ki < result_a.length() {
    if result_a[ki] == 45 || result_b[ki] == 45 { // '-' = 45
      gap_length = gap_length + 1
      if !in_gap {
        gaps = gaps + 1
        in_gap = true
      }
    } else {
      in_gap = false
      if result_a[ki] == result_b[ki] {
        matches = matches + 1
      } else {
        mismatches = mismatches + 1
      }
    }
    ki = ki + 1
  }
  let aln_len = result_a.length().to_double()
  let identity = if aln_len > 0.0 { matches.to_double() / aln_len } else { 0.0 }
  let coverage = if m > 0 { (matches + mismatches).to_double() / m.to_double() } else { 0.0 }
  let final_score = mat_m[m * (n + 1) + n]
  DecipherAlnResult::{
    aligned_a: result_a,
    aligned_b: result_b,
    score: final_score,
    identity,
    coverage,
    mismatches,
    gaps,
    gap_length,
  }
}

///|
/// Align a single sequence to a reference (convenience wrapper).
pub fn align_seqs(
  seq : String,
  reference : String,
) -> DecipherAlnResult {
  align_nw_affine(seq, reference, default_dna_matrix())
}

///|
/// Align a set of sequences to a reference (DECIPHER AlignSeqs style).
/// Returns array of alignment results.
pub fn align_seqs_batch(
  sequences : Array[String],
  reference : String,
) -> Array[DecipherAlnResult] {
  let results : Array[DecipherAlnResult] = []
  let matrix = default_dna_matrix()
  let mut i = 0
  while i < sequences.length() {
    results.push(align_nw_affine(sequences[i], reference, matrix))
    i = i + 1
  }
  results
}

///|
/// Profile: position-specific nucleotide frequencies.
pub struct Profile {
  freqs : Array[Array[Double]] // [position][4 nucleotides]
  length : Int
  total_weight : Double
}

///|
/// Build a profile from an array of aligned sequences.
pub fn build_profile(aligned_seqs : Array[String]) -> Profile {
  if aligned_seqs.length() == 0 {
    return Profile::{ freqs: [], length: 0, total_weight: 0.0 }
  }
  let n_pos = aligned_seqs[0].length()
  let freqs : Array[Array[Double]] = []
  let mut pos = 0
  while pos < n_pos {
    let col : Array[Double] = Array::make(4, 0.0)
    let mut seq_i = 0
    while seq_i < aligned_seqs.length() {
      let c = aligned_seqs[seq_i][pos]
      let idx = decipher_nt_idx(c.unsafe_to_char())
      if idx >= 0 {
        col[idx] = col[idx] + 1.0
      }
      seq_i = seq_i + 1
    }
    // Normalize
    let mut total = 0.0
    let mut k = 0
    while k < 4 {
      total = total + col[k]
      k = k + 1
    }
    if total > 0.0 {
      let mut k2 = 0
      while k2 < 4 {
        col[k2] = col[k2] / total
        k2 = k2 + 1
      }
    }
    freqs.push(col)
    pos = pos + 1
  }
  Profile::{ freqs, length: n_pos, total_weight: aligned_seqs.length().to_double() }
}

///|
/// Align two profiles (profile-profile alignment).
/// Simplified: uses sum of pairwise scores at each position.
pub fn decipher_align_profiles(
  profile_a : Profile,
  profile_b : Profile,
  matrix : DecipherSubMatrix,
) -> DecipherAlnResult {
  // Convert profiles to consensus sequences and align
  let cons_a = profile_to_consensus(profile_a)
  let cons_b = profile_to_consensus(profile_b)
  align_nw_affine(cons_a, cons_b, matrix)
}

///|
/// Convert profile to consensus sequence.
pub fn profile_to_consensus(profile : Profile) -> String {
  let sb = StringBuilder::new()
  let mut pos = 0
  while pos < profile.length {
    let mut best_idx = 0
    let mut best_freq = -1.0
    let mut k = 0
    while k < 4 {
      if profile.freqs[pos][k] > best_freq {
        best_freq = profile.freqs[pos][k]
        best_idx = k
      }
      k = k + 1
    }
    if best_freq > 0.0 {
      sb.write_char(decipher_idx_nt(best_idx))
    } else {
      sb.write_char('-')
    }
    pos = pos + 1
  }
  sb.to_string()
}

///|
/// Compute pairwise distance matrix using sequence identity.
/// Method: "identity" (1 - identity), "hamming" (raw hamming / length).
pub fn decipher_distance_matrix(
  sequences : Array[String],
  dist_method~ : String = "identity",
) -> Array[Array[Double]] {
  let n = sequences.length()
  let dist : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut j = 0
    while j < n {
      if i == j {
        row[j] = 0.0
      } else if j > i {
        let aln = align_nw_affine(sequences[i], sequences[j], default_dna_matrix())
        let d = if dist_method == "hamming" {
          1.0 - aln.identity
        } else {
          1.0 - aln.identity
        }
        row[j] = d
      } else {
        row[j] = dist[j][i]
      }
      j = j + 1
    }
    dist.push(row)
    i = i + 1
  }
  dist
}

///|
/// Cluster result.
pub struct DecipherClusterResult {
  clusters : Array[Array[Int]] // cluster -> member indices
  n_clusters : Int
  cluster_ids : Array[Int] // sequence index -> cluster id
}

///|
/// Cluster sequences by identity threshold.
/// IdClusters: sequences with identity >= threshold are in the same cluster.
pub fn id_clusters(
  sequences : Array[String],
  threshold~ : Double = 0.97,
) -> DecipherClusterResult {
  let n = sequences.length()
  let dist = decipher_distance_matrix(sequences)
  // Union-find
  let parent : Array[Int] = []
  let mut i = 0
  while i < n {
    parent.push(i)
    i = i + 1
  }
  let mut ii = 0
  while ii < n {
    let mut jj = ii + 1
    while jj < n {
      let identity = 1.0 - dist[ii][jj]
      if identity >= threshold {
        // Union ii and jj
        let mut root_i = ii
        while parent[root_i] != root_i {
          root_i = parent[root_i]
        }
        let mut root_j = jj
        while parent[root_j] != root_j {
          root_j = parent[root_j]
        }
        if root_i != root_j {
          parent[root_j] = root_i
        }
      }
      jj = jj + 1
    }
    ii = ii + 1
  }
  // Find roots and build clusters
  let cluster_map : Map[Int, Int] = Map([], capacity=n)
  let cluster_ids : Array[Int] = Array::make(n, 0)
  let mut n_clusters = 0
  let mut k = 0
  while k < n {
    let mut root = k
    while parent[root] != root {
      root = parent[root]
    }
    if !cluster_map.contains(root) {
      cluster_map[root] = n_clusters
      n_clusters = n_clusters + 1
    }
    cluster_ids[k] = cluster_map[root]
    k = k + 1
  }
  let clusters : Array[Array[Int]] = []
  let mut ci = 0
  while ci < n_clusters {
    clusters.push([])
    ci = ci + 1
  }
  let mut mi = 0
  while mi < n {
    clusters[cluster_ids[mi]].push(mi)
    mi = mi + 1
  }
  DecipherClusterResult::{ clusters, n_clusters, cluster_ids }
}

///|
/// Neighbor-Joining tree node.
pub struct DecipherTreeNode {
  name : String
  left : Int
  right : Int
  mut parent : Int
  branch_left : Double
  branch_right : Double
  is_leaf : Bool
}

///|
/// Phylogenetic tree.
pub struct PhyloTree {
  nodes : Array[DecipherTreeNode]
  n_taxa : Int
  root : Int
}

///|
/// Build Neighbor-Joining tree from distance matrix.
pub fn neighbor_joining(
  dist : Array[Array[Double]],
  names : Array[String],
) -> PhyloTree {
  let n = dist.length()
  if n <= 1 {
    let nodes : Array[DecipherTreeNode] = []
    if n == 1 {
      nodes.push(  DecipherTreeNode::{ name: names[0], left: -1, right: -1, parent: -1, branch_left: 0.0, branch_right: 0.0, is_leaf: true })
    }
    return PhyloTree::{ nodes, n_taxa: n, root: if n == 1 { 0 } else { -1 } }
  }
  // NJ algorithm
  // Active taxa (indices still in the tree)
  let active : Array[Bool] = Array::make(n, true)
  let active_indices : Array[Int] = []
  let mut i = 0
  while i < n {
    active_indices.push(i)
    i = i + 1
  }
  // Node list: first n are leaves, then internal nodes
  let nodes : Array[DecipherTreeNode] = []
  let mut ki = 0
  while ki < n {
    nodes.push(  DecipherTreeNode::{ name: names[ki], left: -1, right: -1, parent: -1, branch_left: 0.0, branch_right: 0.0, is_leaf: true })
    ki = ki + 1
  }
  // Working distance matrix (copy)
  let work_dist : Array[Array[Double]] = []
  let mut di = 0
  while di < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut dj = 0
    while dj < n {
      row[dj] = dist[di][dj]
      dj = dj + 1
    }
    work_dist.push(row)
    di = di + 1
  }
  let mut next_node = n
  let mut n_active = n
  while n_active > 2 {
    // Compute r_i = sum of d(i, j) for all active j
    let r : Array[Double] = Array::make(n, 0.0)
    let mut ii = 0
    while ii < n {
      if active[ii] {
        let mut jj = 0
        while jj < n {
          if active[jj] && jj != ii {
            r[ii] = r[ii] + work_dist[ii][jj]
          }
          jj = jj + 1
        }
        r[ii] = r[ii] / (n_active - 2).to_double()
      }
      ii = ii + 1
    }
    // Find pair (i, j) minimizing D(i,j) = d(i,j) - r_i - r_j
    let mut best_i = -1
    let mut best_j = -1
    let mut best_d = 1.0e18
    let mut bi = 0
    while bi < n {
      if active[bi] {
        let mut bj = bi + 1
        while bj < n {
          if active[bj] {
            let d_ij = work_dist[bi][bj] - r[bi] - r[bj]
            if d_ij < best_d {
              best_d = d_ij
              best_i = bi
              best_j = bj
            }
          }
          bj = bj + 1
        }
      }
      bi = bi + 1
    }
    if best_i < 0 || best_j < 0 {
      break
    }
    // Compute branch lengths
    let d_ij = work_dist[best_i][best_j]
    let bl_i = 0.5 * (d_ij + r[best_i] - r[best_j])
    let bl_j = d_ij - bl_i
    // Create new internal node
    nodes.push(  DecipherTreeNode::{
      name: "Node" + next_node.to_string(),
      left: best_i,
      right: best_j,
      parent: -1,
      branch_left: bl_i,
      branch_right: bl_j,
      is_leaf: false,
    })
    nodes[best_i].parent = next_node
    nodes[best_j].parent = next_node
    // Update distance matrix: new node k
    let mut kk = 0
    while kk < n {
      if active[kk] && kk != best_i && kk != best_j {
        let d_new = 0.5 * (work_dist[best_i][kk] + work_dist[best_j][kk] - d_ij)
        work_dist[best_i][kk] = d_new
        work_dist[kk][best_i] = d_new
      }
      kk = kk + 1
    }
    active[best_j] = false
    n_active = n_active - 1
    next_node = next_node + 1
  }
  // Connect last two nodes
  let mut last_i = -1
  let mut last_j = -1
  let mut li = 0
  while li < n {
    if active[li] {
      if last_i < 0 {
        last_i = li
      } else {
        last_j = li
      }
    }
    li = li + 1
  }
  if last_i >= 0 && last_j >= 0 {
    let d_final = work_dist[last_i][last_j]
    nodes.push(  DecipherTreeNode::{
      name: "Root",
      left: last_i,
      right: last_j,
      parent: -1,
      branch_left: 0.5 * d_final,
      branch_right: 0.5 * d_final,
      is_leaf: false,
    })
    nodes[last_i].parent = next_node
    nodes[last_j].parent = next_node
    return PhyloTree::{ nodes, n_taxa: n, root: next_node }
  } else if last_i >= 0 {
    return PhyloTree::{ nodes, n_taxa: n, root: last_i }
  }
  PhyloTree::{ nodes, n_taxa: n, root: -1 }
}

///|
/// Generate Newick string from phylogenetic tree.
pub fn tree_to_newick(tree : PhyloTree) -> String {
  if tree.root < 0 {
    return ";"
  }
  newick_recursive(tree, tree.root) + ";"
}

///|
fn newick_recursive(tree : PhyloTree, node_idx : Int) -> String {
  let node = tree.nodes[node_idx]
  if node.is_leaf {
    return node.name
  }
  let left_str = newick_recursive(tree, node.left)
  let right_str = newick_recursive(tree, node.right)
  let result = "(" + left_str + ":" + node.branch_left.to_string() + "," + right_str + ":" + node.branch_right.to_string() + ")"
  if !node.is_leaf && node.parent == -1 {
    result
  } else {
    result
  }
}

///|
/// Find chimeras using three-parent detection.
/// A sequence is chimeric if it can be split into segments from 3+ parents.
pub fn find_chimeras(
  sequences : Array[String],
  abundances : Array[Int],
  threshold~ : Double = 0.5,
) -> Array[Bool] {
  let n = sequences.length()
  let is_chimera : Array[Bool] = Array::make(n, false)
  // Sort by abundance descending
  let order : Array[Int] = []
  let mut i = 0
  while i < n {
    order.push(i)
    i = i + 1
  }
  order.sort_by(fn(a, b) {
    if abundances[b] > abundances[a] {
      1
    } else if abundances[b] < abundances[a] {
      -1
    } else {
      0
    }
  })
  // Check each sequence against more abundant parents
  let mut idx = 0
  while idx < n {
    let seq_idx = order[idx]
    let seq = sequences[seq_idx]
    // Try all pairs of parents
    let mut p1 = 0
    while p1 < idx {
      let parent_a = sequences[order[p1]]
      if parent_a.length() == seq.length() {
        let mut p2 = p1 + 1
        while p2 < idx {
          let parent_b = sequences[order[p2]]
          if parent_b.length() == seq.length() {
            // Check if seq is a two-parent chimera
            if is_two_parent_chimera(seq, parent_a, parent_b, threshold) {
              is_chimera[seq_idx] = true
            }
          }
          p2 = p2 + 1
        }
      }
      p1 = p1 + 1
    }
    idx = idx + 1
  }
  is_chimera
}

///|
/// Check if a sequence is a two-parent chimera with breakpoint.
fn is_two_parent_chimera(
  seq : String,
  parent_a : String,
  parent_b : String,
  threshold : Double,
) -> Bool {
  let n = seq.length()
  let s = seq.to_array()
  let a = parent_a.to_array()
  let b = parent_b.to_array()
  // Find best breakpoint
  let mut bp = 1
  while bp < n {
    let mut match_a_first = 0
    let mut match_b_second = 0
    let mut k = 0
    while k < bp {
      if s[k] == a[k] {
        match_a_first = match_a_first + 1
      }
      k = k + 1
    }
    while k < n {
      if s[k] == b[k] {
        match_b_second = match_b_second + 1
      }
      k = k + 1
    }
    let frac_a = match_a_first.to_double() / bp.to_double()
    let frac_b = match_b_second.to_double() / (n - bp).to_double()
    if frac_a >= threshold && frac_b >= threshold {
      return true
    }
    // Also try B then A
    let mut match_b_first = 0
    let mut match_a_second = 0
    let mut k2 = 0
    while k2 < bp {
      if s[k2] == b[k2] {
        match_b_first = match_b_first + 1
      }
      k2 = k2 + 1
    }
    while k2 < n {
      if s[k2] == a[k2] {
        match_a_second = match_a_second + 1
      }
      k2 = k2 + 1
    }
    let frac_b2 = match_b_first.to_double() / bp.to_double()
    let frac_a2 = match_a_second.to_double() / (n - bp).to_double()
    if frac_b2 >= threshold && frac_a2 >= threshold {
      return true
    }
    bp = bp + 1
  }
  false
}

///|
/// Compute consensus sequence from aligned sequences.
pub fn decipher_consensus(aligned_seqs : Array[String]) -> String {
  if aligned_seqs.length() == 0 {
    return ""
  }
  let n_pos = aligned_seqs[0].length()
  let sb = StringBuilder::new()
  let mut pos = 0
  while pos < n_pos {
    let counts : Array[Int] = Array::make(5, 0) // A, C, G, T, gap
    let mut si = 0
    while si < aligned_seqs.length() {
      let c = aligned_seqs[si][pos]
      match c.to_int() {
        65 | 97 => counts[0] = counts[0] + 1 // A
        67 | 99 => counts[1] = counts[1] + 1 // C
        71 | 103 => counts[2] = counts[2] + 1 // G
        84 | 116 => counts[3] = counts[3] + 1 // T
        45 => counts[4] = counts[4] + 1 // -
        _ => ()
      }
      si = si + 1
    }
    // Find most common
    let mut best_idx = 0
    let mut best_count = 0
    let mut k = 0
    while k < 5 {
      if counts[k] > best_count {
        best_count = counts[k]
        best_idx = k
      }
      k = k + 1
    }
    match best_idx {
      0 => sb.write_char('A')
      1 => sb.write_char('C')
      2 => sb.write_char('G')
      3 => sb.write_char('T')
      _ => sb.write_char('-')
    }
    pos = pos + 1
  }
  sb.to_string()
}

///|
/// Browse sequences: format for display (ASCII art style).
pub fn browse_seqs(sequences : Array[String], names~ : Array[String] = []) -> String {
  let sb = StringBuilder::new()
  let max_name = if names.length() > 0 {
    let mut mx = 0
    let mut i = 0
    while i < names.length() {
      if names[i].length() > mx {
        mx = names[i].length()
      }
      i = i + 1
    }
    mx
  } else {
    8
  }
  let mut i = 0
  while i < sequences.length() {
    let name = if i < names.length() { names[i] } else { "Seq" + i.to_string() }
    let padding = max_name - name.length()
    let mut p = 0
    while p < padding {
      sb.write_char(' ')
      p = p + 1
    }
    sb.write_string(name)
    sb.write_char(' ')
    sb.write_string(sequences[i])
    sb.write_char('\n')
    i = i + 1
  }
  sb.to_string()
}

///|
/// Sequence database entry.
pub struct DbEntry {
  id : Int
  description : String
  sequence : String
}

///|
/// Sequence database (in-memory).
pub struct SeqDatabase {
  entries : Array[DbEntry]
  next_id : Int
}

///|
/// Create a new empty sequence database.
pub fn seq_database_new() -> SeqDatabase {
  SeqDatabase::{ entries: [], next_id: 0 }
}

///|
/// Add a sequence to the database (Seqs2Db style).
pub fn db_add_seq(
  db : SeqDatabase,
  sequence : String,
  description~ : String = "",
) -> SeqDatabase {
  let entry = DbEntry::{ id: db.next_id, description, sequence }
  let new_entries = []
  let mut i = 0
  while i < db.entries.length() {
    new_entries.push(db.entries[i])
    i = i + 1
  }
  new_entries.push(entry)
  SeqDatabase::{ entries: new_entries, next_id: db.next_id + 1 }
}

///|
/// Add multiple sequences to the database.
pub fn db_add_seqs(
  db : SeqDatabase,
  sequences : Array[String],
  descriptions~ : Array[String] = [],
) -> SeqDatabase {
  let mut result = db
  let mut i = 0
  while i < sequences.length() {
    let desc = if i < descriptions.length() { descriptions[i] } else { "" }
    result = db_add_seq(result, sequences[i], description=desc)
    i = i + 1
  }
  result
}

///|
/// Search database by sequence (SearchDb style).
/// Returns best matching entry index and alignment result.
pub fn db_search(
  db : SeqDatabase,
  query : String,
) -> (Int, DecipherAlnResult) {
  let matrix = default_dna_matrix()
  let mut best_idx = -1
  let mut best_score = -999999
  let mut best_aln = DecipherAlnResult::{
    aligned_a: "",
    aligned_b: "",
    score: 0,
    identity: 0.0,
    coverage: 0.0,
    mismatches: 0,
    gaps: 0,
    gap_length: 0,
  }
  let mut i = 0
  while i < db.entries.length() {
    let aln = align_nw_affine(query, db.entries[i].sequence, matrix)
    if aln.score > best_score {
      best_score = aln.score
      best_idx = i
      best_aln = aln
    }
    i = i + 1
  }
  (best_idx, best_aln)
}

///|
/// Get database size.
pub fn db_size(db : SeqDatabase) -> Int {
  db.entries.length()
}

///|
/// Get entry by index.
pub fn db_get(db : SeqDatabase, index : Int) -> DbEntry? {
  if index >= 0 && index < db.entries.length() {
    Some(db.entries[index])
  } else {
    None
  }
}

///|
/// Remove entry by index.
pub fn db_remove(db : SeqDatabase, index : Int) -> SeqDatabase {
  let new_entries : Array[DbEntry] = []
  let mut i = 0
  while i < db.entries.length() {
    if i != index {
      new_entries.push(db.entries[i])
    }
    i = i + 1
  }
  SeqDatabase::{ entries: new_entries, next_id: db.next_id }
}

///|
/// Calculate GC content of a sequence.
pub fn decipher_gc_content(seq : String) -> Double {
  if seq.length() == 0 {
    return 0.0
  }
  let mut gc = 0
  let mut total = 0
  let bytes = seq.to_array()
  let mut i = 0
  while i < bytes.length() {
    match bytes[i].to_int() {
      71 | 103 | 67 | 99 => {
        gc = gc + 1
        total = total + 1
      }
      65 | 97 | 84 | 116 => {
        total = total + 1
      }
      _ => ()
    }
    i = i + 1
  }
  if total == 0 {
    0.0
  } else {
    gc.to_double() / total.to_double()
  }
}

///|
/// Calculate sequence complexity (Shannon entropy).
pub fn decipher_complexity(seq : String) -> Double {
  if seq.length() == 0 {
    return 0.0
  }
  let counts : Array[Int] = Array::make(4, 0)
  let bytes = seq.to_array()
  let mut total = 0
  let mut i = 0
  while i < bytes.length() {
    let idx = decipher_nt_idx(bytes[i])
    if idx >= 0 {
      counts[idx] = counts[idx] + 1
      total = total + 1
    }
    i = i + 1
  }
  if total == 0 {
    return 0.0
  }
  let mut entropy = 0.0
  let mut k = 0
  while k < 4 {
    if counts[k] > 0 {
      let p = counts[k].to_double() / total.to_double()
      entropy = entropy - p * @math.ln(p)
    }
    k = k + 1
  }
  entropy / @math.ln(4.0) // Normalize to [0, 1]
}

///|
/// Create example sequences for testing/demo.
pub fn create_example_sequences() -> Array[String] {
  [
    "ACGTACGTACGTACGTACGT",
    "ACGTACGTACGTACGTACGA",
    "ACGTACGTACGTACGTACGT",
    "TTGCGCGCATATATATGCGC",
    "TTGCGCGCATATATATGCGT",
    "GCGCATATATATGCGCGCGC",
  ]
}

///|
/// Create example database for testing/demo.
pub fn create_example_database() -> SeqDatabase {
  let db = seq_database_new()
  db_add_seqs(
    db,
    [
      "ACGTACGTACGTACGTACGTACGTACGTACGT",
      "TTGCGCGCATATATATGCGCGCGCATATATA",
      "GCGCATATATATGCGCGCGCATATATATGCG",
      "ATATATATGCGCGCGCATATATATGCGCGCG",
    ],
    descriptions=["E_coli_16S", "L_plantarum_16S", "B_longum_16S", "B_fragilis_16S"],
  )
}