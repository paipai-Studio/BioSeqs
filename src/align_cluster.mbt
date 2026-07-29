///|
/// Progressive Multiple Sequence Alignment (ClustalW-like algorithm).
///
/// Implements progressive alignment via:
/// 1. Computing all pairwise distances using Needleman-Wunsch alignment
/// 2. Building a UPGMA guide tree from the distance matrix
/// 3. Traversing the guide tree from leaves to root
/// 4. At each internal node, aligning the two child profiles
/// 5. Returning the final multiple sequence alignment
///
/// References:
/// - Thompson JD, Higgins DG, Gibson TJ (1994). "CLUSTAL W: improving the
///   sensitivity of progressive multiple sequence alignment through
///   sequence weighting, position-specific gap penalties and weight matrix
///   choice." Nucleic Acids Res. 22(22):4673-80.
/// - Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis.
///   Cambridge University Press, 1998.

///|
/// Number of character categories in the profile frequency matrix.
/// Indices 0-25: A-Z, Index 26: gap '-'.
const PROFILE_ALPHABET_SIZE : Int = 27

///|
/// Map a UInt16 character code to its profile frequency index.
fn char_to_profile_index(c : UInt16) -> Int {
  if c == 45 {
    26
  } else {
    let upper = if c >= 97 && c <= 122 { c - 32 } else { c }
    if upper >= 65 && upper <= 90 {
      (upper - 65).to_int()
    } else {
      26
    }
  }
}

///|
/// Convert a profile index back to a UInt16 character code.
fn profile_index_to_uint16(idx : Int) -> UInt16 {
  if idx >= 0 && idx <= 25 {
    (idx + 65).to_uint16()
  } else {
    45
  }
}

///|
/// Progressive alignment clusterer with customizable scoring parameters.
pub struct AlignClusterer {
  match_score : Double
  mismatch_score : Double
  open_gap_penalty : Double
  extend_gap_penalty : Double
} derive(Eq)

///|
/// Create an AlignClusterer with default scoring parameters.
/// Defaults: match=2.0, mismatch=-1.0, open_gap=-2.0, extend_gap=-0.5.
pub fn align_cluster_new() -> AlignClusterer {
  AlignClusterer::{
    match_score: 2.0,
    mismatch_score: -1.0,
    open_gap_penalty: -2.0,
    extend_gap_penalty: -0.5,
  }
}

///|
/// Create an AlignClusterer with custom scoring parameters.
///
/// `match_score` — score for a matching character pair (positive).
/// `mismatch_score` — score for a mismatching character pair (negative).
/// `open_gap_penalty` — penalty for opening a new gap (negative).
/// `extend_gap_penalty` — penalty for extending an existing gap (negative).
pub fn align_cluster_with_scoring(
  match_score : Double,
  mismatch_score : Double,
  open_gap_penalty : Double,
  extend_gap_penalty : Double,
) -> AlignClusterer {
  AlignClusterer::{
    match_score,
    mismatch_score,
    open_gap_penalty,
    extend_gap_penalty,
  }
}

///|
/// Distance matrix for storing pairwise distances between sequences.
pub struct MSADistanceMatrix {
  matrix : Array[Array[Double]]
  sequence_ids : Array[String]
} derive(Eq)

///|
/// Create a new MSADistanceMatrix.
pub fn MSADistanceMatrix::new(
  matrix : Array[Array[Double]],
  sequence_ids : Array[String],
) -> MSADistanceMatrix {
  MSADistanceMatrix::{ matrix, sequence_ids }
}

///|
/// Get the distance matrix as a 2D array.
pub fn MSADistanceMatrix::get_matrix(self : MSADistanceMatrix) -> Array[Array[Double]] {
  self.matrix
}

///|
/// Get the sequence IDs.
pub fn MSADistanceMatrix::get_sequence_ids(self : MSADistanceMatrix) -> Array[String] {
  self.sequence_ids
}

///|
/// Get the distance between two sequences by index.
pub fn MSADistanceMatrix::get(self : MSADistanceMatrix, i : Int, j : Int) -> Double {
  self.matrix[i][j]
}

///|
/// A node in the guide tree.
pub struct GuideTreeNode {
  node_id : Int
  children : Array[Int]
  height : Double
  is_leaf : Bool
  sequence_index : Int
} derive(Eq)

///|
/// Create a leaf guide tree node.
fn guide_tree_leaf(node_id : Int, sequence_index : Int) -> GuideTreeNode {
  GuideTreeNode::{
    node_id,
    children: [],
    height: 0.0,
    is_leaf: true,
    sequence_index,
  }
}

///|
/// Create an internal guide tree node.
fn guide_tree_internal(
  node_id : Int,
  children : Array[Int],
  height : Double,
) -> GuideTreeNode {
  GuideTreeNode::{
    node_id,
    children,
    height,
    is_leaf: false,
    sequence_index: -1,
  }
}

///|
/// Guide tree for progressive alignment.
pub struct GuideTree {
  nodes : Array[GuideTreeNode]
} derive(Eq)

///|
/// Create a new GuideTree from an array of nodes.
pub fn GuideTree::new(nodes : Array[GuideTreeNode]) -> GuideTree {
  GuideTree::{ nodes }
}

///|
/// Get the nodes of the guide tree.
pub fn GuideTree::get_nodes(self : GuideTree) -> Array[GuideTreeNode] {
  self.nodes
}

///|
/// Get the root node (last node in the array).
pub fn GuideTree::get_root(self : GuideTree) -> GuideTreeNode {
  let n = self.nodes.length()
  if n > 0 {
    self.nodes[n - 1]
  } else {
    self.nodes[0]
  }
}

///|
/// Compute pairwise distance between two sequences using alignment.
///
/// Returns a distance value between 0.0 (identical) and 1.0 (completely different).
/// The distance is derived from the alignment score normalized by sequence length.
pub fn compute_distance(
  clusterer : AlignClusterer,
  seq1 : String,
  seq2 : String,
) -> Double {
  let (_, _, score) = align_pairwise_with_scoring(
    seq1,
    seq2,
    clusterer.match_score,
    clusterer.mismatch_score,
    clusterer.open_gap_penalty,
    clusterer.extend_gap_penalty,
  )
  let len1 = seq1.length()
  let len2 = seq2.length()
  let min_len = if len1 < len2 { len1 } else { len2 }
  if min_len == 0 {
    return 1.0
  }
  let sim = score / (min_len.to_double() * clusterer.match_score)
  if sim > 1.0 {
    0.0
  } else if sim < 0.0 {
    1.0
  } else {
    1.0 - sim
  }
}

///|
/// Build a distance matrix from an array of sequences.
///
/// `clusterer` — AlignClusterer with scoring parameters.
/// `sequences` — array of sequence strings.
/// `sequence_ids` — array of sequence identifiers.
///
/// Returns a MSADistanceMatrix containing all pairwise distances.
pub fn build_distance_matrix(
  clusterer : AlignClusterer,
  sequences : Array[String],
  sequence_ids : Array[String],
) -> MSADistanceMatrix {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else {
        let d = compute_distance(clusterer, sequences[i], sequences[j])
        row.push(d)
      }
    }
    matrix.push(row)
  }
  MSADistanceMatrix::{ matrix, sequence_ids }
}

///|
/// Build a UPGMA guide tree from a distance matrix.
///
/// UPGMA (Unweighted Pair Group Method with Arithmetic mean) builds a
/// hierarchical clustering tree by iteratively merging the two closest
/// clusters.
///
/// Returns a GuideTree with leaf nodes for each sequence and internal
/// nodes representing merge events.
pub fn upgma_tree(distance_matrix : MSADistanceMatrix) -> GuideTree {
  let n = distance_matrix.matrix.length()
  if n == 0 {
    return GuideTree::{ nodes: [] }
  }
  if n == 1 {
    let node = guide_tree_leaf(0, 0)
    return GuideTree::{ nodes: [node] }
  }

  let nodes : Array[GuideTreeNode] = Array::new()
  for i = 0; i < n; i = i + 1 {
    nodes.push(guide_tree_leaf(i, i))
  }

  let mut dist = copy_distance_matrix(distance_matrix.matrix)
  let mut active : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    active.push(i)
  }
  let mut next_id = n

  while active.length() > 1 {
    let m = active.length()

    let mut min_d = 100000000000000000000000000000.0
    let mut min_i = 0
    let mut min_j = 1
    for i = 0; i < m; i = i + 1 {
      for j = i + 1; j < m; j = j + 1 {
        let d = dist[i][j]
        if d < min_d {
          min_d = d
          min_i = i
          min_j = j
        }
      }
    }

    let i_idx = active[min_i]
    let j_idx = active[min_j]
    let height = min_d / 2.0

    let size_i = count_leaves_in_node(nodes, i_idx)
    let size_j = count_leaves_in_node(nodes, j_idx)

    let new_node = guide_tree_internal(
      next_id,
      [i_idx, j_idx],
      height,
    )
    nodes.push(new_node)
    let new_node_idx = next_id
    next_id = next_id + 1

    let new_dists_row : Array[Double] = Array::new()
    for k = 0; k < m; k = k + 1 {
      if k != min_i && k != min_j {
        let d = (
          size_i.to_double() * dist[min_i][k] +
          size_j.to_double() * dist[min_j][k]
        ) / (size_i.to_double() + size_j.to_double())
        new_dists_row.push(d)
      }
    }
    new_dists_row.push(0.0)

    let new_active : Array[Int] = Array::new()
    for k = 0; k < m; k = k + 1 {
      if k != min_i && k != min_j {
        new_active.push(active[k])
      }
    }
    new_active.push(new_node_idx)

    let new_dist : Array[Array[Double]] = Array::new()
    for r = 0; r < m; r = r + 1 {
      if r != min_i && r != min_j {
        let row : Array[Double] = Array::new()
        for c = 0; c < m; c = c + 1 {
          if c != min_i && c != min_j {
            row.push(dist[r][c])
          }
        }
        row.push(new_dists_row[new_dist.length()])
        new_dist.push(row)
      }
    }
    new_dist.push(new_dists_row)

    active = new_active
    dist = new_dist
  }

  GuideTree::{ nodes }
}

///|
/// Count the number of leaf nodes under a given node.
fn count_leaves_in_node(
  nodes : Array[GuideTreeNode],
  node_idx : Int,
) -> Int {
  let node = nodes[node_idx]
  if node.is_leaf {
    1
  } else {
    let mut count = 0
    for child_idx in node.children {
      count = count + count_leaves_in_node(nodes, child_idx)
    }
    count
  }
}

///|
/// Copy a distance matrix.
fn copy_distance_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for row in matrix {
    let new_row : Array[Double] = Array::new()
    for v in row {
      new_row.push(v)
    }
    result.push(new_row)
  }
  result
}

///|
/// Align two sequences using Needleman-Wunsch global alignment.
///
/// `seq1` — first sequence string.
/// `seq2` — second sequence string.
///
/// Returns a tuple of (aligned_seq1, aligned_seq2, alignment_score).
pub fn align_pairwise(
  seq1 : String,
  seq2 : String,
) -> (String, String, Double) {
  align_pairwise_with_scoring(
    seq1,
    seq2,
    2.0,
    -1.0,
    -2.0,
    -0.5,
  )
}

///|
/// Align two sequences with custom scoring using NW algorithm.
///
/// Uses affine gap penalties (open vs extend).
pub fn align_pairwise_with_scoring(
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  open_gap : Double,
  extend_gap : Double,
) -> (String, String, Double) {
  let n = seq1.length()
  let m = seq2.length()

  if n == 0 && m == 0 {
    return ("", "", 0.0)
  }
  if n == 0 {
    let gap_seq = repeat_char('-', m)
    return (gap_seq, seq2, extend_gap * m.to_double())
  }
  if m == 0 {
    let gap_seq = repeat_char('-', n)
    return (seq1, gap_seq, extend_gap * n.to_double())
  }

  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[Int]] = Array::new()
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array::new()
    let trow : Array[Int] = Array::new()
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(0)
    }
    dp.push(drow)
    trace.push(trow)
  }

  dp[0][0] = 0.0
  trace[0][0] = 0
  for i = 1; i <= n; i = i + 1 {
    dp[i][0] = open_gap + (i - 1).to_double() * extend_gap
    trace[i][0] = 2
  }
  for j = 1; j <= m; j = j + 1 {
    dp[0][j] = open_gap + (j - 1).to_double() * extend_gap
    trace[0][j] = 1
  }

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = seq1.unsafe_get(i - 1)
      let c2 = seq2.unsafe_get(j - 1)
      let sub_score = if c1 == c2 {
        match_score
      } else {
        mismatch_score
      }

      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] +
        (if trace[i - 1][j] == 2 { extend_gap } else { open_gap })
      let left = dp[i][j - 1] +
        (if trace[i][j - 1] == 1 { extend_gap } else { open_gap })

      let (val, dir) = if diag >= up && diag >= left {
        (diag, 0)
      } else if up >= left {
        (up, 2)
      } else {
        (left, 1)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }

  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = n
  let mut j = m
  while i > 0 || j > 0 {
    if i > 0 && j > 0 && trace[i][j] == 0 {
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      i = i - 1
      j = j - 1
    } else if i > 0 && trace[i][j] == 2 {
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if j > 0 && trace[i][j] == 1 {
      r1.write_char('-')
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      j = j - 1
    } else if i > 0 {
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if j > 0 {
      r1.write_char('-')
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      j = j - 1
    } else {
      break
    }
  }

  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let score = dp[n][m]
  (aligned1, aligned2, score)
}

///|
/// Build a position weight matrix (frequency profile) from a set of aligned sequences.
///
/// `sequences` — array of aligned sequences (all same length, may contain gaps).
///
/// Returns a 2D array where `result[col][idx]` is the frequency of the character
/// at profile index `idx` (0-25: A-Z, 26: gap '-') at column `col`.
pub fn profile_from_sequences(
  sequences : Array[String],
) -> Array[Array[Double]] {
  if sequences.length() == 0 {
    return []
  }
  let n_cols = sequences[0].length()
  let n_seqs = sequences.length()
  let profile : Array[Array[Double]] = Array::new()

  for col = 0; col < n_cols; col = col + 1 {
    let col_freq : Array[Double] = Array::new()
    for k = 0; k < PROFILE_ALPHABET_SIZE; k = k + 1 {
      col_freq.push(0.0)
    }
    for seq_idx = 0; seq_idx < n_seqs; seq_idx = seq_idx + 1 {
      let c = sequences[seq_idx].unsafe_get(col)
      let idx = char_to_profile_index(c)
      col_freq[idx] = col_freq[idx] + 1.0
    }
    for k = 0; k < PROFILE_ALPHABET_SIZE; k = k + 1 {
      col_freq[k] = col_freq[k] / n_seqs.to_double()
    }
    profile.push(col_freq)
  }
  profile
}

///|
/// Score the alignment of two profile columns using dot product of frequencies.
///
/// `profile1_col` — frequency array for column from first profile.
/// `profile2_col` — frequency array for column from second profile.
///
/// Returns a score representing the similarity of the two columns.
/// A high score indicates similar character distributions.
pub fn score_profile_position(
  profile1_col : Array[Double],
  profile2_col : Array[Double],
) -> Double {
  let len = if profile1_col.length() < profile2_col.length() {
    profile1_col.length()
  } else {
    profile2_col.length()
  }
  let mut score = 0.0
  for k = 0; k < len; k = k + 1 {
    score = score + profile1_col[k] * profile2_col[k]
  }
  score
}

///|
/// Align two sets of aligned sequences (profiles) using profile-profile alignment.
///
/// `profile1` — first set of already-aligned sequences.
/// `profile2` — second set of already-aligned sequences.
///
/// Returns a tuple of (aligned_consensus1, aligned_consensus2, score).
/// The aligned sequences are consensus sequences representing the alignment.
pub fn align_profiles(
  profile1 : Array[String],
  profile2 : Array[String],
) -> (String, String, Double) {
  if profile1.length() == 0 && profile2.length() == 0 {
    return ("", "", 0.0)
  }
  if profile1.length() == 0 {
    let rep2 = if profile2.length() > 0 { profile2[0] } else { "" }
    let gap1 = repeat_char('-', rep2.length())
    return (gap1, rep2, 0.0)
  }
  if profile2.length() == 0 {
    let rep1 = if profile1.length() > 0 { profile1[0] } else { "" }
    let gap2 = repeat_char('-', rep1.length())
    return (rep1, gap2, 0.0)
  }

  let profile1_matrix = profile_from_sequences(profile1)
  let profile2_matrix = profile_from_sequences(profile2)

  let len1 = profile1_matrix.length()
  let len2 = profile2_matrix.length()

  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[Int]] = Array::new()
  for i = 0; i <= len1; i = i + 1 {
    let drow : Array[Double] = Array::new()
    let trow : Array[Int] = Array::new()
    for j = 0; j <= len2; j = j + 1 {
      drow.push(0.0)
      trow.push(0)
    }
    dp.push(drow)
    trace.push(trow)
  }

  let gap_penalty = -1.0
  dp[0][0] = 0.0
  trace[0][0] = 0
  for i = 1; i <= len1; i = i + 1 {
    dp[i][0] = i.to_double() * gap_penalty
    trace[i][0] = 2
  }
  for j = 1; j <= len2; j = j + 1 {
    dp[0][j] = j.to_double() * gap_penalty
    trace[0][j] = 1
  }

  for i = 1; i <= len1; i = i + 1 {
    for j = 1; j <= len2; j = j + 1 {
      let col1 = profile1_matrix[i - 1]
      let col2 = profile2_matrix[j - 1]
      let sub_score = score_profile_position(col1, col2) * 2.0 - 1.0

      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + gap_penalty
      let left = dp[i][j - 1] + gap_penalty

      let (val, dir) = if diag >= up && diag >= left {
        (diag, 0)
      } else if up >= left {
        (up, 2)
      } else {
        (left, 1)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }

  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = len1
  let mut j = len2
  while i > 0 || j > 0 {
    if i > 0 && j > 0 && trace[i][j] == 0 {
      let col1 = profile1_matrix[i - 1]
      let col2 = profile2_matrix[j - 1]
      let consensus_idx1 = argmax(col1)
      let consensus_idx2 = argmax(col2)
      r1.write_char(profile_index_to_uint16(consensus_idx1).unsafe_to_char())
      r2.write_char(profile_index_to_uint16(consensus_idx2).unsafe_to_char())
      i = i - 1
      j = j - 1
    } else if i > 0 && trace[i][j] == 2 {
      let col1 = profile1_matrix[i - 1]
      let consensus_idx = argmax(col1)
      r1.write_char(profile_index_to_uint16(consensus_idx).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if j > 0 && trace[i][j] == 1 {
      let col2 = profile2_matrix[j - 1]
      let consensus_idx = argmax(col2)
      r1.write_char('-')
      r2.write_char(profile_index_to_uint16(consensus_idx).unsafe_to_char())
      j = j - 1
    } else if i > 0 {
      let col1 = profile1_matrix[i - 1]
      let consensus_idx = argmax(col1)
      r1.write_char(profile_index_to_uint16(consensus_idx).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if j > 0 {
      let col2 = profile2_matrix[j - 1]
      let consensus_idx = argmax(col2)
      r1.write_char('-')
      r2.write_char(profile_index_to_uint16(consensus_idx).unsafe_to_char())
      j = j - 1
    } else {
      break
    }
  }

  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let score = dp[len1][len2]
  (aligned1, aligned2, score)
}

///|
/// Return the index of the maximum value in an array.
fn argmax(arr : Array[Double]) -> Int {
  if arr.length() == 0 {
    return 0
  }
  let mut max_idx = 0
  let mut max_val = arr[0]
  for i = 1; i < arr.length(); i = i + 1 {
    if arr[i] > max_val {
      max_val = arr[i]
      max_idx = i
    }
  }
  max_idx
}

///|
/// Progressive multiple sequence alignment.
///
/// `clusterer` — AlignClusterer with scoring parameters.
/// `sequences` — array of sequence strings to align.
/// `sequence_ids` — array of sequence identifiers.
///
/// Returns a MultipleSeqAlignment containing the aligned sequences.
///
/// Algorithm:
/// 1. Compute all pairwise distances
/// 2. Build UPGMA guide tree
/// 3. Follow guide tree from leaves to root
/// 4. At each internal node, align the two child profiles
/// 5. Return final alignment
pub fn progressive_alignment(
  clusterer : AlignClusterer,
  sequences : Array[String],
  sequence_ids : Array[String],
) -> MultipleSeqAlignment {
  let n = sequences.length()
  if n == 0 {
    return create_msa_safe([])
  }
  if n == 1 {
    let rec = SeqRecord::new(Seq::new(sequences[0]), id=sequence_ids[0])
    return create_msa_safe([rec])
  }

  let dist_matrix = build_distance_matrix(clusterer, sequences, sequence_ids)
  let tree = upgma_tree(dist_matrix)

  let alignment_map : Map[Int, Array[String]] = Map([], capacity=n * 2)
  for i = 0; i < n; i = i + 1 {
    alignment_map[i] = [sequences[i]]
  }

  let nodes = tree.nodes
  let processed : Map[Int, Bool] = Map([], capacity=n * 2)

  let sorted_nodes = topological_sort_nodes(nodes)
  for node_idx in sorted_nodes {
    let node = nodes[node_idx]
    if node.is_leaf {
      processed[node_idx] = true
    } else {
      let child_results : Array[Array[String]] = Array::new()
      let mut all_processed = true
      for child_idx in node.children {
        match processed.get(child_idx) {
          Some(true) => {
            match alignment_map.get(child_idx) {
              Some(profile) => child_results.push(profile)
              None => all_processed = false
            }
          }
          _ => all_processed = false
        }
      }
      if all_processed && node.children.length() >= 2 {
        let left_profile = child_results[0]
        let right_profile = child_results[1]
        let (aligned_left, aligned_right, _) = align_profiles(
          left_profile,
          right_profile,
        )
        let left_aligned_profile : Array[String] = Array::new()
        let right_aligned_profile : Array[String] = Array::new()
        for s in left_profile {
          let new_s = apply_alignment_to_sequence(s, aligned_left)
          left_aligned_profile.push(new_s)
        }
        for s in right_profile {
          let new_s = apply_alignment_to_sequence(s, aligned_right)
          right_aligned_profile.push(new_s)
        }
        let combined : Array[String] = Array::new()
        for s in left_aligned_profile {
          combined.push(s)
        }
        for s in right_aligned_profile {
          combined.push(s)
        }
        alignment_map[node_idx] = combined
        processed[node_idx] = true
      } else if node.children.length() == 1 {
        match alignment_map.get(node.children[0]) {
          Some(profile) => {
            alignment_map[node_idx] = profile
            processed[node_idx] = true
          }
          None => ()
        }
      }
    }
  }

  let root_idx = nodes.length() - 1
  let final_profiles = match alignment_map.get(root_idx) {
    Some(p) => p
    None => {
      let profiles : Array[Array[String]] = Array::new()
      for i = 0; i < n; i = i + 1 {
        profiles.push([sequences[i]])
      }
      let mut result = profiles[0]
      for i = 1; i < n; i = i + 1 {
        let (a1, a2, _) = align_profiles(result, profiles[i])
        let new_result : Array[String] = Array::new()
        for s in result {
          new_result.push(apply_alignment_to_sequence(s, a1))
        }
        let new_i : Array[String] = Array::new()
        for s in profiles[i] {
          new_i.push(apply_alignment_to_sequence(s, a2))
        }
        result = new_result
        for s in new_i {
          result.push(s)
        }
      }
      result
    }
  }

  let records : Array[SeqRecord] = Array::new()
  let mut seq_id_idx = 0
  for s in final_profiles {
    let id = if seq_id_idx < sequence_ids.length() {
      sequence_ids[seq_id_idx]
    } else {
      "seq_" + seq_id_idx.to_string()
    }
    records.push(SeqRecord::new(Seq::new(s), id=id))
    seq_id_idx = seq_id_idx + 1
  }
  create_msa_safe(records)
}

///|
/// Create a MultipleSeqAlignment safely, handling AlignIOError.
fn create_msa_safe(records : Array[SeqRecord]) -> MultipleSeqAlignment {
  let empty : MultipleSeqAlignment = MultipleSeqAlignment::{
    records: [],
    annotations: Map([], capacity=0),
    column_annotations: Map([], capacity=0),
  }
  try {
    MultipleSeqAlignment::new(records)
  } catch {
    _ => empty
  }
}

///|
/// Topologically sort guide tree nodes (leaves first, root last).
fn topological_sort_nodes(nodes : Array[GuideTreeNode]) -> Array[Int] {
  let n = nodes.length()
  let visited : Array[Bool] = Array::new()
  for i = 0; i < n; i = i + 1 {
    visited.push(false)
  }
  let result : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    if !visited[i] {
      topological_sort_visit(nodes, i, visited, result)
    }
  }
  result
}

///|
/// DFS helper for topological sorting.
fn topological_sort_visit(
  nodes : Array[GuideTreeNode],
  idx : Int,
  visited : Array[Bool],
  result : Array[Int],
) -> Unit {
  visited[idx] = true
  let node = nodes[idx]
  for child_idx in node.children {
    if !visited[child_idx] {
      topological_sort_visit(nodes, child_idx, visited, result)
    }
  }
  result.push(idx)
}

///|
/// Apply an alignment trace to a sequence, inserting gaps where needed.
///
/// `original_seq` — the original sequence string.
/// `aligned_representative` — the aligned version of the representative,
///   which may contain gap characters.
///
/// Returns the original sequence with gaps inserted at the same positions
/// as in the aligned representative.
fn apply_alignment_to_sequence(
  original_seq : String,
  aligned_representative : String,
) -> String {
  let orig_len = original_seq.length()
  let aligned_len = aligned_representative.length()
  if aligned_len == 0 {
    return ""
  }
  let result = StringBuilder::new()
  let mut orig_pos = 0
  for i = 0; i < aligned_len; i = i + 1 {
    let c = aligned_representative.unsafe_get(i)
    if c == 45 {
      result.write_char('-')
    } else if orig_pos < orig_len {
      result.write_char(original_seq.unsafe_get(orig_pos).unsafe_to_char())
      orig_pos = orig_pos + 1
    } else {
      result.write_char('-')
    }
  }
  result.to_string()
}

///|
/// Get aligned sequences only (array of strings) from progressive alignment.
///
/// `clusterer` — AlignClusterer with scoring parameters.
/// `sequences` — array of sequence strings to align.
///
/// Returns an array of aligned sequence strings.
pub fn get_aligned_sequences(
  clusterer : AlignClusterer,
  sequences : Array[String],
) -> Array[String] {
  let ids : Array[String] = Array::new()
  for i = 0; i < sequences.length(); i = i + 1 {
    ids.push("seq_" + i.to_string())
  }
  let alignment = progressive_alignment(clusterer, sequences, ids)
  let result : Array[String] = Array::new()
  for rec in alignment.iter() {
    result.push(rec.seq.to_string())
  }
  result
}

///|
/// Create a string of repeated characters.
fn repeat_char(c : Char, n : Int) -> String {
  let buf = StringBuilder::new(size_hint=n)
  for i = 0; i < n; i = i + 1 {
    buf.write_char(c)
  }
  buf.to_string()
}