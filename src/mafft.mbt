///|
/// MAFFT-lite: progressive multiple sequence alignment.
///
/// Implements the core pipeline of MAFFT (FFT-NS-2 style):
///   1. K-mer distance for fast pairwise distance estimation
///   2. Guide tree construction (UPGMA)
///   3. Progressive alignment (profile-profile NW with BLOSUM62)
///   4. Iterative refinement (realign each sequence to profile of rest)
///
/// This is a pure-MoonBit reimplementation without FFT optimization.
/// Algorithm follows Katoh & Standley, "MAFFT multiple sequence alignment
/// software version 7" (Mol Biol Evol 2013).

// ============================================================================
// 1. K-mer distance
// ============================================================================

///|
/// Count k-mer occurrences in a sequence.
fn kmer_counts(seq : String, k : Int) -> Map[Int, Int] {
  let n = seq.length()
  let counts : Map[Int, Int] = Map([], capacity=0)
  if n < k {
    return counts
  }
  for i = 0; i <= n - k; i = i + 1 {
    let mut h = 0
    let mut valid = true
    for j = 0; j < k; j = j + 1 {
      let c = seq[i + j].unsafe_to_char()
      let v = match c {
        'A' | 'a' => 0
        'C' | 'c' => 1
        'G' | 'g' => 2
        'T' | 't' | 'U' | 'u' => 3
        _ => {
          valid = false
          0
        }
      }
      h = (h << 2) | v
    }
    if valid {
      match counts.get(h) {
        Some(v) => counts[h] = v + 1
        None => counts[h] = 1
      }
    }
  }
  counts
}

///|
/// K-mer distance between two sequences.
/// Distance = 1 - (shared k-mers / total unique k-mers).
pub fn mafft_kmer_distance(seq1 : String, seq2 : String, k : Int) -> Double {
  let c1 = kmer_counts(seq1, k)
  let c2 = kmer_counts(seq2, k)
  let mut shared = 0
  let mut total = 0
  for key, v1 in c1 {
    total += v1
    match c2.get(key) {
      Some(v2) => {
        if v1 < v2 {
          shared += v1
        } else {
          shared += v2
        }
      }
      None => ()
    }
  }
  for key, v2 in c2 {
    match c1.get(key) {
      Some(_) => ()
      None => total += v2
    }
  }
  if total == 0 {
    return 1.0
  }
  1.0 - shared.to_double() / total.to_double()
}

///|
/// Build a k-mer distance matrix for a set of sequences.
pub fn mafft_distance_matrix(
  sequences : Array[String],
  k : Int,
) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row[j] = 0.0
      } else if j > i {
        row[j] = mafft_kmer_distance(sequences[i], sequences[j], k)
      } else {
        row[j] = matrix[j][i]
      }
    }
    matrix.push(row)
  }
  matrix
}

// ============================================================================
// 2. Guide tree (UPGMA)
// ============================================================================

///|
/// Guide tree node for MSA.
pub struct MsaTreeNode {
  pub children : Array[MsaTreeNode]
  pub is_leaf : Bool
  pub seq_index : Int
  pub height : Double
}

///|
/// Build a UPGMA guide tree from a distance matrix.
pub fn mafft_upgma_tree(
  distances : Array[Array[Double]],
) -> MsaTreeNode {
  let n = distances.length()
  if n == 0 {
    return { children: [], is_leaf: true, seq_index: -1, height: 0.0 }
  }
  if n == 1 {
    return { children: [], is_leaf: true, seq_index: 0, height: 0.0 }
  }
  let clusters : Array[MsaTreeNode] = []
  let weights : Array[Int] = Array::make(n, 1)
  for i = 0; i < n; i = i + 1 {
    clusters.push({ children: [], is_leaf: true, seq_index: i, height: 0.0 })
  }
  let dist : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = distances[i][j]
    }
    dist.push(row)
  }
  let active : Array[Int] = Array::makei(n, fn(i) { i })
  while active.length() > 1 {
    let mut min_dist = 1.0e18
    let mut min_i = 0
    let mut min_j = 1
    for ai = 0; ai < active.length(); ai = ai + 1 {
      for aj = ai + 1; aj < active.length(); aj = aj + 1 {
        let i = active[ai]
        let j = active[aj]
        if dist[i][j] < min_dist {
          min_dist = dist[i][j]
          min_i = ai
          min_j = aj
        }
      }
    }
    let ci = active[min_i]
    let cj = active[min_j]
    let new_height = min_dist / 2.0
    let merged = {
      children: [clusters[ci], clusters[cj]],
      is_leaf: false,
      seq_index: -1,
      height: new_height,
    }
    let wi = weights[ci].to_double()
    let wj = weights[cj].to_double()
    let new_idx = clusters.length()
    clusters.push(merged)
    weights.push(weights[ci] + weights[cj])
    let new_row : Array[Double] = Array::make(clusters.length(), 0.0)
    for k = 0; k < clusters.length(); k = k + 1 {
      if k == new_idx {
        new_row[k] = 0.0
      } else {
        let d_ik = if k < dist[ci].length() { dist[ci][k] } else { 0.0 }
        let d_jk = if k < dist[cj].length() { dist[cj][k] } else { 0.0 }
        new_row[k] = (wi * d_ik + wj * d_jk) / (wi + wj)
      }
    }
    for i = 0; i < dist.length(); i = i + 1 {
      while dist[i].length() < clusters.length() {
        dist[i].push(0.0)
      }
      dist[i][new_idx] = new_row[i]
    }
    dist.push(new_row)
    let _ = active.remove(min_j)
    let _ = active.remove(min_i)
    active.push(new_idx)
  }
  clusters[active[0]]
}

// ============================================================================
// 3. Profile-profile alignment
// ============================================================================

///|
/// Score for aligning two columns (characters, possibly gaps).
fn mafft_col_score(c1 : Char, c2 : Char, matrix : SubstitutionMatrix) -> Double {
  if c1 == '-' && c2 == '-' {
    return 0.0
  }
  if c1 == '-' || c2 == '-' {
    return -4.0
  }
  matrix.get_score_case_insensitive(c1.to_string(), c2.to_string()).to_double()
}

///|
/// Align two profiles using NW DP with substitution matrix.
pub fn mafft_align_profiles(
  profile1 : Array[String],
  profile2 : Array[String],
  matrix : SubstitutionMatrix,
  gap_open : Double,
) -> Array[String] {
  if profile1.length() == 0 {
    return profile2
  }
  if profile2.length() == 0 {
    return profile1
  }
  let m = profile1[0].length()
  let n = profile2[0].length()
  if m == 0 {
    let gap_seq = "-".repeat(n)
    let result : Array[String] = []
    for s in profile1 {
      result.push(gap_seq)
    }
    for s in profile2 {
      result.push(s)
    }
    return result
  }
  if n == 0 {
    let gap_seq = "-".repeat(m)
    let result : Array[String] = []
    for s in profile1 {
      result.push(s)
    }
    for s in profile2 {
      result.push(gap_seq)
    }
    return result
  }
  // Profile-profile DP: score match by summing BLOSUM scores across all
  // non-gap pairs between profile1 column i and profile2 column j
  let dp : Array[Array[Double]] = Array::makei(m + 1, fn(_i) {
    Array::make(n + 1, 0.0)
  })
  let trace : Array[Array[Int]] = Array::makei(m + 1, fn(_i) {
    Array::make(n + 1, 0)
  })
  let n1 = profile1.length()
  let n2 = profile2.length()
  for i = 1; i <= m; i = i + 1 {
    dp[i][0] = dp[i - 1][0] + gap_open * n1.to_double()
    trace[i][0] = 1
  }
  for j = 1; j <= n; j = j + 1 {
    dp[0][j] = dp[0][j - 1] + gap_open * n2.to_double()
    trace[0][j] = 2
  }
  for i = 1; i <= m; i = i + 1 {
    for j = 1; j <= n; j = j + 1 {
      // Match score: sum of all pairwise scores between columns
      let mut match_score = 0.0
      for a = 0; a < n1; a = a + 1 {
        let ca = profile1[a][i - 1].unsafe_to_char()
        for b = 0; b < n2; b = b + 1 {
          let cb = profile2[b][j - 1].unsafe_to_char()
          if ca == '-' || cb == '-' {
            if ca == '-' && cb == '-' {
              // Both gaps, no penalty
              continue
            }
            match_score += gap_open
          } else {
            match_score += matrix.get_score_case_insensitive(
              ca.to_string(),
              cb.to_string(),
            ).to_double()
          }
        }
      }
      let match_val = dp[i - 1][j - 1] + match_score
      let up_val = dp[i - 1][j] + gap_open * n1.to_double()
      let left_val = dp[i][j - 1] + gap_open * n2.to_double()
      if match_val >= up_val && match_val >= left_val {
        dp[i][j] = match_val
        trace[i][j] = 0
      } else if up_val >= left_val {
        dp[i][j] = up_val
        trace[i][j] = 1
      } else {
        dp[i][j] = left_val
        trace[i][j] = 2
      }
    }
  }
  // Traceback: produce ops sequence (0=match, 1=up/gap-in-2, 2=left/gap-in-1)
  let ops : Array[Int] = Array(capacity=0)
  let mut i = m
  let mut j = n
  while i > 0 || j > 0 {
    if i > 0 && j > 0 && trace[i][j] == 0 {
      ops.push(0)
      i -= 1
      j -= 1
    } else if i > 0 && (j == 0 || trace[i][j] == 1) {
      ops.push(1)
      i -= 1
    } else {
      ops.push(2)
      j -= 1
    }
  }
  // Reverse ops
  let ops_len = ops.length()
  for k = 0; k < ops_len / 2; k = k + 1 {
    let tmp = ops[k]
    ops[k] = ops[ops_len - 1 - k]
    ops[ops_len - 1 - k] = tmp
  }
  // Apply ops to all sequences
  let result : Array[String] = []
  for s in profile1 {
    result.push(apply_profile_ops(s, ops, 1))
  }
  for s in profile2 {
    result.push(apply_profile_ops(s, ops, 2))
  }
  result
}

///|
/// Apply profile alignment ops to a sequence.
/// ops: 0=match (take from both), 1=up (take from seq1, gap in seq2), 2=left (gap in seq1, take from seq2)
/// side: 1 for seq from profile1, 2 for seq from profile2
fn apply_profile_ops(seq : String, ops : Array[Int], side : Int) -> String {
  let buf : StringBuilder = StringBuilder()
  let mut si = 0
  for k = 0; k < ops.length(); k = k + 1 {
    let op = ops[k]
    if op == 0 {
      // Match: take from both
      if si < seq.length() {
        buf.write_char(seq[si].unsafe_to_char())
      } else {
        buf.write_char('-')
      }
      si += 1
    } else if op == 1 {
      // Up: take from profile1, gap in profile2
      if side == 1 {
        if si < seq.length() {
          buf.write_char(seq[si].unsafe_to_char())
        } else {
          buf.write_char('-')
        }
        si += 1
      } else {
        buf.write_char('-')
      }
    } else {
      // Left: gap in profile1, take from profile2
      if side == 2 {
        if si < seq.length() {
          buf.write_char(seq[si].unsafe_to_char())
        } else {
          buf.write_char('-')
        }
        si += 1
      } else {
        buf.write_char('-')
      }
    }
  }
  buf.to_string()
}

// ============================================================================
// 4. Progressive alignment
// ============================================================================

///|
/// Progressive alignment following a guide tree.
fn progressive_align(
  tree : MsaTreeNode,
  sequences : Array[String],
  matrix : SubstitutionMatrix,
  gap_open : Double,
) -> (Array[String], Array[Int]) {
  if tree.is_leaf {
    return ([sequences[tree.seq_index]], [tree.seq_index])
  }
  let mut result : Array[String] = []
  let mut order : Array[Int] = []
  for child in tree.children {
    let (child_profile, child_order) = progressive_align(
      child,
      sequences,
      matrix,
      gap_open,
    )
    if result.length() == 0 {
      result = child_profile
      order = child_order
    } else {
      result = mafft_align_profiles(result, child_profile, matrix, gap_open)
      for idx in child_order {
        order.push(idx)
      }
    }
  }
  (result, order)
}

///|
/// Reorder alignment to match original input order.
fn reorder_alignment(
  alignment : Array[String],
  order : Array[Int],
) -> Array[String] {
  let n = alignment.length()
  let result : Array[String] = Array::make(n, "")
  for i = 0; i < n; i = i + 1 {
    result[order[i]] = alignment[i]
  }
  result
}

// ============================================================================
// 5. Iterative refinement
// ============================================================================

///|
/// Compute the SP (sum-of-pairs) score of an alignment.
fn sp_score(alignment : Array[String], matrix : SubstitutionMatrix) -> Double {
  if alignment.length() < 2 {
    return 0.0
  }
  let n = alignment.length()
  let m = alignment[0].length()
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      for k = 0; k < m; k = k + 1 {
        let c1 = alignment[i][k].unsafe_to_char()
        let c2 = alignment[j][k].unsafe_to_char()
        total += mafft_col_score(c1, c2, matrix)
      }
    }
  }
  total
}

///|
/// Remove gaps from a sequence.
fn remove_gaps(seq : String) -> String {
  let buf : StringBuilder = StringBuilder()
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i].unsafe_to_char()
    if c != '-' {
      buf.write_char(c)
    }
  }
  buf.to_string()
}

///|
/// Remove all-gap columns from an alignment.
fn remove_all_gap_columns(alignment : Array[String]) -> Array[String] {
  if alignment.length() == 0 {
    return []
  }
  let m = alignment[0].length()
  let keep : Array[Bool] = Array::make(m, true)
  for col = 0; col < m; col = col + 1 {
    let mut all_gap = true
    for row = 0; row < alignment.length(); row = row + 1 {
      if alignment[row][col].unsafe_to_char() != '-' {
        all_gap = false
        break
      }
    }
    keep[col] = !all_gap
  }
  let result : Array[String] = []
  for row = 0; row < alignment.length(); row = row + 1 {
    let buf : StringBuilder = StringBuilder()
    for col = 0; col < m; col = col + 1 {
      if keep[col] {
        buf.write_char(alignment[row][col].unsafe_to_char())
      }
    }
    result.push(buf.to_string())
  }
  result
}

///|
/// One round of iterative refinement.
fn refine_once(
  alignment : Array[String],
  matrix : SubstitutionMatrix,
  gap_open : Double,
) -> Array[String] {
  let n = alignment.length()
  if n <= 2 {
    return alignment
  }
  let mut best = alignment
  let mut best_score = sp_score(best, matrix)
  for i = 0; i < n; i = i + 1 {
    let seq_i = remove_gaps(alignment[i])
    let remaining : Array[String] = []
    for j = 0; j < n; j = j + 1 {
      if j != i {
        remaining.push(alignment[j])
      }
    }
    let realigned = mafft_align_profiles(remaining, [seq_i], matrix, gap_open)
    let cleaned = remove_all_gap_columns(realigned)
    let new_score = sp_score(cleaned, matrix)
    if new_score > best_score {
      best = cleaned
      best_score = new_score
    }
  }
  best
}

// ============================================================================
// 6. Public API
// ============================================================================

///|
/// MSA result: aligned sequences, names, and score.
pub struct MsaResult {
  pub aligned_sequences : Array[String]
  pub names : Array[String]
  pub score : Double
  pub method : String
}

///|
/// MAFFT-lite alignment (FFT-NS-2 style).
/// Uses k-mer distance, UPGMA guide tree, progressive alignment,
/// and one round of iterative refinement.
pub fn mafft_align(
  sequences : Array[String],
  names : Array[String],
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "mafft" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "mafft" }
  }
  let matrix = blosum62()
  let gap_open = -10.0
  let k = 6
  // Step 1: k-mer distance
  let dist = mafft_distance_matrix(sequences, k)
  // Step 2: UPGMA guide tree
  let tree = mafft_upgma_tree(dist)
  // Step 3: progressive alignment
  let (alignment_raw, order) = progressive_align(tree, sequences, matrix, gap_open)
  let alignment = reorder_alignment(alignment_raw, order)
  // Step 4: iterative refinement
  let refined = refine_once(alignment, matrix, gap_open)
  let cleaned = remove_all_gap_columns(refined)
  let score = sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "mafft" }
}

///|
/// MAFFT-lite alignment with custom parameters.
pub fn mafft_align_with_params(
  sequences : Array[String],
  names : Array[String],
  k : Int,
  gap_open : Double,
  num_refine_rounds : Int,
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "mafft" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "mafft" }
  }
  let matrix = blosum62()
  let dist = mafft_distance_matrix(sequences, k)
  let tree = mafft_upgma_tree(dist)
  let (alignment_raw, order) = progressive_align(tree, sequences, matrix, gap_open)
  let alignment = reorder_alignment(alignment_raw, order)
  let mut refined = alignment
  for _i = 0; _i < num_refine_rounds; _i = _i + 1 {
    refined = refine_once(refined, matrix, gap_open)
  }
  let cleaned = remove_all_gap_columns(refined)
  let score = sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "mafft" }
}

///|
/// Convert MSA result to FASTA format.
pub fn msa_to_fasta(result : MsaResult) -> String {
  let buf : StringBuilder = StringBuilder()
  for i = 0; i < result.aligned_sequences.length(); i = i + 1 {
    let name = if i < result.names.length() { result.names[i] } else { "seq" + i.to_string() }
    buf.write_string(">" + name + "\n")
    buf.write_string(result.aligned_sequences[i] + "\n")
  }
  buf.to_string()
}

///|
/// Convert MSA result to Clustal format (simplified).
pub fn msa_to_clustal(result : MsaResult) -> String {
  let buf : StringBuilder = StringBuilder()
  let n = result.aligned_sequences.length()
  if n == 0 {
    return ""
  }
  let m = result.aligned_sequences[0].length()
  let block_size = 60
  let mut pos = 0
  while pos < m {
    let end = if pos + block_size < m { pos + block_size } else { m }
    for i = 0; i < n; i = i + 1 {
      let name = if i < result.names.length() { result.names[i] } else { "seq" + i.to_string() }
      let seq_part = result.aligned_sequences[i][pos:end].to_string()
      buf.write_string(name + "\t" + seq_part + "\n")
    }
    buf.write_string("\n")
    pos = end
  }
  buf.to_string()
}
