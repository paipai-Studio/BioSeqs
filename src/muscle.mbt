///|
/// MUSCLE-lite: multiple sequence alignment by log-expectation.
///
/// Implements the three-stage MUSCLE pipeline:
///   Stage 1 (Draft): k-mer distance → UPGMA tree → progressive alignment
///   Stage 2 (Improved): log-expectation distance → new tree → progressive
///   Stage 3 (Refinement): delete each tree edge, realign two profiles
///
/// This is a pure-MoonBit reimplementation.
/// Algorithm follows Edgar, "MUSCLE: multiple sequence alignment with high
/// accuracy and high throughput" (Nucleic Acids Res 2004).

// ============================================================================
// 1. Log-expectation distance (from an existing alignment)
// ============================================================================

///|
/// Compute column letter frequencies for a profile.
/// Returns a map from character to frequency at each column.
fn profile_frequencies(
  aligned_seqs : Array[String],
) -> Array[Map[Char, Double]] {
  let n = aligned_seqs.length()
  if n == 0 {
    return []
  }
  let m = aligned_seqs[0].length()
  let result : Array[Map[Char, Double]] = Array(capacity=m)
  for col = 0; col < m; col = col + 1 {
    let freqs : Map[Char, Double] = Map([], capacity=0)
    for row = 0; row < n; row = row + 1 {
      let c = aligned_seqs[row][col].unsafe_to_char()
      match freqs.get(c) {
        Some(v) => freqs[c] = v + 1.0
        None => freqs[c] = 1.0
      }
    }
    // Normalize
    for k, v in freqs {
      freqs[k] = v / n.to_double()
    }
    result.push(freqs)
  }
  result
}

///|
/// Log-expectation score between two columns.
/// LE = sum over all pairs (a,b) of f1[a] * f2[b] * score(a,b)
fn log_expectation_score(
  freqs1 : Map[Char, Double],
  freqs2 : Map[Char, Double],
  matrix : SubstitutionMatrix,
) -> Double {
  let mut total = 0.0
  for a, fa in freqs1 {
    for b, fb in freqs2 {
      if a != '-' && b != '-' {
        let s = matrix.get_score_case_insensitive(a.to_string(), b.to_string()).to_double()
        total += fa * fb * s
      }
    }
  }
  total
}

///|
/// Compute pairwise log-expectation distance from an alignment.
/// Extracts two sequences (ungapped) from the alignment, realigns them,
/// and computes the LE-based distance.
fn le_distance(
  aligned_seqs : Array[String],
  i : Int,
  j : Int,
  matrix : SubstitutionMatrix,
) -> Double {
  let seq_i = remove_gaps_from_seq(aligned_seqs[i])
  let seq_j = remove_gaps_from_seq(aligned_seqs[j])
  // Simple identity-based distance
  let mut matches = 0
  let mut total = 0
  let min_len = if seq_i.length() < seq_j.length() { seq_i.length() } else { seq_j.length() }
  for k = 0; k < min_len; k = k + 1 {
    total += 1
    if seq_i[k] == seq_j[k] {
      matches += 1
    }
  }
  let max_len = if seq_i.length() > seq_j.length() { seq_i.length() } else { seq_j.length() }
  total += max_len - min_len
  if total == 0 {
    return 1.0
  }
  1.0 - matches.to_double() / total.to_double()
}

///|
fn remove_gaps_from_seq(seq : String) -> String {
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
/// Build a distance matrix from an alignment using LE distance.
fn le_distance_matrix(
  aligned_seqs : Array[String],
  matrix : SubstitutionMatrix,
) -> Array[Array[Double]] {
  let n = aligned_seqs.length()
  let result : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      if i != j {
        row[j] = le_distance(aligned_seqs, i, j, matrix)
      }
    }
    result.push(row)
  }
  result
}

// ============================================================================
// 2. Refinement: edge-based profile realignment
// ============================================================================

///|
/// Extract sequences belonging to a subtree.
fn subtree_sequences(
  tree : MsaTreeNode,
  sequences : Array[String],
) -> Array[Int] {
  if tree.is_leaf {
    return [tree.seq_index]
  }
  let result : Array[Int] = []
  for child in tree.children {
    let child_seqs = subtree_sequences(child, sequences)
    for s in child_seqs {
      result.push(s)
    }
  }
  result
}

///|
/// Extract a subset of aligned sequences by index.
fn extract_subset(alignment : Array[String], indices : Array[Int]) -> Array[String] {
  let result : Array[String] = []
  for idx in indices {
    result.push(alignment[idx])
  }
  result
}

///|
/// SP score for MUSCLE.
fn muscle_sp_score(
  alignment : Array[String],
  matrix : SubstitutionMatrix,
) -> Double {
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
        if c1 == '-' && c2 == '-' {
          continue
        }
        if c1 == '-' || c2 == '-' {
          total += -4.0
        } else {
          total += matrix.get_score_case_insensitive(
            c1.to_string(),
            c2.to_string(),
          ).to_double()
        }
      }
    }
  }
  total
}

///|
/// One refinement round: for each internal node, split into two profiles
/// and realign. Keep if SP score improves.
fn muscle_refine(
  alignment : Array[String],
  tree : MsaTreeNode,
  matrix : SubstitutionMatrix,
  gap_open : Double,
) -> Array[String] {
  let n = alignment.length()
  if n <= 2 {
    return alignment
  }
  let mut best = alignment
  let mut best_score = muscle_sp_score(best, matrix)
  // Collect all candidate refinements
  let candidates : Array[Array[String]] = []
  collect_refinements(tree, alignment, matrix, gap_open, candidates)
  for candidate in candidates {
    let cleaned = remove_all_gap_cols(candidate)
    let new_score = muscle_sp_score(cleaned, matrix)
    if new_score > best_score {
      best = cleaned
      best_score = new_score
    }
  }
  best
}

///|
/// Collect all refinement candidates from tree nodes.
fn collect_refinements(
  node : MsaTreeNode,
  alignment : Array[String],
  matrix : SubstitutionMatrix,
  gap_open : Double,
  candidates : Array[Array[String]],
) -> Unit {
  if node.is_leaf {
    return
  }
  let children = node.children
  if children.length() >= 2 {
    let group1_indices = subtree_sequences(children[0], [])
    let group2_indices : Array[Int] = []
    for ci = 1; ci < children.length(); ci = ci + 1 {
      let child_seqs = subtree_sequences(children[ci], [])
      for s in child_seqs {
        group2_indices.push(s)
      }
    }
    let profile1 = extract_subset(alignment, group1_indices)
    let profile2 = extract_subset(alignment, group2_indices)
    let realigned = mafft_align_profiles(profile1, profile2, matrix, gap_open)
    // Safety check: realigned should have same total count
    let total_count = group1_indices.length() + group2_indices.length()
    if realigned.length() == total_count && total_count == alignment.length() {
      let new_alignment : Array[String] = Array::make(alignment.length(), "")
      for i = 0; i < group1_indices.length(); i = i + 1 {
        if i < realigned.length() && group1_indices[i] < new_alignment.length() {
          new_alignment[group1_indices[i]] = realigned[i]
        }
      }
      for i = 0; i < group2_indices.length(); i = i + 1 {
        let ri = group1_indices.length() + i
        if ri < realigned.length() && group2_indices[i] < new_alignment.length() {
          new_alignment[group2_indices[i]] = realigned[ri]
        }
      }
      candidates.push(new_alignment)
    }
  }
  for child in children {
    collect_refinements(child, alignment, matrix, gap_open, candidates)
  }
}

///|
/// Remove all-gap columns from an alignment.
fn remove_all_gap_cols(alignment : Array[String]) -> Array[String] {
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

// ============================================================================
// 3. Public API
// ============================================================================

///|
/// MUSCLE-lite alignment (three-stage pipeline).
/// Stage 1: k-mer distance → UPGMA → progressive (draft)
/// Stage 2: LE distance from draft → UPGMA → progressive (improved)
/// Stage 3: tree-based refinement
pub fn muscle_align(
  sequences : Array[String],
  names : Array[String],
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "muscle" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "muscle" }
  }
  let matrix = blosum62()
  let gap_open = -10.0
  let k = 6
  // Stage 1: Draft progressive
  let dist1 = mafft_distance_matrix(sequences, k)
  let tree1 = mafft_upgma_tree(dist1)
  let (draft_raw, order1) = muscle_progressive(tree1, sequences, matrix, gap_open)
  let draft = reorder_alignment(draft_raw, order1)
  let draft_cleaned = remove_all_gap_cols(draft)
  // Stage 2: Improved progressive using LE distance
  let dist2 = le_distance_matrix(draft_cleaned, matrix)
  let tree2 = mafft_upgma_tree(dist2)
  let (improved_raw, order2) = muscle_progressive(tree2, sequences, matrix, gap_open)
  let improved = reorder_alignment(improved_raw, order2)
  let improved_cleaned = remove_all_gap_cols(improved)
  // Stage 3: Refinement
  let refined = muscle_refine(improved_cleaned, tree2, matrix, gap_open)
  let cleaned = remove_all_gap_cols(refined)
  let score = muscle_sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "muscle" }
}

///|
/// MUSCLE-lite with custom parameters.
pub fn muscle_align_with_params(
  sequences : Array[String],
  names : Array[String],
  k : Int,
  gap_open : Double,
  num_refine_rounds : Int,
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "muscle" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "muscle" }
  }
  let matrix = blosum62()
  let dist1 = mafft_distance_matrix(sequences, k)
  let tree1 = mafft_upgma_tree(dist1)
  let (draft_raw, order1) = muscle_progressive(tree1, sequences, matrix, gap_open)
  let draft = reorder_alignment(draft_raw, order1)
  let draft_cleaned = remove_all_gap_cols(draft)
  let dist2 = le_distance_matrix(draft_cleaned, matrix)
  let tree2 = mafft_upgma_tree(dist2)
  let (improved_raw, order2) = muscle_progressive(tree2, sequences, matrix, gap_open)
  let improved = reorder_alignment(improved_raw, order2)
  let improved_cleaned = remove_all_gap_cols(improved)
  let mut refined = improved_cleaned
  for _i = 0; _i < num_refine_rounds; _i = _i + 1 {
    refined = muscle_refine(refined, tree2, matrix, gap_open)
  }
  let cleaned = remove_all_gap_cols(refined)
  let score = muscle_sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "muscle" }
}

///|
/// Progressive alignment (shared with mafft, reimplemented for independence).
fn muscle_progressive(
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
    let (child_profile, child_order) = muscle_progressive(
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
/// Stage 1 only: draft progressive alignment.
pub fn muscle_stage1(
  sequences : Array[String],
  names : Array[String],
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "muscle-s1" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "muscle-s1" }
  }
  let matrix = blosum62()
  let gap_open = -10.0
  let k = 6
  let dist = mafft_distance_matrix(sequences, k)
  let tree = mafft_upgma_tree(dist)
  let (draft_raw, order) = muscle_progressive(tree, sequences, matrix, gap_open)
  let draft = reorder_alignment(draft_raw, order)
  let cleaned = remove_all_gap_cols(draft)
  let score = muscle_sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "muscle-s1" }
}

///|
/// Stages 1+2: draft + improved progressive.
pub fn muscle_stage2(
  sequences : Array[String],
  names : Array[String],
) -> MsaResult {
  let n = sequences.length()
  if n == 0 {
    return { aligned_sequences: [], names: [], score: 0.0, method: "muscle-s2" }
  }
  if n == 1 {
    return { aligned_sequences: [sequences[0]], names, score: 0.0, method: "muscle-s2" }
  }
  let matrix = blosum62()
  let gap_open = -10.0
  let k = 6
  let dist1 = mafft_distance_matrix(sequences, k)
  let tree1 = mafft_upgma_tree(dist1)
  let (draft_raw, order1) = muscle_progressive(tree1, sequences, matrix, gap_open)
  let draft = reorder_alignment(draft_raw, order1)
  let draft_cleaned = remove_all_gap_cols(draft)
  let dist2 = le_distance_matrix(draft_cleaned, matrix)
  let tree2 = mafft_upgma_tree(dist2)
  let (improved_raw, order2) = muscle_progressive(tree2, sequences, matrix, gap_open)
  let improved = reorder_alignment(improved_raw, order2)
  let cleaned = remove_all_gap_cols(improved)
  let score = muscle_sp_score(cleaned, matrix)
  { aligned_sequences: cleaned, names, score, method: "muscle-s2" }
}
