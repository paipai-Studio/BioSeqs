/// biotite_align_xdrop.mbt — Local gapped alignment with X-drop termination
/// Inspired by biotite.sequence.align.align_local_gapped

///|
/// X-drop alignment result.
pub struct BiotiteXdropResult {
  score : Double
  start1 : Int32 // start position in sequence 1
  end1 : Int32 // end position in sequence 1
  start2 : Int32 // start position in sequence 2
  end2 : Int32 // end position in sequence 2
  alignment : Array[(Int32, Int32)] // aligned position pairs (-1 = gap)
}

///|
/// Perform local gapped alignment with X-drop termination.
///
/// The X-drop algorithm extends an alignment from a seed until the score
/// drops by more than X from the best score seen so far. This is the
/// algorithm used by BLAST's gapped extension step.
///
/// Parameters:
/// - seq1, seq2: sequences as integer code arrays
/// - seed1, seed2: starting position of the seed in seq1 and seq2
/// - seed_len: length of the seed
/// - sub_matrix: substitution matrix [code1][code2] -> score
/// - gap_open: gap opening penalty (positive)
/// - gap_extend: gap extension penalty (positive)
/// - x_drop: X-drop threshold (stop when score drops X below best)
pub fn biotite_align_xdrop(
  seq1 : Array[Int32],
  seq2 : Array[Int32],
  seed1 : Int32,
  seed2 : Int32,
  seed_len : Int32,
  sub_matrix : Array[Array[Double]],
  _gap_open~ : Double = 5.0,
  _gap_extend~ : Double = 1.0,
  x_drop~ : Double = 20.0,
) -> BiotiteXdropResult {
  // Start from the seed and extend in both directions
  // Forward extension
  let (fwd_end1, fwd_end2, fwd_score, fwd_align) = biotite_xdrop_extend(
    seq1,
    seq2,
    seed1 + seed_len,
    seed2 + seed_len,
    1, // forward direction
    sub_matrix,
    x_drop,
  )
  // Backward extension
  let (bwd_end1, bwd_end2, bwd_score, bwd_align) = biotite_xdrop_extend(
    seq1,
    seq2,
    seed1 - 1,
    seed2 - 1,
    -1, // backward direction
    sub_matrix,
    x_drop,
  )
  // Combine results
  let start1 = bwd_end1 + 1
  let start2 = bwd_end2 + 1
  let end1 = fwd_end1
  let end2 = fwd_end2
  // Build alignment: backward (reversed) + seed + forward
  let alignment : Array[(Int32, Int32)] = Array(capacity=0)
  for i = bwd_align.length() - 1; i >= 0; i = i - 1 {
    alignment.push(bwd_align[i])
  }
  for i = 0; i < seed_len; i = i + 1 {
    alignment.push((seed1 + i, seed2 + i))
  }
  for i = 0; i < fwd_align.length(); i = i + 1 {
    alignment.push(fwd_align[i])
  }
  { score: fwd_score + bwd_score, start1, end1, start2, end2, alignment, }
}

///|
/// Extend alignment in one direction with X-drop.
/// Returns (end_pos1, end_pos2, score, alignment_pairs)
fn biotite_xdrop_extend(
  seq1 : Array[Int32],
  seq2 : Array[Int32],
  start1 : Int32,
  start2 : Int32,
  direction : Int32,
  sub_matrix : Array[Array[Double]],
  x_drop : Double,
) -> (Int32, Int32, Double, Array[(Int32, Int32)]) {
  let n1 = seq1.length()
  let n2 = seq2.length()
  let alignment : Array[(Int32, Int32)] = Array(capacity=0)
  let mut i = start1
  let mut j = start2
  let mut best_score = 0.0
  let mut current_score = 0.0
  let mut last_i = start1 - direction
  let mut last_j = start2 - direction
  while true {
    // Check bounds
    let in_bounds = if direction > 0 {
      i < n1 && j < n2
    } else {
      i >= 0 && j >= 0
    }
    if !in_bounds {
      break
    }
    // Diagonal step (match/mismatch)
    let diag_score = if i >= 0 && i < n1 && j >= 0 && j < n2 {
      if seq1[i] < sub_matrix.length() && seq2[j] < sub_matrix[0].length() {
        sub_matrix[seq1[i]][seq2[j]]
      } else {
        0.0
      }
    } else {
      -1000.0
    }
    if diag_score > -100.0 {
      current_score = current_score + diag_score
      alignment.push((i, j))
      last_i = i
      last_j = j
      i = i + direction
      j = j + direction
    } else {
      break
    }
    // Update best score
    if current_score > best_score {
      best_score = current_score
    }
    // X-drop check
    if best_score - current_score > x_drop {
      break
    }
  }
  (last_i, last_j, best_score, alignment)
}
