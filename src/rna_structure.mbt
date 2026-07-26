///|
/// RNA secondary structure prediction using the Nussinov algorithm.
///
/// References:
/// - Nussinov M (1980). "A simple loop minimization algorithm for the
///   prediction of secondary structure in ribonucleic acids."
///   J theor Biol. 86(3):373-378.
/// - Nussinov M, Jacobson AB (1980). "Fast algorithm for the prediction
///   of RNA secondary structures." Proc Natl Acad Sci USA. 77(12):6309-6313.
/// - Zuker M (2003). "Mfold web server for nucleic acid folding and
///   hybridization prediction." Nucleic Acids Res. 31(13):3406-3415.

///|
/// Result of RNA secondary structure prediction.
/// Contains the input sequence, predicted structure in dot-bracket notation,
/// minimum free energy, base pairs, and classified loop regions.
pub struct RNAStructureResult {
  sequence : String
  structure : String
  mfe : Double
  paired_positions : Array[(Int, Int)]
  loop_regions : Array[LoopRegion]
}

///|
/// A loop region in the RNA secondary structure.
/// Classified as hairpin, internal, bulge, or multi-loop.
pub struct LoopRegion {
  start : Int
  end : Int
  loop_type : String
} derive(Eq)

///|
/// Configurable parameters for RNA structure prediction.
///
/// - loop_min_length: minimum number of unpaired bases in a loop (default 4)
/// - stacking_energy: bonus for consecutive base pairs (default -2.0 kcal/mol)
/// - loop_penalty: base loop penalty per unpaired base (default 0.5 kcal/mol)
/// - mismatch_penalty: penalty for non-canonical base pairs (default 1.0)
pub struct RNAStructurePredictor {
  loop_min_length : Int
  stacking_energy : Double
  loop_penalty : Double
  mismatch_penalty : Double
}

///|
/// Create an RNA structure predictor with default parameters.
pub fn rna_predictor_new() -> RNAStructurePredictor {
  RNAStructurePredictor::{
    loop_min_length: 4,
    stacking_energy: -2.0,
    loop_penalty: 0.5,
    mismatch_penalty: 1.0,
  }
}

///|
/// Create an RNA structure predictor with custom parameters.
///
/// - loop_min: minimum loop length (number of unpaired bases)
/// - stacking: stacking energy bonus in kcal/mol (negative = favorable)
/// - loop_penalty: penalty per unpaired base in a loop
/// - mismatch_penalty: penalty for non-canonical base pairs
pub fn rna_predictor_with_params(
  loop_min : Int,
  stacking : Double,
  loop_penalty : Double,
  mismatch_penalty : Double,
) -> RNAStructurePredictor {
  RNAStructurePredictor::{
    loop_min_length: loop_min,
    stacking_energy: stacking,
    loop_penalty,
    mismatch_penalty,
  }
}

///|
/// Check if two RNA bases can form a base pair.
/// Supports A-U, U-A, G-C, C-G, G-U, U-G (wobble) pairs.
/// Both uppercase and lowercase are accepted.
pub fn is_complementary(base1 : UInt16, base2 : UInt16) -> Bool {
  let b1 = rna_normalize_base(base1)
  let b2 = rna_normalize_base(base2)
  if b1 == 65 && b2 == 85 { return true }
  if b1 == 85 && b2 == 65 { return true }
  if b1 == 71 && b2 == 67 { return true }
  if b1 == 67 && b2 == 71 { return true }
  if b1 == 71 && b2 == 85 { return true }
  if b1 == 85 && b2 == 71 { return true }
  false
}

///|
/// Calculate the loop energy penalty for a given loop length.
/// Returns penalty = loop_penalty * length.
pub fn compute_loop_energy(length : Int) -> Double {
  if length <= 0 { return 0.0 }
  0.5 * length.to_double()
}

///|
/// Run the Nussinov algorithm to compute the DP matrix.
/// Returns a 2D array where matrix[i][j] holds the minimum free energy
/// for the subsequence from i to j (inclusive).
///
/// Time complexity: O(n^3), Space complexity: O(n^2).
pub fn nussinov_matrix(
  sequence : Array[UInt16],
  predictor : RNAStructurePredictor,
) -> Array[Array[Double]] {
  let n = sequence.length()
  if n == 0 { return [] }

  let dp : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n; j = j + 1 {
      row.push(0.0)
    }
    dp.push(row)
  }

  for len = 1; len < n; len = len + 1 {
    for i = 0; i < n - len; i = i + 1 {
      let j = i + len
      let min_energy = rna_compute_energy_for_interval(
        sequence, i, j, predictor, dp,
      )
      dp[i][j] = min_energy
    }
  }

  dp
}

///|
/// Traceback through the DP matrix to find the optimal set of base pairs.
/// Returns an array of (i, j) tuples representing paired positions.
pub fn nussinov_traceback(
  matrix : Array[Array[Double]],
  sequence : Array[UInt16],
  predictor : RNAStructurePredictor,
) -> Array[(Int, Int)] {
  let n = sequence.length()
  let base_pairs : Array[(Int, Int)] = Array::new()
  if n == 0 { return base_pairs }

  rna_trace_from_matrix(matrix, sequence, predictor, 0, n - 1, base_pairs)

  base_pairs
}

///|
/// Predict the secondary structure of an RNA sequence using the Nussinov algorithm.
///
/// - predictor: configuration parameters
/// - sequence: RNA sequence string (A, U, G, C characters)
///
/// Returns an RNAStructureResult with the predicted structure.
pub fn rna_predict_secondary_structure(
  predictor : RNAStructurePredictor,
  sequence : String,
) -> RNAStructureResult {
  let n = sequence.length()
  if n == 0 {
    return RNAStructureResult::{
      sequence,
      structure: "",
      mfe: 0.0,
      paired_positions: [],
      loop_regions: [],
    }
  }

  let seq_arr : Array[UInt16] = Array::new()
  for i = 0; i < n; i = i + 1 {
    seq_arr.push(sequence.unsafe_get(i))
  }

  if n <= 2 {
    let structure_str = rna_build_dot_bracket(n, [])
    return RNAStructureResult::{
      sequence,
      structure: structure_str,
      mfe: 0.0,
      paired_positions: [],
      loop_regions: [],
    }
  }

  let matrix = nussinov_matrix(seq_arr, predictor)
  let base_pairs = nussinov_traceback(matrix, seq_arr, predictor)
  let structure_str = mfe_to_dot_bracket(n, base_pairs)
  let loop_regions = identify_loop_regions(seq_arr, base_pairs, predictor)
  let mfe = rna_compute_mfe(seq_arr, base_pairs, predictor, loop_regions)

  RNAStructureResult::{
    sequence,
    structure: structure_str,
    mfe,
    paired_positions: base_pairs,
    loop_regions,
  }
}

///|
/// Alias for rna_predict_secondary_structure.
pub fn rna_predict_with_mfe(
  predictor : RNAStructurePredictor,
  sequence : String,
) -> RNAStructureResult {
  rna_predict_secondary_structure(predictor, sequence)
}

///|
/// Get the dot-bracket notation from an RNA structure result.
pub fn structure_to_dot_bracket(result : RNAStructureResult) -> String {
  result.structure
}

///|
/// Get the list of base pairs from an RNA structure result.
pub fn get_base_pairs(result : RNAStructureResult) -> Array[(Int, Int)] {
  result.paired_positions
}

///|
/// Get the minimum free energy from an RNA structure result.
pub fn get_mfe(result : RNAStructureResult) -> Double {
  result.mfe
}

///|
/// Classify loop regions in the RNA secondary structure.
/// Identifies hairpin, internal, bulge, and multi-loop regions.
pub fn identify_loop_regions(
  sequence : Array[UInt16],
  base_pairs : Array[(Int, Int)],
  predictor : RNAStructurePredictor,
) -> Array[LoopRegion] {
  let loop_regions : Array[LoopRegion] = Array::new()
  let n = sequence.length()
  if n == 0 { return loop_regions }

  let paired = rna_build_paired_array(n, base_pairs)
  rna_find_loops_recursive(
    sequence, base_pairs, 0, n - 1, paired, loop_regions,
  )

  loop_regions
}

///|
/// Convert base pairs to dot-bracket notation.
///
/// - sequence_length: length of the RNA sequence
/// - base_pairs: array of (i, j) base pair positions
///
/// Returns a string where '(' and ')' indicate paired bases
/// and '.' indicates unpaired bases.
pub fn mfe_to_dot_bracket(
  sequence_length : Int,
  base_pairs : Array[(Int, Int)],
) -> String {
  rna_build_dot_bracket(sequence_length, base_pairs)
}

// --- Internal helper functions ---

///|
/// Normalize a base to uppercase UInt16 code.
fn rna_normalize_base(base : UInt16) -> UInt16 {
  if base >= 97 && base <= 122 { base - 32 } else { base }
}

///|
/// Get the base pair energy contribution.
/// Watson-Crick (A-U, G-C): -1.0 (favorable)
/// Wobble (G-U): -0.5 (less favorable)
/// Non-complementary: mismatch_penalty (unfavorable)
fn rna_base_pair_energy(
  b1 : UInt16,
  b2 : UInt16,
  predictor : RNAStructurePredictor,
) -> Double {
  let nb1 = rna_normalize_base(b1)
  let nb2 = rna_normalize_base(b2)

  let is_watson_crick = (nb1 == 65 && nb2 == 85) ||
    (nb1 == 85 && nb2 == 65) ||
    (nb1 == 71 && nb2 == 67) ||
    (nb1 == 67 && nb2 == 71)

  if is_watson_crick { return -1.0 }

  let is_wobble = (nb1 == 71 && nb2 == 85) || (nb1 == 85 && nb2 == 71)
  if is_wobble { return -0.5 }

  predictor.mismatch_penalty
}

///|
/// Compute the minimum energy for the interval [i, j].
/// Used during DP matrix filling.
fn rna_compute_energy_for_interval(
  sequence : Array[UInt16],
  i : Int,
  j : Int,
  predictor : RNAStructurePredictor,
  dp : Array[Array[Double]],
) -> Double {
  let loop_min = predictor.loop_min_length

  let unpaired_i = dp[i + 1][j]
  let unpaired_j = dp[i][j - 1]

  let mut min_energy = if unpaired_i < unpaired_j { unpaired_i } else { unpaired_j }

  let loop_len = j - i - 1
  if loop_len >= loop_min {
    let bp_energy = rna_base_pair_energy(sequence[i], sequence[j], predictor)
    let loop_penalty = 0.6 * loop_len.to_double()
    let paired_energy = dp[i + 1][j - 1] + bp_energy + loop_penalty
    if paired_energy < min_energy {
      min_energy = paired_energy
    }
  }

  let mut k = i
  while k < j {
    let bifurcation_energy = dp[i][k] + dp[k + 1][j]
    if bifurcation_energy < min_energy {
      min_energy = bifurcation_energy
    }
    k = k + 1
  }

  min_energy
}

///|
/// Absolute value helper for Double.
fn rna_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}

///|
/// Traceback through the DP matrix directly to find base pairs.
fn rna_trace_from_matrix(
  matrix : Array[Array[Double]],
  sequence : Array[UInt16],
  predictor : RNAStructurePredictor,
  i : Int,
  j : Int,
  base_pairs : Array[(Int, Int)],
) -> Unit {
  if i >= j { return }

  let n = sequence.length()
  let loop_min = predictor.loop_min_length
  let current = matrix[i][j]
  let eps = 1.0e-9

  if i + 1 < n && rna_abs(matrix[i + 1][j] - current) < eps {
    rna_trace_from_matrix(matrix, sequence, predictor, i + 1, j, base_pairs)
    return
  }

  if j - 1 >= 0 && rna_abs(matrix[i][j - 1] - current) < eps {
    rna_trace_from_matrix(matrix, sequence, predictor, i, j - 1, base_pairs)
    return
  }

  let loop_len = j - i - 1
  if loop_len >= loop_min && i + 1 < n && j - 1 >= 0 {
    let bp_energy = rna_base_pair_energy(sequence[i], sequence[j], predictor)
    let loop_penalty = 0.6 * loop_len.to_double()
    let paired_energy = matrix[i + 1][j - 1] + bp_energy + loop_penalty
    if rna_abs(paired_energy - current) < eps {
      base_pairs.push((i, j))
      rna_trace_from_matrix(
        matrix, sequence, predictor, i + 1, j - 1, base_pairs,
      )
      return
    }
  }

  let mut k = i
  while k < j {
    let bifurcation_energy = matrix[i][k] + matrix[k + 1][j]
    if rna_abs(bifurcation_energy - current) < eps {
      rna_trace_from_matrix(matrix, sequence, predictor, i, k, base_pairs)
      rna_trace_from_matrix(
        matrix, sequence, predictor, k + 1, j, base_pairs,
      )
      return
    }
    k = k + 1
  }
}

///|
/// Build a dot-bracket string from base pairs.
fn rna_build_dot_bracket(n : Int, base_pairs : Array[(Int, Int)]) -> String {
  let buf : FixedArray[UInt16] = FixedArray::make(n, '.'.to_int().to_uint16())
  for idx = 0; idx < base_pairs.length(); idx = idx + 1 {
    let (i, j) = base_pairs[idx]
    buf[i] = '('.to_int().to_uint16()
    buf[j] = ')'.to_int().to_uint16()
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Build an array indicating which positions are paired.
fn rna_build_paired_array(n : Int, base_pairs : Array[(Int, Int)]) -> Array[Int] {
  let paired : Array[Int] = Array::make(n, -1)
  for idx = 0; idx < base_pairs.length(); idx = idx + 1 {
    let (i, j) = base_pairs[idx]
    paired[i] = j
    paired[j] = i
  }
  paired
}

///|
/// Recursively find and classify loop regions.
fn rna_find_loops_recursive(
  sequence : Array[UInt16],
  base_pairs : Array[(Int, Int)],
  i : Int,
  j : Int,
  paired : Array[Int],
  loop_regions : Array[LoopRegion],
) -> Unit {
  if i >= j { return }

  let n = sequence.length()
  let has_pair = if i < n && paired[i] == j { true } else { false }

  if has_pair {
    let loop_len = j - i - 1
    if loop_len > 0 {
      let loop_type = rna_classify_loop(sequence, i, j, paired)
      loop_regions.push(LoopRegion::{
        start: i + 1,
        end: j - 1,
        loop_type,
      })

      let inner_i = i + 1
      let inner_j = j - 1
      if inner_i < inner_j {
        let sub_paired = rna_build_sub_paired(inner_i, inner_j, paired)
        rna_find_loops_recursive(
          sequence, base_pairs, inner_i, inner_j, sub_paired, loop_regions,
        )
      }
    }
  } else {
    let mut k = i
    while k < j {
      let k_pair = if k < n { paired[k] } else { -1 }
      if k_pair > 0 && k_pair <= j {
        let loop_type = rna_classify_loop(sequence, k, k_pair, paired)
        loop_regions.push(LoopRegion::{
          start: k + 1,
          end: k_pair - 1,
          loop_type,
        })
        if k + 1 < k_pair - 1 {
          let sub_paired = rna_build_sub_paired(k + 1, k_pair - 1, paired)
          rna_find_loops_recursive(
            sequence, base_pairs, k + 1, k_pair - 1, sub_paired, loop_regions,
          )
        }
        k = k_pair + 1
      } else {
        k = k + 1
      }
    }
  }
}

///|
/// Classify a loop as hairpin, internal, bulge, or multi-loop.
fn rna_classify_loop(
  sequence : Array[UInt16],
  i : Int,
  j : Int,
  paired : Array[Int],
) -> String {
  let n = sequence.length()

  let has_inner_pairs = rna_has_inner_pairs(i, j, paired)
  let has_left_pair = i > 0 && paired[i - 1] >= 0
  let has_right_pair = j < n - 1 && paired[j + 1] >= 0

  if has_inner_pairs {
    if has_left_pair && has_right_pair {
      return "multi"
    }
    if has_left_pair || has_right_pair {
      return "multi"
    }
    return "hairpin"
  }

  if has_left_pair && has_right_pair {
    return "internal"
  }
  if has_left_pair || has_right_pair {
    return "bulge"
  }

  "hairpin"
}

///|
/// Check if there are base pairs within the loop [i+1, j-1].
fn rna_has_inner_pairs(i : Int, j : Int, paired : Array[Int]) -> Bool {
  let mut k = i + 1
  while k < j {
    if paired[k] > 0 && paired[k] < j && paired[k] > i {
      return true
    }
    k = k + 1
  }
  false
}

///|
/// Build a sub-array of paired information for a given range.
fn rna_build_sub_paired(
  start : Int,
  end : Int,
  paired : Array[Int],
) -> Array[Int] {
  let n = paired.length()
  let sub : Array[Int] = Array::make(n, -1)
  let mut k = start
  while k <= end {
    let p = paired[k]
    if p >= start && p <= end {
      sub[k] = p
      sub[p] = k
    }
    k = k + 1
  }
  sub
}

///|
/// Compute the final minimum free energy.
fn rna_compute_mfe(
  sequence : Array[UInt16],
  base_pairs : Array[(Int, Int)],
  predictor : RNAStructurePredictor,
  loop_regions : Array[LoopRegion],
) -> Double {
  let n = sequence.length()
  if n <= 2 { return 0.0 }

  let mut mfe = 0.0

  for idx = 0; idx < base_pairs.length(); idx = idx + 1 {
    let (i, j) = base_pairs[idx]
    let bp_e = rna_base_pair_energy(sequence[i], sequence[j], predictor)
    mfe = mfe + bp_e
  }

  for idx = 0; idx < loop_regions.length(); idx = idx + 1 {
    let region = loop_regions[idx]
    let loop_len = region.end - region.start + 1
    if loop_len <= 0 { continue }

    let penalty = match region.loop_type {
      "hairpin" => (predictor.loop_penalty + 0.1) * loop_len.to_double()
      "internal" => (predictor.loop_penalty + 0.3) * loop_len.to_double()
      "bulge" => (predictor.loop_penalty + 0.1) * loop_len.to_double()
      "multi" => (predictor.loop_penalty + 0.3) * loop_len.to_double()
      _ => predictor.loop_penalty * loop_len.to_double()
    }
    mfe = mfe + penalty
  }

  let stacking = rna_count_stacking(base_pairs, predictor)
  mfe = mfe + stacking

  mfe
}

///|
/// Count stacking interactions and return the total stacking energy.
fn rna_count_stacking(
  base_pairs : Array[(Int, Int)],
  predictor : RNAStructurePredictor,
) -> Double {
  let mut stacking_count = 0
  let bp_count = base_pairs.length()

  for idx = 0; idx < bp_count; idx = idx + 1 {
    let (i, j) = base_pairs[idx]
    let has_next = i + 1 < j - 1
    if !has_next { continue }

    for idx2 = 0; idx2 < bp_count; idx2 = idx2 + 1 {
      if idx2 == idx { continue }
      let (i2, j2) = base_pairs[idx2]
      if i2 == i + 1 && j2 == j - 1 {
        stacking_count = stacking_count + 1
        break
      }
    }
  }

  stacking_count.to_double() * predictor.stacking_energy
}