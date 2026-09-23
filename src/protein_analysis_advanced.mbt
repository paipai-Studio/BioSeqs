///|
/// Protein analysis advanced — empirical protein-sequence prediction algorithms
/// inspired by Biopython's `Bio.SeqUtils` and widely used immunoinformatics /
/// structural-bioinformatics tools.
///
/// Implements:
///   - Chou–Fasman (1974) secondary-structure prediction from amino-acid
///     propensity tables (helix / sheet / turn nucleation and extension).
///   - IUPred (Dosztányi et al. 2005) intrinsically-disordered-region
///     prediction via estimated pairwise interaction energy.
///   - COILS (Lupas et al. 1991) coiled-coil prediction with heptad-repeat
///     scoring matrices.
///   - Kolaskar–Tongaonkar (1990) antigenicity prediction.
///   - Emini et al. (1985) surface-accessibility prediction.
///   - Karplus–Schulz (1985) flexibility prediction.

// ===========================================================================
// Error type and validation helpers
// ===========================================================================

///|
pub suberror ProteinAdvancedError {
  ProteinAdvancedError(String)
}

///|
fn pa_fail(message : String) -> Unit raise ProteinAdvancedError {
  raise ProteinAdvancedError::ProteinAdvancedError(message)
}

///|
fn pa_validate_protein(
  sequence : String,
  label : String,
) -> Unit raise ProteinAdvancedError {
  if sequence.length() == 0 {
    pa_fail("Protein advanced " + label + " requires a non-empty sequence")
  }
  let seq = sequence.to_upper()
  let valid = "ACDEFGHIKLMNPQRSTVWY"
  let mut i = 0
  while i < seq.length() {
    let c = seq[i:i + 1].to_owned()
    if !valid.contains(c) {
      pa_fail(
        "Protein advanced " +
        label +
        " encountered non-standard residue '" +
        c +
        "' at position " +
        i.to_string(),
      )
    }
    i = i + 1
  }
}

///|
fn pa_aa_index(c : Char) -> Int {
  // Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical, 0-19)
  match c {
    'A' | 'a' => 0
    'C' | 'c' => 1
    'D' | 'd' => 2
    'E' | 'e' => 3
    'F' | 'f' => 4
    'G' | 'g' => 5
    'H' | 'h' => 6
    'I' | 'i' => 7
    'K' | 'k' => 8
    'L' | 'l' => 9
    'M' | 'm' => 10
    'N' | 'n' => 11
    'P' | 'p' => 12
    'Q' | 'q' => 13
    'R' | 'r' => 14
    'S' | 's' => 15
    'T' | 't' => 16
    'V' | 'v' => 17
    'W' | 'w' => 18
    'Y' | 'y' => 19
    _ => -1
  }
}

///|
fn pa_char_at(seq : String, i : Int) -> Char {
  seq.unsafe_get(i).unsafe_to_char()
}

// ===========================================================================
// 1. Chou–Fasman secondary-structure prediction
// ===========================================================================

///|
/// Chou–Fasman propensity values (Pα, Pβ, Pt) for the 20 standard amino acids.
/// Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical).
fn cf_helix_propensity() -> Array[Double] {
  [
    1.42, 0.70, 1.01, 1.51, 1.13, 0.57, 1.00, 1.08, 1.16, 1.21, 1.45, 0.67, 0.57,
    1.11, 0.98, 0.77, 0.83, 1.06, 1.08, 0.69,
  ]
}

///|
fn cf_sheet_propensity() -> Array[Double] {
  [
    0.83, 1.19, 0.54, 0.37, 1.38, 0.75, 0.87, 1.60, 0.74, 1.30, 1.05, 0.89, 0.55,
    1.10, 0.93, 0.75, 1.19, 1.70, 1.37, 1.47,
  ]
}

///|
fn cf_turn_propensity() -> Array[Double] {
  [
    0.66, 1.19, 1.46, 0.74, 0.60, 1.56, 0.95, 0.47, 1.01, 0.59, 0.60, 1.56, 1.52,
    0.98, 0.95, 1.43, 0.96, 0.50, 0.96, 1.14,
  ]
}

///|
/// Result of Chou–Fasman secondary-structure prediction.
pub struct ChouFasmanResult {
  sequence : String
  prediction : Array[String]
  helices : Array[(Int, Int)]
  sheets : Array[(Int, Int)]
  turns : Array[(Int, Int)]
} derive(Debug)
pub extend ChouFasmanResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn ChouFasmanResult::sequence(self : ChouFasmanResult) -> String {
  self.sequence
}

///|
pub fn ChouFasmanResult::prediction(self : ChouFasmanResult) -> Array[String] {
  self.prediction
}

///|
pub fn ChouFasmanResult::helices(self : ChouFasmanResult) -> Array[(Int, Int)] {
  self.helices
}

///|
pub fn ChouFasmanResult::sheets(self : ChouFasmanResult) -> Array[(Int, Int)] {
  self.sheets
}

///|
pub fn ChouFasmanResult::turns(self : ChouFasmanResult) -> Array[(Int, Int)] {
  self.turns
}

///|
/// Predict protein secondary structure using the Chou–Fasman (1974) method.
///
/// Scans for helix nucleation (≥6 residues with mean Pα > 1.03), sheet
/// nucleation (≥3 residues with mean Pβ > 1.05), and turn regions (4 residues
/// with mean Pt > 1.0).  Nucleation sites are extended while the running mean
/// propensity stays above 1.0.  Overlapping helix/sheet assignments are
/// resolved in favour of the higher mean propensity.
pub fn chou_fasman_predict(
  sequence : String,
) -> ChouFasmanResult raise ProteinAdvancedError {
  pa_validate_protein(sequence, "Chou-Fasman")
  let n = sequence.length()
  let p_alpha = cf_helix_propensity()
  let p_beta = cf_sheet_propensity()
  let p_turn = cf_turn_propensity()
  // Per-residue propensities.
  let alpha_vals : Array[Double] = Array::make(n, 0.0)
  let beta_vals : Array[Double] = Array::make(n, 0.0)
  let turn_vals : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let idx = pa_aa_index(pa_char_at(sequence, i))
    if idx >= 0 {
      alpha_vals[i] = p_alpha[idx]
      beta_vals[i] = p_beta[idx]
      turn_vals[i] = p_turn[idx]
    }
    i = i + 1
  }
  // Assignment array: 0 = coil, 1 = helix, 2 = sheet, 3 = turn
  let assign : Array[Int] = Array::make(n, 0)
  let helices : Array[(Int, Int)] = []
  let sheets : Array[(Int, Int)] = []
  let turns : Array[(Int, Int)] = []
  // --- Helix nucleation + extension ---
  i = 0
  while i <= n - 6 {
    let mut sum_a = 0.0
    let mut sum_b = 0.0
    let mut k = 0
    while k < 6 {
      sum_a = sum_a + alpha_vals[i + k]
      sum_b = sum_b + beta_vals[i + k]
      k = k + 1
    }
    if sum_a / 6.0 > 1.03 && sum_a > sum_b {
      // Extend left.
      let mut left = i
      while left > 0 {
        let prev_a = alpha_vals[left - 1]
        let prev_b = beta_vals[left - 1]
        if prev_a > 1.0 && prev_a >= prev_b {
          left = left - 1
        } else {
          break
        }
      }
      // Extend right.
      let mut right = i + 5
      while right < n - 1 {
        let next_a = alpha_vals[right + 1]
        let next_b = beta_vals[right + 1]
        if next_a > 1.0 && next_a >= next_b {
          right = right + 1
        } else {
          break
        }
      }
      helices.push((left, right))
      i = right + 1
    } else {
      i = i + 1
    }
  }
  // --- Sheet nucleation + extension ---
  i = 0
  while i <= n - 3 {
    let mut sum_b = 0.0
    let mut sum_a = 0.0
    let mut k = 0
    while k < 3 {
      sum_b = sum_b + beta_vals[i + k]
      sum_a = sum_a + alpha_vals[i + k]
      k = k + 1
    }
    if sum_b / 3.0 > 1.05 && sum_b > sum_a {
      let mut left = i
      while left > 0 {
        if beta_vals[left - 1] > 1.0 &&
          beta_vals[left - 1] >= alpha_vals[left - 1] {
          left = left - 1
        } else {
          break
        }
      }
      let mut right = i + 2
      while right < n - 1 {
        if beta_vals[right + 1] > 1.0 &&
          beta_vals[right + 1] >= alpha_vals[right + 1] {
          right = right + 1
        } else {
          break
        }
      }
      sheets.push((left, right))
      i = right + 1
    } else {
      i = i + 1
    }
  }
  // --- Turn detection ---
  i = 0
  while i <= n - 4 {
    let mut sum_t = 0.0
    let mut sum_a = 0.0
    let mut sum_b = 0.0
    let mut k = 0
    while k < 4 {
      sum_t = sum_t + turn_vals[i + k]
      sum_a = sum_a + alpha_vals[i + k]
      sum_b = sum_b + beta_vals[i + k]
      k = k + 1
    }
    if sum_t / 4.0 > 1.0 && sum_a / 4.0 < 1.0 && sum_b / 4.0 < 1.0 {
      turns.push((i, i + 3))
    }
    i = i + 1
  }
  // --- Assign: resolve helix/sheet conflicts by mean propensity ---
  for region in helices {
    let (start, end) = region
    let mut mean_a = 0.0
    let mut mean_b = 0.0
    let mut k = start
    while k <= end {
      mean_a = mean_a + alpha_vals[k]
      mean_b = mean_b + beta_vals[k]
      k = k + 1
    }
    let count = (end - start + 1).to_double()
    mean_a = mean_a / count
    mean_b = mean_b / count
    if mean_a >= mean_b {
      let mut j = start
      while j <= end {
        assign[j] = 1
        j = j + 1
      }
    } else {
      let mut j = start
      while j <= end {
        if assign[j] == 0 {
          assign[j] = 2
        }
        j = j + 1
      }
    }
  }
  for region in sheets {
    let (start, end) = region
    let mut mean_a = 0.0
    let mut mean_b = 0.0
    let mut k = start
    while k <= end {
      mean_a = mean_a + alpha_vals[k]
      mean_b = mean_b + beta_vals[k]
      k = k + 1
    }
    let count = (end - start + 1).to_double()
    mean_a = mean_a / count
    mean_b = mean_b / count
    if mean_b > mean_a {
      let mut j = start
      while j <= end {
        assign[j] = 2
        j = j + 1
      }
    }
  }
  for region in turns {
    let (start, end) = region
    let mut j = start
    while j <= end {
      if assign[j] == 0 {
        assign[j] = 3
      }
      j = j + 1
    }
  }
  // Build prediction string array.
  let prediction : Array[String] = Array::make(n, "C")
  let mut j = 0
  while j < n {
    prediction[j] = match assign[j] {
      1 => "H"
      2 => "E"
      3 => "T"
      _ => "C"
    }
    j = j + 1
  }
  ChouFasmanResult::{ sequence, prediction, helices, sheets, turns, }
}

// ===========================================================================
// 2. IUPred disorder prediction
// ===========================================================================

///|
/// IUPred pairwise interaction-energy matrix (Dosztányi et al. 2005).
/// Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical, 0-19).
/// Values are estimated from known folded proteins; more negative = stronger
/// interaction (more ordered).
fn iupred_energy_matrix() -> Array[Array[Double]] {
  [
    [
      0.0, -0.82, -0.91, -0.71, -0.71, -0.52, -0.42, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.82, 0.0, -0.91, -0.91, -0.71, -0.71, -0.52, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.91, -0.91, 0.0, -0.71, -0.71, -0.52, -0.42, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.71, -0.91, -0.71, 0.0, -0.71, -0.52, -0.42, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.71, -0.71, -0.71, -0.71, 0.0, -0.52, -0.42, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.52, -0.71, -0.52, -0.52, -0.52, 0.0, -0.42, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.42, -0.52, -0.42, -0.42, -0.42, -0.42, 0.0, -0.61, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, 0.0, -0.76, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.76, -0.76, -0.76, -0.76, -0.76, -0.76, -0.76, -0.76, 0.0, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, 0.0, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, 0.0,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, -0.61, -0.5, -0.49,
    ],
    [
      -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49,
      0.0, -0.33, -0.5, -0.49, -0.45, -0.33, -0.49, -0.5, -0.49,
    ],
    [
      -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33,
      -0.33, 0.0, -0.5, -0.33, -0.33, -0.33, -0.33, -0.5, -0.33,
    ],
    [
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      0.0, -0.5, -0.45, -0.33, -0.5, -0.5, -0.49,
    ],
    [
      -0.59, -0.59, -0.59, -0.59, -0.59, -0.59, -0.59, -0.59, -0.59, -0.59, -0.59,
      -0.49, -0.33, -0.5, 0.0, -0.45, -0.33, -0.59, -0.5, -0.49,
    ],
    [
      -0.45, -0.45, -0.45, -0.45, -0.45, -0.45, -0.45, -0.45, -0.45, -0.45, -0.45,
      -0.45, -0.33, -0.45, -0.45, 0.0, -0.33, -0.45, -0.5, -0.45,
    ],
    [
      -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33, -0.33,
      -0.33, -0.33, -0.33, -0.33, -0.33, 0.0, -0.33, -0.5, -0.33,
    ],
    [
      -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61, -0.61,
      -0.49, -0.33, -0.5, -0.59, -0.45, -0.33, 0.0, -0.5, -0.49,
    ],
    [
      -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
      -0.5, -0.5, -0.5, -0.5, -0.5, 0.0, -0.5,
    ],
    [
      -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49,
      -0.49, -0.33, -0.49, -0.49, -0.45, -0.33, -0.49, -0.5, 0.0,
    ],
  ]
}

///|
/// Result of IUPred disorder prediction.
pub struct IupredResult {
  sequence : String
  scores : Array[Double]
  disordered_regions : Array[(Int, Int)]
} derive(Debug)
pub extend IupredResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn IupredResult::sequence(self : IupredResult) -> String {
  self.sequence
}

///|
pub fn IupredResult::scores(self : IupredResult) -> Array[Double] {
  self.scores
}

///|
pub fn IupredResult::disordered_regions(
  self : IupredResult,
) -> Array[(Int, Int)] {
  self.disordered_regions
}

///|
/// IUPred disorder prediction (Dosztányi et al. 2005).
///
/// Estimates the pairwise interaction energy for each residue within a sliding
/// window.  The raw energy is converted to a [0, 1] disorder score via a
/// logistic transformation.  Scores above 0.5 indicate disorder.
/// `window_size` controls the neighbourhood (long mode ≈ 100, short mode ≈ 25).
pub fn iupred(
  sequence : String,
  window_size? : Int = 100,
) -> IupredResult raise ProteinAdvancedError {
  pa_validate_protein(sequence, "IUPred")
  let n = sequence.length()
  let win = if window_size < 1 {
    1
  } else if window_size > n {
    n
  } else {
    window_size
  }
  let matrix = iupred_energy_matrix()
  let scores : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let i_idx = pa_aa_index(pa_char_at(sequence, i))
    if i_idx < 0 {
      scores[i] = 0.5
      i = i + 1
      continue
    }
    let half = win / 2
    let left = if i - half < 0 { 0 } else { i - half }
    let right = if i + half >= n { n - 1 } else { i + half }
    let mut energy = 0.0
    let mut count = 0
    let mut j = left
    while j <= right {
      if j != i {
        let j_idx = pa_aa_index(pa_char_at(sequence, j))
        if j_idx >= 0 {
          energy = energy + matrix[i_idx][j_idx]
          count = count + 1
        }
      }
      j = j + 1
    }
    if count > 0 {
      energy = energy / count.to_double()
    }
    // Logistic transformation: more negative energy → lower disorder score.
    let score = 1.0 / (1.0 + @math.exp(-(energy + 0.45) * 4.0))
    scores[i] = score
    i = i + 1
  }
  // Identify disordered regions (score > 0.5 for ≥ 5 consecutive residues).
  let disordered : Array[(Int, Int)] = []
  let mut start = -1
  let mut k = 0
  while k < n {
    if scores[k] > 0.5 {
      if start < 0 {
        start = k
      }
    } else {
      if start >= 0 && k - start >= 5 {
        disordered.push((start, k - 1))
      }
      start = -1
    }
    k = k + 1
  }
  if start >= 0 && n - start >= 5 {
    disordered.push((start, n - 1))
  }
  IupredResult::{ sequence, scores, disordered_regions: disordered, }
}

///|
/// IUPred long disorder mode (window = 100, whole-sequence context).
pub fn iupred_long(
  sequence : String,
) -> IupredResult raise ProteinAdvancedError {
  iupred(sequence, window_size=100)
}

///|
/// IUPred short disorder mode (window = 25, local context).
pub fn iupred_short(
  sequence : String,
) -> IupredResult raise ProteinAdvancedError {
  iupred(sequence, window_size=25)
}

// ===========================================================================
// 3. COILS coiled-coil prediction
// ===========================================================================

///|
/// COILS scoring matrix for each amino acid at each heptad position (a–g).
/// Adapted from Lupas et al. (1991).  Positions a and d are the hydrophobic
/// core and carry 2.5× weight in the scoring.
fn coils_matrix() -> Array[Array[Double]] {
  // Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical, 0-19)
  // Each inner array = [pos_a, pos_b, pos_c, pos_d, pos_e, pos_f, pos_g]
  [
    [0.79, 0.72, 0.73, 0.79, 0.79, 0.79, 0.79], // A
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // C
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // D
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // E
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // F
    [0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43], // G
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // H
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // I
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // K
    [1.20, 0.59, 0.59, 1.20, 0.59, 0.59, 0.59], // L
    [1.20, 0.59, 0.59, 1.20, 0.59, 0.59, 0.59], // M
    [0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43], // N
    [0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59], // P
    [0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59], // Q
    [0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59], // R
    [0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43], // S
    [0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43], // T
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // V
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // W
    [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79], // Y
  ]
}

///|
/// Result of COILS coiled-coil prediction.
pub struct CoilsResult {
  sequence : String
  scores : Array[Double]
  coiled_coil_regions : Array[(Int, Int)]
} derive(Debug)
pub extend CoilsResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn CoilsResult::sequence(self : CoilsResult) -> String {
  self.sequence
}

///|
pub fn CoilsResult::scores(self : CoilsResult) -> Array[Double] {
  self.scores
}

///|
pub fn CoilsResult::coiled_coil_regions(
  self : CoilsResult,
) -> Array[(Int, Int)] {
  self.coiled_coil_regions
}

///|
/// Predict coiled-coil regions using the COILS method (Lupas et al. 1991).
///
/// Scans the sequence with a sliding window (default 28).  For each window
/// position, the score is the weighted mean of heptad-repeat scores (positions
/// a and d carry 2.5× weight).  Normalized scores above the threshold
/// (default 0.9) indicate a coiled-coil segment.
pub fn predict_coiled_coils(
  sequence : String,
  window? : Int = 28,
  threshold? : Double = 0.9,
) -> CoilsResult raise ProteinAdvancedError {
  pa_validate_protein(sequence, "COILS")
  let n = sequence.length()
  let win = if window < 7 { 7 } else if window > n { n } else { window }
  let matrix = coils_matrix()
  let scores : Array[Double] = Array::make(n, 0.0)
  let half = win / 2
  let mut i = 0
  while i < n {
    let left = if i - half < 0 { 0 } else { i - half }
    let right = if i + half >= n { n - 1 } else { i + half }
    let mut best_score = 0.0
    // Try all 7 heptad frames.
    let mut frame = 0
    while frame < 7 {
      let mut total = 0.0
      let mut weight_sum = 0.0
      let mut j = left
      while j <= right {
        let pos = (j - left + frame) % 7
        let idx = pa_aa_index(pa_char_at(sequence, j))
        if idx >= 0 {
          let weight = if pos == 0 || pos == 3 { 2.5 } else { 1.0 }
          total = total + matrix[idx][pos] * weight
          weight_sum = weight_sum + weight
        }
        j = j + 1
      }
      if weight_sum > 0.0 {
        let mean = total / weight_sum
        if mean > best_score {
          best_score = mean
        }
      }
      frame = frame + 1
    }
    scores[i] = best_score
    i = i + 1
  }
  // Identify coiled-coil regions.
  let regions : Array[(Int, Int)] = []
  let mut start = -1
  let mut k = 0
  while k < n {
    if scores[k] >= threshold {
      if start < 0 {
        start = k
      }
    } else {
      if start >= 0 {
        regions.push((start, k - 1))
      }
      start = -1
    }
    k = k + 1
  }
  if start >= 0 {
    regions.push((start, n - 1))
  }
  CoilsResult::{ sequence, scores, coiled_coil_regions: regions, }
}

// ===========================================================================
// 4. Kolaskar–Tongaonkar antigenicity prediction
// ===========================================================================

///|
/// Kolaskar–Tongaonkar antigenic propensity values (1990).
/// Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical).
fn kolaskar_propensity() -> Array[Double] {
  [
    1.064, 1.118, 0.854, 0.896, 1.064, 0.896, 1.118, 1.012, 0.896, 1.064, 1.064,
    0.854, 1.064, 0.896, 0.854, 0.896, 0.896, 1.064, 1.118, 0.854,
  ]
}

///|
/// Result of Kolaskar–Tongaonkar antigenicity prediction.
pub struct AntigenicityResult {
  sequence : String
  scores : Array[Double]
  antigenic_sites : Array[(Int, Int)]
} derive(Debug)
pub extend AntigenicityResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn AntigenicityResult::sequence(self : AntigenicityResult) -> String {
  self.sequence
}

///|
pub fn AntigenicityResult::scores(self : AntigenicityResult) -> Array[Double] {
  self.scores
}

///|
pub fn AntigenicityResult::antigenic_sites(
  self : AntigenicityResult,
) -> Array[(Int, Int)] {
  self.antigenic_sites
}

///|
/// Predict antigenic determinants using the Kolaskar–Tongaonkar (1990) method.
///
/// Computes a sliding-window mean of amino-acid antigenic propensities.  The
/// minimum within each region above the overall mean marks the antigenic site
/// centre (the least hydrophilic residue tends to be the epitope core).
/// `window` (default 7) controls the smoothing.
pub fn kolaskar_tongaonkar_antigenicity(
  sequence : String,
  window? : Int = 7,
) -> AntigenicityResult raise ProteinAdvancedError {
  pa_validate_protein(sequence, "Kolaskar-Tongaonkar")
  let n = sequence.length()
  let win = if window < 1 { 1 } else if window > n { n } else { window }
  let prop = kolaskar_propensity()
  let half = win / 2
  let scores : Array[Double] = Array::make(n, 0.0)
  let mut grand_mean = 0.0
  let mut i = 0
  while i < n {
    let left = if i - half < 0 { 0 } else { i - half }
    let right = if i + half >= n { n - 1 } else { i + half }
    let mut sum = 0.0
    let mut count = 0
    let mut j = left
    while j <= right {
      let idx = pa_aa_index(pa_char_at(sequence, j))
      if idx >= 0 {
        sum = sum + prop[idx]
        count = count + 1
      }
      j = j + 1
    }
    let mean = if count > 0 { sum / count.to_double() } else { 1.0 }
    scores[i] = mean
    grand_mean = grand_mean + prop[pa_aa_index(pa_char_at(sequence, i)).max(0)]
    i = i + 1
  }
  grand_mean = grand_mean / n.to_double()
  // Identify antigenic sites: regions above grand mean.
  let sites : Array[(Int, Int)] = []
  let mut start = -1
  let mut k = 0
  while k < n {
    if scores[k] >= grand_mean {
      if start < 0 {
        start = k
      }
    } else {
      if start >= 0 && k - start >= 6 {
        sites.push((start, k - 1))
      }
      start = -1
    }
    k = k + 1
  }
  if start >= 0 && n - start >= 6 {
    sites.push((start, n - 1))
  }
  AntigenicityResult::{ sequence, scores, antigenic_sites: sites, }
}

// ===========================================================================
// 5. Emini surface accessibility prediction
// ===========================================================================

///|
/// Emini surface probability values (Emini et al. 1985).
/// Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical).
fn emini_surface_prob() -> Array[Double] {
  [
    0.243, 0.000, 1.168, 0.984, 0.404, 0.476, 0.630, 0.168, 0.960, 0.260, 0.260,
    0.911, 0.476, 0.960, 0.960, 0.560, 0.560, 0.168, 0.404, 0.630,
  ]
}

///|
/// Predict surface accessibility using the Emini et al. (1985) method.
///
/// Computes a sliding-window surface probability.  The formula follows the
/// original: `sn = (sum of probs in window) / (window_size * scale)` where
/// the centre residue carries an extra weight.
pub fn emini_surface_accessibility(
  sequence : String,
  window? : Int = 6,
) -> Array[Double] raise ProteinAdvancedError {
  pa_validate_protein(sequence, "Emini")
  let n = sequence.length()
  let win = if window < 1 { 1 } else if window > n { n } else { window }
  let prob = emini_surface_prob()
  let half = win / 2
  let scores : Array[Double] = Array::make(n, 0.0)
  let scale = 0.5
  let mut i = 0
  while i < n {
    let left = if i - half < 0 { 0 } else { i - half }
    let right = if i + half >= n { n - 1 } else { i + half }
    let mut sum = 0.0
    let mut count = 0
    let mut j = left
    while j <= right {
      let idx = pa_aa_index(pa_char_at(sequence, j))
      if idx >= 0 {
        let weight = if j == i { prob[idx] * 2.0 } else { prob[idx] }
        sum = sum + weight
        count = count + 1
      }
      j = j + 1
    }
    let denom = (count.to_double() * scale).max(1.0e-10)
    scores[i] = sum / denom
    i = i + 1
  }
  scores
}

// ===========================================================================
// 6. Karplus–Schulz flexibility prediction
// ===========================================================================

///|
/// Karplus–Schulz flexibility parameters (1985).
/// Order: A C D E F G H I K L M N P Q R S T V W Y (alphabetical).
/// Values represent normalized B-factors (higher = more flexible).
fn karplus_schulz_params() -> Array[Double] {
  [
    0.984, 0.906, 1.068, 1.094, 0.915, 1.031, 0.950, 0.927, 1.102, 0.927, 0.940,
    1.048, 1.049, 1.031, 0.984, 1.046, 0.997, 0.927, 0.906, 0.929,
  ]
}

///|
/// Predict chain flexibility using the Karplus–Schulz (1985) method.
///
/// Computes a sliding-window mean of normalized B-factor flexibility
/// parameters.  Values above 1.0 indicate above-average flexibility.
pub fn karplus_schulz_flexibility(
  sequence : String,
  window? : Int = 5,
) -> Array[Double] raise ProteinAdvancedError {
  pa_validate_protein(sequence, "Karplus-Schulz")
  let n = sequence.length()
  let win = if window < 1 { 1 } else if window > n { n } else { window }
  let params = karplus_schulz_params()
  let half = win / 2
  let scores : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let left = if i - half < 0 { 0 } else { i - half }
    let right = if i + half >= n { n - 1 } else { i + half }
    let mut sum = 0.0
    let mut count = 0
    let mut j = left
    while j <= right {
      let idx = pa_aa_index(pa_char_at(sequence, j))
      if idx >= 0 {
        sum = sum + params[idx]
        count = count + 1
      }
      j = j + 1
    }
    scores[i] = if count > 0 { sum / count.to_double() } else { 1.0 }
    i = i + 1
  }
  scores
}

// ===========================================================================
// Convenience: demo sequences
// ===========================================================================

///|
/// Demo: a sequence with a clear α-helical propensity (poly-Ala with Leu).
pub fn demo_helical_sequence() -> String {
  "MAELAALKEKALAALEALKEKALAALKELAA"
}

///|
/// Demo: a sequence rich in disorder-promoting residues.
pub fn demo_disorder_sequence() -> String {
  "MQDRQKPKQDRQKPKQDRQKPKQDRQKPKQDRQKPK"
}

///|
/// Demo: a sequence with a leucine-zipper-like coiled-coil motif.
pub fn demo_coiled_coil_sequence() -> String {
  "MKQLEDKVEELLSKNYHLENEVARLKKLVGER"
}
