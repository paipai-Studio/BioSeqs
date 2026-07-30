///|
/// Bio.SeqUtils.ProtDao module - Protein Disorder prediction.
/// Inspired by Biopython Bio.SeqUtils.ProtDao (IUPred).
///
/// Implements the IUPred algorithm for predicting intrinsically
/// disordered regions in protein sequences. The method uses
/// position-specific amino acid scores and a smoothing window
/// to classify each residue as ordered or disordered.
///
/// References:
///   - https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtDao.html
///   - Dosztanyi Z, Csizmok V, Tompa P, Simon I (2005)
///     IUPred: web server for the prediction of intrinsically
///     unstructured regions of proteins based on estimated energy
///     content. Bioinformatics 21(16):3433-3434.

// ============================================================================
// Amino Acid Disorder Scores (IUPred-like)
// ============================================================================

///|
/// Get the IUPred disorder score for an amino acid.
/// Higher values indicate higher propensity for disorder.
/// Based on the IUPred algorithm's amino acid propensity scale.
pub fn prot_dao_disorder_score(aa : UInt16) -> Double {
  let c = if aa >= 97 && aa <= 122 { aa - 32 } else { aa }
  match c {
    65 => 0.06  // 'A'
    67 => -0.02  // 'C'
    68 => 0.19  // 'D'
    69 => 0.74  // 'E'
    70 => -0.09  // 'F'
    71 => 0.16  // 'G'
    72 => -0.42  // 'H'
    73 => -0.48  // 'I'
    75 => 0.58  // 'K'
    76 => -0.03  // 'L'
    77 => -0.39  // 'M'
    78 => 0.03  // 'N'
    80 => 0.98  // 'P'
    81 => 0.31  // 'Q'
    82 => 0.41  // 'R'
    83 => 0.34  // 'S'
    84 => 0.59  // 'T'
    86 => -0.26  // 'V'
    87 => -0.19  // 'W'
    89 => -0.49  // 'Y'
    _ => 0.0
  }
}

///|
/// Get the IUPred energy score for an amino acid.
/// Based on the IUPred pairwise energy approximation.
pub fn prot_dao_energy_score(aa : UInt16) -> Double {
  let c = if aa >= 97 && aa <= 122 { aa - 32 } else { aa }
  match c {
    65 => 1.1471  // 'A'
    67 => 2.6497  // 'C'
    68 => 1.9748  // 'D'
    69 => 1.3577  // 'E'
    70 => 2.0289  // 'F'
    71 => 1.8218  // 'G'
    72 => -1.0930  // 'H'
    73 => 0.2643  // 'I'
    75 => 0.6636  // 'K'
    76 => 0.5035  // 'L'
    77 => 2.5789  // 'M'
    78 => 1.8929  // 'N'
    80 => 2.8110  // 'P'
    81 => 1.1740  // 'Q'
    82 => -0.3509  // 'R'
    83 => 2.3317  // 'S'
    84 => 0.7361  // 'T'
    86 => 0.2730  // 'V'
    87 => 1.5749  // 'W'
    89 => 1.6374  // 'Y'
    _ => 0.0
  }
}

// ============================================================================
// Data Structures
// ============================================================================

///|
/// A single residue disorder prediction result.
pub struct DisorderResidue {
  position : Int
  amino_acid : UInt16
  mut disorder_score : Double
  energy_score : Double
  mut is_disordered : Bool
  mut is_disordered_long : Bool
}

///|
/// Create a disorder residue result.
pub fn DisorderResidue::new(
  position : Int,
  amino_acid : UInt16,
  disorder_score : Double,
  energy_score : Double,
) -> DisorderResidue {
  DisorderResidue::{
    position,
    amino_acid,
    disorder_score,
    energy_score,
    is_disordered: false,
    is_disordered_long: false,
  }
}

///|
/// A region of disorder in the protein sequence.
pub struct DisorderRegion {
  start : Int
  end_ : Int
  length : Int
  avg_score : Double
  max_score : Double
  region_type : String
}

///|
/// Create a disorder region.
pub fn DisorderRegion::new(
  start : Int,
  end_ : Int,
  avg_score : Double,
  max_score : Double,
  region_type? : String = "disordered",
) -> DisorderRegion {
  DisorderRegion::{
    start,
    end_,
    length: end_ - start,
    avg_score,
    max_score,
    region_type,
  }
}

///|
/// Full disorder prediction result for a protein.
pub struct DisorderResult {
  sequence : String
  residues : Array[DisorderResidue]
  regions : Array[DisorderRegion]
  n_disordered : Int
  n_disordered_long : Int
  fraction_disordered : Double
  threshold_disordered : Double
  threshold_long : Double
  method : String
}

// ============================================================================
// Core Prediction Functions
// ============================================================================

///|
/// Predict disorder for a protein sequence using IUPred-like algorithm.
///
/// # Parameters
/// - sequence: The protein sequence string
/// - threshold_disordered: Threshold for short disordered regions (default 0.5)
/// - threshold_long: Threshold for long disordered regions (default 0.686)
///
/// # Returns
/// DisorderResult with per-residue scores and predicted regions.
pub fn prot_dao_predict(
  sequence : String,
  threshold_disordered? : Double = 0.5,
  threshold_long? : Double = 0.686,
) -> DisorderResult {
  let seq_len = sequence.length()
  let residues : Array[DisorderResidue] = Array::new()

  // Step 1: Compute raw disorder scores per residue
  let mut i = 0
  while i < seq_len {
    let aa = sequence.unsafe_get(i)
    let ds = prot_dao_disorder_score(aa)
    let es = prot_dao_energy_score(aa)
    residues.push(DisorderResidue::new(i + 1, aa, ds, es))
    i = i + 1
  }

  // Step 2: Smoothing window (IUPred uses a specific window size)
  let window = if seq_len < 21 { seq_len } else { 21 }
  let half_w = window / 2

  // Compute smoothed scores
  let smoothed : Array[Double] = Array::new()
  let mut si = 0
  while si < seq_len {
    let start = if si - half_w < 0 { 0 } else { si - half_w }
    let end_ = if si + half_w >= seq_len { seq_len - 1 } else { si + half_w }
    let mut sum = 0.0
    let mut count = 0
    let mut j = start
    while j <= end_ {
      sum = sum + residues[j].disorder_score
      count = count + 1
      j = j + 1
    }
    smoothed.push(if count > 0 { sum / count.to_double() } else { 0.0 })
    si = si + 1
  }

  // Step 3: Classify each residue
  let mut n_disordered = 0
  let mut n_disordered_long = 0
  let mut ci = 0
  while ci < seq_len {
    residues[ci].disorder_score = smoothed[ci]
    residues[ci].is_disordered = smoothed[ci] > threshold_disordered
    residues[ci].is_disordered_long = smoothed[ci] > threshold_long
    if residues[ci].is_disordered { n_disordered = n_disordered + 1 }
    if residues[ci].is_disordered_long { n_disordered_long = n_disordered_long + 1 }
    ci = ci + 1
  }

  // Step 4: Find contiguous disordered regions
  let regions : Array[DisorderRegion] = Array::new()
  let mut ri = 0
  while ri < seq_len {
    if residues[ri].is_disordered {
      // Start of a disordered region
      let region_start = ri
      let mut region_end = ri
      let mut rs = 0.0
      let mut count = 0
      let mut max_s = 0.0

      while region_end < seq_len && residues[region_end].is_disordered {
        let s = smoothed[region_end]
        rs = rs + s
        count = count + 1
        if s > max_s { max_s = s }
        region_end = region_end + 1
      }

      let avg_s = if count > 0 { rs / count.to_double() } else { 0.0 }
      let is_long = region_end - region_start >= 30
      let type_ = if is_long { "long disordered" } else { "disordered" }

      regions.push(DisorderRegion::new(
        region_start + 1,
        region_end,
        avg_s,
        max_s,
        region_type=type_,
      ))
      ri = region_end
    } else {
      ri = ri + 1
    }
  }

  let frac = if seq_len > 0 {
    n_disordered.to_double() / seq_len.to_double()
  } else {
    0.0
  }

  DisorderResult::{
    sequence,
    residues,
    regions,
    n_disordered,
    n_disordered_long,
    fraction_disordered: frac,
    threshold_disordered,
    threshold_long,
    method: "IUPred-like",
  }
}

// ============================================================================
// DisorderResult Methods
// ============================================================================

///|
/// Get number of disordered residues.
pub fn DisorderResult::get_n_disordered(self : DisorderResult) -> Int {
  self.n_disordered
}

///|
/// Get number of long disordered residues.
pub fn DisorderResult::get_n_disordered_long(self : DisorderResult) -> Int {
  self.n_disordered_long
}

///|
/// Get fraction of disordered residues.
pub fn DisorderResult::get_fraction_disordered(self : DisorderResult) -> Double {
  self.fraction_disordered
}

///|
/// Get disordered regions.
pub fn DisorderResult::get_regions(self : DisorderResult) -> Array[DisorderRegion] {
  self.regions.copy()
}

///|
/// Get number of disordered regions.
pub fn DisorderResult::get_n_regions(self : DisorderResult) -> Int {
  self.regions.length()
}

///|
/// Get the longest disordered region.
pub fn DisorderResult::get_longest_region(self : DisorderResult) -> DisorderRegion {
  if self.regions.length() == 0 {
    return DisorderRegion::new(0, 0, 0.0, 0.0, region_type="none")
  }
  let mut best = self.regions[0]
  let mut i = 1
  while i < self.regions.length() {
    if self.regions[i].length > best.length {
      best = self.regions[i]
    }
    i = i + 1
  }
  best
}

///|
/// Get disordered residues as a string.
pub fn DisorderResult::disordered_sequence(self : DisorderResult) -> String {
  let mut s = ""
  let mut i = 0
  while i < self.residues.length() {
    if self.residues[i].is_disordered {
      s = s + self.residues[i].amino_acid.to_char().to_string()
    } else {
      s = s + "-"
    }
    i = i + 1
  }
  s
}

///|
/// Get disorder score array.
pub fn DisorderResult::get_scores(self : DisorderResult) -> Array[Double] {
  let scores : Array[Double] = Array::new()
  let mut i = 0
  while i < self.residues.length() {
    scores.push(self.residues[i].disorder_score)
    i = i + 1
  }
  scores
}

///|
/// Summary of disorder prediction.
pub fn DisorderResult::summary(self : DisorderResult) -> String {
  let mut s = "Disorder Prediction Summary\n"
  s = s + "Method: " + self.method + "\n"
  s = s + "Sequence length: " + self.sequence.length().to_string() + "\n"
  s = s + "Disordered residues: " + self.n_disordered.to_string() + " (" + self.fraction_disordered.to_string() + ")\n"
  s = s + "Long disordered residues: " + self.n_disordered_long.to_string() + "\n"
  s = s + "Disordered regions: " + self.regions.length().to_string() + "\n"
  s = s + "Threshold (disordered): " + self.threshold_disordered.to_string() + "\n"
  s = s + "Threshold (long): " + self.threshold_long.to_string() + "\n"
  s
}

///|
/// Format disorder result as ASCII visualization.
pub fn DisorderResult::to_ascii(self : DisorderResult, width? : Int = 80) -> String {
  let w = if width > 0 { width } else { 80 }
  let seq_len = self.sequence.length()
  if seq_len == 0 { return "Empty sequence" }

  // Create the score track
  let mut s = "Score: "
  let mut i = 0
  while i < seq_len {
    let score = self.residues[i].disorder_score
    let bar = if score > self.threshold_long { "##" } else if score > self.threshold_disordered { "##" } else { "--" }
    s = s + bar
    i = i + 1
  }
  s = s + "\n"

  // Position ruler
  let mut ruler = "Pos:   "
  let mut pos = 1
  while pos <= seq_len {
    let label = pos.to_string()
    let label_len = label.length()
    if pos == 1 || pos % 10 == 0 {
      if pos > 1 {
        let spaces = "          "
        let take = if spaces.length() >= label_len { label_len } else { spaces.length() }
        ruler = ruler + spaces.substring(start=0, end=take)
      }
      ruler = ruler + label
    } else {
      ruler = ruler + " "
    }
    pos = pos + 1
  }
  s = s + ruler + "\n"

  // Sequence
  s = s + "Seq:   " + self.sequence + "\n"

  // Disorder track
  s = s + "Dis:   "
  let mut j = 0
  while j < seq_len {
    if self.residues[j].is_disordered {
      s = s + "*"
    } else {
      s = s + " "
    }
    j = j + 1
  }
  s = s + "\n"

  // Region summary
  s = s + "\nDisordered Regions:\n"
  let mut ri = 0
  while ri < self.regions.length() {
    let reg = self.regions[ri]
    s = s + "  Region " + (ri + 1).to_string() + ": " + reg.start.to_string() + "-" + reg.end_.to_string() + " (length=" + reg.length.to_string() + ", avg_score=" + reg.avg_score.to_string() + ", type=" + reg.region_type + ")\n"
    ri = ri + 1
  }

  s
}

// ============================================================================
// Sample Data
// ============================================================================

///|
/// Generate a sample protein sequence for disorder prediction.
pub fn prot_dao_sample_sequence() -> String {
  // A protein with mixed ordered and disordered regions
  // N-terminal disordered region, central ordered region, C-terminal disordered
  let n_term_disordered = "MSEGKSSGGSSSGRRSSSGGGSSGGKSGKSSGR"  // disordered
  let ordered_core = "PKEIVLDSNGSLEVVKNGKVVVESTGNKVVLTSGQE"  // ordered
  let linker = "SSSSSS"  // disordered linker
  let second_domain = "KLPNLVISGDAGHVVVDAESNGKIVLLPSN"  // ordered
  let c_term_disordered = "GSSGSSGSGGKSGSSGSGRGSSGGNSSGKSSG"  // disordered
  n_term_disordered + ordered_core + linker + second_domain + c_term_disordered
}

///|
/// Generate sample disorder prediction result.
pub fn prot_dao_sample() -> DisorderResult {
  let seq = prot_dao_sample_sequence()
  prot_dao_predict(seq)
}
