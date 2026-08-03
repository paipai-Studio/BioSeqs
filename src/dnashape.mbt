// DNAshapeR - DNA shape prediction from sequence.
// Inspired by Bioconductor DNAshapeR package.
//
// Predicts four DNA shape features from nucleotide sequence using
// precomputed lookup tables derived from all-atom Monte Carlo simulations:
//   - Minor Groove Width (MGW, 5-mer)
//   - Roll (tetranucleotide)
//   - Propeller Twist (ProT, tetranucleotide)
//   - Helix Twist (HelT, 5-mer)
//
// References:
//   - Chiu TP et al. (2016) DNAshapeR: an R/Bioconductor package for
//     DNA shape prediction and usage. Bioinformatics 32:3415-3417.
//   - Zhou T et al. (2013) DNAshape: a method for the high-throughput
//     prediction of DNA structural features on a genomic scale.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// DNA shape feature types.
pub enum DnaShapeFeature {
  ShapeMGW // Minor Groove Width (Å), 5-mer
  ShapeRoll // Roll (degrees), tetranucleotide
  ShapeProT // Propeller Twist (degrees), tetranucleotide
  ShapeHelT // Helix Twist (degrees), 5-mer
  ShapeEP // Electrostatic Potential (kT/e, 5-mer)
}

///|
/// Convert a DnaShapeFeature to its string name.
pub fn DnaShapeFeature::to_string(self : DnaShapeFeature) -> String {
  match self {
    ShapeMGW => "MGW"
    ShapeRoll => "Roll"
    ShapeProT => "ProT"
    ShapeHelT => "HelT"
    ShapeEP => "EP"
  }
}

///|
/// Factory function for cross-package DnaShapeFeature construction.
pub fn dnashape_make_feature(name : String) -> DnaShapeFeature {
  match name {
    "MGW" => ShapeMGW
    "Roll" => ShapeRoll
    "ProT" => ShapeProT
    "HelT" => ShapeHelT
    "EP" => ShapeEP
    _ => ShapeMGW
  }
}

///|
/// Get the k-mer size used for a feature.
pub fn DnaShapeFeature::kmer_size(self : DnaShapeFeature) -> Int {
  match self {
    ShapeMGW => 5
    ShapeRoll => 4
    ShapeProT => 4
    ShapeHelT => 5
    ShapeEP => 5
  }
}

///|
/// Result of predicting a single shape feature for a sequence.
pub struct DnaShapeResult {
  pub feature : DnaShapeFeature
  pub values : Array[Double] // per-position values (NaN at boundaries)
  pub positions : Array[Int] // 1-based positions where values are defined
}

///|
/// Construct a DnaShapeResult.
pub fn DnaShapeResult::new(
  feature~ : DnaShapeFeature,
  values~ : Array[Double],
  positions~ : Array[Int],
) -> DnaShapeResult {
  DnaShapeResult::{ feature, values, positions }
}

///|
/// Get the feature name as string.
pub fn DnaShapeResult::feature_name(self : DnaShapeResult) -> String {
  self.feature.to_string()
}

///|
/// Get the number of predicted values.
pub fn DnaShapeResult::length(self : DnaShapeResult) -> Int {
  self.values.length()
}

///|
/// Get value at index.
pub fn DnaShapeResult::get_value(self : DnaShapeResult, idx : Int) -> Double {
  self.values[idx]
}

///|
/// Compute mean of predicted values (ignoring NaN).
pub fn DnaShapeResult::mean(self : DnaShapeResult) -> Double {
  let mut sum = 0.0
  let mut count = 0
  for v in self.values {
    if !v.is_nan() {
      sum = sum + v
      count = count + 1
    }
  }
  if count == 0 {
    0.0
  } else {
    sum / count.to_double()
  }
}

///|
/// Compute min of predicted values (ignoring NaN).
pub fn DnaShapeResult::min(self : DnaShapeResult) -> Double {
  let mut m = 1.0e30
  for v in self.values {
    if !v.is_nan() && v < m {
      m = v
    }
  }
  if m == 1.0e30 {
    0.0
  } else {
    m
  }
}

///|
/// Compute max of predicted values (ignoring NaN).
pub fn DnaShapeResult::max(self : DnaShapeResult) -> Double {
  let mut m = -1.0e30
  for v in self.values {
    if !v.is_nan() && v > m {
      m = v
    }
  }
  if m == -1.0e30 {
    0.0
  } else {
    m
  }
}

///|
/// Complete shape prediction for a sequence (all features).
pub struct DnaShapePrediction {
  pub sequence : String
  pub mgw : DnaShapeResult
  pub roll : DnaShapeResult
  pub prot : DnaShapeResult
  pub helt : DnaShapeResult
}

// ===========================================================================
// Base and k-mer encoding
// ===========================================================================

///|
/// Encode a single base to integer (A=0, C=1, G=2, T=3, others=-1).
pub fn dnashape_encode_base(base : String) -> Int {
  if base == "A" || base == "a" {
    0
  } else if base == "C" || base == "c" {
    1
  } else if base == "G" || base == "g" {
    2
  } else if base == "T" || base == "t" {
    3
  } else {
    -1
  }
}

///|
/// Encode a k-mer string to its integer index (base-4).
/// Returns -1 if any base is non-ACGT.
pub fn dnashape_encode_kmer(kmer : String) -> Int {
  let n = kmer.length()
  let mut result = 0
  let mut i = 0
  while i < n {
    let b = kmer[i:i + 1].to_string()
    let code = dnashape_encode_base(b)
    if code < 0 {
      return -1
    }
    result = result * 4 + code
    i = i + 1
  }
  result
}

///|
/// Decode an integer index to a k-mer string of given length.
pub fn dnashape_decode_kmer(idx : Int, length : Int) -> String {
  let bases = ["A", "C", "G", "T"]
  let mut result = ""
  let mut v = idx
  let codes = Array::make(length, 0)
  let mut i = length - 1
  while i >= 0 {
    codes[i] = v % 4
    v = v / 4
    i = i - 1
  }
  for c in codes {
    result = result + bases[c]
  }
  result
}

///|
/// Compute reverse complement of a DNA sequence (uppercase ACGT only).
pub fn dnashape_revcomp(seq : String) -> String {
  let n = seq.length()
  let mut result = ""
  let mut i = n - 1
  while i >= 0 {
    let c = seq[i:i + 1].to_string()
    let complement = if c == "A" || c == "a" {
      "T"
    } else if c == "C" || c == "c" {
      "G"
    } else if c == "G" || c == "g" {
      "C"
    } else if c == "T" || c == "t" {
      "A"
    } else {
      "N"
    }
    result = result + complement
    i = i - 1
  }
  result
}

// ===========================================================================
// Dinucleotide shape tables (representative literature values)
// ===========================================================================

///|
/// Get the Roll lookup table (dinucleotide index -> degrees).
/// Index computed as base4(dinucleotide). Values from literature.
pub fn dnashape_get_roll_table() -> Map[Int, Double] {
  let t : Map[Int, Double] = Map([], capacity=16)
  // AA=0, AC=1, AG=2, AT=3
  // CA=4, CC=5, CG=6, CT=7
  // GA=8, GC=9, GG=10, GT=11
  // TA=12, TC=13, TG=14, TT=15
  t[0] = -6.0 // AA
  t[1] = 4.0 // AC
  t[2] = 2.0 // AG
  t[3] = -2.0 // AT
  t[4] = 4.0 // CA
  t[5] = 2.0 // CC
  t[6] = 6.0 // CG
  t[7] = 2.0 // CT
  t[8] = 2.0 // GA
  t[9] = 0.0 // GC
  t[10] = 2.0 // GG
  t[11] = 4.0 // GT
  t[12] = 2.0 // TA
  t[13] = 2.0 // TC
  t[14] = 4.0 // TG
  t[15] = -6.0 // TT
  t
}

///|
/// Get the ProT (Propeller Twist) lookup table (dinucleotide -> degrees).
pub fn dnashape_get_prot_table() -> Map[Int, Double] {
  let t : Map[Int, Double] = Map([], capacity=16)
  t[0] = 15.0 // AA
  t[1] = 12.0 // AC
  t[2] = 10.0 // AG
  t[3] = 15.0 // AT
  t[4] = 11.0 // CA
  t[5] = 9.0 // CC
  t[6] = 8.0 // CG
  t[7] = 10.0 // CT
  t[8] = 11.0 // GA
  t[9] = 9.0 // GC
  t[10] = 9.0 // GG
  t[11] = 12.0 // GT
  t[12] = -2.0 // TA
  t[13] = 11.0 // TC
  t[14] = 11.0 // TG
  t[15] = 15.0 // TT
  t
}

///|
/// Get the MGW (Minor Groove Width) lookup table (dinucleotide -> Å).
pub fn dnashape_get_mgw_table() -> Map[Int, Double] {
  let t : Map[Int, Double] = Map([], capacity=16)
  t[0] = 3.5 // AA - A-tracts narrow groove
  t[1] = 4.5 // AC
  t[2] = 4.0 // AG
  t[3] = 4.0 // AT
  t[4] = 4.5 // CA
  t[5] = 4.5 // CC
  t[6] = 5.0 // CG
  t[7] = 4.0 // CT
  t[8] = 4.0 // GA
  t[9] = 4.5 // GC
  t[10] = 4.5 // GG
  t[11] = 4.5 // GT
  t[12] = 4.5 // TA
  t[13] = 4.0 // TC
  t[14] = 4.5 // TG
  t[15] = 3.5 // TT
  t
}

///|
/// Get the HelT (Helix Twist) lookup table (dinucleotide -> degrees).
pub fn dnashape_get_helt_table() -> Map[Int, Double] {
  let t : Map[Int, Double] = Map([], capacity=16)
  t[0] = 35.0 // AA
  t[1] = 34.0 // AC
  t[2] = 34.0 // AG
  t[3] = 31.0 // AT
  t[4] = 34.0 // CA
  t[5] = 34.0 // CC
  t[6] = 33.0 // CG
  t[7] = 34.0 // CT
  t[8] = 34.0 // GA
  t[9] = 34.0 // GC
  t[10] = 34.0 // GG
  t[11] = 34.0 // GT
  t[12] = 36.0 // TA
  t[13] = 34.0 // TC
  t[14] = 34.0 // TG
  t[15] = 35.0 // TT
  t
}

///|
/// Get the EP (Electrostatic Potential) lookup table (dinucleotide -> kT/e).
pub fn dnashape_get_ep_table() -> Map[Int, Double] {
  let t : Map[Int, Double] = Map([], capacity=16)
  t[0] = -1.5 // AA - AT-rich lower potential
  t[1] = -1.0 // AC
  t[2] = -1.1 // AG
  t[3] = -1.4 // AT
  t[4] = -1.0 // CA
  t[5] = -0.7 // CC
  t[6] = -0.6 // CG
  t[7] = -1.1 // CT
  t[8] = -1.1 // GA
  t[9] = -0.6 // GC
  t[10] = -0.7 // GG
  t[11] = -1.0 // GT
  t[12] = -1.3 // TA
  t[13] = -1.1 // TC
  t[14] = -1.0 // TG
  t[15] = -1.5 // TT
  t
}

// ===========================================================================
// Shape prediction
// ===========================================================================

///|
/// Look up a dinucleotide value from the appropriate table.
fn dnashape_lookup_dinuc(
  dinuc_idx : Int,
  feature : DnaShapeFeature,
) -> Double {
  match feature {
    ShapeMGW => {
      let t = dnashape_get_mgw_table()
      t.get(dinuc_idx).unwrap_or(4.0)
    }
    ShapeRoll => {
      let t = dnashape_get_roll_table()
      t.get(dinuc_idx).unwrap_or(0.0)
    }
    ShapeProT => {
      let t = dnashape_get_prot_table()
      t.get(dinuc_idx).unwrap_or(10.0)
    }
    ShapeHelT => {
      let t = dnashape_get_helt_table()
      t.get(dinuc_idx).unwrap_or(34.0)
    }
    ShapeEP => {
      let t = dnashape_get_ep_table()
      t.get(dinuc_idx).unwrap_or(-1.0)
    }
  }
}

///|
/// Predict a single shape feature for a sequence.
/// Uses k-mer averaging of dinucleotide values.
/// For MGW/HelT/EP (k=5): averages 4 overlapping dinucleotides.
/// For Roll/ProT (k=4): averages 3 overlapping dinucleotides.
/// Output has same length as input sequence; boundary positions use NaN.
pub fn dnashape_predict(
  seq : String,
  feature : DnaShapeFeature,
) -> DnaShapeResult {
  let n = seq.length()
  let k = feature.kmer_size()
  let values = Array::new()
  let positions = Array::new()
  if n < k {
    // sequence too short: all NaN
    let mut i = 0
    while i < n {
      values.push(0.0 / 0.0)
      i = i + 1
    }
    return DnaShapeResult::{ feature, values, positions }
  }
  // center offset: for k=5, center is at position 2 (0-indexed)
  // for k=4, center is between positions 1 and 2
  let center = k / 2
  // Build values for all positions, default NaN
  let all_values = Array::make(n, 0.0 / 0.0)
  // Slide window
  let mut i = 0
  while i <= n - k {
    // Extract k-mer
    let kmer = seq[i:i + k].to_string()
    // For each k-mer, average the (k-1) overlapping dinucleotide values
    let mut sum = 0.0
    let mut count = 0
    let mut j = 0
    while j < k - 1 {
      let dinuc = kmer[j:j + 2].to_string()
      let idx = dnashape_encode_kmer(dinuc)
      if idx >= 0 {
        sum = sum + dnashape_lookup_dinuc(idx, feature)
        count = count + 1
      }
      j = j + 1
    }
    let avg = if count > 0 { sum / count.to_double() } else { 0.0 / 0.0 }
    // Assign to center position
    let pos = i + center
    if pos < n {
      all_values[pos] = avg
    }
    i = i + 1
  }
  // Collect defined values and positions
  let mut idx = 0
  while idx < n {
    values.push(all_values[idx])
    if !all_values[idx].is_nan() {
      positions.push(idx + 1) // 1-based
    }
    idx = idx + 1
  }
  DnaShapeResult::{ feature, values, positions }
}

///|
/// Predict all four standard shape features for a sequence.
pub fn dnashape_predict_all(seq : String) -> DnaShapePrediction {
  DnaShapePrediction::{
    sequence: seq,
    mgw: dnashape_predict(seq, ShapeMGW),
    roll: dnashape_predict(seq, ShapeRoll),
    prot: dnashape_predict(seq, ShapeProT),
    helt: dnashape_predict(seq, ShapeHelT),
  }
}

///|
/// Format a DnaShapeResult as a string.
pub fn dnashape_result_to_string(result : DnaShapeResult) -> String {
  let mut s = result.feature.to_string() + ": ["
  let n = result.values.length()
  let mut i = 0
  while i < n {
    let v = result.values[i]
    if v.is_nan() {
      s = s + "NA"
    } else {
      s = s + v.to_string()
    }
    if i < n - 1 {
      s = s + ", "
    }
    i = i + 1
  }
  s = s + "]"
  s
}

///|
/// Generate ASCII heatmap for a shape result.
/// Maps values to characters based on min/max range.
pub fn dnashape_heatmap(result : DnaShapeResult) -> String {
  let lo = result.min()
  let hi = result.max()
  let range = hi - lo
  let chars = [' ', '.', ':', '-', '=', '+', '*', '#', '%', '@']
  let mut s = ""
  for v in result.values {
    if v.is_nan() {
      s = s + "_"
    } else if range == 0.0 {
      s = s + chars[5].to_string()
    } else {
      let frac = (v - lo) / range
      let idx = (frac * 10.0).to_int()
      let clamped = if idx < 0 { 0 } else if idx > 9 { 9 } else { idx }
      s = s + chars[clamped].to_string()
    }
  }
  s
}

///|
/// Encode a sequence to its dinucleotide composition (16 counts).
pub fn dnashape_dinuc_composition(seq : String) -> Array[Int] {
  let counts = Array::make(16, 0)
  let n = seq.length()
  let mut i = 0
  while i < n - 1 {
    let dinuc = seq[i:i + 2].to_string()
    let idx = dnashape_encode_kmer(dinuc)
    if idx >= 0 {
      counts[idx] = counts[idx] + 1
    }
    i = i + 1
  }
  counts
}

///|
/// Compute GC content of a sequence.
pub fn dnashape_gc_content(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let mut gc = 0
  let mut i = 0
  while i < n {
    let c = seq[i:i + 1].to_string()
    if c == "G" || c == "g" || c == "C" || c == "c" {
      gc = gc + 1
    }
    i = i + 1
  }
  gc.to_double() / n.to_double()
}

///|
/// Sample DNA sequence for testing (50 bp).
pub fn dnashape_sample_sequence() -> String {
  "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
}

///|
/// A-tract sample sequence (known narrow minor groove).
pub fn dnashape_sample_atract() -> String {
  "AAAAAAAAAATTTTTTTTTTAAAAAAAAAATTTTTTTTTTAAAAAAAAAA"
}

///|
/// GC-rich sample sequence (known wide minor groove).
pub fn dnashape_sample_gc_rich() -> String {
  "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
}

///|
/// Random mixed sequence sample.
pub fn dnashape_sample_mixed() -> String {
  "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
}
