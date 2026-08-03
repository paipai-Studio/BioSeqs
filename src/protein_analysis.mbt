// ============================================================================
// Protein Analysis Module
// ============================================================================
// Comprehensive protein sequence analysis tools for:
//   - Conservation scoring
//   - Hydrophobicity calculation (Kyte-Doolittle)
//   - Secondary structure prediction (GOR method)
//   - Antigenicity prediction
//   - Transmembrane domain prediction
//   - Dipeptide/tripeptide composition
// ============================================================================

// ============================================================================
// Kyte-Doolittle Hydrophobicity Scale
// ============================================================================

fn protein_get_kd(c : Char) -> Double {
  match c {
    'A' => 1.8
    'R' => -4.5
    'N' => -3.5
    'D' => -3.5
    'C' => 2.5
    'Q' => -3.5
    'E' => -3.5
    'G' => -0.4
    'H' => -3.2
    'I' => 4.5
    'L' => 3.8
    'K' => -3.9
    'M' => 1.9
    'F' => 2.8
    'P' => -1.6
    'S' => -0.8
    'T' => -0.7
    'W' => -0.9
    'Y' => -1.3
    'V' => 4.2
    _ => 0.0
  }
}

// ============================================================================
// GOR Secondary Structure Parameters
// ============================================================================
// GOR method (Garnier-Osguthorpe-Robson) for secondary structure prediction.
// Each amino acid has parameters for:
//   - P_alpha: probability of being in alpha-helix
//   - P_beta: probability of being in beta-sheet
//   - P_turn: probability of being in turn/coil
// ============================================================================

struct GorParams {
  pa : Double
  pb : Double
  pt : Double
}

fn protein_get_gor(c : Char) -> GorParams {
  match c {
    'A' => GorParams::{ pa: 1.42, pb: 0.83, pt: 0.66 }
    'R' => GorParams::{ pa: 0.98, pb: 0.93, pt: 0.95 }
    'N' => GorParams::{ pa: 0.67, pb: 0.89, pt: 1.56 }
    'D' => GorParams::{ pa: 1.01, pb: 0.54, pt: 1.46 }
    'C' => GorParams::{ pa: 0.70, pb: 1.19, pt: 1.19 }
    'Q' => GorParams::{ pa: 1.11, pb: 1.10, pt: 0.98 }
    'E' => GorParams::{ pa: 1.51, pb: 0.37, pt: 0.74 }
    'G' => GorParams::{ pa: 0.57, pb: 0.75, pt: 1.56 }
    'H' => GorParams::{ pa: 1.00, pb: 0.87, pt: 0.95 }
    'I' => GorParams::{ pa: 1.08, pb: 1.60, pt: 0.47 }
    'L' => GorParams::{ pa: 1.21, pb: 1.30, pt: 0.59 }
    'K' => GorParams::{ pa: 1.16, pb: 0.74, pt: 1.01 }
    'M' => GorParams::{ pa: 1.45, pb: 1.05, pt: 0.60 }
    'F' => GorParams::{ pa: 1.13, pb: 1.38, pt: 0.60 }
    'P' => GorParams::{ pa: 0.57, pb: 0.55, pt: 1.52 }
    'S' => GorParams::{ pa: 0.77, pb: 0.75, pt: 1.43 }
    'T' => GorParams::{ pa: 0.83, pb: 1.19, pt: 0.96 }
    'W' => GorParams::{ pa: 1.08, pb: 1.37, pt: 0.96 }
    'Y' => GorParams::{ pa: 0.69, pb: 1.47, pt: 1.14 }
    'V' => GorParams::{ pa: 1.06, pb: 1.70, pt: 0.50 }
    _ => GorParams::{ pa: 1.00, pb: 1.00, pt: 1.00 }
  }
}

// ============================================================================
// TM Scale (transmembrane propensity)
// ============================================================================

struct TmParams {
  inside : Double
  outside : Double
  tm : Double
}

fn protein_get_tm(c : Char) -> TmParams {
  match c {
    'A' => TmParams::{ inside: 1.1, outside: 0.56, tm: 1.32 }
    'R' => TmParams::{ inside: 0.0, outside: 1.0, tm: 0.0 }
    'N' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.6 }
    'D' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.7 }
    'C' => TmParams::{ inside: 1.1, outside: 0.56, tm: 0.62 }
    'Q' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.2 }
    'E' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.3 }
    'G' => TmParams::{ inside: 0.5, outside: 0.2, tm: 0.48 }
    'H' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.5 }
    'I' => TmParams::{ inside: 1.0, outside: 0.32, tm: 1.0 }
    'L' => TmParams::{ inside: 1.0, outside: 0.36, tm: 1.24 }
    'K' => TmParams::{ inside: 0.1, outside: 1.0, tm: 0.0 }
    'M' => TmParams::{ inside: 1.0, outside: 0.36, tm: 1.11 }
    'F' => TmParams::{ inside: 1.0, outside: 0.42, tm: 1.02 }
    'P' => TmParams::{ inside: 0.5, outside: 0.2, tm: -0.3 }
    'S' => TmParams::{ inside: 0.5, outside: 0.2, tm: 0.38 }
    'T' => TmParams::{ inside: 0.5, outside: 0.2, tm: 0.38 }
    'W' => TmParams::{ inside: 1.0, outside: 0.44, tm: 1.0 }
    'Y' => TmParams::{ inside: 0.5, outside: 0.2, tm: 0.69 }
    'V' => TmParams::{ inside: 1.0, outside: 0.32, tm: 1.19 }
    _ => TmParams::{ inside: 0.5, outside: 0.2, tm: 0.0 }
  }
}

// ============================================================================
// Shannon Entropy for Conservation
// ============================================================================

fn protein_shannon_entropy(freqs : Array[Double]) -> Double {
  let mut h = 0.0
  let mut i = 0
  while i < freqs.length() {
    let f = freqs[i]
    if f > 0.0 {
      h = h - f * @math.log2(f)
    }
    i = i + 1
  }
  h
}

// ============================================================================
// Amino Acid Frequency at Position
// ============================================================================

fn protein_aa_frequency_at(alignment : Array[String], pos : Int) -> Array[Double] {
  let freqs : Array[Double] = Array::new()
  let n = alignment.length().to_double()
  if n == 0.0 { return freqs }
  let aa_map = Map([], capacity=26)
  let mut i = 0
  while i < alignment.length() {
    let seq = alignment[i]
    if pos < seq.length() {
      let char_code = seq.unsafe_get(pos).to_int()
      let c = if char_code >= 97 && char_code <= 122 {
        Char::from_int(char_code - 32)
      } else {
        Char::from_int(char_code)
      }
      let count = match aa_map.get(c) {
        Some(v) => v + 1.0
        None => 1.0
      }
      aa_map[c] = count
    }
    i = i + 1
  }
  let keys = aa_map.keys().collect()
  let mut j = 0
  while j < keys.length() {
    let key = keys[j]
    let freq = aa_map[key] / n
    freqs.push(freq)
    j = j + 1
  }
  freqs
}

// ============================================================================
// Public API
// ============================================================================

///|
/// Calculate per-position conservation score from a multiple sequence alignment.
/// Returns scores from 0 (no conservation) to 1 (perfect conservation).
pub fn protein_conservation(alignment : Array[String]) -> Array[Double] {
  if alignment.length() == 0 { return [] }
  let aln_len = alignment[0].length()
  let alphabet_size = 20.0
  let conservation : Array[Double] = Array::new()
  let mut pos = 0
  while pos < aln_len {
    let aa_freq = protein_aa_frequency_at(alignment, pos)
    let entropy = protein_shannon_entropy(aa_freq)
    let score = 1.0 - entropy / @math.log2(alphabet_size)
    conservation.push(score)
    pos = pos + 1
  }
  conservation
}

///|
/// Calculate Kyte-Doolittle hydropathy values for a protein sequence.
/// Returns array of hydropathy values averaged over a sliding window.
pub fn protein_hydropathy(sequence : String, window? : Int = 7) -> Array[Double] {
  let seq = sequence.to_upper()
  let w = if window <= 0 { 7 } else { window }
  let result : Array[Double] = Array::new()
  if seq.length() < w { return result }
  let mut i = 0
  while i <= seq.length() - w {
    let mut sum = 0.0
    let mut j = 0
    while j < w {
      let c = seq.unsafe_get(i + j).unsafe_to_char()
      sum = sum + protein_get_kd(c)
      j = j + 1
    }
    result.push(sum / w.to_double())
    i = i + 1
  }
  result
}

///|
/// Predict secondary structure using the GOR method.
/// Returns array of ('H' for helix, 'E' for sheet, 'T' for turn, 'C' for coil).
pub fn protein_secondary_structure(sequence : String) -> Array[String] {
  let seq = sequence.to_upper()
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let params = protein_get_gor(c)
    let ss = if params.pa >= params.pb && params.pa >= params.pt {
      "H"
    } else if params.pb >= params.pt {
      "E"
    } else if params.pt > 1.0 {
      "T"
    } else {
      "C"
    }
    result.push(ss)
    i = i + 1
  }
  result
}

///|
/// Calculate antigenicity scale (Hopp-Woods) for each residue.
/// Positive values indicate potential antigenic sites.
pub fn protein_antigenicity(sequence : String) -> Array[Double] {
  let seq = sequence.to_upper()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let val = match c {
      'A' => 0.06
      'R' => 0.37
      'N' => 0.07
      'D' => 0.16
      'C' => 0.08
      'Q' => 0.03
      'E' => 0.03
      'G' => 0.16
      'H' => -0.01
      'I' => -0.06
      'L' => -0.02
      'K' => 0.16
      'M' => -0.03
      'F' => -0.05
      'P' => 0.18
      'S' => 0.03
      'T' => -0.01
      'W' => -0.02
      'Y' => -0.06
      'V' => -0.07
      _ => 0.0
    }
    result.push(val)
    i = i + 1
  }
  result
}

///|
/// Predict transmembrane domains using the TM scale.
/// Returns array of scores: positive = transmembrane, negative = inside/outside.
pub fn protein_transmembrane(sequence : String, window? : Int = 19, threshold? : Double = 1.6) -> Array[Double] {
  let seq = sequence.to_upper()
  let w = if window <= 0 { 19 } else { window }
  let result : Array[Double] = Array::new()
  if seq.length() < w { return result }
  let mut i = 0
  while i <= seq.length() - w {
    let mut sum = 0.0
    let mut j = 0
    while j < w {
      let c = seq.unsafe_get(i + j).unsafe_to_char()
      let tm = protein_get_tm(c)
      sum = sum + tm.tm
      j = j + 1
    }
    let score = sum / w.to_double()
    result.push(score)
    i = i + 1
  }
  result
}

///|
/// Get the transmembrane propensity score for a specific amino acid.
pub fn protein_transmembrane_propensity(c : Char) -> Double {
  let tm = protein_get_tm(c)
  tm.tm
}

///|
/// Calculate amino acid composition (frequency of each amino acid).
pub fn protein_aa_composition(sequence : String) -> Map[Char, Double] {
  let seq = sequence.to_upper()
  let n = seq.length().to_double()
  let comp = Map([], capacity=20)
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let frac = 1.0 / n
    match comp.get(c) {
      Some(v) => comp[c] = v + frac
      None => comp[c] = frac
    }
    i = i + 1
  }
  comp
}

///|
/// Calculate dipeptide frequency composition.
pub fn protein_dipeptide_composition(sequence : String) -> Map[String, Double] {
  let seq = sequence.to_upper()
  let n = (seq.length() - 1).to_double()
  let comp = Map([], capacity=400)
  if seq.length() < 2 { return comp }
  let mut i = 0
  while i < seq.length() - 1 {
    let di = seq[i:i + 2].to_owned()
    let frac = 1.0 / n
    match comp.get(di) {
      Some(v) => comp[di] = v + frac
      None => comp[di] = frac
    }
    i = i + 1
  }
  comp
}

///|
/// Calculate tripeptide frequency composition.
pub fn protein_tripeptide_composition(sequence : String) -> Map[String, Double] {
  let seq = sequence.to_upper()
  let n = (seq.length() - 2).to_double()
  let comp = Map([], capacity=8000)
  if seq.length() < 3 { return comp }
  let mut i = 0
  while i < seq.length() - 2 {
    let tri = seq[i:i + 3].to_owned()
    let frac = 1.0 / n
    match comp.get(tri) {
      Some(v) => comp[tri] = v + frac
      None => comp[tri] = frac
    }
    i = i + 1
  }
  comp
}

///|
/// Generate a simple summary report for a protein sequence.
pub fn protein_summary(sequence : String) -> String {
  let seq = sequence.to_upper()
  let mut out = "Protein Sequence Analysis Report\n"
  out = out + "================================\n"
  out = out + "Length: " + seq.length().to_string() + " residues\n"
  let comp = protein_aa_composition(seq)
  let keys = comp.keys().collect()
  let sorted_keys : Array[Char] = Array::new()
  let mut ki = 0
  while ki < keys.length() {
    sorted_keys.push(keys[ki])
    ki = ki + 1
  }
  sorted_keys.sort_by(fn(a : Char, b : Char) -> Int {
    a.to_int() - b.to_int()
  })
  out = out + "\nAmino Acid Composition:\n"
  let mut si = 0
  while si < sorted_keys.length() {
    let aa = sorted_keys[si]
    let freq = comp[aa] * 100.0
    out = out + "  " + aa.to_string() + ": " + freq.to_string() + "%\n"
    si = si + 1
  }
  let hydro = protein_hydropathy(seq, window=7)
  let mut total_hydro = 0.0
  let mut hi = 0
  while hi < hydro.length() {
    total_hydro = total_hydro + hydro[hi]
    hi = hi + 1
  }
  let avg_hydro = if hydro.length() > 0 { total_hydro / hydro.length().to_double() } else { 0.0 }
  out = out + "\nAverage Hydropathy (window=7): " + avg_hydro.to_string() + "\n"
  let tm = protein_transmembrane(seq, window=19)
  let mut has_tm = false
  let mut ti = 0
  while ti < tm.length() {
    if tm[ti] > 1.6 { has_tm = true }
    ti = ti + 1
  }
  out = out + "Transmembrane domain detected: " + (if has_tm { "Yes" } else { "No" }) + "\n"
  let ss = protein_secondary_structure(seq)
  let mut h_count = 0
  let mut e_count = 0
  let mut t_count = 0
  let mut c_count = 0
  let mut ssi = 0
  while ssi < ss.length() {
    match ss[ssi] {
      "H" => h_count = h_count + 1
      "E" => e_count = e_count + 1
      "T" => t_count = t_count + 1
      _ => c_count = c_count + 1
    }
    ssi = ssi + 1
  }
  let total = seq.length().to_double()
  out = out + "\nSecondary Structure Prediction (GOR):\n"
  out = out + "  Helix (H): " + (h_count.to_double() / total * 100.0).to_string() + "%\n"
  out = out + "  Sheet (E): " + (e_count.to_double() / total * 100.0).to_string() + "%\n"
  out = out + "  Turn (T): " + (t_count.to_double() / total * 100.0).to_string() + "%\n"
  out = out + "  Coil (C): " + (c_count.to_double() / total * 100.0).to_string() + "%\n"
  out
}
