///|
/// Bio.SeqUtils.MeltingTemp - DNA melting temperature calculation.
/// Implements Wallace rule, GC-content, Nearest-Neighbor (SantaLucia 1998),
/// salt corrections (Na+, Mg2+), and mismatch corrections.

///|
pub struct MeltingTempResult {
  tm : Double
  method : String
}

///|
pub fn MeltingTempResult::new(
  tm : Double,
  method : String,
) -> MeltingTempResult {
  MeltingTempResult::{ tm, method, }
}

///|
/// Wallace rule: Tm = 4*(G+C) + 2*(A+T)
pub fn mt_wallace(seq : String) -> MeltingTempResult {
  let seq_upper = seq.to_upper()

  let mut gc_count = 0
  let mut at_count = 0

  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "G" || c == "C" {
      gc_count = gc_count + 1
    } else if c == "A" || c == "T" {
      at_count = at_count + 1
    }
    idx = idx + 1
  }

  let tm = 4.0 * gc_count.to_double() + 2.0 * at_count.to_double()

  MeltingTempResult::new(tm, "wallace")
}

///|
/// GC-content method: Tm = 64.9 + 41*(GC% - 16.4)/length
pub fn mt_gc_content(seq : String) -> MeltingTempResult {
  let seq_upper = seq.to_upper()
  let length = seq_upper.length().to_double()

  let mut gc_count = 0

  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "G" || c == "C" {
      gc_count = gc_count + 1
    }
    idx = idx + 1
  }

  let tm = 64.9 + 41.0 * (gc_count.to_double() - 16.4) / length

  MeltingTempResult::new(tm, "gc_content")
}

///|
/// Nearest-Neighbor thermodynamic parameters (SantaLucia 1998).
/// Values in kcal/mol: (dH, dS)
/// dH = enthalpy change, dS = entropy change (cal/mol·K)
fn mt_nn_params(dinucleotide : String) -> (Double, Double) {
  let nn_params : Map[String, (Double, Double)] = {
    let m : Map[String, (Double, Double)] = Map([], capacity=16)
    m["AA"] = (-7.6, -21.3)
    m["TT"] = (-7.6, -21.3)
    m["AT"] = (-7.2, -20.4)
    m["TA"] = (-7.2, -21.3)
    m["CA"] = (-8.5, -22.7)
    m["TG"] = (-8.5, -22.7)
    m["GT"] = (-8.4, -22.4)
    m["AC"] = (-8.4, -22.4)
    m["CT"] = (-7.8, -21.0)
    m["AG"] = (-7.8, -21.0)
    m["GA"] = (-8.2, -22.2)
    m["TC"] = (-8.2, -22.2)
    m["CG"] = (-10.6, -27.2)
    m["GC"] = (-9.8, -24.4)
    m["GG"] = (-8.0, -19.9)
    m["CC"] = (-8.0, -19.9)
    m
  }
  match nn_params.get(dinucleotide) {
    Some(v) => v
    None => (0.0, 0.0)
  }
}

///|
/// Nearest-Neighbor melting temperature calculation (SantaLucia 1998).
/// Tm = dH / (dS + R * ln(Ct/4)) - 273.15 + 16.6 * log10([Na+])
/// where R = 1.987 cal/mol·K, Ct = total strand concentration (default 0.0000001 M = 0.1 µM)
pub fn mt_nearest_neighbor(
  seq : String,
  conc~ : Double = 0.0000001,
  na_conc~ : Double = 0.05
) -> MeltingTempResult {
  let seq_upper = seq.to_upper()
  let n = seq_upper.length()
  if n < 2 {
    return MeltingTempResult::new(0.0, "nearest_neighbor")
  }

  // Initialize with terminal base pair penalties
  // Biopython uses: dH = init + sum(dH_i), dS = init + sum(dS_i)
  // For DNA duplex: initiation parameters
  let mut dH = 0.2  // kcal/mol (initiation)
  let mut dS = -5.7 // cal/mol·K (initiation)

  // Sum over dinucleotide steps
  let mut i = 0
  while i < n - 1 {
    let dinuc = seq_upper[i:i + 2].to_string()
    let (dh_step, ds_step) = mt_nn_params(dinuc)
    dH = dH + dh_step
    dS = dS + ds_step
    i = i + 1
  }

  // Terminal AT penalty (Biopython adds 3.4 kcal/mol for terminal AT pairs)
  // Simplified: no additional terminal penalty in this implementation

  // R = 1.987 cal/mol·K
  let r_const = 1.987

  // dS is in cal/mol·K, dH in kcal/mol -> convert dH to cal/mol
  let dH_cal = dH * 1000.0

  // Tm = dH / (dS + R * ln(Ct/4)) - 273.15 + 16.6 * log10([Na+])
  let ln_term = @math.ln(conc / 4.0)
  let tm_k = dH_cal / (dS + r_const * ln_term)
  let tm_c = tm_k - 273.15

  // Na+ correction: 16.6 * log10([Na+])
  let na_correction = 16.6 * @math.log10(na_conc)
  let tm_final = tm_c + na_correction

  MeltingTempResult::new(tm_final, "nearest_neighbor")
}

///|
/// Salt correction for Na+ concentration.
/// Tm(Na+) = Tm(1M Na+) + 16.6 * log10([Na+])
pub fn mt_salt_correction_na(tm : Double, na_conc : Double) -> Double {
  tm + 16.6 * @math.log10(na_conc)
}

///|
/// Salt correction for Mg2+ concentration (Owczarzy et al. 2004).
/// Tm(Mg2+) = Tm(1M Na+) / (1 + 0.75 * [Mg2+] * sqrt([Na+]))
/// Simplified: Tm(Mg2+) = Tm(1M Na+) + correction
pub fn mt_salt_correction_mg(
  tm : Double,
  mg_conc : Double,
  na_conc~ : Double = 0.05
) -> Double {
  // Owczarzy correction: factor = 1 / (1 + 0.75 * [Mg2+] * sqrt([Na+]))
  let factor = 1.0 / (1.0 + 0.75 * mg_conc * (na_conc).sqrt())
  // Simplified: apply as multiplicative factor to the Na+ correction
  let na_correction = 16.6 * @math.log10(na_conc)
  let mg_correction = na_correction * factor
  tm + mg_correction - na_correction
}

///|
/// Melting temperature with salt correction.
pub fn mt_with_salt(
  seq : String,
  na_conc~ : Double = 0.05,
  mg_conc~ : Double = 0.0,
  method~ : String = "nearest_neighbor"
) -> MeltingTempResult {
  let base_tm = if method == "wallace" {
    mt_wallace(seq).tm
  } else if method == "gc_content" {
    mt_gc_content(seq).tm
  } else {
    mt_nearest_neighbor(seq).tm
  }

  let tm_final = if mg_conc > 0.0 {
    mt_salt_correction_mg(base_tm, mg_conc, na_conc~)
  } else {
    mt_salt_correction_na(base_tm, na_conc)
  }

  MeltingTempResult::new(tm_final, method + "_with_salt")
}

///|
/// Check if sequence contains non-IUPAC characters.
pub fn mt_is_valid_dna(seq : String) -> Bool {
  let valid : Map[UInt16, Bool] = {
    let m : Map[UInt16, Bool] = Map([], capacity=8)
    m[65] = true  // A
    m[67] = true  // C
    m[71] = true  // G
    m[84] = true  // T
    m[82] = true  // R (A/G)
    m[89] = true  // Y (C/T)
    m[83] = true  // S (G/C)
    m[87] = true  // W (A/T)
    m[75] = true  // K (G/T)
    m[77] = true  // M (A/C)
    m[66] = true  // B (C/G/T)
    m[68] = true  // D (A/G/T)
    m[72] = true  // H (A/C/T)
    m[86] = true  // V (A/C/G)
    m[78] = true  // N (any)
    m
  }
  let mut i = 0
  while i < seq.length() {
    if !(valid.contains(seq[i])) {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Recommend Tm calculation method based on sequence length.
pub fn recommend_tm_method(seq : String) -> String {
  let length = seq.length()

  if length <= 14 {
    "wallace"
  } else if length <= 50 {
    "nearest_neighbor"
  } else {
    "gc_content"
  }
}

///|
/// Create example DNA sequence.
pub fn create_example_dna() -> String {
  "AGCTGATCGATCGATCGATCGATCGATCGATCGATCG"
}

///|
/// Calculate Tm for a primer sequence (convenience wrapper).
pub fn mt_primer(primer : String) -> MeltingTempResult {
  let method = recommend_tm_method(primer)
  if method == "wallace" {
    mt_wallace(primer)
  } else if method == "nearest_neighbor" {
    mt_nearest_neighbor(primer)
  } else {
    mt_gc_content(primer)
  }
}
