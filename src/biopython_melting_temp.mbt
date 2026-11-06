/// biopython_melting_temp.mbt — DNA melting temperature calculation
/// Inspired by Biopython Bio.SeqUtils.MeltingTemp: Tm Wallace, GC, NN

///|
/// Wallace rule: Tm = 2*(A+T) + 4*(G+C)
/// Simple rule, valid for short oligos (< 14 bp).
pub fn biopython_tm_wallace(sequence : String) -> Double {
  let mut at = 0.0
  let mut gc = 0.0
  for i = 0; i < sequence.length(); i = i + 1 {
    match sequence[i].to_int() {
      65 | 84 | 97 | 116 => at = at + 1.0 // A, T
      71 | 67 | 103 | 99 => gc = gc + 1.0 // G, C
      _ => ()
    }
  }
  2.0 * at + 4.0 * gc
}

///|
/// GC-content based Tm.
/// Tm = 64.9 + 41 * (G+C - 16.4) / (A+T+G+C)
pub fn biopython_tm_gc(
  sequence : String,
  na_conc~ : Double = 50.0, // mM
) -> Double {
  let mut gc = 0.0
  let mut n = 0.0
  for i = 0; i < sequence.length(); i = i + 1 {
    match sequence[i].to_int() {
      71 | 67 | 103 | 99 => {
        gc = gc + 1.0
        n = n + 1.0
      }
      65 | 84 | 97 | 116 => n = n + 1.0
      _ => ()
    }
  }
  if n == 0.0 {
    return 0.0
  }
  let tm = 64.9 + 41.0 * (gc - 16.4) / n
  // Salt correction
  let salt_corr = 16.6 * @math.log10(1.0 + na_conc / 1000.0)
  tm + salt_corr
}

///|
/// Nearest-neighbor thermodynamics (simplified).
/// Uses unified NN parameters (SantaLucia 1998).
pub fn biopython_tm_nearest_neighbor(
  sequence : String,
  na_conc~ : Double = 50.0, // mM
  primer_conc~ : Double = 250.0, // nM
) -> Double {
  let n = sequence.length()
  if n < 2 {
    return 0.0
  }
  let seq = sequence.to_upper()
  // NN parameters: (dH kcal/mol, dS cal/mol·K) for each dinucleotide
  let (dH, dS) = biopython_nn_thermo(seq)
  // Salt correction
  let na = na_conc / 1000.0 // M
  let ct = primer_conc / 1000000000.0 // M
  // Tm = dH*1000 / (dS + R * ln(Ct/x)) - salt correction
  let r = 1.987 // cal/(mol·K)
  let x = 4.0 // for non-self-complementary
  let tm = dH * 1000.0 / (dS + r * @math.ln(ct / x))
  // Salt correction (SantaLucia 1998)
  let salt_corr = 0.368 * (n - 1).to_double() * @math.ln(na)
  tm + salt_corr
}

///|
/// Get nearest-neighbor thermodynamics for a sequence.
/// Returns (total_dH, total_dS).
fn biopython_nn_thermo(seq : String) -> (Double, Double) {
  // SantaLucia (1998) unified parameters
  // dH in kcal/mol, dS in cal/(mol·K)
  let mut total_dH = 0.0
  let mut total_dS = 0.0
  for i = 0; i < seq.length() - 1; i = i + 1 {
    let nn = seq[i:i + 2].to_owned()
    let (dh, ds) = match nn {
      "AA" | "TT" => (-7.9, -22.2)
      "AT" => (-7.2, -20.4)
      "TA" => (-7.2, -21.3)
      "CA" | "TG" => (-8.5, -22.7)
      "GT" | "AC" => (-8.4, -22.4)
      "CT" | "AG" => (-7.8, -21.0)
      "GA" | "TC" => (-8.2, -22.2)
      "CG" => (-10.6, -27.2)
      "GC" => (-9.8, -24.4)
      "GG" | "CC" => (-8.0, -19.9)
      _ => (0.0, 0.0)
    }
    total_dH = total_dH + dh
    total_dS = total_dS + ds
  }
  // Initiation parameters
  // Terminal AT: dH = +2.3, dS = +4.1
  // Terminal GC: dH = +0.1, dS = -2.8
  let first = seq[0:1].to_owned()
  let last = seq[seq.length() - 1:].to_owned()
  for end in [first, last] {
    if end == "A" || end == "T" {
      total_dH = total_dH + 2.3
      total_dS = total_dS + 4.1
    } else {
      total_dH = total_dH + 0.1
      total_dS = total_dS + -2.8
    }
  }
  (total_dH, total_dS)
}

///|
/// Chemically modified bases Tm correction (simplified).
pub fn biopython_tm_chem_correction(
  tm : Double,
  dms_conc~ : Double = 0.0, // DMSO %
  formamide_conc~ : Double = 0.0, // formamide %
) -> Double {
  // DMSO: -0.6°C per 1% DMSO
  // Formamide: -0.5°C per 1% formamide
  tm - 0.6 * dms_conc - 0.5 * formamide_conc
}
