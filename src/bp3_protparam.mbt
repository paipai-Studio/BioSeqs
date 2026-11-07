///|
/// bp3_protparam.mbt — Biopython round 3: ProtParam + seq1/seq3 + ProteinData
/// #8 seq1/seq3 — 1-letter ↔ 3-letter amino acid conversion
/// #9 ProtParam — unified protein parameter calculator
/// #10 ProteinData — amino acid property data tables

///|
/// #8 seq1/seq3 conversion
let bp3_aa_1to3 : Map[String, String] = Map([
  ("A", "ALA"),
  ("R", "ARG"),
  ("N", "ASN"),
  ("D", "ASP"),
  ("C", "CYS"),
  ("E", "GLU"),
  ("Q", "GLN"),
  ("G", "GLY"),
  ("H", "HIS"),
  ("I", "ILE"),
  ("L", "LEU"),
  ("K", "LYS"),
  ("M", "MET"),
  ("F", "PHE"),
  ("P", "PRO"),
  ("S", "SER"),
  ("T", "THR"),
  ("W", "TRP"),
  ("Y", "TYR"),
  ("V", "VAL"),
])

///|
let bp3_aa_3to1 : Map[String, String] = Map([
  ("ALA", "A"),
  ("ARG", "R"),
  ("ASN", "N"),
  ("ASP", "D"),
  ("CYS", "C"),
  ("GLU", "E"),
  ("GLN", "Q"),
  ("GLY", "G"),
  ("HIS", "H"),
  ("ILE", "I"),
  ("LEU", "L"),
  ("LYS", "K"),
  ("MET", "M"),
  ("PHE", "F"),
  ("PRO", "P"),
  ("SER", "S"),
  ("THR", "T"),
  ("TRP", "W"),
  ("TYR", "Y"),
  ("VAL", "V"),
])

///|
pub fn bp3_seq1(seq3 : String) -> String {
  let r = bp3_aa_3to1.get(seq3.to_upper())
  match r {
    Some(v) => v
    None => "X"
  }
}

///|
pub fn bp3_seq3(seq1 : String) -> String {
  let r = bp3_aa_1to3.get(seq1.to_upper())
  match r {
    Some(v) => v
    None => "XXX"
  }
}

///|
pub fn bp3_seq1_to_string(three_letter_seq : String) -> String {
  let sb = StringBuilder()
  let n = three_letter_seq.length()
  let mut i = 0
  while i + 3 <= n {
    let codon = three_letter_seq[i:i + 3].to_owned()
    let _ = sb.write(bp3_seq1(codon))
    i = i + 3
  }
  sb.to_string()
}

///|
pub fn bp3_seq3_to_string(one_letter_seq : String) -> String {
  let sb = StringBuilder()
  for i in 0..<one_letter_seq.length() {
    let ch = one_letter_seq[i:i + 1].to_owned()
    let _ = sb.write(bp3_seq3(ch))
  }
  sb.to_string()
}

///|
/// #10 ProteinData — amino acid property data tables
let bp3_hydropathy : Map[String, Double] = Map([
  ("I", 4.5),
  ("V", 4.2),
  ("L", 3.8),
  ("F", 2.8),
  ("C", 2.5),
  ("M", 1.9),
  ("A", 1.8),
  ("G", -0.4),
  ("T", -0.7),
  ("S", -0.8),
  ("W", -0.9),
  ("Y", -1.3),
  ("P", -1.6),
  ("H", -3.2),
  ("E", -3.5),
  ("Q", -3.5),
  ("D", -3.5),
  ("N", -3.5),
  ("K", -3.9),
  ("R", -4.5),
])

///|
let bp3_helix_propensity : Map[String, Double] = Map([
  ("A", 1.42),
  ("R", 0.98),
  ("N", 0.67),
  ("D", 1.01),
  ("C", 0.70),
  ("E", 1.51),
  ("Q", 1.11),
  ("G", 0.57),
  ("H", 1.00),
  ("I", 1.08),
  ("L", 1.21),
  ("K", 1.16),
  ("M", 1.45),
  ("F", 1.13),
  ("P", 0.57),
  ("S", 0.77),
  ("T", 0.83),
  ("W", 1.08),
  ("Y", 0.69),
  ("V", 1.06),
])

///|
let bp3_sheet_propensity : Map[String, Double] = Map([
  ("A", 0.83),
  ("R", 0.93),
  ("N", 0.89),
  ("D", 0.54),
  ("C", 1.19),
  ("E", 0.37),
  ("Q", 1.10),
  ("G", 0.75),
  ("H", 0.87),
  ("I", 1.60),
  ("L", 1.30),
  ("K", 0.74),
  ("M", 1.05),
  ("F", 1.38),
  ("P", 0.55),
  ("S", 0.75),
  ("T", 1.19),
  ("W", 1.37),
  ("Y", 1.47),
  ("V", 1.70),
])

///|
let bp3_turn_propensity : Map[String, Double] = Map([
  ("A", 0.66),
  ("R", 0.95),
  ("N", 1.56),
  ("D", 1.46),
  ("C", 0.54),
  ("E", 0.74),
  ("Q", 0.98),
  ("G", 1.56),
  ("H", 0.95),
  ("I", 0.47),
  ("L", 0.59),
  ("K", 1.01),
  ("M", 0.60),
  ("F", 0.60),
  ("P", 1.52),
  ("S", 1.43),
  ("T", 0.96),
  ("W", 0.96),
  ("Y", 1.14),
  ("V", 0.50),
])

///|
let bp3_flexibility : Map[String, Double] = Map([
  ("A", 0.357),
  ("R", 0.529),
  ("N", 0.463),
  ("D", 0.511),
  ("C", 0.346),
  ("E", 0.502),
  ("Q", 0.493),
  ("G", 0.544),
  ("H", 0.323),
  ("I", 0.376),
  ("L", 0.368),
  ("K", 0.529),
  ("M", 0.368),
  ("F", 0.280),
  ("P", 0.509),
  ("S", 0.493),
  ("T", 0.442),
  ("W", 0.278),
  ("Y", 0.294),
  ("V", 0.395),
])

///|
pub fn bp3_aa_hydropathy(aa : String) -> Double {
  bp3_hydropathy.get(aa).unwrap_or(0.0)
}

///|
pub fn bp3_aa_helix_propensity(aa : String) -> Double {
  bp3_helix_propensity.get(aa).unwrap_or(1.0)
}

///|
pub fn bp3_aa_sheet_propensity(aa : String) -> Double {
  bp3_sheet_propensity.get(aa).unwrap_or(1.0)
}

///|
pub fn bp3_aa_turn_propensity(aa : String) -> Double {
  bp3_turn_propensity.get(aa).unwrap_or(1.0)
}

///|
pub fn bp3_aa_flexibility(aa : String) -> Double {
  bp3_flexibility.get(aa).unwrap_or(0.4)
}

///|
/// #9 ProtParam — unified protein parameter calculator
pub struct Bp3ProtParam {
  sequence : String
  mw : Double
  gravy : Double
  aromaticity : Double
  instability_index : Double
  pI : Double
  extinction_coeff : Double
  helix_fraction : Double
  sheet_fraction : Double
  turn_fraction : Double
}

///|
let bp3_aa_mw : Map[String, Double] = Map([
  ("A", 89.09),
  ("R", 174.20),
  ("N", 132.12),
  ("D", 133.10),
  ("C", 121.16),
  ("E", 147.13),
  ("Q", 146.15),
  ("G", 75.07),
  ("H", 155.16),
  ("I", 131.17),
  ("L", 131.17),
  ("K", 146.19),
  ("M", 149.21),
  ("F", 165.19),
  ("P", 115.13),
  ("S", 105.09),
  ("T", 119.12),
  ("W", 204.23),
  ("Y", 181.19),
  ("V", 117.15),
])

///|
pub fn bp3_protparam(sequence : String) -> Bp3ProtParam {
  let n = sequence.length()
  let mut mw = 18.02
  let mut gravy_sum = 0.0
  let mut aromatic = 0
  let mut instability = 0.0
  let mut helix_sum = 0.0
  let mut sheet_sum = 0.0
  let mut turn_sum = 0.0
  let mw_map = bp3_aa_mw
  let diw : Array[Double] = []
  for i in 0..<n {
    let aa = sequence[i:i + 1].to_owned()
    let m = mw_map.get(aa).unwrap_or(110.0)
    mw = mw + m
    gravy_sum = gravy_sum + bp3_aa_hydropathy(aa)
    if aa == "W" || aa == "Y" || aa == "F" {
      aromatic = aromatic + 1
    }
    helix_sum = helix_sum + bp3_aa_helix_propensity(aa)
    sheet_sum = sheet_sum + bp3_aa_sheet_propensity(aa)
    turn_sum = turn_sum + bp3_aa_turn_propensity(aa)
    diw.push(bp3_aa_hydropathy(aa))
  }
  let _nm1 = n - 1
  for i in 0..<_nm1 {
    let aa1 = sequence[i:i + 1].to_owned()
    let aa2 = sequence[i + 1:i + 2].to_owned()
    instability = instability + diw[i] * diw[i + 1]
    ignore(aa1)
    ignore(aa2)
  }
  let gravy = if n > 0 { gravy_sum / n.to_double() } else { 0.0 }
  let aromaticity = if n > 0 {
    aromatic.to_double() / n.to_double()
  } else {
    0.0
  }
  let instability_index = if n > 1 {
    10.0 / (n.to_double() - 1.0) * instability
  } else {
    0.0
  }
  let helix_fraction = if n > 0 { helix_sum / n.to_double() / 3.0 } else { 0.0 }
  let sheet_fraction = if n > 0 { sheet_sum / n.to_double() / 3.0 } else { 0.0 }
  let turn_fraction = if n > 0 { turn_sum / n.to_double() / 3.0 } else { 0.0 }
  let mut n_trp = 0
  let mut n_tyr = 0
  let mut n_cys = 0
  for i in 0..<n {
    let aa = sequence[i:i + 1].to_owned()
    if aa == "W" {
      n_trp = n_trp + 1
    }
    if aa == "Y" {
      n_tyr = n_tyr + 1
    }
    if aa == "C" {
      n_cys = n_cys + 1
    }
  }
  let extinction_coeff = n_trp.to_double() * 5500.0 +
    n_tyr.to_double() * 1490.0 +
    n_cys.to_double() * 125.0
  let pI = bp3_protparam_pi(sequence)
  {
    sequence,
    mw,
    gravy,
    aromaticity,
    instability_index,
    pI,
    extinction_coeff,
    helix_fraction,
    sheet_fraction,
    turn_fraction,
  }
}

///|
fn bp3_protparam_pi(sequence : String) -> Double {
  let mut n_pos = 0.0
  let mut n_neg = 0.0
  for i in 0..<sequence.length() {
    let aa = sequence[i:i + 1].to_owned()
    match aa {
      "K" => n_pos = n_pos + 1.0
      "R" => n_pos = n_pos + 1.0
      "H" => n_pos = n_pos + 1.0
      "D" => n_neg = n_neg + 1.0
      "E" => n_neg = n_neg + 1.0
      _ => ()
    }
  }
  if n_pos == 0.0 && n_neg == 0.0 {
    return 7.0
  }
  if n_pos == 0.0 {
    return 3.0
  }
  if n_neg == 0.0 {
    return 11.0
  }
  let mut best_ph = 7.0
  let mut best_charge = 999999.0
  let mut ph = 0.0
  while ph <= 14.0 {
    let charge = n_pos / (1.0 + @math.pow(10.0, ph - 10.5)) -
      n_neg / (1.0 + @math.pow(10.0, 3.65 - ph))
    let abs_charge = if charge < 0.0 { -charge } else { charge }
    if abs_charge < best_charge {
      best_charge = abs_charge
      best_ph = ph
    }
    ph = ph + 0.1
  }
  best_ph
}

///|
pub fn bp3_protparam_charge_at_ph(sequence : String, ph : Double) -> Double {
  let mut charge = 0.0
  for i in 0..<sequence.length() {
    let aa = sequence[i:i + 1].to_owned()
    match aa {
      "K" => charge = charge + 1.0 / (1.0 + @math.pow(10.0, ph - 10.5))
      "R" => charge = charge + 1.0 / (1.0 + @math.pow(10.0, ph - 12.5))
      "H" => charge = charge + 1.0 / (1.0 + @math.pow(10.0, ph - 6.0))
      "D" => charge = charge - 1.0 / (1.0 + @math.pow(10.0, 3.65 - ph))
      "E" => charge = charge - 1.0 / (1.0 + @math.pow(10.0, 4.25 - ph))
      _ => ()
    }
  }
  charge
}

///|
pub fn _bp3_protparam_dummy() -> Unit {
  let _ = Bp3ProtParam::{
    sequence: "",
    mw: 0.0,
    gravy: 0.0,
    aromaticity: 0.0,
    instability_index: 0.0,
    pI: 0.0,
    extinction_coeff: 0.0,
    helix_fraction: 0.0,
    sheet_fraction: 0.0,
    turn_fraction: 0.0,
  }
}
