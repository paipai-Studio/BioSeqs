///|
/// ProtParam - Protein parameter analysis module inspired by Biopython Bio.SeqUtils.ProtParam.
/// Extends the existing ProteinAnalysis type with additional methods.

///|
pub fn ProteinAnalysis::instability_index(self : ProteinAnalysis) -> Double {
  if self.sequence.length() < 2 {
    return 0.0
  }
  let seq = self.sequence.to_upper()
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < seq.length() - 1 {
    let aa1 = seq.unsafe_get(i).unsafe_to_char()
    let aa2 = seq.unsafe_get(i + 1).unsafe_to_char()
    sum = sum + instability_score(aa1, aa2)
    count = count + 1
    i = i + 1
  }
  sum / count.to_double() * 10.0
}

///|
fn instability_score(a1 : Char, a2 : Char) -> Double {
  let aa = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
  let scores : Array[Array[Double]] = Array::new()
  
  scores.push([4.9, -1.8, -4.8, -3.5, 0.2, -8.0, -4.1, 0.0, -4.6, 4.2, 4.9, -3.9, 1.9, 2.8, -4.6, 0.4, 0.5, -0.9, -1.3, 4.2])
  scores.push([-1.8, 10.5, -2.8, -4.0, -5.5, -3.5, -1.0, -0.8, -3.2, 5.2, 5.9, 10.8, -2.5, 2.5, -4.0, 0.4, -1.5, -1.6, -2.8, -4.8])
  scores.push([-1.6, -6.1, -5.1, -1.0, -5.6, -3.0, -4.6, 0.2, -0.6, -4.7, -0.2, -1.8, 2.0, -3.5, -4.0, -6.2, -4.0, -0.8, -2.0, 0.3])
  scores.push([-3.6, -3.5, -6.1, -2.4, -4.0, -5.5, -0.8, -3.0, -4.1, -0.5, -5.5, -5.5, -2.8, -1.0, 0.4, -4.0, 0.8, -6.0, -1.0, -4.5])
  scores.push([-8.0, -3.5, -6.2, -7.0, -7.0, -0.8, -1.0, -5.0, -3.0, -3.5, -5.5, -5.0, -6.0, -4.5, -5.0, -4.1, -1.0, -6.2, -4.0, -5.4])
  scores.push([-6.0, -6.0, -2.5, -2.0, -3.0, -5.0, -4.0, -4.5, -5.5, -4.0, -5.0, 0.0, -0.8, -2.8, -3.0, -6.1, -5.9, -0.4, -1.8, -2.0])
  scores.push([-4.0, -3.0, -3.5, -5.0, 0.0, -4.6, -3.2, -5.4, -3.0, -2.5, -2.5, -3.2, -1.5, -1.0, -3.0, -4.0, -5.0, -4.5, -6.0, -4.2])
  scores.push([5.2, 4.2, -0.6, -2.5, 4.3, 4.9, 5.7, 2.5, 3.8, -1.6, -0.6, 0.0, -0.9, 3.8, 4.2, 4.9, 5.9, 4.9, -1.0, -2.5])
  scores.push([4.9, 4.9, 5.5, 3.8, 3.8, -1.6, -1.0, 0.0, -0.9, 3.8, 4.2, -3.9, 10.8, 5.5, -3.5, -3.2, 5.7, 5.5, 10.0, 2.8])
  scores.push([3.0, -4.5, -1.5, -2.0, -3.0, 2.5, 3.0, 1.9, -2.5, 3.8, -1.8, -1.5, 2.5, 3.8, 2.8, 5.2, 3.8, 3.8, -2.0, -1.8])
  scores.push([-1.0, 1.9, 2.8, 4.0, 2.8, 2.5, 3.8, -2.0, -1.0, 3.8, 3.8, 3.0, 3.8, 6.0, -1.0, -1.0, 0.0, 2.0, -4.6, -4.0])
  scores.push([-4.0, -5.5, -5.0, -4.0, -3.0, -1.6, -1.6, -4.5, -2.0, -1.0, -4.0, -3.5, -4.0, -2.0, -3.0, 0.4, 0.2, -0.6, 1.7, -1.5])
  scores.push([-0.6, -1.0, -1.5, -1.8, -1.0, -3.5, -2.5, -3.5, -3.0, -3.5, 0.5, -6.0, 0.0, -1.0, -1.0, 0.0, 0.0, -2.0, -1.0, -1.0])
  scores.push([-3.0, -2.0, -3.0, -3.5, -3.5, -0.9, -1.6, -0.9, -1.5, -3.0, -0.9, -0.9, -3.0, 1.9, 4.4, -4.0, -3.5, -2.5, -1.0, -2.0])
  scores.push([-1.3, -2.4, 3.8, -3.0, -1.0, 3.8, 3.8, -2.5, 2.0, 4.0, -3.0, -2.0, -3.0, -1.0, -1.5, 4.2, -2.8, -5.5, -0.4, -5.0])
  scores.push([4.2, 4.2, -2.8, -1.0, 3.8, -1.6, -0.4, -0.5, 3.8, 4.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
  scores.push([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
  scores.push([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
  scores.push([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
  scores.push([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
  
  let mut idx1 = -1
  let mut idx2 = -1
  let mut j = 0
  while j < 20 {
    if aa[j] == a1 { idx1 = j }
    if aa[j] == a2 { idx2 = j }
    j = j + 1
  }
  
  if idx1 >= 0 && idx2 >= 0 && idx1 < scores.length() {
    let row = scores[idx1]
    if idx2 < row.length() { row[idx2] } else { 0.0 }
  } else {
    0.0
  }
}

///|
let pka_values : Map[Char, (Double, Double, Double)] = Map(
  [
    ('A', (2.34, 9.69, 0.0)), ('R', (2.17, 9.04, 12.48)),
    ('N', (2.02, 8.80, 0.0)), ('D', (1.88, 9.60, 3.65)),
    ('C', (1.71, 10.28, 8.33)), ('E', (2.19, 9.67, 4.25)),
    ('Q', (2.17, 9.13, 0.0)), ('G', (2.34, 9.60, 0.0)),
    ('H', (1.82, 9.17, 6.00)), ('I', (2.36, 9.60, 0.0)),
    ('L', (2.36, 9.60, 0.0)), ('K', (2.18, 8.95, 10.53)),
    ('M', (2.28, 9.21, 0.0)), ('F', (1.83, 9.13, 0.0)),
    ('P', (1.99, 10.60, 0.0)), ('S', (2.21, 9.15, 0.0)),
    ('T', (2.09, 9.10, 0.0)), ('W', (2.83, 9.39, 0.0)),
    ('Y', (2.32, 9.62, 10.07)), ('V', (2.32, 9.62, 0.0)),
  ],
  capacity=20,
)

///|
pub fn ProteinAnalysis::isoelectric_point(self : ProteinAnalysis) -> Double {
  let mut ph_min = 0.0
  let mut ph_max = 14.0
  let mut ph_mid = 7.0
  let mut iter = 0
  while iter < 50 {
    ph_mid = (ph_min + ph_max) / 2.0
    let charge = self.charge_at_ph(ph_mid)
    if charge > 0.0 {
      ph_min = ph_mid
    } else {
      ph_max = ph_mid
    }
    iter = iter + 1
  }
  ph_mid
}

///|
pub fn ProteinAnalysis::charge_at_ph(self : ProteinAnalysis, ph : Double) -> Double {
  let seq = self.sequence.to_upper()
  let mut charge = 0.0
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).unsafe_to_char()
    match pka_values.get(c) {
      Some(p) => {
        let n_term = p.0
        let c_term = p.1
        let side = p.2
        if i == 0 {
          charge = charge + 1.0 / (1.0 + exp((ph - n_term) * 2.3026))
        }
        if i == seq.length() - 1 {
          charge = charge - 1.0 / (1.0 + exp((c_term - ph) * 2.3026))
        }
        if side > 0.0 {
          if side > 7.0 {
            charge = charge + 1.0 / (1.0 + exp((ph - side) * 2.3026))
          } else {
            charge = charge - 1.0 / (1.0 + exp((side - ph) * 2.3026))
          }
        }
      }
      None => ()
    }
    i = i + 1
  }
  charge
}

///|
let helix_p : Map[Char, Double] = Map(
  [
    ('A', 1.42), ('R', 1.05), ('N', 0.67), ('D', 1.01), ('C', 0.70),
    ('E', 1.51), ('Q', 1.11), ('G', 0.57), ('H', 1.00), ('I', 1.08),
    ('L', 1.21), ('K', 1.16), ('M', 1.45), ('F', 1.13), ('P', 0.57),
    ('S', 0.77), ('T', 0.83), ('W', 1.08), ('Y', 1.07), ('V', 1.06),
  ],
  capacity=20,
)

///|
let sheet_p : Map[Char, Double] = Map(
  [
    ('A', 0.83), ('R', 0.74), ('N', 0.89), ('D', 0.54), ('C', 1.19),
    ('E', 0.37), ('Q', 1.11), ('G', 0.75), ('H', 0.87), ('I', 1.60),
    ('L', 1.30), ('K', 0.74), ('M', 1.05), ('F', 1.38), ('P', 0.55),
    ('S', 0.75), ('T', 1.19), ('W', 1.37), ('Y', 1.47), ('V', 1.70),
  ],
  capacity=20,
)

///|
pub fn ProteinAnalysis::secondary_structure(self : ProteinAnalysis) -> (Double, Double, Double) {
  let seq = self.sequence.to_upper()
  let n = seq.length()
  if n == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut h = 0.0
  let mut s = 0.0
  let mut i = 0
  while i < n {
    let aa = seq.unsafe_get(i).unsafe_to_char()
    h = h + helix_p.get(aa).unwrap_or(1.0)
    s = s + sheet_p.get(aa).unwrap_or(1.0)
    i = i + 1
  }
  let len = n.to_double()
  let c = n.to_double() - h - s
  (h / len, s / len, c / len)
}

///|
pub fn ProteinAnalysis::amino_acid_composition(self : ProteinAnalysis) -> Map[Char, Double] {
  let counts = self.count_amino_acids()
  let total = self.sequence.length().to_double()
  let comp = Map([], capacity=20)
  for aa, count in counts {
    comp[aa] = count.to_double() / total
  }
  comp
}

///|
pub struct SignalPeptideResult {
  has_signal_peptide : Bool
  cleavage_site : Int?
  score : Double
}

///|
pub fn ProteinAnalysis::predict_signal_peptide(self : ProteinAnalysis) -> SignalPeptideResult {
  let seq = self.sequence.to_upper()
  let n = seq.length()
  if n < 20 {
    return SignalPeptideResult::{ has_signal_peptide: false, cleavage_site: None, score: 0.0 }
  }
  let mut score = 0.0
  let mut hydrophobic_score = 0.0
  let mut n_terminal_positive = 0
  let mut i = 0
  while i < min(15, n) {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if c == 'K' || c == 'R' {
      n_terminal_positive = n_terminal_positive + 1
    }
    i = i + 1
  }
  i = 8
  while i < min(25, n) {
    let c = seq.unsafe_get(i).unsafe_to_char()
    match gravy_values.get(c) {
      Some(v) => if v > 0.0 { hydrophobic_score = hydrophobic_score + v }
      None => ()
    }
    i = i + 1
  }
  if n_terminal_positive >= 1 {
    score = score + 1.0
  }
  if hydrophobic_score > 2.5 {
    score = score + hydrophobic_score / 5.0
  }
  let mut cleavage_score = 0.0
  let mut cleavage_pos = -1
  i = 15
  while i < min(n - 1, 30) {
    let c1 = seq.unsafe_get(i).unsafe_to_char()
    let c2 = seq.unsafe_get(i + 1).unsafe_to_char()
    if (c1 == 'A' || c1 == 'G' || c1 == 'S') && c2 == 'A' {
      cleavage_score = cleavage_score + 1.0
      cleavage_pos = i
    }
    i = i + 1
  }
  if cleavage_score > 0.0 {
    score = score + cleavage_score / 3.0
  }
  let has_signal = score > 1.5
  SignalPeptideResult::{
    has_signal_peptide: has_signal,
    cleavage_site: if cleavage_pos >= 0 && has_signal { Some(cleavage_pos) } else { None },
    score: score,
  }
}

///|
fn min(a : Int, b : Int) -> Int {
  if a < b { a } else { b }
}

///|
fn exp(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 50 {
    term = term * x / i.to_double()
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
pub fn ProteinAnalysis::length(self : ProteinAnalysis) -> Int {
  self.sequence.length()
}