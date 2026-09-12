///|
/// Machine learning feature extraction module for biological sequences.

///|
pub struct NucleotideComposition {
  a : Double
  t : Double
  g : Double
  c : Double
  other : Double
}

///|
pub fn kmer_frequency(
  seq : String,
  k : Int,
  normalize : Bool,
) -> Map[String, Double] {
  let n = seq.length()
  if k <= 0 || k > n {
    return Map([], capacity=0)
  }

  let total = n - k + 1
  let freq : Map[String, Double] = Map([], capacity=total)

  for i in 0..<total {
    let kmer = seq[i:i + k].to_owned()
    match freq.get(kmer) {
      Some(v) => freq[kmer] = v + 1.0
      None => freq[kmer] = 1.0
    }
  }

  if normalize && total > 0 {
    let keys = freq.keys().collect()
    for i in 0..<keys.length() {
      freq[keys[i]] = freq[keys[i]] / total.to_double()
    }
  }

  freq
}

///|
pub fn dna_canonical_kmer_frequency(
  seq : String,
  k : Int,
  normalize : Bool,
) -> Map[String, Double] {
  let n = seq.length()
  if k <= 0 || k > n {
    return Map([], capacity=0)
  }

  let total = n - k + 1
  let freq : Map[String, Double] = Map([], capacity=total)

  for i in 0..<total {
    let kmer = seq[i:i + k].to_owned()
    let rc = reverse_complement_str(kmer)

    let canonical = if kmer <= rc { kmer } else { rc }

    match freq.get(canonical) {
      Some(v) => freq[canonical] = v + 1.0
      None => freq[canonical] = 1.0
    }
  }

  if normalize && total > 0 {
    let keys = freq.keys().collect()
    for i in 0..<keys.length() {
      freq[keys[i]] = freq[keys[i]] / total.to_double()
    }
  }

  freq
}

///|
fn reverse_complement_str(s : String) -> String {
  let n = s.length()
  let mut result = ""
  for i in 0..<n {
    let idx = n - 1 - i
    let c = s.unsafe_get(idx)
    let complement = match c {
      65 => 'T'
      84 => 'A'
      71 => 'C'
      67 => 'G'
      97 => 't'
      116 => 'a'
      103 => 'c'
      99 => 'g'
      _ => c.to_char().unwrap_or('?')
    }
    result = "\{result}\{complement}"
  }
  result
}

///|
pub fn nucleotide_composition(seq : String) -> NucleotideComposition {
  let n = seq.length()
  if n == 0 {
    return NucleotideComposition::{
      a: 0.0,
      t: 0.0,
      g: 0.0,
      c: 0.0,
      other: 0.0,
    }
  }

  let mut a_count = 0.0
  let mut t_count = 0.0
  let mut g_count = 0.0
  let mut c_count = 0.0
  let mut o_count = 0.0

  for i in 0..<n {
    let ch = seq.unsafe_get(i)
    if ch == 65 || ch == 97 {
      a_count = a_count + 1.0
    } else if ch == 84 || ch == 116 {
      t_count = t_count + 1.0
    } else if ch == 71 || ch == 103 {
      g_count = g_count + 1.0
    } else if ch == 67 || ch == 99 {
      c_count = c_count + 1.0
    } else {
      o_count = o_count + 1.0
    }
  }

  let len = n.to_double()
  NucleotideComposition::{
    a: a_count / len,
    t: t_count / len,
    g: g_count / len,
    c: c_count / len,
    other: o_count / len,
  }
}

///|
pub fn amino_acid_composition(seq : String) -> Map[UInt16, Double] {
  let n = seq.length()
  if n == 0 {
    return Map([], capacity=0)
  }

  let freq : Map[UInt16, Double] = Map([], capacity=20)
  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    match freq.get(c) {
      Some(v) => freq[c] = v + 1.0
      None => freq[c] = 1.0
    }
  }

  let len = n.to_double()
  let keys = freq.keys().collect()
  for i in 0..<keys.length() {
    freq[keys[i]] = freq[keys[i]] / len
  }

  freq
}

///|
pub fn dipeptide_composition(seq : String) -> Map[String, Double] {
  let n = seq.length()
  if n < 2 {
    return Map([], capacity=0)
  }

  let total = n - 1
  let freq : Map[String, Double] = Map([], capacity=total)
  let upper_seq = seq.to_upper()

  for i in 0..<total {
    let aa1 = upper_seq.unsafe_get(i).to_char().unwrap_or('?')
    let aa2 = upper_seq.unsafe_get(i + 1).to_char().unwrap_or('?')
    let dipeptide = "\{aa1}\{aa2}"

    match freq.get(dipeptide) {
      Some(v) => freq[dipeptide] = v + 1.0
      None => freq[dipeptide] = 1.0
    }
  }

  let len = total.to_double()
  let keys = freq.keys().collect()
  for i in 0..<keys.length() {
    freq[keys[i]] = freq[keys[i]] / len
  }

  freq
}

///|
pub fn tripeptide_composition(seq : String) -> Map[String, Double] {
  let n = seq.length()
  if n < 3 {
    return Map([], capacity=0)
  }

  let total = n - 2
  let freq : Map[String, Double] = Map([], capacity=total)
  let upper_seq = seq.to_upper()

  for i in 0..<total {
    let aa1 = upper_seq.unsafe_get(i).to_char().unwrap_or('?')
    let aa2 = upper_seq.unsafe_get(i + 1).to_char().unwrap_or('?')
    let aa3 = upper_seq.unsafe_get(i + 2).to_char().unwrap_or('?')
    let tripeptide = "\{aa1}\{aa2}\{aa3}"

    match freq.get(tripeptide) {
      Some(v) => freq[tripeptide] = v + 1.0
      None => freq[tripeptide] = 1.0
    }
  }

  let len = total.to_double()
  let keys = freq.keys().collect()
  for i in 0..<keys.length() {
    freq[keys[i]] = freq[keys[i]] / len
  }

  freq
}

///|
pub let aa_hydrophobicity : Map[UInt16, Double] = {
  let m : Map[UInt16, Double] = Map([], capacity=20)
  m['A'] = 1.8
  m['R'] = -4.5
  m['N'] = -3.5
  m['D'] = -3.5
  m['C'] = 2.5
  m['E'] = -3.5
  m['Q'] = -3.5
  m['G'] = -0.4
  m['H'] = -3.2
  m['I'] = 4.5
  m['L'] = 3.8
  m['K'] = -3.9
  m['M'] = 1.9
  m['F'] = 2.8
  m['P'] = -1.6
  m['S'] = -0.8
  m['T'] = -0.7
  m['W'] = -0.9
  m['Y'] = -1.3
  m['V'] = 4.2
  m
}

///|
pub let aa_charge : Map[UInt16, Double] = {
  let m : Map[UInt16, Double] = Map([], capacity=20)
  m['A'] = 0.0
  m['R'] = 1.0
  m['N'] = 0.0
  m['D'] = -1.0
  m['C'] = 0.0
  m['E'] = -1.0
  m['Q'] = 0.0
  m['G'] = 0.0
  m['H'] = 0.0
  m['I'] = 0.0
  m['L'] = 0.0
  m['K'] = 1.0
  m['M'] = 0.0
  m['F'] = 0.0
  m['P'] = 0.0
  m['S'] = 0.0
  m['T'] = 0.0
  m['W'] = 0.0
  m['Y'] = 0.0
  m['V'] = 0.0
  m
}

///|
pub let aa_polarity : Map[UInt16, Double] = {
  let m : Map[UInt16, Double] = Map([], capacity=20)
  m['A'] = 0.0
  m['R'] = 1.0
  m['N'] = 1.0
  m['D'] = 1.0
  m['C'] = 0.5
  m['E'] = 1.0
  m['Q'] = 1.0
  m['G'] = 0.0
  m['H'] = 1.0
  m['I'] = 0.0
  m['L'] = 0.0
  m['K'] = 1.0
  m['M'] = 0.0
  m['F'] = 0.0
  m['P'] = 0.0
  m['S'] = 1.0
  m['T'] = 1.0
  m['W'] = 0.5
  m['Y'] = 0.5
  m['V'] = 0.0
  m
}

///|
pub fn avg_hydrophobicity(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let mut total = 0.0
  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    total = total + aa_hydrophobicity.get(c).unwrap_or(0.0)
  }

  total / n.to_double()
}

///|
pub fn avg_charge(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let mut total = 0.0
  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    total = total + aa_charge.get(c).unwrap_or(0.0)
  }

  total / n.to_double()
}

///|
pub fn avg_polarity(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let mut total = 0.0
  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    total = total + aa_polarity.get(c).unwrap_or(0.0)
  }

  total / n.to_double()
}

///|
pub fn molecular_weight(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let mut weight = 18.0

  let aa_weights : Map[UInt16, Double] = Map([], capacity=20)
  aa_weights['A'] = 89.09
  aa_weights['R'] = 174.20
  aa_weights['N'] = 132.12
  aa_weights['D'] = 133.10
  aa_weights['C'] = 121.15
  aa_weights['E'] = 147.13
  aa_weights['Q'] = 146.15
  aa_weights['G'] = 75.07
  aa_weights['H'] = 155.16
  aa_weights['I'] = 131.17
  aa_weights['L'] = 131.17
  aa_weights['K'] = 146.19
  aa_weights['M'] = 149.21
  aa_weights['F'] = 165.19
  aa_weights['P'] = 115.13
  aa_weights['S'] = 105.09
  aa_weights['T'] = 119.12
  aa_weights['W'] = 204.23
  aa_weights['Y'] = 181.19
  aa_weights['V'] = 117.15

  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    weight = weight + aa_weights.get(c).unwrap_or(0.0) - 18.0
  }

  weight
}

///|
pub fn secondary_structure_propensity(seq : String) -> (Double, Double, Double) {
  let n = seq.length()
  if n == 0 {
    return (0.0, 0.0, 0.0)
  }

  let mut helix = 0.0
  let mut sheet = 0.0
  let mut coil = 0.0

  let helix_p : Map[UInt16, Double] = Map([], capacity=20)
  helix_p['A'] = 1.42
  helix_p['R'] = 1.05
  helix_p['N'] = 0.67
  helix_p['D'] = 1.01
  helix_p['C'] = 0.70
  helix_p['E'] = 1.51
  helix_p['Q'] = 1.11
  helix_p['G'] = 0.57
  helix_p['H'] = 1.00
  helix_p['I'] = 1.08
  helix_p['L'] = 1.21
  helix_p['K'] = 1.16
  helix_p['M'] = 1.45
  helix_p['F'] = 1.13
  helix_p['P'] = 0.57
  helix_p['S'] = 0.77
  helix_p['T'] = 0.83
  helix_p['W'] = 1.08
  helix_p['Y'] = 1.07
  helix_p['V'] = 1.06

  let sheet_p : Map[UInt16, Double] = Map([], capacity=20)
  sheet_p['A'] = 0.83
  sheet_p['R'] = 0.74
  sheet_p['N'] = 0.89
  sheet_p['D'] = 0.54
  sheet_p['C'] = 1.19
  sheet_p['E'] = 0.37
  sheet_p['Q'] = 1.11
  sheet_p['G'] = 0.75
  sheet_p['H'] = 0.87
  sheet_p['I'] = 1.60
  sheet_p['L'] = 1.30
  sheet_p['K'] = 0.74
  sheet_p['M'] = 1.05
  sheet_p['F'] = 1.38
  sheet_p['P'] = 0.55
  sheet_p['S'] = 0.75
  sheet_p['T'] = 1.19
  sheet_p['W'] = 1.37
  sheet_p['Y'] = 1.47
  sheet_p['V'] = 1.70

  let coil_p : Map[UInt16, Double] = Map([], capacity=20)
  coil_p['A'] = 1.00
  coil_p['R'] = 1.00
  coil_p['N'] = 1.00
  coil_p['D'] = 1.00
  coil_p['C'] = 1.00
  coil_p['E'] = 1.00
  coil_p['Q'] = 1.00
  coil_p['G'] = 1.00
  coil_p['H'] = 1.00
  coil_p['I'] = 1.00
  coil_p['L'] = 1.00
  coil_p['K'] = 1.00
  coil_p['M'] = 1.00
  coil_p['F'] = 1.00
  coil_p['P'] = 1.00
  coil_p['S'] = 1.00
  coil_p['T'] = 1.00
  coil_p['W'] = 1.00
  coil_p['Y'] = 1.00
  coil_p['V'] = 1.00

  let upper_seq = seq.to_upper()

  for i in 0..<n {
    let c = upper_seq.unsafe_get(i)
    helix = helix + helix_p.get(c).unwrap_or(0.0)
    sheet = sheet + sheet_p.get(c).unwrap_or(0.0)
    coil = coil + coil_p.get(c).unwrap_or(0.0)
  }

  let len = n.to_double()
  (helix / len, sheet / len, coil / len)
}

///|
pub fn position_specific_features(
  seq : String,
  n_positions : Int,
) -> Array[Double] {
  let len = seq.length()
  if len == 0 || n_positions <= 0 {
    return Array::new()
  }

  let features : Array[Double] = Array::new()
  let upper_seq = seq.to_upper()

  for i in 0..<n_positions {
    let idx = i * len / n_positions
    let actual_idx = if idx >= len { len - 1 } else { idx }
    let c = upper_seq.unsafe_get(actual_idx)

    features.push(aa_hydrophobicity.get(c).unwrap_or(0.0))
    features.push(aa_charge.get(c).unwrap_or(0.0))
    features.push(aa_polarity.get(c).unwrap_or(0.0))
  }

  features
}

///|
pub fn protein_feature_vector(seq : String) -> Array[Double] {
  let features : Array[Double] = Array::new()

  let aa_order : Array[UInt16] = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
    'T', 'W', 'Y', 'V',
  ]
  let aa_comp = amino_acid_composition(seq)
  for i in 0..<aa_order.length() {
    features.push(aa_comp.get(aa_order[i]).unwrap_or(0.0))
  }

  let dipep = dipeptide_composition(seq)
  let dipep_keys = dipep.keys().collect()
  let dipep_limit = min_int(dipep_keys.length(), 20)
  for i in 0..<dipep_limit {
    features.push(dipep.get(dipep_keys[i]).unwrap_or(0.0))
  }

  features.push(avg_hydrophobicity(seq))
  features.push(avg_charge(seq))
  features.push(avg_polarity(seq))
  features.push(molecular_weight(seq))

  let ss = secondary_structure_propensity(seq)
  features.push(ss.0)
  features.push(ss.1)
  features.push(ss.2)

  let pos_feat = position_specific_features(seq, 10)
  for i in 0..<pos_feat.length() {
    features.push(pos_feat[i])
  }

  features
}

///|
pub fn dna_feature_vector(seq : String) -> Array[Double] {
  let features : Array[Double] = Array::new()

  let comp = nucleotide_composition(seq)
  features.push(comp.a)
  features.push(comp.t)
  features.push(comp.g)
  features.push(comp.c)

  let dinuc = kmer_frequency(seq, 2, true)
  let dinuc_order = [
    "AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA",
    "CT", "CG", "CC",
  ]
  for i in 0..<dinuc_order.length() {
    features.push(dinuc.get(dinuc_order[i]).unwrap_or(0.0))
  }

  let trinuc = kmer_frequency(seq, 3, true)
  let trinuc_keys = trinuc.keys().collect()
  let trinuc_limit = min_int(trinuc_keys.length(), 20)
  for i in 0..<trinuc_limit {
    features.push(trinuc.get(trinuc_keys[i]).unwrap_or(0.0))
  }

  features.push(comp.g + comp.c)
  features.push(seq.length().to_double())

  features
}

///|
pub fn protein_feature_names() -> Array[String] {
  let names : Array[String] = Array::new()

  let aa_order : Array[UInt16] = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
    'T', 'W', 'Y', 'V',
  ]
  for i in 0..<aa_order.length() {
    names.push("aa_\{aa_order[i].to_char().unwrap_or('?')}")
  }

  for i in 0..<20 {
    names.push("dipep_" + i.to_string())
  }

  names.push("avg_hydrophobicity")
  names.push("avg_charge")
  names.push("avg_polarity")
  names.push("molecular_weight")
  names.push("helix_propensity")
  names.push("sheet_propensity")
  names.push("coil_propensity")

  for i in 0..<10 {
    names.push("pos_" + i.to_string() + "_hydrophobicity")
    names.push("pos_" + i.to_string() + "_charge")
    names.push("pos_" + i.to_string() + "_polarity")
  }

  names
}

///|
pub fn dna_feature_names() -> Array[String] {
  let names : Array[String] = Array::new()

  names.push("A")
  names.push("T")
  names.push("G")
  names.push("C")

  let dinuc_order = [
    "AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC", "GA", "GT", "GG", "GC", "CA",
    "CT", "CG", "CC",
  ]
  for i in 0..<dinuc_order.length() {
    names.push("dinuc_" + dinuc_order[i])
  }

  for i in 0..<20 {
    names.push("trinuc_" + i.to_string())
  }

  names.push("GC_content")
  names.push("length")

  names
}

///|
fn min_int(a : Int, b : Int) -> Int {
  if a < b {
    a
  } else {
    b
  }
}
