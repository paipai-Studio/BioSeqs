///|
/// Biostrings - Sequence analysis functions inspired by R/Bioconductor Biostrings.

///|
fn log_nat(x : Double) -> Double {
  if x <= 0.0 {
    @double.not_a_number
  } else if x == 1.0 {
    0.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let y2 = y * y
    let mut sum = y
    let mut term = y
    let mut i = 1
    while i <= 20 {
      term = term * y2
      sum = sum + term / (2 * i + 1).to_double()
      i = i + 1
    }
    2.0 * sum
  }
}

///|
let iupac_nucleotide_map : Map[UInt16, Array[UInt16]] = Map([
  (65, [65]),
  (67, [67]),
  (71, [71]),
  (84, [84]),
  (85, [85]),
  (77, [65, 67]),
  (82, [65, 71]),
  (87, [65, 84]),
  (83, [67, 71]),
  (89, [67, 84]),
  (75, [71, 84]),
  (86, [65, 67, 71]),
  (72, [65, 67, 84]),
  (68, [65, 71, 84]),
  (66, [67, 71, 84]),
  (78, [65, 67, 71, 84]),
])

///|
pub struct NucleotideFreq {
  a : Int
  c : Int
  g : Int
  t : Int
  other : Int
  n : Int
}

///|
pub struct AminoAcidFreq {
  a : Int
  c : Int
  d : Int
  e : Int
  f : Int
  g : Int
  h : Int
  i : Int
  k : Int
  l : Int
  m : Int
  n : Int
  p : Int
  q : Int
  r : Int
  s : Int
  t : Int
  v : Int
  w : Int
  y : Int
  other : Int
}

///|
pub fn letter_frequency_dna(seq : String) -> NucleotideFreq {
  let mut a = 0
  let mut c = 0
  let mut g = 0
  let mut t = 0
  let mut n = 0
  let mut other = 0

  let mut i = 0
  while i < seq.length() {
    let ch = seq.unsafe_get(i)
    let upper = if ch >= 97 && ch <= 122 { ch - 32 } else { ch }
    if upper == 65 {
      a = a + 1
    } else if upper == 67 {
      c = c + 1
    } else if upper == 71 {
      g = g + 1
    } else if upper == 84 {
      t = t + 1
    } else if upper == 78 {
      n = n + 1
    } else {
      other = other + 1
    }
    i = i + 1
  }

  NucleotideFreq::{ a, c, g, t, n, other }
}

///|
pub fn letter_frequency_aa(seq : String) -> AminoAcidFreq {
  let counts : FixedArray[Int] = FixedArray::make(21, 0)

  let mut i = 0
  while i < seq.length() {
    let ch = seq.unsafe_get(i)
    let upper = if ch >= 97 && ch <= 122 { ch - 32 } else { ch }
    let idx = if upper == 65 {
      0
    } else if upper == 67 {
      1
    } else if upper == 68 {
      2
    } else if upper == 69 {
      3
    } else if upper == 70 {
      4
    } else if upper == 71 {
      5
    } else if upper == 72 {
      6
    } else if upper == 73 {
      7
    } else if upper == 75 {
      8
    } else if upper == 76 {
      9
    } else if upper == 77 {
      10
    } else if upper == 78 {
      11
    } else if upper == 80 {
      12
    } else if upper == 81 {
      13
    } else if upper == 82 {
      14
    } else if upper == 83 {
      15
    } else if upper == 84 {
      16
    } else if upper == 86 {
      17
    } else if upper == 87 {
      18
    } else if upper == 89 {
      19
    } else {
      20
    }
    counts[idx] = counts[idx] + 1
    i = i + 1
  }

  AminoAcidFreq::{
    a: counts[0],
    c: counts[1],
    d: counts[2],
    e: counts[3],
    f: counts[4],
    g: counts[5],
    h: counts[6],
    i: counts[7],
    k: counts[8],
    l: counts[9],
    m: counts[10],
    n: counts[11],
    p: counts[12],
    q: counts[13],
    r: counts[14],
    s: counts[15],
    t: counts[16],
    v: counts[17],
    w: counts[18],
    y: counts[19],
    other: counts[20],
  }
}

///|
pub fn oligonucleotide_frequency(
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

  let mut i = 0
  while i < total {
    let kmer = seq[i:i + k].to_owned()
    match freq.get(kmer) {
      Some(v) => freq[kmer] = v + 1.0
      None => freq[kmer] = 1.0
    }
    i = i + 1
  }

  if normalize && total > 0 {
    let keys = freq.keys().collect()
    let mut j = 0
    while j < keys.length() {
      freq[keys[j]] = freq[keys[j]] / total.to_double()
      j = j + 1
    }
  }

  freq
}

///|
pub fn dinucleotide_frequency(
  seq : String,
  normalize : Bool,
) -> Map[String, Double] {
  oligonucleotide_frequency(seq, 2, normalize)
}

///|
pub fn trinucleotide_frequency(
  seq : String,
  normalize : Bool,
) -> Map[String, Double] {
  oligonucleotide_frequency(seq, 3, normalize)
}

///|
pub fn rscu(seq : String) -> Map[String, Double] {
  let codon_freq = trinucleotide_frequency(seq, false)
  let aa_codons : Map[UInt16, Array[String]] = Map([
    (70, ["TTT", "TTC"]),
    (76, ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
    (73, ["ATT", "ATC", "ATA"]),
    (77, ["ATG"]),
    (86, ["GTT", "GTC", "GTA", "GTG"]),
    (83, ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
    (80, ["CCT", "CCC", "CCA", "CCG"]),
    (84, ["ACT", "ACC", "ACA", "ACG"]),
    (65, ["GCT", "GCC", "GCA", "GCG"]),
    (89, ["TAT", "TAC"]),
    (72, ["CAT", "CAC"]),
    (81, ["CAA", "CAG"]),
    (78, ["AAT", "AAC"]),
    (75, ["AAA", "AAG"]),
    (68, ["GAT", "GAC"]),
    (69, ["GAA", "GAG"]),
    (67, ["TGT", "TGC"]),
    (87, ["TGG"]),
    (82, ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
    (71, ["GGT", "GGC", "GGA", "GGG"]),
  ])

  let rscu_result : Map[String, Double] = Map([], capacity=64)

  let aa_keys = aa_codons.keys().collect()
  let mut i = 0
  while i < aa_keys.length() {
    let aa = aa_keys[i]
    let codons = aa_codons[aa]
    let n_codons = codons.length().to_double()

    let mut total = 0.0
    let mut j = 0
    while j < codons.length() {
      total = total + codon_freq.get(codons[j]).unwrap_or(0.0)
      j = j + 1
    }

    let mut k = 0
    while k < codons.length() {
      let codon = codons[k]
      let observed = codon_freq.get(codon).unwrap_or(0.0)
      let expected = if total > 0.0 { total / n_codons } else { 1.0 }
      rscu_result[codon] = observed / expected
      k = k + 1
    }
    i = i + 1
  }

  rscu_result
}

///|
pub fn gc_content_by_position(seq : String, window_size : Int) -> Array[Double] {
  let n = seq.length()
  if window_size <= 0 || n == 0 {
    return Array::new()
  }

  let result = Array::new()
  let step = window_size
  let mut start = 0

  while start < n {
    let end = if start + window_size <= n { start + window_size } else { n }
    let window = seq[start:end]

    let mut gc = 0
    let mut total = 0
    let mut i = 0
    while i < window.length() {
      let ch = window.unsafe_get(i)
      let upper = if ch >= 97 && ch <= 122 { ch - 32 } else { ch }
      if upper == 71 || upper == 67 {
        gc = gc + 1
      }
      total = total + 1
      i = i + 1
    }

    result.push(
      if total > 0 {
        gc.to_double() / total.to_double()
      } else {
        0.0
      },
    )
    start = start + step
  }

  result
}

///|
pub fn shannon_entropy(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let freq : Map[UInt16, Int] = Map([], capacity=4)

  let mut i = 0
  while i < n {
    let ch = seq.unsafe_get(i)
    let upper = if ch >= 97 && ch <= 122 { ch - 32 } else { ch }
    freq[upper] = freq.get(upper).unwrap_or(0) + 1
    i = i + 1
  }

  let mut entropy = 0.0
  let keys = freq.keys().collect()
  let mut j = 0
  while j < keys.length() {
    let prob = freq[keys[j]].to_double() / n.to_double()
    entropy = entropy - prob * (log_nat(prob) / log_nat(2.0))
    j = j + 1
  }

  entropy
}

///|
pub fn dust_complexity(seq : String) -> Double {
  let n = seq.length()
  if n < 3 {
    return 0.0
  }

  let mut score = 0.0
  let window_size = 64

  let mut i = 0
  while i <= n - window_size {
    let window = seq[i:i + window_size]

    let triples : Map[String, Int] = Map([], capacity=64)
    let mut j = 0
    while j <= window_size - 3 {
      let triple = window[j:j + 3].to_owned()
      triples[triple] = triples.get(triple).unwrap_or(0) + 1
      j = j + 1
    }

    let mut window_score = 0.0
    let keys = triples.keys().collect()
    let mut k = 0
    while k < keys.length() {
      let count = triples[keys[k]]
      if count > 1 {
        window_score = window_score +
          count.to_double() * (count.to_double() - 1.0) / 2.0
      }
      k = k + 1
    }

    score = score + window_score / ((window_size - 2).to_double() / 3.0)
    i = i + window_size
  }

  if score > 0.0 {
    score = score * (100.0 / (n.to_double() / window_size.to_double()).ceil())
  }

  score
}

///|
pub fn iupac_match(pattern : String, subject : String) -> Array[(Int, Int)] {
  let result = Array::new()
  let n = subject.length()
  let m = pattern.length()

  if m == 0 || n < m {
    return result
  }

  let mut i = 0
  while i <= n - m {
    let mut matches = true
    let mut j = 0
    while j < m {
      let p_ch = pattern.unsafe_get(j)
      let p_upper = if p_ch >= 97 && p_ch <= 122 { p_ch - 32 } else { p_ch }
      let s_ch = subject.unsafe_get(i + j)
      let s_upper = if s_ch >= 97 && s_ch <= 122 { s_ch - 32 } else { s_ch }

      if !iupac_compatible(p_upper, s_upper) {
        matches = false
        break
      }
      j = j + 1
    }

    if matches {
      result.push((i, i + m))
    }
    i = i + 1
  }

  result
}

///|
fn iupac_compatible(c1 : UInt16, c2 : UInt16) -> Bool {
  let bases1 = iupac_nucleotide_map.get(c1).unwrap_or([0])
  let bases2 = iupac_nucleotide_map.get(c2).unwrap_or([0])

  let mut i = 0
  while i < bases1.length() {
    let mut j = 0
    while j < bases2.length() {
      if bases1[i] == bases2[j] {
        return true
      }
      j = j + 1
    }
    i = i + 1
  }

  false
}

///|
pub fn tm_nn(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }

  let mut dH = 0.0
  let mut dS = 0.0

  let nn_params : Map[String, (Double, Double)] = Map([
    ("AA", (-7.9, -22.2)),
    ("TT", (-7.9, -22.2)),
    ("AT", (-7.2, -20.4)),
    ("TA", (-7.2, -20.4)),
    ("CA", (-8.5, -22.7)),
    ("TG", (-8.5, -22.7)),
    ("GT", (-8.4, -22.4)),
    ("AC", (-8.4, -22.4)),
    ("CT", (-7.8, -21.0)),
    ("AG", (-7.8, -21.0)),
    ("GA", (-8.2, -22.2)),
    ("TC", (-8.2, -22.2)),
    ("CG", (-10.6, -27.2)),
    ("GC", (-9.8, -24.4)),
    ("GG", (-8.0, -19.9)),
    ("CC", (-8.0, -19.9)),
  ])

  let mut i = 0
  while i < n - 1 {
    let dinucleotide = seq[i:i + 2].to_owned().to_upper()
    match nn_params.get(dinucleotide) {
      Some((h, s)) => {
        dH = dH + h
        dS = dS + s
      }
      None => ()
    }
    i = i + 1
  }

  let tm = 1000.0 *
    dH /
    (dS + 1.987 * (4.6 + log_nat(1.0 / n.to_double()) / log_nat(10.0))) -
    273.15

  tm.max(0.0)
}

///|
pub fn nedit_at(pattern : String, subject : String, fixed : Bool) -> Int {
  let n = pattern.length()
  let m = subject.length()
  let min_len = if n < m { n } else { m }

  let mut edits = 0

  let mut i = 0
  while i < min_len {
    let p_ch = pattern.unsafe_get(i)
    let p_upper = if p_ch >= 97 && p_ch <= 122 { p_ch - 32 } else { p_ch }
    let s_ch = subject.unsafe_get(i)
    let s_upper = if s_ch >= 97 && s_ch <= 122 { s_ch - 32 } else { s_ch }

    if fixed {
      if p_upper != s_upper {
        edits = edits + 1
      }
    } else if !iupac_compatible(p_upper, s_upper) {
      edits = edits + 1
    }
    i = i + 1
  }

  edits + (n - min_len).max(m - min_len)
}

///|
pub fn percent_identity(seq1 : String, seq2 : String) -> Double {
  let n = seq1.length()
  let m = seq2.length()
  let min_len = if n < m { n } else { m }

  if min_len == 0 {
    return 0.0
  }

  let mut matches = 0

  let mut i = 0
  while i < min_len {
    let c1 = seq1.unsafe_get(i)
    let c1_upper = if c1 >= 97 && c1 <= 122 { c1 - 32 } else { c1 }
    let c2 = seq2.unsafe_get(i)
    let c2_upper = if c2 >= 97 && c2 <= 122 { c2 - 32 } else { c2 }

    if c1_upper == c2_upper {
      matches = matches + 1
    }
    i = i + 1
  }

  matches.to_double() / min_len.to_double() * 100.0
}

///|
pub fn iupac_reverse_complement(seq : String) -> String {
  let n = seq.length()
  let result : FixedArray[UInt16] = FixedArray::make(n, 63)

  let mut i = 0
  while i < n {
    let c = seq.unsafe_get(i)
    let complement = if c == 65 {
      84
    } else if c == 84 {
      65
    } else if c == 71 {
      67
    } else if c == 67 {
      71
    } else if c == 85 {
      65
    } else if c == 77 {
      75
    } else if c == 82 {
      89
    } else if c == 87 {
      87
    } else if c == 83 {
      83
    } else if c == 89 {
      82
    } else if c == 75 {
      77
    } else if c == 86 {
      66
    } else if c == 72 {
      68
    } else if c == 68 {
      72
    } else if c == 66 {
      86
    } else if c == 78 {
      78
    } else if c == 97 {
      116
    } else if c == 116 {
      97
    } else if c == 103 {
      99
    } else if c == 99 {
      103
    } else if c == 117 {
      97
    } else if c == 109 {
      107
    } else if c == 114 {
      121
    } else if c == 119 {
      119
    } else if c == 115 {
      115
    } else if c == 121 {
      114
    } else if c == 107 {
      109
    } else if c == 118 {
      98
    } else if c == 104 {
      100
    } else if c == 100 {
      104
    } else if c == 98 {
      118
    } else if c == 110 {
      110
    } else {
      63
    }
    result[n - 1 - i] = complement.to_uint16()
    i = i + 1
  }

  unsafe_from_uint16_fixedarray(result)
}
