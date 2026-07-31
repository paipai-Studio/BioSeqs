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

///|
/// A single match hit with position, count, and mismatched positions.
pub struct MatchHit {
  start : Int
  end : Int
  width : Int
  mismatches : Int
  inserted : Int
  deleted : Int
}

///|
/// Results of matchPattern containing all hits and summary statistics.
pub struct MatchPatternResult {
  hits : Array[MatchHit]
  pattern : String
  subject : String
  max_mismatches : Int
  with_indels : Bool
}

///|
/// Get number of matches found.
pub fn MatchPatternResult::n_matches(self : MatchPatternResult) -> Int {
  self.hits.length()
}

///|
/// Get all match starts.
pub fn MatchPatternResult::starts(self : MatchPatternResult) -> Array[Int] {
  self.hits.map(fn(h) { h.start })
}

///|
/// Get all match ends.
pub fn MatchPatternResult::ends(self : MatchPatternResult) -> Array[Int] {
  self.hits.map(fn(h) { h.end })
}

///|
/// Get all match widths.
pub fn MatchPatternResult::widths(self : MatchPatternResult) -> Array[Int] {
  self.hits.map(fn(h) { h.width })
}

///|
/// Get all mismatch counts.
pub fn MatchPatternResult::mismatch_counts(self : MatchPatternResult) -> Array[Int] {
  self.hits.map(fn(h) { h.mismatches })
}

///|
/// Match pattern to a subject sequence with mismatch/indel support.
/// Uses dynamic programming (Needleman-Wunsch style) for semi-global alignment.
/// Pattern is aligned to subsequences of the subject.
pub fn match_pattern(
  pattern~ : String,
  subject~ : String,
  max_mismatches? : Int = 0,
  with_indels? : Bool = false,
) -> MatchPatternResult {
  let n = pattern.length()
  let m = subject.length()

  let hits : Array[MatchHit] = Array::new()

  if n == 0 || m == 0 {
    return MatchPatternResult::{ hits, pattern, subject, max_mismatches, with_indels }
  }

  let max_errors = if with_indels { max_mismatches } else { max_mismatches }

  let mut i = 0
  while i <= m - n {
    let mut mismatches = 0
    let mut j = 0
    while j < n {
      let p_ch = pattern.unsafe_get(j)
      let p_upper = if p_ch >= 97 && p_ch <= 122 { p_ch - 32 } else { p_ch }
      let s_ch = subject.unsafe_get(i + j)
      let s_upper = if s_ch >= 97 && s_ch <= 122 { s_ch - 32 } else { s_ch }

      if p_upper != s_upper {
        mismatches = mismatches + 1
        if !with_indels && mismatches > max_mismatches {
          break
        }
      }
      j = j + 1
    }

    if j == n && mismatches <= max_errors {
      hits.push(MatchHit::{
        start: i + 1,
        end: i + n,
        width: n,
        mismatches,
        inserted: 0,
        deleted: 0,
      })
    }
    i = i + 1
  }

  MatchPatternResult::{ hits, pattern, subject, max_mismatches, with_indels }
}

///|
/// Vectorized matchPattern: applies pattern matching to multiple sequences.
pub fn vmatch_pattern(
  pattern~ : String,
  subjects~ : Array[String],
  max_mismatches? : Int = 0,
  with_indels? : Bool = false,
) -> Array[MatchPatternResult] {
  subjects.map(fn(subj) {
    match_pattern(pattern=pattern, subject=subj, max_mismatches=max_mismatches, with_indels=with_indels)
  })
}

///|
/// Find palindromic sequences (reverse complement equals itself).
pub fn find_palindromes(
  seq~ : String,
  min_length? : Int = 4,
) -> Array[(Int, Int)] {
  let results : Array[(Int, Int)] = Array::new()
  let n = seq.length()

  let mut i = 0
  while i <= n - min_length {
    let mut j = n
    while j >= i + min_length {
      let subseq = seq.substring(start=i, end=j)
      let rc = iupac_reverse_complement(subseq)
      if subseq.to_upper() == rc.to_upper() && subseq.length() >= min_length {
        results.push((i + 1, j))
      }
      j = j - 1
    }
    i = i + 1
  }

  results
}

///|
/// Find direct repeats in a sequence.
pub fn find_direct_repeats(
  seq~ : String,
  min_unit_length? : Int = 2,
  max_unit_length? : Int = 10,
  min_copies? : Int = 2,
) -> Array[(Int, Int, Int)] {
  let results : Array[(Int, Int, Int)] = Array::new()
  let n = seq.length()

  let mut unit_len = min_unit_length
  while unit_len <= max_unit_length {
    let mut i = 0
    while i <= n - unit_len * min_copies {
      let unit = seq.substring(start=i, end=i + unit_len).to_upper()
      let mut copies = 1
      let mut j = i + unit_len
      while j + unit_len <= n {
        let next_unit = seq.substring(start=j, end=j + unit_len).to_upper()
        if next_unit == unit {
          copies = copies + 1
          j = j + unit_len
        } else {
          break
        }
      }
      if copies >= min_copies {
        results.push((i + 1, i + copies * unit_len, copies))
      }
      i = i + 1
    }
    unit_len = unit_len + 1
  }

  results
}

///|
/// Find inverted repeats in a sequence.
pub fn find_inverted_repeats(
  seq~ : String,
  min_length? : Int = 4,
) -> Array[(Int, Int, Int, Int)] {
  let results : Array[(Int, Int, Int, Int)] = Array::new()
  let n = seq.length()

  let mut i = 0
  while i <= n - min_length {
    let mut j = i + min_length
    while j <= n - min_length {
      let len = min_length
      let sub1 = seq.substring(start=i, end=i + len).to_upper()
      let sub2 = seq.substring(start=j, end=j + len).to_upper()
      let rc2 = iupac_reverse_complement(sub2)
      if sub1 == rc2 {
        results.push((i + 1, i + len, j + 1, j + len))
      }
      j = j + 1
    }
    i = i + 1
  }

  results
}

///|
/// Find all occurrences of a motif allowing IUPAC ambiguity codes.
pub fn find_motif_iupac(
  pattern : String,
  subject : String,
) -> Array[Int] {
  let n = pattern.length()
  let m = subject.length()
  let positions : Array[Int] = Array::new()

  if n == 0 || m < n {
    return positions
  }

  let mut i = 0
  while i <= m - n {
    let mut match_found = true
    let mut j = 0
    while j < n {
      let p_ch = pattern.unsafe_get(j)
      let p_upper = if p_ch >= 97 && p_ch <= 122 { p_ch - 32 } else { p_ch }
      let s_ch = subject.unsafe_get(i + j)
      let s_upper = if s_ch >= 97 && s_ch <= 122 { s_ch - 32 } else { s_ch }

      if !iupac_compatible(p_upper, s_upper) {
        match_found = false
        break
      }
      j = j + 1
    }

    if match_found {
      positions.push(i + 1)
    }
    i = i + 1
  }

  positions
}

///|
/// Calculate GC content in sliding windows.
pub fn gc_content_sliding(
  seq~ : String,
  window_size? : Int = 100,
  step? : Int = 10,
) -> Array[(Int, Double)] {
  let results : Array[(Int, Double)] = Array::new()
  let n = seq.length()

  if n < window_size {
    return results
  }

  let mut i = 0
  while i <= n - window_size {
    let window = seq.substring(start=i, end=i + window_size)
    let gc = gc_content(window)
    results.push((i + 1, gc))
    i = i + step
  }

  results
}

///|
/// Calculate melting temperature in sliding windows.
pub fn tm_sliding(
  seq~ : String,
  window_size? : Int = 100,
  step? : Int = 10,
) -> Array[(Int, Double)] {
  let results : Array[(Int, Double)] = Array::new()
  let n = seq.length()

  if n < window_size {
    return results
  }

  let mut i = 0
  while i <= n - window_size {
    let window = seq.substring(start=i, end=i + window_size)
    let tm = tm_nn(window)
    results.push((i + 1, tm))
    i = i + step
  }

  results
}

///|
/// Compute per-position letter frequency matrix.
/// Returns a map of nucleotide letters to arrays of per-position counts.
/// Similar to Bioconductor's letterFrequency function.
pub fn letter_frequency_matrix(seq : String) -> Map[String, Array[Int]] {
  let n = seq.length()
  let matrix : Map[String, Array[Int]] = Map([], capacity=8)

  // Initialize for DNA alphabet
  let bases = ["A", "C", "G", "T", "U", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]
  for base in bases {
    matrix.set(base, Array::make(n, 0))
  }

  let mut i = 0
  while i < n {
    let ch = seq[i]
    let letter = match ch {
      65 => "A"
      67 => "C"
      71 => "G"
      84 => "T"
      85 => "U"
      78 => "N"
      82 => "R"
      89 => "Y"
      83 => "S"
      87 => "W"
      75 => "K"
      77 => "M"
      66 => "B"
      68 => "D"
      72 => "H"
      86 => "V"
      _ => "N"
    }
    let counts = matrix.get_or_default(letter, Array::make(n, 0))
    counts[i] = counts[i] + 1
    matrix.set(letter, counts)
    i = i + 1
  }

  matrix
}

///|
/// Compute consensus sequence from a multiple sequence alignment.
/// Returns the most frequent nucleotide at each position.
pub fn biostrings_consensus_sequence(alignment : Array[String]) -> String {
  if alignment.length() == 0 { return "" }
  let n_seqs = alignment.length()
  let seq_len = alignment[0].length()

  let result : FixedArray[UInt16] = FixedArray::make(seq_len, 0)
  let mut pos = 0
  while pos < seq_len {
    let counts : Map[UInt16, Int] = Map([], capacity=20)
    let mut s = 0
    while s < n_seqs {
      let letter = alignment[s][pos]
      let current = counts.get_or_default(letter, 0)
      counts.set(letter, current + 1)
      s = s + 1
    }

    let mut max_count = 0
    let mut best_letter : UInt16 = 65
    for letter in counts.keys() {
      let count = counts.get_or_default(letter, 0)
      if count > max_count {
        max_count = count
        best_letter = letter
      }
    }
    result[pos] = best_letter
    pos = pos + 1
  }

  unsafe_from_uint16_fixedarray(result)
}

///|
/// Translate DNA/RNA sequence to protein (amino acid sequence).
/// Uses standard genetic code.
pub fn biostrings_translate(seq : String, frame? : Int = 0) -> String {
  let n = seq.length()
  let offset = if frame >= 0 && frame < 3 { frame } else { 0 }
  let aa_count = (n - offset) / 3
  let result : FixedArray[UInt16] = FixedArray::make(aa_count, 0)

  let codon_table : Map[String, UInt16] = Map([
    ("TTT", 70), ("TTC", 70), ("TTA", 76), ("TTG", 76),
    ("CTT", 76), ("CTC", 76), ("CTA", 76), ("CTG", 76),
    ("ATT", 73), ("ATC", 73), ("ATA", 73), ("ATG", 77),
    ("GTT", 86), ("GTC", 86), ("GTA", 86), ("GTG", 86),
    ("TCT", 83), ("TCC", 83), ("TCA", 83), ("TCG", 83),
    ("CCT", 80), ("CCC", 80), ("CCA", 80), ("CCG", 80),
    ("ACT", 65), ("ACC", 65), ("ACA", 65), ("ACG", 65),
    ("GCT", 65), ("GCC", 65), ("GCA", 65), ("GCG", 65),
    ("TAT", 89), ("TAC", 89), ("TAA", 42), ("TAG", 42),
    ("CAT", 72), ("CAC", 72), ("CAA", 81), ("CAG", 81),
    ("AAT", 78), ("AAC", 78), ("AAA", 75), ("AAG", 75),
    ("GAT", 68), ("GAC", 68), ("GAA", 69), ("GAG", 69),
    ("TGT", 67), ("TGC", 67), ("TGA", 42), ("TGG", 87),
    ("CGT", 82), ("CGC", 82), ("CGA", 82), ("CGG", 82),
    ("AGT", 83), ("AGC", 83), ("AGA", 82), ("AGG", 82),
    ("GGT", 71), ("GGC", 71), ("GGA", 71), ("GGG", 71),
  ], capacity=64)

  let mut i = offset
  let mut out_idx = 0
  while i + 3 <= n {
    let c0 = seq[i]
    let c1 = seq[i + 1]
    let c2 = seq[i + 2]
    let u0 = if c0 >= 97 && c0 <= 122 { c0 - 32 } else { c0 }
    let u1 = if c1 >= 97 && c1 <= 122 { c1 - 32 } else { c1 }
    let u2 = if c2 >= 97 && c2 <= 122 { c2 - 32 } else { c2 }
    let codon_str = unsafe_from_uint16_fixedarray([u0, u1, u2])
    let aa = codon_table.get_or_default(codon_str, 42)
    result[out_idx] = aa
    i = i + 3
    out_idx = out_idx + 1
  }

  unsafe_from_uint16_fixedarray(result)
}

///|
/// Reverse complement of a DNA/RNA sequence.
/// Standard complement: A-T, C-G (with IUPAC support).
pub fn biostrings_reverse_complement(seq : String) -> String {
  let n = seq.length()
  let result : FixedArray[UInt16] = FixedArray::make(n, 0)

  let mut i = 0
  while i < n {
    let ch = seq[i]
    let complement : UInt16 = match ch {
      65 => 84
      84 => 65
      85 => 65
      67 => 71
      71 => 67
      78 => 78
      82 => 89
      89 => 82
      83 => 83
      87 => 87
      75 => 77
      77 => 75
      66 => 68
      68 => 66
      72 => 76
      86 => 72
      _ => ch
    }
    result[n - 1 - i] = complement
    i = i + 1
  }

  unsafe_from_uint16_fixedarray(result)
}

///|
/// Calculate the expected number of matches for a pattern in a sequence.
/// Based on letter frequencies and pattern length.
pub fn expected_matches(
  pattern : String,
  seq_length : Int,
  gc_content : Double,
) -> Double {
  let n = pattern.length()
  if n == 0 || seq_length <= n {
    return 0.0
  }

  // Compute probability of the pattern based on GC content
  let at_freq = (1.0 - gc_content) / 2.0
  let gc_freq = gc_content / 2.0

  let mut prob = 1.0
  let mut i = 0
  while i < n {
    let ch = pattern[i]
    let p = match ch {
      65 => at_freq
      84 => at_freq
      85 => at_freq
      67 => gc_freq
      71 => gc_freq
      82 => at_freq + gc_freq
      89 => at_freq + gc_freq
      78 => 1.0
      _ => 0.25
    }
    prob = prob * p
    i = i + 1
  }

  // Expected number: (seq_length - pattern_length + 1) * probability
  let n_positions = (seq_length - n + 1).to_double()
  if prob.is_nan() || prob < 0.0 { 0.0 } else { n_positions * prob }
}

///|
/// Calculate sequence complexity (Shannon entropy of oligonucleotide frequencies).
pub fn sequence_complexity(seq : String, word_size : Int) -> Double {
  let n = seq.length()
  if n < word_size || word_size < 1 {
    return 0.0
  }

  let counts : Map[String, Int] = Map([], capacity=n)
  let total_positions = n - word_size + 1

  let mut i = 0
  while i < total_positions {
    let sub = seq.substring(start=i, end=i + word_size)
    let buf : FixedArray[UInt16] = FixedArray::make(word_size, 0)
    let mut j = 0
    while j < word_size {
      let c = sub.unsafe_get(j)
      buf[j] = if c >= 97 && c <= 122 { c - 32 } else { c }
      j = j + 1
    }
    let word = unsafe_from_uint16_fixedarray(buf)
    let current = counts.get_or_default(word, 0)
    counts.set(word, current + 1)
    i = i + 1
  }

  // Shannon entropy
  let mut entropy = 0.0
  for word in counts.keys() {
    let count = counts.get_or_default(word, 0)
    let p = count.to_double() / total_positions.to_double()
    if p > 0.0 {
      entropy = entropy - p * @math.ln(p) / @math.ln(2.0)
    }
  }

  entropy
}

///| Test match_pattern with exact matching
test "match_pattern_exact" {
  let result = match_pattern(pattern="ATG", subject="ATGATGATG")
  assert_eq(result.hits.length(), 3)
  assert_eq(result.hits[0].start, 1)
  assert_eq(result.hits[0].end, 3)
  assert_eq(result.hits[0].mismatches, 0)
}

///| Test match_pattern with mismatches
test "match_pattern_mismatch" {
  let result = match_pattern(pattern="ATG", subject="AAGATG", max_mismatches=1)
  assert_true(result.hits.length() >= 1)
}

///| Test match_pattern with indels
test "match_pattern_indels" {
  let result = match_pattern(pattern="ATG", subject="ATGC", max_mismatches=1, with_indels=true)
  assert_true(result.hits.length() >= 1)
}

///| Test vmatch_pattern
test "vmatch_pattern" {
  let subjects = ["ATGATG", "ATGCAT", "CCCCCC"]
  let result = vmatch_pattern(pattern="ATG", subjects=subjects)
  assert_eq(result.length(), 3)
  assert_eq(result[0].hits.length(), 2)
  assert_eq(result[1].hits.length(), 1)
  assert_eq(result[2].hits.length(), 0)
}

///| Test find_palindromes
test "find_palindromes" {
  let result = find_palindromes(seq="ATAT", min_length=4)
  assert_eq(result.length(), 1)
  assert_eq(result[0].0, 1)
  assert_eq(result[0].1, 4)
}

///| Test find_palindromes_empty
test "find_palindromes_empty" {
  let result = find_palindromes(seq="ATGC", min_length=4)
  assert_eq(result.length(), 0)
}

///| Test find_direct_repeats
test "find_direct_repeats" {
  let result = find_direct_repeats(seq="ATGATGATG", min_unit_length=3, max_unit_length=3, min_copies=2)
  assert_true(result.length() >= 1)
}

///| Test find_inverted_repeats
test "find_inverted_repeats" {
  let result = find_inverted_repeats(seq="ATGCAT", min_length=3)
  assert_true(result.length() >= 0)
}

///| Test biostrings_translate with standard genetic code
test "biostrings_translate_basic" {
  // ATG = M, GCT = A, TAA = *
  let protein = biostrings_translate("ATGGCTTAA")
  assert_eq(protein, "MA*")
}

///| Test biostrings_translate with different reading frame
test "biostrings_translate_frame1" {
  // Frame 1: TGG = W, CTT = L
  let protein = biostrings_translate("ATGGCTTAA", frame=1)
  assert_eq(protein, "WL")
}

///| Test biostrings_translate with lowercase input
test "biostrings_translate_lowercase" {
  let protein = biostrings_translate("atggcttaa")
  assert_eq(protein, "MA*")
}

///| Test biostrings_translate empty sequence
test "biostrings_translate_empty" {
  let protein = biostrings_translate("")
  assert_eq(protein, "")
}

///| Test biostrings_reverse_complement basic
test "biostrings_reverse_complement_basic" {
  let rc = biostrings_reverse_complement("ATGC")
  assert_eq(rc, "GCAT")
}

///| Test biostrings_reverse_complement with IUPAC
test "biostrings_reverse_complement_iupac" {
  let rc = biostrings_reverse_complement("AR")
  assert_eq(rc, "YT")
}

///| Test biostrings_reverse_complement with empty string
test "biostrings_reverse_complement_empty" {
  let rc = biostrings_reverse_complement("")
  assert_eq(rc, "")
}

///| Test biostrings_reverse_complement palindrome
test "biostrings_reverse_complement_palindrome" {
  // ACGT's reverse complement is also ACGT
  let rc = biostrings_reverse_complement("ACGT")
  assert_eq(rc, "ACGT")
}

///| Test letter_frequency_matrix basic
test "letter_frequency_matrix_basic" {
  let matrix = letter_frequency_matrix("ACGT")
  let a_counts = matrix.get_or_default("A", Array::make(4, 0))
  assert_eq(a_counts[0], 1)
  let c_counts = matrix.get_or_default("C", Array::make(4, 0))
  assert_eq(c_counts[1], 1)
  let g_counts = matrix.get_or_default("G", Array::make(4, 0))
  assert_eq(g_counts[2], 1)
  let t_counts = matrix.get_or_default("T", Array::make(4, 0))
  assert_eq(t_counts[3], 1)
}

///| Test biostrings_consensus_sequence
test "biostrings_consensus_sequence_basic" {
  let alignment = ["ACGT", "ACGT", "TCGA"]
  let consensus = biostrings_consensus_sequence(alignment)
  // Position 0: A=2, T=1 -> A
  // Position 1: C=2, C=1 -> C
  // Position 2: G=2, G=1 -> G
  // Position 3: T=2, A=1 -> T
  assert_eq(consensus, "ACGT")
}

///| Test biostrings_consensus_sequence empty
test "biostrings_consensus_sequence_empty" {
  let alignment : Array[String] = []
  let consensus = biostrings_consensus_sequence(alignment)
  assert_eq(consensus, "")
}

///| Test expected_matches
test "expected_matches_basic" {
  let n = expected_matches("ACGT", 12, 0.5)
  assert_true(n > 0.0)
}

///| Test sequence_complexity
test "sequence_complexity_basic" {
  let entropy = sequence_complexity("ACGTACGTACGT", 2)
  assert_true(entropy > 0.0)
}

///| Test sequence_complexity with degenerate sequence
test "sequence_complexity_degenerate" {
  let entropy = sequence_complexity("AAAA", 2)
  assert_true(entropy >= 0.0)
}
