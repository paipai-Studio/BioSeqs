///|
/// EMBOSS tools module inspired by Biopython Bio.Emboss.
/// Implements core sequence analysis algorithms from the EMBOSS suite.

///|
/// Calculate the GC skew of a sequence.
pub fn emboss_gc_skew(seq : String) -> Double {
  let s = seq.to_upper()
  let mut g = 0
  let mut c = 0
  let mut i = 0
  while i < s.length() {
    let ch = s.unsafe_get(i).unsafe_to_char()
    if ch == 'G' {
      g = g + 1
    } else if ch == 'C' {
      c = c + 1
    }
    i = i + 1
  }
  if g + c == 0 {
    0.0
  } else {
    (g.to_double() - c.to_double()) / (g.to_double() + c.to_double())
  }
}

///|
/// Calculate the AT skew of a sequence.
pub fn emboss_at_skew(seq : String) -> Double {
  let s = seq.to_upper()
  let mut a = 0
  let mut t = 0
  let mut i = 0
  while i < s.length() {
    let ch = s.unsafe_get(i).unsafe_to_char()
    if ch == 'A' {
      a = a + 1
    } else if ch == 'T' {
      t = t + 1
    }
    i = i + 1
  }
  if a + t == 0 {
    0.0
  } else {
    (a.to_double() - t.to_double()) / (a.to_double() + t.to_double())
  }
}

///|
/// Calculate the molecular weight of a DNA/RNA sequence.
pub fn dna_molecular_weight(seq : String) -> Double {
  let s = seq.to_upper()
  let mut weight = 0.0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    let w = if c == 'A' {
      313.21
    } else if c == 'T' {
      304.2
    } else if c == 'C' {
      289.18
    } else if c == 'G' {
      329.21
    } else if c == 'U' {
      290.17
    } else {
      300.0
    }
    weight = weight + w
    i = i + 1
  }
  weight + 18.0
}

///|
/// Calculate the melting temperature (Tm) of a DNA sequence using the nearest neighbor method.
pub fn tm_nearest_neighbor(seq : String) -> Double {
  let s = seq.to_upper()
  let mut sum = 0.0
  let mut gc_count = 0
  let n = s.length()
  if n < 2 {
    return 0.0
  }

  let mut i = 0
  while i < n - 1 {
    let dinuc = substring(s, i, 2)
    let score = if dinuc == "AA" || dinuc == "TT" {
      -9.1
    } else if dinuc == "AT" {
      -8.6
    } else if dinuc == "TA" {
      -6.0
    } else if dinuc == "CA" || dinuc == "TG" {
      -5.8
    } else if dinuc == "GT" || dinuc == "AC" {
      -6.5
    } else if dinuc == "CT" || dinuc == "AG" {
      -7.8
    } else if dinuc == "GA" || dinuc == "TC" {
      -5.6
    } else if dinuc == "CG" {
      -11.9
    } else if dinuc == "GC" {
      -11.1
    } else if dinuc == "GG" || dinuc == "CC" {
      -11.0
    } else {
      -8.0
    }
    sum = sum + score
    if dinuc.unsafe_get(0).unsafe_to_char() == 'G' ||
      dinuc.unsafe_get(0).unsafe_to_char() == 'C' ||
      dinuc.unsafe_get(1).unsafe_to_char() == 'G' ||
      dinuc.unsafe_get(1).unsafe_to_char() == 'C' {
      gc_count = gc_count + 1
    }
    i = i + 1
  }

  let tm = 100.5 +
    41.0 * gc_count.to_double() / n.to_double() -
    820.0 / n.to_double() +
    sum
  tm
}

///|
/// Calculate the melting temperature using the simplified formula.
pub fn tm_simple(seq : String) -> Double {
  let s = seq.to_upper()
  let mut gc = 0
  let mut at = 0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c == 'G' || c == 'C' {
      gc = gc + 1
    } else if c == 'A' || c == 'T' {
      at = at + 1
    }
    i = i + 1
  }
  2.0 * at.to_double() + 4.0 * gc.to_double()
}

///|
/// Reverse complement of a DNA sequence.
pub fn emboss_reverse_complement(seq : String) -> String {
  let s = seq.to_upper()
  let n = s.length()
  let mut result = ""
  let mut i = n - 1
  while i >= 0 {
    let c = s.unsafe_get(i).unsafe_to_char()
    let comp = if c == 'A' {
      "T"
    } else if c == 'T' || c == 'U' {
      "A"
    } else if c == 'G' {
      "C"
    } else if c == 'C' {
      "G"
    } else if c == 'N' || c == 'X' {
      c.to_string()
    } else {
      "N"
    }
    result = result + comp
    i = i - 1
  }
  result
}

///|
/// Translate DNA to protein using the standard genetic code.
pub fn translate(seq : String) -> String {
  let s = seq.to_upper()
  let codon_table = codon_table_std()
  let mut result = ""
  let mut i = 0
  while i < s.length() - 2 {
    let codon = substring(s, i, 3)
    let aa_opt = codon_table.get(codon)
    let aa = if aa_opt is Some(_) {
      match aa_opt {
        Some(a) => a
        None => "X"
      }
    } else {
      "X"
    }
    result = result + aa
    i = i + 3
  }
  result
}

///|
/// Get the standard codon table.
fn codon_table_std() -> Map[String, String] {
  let table = Map([], capacity=64)
  table["TTT"] = "F"
  table["TTC"] = "F"
  table["TTA"] = "L"
  table["TTG"] = "L"
  table["CTT"] = "L"
  table["CTC"] = "L"
  table["CTA"] = "L"
  table["CTG"] = "L"
  table["ATT"] = "I"
  table["ATC"] = "I"
  table["ATA"] = "I"
  table["ATG"] = "M"
  table["GTT"] = "V"
  table["GTC"] = "V"
  table["GTA"] = "V"
  table["GTG"] = "V"
  table["TCT"] = "S"
  table["TCC"] = "S"
  table["TCA"] = "S"
  table["TCG"] = "S"
  table["CCT"] = "P"
  table["CCC"] = "P"
  table["CCA"] = "P"
  table["CCG"] = "P"
  table["ACT"] = "T"
  table["ACC"] = "T"
  table["ACA"] = "T"
  table["ACG"] = "T"
  table["GCT"] = "A"
  table["GCC"] = "A"
  table["GCA"] = "A"
  table["GCG"] = "A"
  table["TAT"] = "Y"
  table["TAC"] = "Y"
  table["TAA"] = "*"
  table["TAG"] = "*"
  table["CAT"] = "H"
  table["CAC"] = "H"
  table["CAA"] = "Q"
  table["CAG"] = "Q"
  table["AAT"] = "N"
  table["AAC"] = "N"
  table["AAA"] = "K"
  table["AAG"] = "K"
  table["GAT"] = "D"
  table["GAC"] = "D"
  table["GAA"] = "E"
  table["GAG"] = "E"
  table["TGT"] = "C"
  table["TGC"] = "C"
  table["TGA"] = "*"
  table["TGG"] = "W"
  table["CGT"] = "R"
  table["CGC"] = "R"
  table["CGA"] = "R"
  table["CGG"] = "R"
  table["AGT"] = "S"
  table["AGC"] = "S"
  table["AGA"] = "R"
  table["AGG"] = "R"
  table["GGT"] = "G"
  table["GGC"] = "G"
  table["GGA"] = "G"
  table["GGG"] = "G"
  table
}

///|
/// Find open reading frames (ORFs) in a DNA sequence.
pub fn find_orfs(
  seq : String,
  min_length? : Int = 30,
) -> Array[(Int, Int, String)] {
  let s = seq.to_upper()
  let result : Array[(Int, Int, String)] = Array::new()
  let start_codons = ["ATG"]
  let stop_codons = ["TAA", "TAG", "TGA"]
  let n = s.length()

  for frame in [0, 1, 2] {
    let mut i = frame
    while i < n - 2 {
      let codon = substring(s, i, 3)
      let is_start = {
        let mut found = false
        for sc in start_codons {
          if codon == sc {
            found = true
            break
          }
        }
        found
      }
      if is_start {
        let mut j = i + 3
        while j < n - 2 {
          let stop_codon = substring(s, j, 3)
          let is_stop = {
            let mut found = false
            for stc in stop_codons {
              if stop_codon == stc {
                found = true
                break
              }
            }
            found
          }
          if is_stop {
            let length = j - i + 3
            if length >= min_length {
              let orf_seq = substring(s, i, length)
              let protein = translate(orf_seq)
              result.push((i, j + 3, protein))
            }
            break
          }
          j = j + 3
        }
        i = j
      } else {
        i = i + 3
      }
    }
  }

  result
}

///|
/// Calculate the Hamming distance between two sequences.
pub fn hamming_distance(seq1 : String, seq2 : String) -> Int {
  let s1 = seq1.to_upper()
  let s2 = seq2.to_upper()
  let mut dist = 0
  let min_len = if s1.length() < s2.length() {
    s1.length()
  } else {
    s2.length()
  }
  let mut i = 0
  while i < min_len {
    if s1.unsafe_get(i).unsafe_to_char() != s2.unsafe_get(i).unsafe_to_char() {
      dist = dist + 1
    }
    i = i + 1
  }
  dist + (s1.length() - min_len).abs()
}

///|
/// Calculate the Levenshtein distance between two sequences.
pub fn levenshtein_distance(seq1 : String, seq2 : String) -> Int {
  let s1 = seq1.to_upper()
  let s2 = seq2.to_upper()
  let m = s1.length()
  let n = s2.length()

  let matrix : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i <= m {
    let row : Array[Int] = Array::new()
    let mut j = 0
    while j <= n {
      row.push(0)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }

  i = 0
  while i <= m {
    matrix[i][0] = i
    i = i + 1
  }
  let mut j = 0
  while j <= n {
    matrix[0][j] = j
    j = j + 1
  }

  i = 1
  while i <= m {
    j = 1
    while j <= n {
      let cost = if s1.unsafe_get(i - 1).unsafe_to_char() ==
        s2.unsafe_get(j - 1).unsafe_to_char() {
        0
      } else {
        1
      }
      let del = matrix[i - 1][j] + 1
      let ins = matrix[i][j - 1] + 1
      let sub = matrix[i - 1][j - 1] + cost
      matrix[i][j] = min3(del, ins, sub)
      j = j + 1
    }
    i = i + 1
  }

  matrix[m][n]
}

///|
/// Helper function: minimum of three integers.
fn min3(a : Int, b : Int, c : Int) -> Int {
  let m = if a < b { a } else { b }
  if m < c {
    m
  } else {
    c
  }
}

///|
/// Count k-mers in a sequence.
pub fn emboss_count_kmers(seq : String, k : Int) -> Map[String, Int] {
  let s = seq.to_upper()
  let result = Map([], capacity=128)
  let n = s.length()
  if k <= 0 || k > n {
    return result
  }

  let mut i = 0
  while i <= n - k {
    let kmer = substring(s, i, k)
    let valid = {
      let mut v = true
      let mut j = 0
      while j < kmer.length() {
        let c = kmer.unsafe_get(j).unsafe_to_char()
        if c != 'A' &&
          c != 'T' &&
          c != 'C' &&
          c != 'G' &&
          c != 'N' &&
          c != 'U' &&
          c != 'R' &&
          c != 'Y' &&
          c != 'S' &&
          c != 'W' &&
          c != 'K' &&
          c != 'M' &&
          c != 'B' &&
          c != 'D' &&
          c != 'H' &&
          c != 'V' {
          v = false
          break
        }
        j = j + 1
      }
      v
    }
    if valid {
      match result.get(kmer) {
        Some(count) => result[kmer] = count + 1
        None => result[kmer] = 1
      }
    }
    i = i + 1
  }

  result
}

///|
/// Calculate the isoelectric point (pI) of a protein sequence.
pub fn isoelectric_point(seq : String) -> Double {
  let s = seq.to_upper()
  let mut pos_charged = 0
  let mut neg_charged = 0
  let mut histidine = 0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    match c {
      'K' => pos_charged = pos_charged + 1
      'R' => pos_charged = pos_charged + 1
      'H' => histidine = histidine + 1
      'D' => neg_charged = neg_charged + 1
      'E' => neg_charged = neg_charged + 1
      _ => ()
    }
    i = i + 1
  }

  let total = s.length()
  if total == 0 {
    return 7.0
  }

  let pos_ratio = pos_charged.to_double() / total.to_double()
  let neg_ratio = neg_charged.to_double() / total.to_double()

  if pos_ratio > neg_ratio {
    7.0 + (pos_ratio - neg_ratio) * 5.0
  } else if neg_ratio > pos_ratio {
    7.0 - (neg_ratio - pos_ratio) * 5.0
  } else {
    7.0
  }
}

///|
/// Calculate the aliphatic index of a protein sequence.
pub fn aliphatic_index(seq : String) -> Double {
  let s = seq.to_upper()
  let mut a = 0
  let mut v = 0
  let mut i = 0
  let mut l = 0
  let mut total = 0
  let mut idx = 0
  while idx < s.length() {
    let c = s.unsafe_get(idx).unsafe_to_char()
    match c {
      'A' => a = a + 1
      'V' => v = v + 1
      'I' => i = i + 1
      'L' => l = l + 1
      _ => ()
    }
    total = total + 1
    idx = idx + 1
  }
  if total == 0 {
    0.0
  } else {
    (
      a.to_double() * 100.0 +
      v.to_double() * 86.0 +
      i.to_double() * 108.0 +
      l.to_double() * 197.0
    ) /
    total.to_double()
  }
}

///|
/// Calculate the grand average of hydropathicity (GRAVY) score.
pub fn gravy_score(seq : String) -> Double {
  let s = seq.to_upper()
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    let score = if c == 'A' {
      1.8
    } else if c == 'R' {
      -4.5
    } else if c == 'N' || c == 'D' || c == 'E' || c == 'Q' {
      -3.5
    } else if c == 'C' {
      2.5
    } else if c == 'G' {
      -0.4
    } else if c == 'H' {
      -3.2
    } else if c == 'I' {
      4.5
    } else if c == 'L' {
      3.8
    } else if c == 'K' {
      -3.9
    } else if c == 'M' {
      1.9
    } else if c == 'F' {
      2.8
    } else if c == 'P' {
      -1.6
    } else if c == 'S' {
      -0.8
    } else if c == 'T' {
      -0.7
    } else if c == 'W' {
      -0.9
    } else if c == 'Y' {
      -1.3
    } else if c == 'V' {
      4.2
    } else {
      0.0
    }
    sum = sum + score
    count = count + 1
    i = i + 1
  }
  if count == 0 {
    0.0
  } else {
    sum / count.to_double()
  }
}

///|
/// Calculate the instability index of a protein sequence.
pub fn instability_index(seq : String) -> Double {
  let s = seq.to_upper()
  let n = s.length()
  if n < 2 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < n - 1 {
    let aa1 = s.unsafe_get(i).unsafe_to_char()
    let aa2 = s.unsafe_get(i + 1).unsafe_to_char()
    sum = sum + instability_pair_score(aa1, aa2)
    i = i + 1
  }
  sum / (n - 1).to_double() * 10.0
}

///|
/// Get the instability score for an amino acid pair.
fn instability_pair_score(a1 : Char, a2 : Char) -> Double {
  let aa = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
    'T', 'W', 'Y', 'V',
  ]
  let scores : Array[Array[Double]] = Array::new()
  scores.push([
    4.9, -1.8, -4.8, -3.5, 0.2, -8.0, -4.1, 0.0, -4.6, 4.2, 4.9, -3.9, 1.9, 2.8,
    -4.6, 0.4, 0.5, -0.9, -1.3, 4.2,
  ])
  scores.push([
    -1.8, 10.5, -2.8, -4.0, -5.5, -3.5, -1.0, -0.8, -3.2, 5.2, 5.9, 10.8, -2.5, 2.5,
    -4.0, 0.4, -1.5, -1.6, -2.8, -4.8,
  ])
  scores.push([
    -1.6, -6.1, -5.1, -1.0, -5.6, -3.0, -4.6, 0.2, -0.6, -4.7, -0.2, -1.8, 2.0, -3.5,
    -4.0, -6.2, -4.0, -0.8, -2.0, 0.3,
  ])
  scores.push([
    -3.6, -3.5, -6.1, -2.4, -4.0, -5.5, -0.8, -3.0, -4.1, -0.5, -5.5, -5.5, -2.8,
    -1.0, 0.4, -4.0, 0.8, -6.0, -1.0, -4.5,
  ])
  scores.push([
    -8.0, -3.5, -6.2, -7.0, -7.0, -0.8, -1.0, -5.0, -3.0, -3.5, -5.5, -5.0, -6.0,
    -4.5, -5.0, -4.1, -1.0, -6.2, -4.0, -5.4,
  ])
  scores.push([
    -6.0, -6.0, -2.5, -2.0, -3.0, -5.0, -4.0, -4.5, -5.5, -4.0, -5.0, 0.0, -0.8,
    -2.8, -3.0, -6.1, -5.9, -0.4, -1.8, -2.0,
  ])
  scores.push([
    -4.0, -3.0, -3.5, -5.0, 0.0, -4.6, -3.2, -5.4, -3.0, -2.5, -2.5, -3.2, -1.5,
    -1.0, -3.0, -4.0, -5.0, -4.5, -6.0, -4.2,
  ])
  scores.push([
    5.2, 4.2, -0.6, -2.5, 4.3, 4.9, 5.7, 2.5, 3.8, -1.6, -0.6, 0.0, -0.9, 3.8, 4.2,
    4.9, 5.9, 4.9, -1.0, -2.5,
  ])
  scores.push([
    4.9, 4.9, 5.5, 3.8, 3.8, -1.6, -1.0, 0.0, -0.9, 3.8, 4.2, -3.9, 10.8, 5.5, -3.5,
    -3.2, 5.7, 5.5, 10.0, 2.8,
  ])
  scores.push([
    3.0, -4.5, -1.5, -2.0, -3.0, 2.5, 3.0, 1.9, -2.5, 3.8, -1.8, -1.5, 2.5, 3.8,
    2.8, 5.2, 3.8, 3.8, -2.0, -1.8,
  ])
  scores.push([
    -1.0, 1.9, 2.8, 4.0, 2.8, 2.5, 3.8, -2.0, -1.0, 3.8, 3.8, 3.0, 3.8, 6.0, -1.0,
    -1.0, 0.0, 2.0, -4.6, -4.0,
  ])
  scores.push([
    -4.0, -5.5, -5.0, -4.0, -3.0, -1.6, -1.6, -4.5, -2.0, -1.0, -4.0, -3.5, -4.0,
    -2.0, -3.0, 0.4, 0.2, -0.6, 1.7, -1.5,
  ])
  scores.push([
    -4.0, -5.5, -5.0, -4.0, -3.0, -1.6, -1.6, -4.5, -2.0, -1.0, -4.0, -3.5, -4.0,
    -2.0, -3.0, 0.4, 0.2, -0.6, 1.7, -1.5,
  ])
  scores.push([
    -4.5, -2.5, -3.5, -1.0, -4.5, -2.8, -1.0, 3.8, 5.5, 3.8, 6.0, -2.0, -2.0, 6.0,
    -2.5, -1.0, -1.0, 4.9, 3.7, -0.5,
  ])
  scores.push([
    -4.6, -4.0, -4.0, 0.4, -5.0, -3.0, -3.0, 4.2, -3.5, 2.8, -1.0, -3.0, -3.0, -2.5,
    -1.6, 4.6, -0.2, -1.8, -3.0, -1.6,
  ])
  scores.push([
    0.4, 0.4, -6.2, -4.0, -4.1, -6.1, -4.0, 4.9, -3.2, 5.2, -1.0, 0.4, 0.4, -1.0,
    4.6, 0.8, 0.5, -0.7, -0.2, 5.5,
  ])
  scores.push([
    0.5, -1.5, -4.0, 0.8, -1.0, -5.9, -5.0, 5.9, 5.7, 3.8, 0.0, 0.2, 0.2, -1.0, -0.2,
    0.5, 0.3, -1.1, 0.4, 4.3,
  ])
  scores.push([
    -0.9, -1.6, -0.8, -6.0, -6.2, -0.4, -4.5, 4.9, 5.5, 3.8, 2.0, -0.6, -0.6, 4.9,
    -1.8, -0.7, -1.1, -0.9, 1.9, 3.6,
  ])
  scores.push([
    -1.3, -2.8, -2.0, -1.0, -4.0, -1.8, -6.0, -1.0, 10.0, -2.0, -4.6, 1.7, 1.7, 3.7,
    -3.0, -0.2, 0.4, 1.9, -0.7, -1.0,
  ])
  scores.push([
    4.2, -4.8, 0.3, -4.5, -5.4, -2.0, -4.2, -2.5, 2.8, -1.8, -4.0, -1.5, -1.5, -0.5,
    -1.6, 5.5, 4.3, 3.6, -1.0, 4.2,
  ])

  let mut idx1 = -1
  let mut j = 0
  while j < aa.length() {
    if aa[j] == a1 {
      idx1 = j
      break
    }
    j = j + 1
  }
  let mut idx2 = -1
  j = 0
  while j < aa.length() {
    if aa[j] == a2 {
      idx2 = j
      break
    }
    j = j + 1
  }
  if idx1 >= 0 &&
    idx2 >= 0 &&
    idx1 < scores.length() &&
    idx2 < scores[0].length() {
    scores[idx1][idx2]
  } else {
    0.0
  }
}

///|
/// Calculate the amino acid composition of a protein sequence.
pub fn emboss_amino_acid_composition(seq : String) -> Map[String, Double] {
  let s = seq.to_upper()
  let result = Map([], capacity=20)
  let aa_list = [
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S",
    "T", "W", "Y", "V",
  ]
  for aa in aa_list {
    result[aa] = 0.0
  }
  let mut total = 0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    let aa_str = c.to_string()
    match result.get(aa_str) {
      Some(_) => {
        result[aa_str] = result[aa_str] + 1.0
        total = total + 1
      }
      None => ()
    }
    i = i + 1
  }
  if total > 0 {
    for aa in aa_list {
      result[aa] = result[aa] / total.to_double() * 100.0
    }
  }
  result
}

///|
/// Sample DNA sequence for testing.
pub fn sample_dna_sequence() -> String {
  "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA"
}

///|
/// Sample protein sequence for testing.
pub fn sample_protein_sequence() -> String {
  "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
}
