///|
/// biotite_align_utils3.mbt — Biotite round 3: alignment utility functions
/// 7 features: KmerAlphabet, SimilarityRule/ScoreThresholdRule, MincodeSelector,
/// Permutation/FrequencyPermutation, pairwise identity matrix, gap utilities,
/// alignment extraction functions.

// #7 KmerAlphabet — k-mers as alphabet symbols
pub struct Biotite3KmerAlphabet {
  base_alphabet : Array[String]
  k : Int32
  kmer_symbols : Array[String]
}

///|
pub fn biotite3_kmer_alphabet(
  base : Array[String],
  k : Int32,
) -> Biotite3KmerAlphabet {
  let kmer_symbols = Array(capacity=0)
  biotite3_kmer_alphabet_generate(base, k, "", kmer_symbols)
  { base_alphabet: base, k, kmer_symbols, }
}

///|
fn biotite3_kmer_alphabet_generate(
  base : Array[String],
  k : Int32,
  prefix : String,
  result : Array[String],
) -> Unit {
  if k == 0 {
    result.push(prefix)
    return
  }
  for i in 0..<base.length() {
    biotite3_kmer_alphabet_generate(base, k - 1, prefix + base[i], result)
  }
}

///|
pub fn biotite3_kmer_alphabet_encode(
  kmer_alpha : Biotite3KmerAlphabet,
  symbols : Array[String],
) -> Int32 {
  let kmer = {
    let sb = StringBuilder()
    for i in 0..<symbols.length() {
      let _ = sb.write(symbols[i])
    }
    sb.to_string()
  }
  let mut idx = -1
  for i in 0..<kmer_alpha.kmer_symbols.length() {
    if kmer_alpha.kmer_symbols[i] == kmer {
      idx = i
      break
    }
  }
  idx
}

///|
pub fn biotite3_kmer_alphabet_decode(
  kmer_alpha : Biotite3KmerAlphabet,
  code : Int32,
) -> Array[String] {
  if code < 0 || code >= kmer_alpha.kmer_symbols.length() {
    return []
  }
  let kmer = kmer_alpha.kmer_symbols[code]
  let result = Array(capacity=kmer_alpha.k)
  for i in 0..<kmer.length() {
    result.push(kmer[i:i + 1].to_owned())
  }
  result
}

// #8 SimilarityRule / ScoreThresholdRule

///|
pub struct Biotite3ScoreThresholdRule {
  threshold : Double
  alphabet : Array[String]
  scores : Array[Array[Double]]
}

///|
pub fn biotite3_score_threshold_rule(
  alphabet : Array[String],
  scores : Array[Array[Double]],
  threshold : Double,
) -> Biotite3ScoreThresholdRule {
  { threshold, alphabet, scores, }
}

///|
pub fn biotite3_score_threshold_find_similar(
  rule : Biotite3ScoreThresholdRule,
  query_code : Int32,
) -> Array[Int32] {
  let result = Array(capacity=0)
  if query_code < 0 || query_code >= rule.alphabet.length() {
    return result
  }
  for i in 0..<rule.alphabet.length() {
    if i != query_code &&
      query_code < rule.scores.length() &&
      i < rule.scores[query_code].length() {
      if rule.scores[query_code][i] >= rule.threshold {
        result.push(i)
      }
    }
  }
  result
}

// #9 MincodeSelector — select 1/compression smallest k-mers

///|
pub fn biotite3_mincode_select(
  symbols : Array[String],
  k : Int32,
  n_select : Int32,
) -> Array[Int32] {
  let n_kmers = symbols.length() - k + 1
  if n_kmers <= 0 {
    return []
  }
  let kmer_hashes = Array(capacity=n_kmers)
  for i in 0..<n_kmers {
    let mut h = 0
    for j in 0..<k {
      let c = symbols[i + j]
      for c_idx in 0..<c.length() {
        h = h * 31 + c[c_idx].to_int()
      }
    }
    kmer_hashes.push((h, i))
  }
  kmer_hashes.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  let result = Array(capacity=n_select)
  let count = if n_select < kmer_hashes.length() {
    n_select
  } else {
    kmer_hashes.length()
  }
  for i in 0..<count {
    result.push(kmer_hashes[i].1)
  }
  result.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  result
}

// #10 Permutation / FrequencyPermutation — k-mer ordering

///|
pub struct Biotite3FrequencyPermutation {
  kmer_frequencies : Array[(String, Int32)]
}

///|
pub fn biotite3_frequency_permutation(
  kmers : Array[String],
) -> Biotite3FrequencyPermutation {
  let freq : Map[String, Int32] = Map([])
  for i in 0..<kmers.length() {
    freq[kmers[i]] = freq.get(kmers[i]).unwrap_or(0) + 1
  }
  let pairs = Array(capacity=0)
  for k in freq.keys() {
    pairs.push((k, freq[k]))
  }
  pairs.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  { kmer_frequencies: pairs, }
}

///|
pub fn biotite3_frequency_permutation_rank(
  perm : Biotite3FrequencyPermutation,
  kmer : String,
) -> Int32 {
  let mut rank = -1
  for i in 0..<perm.kmer_frequencies.length() {
    if perm.kmer_frequencies[i].0 == kmer {
      rank = i
      break
    }
  }
  rank
}

///|
pub fn biotite3_random_permutation(
  kmers : Array[String],
  seed : Int32,
) -> Array[String] {
  let result = Array(capacity=kmers.length())
  for i in 0..<kmers.length() {
    result.push(kmers[i])
  }
  let mut s = seed
  let _upper_i = result.length() - 1
  for i in 0..<_upper_i {
    s = (s * 1103515245 + 12345) % 2147483647
    let j = i + s % (result.length() - i)
    let tmp = result[i]
    result[i] = result[j]
    result[j] = tmp
  }
  result
}

// #11 get_pairwise_sequence_identity — identity matrix for MSA

///|
pub fn biotite3_pairwise_identity_matrix(
  sequences : Array[String],
) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      if i == j {
        row.push(1.0)
      } else {
        let s1 = sequences[i]
        let s2 = sequences[j]
        let min_len = if s1.length() < s2.length() {
          s1.length()
        } else {
          s2.length()
        }
        let mut matches = 0
        for k in 0..<min_len {
          if s1[k] == s2[k] {
            matches = matches + 1
          }
        }
        row.push(matches.to_double() / min_len.to_double())
      }
    }
    matrix.push(row)
  }
  matrix
}

// #12 find_terminal_gaps / remove_gaps / remove_terminal_gaps

///|
pub fn biotite3_find_terminal_gaps(
  aligned : Array[String],
  gap_char : String,
) -> (Int32, Int32) {
  let mut left = 0
  let mut right = 0
  if aligned.length() == 0 {
    return (0, 0)
  }
  let seq_len = aligned[0].length()
  let mut still_gap = true
  for col in 0..<seq_len {
    if !still_gap {
      break
    }
    let mut all_gap = true
    for row in 0..<aligned.length() {
      if col < aligned[row].length() {
        let ch = aligned[row][col:col + 1].to_owned()
        if ch != gap_char {
          all_gap = false
          break
        }
      }
    }
    if all_gap {
      left = left + 1
    } else {
      still_gap = false
    }
  }
  still_gap = true
  for col_idx in 0..<seq_len {
    let col = seq_len - 1 - col_idx
    if !still_gap {
      break
    }
    let mut all_gap = true
    for row in 0..<aligned.length() {
      if col < aligned[row].length() {
        let ch = aligned[row][col:col + 1].to_owned()
        if ch != gap_char {
          all_gap = false
          break
        }
      }
    }
    if all_gap {
      right = right + 1
    } else {
      still_gap = false
    }
  }
  (left, right)
}

///|
pub fn biotite3_remove_terminal_gaps(
  aligned : Array[String],
  gap_char : String,
) -> Array[String] {
  let (left, right) = biotite3_find_terminal_gaps(aligned, gap_char)
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  let result = Array(capacity=aligned.length())
  for i in 0..<aligned.length() {
    let end = seq_len - right
    if left < end && end <= aligned[i].length() {
      result.push(aligned[i][left:end].to_owned())
    } else {
      result.push("")
    }
  }
  result
}

///|
pub fn biotite3_remove_all_gaps(
  aligned : Array[String],
  gap_char : String,
) -> Array[String] {
  let result = Array(capacity=aligned.length())
  for i in 0..<aligned.length() {
    let sb = StringBuilder()
    for j in 0..<aligned[i].length() {
      let ch = aligned[i][j:j + 1].to_owned()
      if ch != gap_char {
        let _ = sb.write(ch)
      }
    }
    result.push(sb.to_string())
  }
  result
}

// #13 score / get_codes / get_symbols from alignment

///|
pub fn biotite3_alignment_score(
  aligned : Array[String],
  matrix : Array[Array[Double]],
  alphabet : Array[String],
) -> Double {
  let mut total = 0.0
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  for col in 0..<seq_len {
    for i in 0..<aligned.length() {
      let _start_j = i + 1
      for j in _start_j..<aligned.length() {
        if col < aligned[i].length() && col < aligned[j].length() {
          let c1 = aligned[i][col:col + 1].to_owned()
          let c2 = aligned[j][col:col + 1].to_owned()
          let mut idx1 = -1
          let mut idx2 = -1
          for a in 0..<alphabet.length() {
            if alphabet[a] == c1 {
              idx1 = a
            }
            if alphabet[a] == c2 {
              idx2 = a
            }
          }
          if idx1 >= 0 &&
            idx2 >= 0 &&
            idx1 < matrix.length() &&
            idx2 < matrix[idx1].length() {
            total = total + matrix[idx1][idx2]
          }
        }
      }
    }
  }
  total
}

///|
pub fn biotite3_alignment_get_codes(
  aligned : Array[String],
  alphabet : Array[String],
) -> Array[Array[Int32]] {
  let result = Array(capacity=aligned.length())
  for i in 0..<aligned.length() {
    let codes = Array(capacity=aligned[i].length())
    for j in 0..<aligned[i].length() {
      let ch = aligned[i][j:j + 1].to_owned()
      let mut idx = -1
      for a in 0..<alphabet.length() {
        if alphabet[a] == ch {
          idx = a
          break
        }
      }
      codes.push(idx)
    }
    result.push(codes)
  }
  result
}

///|
pub fn biotite3_alignment_get_symbols(
  aligned : Array[String],
) -> Array[Array[String]] {
  let result = Array(capacity=aligned.length())
  for i in 0..<aligned.length() {
    let symbols = Array(capacity=aligned[i].length())
    for j in 0..<aligned[i].length() {
      symbols.push(aligned[i][j:j + 1].to_owned())
    }
    result.push(symbols)
  }
  result
}

///|
pub fn _biotite_align_utils3_dummy() -> Unit {
  let _ = Biotite3KmerAlphabet::{ base_alphabet: [], k: 0, kmer_symbols: [], }
  let _ = Biotite3ScoreThresholdRule::{
    threshold: 0.0,
    alphabet: [],
    scores: [],
  }
  let _ = Biotite3FrequencyPermutation::{ kmer_frequencies: [], }
}
