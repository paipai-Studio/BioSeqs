///|
/// Sequence utilities module.
let dna_molecular_weights : Map[String, Double] = Map([
  ("A", 313.2),
  ("T", 304.2),
  ("C", 289.2),
  ("G", 329.2),
  ("N", 0.0),
])

///|
let protein_molecular_weights : Map[String, Double] = Map([
  ("A", 89.1),
  ("R", 174.2),
  ("N", 132.1),
  ("D", 133.1),
  ("C", 121.2),
  ("E", 147.1),
  ("Q", 146.2),
  ("G", 75.1),
  ("H", 155.2),
  ("I", 131.2),
  ("L", 131.2),
  ("K", 146.2),
  ("M", 149.2),
  ("F", 165.2),
  ("P", 115.1),
  ("S", 105.1),
  ("T", 119.1),
  ("W", 204.2),
  ("Y", 181.2),
  ("V", 117.1),
])

///|
let protein_one_letter : Map[String, String] = Map([
  ("Ala", "A"),
  ("Arg", "R"),
  ("Asn", "N"),
  ("Asp", "D"),
  ("Cys", "C"),
  ("Glu", "E"),
  ("Gln", "Q"),
  ("Gly", "G"),
  ("His", "H"),
  ("Ile", "I"),
  ("Leu", "L"),
  ("Lys", "K"),
  ("Met", "M"),
  ("Phe", "F"),
  ("Pro", "P"),
  ("Ser", "S"),
  ("Thr", "T"),
  ("Trp", "W"),
  ("Tyr", "Y"),
  ("Val", "V"),
])

///|
let protein_three_letter : Map[String, String] = Map([
  ("A", "Ala"),
  ("R", "Arg"),
  ("N", "Asn"),
  ("D", "Asp"),
  ("C", "Cys"),
  ("E", "Glu"),
  ("Q", "Gln"),
  ("G", "Gly"),
  ("H", "His"),
  ("I", "Ile"),
  ("L", "Leu"),
  ("K", "Lys"),
  ("M", "Met"),
  ("F", "Phe"),
  ("P", "Pro"),
  ("S", "Ser"),
  ("T", "Thr"),
  ("W", "Trp"),
  ("Y", "Tyr"),
  ("V", "Val"),
])

///|
pub fn sequtils_molecular_weight(seq : String, seq_type : String) -> Double {
  let weights = match seq_type {
    "dna" => dna_molecular_weights
    "protein" => protein_molecular_weights
    _ => dna_molecular_weights
  }

  let mut weight = 0.0
  for i = 0; i < seq.length(); i = i + 1 {
    let char = seq[i:i + 1].to_owned()
    let added = match weights.get(char) {
      Some(w) => w
      None => 0.0
    }
    weight = weight + added
  }

  weight - (seq.length() - 1).to_double() * 18.0
}

///|
pub fn gc_content(seq : String) -> Double {
  if seq.length() == 0 {
    return 0.0
  }

  let mut gc = 0
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned().to_upper()
    if c == "G" || c == "C" {
      gc = gc + 1
    }
  }

  gc.to_double() / seq.length().to_double()
}

///|
pub fn melting_temperature(seq : String) -> Double {
  let len = seq.length()
  let gc = gc_content(seq)

  if len <= 14 {
    2.0 * (len.to_double() - gc * len.to_double()) + 4.0 * gc * len.to_double()
  } else {
    64.9 + 41.0 * (gc * len.to_double() - 16.4) / len.to_double()
  }
}

///|
pub fn seq1(seq : String) -> String {
  let mut result = ""
  let mut i = 0
  while i + 3 <= seq.length() {
    let triplet = seq[i:i + 3].to_owned()
    let aa = match protein_one_letter.get(triplet) {
      Some(a) => a
      None => "X"
    }
    result = result + aa
    i = i + 3
  }
  result
}

///|
pub fn seq3(seq : String) -> String {
  let mut result = ""
  for i = 0; i < seq.length(); i = i + 1 {
    let aa = seq[i:i + 1].to_owned()
    let triplet = match protein_three_letter.get(aa) {
      Some(t) => t
      None => "Xaa"
    }
    result = result + triplet + " "
  }
  if result.length() > 0 {
    let len = result.length()
    if result[len - 1:len] == " " {
      result[0:len - 1].to_owned()
    } else {
      result
    }
  } else {
    result
  }
}

///|
pub fn is_dna(seq : String) -> Bool {
  let allowed = ["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"]
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned()
    let mut is_allowed = false
    for j = 0; j < allowed.length() && !is_allowed; j = j + 1 {
      if c == allowed[j] {
        is_allowed = true
      }
    }
    if !is_allowed {
      return false
    }
  }
  true
}

///|
pub fn is_rna(seq : String) -> Bool {
  let allowed = ["A", "U", "C", "G", "N", "a", "u", "c", "g", "n"]
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned()
    let mut is_allowed = false
    for j = 0; j < allowed.length() && !is_allowed; j = j + 1 {
      if c == allowed[j] {
        is_allowed = true
      }
    }
    if !is_allowed {
      return false
    }
  }
  true
}

///|
pub fn is_protein(seq : String) -> Bool {
  let allowed = [
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S",
    "T", "W", "Y", "V", "X", "*",
  ]
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned().to_upper()
    let mut is_allowed = false
    for j = 0; j < allowed.length() && !is_allowed; j = j + 1 {
      if c == allowed[j] {
        is_allowed = true
      }
    }
    if !is_allowed {
      return false
    }
  }
  true
}

///|
pub fn count_kmers(seq : String, k : Int) -> Map[String, Int] {
  if k <= 0 || seq.length() < k {
    return Map([], capacity=0)
  }

  let total = seq.length() - k + 1
  let counts : Map[String, Int] = Map([], capacity=total)

  for i = 0; i <= seq.length() - k; i = i + 1 {
    let kmer = seq[i:i + k].to_owned()
    match counts.get(kmer) {
      Some(c) => counts[kmer] = c + 1
      None => counts[kmer] = 1
    }
  }
  counts
}

///|
pub fn sequtils_shannon_entropy(seq : String) -> Double {
  if seq.length() == 0 {
    return 0.0
  }

  let counts : Map[String, Int] = Map([], capacity=seq.length())
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned()
    match counts.get(c) {
      Some(count) => counts[c] = count + 1
      None => counts[c] = 1
    }
  }

  let len = seq.length().to_double()
  let keys = counts.keys().collect()
  let mut entropy = 0.0

  for idx = 0; idx < keys.length(); idx = idx + 1 {
    let key = keys[idx]
    let prob = counts[key].to_double() / len
    if prob > 0.0 {
      entropy = entropy - prob * sequtils_log_nat(prob)
    }
  }

  entropy
}

///|
fn sequtils_log_nat(x : Double) -> Double {
  if x <= 0.0 {
    return -1000000000000000000000000000000000.0
  }
  let term = (x - 1.0) / (x + 1.0)
  let term_squared = term * term

  let mut sum = term
  let mut current_term = term
  let mut i = 1
  while i <= 50 {
    current_term = current_term * term_squared
    sum = sum + current_term / (2 * i + 1).to_double()
    i = i + 1
  }

  2.0 * sum
}

///|
pub fn ambiguity_codes(seq : String) -> String {
  let mut result = ""
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i:i + 1].to_owned().to_upper()
    let expanded = match c {
      "A" => "A"
      "T" => "T"
      "C" => "C"
      "G" => "G"
      "N" => "ATCG"
      "R" => "AG"
      "Y" => "CT"
      "S" => "CG"
      "W" => "AT"
      "K" => "GT"
      "M" => "AC"
      "B" => "CGT"
      "D" => "AGT"
      "H" => "ACT"
      "V" => "ACG"
      _ => c
    }
    result = result + expanded
  }
  result
}

///|
pub fn codon_usage(seq : String) -> Map[String, Int] {
  if seq.length() < 3 {
    return Map([], capacity=0)
  }

  let total = seq.length() / 3
  let usage : Map[String, Int] = Map([], capacity=total)
  let mut i = 0
  while i + 3 <= seq.length() {
    let codon = seq[i:i + 3].to_owned()
    match usage.get(codon) {
      Some(c) => usage[codon] = c + 1
      None => usage[codon] = 1
    }
    i = i + 3
  }
  usage
}

///|
pub fn reverse_complement(seq : String) -> String {
  let mut result = ""
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i:i + 1].to_owned()
    let comp = match c.to_upper() {
      "A" => "T"
      "T" => "A"
      "C" => "G"
      "G" => "C"
      "U" => "A"
      _ => c
    }
    result = result + comp
  }
  result
}

///|
/// Calculate GC skew of a DNA sequence in a sliding window.
/// GC skew = (G - C) / (G + C)
pub fn gc_skew(seq : String, window_size : Int, step_size : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let len = seq.length()
  if len < window_size { return result }
  
  let mut i = 0
  while i <= len - window_size {
    let window = seq[i:i + window_size].to_owned().to_upper()
    let mut g_count = 0
    let mut c_count = 0
    let mut j = 0
    while j < window.length() {
      let c = window[j:j + 1]
      if c == "G" { g_count = g_count + 1 }
      if c == "C" { c_count = c_count + 1 }
      j = j + 1
    }
    let gc_total = g_count.to_double() + c_count.to_double()
    if gc_total > 0.0 {
      result.push((g_count.to_double() - c_count.to_double()) / gc_total)
    } else {
      result.push(0.0)
    }
    i = i + step_size
  }
  result
}

///|
/// Calculate AT skew of a DNA sequence in a sliding window.
/// AT skew = (A - T) / (A + T)
pub fn at_skew(seq : String, window_size : Int, step_size : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let len = seq.length()
  if len < window_size { return result }
  
  let mut i = 0
  while i <= len - window_size {
    let window = seq[i:i + window_size].to_owned().to_upper()
    let mut a_count = 0
    let mut t_count = 0
    let mut j = 0
    while j < window.length() {
      let c = window[j:j + 1]
      if c == "A" { a_count = a_count + 1 }
      if c == "T" { t_count = t_count + 1 }
      j = j + 1
    }
    let at_total = a_count.to_double() + t_count.to_double()
    if at_total > 0.0 {
      result.push((a_count.to_double() - t_count.to_double()) / at_total)
    } else {
      result.push(0.0)
    }
    i = i + step_size
  }
  result
}

///|
/// Calculate sequence similarity between two sequences.
pub fn sequence_similarity(seq1 : String, seq2 : String) -> Double {
  let len = min_val(seq1.length().to_double(), seq2.length().to_double())
  if len == 0.0 { return 0.0 }
  
  let mut matches = 0.0
  let mut i = 0
  while i < len.to_int() {
    if seq1[i:i + 1] == seq2[i:i + 1] {
      matches = matches + 1.0
    }
    i = i + 1
  }
  matches / len
}

///|
/// Calculate the edit distance (Levenshtein distance) between two sequences.
pub fn edit_distance(seq1 : String, seq2 : String) -> Int {
  let m = seq1.length()
  let n = seq2.length()
  
  if m == 0 { return n }
  if n == 0 { return m }
  
  // Dynamic programming approach
  let dp : Array[Array[Int]] = Array::make(m + 1, Array::make(n + 1, 0))
  
  let mut i = 0
  while i <= m {
    dp[i][0] = i
    i = i + 1
  }
  
  let mut j = 0
  while j <= n {
    dp[0][j] = j
    j = j + 1
  }
  
  i = 1
  while i <= m {
    j = 1
    while j <= n {
      if seq1[i - 1:i] == seq2[j - 1:j] {
        dp[i][j] = dp[i - 1][j - 1]
      } else {
        dp[i][j] = min_val3(
          dp[i - 1][j] + 1,
          dp[i][j - 1] + 1,
          dp[i - 1][j - 1] + 1,
        )
      }
      j = j + 1
    }
    i = i + 1
  }
  
  dp[m][n]
}

fn min_val(a : Double, b : Double) -> Double {
  if a < b { a } else { b }
}

fn min_val_int(a : Int, b : Int) -> Int {
  if a < b { a } else { b }
}

fn min_val3(a : Int, b : Int, c : Int) -> Int {
  min_val_int(min_val_int(a, b), c)
}
