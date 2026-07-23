///|
/// Sequence motifs analysis module.
pub struct Motif {
  name : String
  pwm : Array[Array[Double]]
  alphabet : Array[String]
  sites : Array[String]
}

///|
pub struct MotifHit {
  position : Int
  score : Double
  strand : String
  sequence : String
}

///|
pub type PSSM = Array[Array[Double]]

///|
let dna_alphabet : Array[String] = ["A", "C", "G", "T"]

///|
let protein_alphabet : Array[String] = [
  "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S",
  "T", "W", "Y", "V",
]

///|
pub fn Motif::new(
  name : String,
  pwm : Array[Array[Double]],
  alphabet : Array[String],
  sites : Array[String],
) -> Motif {
  Motif::{ name, pwm, alphabet, sites }
}

///|
pub fn MotifHit::new(
  position : Int,
  score : Double,
  strand : String,
  sequence : String,
) -> MotifHit {
  MotifHit::{ position, score, strand, sequence }
}

///|
pub fn count_matrix(sites : Array[String]) -> Array[Array[Int]] {
  if sites.length() == 0 {
    return []
  }
  let width = sites[0].length()
  let alphabet_size = dna_alphabet.length()
  let counts = Array::make(width, Array::make(alphabet_size, 0))

  for idx = 0; idx < sites.length(); idx = idx + 1 {
    let site = sites[idx]
    for pos = 0; pos < width; pos = pos + 1 {
      let char = site[pos:pos + 1].to_owned()
      for alpha_idx = 0; alpha_idx < alphabet_size; alpha_idx = alpha_idx + 1 {
        if char == dna_alphabet[alpha_idx] {
          counts[pos][alpha_idx] = counts[pos][alpha_idx] + 1
        }
      }
    }
  }

  counts
}

///|
pub fn pwm_from_counts(
  counts : Array[Array[Int]],
  pseudocount : Double,
) -> Array[Array[Double]] {
  if counts.length() == 0 {
    return []
  }
  let width = counts.length()
  let alphabet_size = counts[0].length()
  let pwm = Array::make(width, Array::make(alphabet_size, 0.0))

  for pos = 0; pos < width; pos = pos + 1 {
    let row = counts[pos]
    let mut total = 0.0
    for i = 0; i < row.length(); i = i + 1 {
      total = total + row[i].to_double()
    }
    let norm_total = total + pseudocount * alphabet_size.to_double()

    for alpha_idx = 0; alpha_idx < alphabet_size; alpha_idx = alpha_idx + 1 {
      pwm[pos][alpha_idx] = (row[alpha_idx].to_double() + pseudocount) /
        norm_total
    }
  }

  pwm
}

///|
pub fn pssm_from_pwm(
  pwm : Array[Array[Double]],
  background : Array[Double],
) -> Array[Array[Double]] {
  if pwm.length() == 0 {
    return []
  }
  let width = pwm.length()
  let alphabet_size = pwm[0].length()
  let pssm = Array::make(width, Array::make(alphabet_size, 0.0))

  for pos = 0; pos < width; pos = pos + 1 {
    let row = pwm[pos]
    for alpha_idx = 0; alpha_idx < alphabet_size; alpha_idx = alpha_idx + 1 {
      if row[alpha_idx] > 0.0 && background[alpha_idx] > 0.0 {
        pssm[pos][alpha_idx] = motif_log_nat(
          row[alpha_idx] / background[alpha_idx],
        )
      } else {
        pssm[pos][alpha_idx] = 0.0
      }
    }
  }

  pssm
}

///|
fn motif_log_nat(x : Double) -> Double {
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
pub fn create_motif(name : String, sites : Array[String]) -> Motif {
  let counts = count_matrix(sites)
  let pwm = pwm_from_counts(counts, 0.5)
  Motif::new(name, pwm, dna_alphabet, sites)
}

///|
pub fn score_sequence(
  pssm : Array[Array[Double]],
  sequence : String,
  alphabet : Array[String],
) -> Double {
  let mut score = 0.0
  for pos = 0; pos < pssm.length() && pos < sequence.length(); pos = pos + 1 {
    let char = sequence[pos:pos + 1].to_owned()
    let mut found = false
    for alpha_idx = 0
        alpha_idx < alphabet.length() && !found
        alpha_idx = alpha_idx + 1 {
      if char == alphabet[alpha_idx] {
        score = score + pssm[pos][alpha_idx]
        found = true
      }
    }
    if !found {
      score = score - 1.0
    }
  }
  score
}

///|
pub fn search_motif(
  motif : Motif,
  sequence : String,
  threshold : Double,
) -> Array[MotifHit] {
  let background = [0.25, 0.25, 0.25, 0.25]
  let pssm = pssm_from_pwm(motif.pwm, background)
  let width = pssm.length()
  let hits : Array[MotifHit] = Array::new()

  for pos = 0; pos <= sequence.length() - width; pos = pos + 1 {
    let subseq = sequence[pos:pos + width].to_owned()
    let score = score_sequence(pssm, subseq, motif.alphabet)
    if score >= threshold {
      hits.push(MotifHit::new(pos, score, "+", subseq))
    }
  }

  let rev_seq = motif_reverse_complement(sequence)

  for pos = 0; pos <= rev_seq.length() - width; pos = pos + 1 {
    let subseq = rev_seq[pos:pos + width].to_owned()
    let score = score_sequence(pssm, subseq, motif.alphabet)
    if score >= threshold {
      let orig_pos = sequence.length() - pos - width
      hits.push(MotifHit::new(orig_pos, score, "-", subseq))
    }
  }

  hits
}

///|
fn motif_reverse_complement(seq : String) -> String {
  let mut result = ""
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i:i + 1].to_owned()
    let comp = match c {
      "A" => "T"
      "T" => "A"
      "C" => "G"
      "G" => "C"
      _ => c
    }
    result = result + comp
  }
  result
}

///|
pub fn motif_consensus(motif : Motif) -> String {
  if motif.pwm.length() == 0 {
    return ""
  }

  let mut result = ""
  for pos = 0; pos < motif.pwm.length(); pos = pos + 1 {
    let row = motif.pwm[pos]
    let mut max_prob = 0.0
    let mut max_char = ""
    for alpha_idx = 0
        alpha_idx < motif.alphabet.length()
        alpha_idx = alpha_idx + 1 {
      if row[alpha_idx] > max_prob {
        max_prob = row[alpha_idx]
        max_char = motif.alphabet[alpha_idx]
      }
    }
    result = result + max_char
  }

  result
}

///|
pub fn motif_ic(motif : Motif) -> Double {
  if motif.pwm.length() == 0 {
    return 0.0
  }
  let background = 0.25
  let mut ic = 0.0

  for pos = 0; pos < motif.pwm.length(); pos = pos + 1 {
    let row = motif.pwm[pos]
    for alpha_idx = 0
        alpha_idx < motif.alphabet.length()
        alpha_idx = alpha_idx + 1 {
      let p = row[alpha_idx]
      if p > 0.0 {
        ic = ic + p * motif_log_nat(p / background)
      }
    }
  }

  ic
}

///|
pub fn meme_to_motif(meme_content : String) -> Motif {
  let lines = meme_content.split("\n").to_array()
  let mut name = ""
  let sites : Array[String] = Array::new()

  let mut in_sites = false
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.has_prefix("MOTIF") {
      let parts = line.split(" ").to_array()
      if parts.length() >= 2 {
        name = parts[1].to_owned()
      }
    } else if line.has_prefix("letter-probability") {
      in_sites = true
    } else if line == "" && in_sites {
      in_sites = false
    } else if in_sites && !line.has_prefix("#") {
      sites.push(line)
    }
  }

  let width = sites.length()
  let pwm = Array::make(width, Array::make(4, 0.0))

  for pos = 0; pos < width; pos = pos + 1 {
    let probs = sites[pos].split(" ").to_array()
    for alpha_idx = 0; alpha_idx < 4; alpha_idx = alpha_idx + 1 {
      if alpha_idx < probs.length() {
        pwm[pos][alpha_idx] = match parse_float_str(probs[alpha_idx].to_owned()) {
          Some(v) => v
          None => 0.0
        }
      }
    }
  }

  Motif::new(name, pwm, dna_alphabet, sites)
}

///|
pub fn motif_to_meme(motif : Motif) -> String {
  let header = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies:\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF " +
    motif.name +
    "\n"
  let matrix_header = "letter-probability matrix: alength= 4 w= " +
    motif.pwm.length().to_string() +
    "\n"

  let mut matrix_rows = ""
  for pos = 0; pos < motif.pwm.length(); pos = pos + 1 {
    let row = motif.pwm[pos]
    let mut row_str = ""
    for alpha_idx = 0
        alpha_idx < motif.alphabet.length()
        alpha_idx = alpha_idx + 1 {
      row_str = row_str + row[alpha_idx].to_string() + " "
    }
    matrix_rows = matrix_rows + row_str + "\n"
  }

  header + matrix_header + matrix_rows + "\n"
}

///|
pub fn align_motifs(motif1 : Motif, motif2 : Motif) -> Double {
  if motif1.pwm.length() == 0 || motif2.pwm.length() == 0 {
    return 0.0
  }
  let min_len = if motif1.pwm.length() < motif2.pwm.length() {
    motif1.pwm.length()
  } else {
    motif2.pwm.length()
  }
  let mut score = 0.0

  for pos = 0; pos < min_len; pos = pos + 1 {
    let row1 = motif1.pwm[pos]
    let row2 = motif2.pwm[pos]
    for alpha_idx = 0
        alpha_idx < motif1.alphabet.length()
        alpha_idx = alpha_idx + 1 {
      score = score + row1[alpha_idx] * row2[alpha_idx]
    }
  }

  score / min_len.to_double()
}

///|
/// Calculate per-position information content for a motif.
/// Returns an array of information content values for each position.
pub fn motif_per_position_ic(motif : Motif) -> Array[Double] {
  let result : Array[Double] = Array::new()
  if motif.pwm.length() == 0 { return result }
  
  let background = 0.25
  let mut pos = 0
  while pos < motif.pwm.length() {
    let row = motif.pwm[pos]
    let mut ic = 0.0
    let mut alpha_idx = 0
    while alpha_idx < motif.alphabet.length() {
      let p = row[alpha_idx]
      if p > 0.0 {
        ic = ic + p * log2_val(p / background)
      }
      alpha_idx = alpha_idx + 1
    }
    result.push(ic)
    pos = pos + 1
  }
  result
}

///|
/// Calculate the total information content of a motif in bits.
pub fn motif_total_ic(motif : Motif) -> Double {
  let per_pos = motif_per_position_ic(motif)
  let mut total = 0.0
  let mut i = 0
  while i < per_pos.length() {
    total = total + per_pos[i]
    i = i + 1
  }
  total
}

///|
/// Generate a sequence logo data from aligned sequences.
/// Returns an array of (position, letter, frequency, information_content).
pub fn generate_logo_data(sequences : Array[String]) -> Array[(Int, String, Double, Double)] {
  if sequences.length() == 0 { return Array::new() }
  
  let width = sequences[0].length()
  let logo_data : Array[(Int, String, Double, Double)] = Array::new()
  
  let alphabet = ["A", "C", "G", "T"]
  let n = sequences.length().to_double()
  
  let mut pos = 0
  while pos < width {
    let counts = Map([], capacity=4)
    let mut i = 0
    while i < sequences.length() {
      let base = sequences[i][pos:pos + 1].to_owned().to_upper()
      counts.set(base, counts.get(base).unwrap_or(0) + 1)
      i = i + 1
    }
    
    let background = 0.25
    let mut ic = 0.0
    let mut max_freq = 0.0
    let mut max_base = ""
    
    let mut ai = 0
    while ai < alphabet.length() {
      let base = alphabet[ai]
      let count = counts.get(base).unwrap_or(0).to_double()
      let freq = count / n
      
      if freq > max_freq {
        max_freq = freq
        max_base = base
      }
      
      if freq > 0.0 {
        ic = ic + freq * (log2_val(freq / background))
      }
      ai = ai + 1
    }
    
    logo_data.push((pos, max_base, max_freq, ic))
    pos = pos + 1
  }
  
  logo_data
}

///|
/// Calculate the similarity between two motifs using Pearson correlation.
pub fn motif_pearson_correlation(motif1 : Motif, motif2 : Motif) -> Double {
  if motif1.pwm.length() == 0 || motif2.pwm.length() == 0 {
    return 0.0
  }
  
  let min_len = if motif1.pwm.length() < motif2.pwm.length() {
    motif1.pwm.length()
  } else {
    motif2.pwm.length()
  }
  
  let n = (min_len * motif1.alphabet.length()).to_double()
  
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  
  let mut pos = 0
  while pos < min_len {
    let row1 = motif1.pwm[pos]
    let row2 = motif2.pwm[pos]
    let mut alpha_idx = 0
    while alpha_idx < motif1.alphabet.length() {
      let x = row1[alpha_idx]
      let y = row2[alpha_idx]
      sum_x = sum_x + x
      sum_y = sum_y + y
      sum_xy = sum_xy + x * y
      sum_x2 = sum_x2 + x * x
      sum_y2 = sum_y2 + y * y
      alpha_idx = alpha_idx + 1
    }
    pos = pos + 1
  }
  
  let numerator = n * sum_xy - sum_x * sum_x
  let denom1 = sqrt_val(n * sum_x2 - sum_x * sum_x)
  let denom2 = sqrt_val(n * sum_y2 - sum_y * sum_y)
  
  if denom1 * denom2 > 0.0 {
    numerator / (denom1 * denom2)
  } else {
    0.0
  }
}

///|
/// Perform motif enrichment analysis.
/// Returns a list of motifs found in a set of sequences with their enrichment scores.
pub fn motif_enrichment(
  motifs : Array[Motif],
  sequences : Array[String],
  background? : Array[String] = [],
) -> Array[(String, Double, Double)] {
  let results : Array[(String, Double, Double)] = Array::new()
  
  // Count motif occurrences in target sequences
  let motif_counts : Map[String, Int] = Map([], capacity=motifs.length())
  let seq_lengths : Map[String, Double] = Map([], capacity=motifs.length())
  
  let mut mi = 0
  while mi < motifs.length() {
    let motif = motifs[mi]
    let width = motif.pwm.length()
    let mut total_sites = 0
    
    let mut si = 0
    while si < sequences.length() {
      let seq = sequences[si]
      let mut pos = 0
      while pos <= seq.length() - width {
        let subseq = seq[pos:pos + width].to_owned()
        let score = score_sequence(pssm_from_pwm(motif.pwm, [0.25, 0.25, 0.25, 0.25]), subseq, motif.alphabet)
        if score > 0.0 {
          total_sites = total_sites + 1
        }
        pos = pos + 1
      }
      si = si + 1
    }
    
    motif_counts.set(motif.name, total_sites)
    seq_lengths.set(motif.name, sequences.length().to_double())
    mi = mi + 1
  }
  
  // Calculate enrichment score (simplified)
  let names = motif_counts.keys().collect()
  let mut ni = 0
  while ni < names.length() {
    let name = names[ni]
    let count = motif_counts[name]
    let total_seqs = sequences.length().to_double()
    
    let frequency = if total_seqs > 0.0 {
      count.to_double() / total_seqs
    } else {
      0.0
    }
    
    let enrichment = if frequency > 0.0 && total_seqs > 0.0 {
      log2_val(frequency * total_seqs)
    } else {
      0.0
    }
    
    results.push((name, count.to_double(), enrichment))
    ni = ni + 1
  }
  
  results
}

///|
/// Calculate the GC content of a motif (average GC frequency).
pub fn motif_gc_content(motif : Motif) -> Double {
  if motif.pwm.length() == 0 { return 0.0 }
  
  let mut gc_sum = 0.0
  let mut count = 0
  let mut pos = 0
  while pos < motif.pwm.length() {
    let row = motif.pwm[pos]
    // Assuming alphabet is [A, C, G, T]
    // C is index 1, G is index 2
    gc_sum = gc_sum + row[1] + row[2]
    count = count + 1
    pos = pos + 1
  }
  
  if count > 0 {
    gc_sum / count.to_double()
  } else {
    0.0
  }
}

///|
/// Generate a random motif with given properties.
pub fn generate_random_motif(
  name : String,
  length : Int,
  gc_content? : Double = 0.5,
) -> Motif {
  let pwm = Array::make(length, Array::make(4, 0.0))
  
  let at_content = 1.0 - gc_content
  
  let mut pos = 0
  while pos < length {
    let row = Array::make(4, 0.0)
    // A, T get equal share of AT content
    // C, G get equal share of GC content
    let at_share = at_content / 2.0
    let gc_share = gc_content / 2.0
    
    row[0] = at_share // A
    row[1] = gc_share // C
    row[2] = gc_share // G
    row[3] = at_share // T
    
    pwm[pos] = row
    pos = pos + 1
  }
  
  Motif::new(name, pwm, dna_alphabet, [])
}

///|
/// Parse a float string to Double.
fn parse_float_str(s : String) -> Double? {
  if s == "." || s.length() == 0 {
    return None
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  
  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  
  let mut has_dot = false
  let mut decimal_place = 1.0
  
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      has_dot = true
      i = i + 1
      continue
    }
    if c >= zero && c <= zero + 9 {
      let digit = (c - zero).to_double()
      result = result * 10.0 + digit
      if has_dot {
        decimal_place = decimal_place * 10.0
      }
    }
    i = i + 1
  }
  
  Some(result * sign / decimal_place)
}

///|
/// Calculate the square root of a number.
fn sqrt_val(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  let mut guess = x / 2.0
  let mut i = 0
  while i < 100 {
    guess = (guess + x / guess) / 2.0
    i = i + 1
  }
  guess
}

///|
/// Calculate log2 of a number.
fn log2_val(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  // Use natural log conversion: log2(x) = ln(x) / ln(2)
  let ln2 = 0.6931471805599453
  let mut result = 0.0
  let term = (x - 1.0) / (x + 1.0)
  let term2 = term * term
  let mut sum = term
  let mut cur_term = term
  let mut i = 1
  while i <= 50 {
    cur_term = cur_term * term2
    sum = sum + cur_term / (2 * i + 1).to_double()
    i = i + 1
  }
  let ln_x = 2.0 * sum
  ln_x / ln2
}
