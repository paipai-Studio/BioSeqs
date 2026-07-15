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
        pwm[pos][alpha_idx] = match parse_double(probs[alpha_idx].to_owned()) {
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
