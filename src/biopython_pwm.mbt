/// biopython_pwm.mbt — Position Weight Matrix
/// Inspired by Biopython Bio.motifs: PWM, PSSM, log-odds scoring

///|
/// Position Weight Matrix (PWM) / Position Specific Scoring Matrix.
/// Built from a set of aligned sequences of equal length.
pub struct BiopythonPWM {
  alphabet : Array[String]
  matrix : Array[Array[Double]] // position -> alphabet_index -> log-odds score
  background : Array[Double] // background frequencies
  pseudocount : Double
}

///|
/// Build a PWM from aligned sequences.
pub fn BiopythonPWM::from_sequences(
  sequences : Array[String],
  alphabet~ : Array[String] = ["A", "C", "G", "T"],
  background~ : Array[Double] = [0.25, 0.25, 0.25, 0.25],
  pseudocount~ : Double = 0.1,
) -> BiopythonPWM {
  if sequences.length() == 0 {
    return { alphabet, matrix: [], background, pseudocount, }
  }
  let n_seqs = sequences.length()
  let aln_len = sequences[0].length()
  let n_sym = alphabet.length()
  let matrix : Array[Array[Double]] = Array::make(
    aln_len,
    Array::make(n_sym, 0.0),
  )
  for pos = 0; pos < aln_len; pos = pos + 1 {
    let counts : Array[Double] = Array::make(n_sym, pseudocount)
    for seq_idx = 0; seq_idx < n_seqs; seq_idx = seq_idx + 1 {
      let ch = sequences[seq_idx][pos:pos + 1].to_owned()
      for sym_idx = 0; sym_idx < n_sym; sym_idx = sym_idx + 1 {
        if ch == alphabet[sym_idx] {
          counts[sym_idx] = counts[sym_idx] + 1.0
          break
        }
      }
    }
    // Convert to log-odds
    let total = biopython_pwm_sum(counts)
    for sym_idx = 0; sym_idx < n_sym; sym_idx = sym_idx + 1 {
      let freq = counts[sym_idx] / total
      if freq > 0.0 && background[sym_idx] > 0.0 {
        matrix[pos][sym_idx] = @math.log2(freq / background[sym_idx])
      } else {
        matrix[pos][sym_idx] = -10.0
      }
    }
  }
  { alphabet, matrix, background, pseudocount, }
}

///|
/// Score a sequence against the PWM.
pub fn BiopythonPWM::score(self : BiopythonPWM, sequence : String) -> Double {
  let n = sequence.length()
  let pwm_len = self.matrix.length()
  if n < pwm_len {
    return -1000000000000000.0
  }
  let mut best_score = -1000000000000000.0
  // Slide PWM across sequence
  for start = 0; start <= n - pwm_len; start = start + 1 {
    let mut score = 0.0
    for pos = 0; pos < pwm_len; pos = pos + 1 {
      let ch = sequence[start + pos:start + pos + 1].to_owned()
      for sym_idx = 0; sym_idx < self.alphabet.length(); sym_idx = sym_idx + 1 {
        if ch == self.alphabet[sym_idx] {
          score = score + self.matrix[pos][sym_idx]
          break
        }
      }
    }
    if score > best_score {
      best_score = score
    }
  }
  best_score
}

///|
/// Scan a sequence and return all hit positions above threshold.
pub fn BiopythonPWM::scan(
  self : BiopythonPWM,
  sequence : String,
  threshold~ : Double = 0.0,
) -> Array[(Int32, Double)] {
  let n = sequence.length()
  let pwm_len = self.matrix.length()
  let hits : Array[(Int32, Double)] = Array(capacity=0)
  if n < pwm_len {
    return hits
  }
  for start = 0; start <= n - pwm_len; start = start + 1 {
    let mut score = 0.0
    for pos = 0; pos < pwm_len; pos = pos + 1 {
      let ch = sequence[start + pos:start + pos + 1].to_owned()
      for sym_idx = 0; sym_idx < self.alphabet.length(); sym_idx = sym_idx + 1 {
        if ch == self.alphabet[sym_idx] {
          score = score + self.matrix[pos][sym_idx]
          break
        }
      }
    }
    if score >= threshold {
      hits.push((start, score))
    }
  }
  hits
}

///|
/// Get the consensus sequence from the PWM.
pub fn BiopythonPWM::consensus(self : BiopythonPWM) -> String {
  let sb = StringBuilder()
  for pos = 0; pos < self.matrix.length(); pos = pos + 1 {
    let mut best_idx = 0
    let mut best_score = -1000000000000000.0
    for sym_idx = 0; sym_idx < self.alphabet.length(); sym_idx = sym_idx + 1 {
      if self.matrix[pos][sym_idx] > best_score {
        best_score = self.matrix[pos][sym_idx]
        best_idx = sym_idx
      }
    }
    let _ = sb.write(self.alphabet[best_idx])
  }
  sb.to_string()
}

///|
/// Get the length of the PWM.
pub fn BiopythonPWM::length(self : BiopythonPWM) -> Int32 {
  self.matrix.length()
}

///|
/// Helper: sum array of doubles
fn biopython_pwm_sum(arr : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < arr.length(); i = i + 1 {
    total = total + arr[i]
  }
  total
}
