/// biotite_infotheory3.mbt — Biotite round 3: information theory
/// 3 features: BLOSUM/PAM construction from alignments, mutual information,
/// Shannon entropy.

// #28 BLOSUM/PAM matrix construction from multiple sequence alignment
pub struct Biotite3SubstMatrix {
  alphabet : Array[String]
  scores : Array[Array[Double]]
  name : String
}

pub fn biotite3_build_blosum(aligned : Array[String], alphabet : Array[String], scale : Double) -> Biotite3SubstMatrix {
  let n = alphabet.length()
  let counts = Array(capacity=n)
  for _i in 0..<n {
    let row = Array(capacity=n)
    for _j in 0..<n {
      row.push(0.0)
    }
    counts.push(row)
  }
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  let mut total_pairs = 0.0
  for col in 0..<seq_len {
    for i in 0..<aligned.length() {
      let _start_j = i + 1
      for j in _start_j..<aligned.length() {
        if col < aligned[i].length() && col < aligned[j].length() {
          let c1 = aligned[i][col:col+1].to_owned()
          let c2 = aligned[j][col:col+1].to_owned()
          let mut idx1 = -1
          let mut idx2 = -1
          for a in 0..<n {
            if alphabet[a] == c1 { idx1 = a }
            if alphabet[a] == c2 { idx2 = a }
          }
          if idx1 >= 0 && idx2 >= 0 {
            counts[idx1][idx2] = counts[idx1][idx2] + 1.0
            counts[idx2][idx1] = counts[idx2][idx1] + 1.0
            total_pairs = total_pairs + 1.0
          }
        }
      }
    }
  }
  let freqs = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      if total_pairs > 0.0 {
        row.push(counts[i][j] / total_pairs)
      } else {
        row.push(0.0)
      }
    }
    freqs.push(row)
  }
  let marginals = Array(capacity=n)
  for i in 0..<n {
    let mut sum = 0.0
    for j in 0..<n {
      sum = sum + freqs[i][j]
    }
    marginals.push(sum)
  }
  let scores = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      let expected = marginals[i] * marginals[j]
      if expected > 0.0 && freqs[i][j] > 0.0 {
        let log_odds = @math.ln(freqs[i][j] / expected) / @math.ln(2.0) * scale
        row.push(log_odds.floor().to_int().to_double())
      } else {
        row.push(-99.0)
      }
    }
    scores.push(row)
  }
  { alphabet, scores, name: "BLOSUM" }
}

pub fn biotite3_build_pam(aligned : Array[String], alphabet : Array[String], scale : Double) -> Biotite3SubstMatrix {
  let blosum = biotite3_build_blosum(aligned, alphabet, scale)
  { ..blosum, name: "PAM" }
}

pub fn biotite3_subst_matrix_get_score(matrix : Biotite3SubstMatrix, sym1 : String, sym2 : String) -> Double {
  let mut idx1 = -1
  let mut idx2 = -1
  for i in 0..<matrix.alphabet.length() {
    if matrix.alphabet[i] == sym1 { idx1 = i }
    if matrix.alphabet[i] == sym2 { idx2 = i }
  }
  if idx1 >= 0 && idx2 >= 0 && idx1 < matrix.scores.length() && idx2 < matrix.scores[idx1].length() {
    matrix.scores[idx1][idx2]
  } else {
    -99.0
  }
}

pub fn biotite3_subst_matrix_to_string(matrix : Biotite3SubstMatrix) -> String {
  let sb = StringBuilder()
  let _ = sb.write("  ")
  for i in 0..<matrix.alphabet.length() {
    let _ = sb.write(matrix.alphabet[i] + " ")
  }
  let _ = sb.write("\n")
  for i in 0..<matrix.scores.length() {
    let _ = sb.write(matrix.alphabet[i] + " ")
    for j in 0..<matrix.scores[i].length() {
      let v = matrix.scores[i][j]
      let s = if v < -50.0 { "*" } else { v.to_string() }
      let _ = sb.write(s + " ")
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

// #29 mutual_information — column co-evolution
pub fn biotite3_mutual_information(aligned : Array[String], col1 : Int, col2 : Int, alphabet : Array[String]) -> Double {
  let n = alphabet.length()
  let joint : Array[Array[Double]] = Array(capacity=n)
  for _i in 0..<n {
    let row = Array(capacity=n)
    for _j in 0..<n {
      row.push(0.0)
    }
    joint.push(row)
  }
  let mut total = 0.0
  for seq in 0..<aligned.length() {
    if col1 < aligned[seq].length() && col2 < aligned[seq].length() {
      let c1 = aligned[seq][col1:col1+1].to_owned()
      let c2 = aligned[seq][col2:col2+1].to_owned()
      let mut idx1 = -1
      let mut idx2 = -1
      for a in 0..<n {
        if alphabet[a] == c1 { idx1 = a }
        if alphabet[a] == c2 { idx2 = a }
      }
      if idx1 >= 0 && idx2 >= 0 {
        joint[idx1][idx2] = joint[idx1][idx2] + 1.0
        total = total + 1.0
      }
    }
  }
  if total == 0.0 {
    return 0.0
  }
  for i in 0..<n {
    for j in 0..<n {
      joint[i][j] = joint[i][j] / total
    }
  }
  let p1 = Array(capacity=n)
  let p2 = Array(capacity=n)
  for i in 0..<n {
    let mut s1 = 0.0
    let mut s2 = 0.0
    for j in 0..<n {
      s1 = s1 + joint[i][j]
      s2 = s2 + joint[j][i]
    }
    p1.push(s1)
    p2.push(s2)
  }
  let mut mi = 0.0
  for i in 0..<n {
    for j in 0..<n {
      if joint[i][j] > 0.0 && p1[i] > 0.0 && p2[j] > 0.0 {
        mi = mi + joint[i][j] * @math.ln(joint[i][j] / (p1[i] * p2[j]))
      }
    }
  }
  mi
}

pub fn biotite3_mutual_information_matrix(aligned : Array[String], alphabet : Array[String]) -> Array[Array[Double]] {
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  let result = Array(capacity=seq_len)
  for i in 0..<seq_len {
    let row = Array(capacity=seq_len)
    for j in 0..<seq_len {
      if i == j {
        row.push(0.0)
      } else {
        row.push(biotite3_mutual_information(aligned, i, j, alphabet))
      }
    }
    result.push(row)
  }
  result
}

// #30 shannon_entropy — per-column conservation
pub fn biotite3_shannon_entropy(aligned : Array[String], col : Int, alphabet : Array[String]) -> Double {
  let n = alphabet.length()
  let counts = Array(capacity=n)
  for _i in 0..<n {
    counts.push(0.0)
  }
  let mut total = 0.0
  for seq in 0..<aligned.length() {
    if col < aligned[seq].length() {
      let ch = aligned[seq][col:col+1].to_owned()
      for a in 0..<n {
        if alphabet[a] == ch {
          counts[a] = counts[a] + 1.0
          total = total + 1.0
          break
        }
      }
    }
  }
  if total == 0.0 {
    return 0.0
  }
  let mut h = 0.0
  for i in 0..<n {
    let p = counts[i] / total
    if p > 0.0 {
      h = h - p * @math.ln(p)
    }
  }
  h
}

pub fn biotite3_shannon_entropy_profile(aligned : Array[String], alphabet : Array[String]) -> Array[Double] {
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  let result = Array(capacity=seq_len)
  for col in 0..<seq_len {
    result.push(biotite3_shannon_entropy(aligned, col, alphabet))
  }
  result
}

pub fn biotite3_relative_entropy(aligned : Array[String], col : Int, alphabet : Array[String], bg_freq : Array[Double]) -> Double {
  let n = alphabet.length()
  let counts = Array(capacity=n)
  for _i in 0..<n {
    counts.push(0.0)
  }
  let mut total = 0.0
  for seq in 0..<aligned.length() {
    if col < aligned[seq].length() {
      let ch = aligned[seq][col:col+1].to_owned()
      for a in 0..<n {
        if alphabet[a] == ch {
          counts[a] = counts[a] + 1.0
          total = total + 1.0
          break
        }
      }
    }
  }
  if total == 0.0 {
    return 0.0
  }
  let mut re = 0.0
  for i in 0..<n {
    let p = counts[i] / total
    let q = if i < bg_freq.length() { bg_freq[i] } else { 1.0 / n.to_double() }
    if p > 0.0 && q > 0.0 {
      re = re + p * @math.ln(p / q)
    }
  }
  re
}

pub fn _biotite_infotheory3_dummy() -> Unit {
  let _ = Biotite3SubstMatrix::{ alphabet: [], scores: [], name: "" }
}
