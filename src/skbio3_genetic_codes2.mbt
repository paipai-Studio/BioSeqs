// Scikit-bio third-round: NCBI genetic code tables (items 30-32) + alignment/composition (items 33-34)
// genetic_code_bacterial, genetic_code_alt_yeast, genetic_code_ciliate
// alignment_identity_matrix, multiplicative_replacement

///|
pub fn skbio3_genetic_code_bacterial() -> Map[String, String] {
  // NCBI genetic code table 11: Bacterial/Archaeal/Plant Plastid
  // Differences from standard: TGA = W (not stop), TAG = stop, ATA = I, AGA/AGG = R
  let table : Map[String, String] = Map([])
  let standard_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIITTTTNNKKSSRRVVVVAAAADDEEGGGG"
  // Build codon table: TGA -> W instead of *
  let bases = "TCAG"
  for i in 0..<4 {
    for j in 0..<4 {
      for k in 0..<4 {
        let codon = {
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(bases[i].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[j].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[k].to_int()))
          sb.to_string()
        }
        let idx = i * 16 + j * 4 + k
        let mut aa = {
          let c = standard_aas[idx].to_int()
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(c))
          sb.to_string()
        }
        // Bacterial code differences:
        // TGA -> W (opal)
        if codon == "TGA" {
          aa = "W"
        }
        // AGA -> R, AGG -> R (same as standard, no change)
        // ATA -> I (same as standard, no change)
        table[codon] = aa
      }
    }
  }
  table
}

///|
pub fn skbio3_genetic_code_alt_yeast() -> Map[String, String] {
  // NCBI genetic code table 12: Alternative Yeast Nuclear
  // Differences from standard: CTG -> S (Ser) instead of L (Leu)
  let table : Map[String, String] = Map([])
  let standard_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIITTTTNNKKSSRRVVVVAAAADDEEGGGG"
  let bases = "TCAG"
  for i in 0..<4 {
    for j in 0..<4 {
      for k in 0..<4 {
        let codon = {
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(bases[i].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[j].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[k].to_int()))
          sb.to_string()
        }
        let idx = i * 16 + j * 4 + k
        let mut aa = {
          let c = standard_aas[idx].to_int()
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(c))
          sb.to_string()
        }
        // Alternative yeast: CTG -> S
        if codon == "CTG" {
          aa = "S"
        }
        table[codon] = aa
      }
    }
  }
  table
}

///|
pub fn skbio3_genetic_code_ciliate() -> Map[String, String] {
  // NCBI genetic code table 6: Ciliate Nuclear
  // Differences: TAA -> Q, TAG -> Q (not stop)
  let table : Map[String, String] = Map([])
  let standard_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIITTTTNNKKSSRRVVVVAAAADDEEGGGG"
  let bases = "TCAG"
  for i in 0..<4 {
    for j in 0..<4 {
      for k in 0..<4 {
        let codon = {
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(bases[i].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[j].to_int()))
          let _ = sb.write_char(Int32::unsafe_to_char(bases[k].to_int()))
          sb.to_string()
        }
        let idx = i * 16 + j * 4 + k
        let mut aa = {
          let c = standard_aas[idx].to_int()
          let sb = StringBuilder()
          let _ = sb.write_char(Int32::unsafe_to_char(c))
          sb.to_string()
        }
        // Ciliate: TAA -> Q, TAG -> Q
        if codon == "TAA" || codon == "TAG" {
          aa = "Q"
        }
        table[codon] = aa
      }
    }
  }
  table
}

///|
pub fn skbio3_translate_with_table(
  seq : String,
  table : Map[String, String],
) -> String {
  // Translate a DNA sequence using a custom genetic code table
  let n = seq.length()
  let sb = StringBuilder()
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned()
    let aa = table.get(codon).unwrap_or("X")
    let _ = sb.write(aa)
    i = i + 3
  }
  sb.to_string()
}

///|
pub fn skbio3_alignment_identity_matrix(
  sequences : Array[String],
) -> Array[Array[Double]] {
  // Pairwise identity matrix for sequences in an alignment
  let n = sequences.length()
  let matrix : Array[Array[Double]] = []
  for _ in 0..<n {
    let row : Array[Double] = []
    for _ in 0..<n {
      row.push(0.0)
    }
    matrix.push(row)
  }
  for i in 0..<n {
    matrix[i][i] = 1.0
    let _next = i + 1
    for j in _next..<n {
      let seq_i = sequences[i]
      let seq_j = sequences[j]
      let len = if seq_i.length() < seq_j.length() {
        seq_i.length()
      } else {
        seq_j.length()
      }
      if len == 0 {
        matrix[i][j] = 0.0
        matrix[j][i] = 0.0
      } else {
        let mut matches = 0
        for k in 0..<len {
          if seq_i[k] == seq_j[k] {
            matches = matches + 1
          }
        }
        let identity = matches.to_double() / len.to_double()
        matrix[i][j] = identity
        matrix[j][i] = identity
      }
    }
  }
  matrix
}

///|
pub fn skbio3_multiplicative_replacement(
  composition : Array[Double],
  count : Double,
) -> Array[Double] {
  // Multiplicative replacement for compositional data with zeros
  // Replaces zeros with a small value delta, renormalizes
  // delta = (count of zeros) / (total count * D) where D is dimension
  let d = composition.length()
  if d == 0 {
    return []
  }
  // Count zeros
  let mut n_zeros = 0
  let mut sum_nonzero = 0.0
  for i in 0..<d {
    if composition[i] <= 0.0 {
      n_zeros = n_zeros + 1
    } else {
      sum_nonzero = sum_nonzero + composition[i]
    }
  }
  if n_zeros == 0 {
    return composition.copy()
  }
  // delta = (n_zeros) / (count * d)
  let delta = if count > 0.0 {
    n_zeros.to_double() / (count * d.to_double())
  } else {
    0.01 / d.to_double()
  }
  // Ensure delta is small enough
  let max_delta = if sum_nonzero > 0.0 {
    sum_nonzero / d.to_double() * 0.5
  } else {
    0.001
  }
  let actual_delta = if delta > max_delta { max_delta } else { delta }
  // Replace zeros and rescale non-zeros
  let result : Array[Double] = []
  let scale = 1.0 - n_zeros.to_double() * actual_delta
  for i in 0..<d {
    if composition[i] <= 0.0 {
      result.push(actual_delta)
    } else {
      result.push(composition[i] * scale / sum_nonzero)
    }
  }
  result
}
