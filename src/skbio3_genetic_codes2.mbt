// Scikit-bio third-round: NCBI genetic code tables (items 30-32) + alignment/composition (items 33-34)
// genetic_code_bacterial, genetic_code_alt_yeast, genetic_code_ciliate
// alignment_identity_matrix, multiplicative_replacement

///|
/// NCBI genetic code table 11: Bacterial/Archaeal/Plant Plastid。
///
/// 数据来源: 全库权威密码子表 codon_table(11) (见 src/codon_utils.mbt)。
/// 相对标准表的差异: TGA = * (终止), ATA = I, AGA/AGG = R。
pub fn skbio3_genetic_code_bacterial() -> Map[String, String] {
  codon_table_bacterial().forward_table()
}

///|
/// NCBI genetic code table 12: Alternative Yeast Nuclear。
///
/// 数据来源: 全库权威密码子表 codon_table(12)。
/// 相对标准表的差异: CTG = S (Ser) 而非 L (Leu)。
pub fn skbio3_genetic_code_alt_yeast() -> Map[String, String] {
  codon_table_alt_yeast().forward_table()
}

///|
/// NCBI genetic code table 6: Ciliate Nuclear。
///
/// 数据来源: 全库权威密码子表 codon_table(6)。
/// 相对标准表的差异: TAA = Q, TAG = Q (不再是终止密码子)。
pub fn skbio3_genetic_code_ciliate() -> Map[String, String] {
  codon_table_ciliate().forward_table()
}

///|
/// 使用自定义遗传密码表翻译 DNA 序列 (兼容入口)。
///
/// 实现已收敛到统一翻译入口 translate_with_codon_map,
/// 因此自动支持 IUPAC 简并密码子与大小写归一。
///
/// 参数:
///   - seq: DNA 序列
///   - table: 密码子 → 氨基酸 映射
pub fn skbio3_translate_with_table(
  seq : String,
  table : Map[String, String],
) -> String {
  translate_with_codon_map(seq, table)
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
