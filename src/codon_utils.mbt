///|
/// Bio.SeqUtils Codon Analysis — MoonBit 实现
/// 对应 Biopython Bio.SeqUtils 核心算法（缺失功能补全）
///   - GC123: 密码子三个位置的 GC 含量 (返回 Array[Double])
///   - six_frame_translations: 六框翻译 (返回结构化结果，支持多密码子表)
///   - CodonAdaptationIndex (CAI): 密码子适应指数
///   - Effective Number of Codons (Nc): Wright 1990 精确版
///
/// 注意: 与 sequtils_advanced.mbt / codon_align.mbt 重名的函数已加 seq_ 前缀区分

// ============ NCBI 密码子表 ============

/// 标准密码子表 (NCBI Table 1)
let _codon_table_standard : Map[String, String] = Map(
  [
    ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
    ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
    ("ATT", "I"), ("ATC", "I"), ("ATA", "I"), ("ATG", "M"),
    ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
    ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
    ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
    ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
    ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
    ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
    ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
    ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
    ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
    ("TGT", "C"), ("TGC", "C"), ("TGA", "*"), ("TGG", "W"),
    ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
    ("AGT", "S"), ("AGC", "S"), ("AGA", "R"), ("AGG", "R"),
    ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
  ],
  capacity=64,
)

/// Vertebrate Mitochondrial (NCBI Table 2)
let _codon_table_vertebrate_mito : Map[String, String] = Map(
  [
    ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
    ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
    ("ATT", "I"), ("ATC", "I"), ("ATA", "M"), ("ATG", "M"),
    ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
    ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
    ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
    ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
    ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
    ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
    ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
    ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
    ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
    ("TGT", "C"), ("TGC", "C"), ("TGA", "W"), ("TGG", "W"),
    ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
    ("AGT", "S"), ("AGC", "S"), ("AGA", "*"), ("AGG", "*"),
    ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
  ],
  capacity=64,
)

/// Yeast Mitochondrial (NCBI Table 3)
let _codon_table_yeast_mito : Map[String, String] = Map(
  [
    ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
    ("CTT", "T"), ("CTC", "T"), ("CTA", "T"), ("CTG", "T"),
    ("ATT", "I"), ("ATC", "I"), ("ATA", "M"), ("ATG", "M"),
    ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
    ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
    ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
    ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
    ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
    ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
    ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
    ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
    ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
    ("TGT", "C"), ("TGC", "C"), ("TGA", "W"), ("TGG", "W"),
    ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
    ("AGT", "S"), ("AGC", "S"), ("AGA", "R"), ("AGG", "R"),
    ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
  ],
  capacity=64,
)

///|
/// 获取 NCBI 密码子表
pub fn codon_table(table_id : Int) -> Map[String, String] {
  match table_id {
    2 => _codon_table_vertebrate_mito
    3 => _codon_table_yeast_mito
    _ => _codon_table_standard
  }
}

///|
/// 六框翻译结果
pub struct SixFrameTranslation {
  frame_id : String      // "Fwd-1", "Fwd-2", "Fwd-3", "Rev-1", ...
  strand : String        // "+" 或 "-"
  frame_offset : Int      // 0, 1, 2
  protein : String       // 翻译后的氨基酸序列
}

// ============ 内部辅助函数 ============

///|
/// 单框翻译（内部）
fn _cu_translate_frame(seq : String, frame : Int, codon_map : Map[String, String]) -> String {
  let n = seq.length()
  let buf = StringBuilder()
  let mut i = frame
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    let aa = match codon_map.get(codon) {
      Some(a) => a
      None => "X"
    }
    buf.write_string(aa)
    i = i + 3
  }
  buf.to_string()
}

///|
/// 密码子频率统计（内部）
fn _cu_codon_counts(sequences : Array[String], codon_map : Map[String, String]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=64)
  let codon_keys = codon_map.keys().to_array()
  let n = codon_keys.length()
  let mut c = 0
  while c < n {
    counts[codon_keys[c]] = 0
    c = c + 1
  }

  let n_seqs = sequences.length()
  let mut s = 0
  while s < n_seqs {
    let seq = sequences[s].to_upper()
    let mut i = 0
    while i + 3 <= seq.length() {
      let codon = seq[i:i + 3].to_owned().to_upper()
      if counts.contains(codon) {
        counts[codon] = counts[codon] + 1
      }
      i = i + 3
    }
    s = s + 1
  }
  counts
}

///|
/// 按氨基酸分组密码子（内部）
fn _cu_aa_groups(codon_map : Map[String, String]) -> Map[String, Array[String]] {
  let groups : Map[String, Array[String]] = Map([], capacity=20)
  let codon_keys = codon_map.keys().to_array()
  let n = codon_keys.length()
  let mut k = 0
  while k < n {
    let codon = codon_keys[k]
    match codon_map.get(codon) {
      Some(aa) => {
        if aa != "*" {
          if !groups.contains(aa) {
            groups[aa] = Array::new()
          }
          groups[aa].push(codon)
        }
      }
      None => ()
    }
    k = k + 1
  }
  groups
}

// ============ 公开 API ============

///|
/// **GC123** — 密码子三个位置的 GC 含量
/// 对应 Biopython Bio.SeqUtils.GC123
/// 返回 [GC_pos1, GC_pos2, GC_pos3]（比例值 0.0-1.0）
///
/// 注意: 本函数返回 Array[Double]，与 sequtils_advanced.mbt 的 gc123 返回 GC123Result 不同
///
/// 例: seq_gc123("ATGGCCATTGTA") => [0.5, 0.5, 0.333333]
pub fn seq_gc123(seq : String) -> Array[Double] {
  let seq_upper = seq.to_upper()
  let gc_counts = [0, 0, 0]
  let total_counts = [0, 0, 0]

  let mut i = 0
  while i + 3 <= seq_upper.length() {
    let pos0 = seq_upper[i:i + 1].to_owned()
    let pos1 = seq_upper[i + 1:i + 2].to_owned()
    let pos2 = seq_upper[i + 2:i + 3].to_owned()
    let bases = [pos0, pos1, pos2]

    let mut j = 0
    while j < 3 {
      let b = bases[j]
      if b == "G" || b == "C" {
        gc_counts[j] = gc_counts[j] + 1
      }
      if b == "A" || b == "T" || b == "G" || b == "C" {
        total_counts[j] = total_counts[j] + 1
      }
      j = j + 1
    }
    i = i + 3
  }

  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < 3 {
    if total_counts[j] == 0 {
      result.push(0.0)
    } else {
      result.push(gc_counts[j].to_double() / total_counts[j].to_double())
    }
    j = j + 1
  }
  result
}

///|
/// **六框翻译** — 对 DNA 序列进行所有 6 种阅读框的翻译
/// 对应 Biopython Bio.SeqUtils.six_frame_translations
///
/// 支持多密码子表 (table 1/2/3)，返回结构化结果
///
/// 注意: 本函数返回 Array[SixFrameTranslation]，与 sequtils_advanced.mbt 的
/// six_frame_translations 返回 String raise SeqError 不同
///
/// 参数:
///   - seq: DNA 序列
///   - table: NCBI 密码子表编号 (1=标准, 2=Vertebrate Mito, 3=Yeast Mito)
///
/// 例:
///   seq_six_frame_translations("ATGGCCATTGTAATGGGCCGCTGAA", table=1)
///   返回 Fwd-1/2/3 + Rev-1/2/3 共 6 个框架
pub fn seq_six_frame_translations(seq : String, table? : Int = 1) -> Array[SixFrameTranslation] {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)
  let result : Array[SixFrameTranslation] = Array::new()

  let mut frame = 0
  while frame < 3 {
    let protein = _cu_translate_frame(seq_upper, frame, codon_map)
    result.push(SixFrameTranslation::{
      frame_id: "Fwd-" + (frame + 1).to_string(),
      strand: "+",
      frame_offset: frame,
      protein: protein,
    })
    frame = frame + 1
  }

  let rc_seq = reverse_complement(seq_upper)
  let mut frame2 = 0
  while frame2 < 3 {
    let protein = _cu_translate_frame(rc_seq, frame2, codon_map)
    result.push(SixFrameTranslation::{
      frame_id: "Rev-" + (frame2 + 1).to_string(),
      strand: "-",
      frame_offset: frame2,
      protein: protein,
    })
    frame2 = frame2 + 1
  }

  result
}

///|
/// **Codon Adaptation Index (CAI) 权重** — 计算密码子使用偏倚度权重
/// 对应 Biopython Bio.SeqUtils.CodonAdaptationIndex
///
/// 算法 (Sharp & Li 1987):
///   w(codon) = freq(codon) / max_freq(同氨基酸密码子组)
///
/// 参数:
///   - sequences: 参考序列集合（通常是高表达基因）
///   - table: NCBI 密码子表编号
///
/// 返回: Map[String, Double] — 每个密码子的 CAI 权重 (0.0-1.0)
///   高 w 值表示该密码子在参考集中是"偏好密码子"
///
/// 例:
///   let w = cai_weights(["ATGGCCAAGG...", "ATGGCCTACAT..."])
///   // w["ATG"] = 1.0 (Met 只有一个密码子)
///   // w["CTG"] = 1.0, w["CTT"] = 0.45 (若 CTG 是最常用 Leu 密码子)
pub fn cai_weights(sequences : Array[String], table? : Int = 1) -> Map[String, Double] {
  let codon_map = codon_table(table)
  let counts = _cu_codon_counts(sequences, codon_map)
  let codon_keys = codon_map.keys().to_array()

  // 计算总密码子数
  let mut total = 0
  let n_codons = codon_keys.length()
  let mut k = 0
  while k < n_codons {
    match counts.get(codon_keys[k]) {
      Some(c) => total = total + c
      None => ()
    }
    k = k + 1
  }

  // 计算每个密码子的频率
  let freqs : Map[String, Double] = Map([], capacity=64)
  let mut k2 = 0
  while k2 < n_codons {
    let codon = codon_keys[k2]
    if total > 0 {
      match counts.get(codon) {
        Some(c) => freqs[codon] = c.to_double() / total.to_double()
        None => freqs[codon] = 0.0
      }
    } else {
      freqs[codon] = 0.0
    }
    k2 = k2 + 1
  }

  // 按氨基酸分组
  let aa_groups = _cu_aa_groups(codon_map)
  let aa_keys = aa_groups.keys().to_array()
  let n_aa = aa_keys.length()

  // 计算权重
  let weights : Map[String, Double] = Map([], capacity=64)
  let mut a = 0
  while a < n_aa {
    let aa = aa_keys[a]
    let aa_codons = match aa_groups.get(aa) {
      Some(c) => c
      None => Array::new()
    }

    // 找该氨基酸密码子组的最大频率
    let mut max_freq = 0.0
    let n_aa_codons = aa_codons.length()
    let mut ac = 0
    while ac < n_aa_codons {
      let codon = aa_codons[ac]
      match freqs.get(codon) {
        Some(f) => if f > max_freq { max_freq = f }
        None => ()
      }
      ac = ac + 1
    }

    // w = freq / max_freq，若 max_freq=0 则全部设为 1.0 (无偏倚)
    let mut ac2 = 0
    while ac2 < n_aa_codons {
      let codon = aa_codons[ac2]
      if max_freq > 0.0 {
        match freqs.get(codon) {
          Some(f) => weights[codon] = f / max_freq
          None => weights[codon] = 1.0
        }
      } else {
        weights[codon] = 1.0
      }
      ac2 = ac2 + 1
    }
    a = a + 1
  }

  weights
}

///|
/// **Codon Adaptation Index (CAI)** — 计算单条序列的密码子适应指数
/// 对应 Biopython Bio.SeqUtils.CodonAdaptationIndex
///
/// 算法 (Sharp & Li 1987):
///   CAI = exp(mean(ln(w_i)))
///   其中 w_i 是每个密码子的 CAI 权重
///
/// 参数:
///   - seq: 编码 DNA 序列 (必须是 3 的倍数)
///   - weights: 由 cai_weights() 计算的密码子权重
///   - table: NCBI 密码子表编号
///
/// 返回: CAI 值 (0.0 - 1.0)
///   - 接近 1.0: 密码子使用高度优化 (偏向参考集的偏好密码子)
///   - 低值: 密码子使用与参考集差异大
///
/// 例:
///   let ref_set = ["ATGGCCAAG...", "ATGGCCTAC...", ...]
///   let w = cai_weights(ref_set)
///   let cai = codon_adaptation_index("ATGGCCAAGATG...", w)
pub fn codon_adaptation_index(
  seq : String,
  weights : Map[String, Double],
  table? : Int = 1,
) -> Double {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)

  let ln_w_values : Array[Double] = Array::new()
  let mut i = 0
  while i + 3 <= seq_upper.length() {
    let codon = seq_upper[i:i + 3].to_owned().to_upper()
    // 必须是有效密码子且不是终止密码
    let is_valid = match codon_map.get(codon) {
      Some(aa) => aa != "*"
      None => false
    }
    if is_valid {
      let w_opt = match weights.get(codon) {
        Some(w) => w
        None => 0.0
      }
      if w_opt > 0.0 {
        ln_w_values.push(@math.ln(w_opt))
      }
    }
    i = i + 3
  }

  if ln_w_values.length() == 0 {
    return 0.0
  }

  let mut sum = 0.0
  let n = ln_w_values.length()
  let mut j = 0
  while j < n {
    sum = sum + ln_w_values[j]
    j = j + 1
  }

  @math.exp(sum / n.to_double())
}

///|
/// **Effective Number of Codons (Nc/ENC)** — Wright 1990 精确算法
/// 对应 Biopython Bio.SeqUtils (codon_align.mbt 有简化版)
///
/// 算法 (Wright 1990):
///   F-hat = (n * sum(p_i^2) - 1) / (n - 1)
///   Nc_aa = 1 / F-hat
///   最终 Nc = sum(Nc_aa for 每个氨基酸)
///   Met 和 Trp 各只有 1 密码子，跳过不计入
///
/// 参数:
///   - seq: 编码 DNA 序列 (3 的倍数)
///   - table: NCBI 密码子表编号
///
/// 返回: Nc 值 (约 20.0 - 61.0)
///   - ~20: 密码子使用极度偏倚 (每个氨基酸只用 1-2 个密码子)
///   - ~61: 密码子使用完全均匀 (61 个有义密码子均匀使用)
///
/// 例:
///   let nc = seq_enc_wright("ATGGCCAAGATGCT...")
pub fn seq_enc_wright(seq : String, table? : Int = 1) -> Double {
  let seq_upper = seq.to_upper()
  let codon_map = codon_table(table)

  // 统计密码子频率
  let codon_counts = _cu_codon_counts([seq_upper], codon_map)

  // 按氨基酸分组
  let aa_groups = _cu_aa_groups(codon_map)
  let aa_keys = aa_groups.keys().to_array()
  let n_aa = aa_keys.length()

  let mut total_nc = 0.0
  let mut skipped_met_trp = 0

  let mut a = 0
  while a < n_aa {
    let aa = aa_keys[a]
    let aa_codons = match aa_groups.get(aa) {
      Some(c) => c
      None => Array::new()
    }

    // Met (M) 和 Trp (W) 只有 1 个密码子，跳过
    if aa_codons.length() <= 1 {
      skipped_met_trp = skipped_met_trp + 1
      a = a + 1
      continue
    }

    // 该氨基酸的密码子总数
    let mut total_c = 0
    let n_aa_c = aa_codons.length()
    let mut ac = 0
    while ac < n_aa_c {
      match codon_counts.get(aa_codons[ac]) {
        Some(c) => total_c = total_c + c
        None => ()
      }
      ac = ac + 1
    }

    if total_c <= 1 {
      // 计数太少，Nc_aa 取密码子种类数
      total_nc = total_nc + aa_codons.length().to_double()
    } else {
      // 计算 sum(p_i^2)
      let mut sum_p2 = 0.0
      let mut ac2 = 0
      while ac2 < n_aa_c {
        let codon = aa_codons[ac2]
        let count = match codon_counts.get(codon) {
          Some(c) => c
          None => 0
        }
        let p = count.to_double() / total_c.to_double()
        sum_p2 = sum_p2 + p * p
        ac2 = ac2 + 1
      }

      // F-hat
      let f_hat = (total_c.to_double() * sum_p2 - 1.0) / (total_c.to_double() - 1.0)
      if f_hat > 0.0 {
        total_nc = total_nc + 1.0 / f_hat
      } else {
        total_nc = total_nc + aa_codons.length().to_double()
      }
    }

    a = a + 1
  }

  total_nc
}
