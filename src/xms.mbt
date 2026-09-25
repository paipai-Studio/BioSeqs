///|
/// Bio.motifs.xms：XMS 加权位点模体格式（转录因子结合位点）。
/// Bio.motifs.microrna：microRNA 结合位点 motif（seed 分类与 PFM 构建）。
///
/// XMS 格式：可选 `>name` 头部，之后每行一条位点序列（可带权重）。
/// microRNA：fasta 风格 `>id comment` + 成熟序列，支持 6mer / 7mer-A1 /
/// 7mer-m8 / 8mer seed 分类与反向互补靶位点展示。

pub(all) struct XmsSite {
  sequence : String
  weight : Double
}

pub(all) struct XmsMotif {
  name : String
  sites : Array[XmsSite]
}

/// microRNA motif 记录。
pub(all) struct MirsRecord {
  id : String
  comment : String
  mature : String
}

/// 构造一条 XMS 位点。
pub fn XmsSite::new(sequence : String, weight? : Double = 1.0) -> XmsSite {
  XmsSite::{ sequence: sequence, weight: weight }
}

/// 去除首尾空白。
fn xms_trim(s : String) -> String {
  let mut i = 0
  let mut j = s.length()
  while i < j && (s[i] == ' ' || s[i] == '\t') {
    i = i + 1
  }
  while j > i && (s[j - 1] == ' ' || s[j - 1] == '\t' || s[j - 1] == '\r') {
    j = j - 1
  }
  if i >= j {
    ""
  } else {
    s[i:j].to_owned()
  }
}


/// 字符串全大写（逐字符）。
fn xms_upper(s : String) -> String {
  let buf = StringBuilder()
  for i in 0..<s.length() {
    buf.write_char(s.unsafe_get(i).unsafe_to_char().to_ascii_uppercase())
  }
  buf.to_string()
}
/// 解析 XMS 文本，返回 motif 列表（支持多个 `>name` 块）。
pub fn xms_parse(text : String) -> Array[XmsMotif] {
  let lines = text.split("\n")
  let motifs : Array[XmsMotif] = Array(capacity=0)
  let mut name = ""
  let mut sites : Array[XmsSite] = Array(capacity=0)
  for raw_view in lines {
    let line = xms_trim(raw_view.to_owned())
    if line.length() == 0 {
      continue
    }
    if line[0] == '>' {
      if sites.length() > 0 || name != "" {
        motifs.push(XmsMotif::{ name: name, sites: sites })
        sites = Array(capacity=0)
      }
      name = line[1:].to_owned()
      continue
    }
    let parts = line.split(" ").to_array()
    let seq_view = parts[0]
    let mut seq = seq_view.to_owned()
    let mut weight = 1.0
    if parts.length() > 1 {
      let w_str = xms_trim(parts[1].to_owned())
      if w_str.length() > 0 {
        weight = parse_double(w_str).unwrap_or(1.0)
      }
    }
    seq = xms_upper(seq)
    sites.push(XmsSite::{ sequence: seq, weight: weight })
  }
  if sites.length() > 0 || name != "" {
    motifs.push(XmsMotif::{ name: name, sites: sites })
  }
  motifs
}

/// microRNA 记录解析（fasta 风格）。
pub fn mirs_parse(text : String) -> Array[MirsRecord] {
  let lines = text.split("\n")
  let records : Array[MirsRecord] = Array(capacity=0)
  let mut id = ""
  let mut comment = ""
  for raw_view in lines {
    let line = xms_trim(raw_view.to_owned())
    if line.length() == 0 {
      continue
    }
    if line[0] == '>' {
      let sp = line[1:].split(" ").to_array()
      id = xms_trim(sp[0].to_owned())
      let mut comment_part = ""
      if sp.length() > 1 {
        let mut buf = ""
        for i in 1..<sp.length() {
          if i > 1 {
            buf = buf + " "
          }
          buf = buf + sp[i].to_owned()
        }
        comment_part = buf
      }
      comment = comment_part
    } else {
      records.push(MirsRecord::{ id: id, comment: comment, mature: xms_upper(line) })
    }
  }
  records
}

/// 加权计数矩阵（L x 4，列序 A/C/G/T）。
pub fn xms_counts(motif : XmsMotif) -> Array[Array[Double]] {
  let mut length = 0
  for s in motif.sites {
    if s.sequence.length() > length {
      length = s.sequence.length()
    }
  }
  let mat : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<length {
    mat.push([0.0, 0.0, 0.0, 0.0])
  }
  for s in motif.sites {
    for i in 0..<s.sequence.length() {
      let idx = xms_base_index(s.sequence.unsafe_get(i).unsafe_to_char())
      if idx >= 0 {
        let row = mat[i]
        row[idx] = row[idx] + s.weight
      }
    }
  }
  mat
}

/// 碱基 A/C/G/T 到列索引。
fn xms_base_index(b : Char) -> Int {
  if b == 'A' {
    return 0
  }
  if b == 'C' {
    return 1
  }
  if b == 'G' {
    return 2
  }
  if b == 'T' {
    return 3
  }
  -1
}

/// 计数矩阵行归一化 → 频率矩阵（PFM）。
pub fn xms_pfm(counts : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for row in counts {
    let mut total = 0.0
    for v in row {
      total = total + v
    }
    let new_row : Array[Double] = Array(capacity=0)
    if total > 0.0 {
      for v in row {
        new_row.push(v / total)
      }
    } else {
      for _ in row {
        new_row.push(0.25)
      }
    }
    out.push(new_row)
  }
  out
}

/// 频率矩阵 → log-odds PWM（背景频率数组，默认 0.25 均匀）。
pub fn xms_pwm(
  pfm : Array[Array[Double]],
  background? : Array[Double],
) -> Array[Array[Double]] {
  let bg = if background is Some(_) { background.unwrap() } else {
    [0.25, 0.25, 0.25, 0.25]
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  for row in pfm {
    let new_row : Array[Double] = Array(capacity=0)
    for i in 0..<4 {
      let ratio = if bg[i] > 0.0 {
        max_double(row[i], 1.0e-9) / bg[i]
      } else {
        1.0
      }
      new_row.push(xms_log2(ratio))
    }
    out.push(new_row)
  }
  out
}

fn max_double(a : Double, b : Double) -> Double {
  if a > b { a } else { b }
}

fn xms_log2(v : Double) -> Double {
  let ln2 = 0.6931471805599453
  if v <= 0.0 {
    0.0
  } else {
    @math.ln(v) / ln2
  }
}

/// 共识序列（每列 argmax）。
pub fn xms_consensus(pfm : Array[Array[Double]]) -> String {
  let sb = StringBuilder()
  for row in pfm {
    let mut best_idx = 0
    let mut best_val = row[0]
    for i in 1..<4 {
      if row[i] > best_val {
        best_val = row[i]
        best_idx = i
      }
    }
    sb.write_char(xms_base_char(best_idx))
  }
  sb.to_string()
}

fn xms_base_char(idx : Int) -> Char {
  if idx == 0 {
    return 'A'
  }
  if idx == 1 {
    return 'C'
  }
  if idx == 2 {
    return 'G'
  }
  'T'
}

/// 每位置信息含量（bits，最大 2.0）与总信息含量。
pub fn xms_info_content(pfm : Array[Array[Double]]) -> (Array[Double], Double) {
  let per_pos : Array[Double] = Array(capacity=0)
  let mut total = 0.0
  for row in pfm {
    let mut ic = 0.0
    for v in row {
      if v > 1.0e-9 {
        ic = ic + v * xms_log2(v)
      }
    }
    let bits = 2.0 + ic
    per_pos.push(bits)
    total = total + bits
  }
  (per_pos, total)
}

/// XMS 序列化为文本（含 `>name` 头，往返兼容）。
pub fn xms_write(motif : XmsMotif) -> String {
  let sb = StringBuilder()
  if motif.name.length() > 0 {
    sb.write_string(">" + motif.name + "\n")
  }
  for s in motif.sites {
    sb.write_string(s.sequence)
    if s.weight != 1.0 {
      sb.write_string(" " + s.weight.to_string())
    }
    sb.write_string("\n")
  }
  sb.to_string()
}

/// 全部位点序列。
pub fn xms_site_sequences(motif : XmsMotif) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for s in motif.sites {
    out.push(s.sequence)
  }
  out
}

/// microRNA：提取成熟序列的 seed（1-based 闭区间 [start, end]）。
pub fn mirs_seed(mature : String, start : Int, end : Int) -> String? {
  if start < 1 || end > mature.length() || start > end {
    return None
  }
  Some(mature[start - 1:end].to_owned())
}

/// 6mer seed：成熟序列第 2-7 位。
pub fn mirs_seed_6mer(mature : String) -> String? {
  mirs_seed(mature, 2, 7)
}

/// 7mer-A1 seed：第 1-7 位。
pub fn mirs_seed_7mer_a1(mature : String) -> String? {
  mirs_seed(mature, 1, 7)
}

/// 7mer-m8 seed：第 2-8 位。
pub fn mirs_seed_7mer_m8(mature : String) -> String? {
  mirs_seed(mature, 2, 8)
}

/// 8mer seed：第 1-8 位。
pub fn mirs_seed_8mer(mature : String) -> String? {
  mirs_seed(mature, 1, 8)
}

/// 按顺序列出该成熟序列的全部 seed 类型与序列。
pub fn mirs_classify_seeds(mature : String) -> Array[(String, String)] {
  let out : Array[(String, String)] = Array(capacity=0)
  let six = mirs_seed_6mer(mature)
  let a1 = mirs_seed_7mer_a1(mature)
  let m8 = mirs_seed_7mer_m8(mature)
  let eight = mirs_seed_8mer(mature)
  if eight is Some(_) {
    out.push(("8mer", eight.unwrap()))
  }
  if m8 is Some(_) {
    out.push(("7mer-m8", m8.unwrap()))
  }
  if a1 is Some(_) {
    out.push(("7mer-A1", a1.unwrap()))
  }
  if six is Some(_) {
    out.push(("6mer", six.unwrap()))
  }
  out
}

/// 反向互补（DNA）。
pub fn mirs_reverse_complement(seq : String) -> String {
  let sb = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    sb.write_char(mirs_complement(seq.unsafe_get(i).unsafe_to_char()))
    i = i - 1
  }
  sb.to_string()
}

fn mirs_complement(b : Char) -> Char {
  if b == 'A' {
    return 'T'
  }
  if b == 'C' {
    return 'G'
  }
  if b == 'G' {
    return 'C'
  }
  if b == 'T' {
    return 'A'
  }
  if b == 'U' {
    return 'A'
  }
  'N'
}

/// microRNA 位点计数矩阵（等价 xms_counts，接受序列数组）。
pub fn mirs_pfm(seqs : Array[String]) -> Array[Array[Double]] {
  let sites : Array[XmsSite] = Array(capacity=0)
  for s in seqs {
    sites.push(XmsSite::{ sequence: s, weight: 1.0 })
  }
  let motif = XmsMotif::{ name: "mirna", sites: sites }
  xms_pfm(xms_counts(motif))
}