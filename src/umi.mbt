///|
/// umi.mbt — UMI(分子标签)提取、聚类、去重与共识.
///
/// 对齐 umi_tools / fgbio 的核心语义(本地确定性实现, 数据不出本机):
/// - read header UMI 提取: 下划线后缀(`read1_ACGT`)、`#`、Illumina 冒号、`|`
/// - 从 read 5' 端剪除 UMI 碱基(序列与质量同步裁剪)
/// - 按 `(UMI, chrom, pos, strand)` 分组为 read family
/// - directional UMI 纠错: 邻近一个编辑距离的低频 UMI 并入同 locus 的
///   高频代表 UMI
/// - 质量加权逐位多数投票生成 family consensus
/// - 汇总统计与 JSON 输出(衔接 Web/报告架构)

///|
/// UMI 处理错误.
pub suberror UmiError {
  UmiError(String)
}

///|
/// read header 中 UMI 的常见分隔约定.
pub(all) enum UmiHeaderFormat {
  /// `@read1_ACGTACGT`
  Underscore
  /// `@read1#ACGTACGT`
  Hash
  /// `@MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y:ACGTACGT`
  Illumina
  /// `@read1|ACGTACGT`
  Pipe
}

///|
/// 解析后的 header.
pub(all) struct UmiHeader {
  /// 去掉 UMI 后的 read 名称
  name : String
  /// 提取到的 UMI(大写)
  umi : String
}

///|
/// 分隔符.
fn umi_separator(format : UmiHeaderFormat) -> String {
  match format {
    Underscore => "_"
    Hash => "#"
    Illumina => ":"
    Pipe => "|"
  }
}

///|
/// 是否为合法的 UMI 碱基(A/C/G/T/N, 大小写均可).
fn umi_is_valid_base(code : Int32) -> Bool {
  match code {
    65 | 67 | 71 | 84 | 78 => true
    97 | 99 | 103 | 116 | 110 => true
    _ => false
  }
}

///|
/// 校验 UMI 字符串.
pub fn umi_is_valid(s : String) -> Bool {
  if s.length() == 0 {
    return false
  }
  for i = 0; i < s.length(); i = i + 1 {
    if !umi_is_valid_base(s.unsafe_get(i).to_int()) {
      return false
    }
  }
  true
}

///|
/// 从 header 解析 UMI 与 read 名称.
pub fn umi_parse_header(
  header : String,
  format : UmiHeaderFormat,
) -> UmiHeader raise UmiError {
  // 只取第一段(忽略 description), 避免把描述中的分隔符当成 UMI.
  let trimmed = header.trim()
  let mut first = trimmed.split(" ").to_array()[0].to_owned()
  // FASTQ/FASTA 记录前缀不属于 ID.
  if first.has_prefix("@") || first.has_prefix(">") {
    first = first[1:].to_owned()
  }
  let sep = umi_separator(format)
  let parts = first.split(sep).to_array()
  if parts.length() < 2 {
    raise UmiError("header `\{first}` has no `\{sep}` separated UMI")
  }
  let umi = parts[parts.length() - 1].to_owned().to_upper()
  if !umi_is_valid(umi) {
    raise UmiError("header `\{first}` has invalid UMI `\{umi}`")
  }
  let sb = StringBuilder()
  for i = 0; i < parts.length() - 1; i = i + 1 {
    if i > 0 {
      sb.write_string(sep)
    }
    sb.write_string(parts[i].to_owned())
  }
  { name: sb.to_string(), umi, }
}

///|
/// 只提取 UMI.
pub fn umi_extract_umi(
  header : String,
  format : UmiHeaderFormat,
) -> String raise UmiError {
  umi_parse_header(header, format).umi
}

///|
/// 只取出去掉 UMI 后的 read 名称.
pub fn umi_strip_header(
  header : String,
  format : UmiHeaderFormat,
) -> String raise UmiError {
  umi_parse_header(header, format).name
}

///|
/// FASTQ 质量串(Phred+33)转质量数组.
pub fn umi_qualities_from_ascii(quality : String) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=quality.length())
  for i = 0; i < quality.length(); i = i + 1 {
    out.push(quality.unsafe_get(i).to_int() - 33)
  }
  out
}

///|
/// 质量数组转 FASTQ 质量串(Phred+33).
pub fn umi_qualities_to_ascii(qualities : Array[Int32]) -> String {
  let sb = StringBuilder()
  for q in qualities {
    sb.write_char((q + 33).to_uint16().unsafe_to_char())
  }
  sb.to_string()
}

///|
/// 单条带 UMI 的 read.
pub(all) struct UmiRead {
  id : String
  sequence : String
  qualities : Array[Int32]
  umi : String
  chrom : String
  pos : Int32
  strand : Int32
}

///|
/// 构造并校验一条 UMI read.
pub fn umi_read(
  id : String,
  sequence : String,
  qualities : Array[Int32],
  umi : String,
  chrom : String,
  pos : Int32,
  strand~ : Int32 = 1,
) -> UmiRead raise UmiError {
  if sequence.length() == 0 {
    raise UmiError("read `\{id}` has empty sequence")
  }
  if qualities.length() != sequence.length() {
    raise UmiError(
      "read `\{id}` quality length \{qualities.length()} != sequence length \{sequence.length()}",
    )
  }
  if !umi_is_valid(umi) {
    raise UmiError("read `\{id}` has invalid UMI `\{umi}`")
  }
  { id, sequence, qualities, umi: umi.to_upper(), chrom, pos, strand, }
}

///|
/// 从 read 5' 端剪除 `umi_len` 个碱基.
pub fn umi_clip_read(
  sequence : String,
  qualities : Array[Int32],
  umi_len : Int32,
) -> (String, Array[Int32]) raise UmiError {
  if umi_len < 0 || umi_len > sequence.length() {
    raise UmiError(
      "cannot clip \{umi_len} bases from read of length \{sequence.length()}",
    )
  }
  if qualities.length() != sequence.length() {
    raise UmiError("quality and sequence length mismatch")
  }
  let seq = sequence[umi_len:].to_owned()
  let qual : Array[Int32] = Array(capacity=sequence.length() - umi_len)
  for i = umi_len; i < sequence.length(); i = i + 1 {
    qual.push(qualities[i])
  }
  (seq, qual)
}

///|
/// UMI 编辑距离(Hamming, 要求等长; 非等长返回较大值).
pub fn umi_edit_distance(a : String, b : String) -> Int32 {
  if a.length() != b.length() {
    return 99
  }
  let mut d = 0
  for i = 0; i < a.length(); i = i + 1 {
    if a.unsafe_get(i) != b.unsafe_get(i) {
      d = d + 1
    }
  }
  d
}

///|
/// 碱基编码: A=0 C=1 G=2 T=3 N=4 other=5.
fn umi_base_index(code : Int32) -> Int32 {
  match code {
    65 | 97 => 0
    67 | 99 => 1
    71 | 103 => 2
    84 | 116 => 3
    78 | 110 => 4
    _ => 5
  }
}

///|
/// 编码转碱基字符.
fn umi_base_char(idx : Int32) -> Char {
  match idx {
    0 => 'A'
    1 => 'C'
    2 => 'G'
    3 => 'T'
    4 => 'N'
    _ => 'N'
  }
}

///|
/// 质量加权逐位多数投票.
///
/// 权重为 `phred + 1`(负质量按 0 处理), 平票按 A/C/G/T/N/other 的
/// 固定顺序决定, 保证结果可复现.
pub fn umi_consensus(
  sequences : Array[String],
  qualities : Array[Array[Int32]],
) -> String raise UmiError {
  if sequences.length() == 0 {
    raise UmiError("cannot build consensus from zero reads")
  }
  if sequences.length() != qualities.length() {
    raise UmiError("sequence and quality collection size mismatch")
  }
  let len = sequences[0].length()
  if len == 0 {
    raise UmiError("cannot build consensus from empty sequences")
  }
  for s in sequences {
    if s.length() != len {
      raise UmiError("consensus input sequences have inconsistent lengths")
    }
  }
  for q in qualities {
    if q.length() != len {
      raise UmiError("consensus input quality length mismatch")
    }
  }
  let sb = StringBuilder()
  for pos = 0; pos < len; pos = pos + 1 {
    let weights : Array[Double] = Array::make(6, 0.0)
    for ri = 0; ri < sequences.length(); ri = ri + 1 {
      let idx = umi_base_index(sequences[ri].unsafe_get(pos).to_int())
      let q = if qualities[ri][pos] > 0 { qualities[ri][pos] } else { 0 }
      weights[idx] = weights[idx] + q.to_double() + 1.0
    }
    let mut best = 0
    for b = 1; b < 6; b = b + 1 {
      if weights[b] > weights[best] {
        best = b
      }
    }
    sb.write_char(umi_base_char(best))
  }
  sb.to_string()
}

///|
/// 一个 UMI family.
pub(all) struct UmiFamily {
  umi : String
  chrom : String
  pos : Int32
  strand : Int32
  mut size : Int32
  consensus : String
  mut mean_quality : Double
}

///|
/// 按 `(UMI, chrom, pos, strand)` 分组并生成共识.
pub fn umi_dedup(
  reads : Array[UmiRead],
  min_family_size~ : Int32 = 1,
) -> Array[UmiFamily] raise UmiError {
  let groups : Map[String, Array[Int32]] = Map([], capacity=reads.length())
  for i = 0; i < reads.length(); i = i + 1 {
    let r = reads[i]
    let key = "\{r.umi}|\{r.chrom}|\{r.pos.to_string()}|\{r.strand.to_string()}"
    match groups.get(key) {
      Some(list) => list.push(i)
      None => groups[key] = [i]
    }
  }
  let out : Array[UmiFamily] = Array(capacity=0)
  for _, indices in groups {
    let seqs : Array[String] = Array(capacity=indices.length())
    let quals : Array[Array[Int32]] = Array(capacity=indices.length())
    let mut q_sum = 0.0
    let mut q_n = 0
    for i in indices {
      seqs.push(reads[i].sequence)
      quals.push(reads[i].qualities)
      for q in reads[i].qualities {
        q_sum = q_sum + q.to_double()
        q_n = q_n + 1
      }
    }
    if indices.length() >= min_family_size {
      let first = reads[indices[0]]
      out.push({
        umi: first.umi,
        chrom: first.chrom,
        pos: first.pos,
        strand: first.strand,
        size: indices.length(),
        consensus: umi_consensus(seqs, quals),
        mean_quality: if q_n > 0 {
          q_sum / q_n.to_double()
        } else {
          0.0
        },
      })
    }
  }
  out.sort_by(fn(a : UmiFamily, b : UmiFamily) -> Int32 {
    let c = a.chrom.compare(b.chrom)
    if c != 0 {
      return c
    }
    if a.pos != b.pos {
      return a.pos - b.pos
    }
    if a.strand != b.strand {
      return a.strand - b.strand
    }
    a.umi.compare(b.umi)
  })
  out
}

///|
/// directional UMI 纠错: 同一 locus 内, 与更强代表 UMI 编辑距离
/// `<= max_edit_distance` 的 family 并入该代表.
pub fn umi_directional_correct(
  families : Array[UmiFamily],
  max_edit_distance~ : Int32 = 1,
) -> Array[UmiFamily] {
  let sorted = families.copy()
  sorted.sort_by(fn(a : UmiFamily, b : UmiFamily) -> Int32 {
    let c = a.chrom.compare(b.chrom)
    if c != 0 {
      return c
    }
    if a.pos != b.pos {
      return a.pos - b.pos
    }
    if a.strand != b.strand {
      return a.strand - b.strand
    }
    if a.size != b.size {
      return b.size - a.size
    }
    a.umi.compare(b.umi)
  })
  let reps : Array[UmiFamily] = Array(capacity=sorted.length())
  for fam in sorted {
    let mut merged = false
    for rep in reps {
      if rep.chrom == fam.chrom &&
        rep.pos == fam.pos &&
        rep.strand == fam.strand &&
        umi_edit_distance(rep.umi, fam.umi) <= max_edit_distance {
        let total = rep.size + fam.size
        rep.mean_quality = if total > 0 {
          (
            rep.mean_quality * rep.size.to_double() +
            fam.mean_quality * fam.size.to_double()
          ) /
          total.to_double()
        } else {
          0.0
        }
        rep.size = total
        merged = true
        break
      }
    }
    if !merged {
      reps.push({
        umi: fam.umi,
        chrom: fam.chrom,
        pos: fam.pos,
        strand: fam.strand,
        size: fam.size,
        consensus: fam.consensus,
        mean_quality: fam.mean_quality,
      })
    }
  }
  reps
}

///|
/// UMI 汇总统计.
pub(all) struct UmiStats {
  n_families : Int32
  n_reads : Int32
  singleton_families : Int32
  max_family_size : Int32
  mean_family_size : Double
  mean_quality : Double
}

///|
/// 统计 families.
pub fn umi_stats(families : Array[UmiFamily]) -> UmiStats {
  let mut reads = 0
  let mut singletons = 0
  let mut max_size = 0
  let mut q_sum = 0.0
  let mut q_n = 0
  for f in families {
    reads = reads + f.size
    if f.size == 1 {
      singletons = singletons + 1
    }
    if f.size > max_size {
      max_size = f.size
    }
    q_sum = q_sum + f.mean_quality
    q_n = q_n + 1
  }
  {
    n_families: families.length(),
    n_reads: reads,
    singleton_families: singletons,
    max_family_size: max_size,
    mean_family_size: if q_n > 0 {
      reads.to_double() / q_n.to_double()
    } else {
      0.0
    },
    mean_quality: if q_n > 0 {
      q_sum / q_n.to_double()
    } else {
      0.0
    },
  }
}

///|
/// families 序列化为 JSON.
pub fn umi_families_to_json(families : Array[UmiFamily]) -> String {
  let arr : Array[JsonValue] = Array(capacity=families.length())
  for f in families {
    let m : Map[String, JsonValue] = Map([], capacity=7)
    m["umi"] = JsonStr(f.umi)
    m["chrom"] = JsonStr(f.chrom)
    m["pos"] = JsonNum(f.pos.to_double())
    m["strand"] = JsonNum(f.strand.to_double())
    m["size"] = JsonNum(f.size.to_double())
    m["consensus"] = JsonStr(f.consensus)
    m["mean_quality"] = JsonNum(f.mean_quality)
    arr.push(JsonObj(m))
  }
  let root : Map[String, JsonValue] = Map([], capacity=2)
  root["n_families"] = JsonNum(families.length().to_double())
  root["families"] = JsonArr(arr)
  json_to_string(JsonObj(root))
}
