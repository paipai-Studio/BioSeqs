///|
/// BEDPE 交互对格式读写与分析.
///
/// BEDPE(UCSC/bedtools)描述两个基因组区间之间的配对关系, 每行 10+ 列:
/// `chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 [fields...]`
/// 坐标为 0-based 半开区间(与 BED 一致); 列分隔符兼容 tab 与连续空格.
/// 典型来源: Hi-C loops、结构变异(SV)、增强子-启动子配对、ChIA-PET.
pub suberror BedpeError {
  BedpeError(String)
}

///|
/// 一条 BEDPE 交互记录.
pub struct BedpePair {
  chrom1 : String
  start1 : Int32
  end1 : Int32
  chrom2 : String
  start2 : Int32
  end2 : Int32
  name : String
  score : Int32
  strand1 : String
  strand2 : String
  fields : Array[String]
}

///|
/// 构造 BEDPE 记录(附加列默认为空).
pub fn BedpePair::new(
  chrom1 : String,
  start1 : Int32,
  end1 : Int32,
  chrom2 : String,
  start2 : Int32,
  end2 : Int32,
  name~ : String = ".",
  score~ : Int32 = 0,
  strand1~ : String = ".",
  strand2~ : String = ".",
  fields~ : Array[String] = [],
) -> BedpePair {
  {
    chrom1,
    start1,
    end1,
    chrom2,
    start2,
    end2,
    name,
    score,
    strand1,
    strand2,
    fields,
  }
}

///|
/// 第一个锚点长度.
pub fn BedpePair::anchor_length1(self : BedpePair) -> Int32 {
  self.end1 - self.start1
}

///|
/// 第二个锚点长度.
pub fn BedpePair::anchor_length2(self : BedpePair) -> Int32 {
  self.end2 - self.start2
}

///|
/// 两锚点是否位于同一染色体.
pub fn BedpePair::is_same_chrom(self : BedpePair) -> Bool {
  self.chrom1 == self.chrom2
}

///|
/// Loop 长度(两锚点内侧距离): 同染色体时为
/// max(start1, start2) - min(end1, end2), 锚点重叠时为 0;
/// 跨染色体返回 -1.
pub fn BedpePair::loop_length(self : BedpePair) -> Int32 {
  if self.chrom1 != self.chrom2 {
    return -1
  }
  let inner = if self.start1 > self.start2 { self.start1 } else { self.start2 }
  let outer = if self.end1 < self.end2 { self.end1 } else { self.end2 }
  if inner > outer {
    inner - outer
  } else {
    0
  }
}

///|
/// 转为一对 BED 区间(name 各带 _1/_2 后缀以保持关联).
pub fn BedpePair::to_bed(self : BedpePair) -> (BedInterval, BedInterval) {
  (
    BedInterval::new(self.chrom1, self.start1, self.end1),
    BedInterval::new(self.chrom2, self.start2, self.end2),
  )
}

///|
/// 任一锚点与区间 [start, end) 在染色体 chrom 上相交.
pub fn BedpePair::overlaps_region(
  self : BedpePair,
  chrom : String,
  start : Int32,
  end : Int32,
) -> Bool {
  if self.chrom1 == chrom && self.start1 < end && self.end1 > start {
    return true
  }
  if self.chrom2 == chrom && self.start2 < end && self.end2 > start {
    return true
  }
  false
}

///|
/// 按连续空白(tab/空格)切分字段, 丢弃空字段.
fn bedpe_split_fields(line : StringView) -> Array[String] {
  let fields : Array[String] = []
  let mut sb = StringBuilder()
  for c in line {
    if c == '\t' || c == ' ' {
      let cur = sb.to_string()
      if cur.length() > 0 {
        fields.push(cur)
      }
      sb = StringBuilder()
    } else {
      sb.write_char(c)
    }
  }
  let cur = sb.to_string()
  if cur.length() > 0 {
    fields.push(cur)
  }
  fields
}

///|
/// 严格解析 Int32(全数字, 允许负号), 失败抛 BedpeError.
fn bedpe_parse_int(
  s : String,
  field : String,
  line_no : Int32,
) -> Int32 raise BedpeError {
  let n = s.length()
  if n == 0 {
    raise BedpeError("Empty " + field + " at BEDPE line " + line_no.to_string())
  }
  let mut i = 0
  if s.unsafe_get(0) == '-' {
    i = 1
    if n == 1 {
      raise BedpeError(
        "Invalid " + field + " '" + s + "' at BEDPE line " + line_no.to_string(),
      )
    }
  }
  let zero = '0'.to_int()
  let nine = '9'.to_int()
  let mut value : Int64 = 0
  let neg = s.unsafe_get(0) == '-'
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > nine {
      raise BedpeError(
        "Invalid " + field + " '" + s + "' at BEDPE line " + line_no.to_string(),
      )
    }
    value = value * 10L + (c - zero).to_int64()
    i = i + 1
  }
  if neg {
    value = -value
  }
  if value > 2147483647L || value < -2147483648L {
    raise BedpeError(
      field + " out of range '" + s + "' at BEDPE line " + line_no.to_string(),
    )
  }
  value.to_int()
}

///|
/// 校验链向字段(+/-.).
fn bedpe_validate_strand(
  s : String,
  field : String,
  line_no : Int32,
) -> Unit raise BedpeError {
  if s != "+" && s != "-" && s != "." {
    raise BedpeError(
      "Invalid " +
      field +
      " '" +
      s +
      "' (expected +/-/.) at BEDPE line " +
      line_no.to_string(),
    )
  }
}

///|
/// 解析 BEDPE 文本. 跳过空行与 '#'/'track'/'browser' 注释行;
/// 至少需要 6 列(两个区间), name/score/strand 缺省为 '.'/0/'.',
/// 第 11 列起原样保留在 fields. score 要求 0-1000.
pub fn parse_bedpe(content : String) -> Array[BedpePair] raise BedpeError {
  let pairs : Array[BedpePair] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let mut line_no = 0
  for raw in lines {
    line_no = line_no + 1
    let line = raw.trim()
    if line.length() == 0 {
      continue
    }
    if line.unsafe_get(0) == '#' {
      continue
    }
    if line.has_prefix("track") || line.has_prefix("browser") {
      continue
    }
    let fields = bedpe_split_fields(line)
    if fields.length() < 6 {
      raise BedpeError(
        "Expected at least 6 columns, got " +
        fields.length().to_string() +
        " at BEDPE line " +
        line_no.to_string(),
      )
    }
    let chrom1 = fields[0]
    let start1 = bedpe_parse_int(fields[1], "start1", line_no)
    let end1 = bedpe_parse_int(fields[2], "end1", line_no)
    let chrom2 = fields[3]
    let start2 = bedpe_parse_int(fields[4], "start2", line_no)
    let end2 = bedpe_parse_int(fields[5], "end2", line_no)
    if end1 < start1 || end2 < start2 {
      raise BedpeError(
        "Interval end before start at BEDPE line " + line_no.to_string(),
      )
    }
    let name = if fields.length() > 6 { fields[6] } else { "." }
    let score = if fields.length() > 7 {
      let v = bedpe_parse_int(fields[7], "score", line_no)
      if v < 0 || v > 1000 {
        raise BedpeError(
          "Score " +
          v.to_string() +
          " out of range 0-1000 at BEDPE line " +
          line_no.to_string(),
        )
      }
      v
    } else {
      0
    }
    let strand1 = if fields.length() > 8 {
      let s = fields[8]
      bedpe_validate_strand(s, "strand1", line_no)
      s
    } else {
      "."
    }
    let strand2 = if fields.length() > 9 {
      let s = fields[9]
      bedpe_validate_strand(s, "strand2", line_no)
      s
    } else {
      "."
    }
    let extra : Array[String] = []
    if fields.length() > 10 {
      for i in 10..<fields.length() {
        extra.push(fields[i])
      }
    }
    pairs.push(
      BedpePair::new(
        chrom1,
        start1,
        end1,
        chrom2,
        start2,
        end2,
        name~,
        score~,
        strand1~,
        strand2~,
        fields=extra,
      ),
    )
  }
  pairs
}

///|
/// 写出 BEDPE 文本(10 标准列 + 附加列, tab 分隔).
pub fn write_bedpe(pairs : Array[BedpePair]) -> String {
  let sb = StringBuilder()
  for p in pairs {
    sb.write_string(
      p.chrom1 +
      "\t" +
      p.start1.to_string() +
      "\t" +
      p.end1.to_string() +
      "\t" +
      p.chrom2 +
      "\t" +
      p.start2.to_string() +
      "\t" +
      p.end2.to_string() +
      "\t" +
      p.name +
      "\t" +
      p.score.to_string() +
      "\t" +
      p.strand1 +
      "\t" +
      p.strand2,
    )
    for f in p.fields {
      sb.write_char('\t')
      sb.write_string(f)
    }
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 过滤同染色体内的交互对.
pub fn bedpe_filter_same_chrom(pairs : Array[BedpePair]) -> Array[BedpePair] {
  let out : Array[BedpePair] = Array(capacity=0)
  for p in pairs {
    if p.is_same_chrom() {
      out.push(p)
    }
  }
  out
}

///|
/// 按染色体对计数交互(键为 "chrom1|chrom2", 保持输入顺序首次出现).
pub fn bedpe_count_chrom_pairs(pairs : Array[BedpePair]) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([], capacity=0)
  for p in pairs {
    let key = p.chrom1 + "|" + p.chrom2
    match counts.get(key) {
      Some(v) => counts[key] = v + 1
      None => counts[key] = 1
    }
  }
  counts
}

///|
/// Loop 长度分布直方图(仅同染色体对): [bin_start, bin_start+bin_size)
/// 计数, 长度超出 (bin_start + bin_size*nbins) 的计入最后一个桶;
/// 返回每桶计数.
pub fn bedpe_distance_histogram(
  pairs : Array[BedpePair],
  bin_start : Int32,
  bin_size : Int32,
  nbins : Int32,
) -> Array[Int32] {
  let hist = Array::make(nbins, 0)
  if bin_size <= 0 || nbins <= 0 {
    return hist
  }
  for p in pairs {
    let d = p.loop_length()
    if d < 0 {
      continue
    }
    if d < bin_start {
      continue
    }
    let q = (d - bin_start) / bin_size
    let idx : Int32 = if q < 0 { 0 } else if q >= nbins { nbins - 1 } else { q }
    hist[idx] = hist[idx] + 1
  }
  hist
}

///|
/// 判断两条交互对是否存在任一锚点组合相交(锚点区间半开).
pub fn bedpe_pair_overlaps(a : BedpePair, b : BedpePair) -> Bool {
  let aa1 = (a.chrom1, a.start1, a.end1)
  let aa2 = (a.chrom2, a.start2, a.end2)
  let ab1 = (b.chrom1, b.start1, b.end1)
  let ab2 = (b.chrom2, b.start2, b.end2)
  let overlap = (x : (String, Int32, Int32), y : (String, Int32, Int32)) => {
    let (cx, sx, ex) = x
    let (cy, sy, ey) = y
    cx == cy && sx < ey && ex > sy
  }
  overlap(aa1, ab1) ||
  overlap(aa1, ab2) ||
  overlap(aa2, ab1) ||
  overlap(aa2, ab2)
}
