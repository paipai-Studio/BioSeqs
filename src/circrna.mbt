///|
/// circrna.mbt — find_circ / CIRI 风格环状 RNA 检测
/// (反向剪接 junction, 本地确定性实现).
///
/// 语义(简化 find_circ 锚点对):
/// - 输入 SAM 的 split-read: 主定位 + SA 补充定位(标准 SA 标签
///   "rname,pos,strand,CIGAR,mapQ,NM");
/// - 反向剪接(BSJ)读段判定: SA 补充定位整体位于主定位上游
///   (SA 参考终点 <= 主参考起点), 即读段跨过 junction 后基因组
///   坐标"倒退"——环状 RNA 的拼接特征(线性融合/剪接为正向);
/// - circRNA 坐标(0-based 半开): 环 = [sa_start, main_end)
///   (SA 定位起点 = acceptor 处环头第一碱基, 主定位终点 = donor
///   位点), BSJ 断点 = (main_end, sa_start), 环长 = main_end -
///   sa_start;
/// - 聚合同 (chrom, sa_start, main_end) 的读段数为 BSJ 支持数;
/// - 过滤: 支持数 < min_support 或环长 < min_length 的丢弃;
/// - 链方向: 主读段与 SA 链一致时标该链, 不一致标 ".".
///
/// 已知简化: 不做 CIRI 的穷举种子扩展(仅 SA split-read),
/// 不区分 Alu 介导与非 Alu 环.
/// 数据不出本机.

///|
/// 一条检测到的 circRNA.
pub(all) struct CircRnaRecord {
  chrom : String
  /// 0-based 半开环区间.
  start : Int32
  end_ : Int32
  /// 环长(bp).
  length : Int32
  /// BSJ 支持读段数.
  support : Int32
  /// "+"/"-"/"."(主与 SA 链不一致).
  strand : String
  /// "circ_\{chrom}_\{start}_\{end}" 风格名称.
  name : String
}

///|
/// 检测参数.
pub(all) struct CircRnaParams {
  /// 最小 BSJ 支持读段数.
  min_support : Int32
  /// 最小环长(bp).
  min_length : Int32
}

///|
/// 默认参数.
pub fn circrna_params() -> CircRnaParams {
  { min_support: 2, min_length: 50, }
}

// ---------------------------------------------------------------------------
// SA 解析与坐标
// ---------------------------------------------------------------------------

///|
/// 解析 SA 标签: "rname,pos,strand,CIGAR,mapQ,NM;..." ->
/// (rname, pos(1-based), strand("+"/"-"), cigar) 列表.
fn circrna_parse_sa(
  sa_value : String,
) -> Array[(String, Int32, String, String)] {
  let entries : Array[(String, Int32, String, String)] = Array(capacity=0)
  let parts = sa_value.split(";").to_array()
  for part in parts {
    let p = part.to_owned()
    if p.length() == 0 {
      continue
    }
    let fields = p.split(",").to_array()
    if fields.length() < 4 {
      continue
    }
    let rname = fields[0].to_owned()
    let pos = parse_int(fields[1].to_owned())
    let strand = fields[2].to_owned()
    let cigar = fields[3].to_owned()
    entries.push((rname, pos, strand, cigar))
  }
  entries
}

///|
/// CIGAR 参考跨度(M/D/N/=/X 消耗参考).
fn circrna_cigar_span(cigar : String) -> Int32 {
  let mut span = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar.unsafe_get(i).unsafe_to_char()
    if c >= '0' && c <= '9' {
      num = num * 10 + (c.to_int() - '0'.to_int())
    } else {
      if c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X' {
        span = span + num
      }
      num = 0
    }
  }
  span
}

// ---------------------------------------------------------------------------
// 检测
// ---------------------------------------------------------------------------

///|
/// 从 SAM 文本检测 circRNA(反向剪接 split-read 聚合).
pub fn circrna_detect(
  sam_text : String,
  params~ : CircRnaParams = circrna_params(),
) -> Array[CircRnaRecord] {
  let sam = parse_sam(sam_text)
  // 聚合键 "chrom:start:end" -> (support, strand)
  let table : Map[String, (Int32, String)] = Map([], capacity=0)
  let coords : Map[String, (String, Int32, Int32)] = Map([], capacity=0)
  for rec in sam.records {
    // 仅主比对(secondary/supplementary 跳过, 避免双计)
    if (rec.flag & 0x100) != 0 || (rec.flag & 0x800) != 0 {
      continue
    }
    let sa = rec.tags.get("SA")
    match sa {
      None => continue
      Some(sa_value) => {
        let supp = circrna_parse_sa(sa_value)
        if supp.length() == 0 {
          continue
        }
        for entry in supp {
          let (sa_chrom, sa_pos, sa_strand, sa_cigar) = entry
          if sa_chrom != rec.rname {
            continue
          }
          // 主定位参考区间(0-based 半开)
          let main_start = rec.pos - 1
          let main_end = main_start + circrna_cigar_span_cigar(rec.cigar)
          // SA 定位参考区间
          let sa_start = sa_pos - 1
          let sa_end = sa_start + circrna_cigar_span(sa_cigar)
          // BSJ 判定: SA 整体在主定位上游(坐标倒退)
          if sa_end > main_start {
            continue
          }
          // 环区间 [sa_start, main_end)
          let circ_start = sa_start
          let circ_end = main_end
          let circ_len = circ_end - circ_start
          if circ_len < params.min_length {
            continue
          }
          // 链方向: 主 flag 0x10 与 SA strand
          let main_strand = if (rec.flag & 0x10) != 0 { "-" } else { "+" }
          let strand = if main_strand == sa_strand { main_strand } else { "." }
          let key = "\{sa_chrom}:\{circ_start}:\{circ_end}"
          let (cnt0, _) = table.get(key).unwrap_or((0, strand))
          table.set(key, (cnt0 + 1, strand))
          coords.set(key, (sa_chrom, circ_start, circ_end))
        }
      }
    }
  }
  // 输出(过滤 + 按 (chrom, start) 排序)
  let out : Array[CircRnaRecord] = Array(capacity=table.length())
  for key in table.keys() {
    let (support, strand) = table.get(key).unwrap()
    if support < params.min_support {
      continue
    }
    let (chrom, start, end_) = coords.get(key).unwrap()
    out.push(CircRnaRecord::{
      chrom,
      start,
      end_,
      length: end_ - start,
      support,
      strand,
      name: "circ_\{chrom}_\{start}_\{end_}",
    })
  }
  out.sort_by(fn(a, b) {
    if a.chrom != b.chrom {
      if a.chrom < b.chrom {
        -1
      } else {
        1
      }
    } else if a.start != b.start {
      a.start - b.start
    } else {
      0
    }
  })
  out
}

///|
/// 从 CigarElement 数组计算参考跨度.
fn circrna_cigar_span_cigar(cigar : Array[CigarElement]) -> Int32 {
  let mut span = 0
  for elem in cigar {
    let op = elem.op
    if op is Match || op is Delete || op is Skip || op is Equal || op is Diff {
      span = span + elem.length
    }
  }
  span
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// circRNA 记录 -> TSV(find_circ 风格: chrom start end name
/// support strand length).
pub fn circrna_tsv(records : Array[CircRnaRecord]) -> String {
  let sb = StringBuilder()
  sb.write_string("chrom\tstart\tend\tname\tsupport\tstrand\tlength\n")
  for r in records {
    sb.write_string(r.chrom)
    sb.write_char('\t')
    sb.write_string(r.start.to_string())
    sb.write_char('\t')
    sb.write_string(r.end_.to_string())
    sb.write_char('\t')
    sb.write_string(r.name)
    sb.write_char('\t')
    sb.write_string(r.support.to_string())
    sb.write_char('\t')
    sb.write_string(r.strand)
    sb.write_char('\t')
    sb.write_string(r.length.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// circRNA 记录 -> JSON 值.
pub fn circrna_value(records : Array[CircRnaRecord]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=records.length())
  for r in records {
    arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(r.chrom)),
          ("start", JsonNum(r.start.to_double())),
          ("end", JsonNum(r.end_.to_double())),
          ("name", JsonStr(r.name)),
          ("support", JsonNum(r.support.to_double())),
          ("strand", JsonStr(r.strand)),
          ("length", JsonNum(r.length.to_double())),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// circRNA 记录 -> JSON 文本.
pub fn circrna_json(records : Array[CircRnaRecord]) -> String {
  json_to_string(circrna_value(records))
}

///|
/// 摘要行.
pub fn circrna_summary_line(records : Array[CircRnaRecord]) -> String {
  let mut total_support = 0
  for r in records {
    total_support = total_support + r.support
  }
  "circrna: records=\{records.length()}, total_bsj_support=\{total_support}"
}

///|
/// 单条记录摘要.
pub fn CircRnaRecord::summary(self : CircRnaRecord) -> String {
  "\{self.name} [\{self.chrom}:\{self.start}-\{self.end_}] len=\{self.length} support=\{self.support} strand=\{self.strand}"
}
