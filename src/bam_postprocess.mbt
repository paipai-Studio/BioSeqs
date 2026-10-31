///|
/// bam_postprocess.mbt — BAM 比对后处理(本地确定性实现):
/// - sam_to_bam: SAM 文本 -> 内存 BamFile 桥接(供后处理与测试构造);
/// - bam_fixmate: samtools fixmate 语义(qname 相邻配对互填 mate 坐标/
///   TLEN 外坐标跨度/镜像 mate flag);
/// - bam_markdup_paired: Picard MarkDuplicates 语义(配对 5'-coordinate
///   mate-key 分组, mapq 最优保留其余打 0x400, 光学重复坐标判别);
/// - bam_bqsr: GATK BaseRecalibrator/ApplyBQSR 简化(covariate 表
///   (quality x 三碱基上下文) 经验错误率重估, known-sites 位点豁免,
///   重写碱基质量).
///
/// 数据不出本机.

// ---------------------------------------------------------------------------
// SAM -> BAM 桥接
// ---------------------------------------------------------------------------

///|
/// SAM CIGAR 操作符 -> BAM CIGAR 操作符(两套同构枚举).
fn bam_pp_map_cigar_op(op : CigarOp) -> BamCigarOp {
  match op {
    Match => BamCigarOp::Match
    Insert => BamCigarOp::Insert
    Delete => BamCigarOp::Delete
    Skip => BamCigarOp::Skip
    SoftClip => BamCigarOp::SoftClip
    HardClip => BamCigarOp::HardClip
    Pad => BamCigarOp::Pad
    Equal => BamCigarOp::Equal
    Diff => BamCigarOp::Diff
  }
}

///|
/// htslib reg2bin 标准算法(0-based 半开区间).
fn bam_pp_reg2bin(beg : Int32, end : Int32) -> Int32 {
  let e = end - 1
  if beg >> 14 == e >> 14 {
    return 4681 + (beg >> 14)
  }
  if beg >> 17 == e >> 17 {
    return 585 + (beg >> 17)
  }
  if beg >> 20 == e >> 20 {
    return 73 + (beg >> 20)
  }
  if beg >> 23 == e >> 23 {
    return 9 + (beg >> 23)
  }
  if beg >> 26 == e >> 26 {
    return 1 + (beg >> 26)
  }
  0
}

///|
/// 数字字符判定.
fn bam_pp_is_digit(c : Char) -> Bool {
  c == '0' ||
  c == '1' ||
  c == '2' ||
  c == '3' ||
  c == '4' ||
  c == '5' ||
  c == '6' ||
  c == '7' ||
  c == '8' ||
  c == '9'
}

///|
/// 整数字符串判定(可带负号).
fn bam_pp_is_int_str(s : String) -> Bool {
  if s.length() == 0 {
    return false
  }
  let start : Int32 = if s.unsafe_get(0).unsafe_to_char() == '-' {
    1
  } else {
    0
  }
  if s.length() <= start {
    return false
  }
  for i = start; i < s.length(); i = i + 1 {
    if bam_pp_is_digit(s.unsafe_get(i).unsafe_to_char()) == false {
      return false
    }
  }
  true
}

///|
/// 浮点字符串判定(整数部分 + 恰一个小数点).
fn bam_pp_is_float_str(s : String) -> Bool {
  let mut has_dot = false
  let mut digits = 0
  let start : Int32 = if s.length() > 0 &&
    s.unsafe_get(0).unsafe_to_char() == '-' {
    1
  } else {
    0
  }
  if s.length() <= start {
    return false
  }
  for i = start; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c == '.' {
      if has_dot {
        return false
      }
      has_dot = true
    } else if bam_pp_is_digit(c) {
      digits = digits + 1
    } else {
      return false
    }
  }
  has_dot && digits > 0
}

///|
/// SAM tag 值(纯值字符串) -> BamTagValue(类型推断: 整数/浮点/单字符/串).
fn bam_pp_tag_value(value : String) -> (Char, BamTagValue) {
  if bam_pp_is_int_str(value) {
    return ('i', BamTagValue::Int(parse_int(value)))
  }
  if bam_pp_is_float_str(value) {
    return ('f', BamTagValue::Float(parse_float_str(value).unwrap_or(0.0)))
  }
  if value.length() == 1 {
    return ('A', BamTagValue::Char(value.unsafe_get(0).unsafe_to_char()))
  }
  ('Z', BamTagValue::Str(value))
}

///|
fn bam_pp_header_type_name(t : SamHeaderType) -> String {
  match t {
    HD => "HD"
    SQ => "SQ"
    RG => "RG"
    PG => "PG"
    CO => "CO"
  }
}

///|
/// SAM 文本 -> 内存 BamFile.
/// SQ 头行建立参考索引; 坐标 1-based -> 0-based; rnext "=" 复用本条
/// ref_id; qual "*" 转空串; tag 值按内容推断类型.
pub fn sam_to_bam(sam_text : String) -> BamFile {
  let sam = parse_sam(sam_text)
  // SQ 头 -> 参考表
  let ref_names : Array[String] = Array(capacity=0)
  let ref_lengths : Array[Int32] = Array(capacity=0)
  let header_text = StringBuilder()
  for h in sam.headers {
    if h.header_type is SQ {
      let name = h.tags.get("SN").unwrap_or("")
      let ln = parse_int(h.tags.get("LN").unwrap_or("0"))
      if name != "" {
        ref_names.push(name)
        ref_lengths.push(ln)
      }
    }
    header_text.write_string("@")
    header_text.write_string(bam_pp_header_type_name(h.header_type))
    header_text.write_char('\n')
  }
  let refs : FixedArray[BamReference] = FixedArray::make(ref_names.length(), BamReference::{
    name: "",
    length: 0,
  })
  for i = 0; i < ref_names.length(); i = i + 1 {
    refs[i] = BamReference::{ name: ref_names[i], length: ref_lengths[i], }
  }
  let header = BamHeader::{ text: header_text.to_string(), references: refs, }
  // 记录转换
  let records : Array[BamRecord] = Array(capacity=sam.records.length())
  for rec in sam.records {
    let ref_id : Int32 = {
      let mut idx = -1
      for i = 0; i < ref_names.length(); i = i + 1 {
        if ref_names[i] == rec.rname {
          idx = i
          break
        }
      }
      idx
    }
    let pos : Int32 = if rec.pos > 0 { rec.pos - 1 } else { -1 }
    let next_ref_id : Int32 = if rec.rnext == "=" {
      ref_id
    } else if rec.rnext == "*" || rec.rnext == "" {
      -1
    } else {
      let mut idx = -1
      for i = 0; i < ref_names.length(); i = i + 1 {
        if ref_names[i] == rec.rnext {
          idx = i
          break
        }
      }
      idx
    }
    let next_pos : Int32 = if rec.pnext > 0 { rec.pnext - 1 } else { -1 }
    let cigar : FixedArray[BamCigarElement] = FixedArray::make(
      rec.cigar.length(),
      BamCigarElement::{ op: BamCigarOp::Match, length: 0, },
    )
    for ci = 0; ci < rec.cigar.length(); ci = ci + 1 {
      cigar[ci] = BamCigarElement::{
        op: bam_pp_map_cigar_op(rec.cigar[ci].op),
        length: rec.cigar[ci].length,
      }
    }
    let tags : Array[BamTag] = Array(capacity=0)
    for key in rec.tags.keys() {
      let value = rec.tags.get(key).unwrap_or("")
      let (code, v) = bam_pp_tag_value(value)
      tags.push(BamTag::{ name: key, type_code: code, value: v, })
    }
    let end : Int32 = if pos >= 0 && cigar.length() > 0 {
      let mut span = 0
      for ci = 0; ci < cigar.length(); ci = ci + 1 {
        let op = cigar[ci].op
        if op is Match ||
          op is Delete ||
          op is Skip ||
          op is Equal ||
          op is Diff {
          span = span + cigar[ci].length
        }
      }
      pos + span
    } else {
      0
    }
    let bin : Int32 = if pos >= 0 { bam_pp_reg2bin(pos, end) } else { 4680 }
    let qual : String = if rec.qual == "*" { "" } else { rec.qual }
    records.push(BamRecord::{
      ref_id,
      pos,
      bin,
      mapq: rec.mapq,
      flag: rec.flag,
      l_read_name: rec.qname.length() + 1,
      cigar_len: cigar.length(),
      l_seq: rec.seq.length(),
      next_ref_id,
      next_pos,
      tlen: rec.tlen,
      qname: rec.qname,
      cigar,
      seq: rec.seq,
      qual,
      tags,
    })
  }
  BamFile::{ header, records, }
}

// ---------------------------------------------------------------------------
// fixmate
// ---------------------------------------------------------------------------

///|
/// 构造替换部分字段的新记录.
fn bam_pp_with_mate(
  rec : BamRecord,
  next_ref_id~ : Int32? = None,
  next_pos~ : Int32? = None,
  tlen~ : Int32? = None,
  flag~ : Int32? = None,
  qual~ : String? = None,
) -> BamRecord {
  BamRecord::{
    ref_id: rec.ref_id,
    pos: rec.pos,
    bin: rec.bin,
    mapq: rec.mapq,
    flag: flag.unwrap_or(rec.flag),
    l_read_name: rec.l_read_name,
    cigar_len: rec.cigar_len,
    l_seq: rec.l_seq,
    next_ref_id: next_ref_id.unwrap_or(rec.next_ref_id),
    next_pos: next_pos.unwrap_or(rec.next_pos),
    tlen: tlen.unwrap_or(rec.tlen),
    qname: rec.qname,
    cigar: rec.cigar,
    seq: rec.seq,
    qual: qual.unwrap_or(rec.qual),
    tags: rec.tags,
  }
}

///|
/// 读段参考跨度终点(0-based 半开).
fn bam_pp_end(rec : BamRecord) -> Int32 {
  let mut span = 0
  for ci = 0; ci < rec.cigar.length(); ci = ci + 1 {
    let op = rec.cigar[ci].op
    if op is Match || op is Delete || op is Skip || op is Equal || op is Diff {
      span = span + rec.cigar[ci].length
    }
  }
  rec.pos + span
}

///|
/// samtools fixmate 语义: 输入按 qname 相邻排序的配对读段,
/// 互填 mate 坐标(next_ref_id/next_pos), 镜像 mate flag(0x8 mate
/// unmapped / 0x20 mate reverse), 并按外坐标跨度重算 TLEN:
/// 双端同参考映射时 TLEN = 覆盖外坐标宽度(0-based 半开), 左端读段
/// 为正、右端为负, 坐标相同按 first(0x40) 为正; 其余情形置 0.
/// 非配对读段原样保留.
pub fn bam_fixmate(bam : BamFile) -> BamFile {
  let records = bam.records
  let out : Array[BamRecord] = Array(capacity=records.length())
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    let paired = (r.flag & 0x1) != 0
    let mate_here = paired &&
      i + 1 < records.length() &&
      records[i + 1].qname == r.qname &&
      (records[i + 1].flag & 0x1) != 0
    if mate_here == false {
      out.push(bam_pp_with_mate(r, tlen=Some(0)))
      i = i + 1
      continue
    }
    let m = records[i + 1]
    // 互填 mate 坐标
    let r2 = bam_pp_with_mate(
      r,
      next_ref_id=Some(m.ref_id),
      next_pos=Some(m.pos),
      flag=Some(bam_pp_mirror_flags(r, m)),
      tlen=Some(bam_pp_pair_tlen(r, m)),
    )
    let m2 = bam_pp_with_mate(
      m,
      next_ref_id=Some(r.ref_id),
      next_pos=Some(r.pos),
      flag=Some(bam_pp_mirror_flags(m, r)),
      tlen=Some(bam_pp_pair_tlen(m, r)),
    )
    out.push(r2)
    out.push(m2)
    i = i + 2
  }
  BamFile::{ header: bam.header, records: out, }
}

///|
/// 镜像 mate flag: 0x8 随 mate 是否 unmapped, 0x20 随 mate 0x10.
fn bam_pp_mirror_flags(rec : BamRecord, mate : BamRecord) -> Int32 {
  let mut flag = rec.flag
  if mate.is_unmapped() {
    flag = flag | 0x8
  } else {
    flag = flag - (flag & 0x8)
  }
  if (mate.flag & 0x10) != 0 {
    flag = flag | 0x20
  } else {
    flag = flag - (flag & 0x20)
  }
  flag
}

///|
/// 外坐标 TLEN(rec 视角).
fn bam_pp_pair_tlen(rec : BamRecord, mate : BamRecord) -> Int32 {
  if rec.is_unmapped() || mate.is_unmapped() || rec.ref_id != mate.ref_id {
    return 0
  }
  let r_end = bam_pp_end(rec)
  let m_end = bam_pp_end(mate)
  let lo = if rec.pos < mate.pos { rec.pos } else { mate.pos }
  let hi = if r_end > m_end { r_end } else { m_end }
  let span = hi - lo
  if rec.pos != mate.pos {
    if rec.pos < mate.pos {
      return span
    }
    return -span
  }
  // 坐标相同: first-in-pair 为正
  if (rec.flag & 0x40) != 0 && (rec.flag & 0x80) == 0 {
    return span
  }
  -span
}

// ---------------------------------------------------------------------------
// 配对感知 markdup
// ---------------------------------------------------------------------------

///|
/// markdup 统计.
pub(all) struct BamMarkdupStats {
  /// 参与判重的配对读段组数.
  n_pair_groups : Int32
  n_dup_pairs : Int32
  /// 其中光学重复对数.
  n_optical : Int32
  /// 单端读段组数.
  n_single_groups : Int32
  n_dup_single : Int32
  /// 重复率(被标记读段 / 参与判重读段).
  dup_rate : Double
}

///|
/// 读段 5' 坐标(reverse 取末端).
fn bam_pp_five_prime(rec : BamRecord) -> Int32 {
  if (rec.flag & 0x10) != 0 {
    bam_pp_end(rec) - 1
  } else {
    rec.pos
  }
}

///|
/// qname 光学坐标(Picard 风格 "prefix:lane:tile:x:y" 的 x/y),
/// 无坐标格式返回 None.
fn bam_pp_optical_coords(qname : String) -> (Int32, Int32)? {
  let parts = qname.split(":").to_array()
  if parts.length() < 5 {
    return None
  }
  let x = parts[parts.length() - 2].to_owned()
  let y = parts[parts.length() - 1].to_owned()
  let xi = bam_pp_parse_int_or(x, -1)
  let yi = bam_pp_parse_int_or(y, -1)
  if xi < 0 || yi < 0 {
    return None
  }
  Some((xi, yi))
}

///|
fn bam_pp_parse_int_or(s : String, fallback : Int32) -> Int32 {
  if bam_pp_is_int_str(s) {
    parse_int(s)
  } else {
    fallback
  }
}

///|
/// 配对 mate-key(规范化: 两端 (ref_id, 5', strand) 按字典序).
fn bam_pp_pair_key(rec : BamRecord, mate : BamRecord) -> String {
  let a = "\{rec.ref_id}:\{bam_pp_five_prime(rec)}:\{(rec.flag & 0x10) != 0}"
  let b = "\{mate.ref_id}:\{bam_pp_five_prime(mate)}:\{(mate.flag & 0x10) != 0}"
  if a <= b {
    a + "|" + b
  } else {
    b + "|" + a
  }
}

///|
/// 读段是否参与 markdup(主比对: 非重复标记、非 secondary/supplementary).
fn bam_pp_markdup_eligible(rec : BamRecord) -> Bool {
  rec.is_unmapped() == false &&
  (rec.flag & 0x100) == 0 &&
  (rec.flag & 0x800) == 0
}

///|
/// Picard MarkDuplicates 语义(简化): 配对按规范化 5'-coordinate
/// mate-key 分组(以 qname 索引定位 first/second, 不依赖输入排序),
/// 组内保留双端 mapq 之和最高(tie 取 qname 字典序最小)不打标,
/// 其余两条均打 0x400; 单端读段按 (ref_id, 5', strand) 独立分组
/// 同法处理. 光学重复: qname 含 Picard 坐标后缀(":...:x:y")且组内
/// 与保留者切比雪夫距离 <= optical_distance 时计入 n_optical(默认 100).
pub fn bam_markdup_paired(
  bam : BamFile,
  optical_distance~ : Int32 = 100,
) -> (BamFile, BamMarkdupStats) {
  let records = bam.records
  let marked : Array[Bool] = Array::make(records.length(), false)
  // qname -> first/second 读段索引(不依赖输入排序)
  let first_idx : Map[String, Int32] = Map([], capacity=0)
  let second_idx : Map[String, Int32] = Map([], capacity=0)
  let single_list : Array[Int32] = Array(capacity=0)
  let mut n_total_reads = 0
  for i = 0; i < records.length(); i = i + 1 {
    let r = records[i]
    if bam_pp_markdup_eligible(r) == false {
      continue
    }
    n_total_reads = n_total_reads + 1
    if (r.flag & 0x1) != 0 {
      if (r.flag & 0x80) != 0 {
        if second_idx.contains(r.qname) == false {
          second_idx.set(r.qname, i)
        }
      } else if first_idx.contains(r.qname) == false {
        first_idx.set(r.qname, i)
      }
    } else {
      single_list.push(i)
    }
  }
  // 配对分桶: mate-key -> pair 索引对列表
  let pair_groups : Map[String, Array[(Int32, Int32)]] = Map([], capacity=0)
  for qname in first_idx.keys() {
    match second_idx.get(qname) {
      Some(si) => {
        let fi = first_idx.get(qname).unwrap()
        let gkey = bam_pp_pair_key(records[fi], records[si])
        match pair_groups.get(gkey) {
          Some(arr) => arr.push((fi, si))
          None => pair_groups.set(gkey, [(fi, si)])
        }
      }
      None => ()
    }
  }
  let mut n_pair_groups = 0
  let mut n_dup_pairs = 0
  let mut n_optical = 0
  let mut n_marked_reads = 0
  for gkey in pair_groups.keys() {
    let group = pair_groups.get(gkey).unwrap()
    n_pair_groups = n_pair_groups + 1
    if group.length() <= 1 {
      continue
    }
    // 保留双端 mapq 之和最高(tie: qname 字典序最小)
    let mut best = group[0]
    let mut best_score = records[group[0].0].mapq + records[group[0].1].mapq
    let mut best_qname = records[group[0].0].qname
    for gi = 1; gi < group.length(); gi = gi + 1 {
      let (fi, si) = group[gi]
      let score = records[fi].mapq + records[si].mapq
      let better = score > best_score ||
        (score == best_score && records[fi].qname < best_qname)
      if better {
        best = (fi, si)
        best_score = score
        best_qname = records[fi].qname
      }
    }
    let keep_xy = bam_pp_optical_coords(records[best.0].qname)
    for gi = 0; gi < group.length(); gi = gi + 1 {
      if group[gi].0 == best.0 && group[gi].1 == best.1 {
        continue
      }
      let (fi, si) = group[gi]
      marked[fi] = true
      marked[si] = true
      n_marked_reads = n_marked_reads + 2
      n_dup_pairs = n_dup_pairs + 1
      match (keep_xy, bam_pp_optical_coords(records[fi].qname)) {
        (Some(k), Some(d)) =>
          if optical_distance > 0 &&
            (k.0 - d.0).abs() <= optical_distance &&
            (k.1 - d.1).abs() <= optical_distance {
            n_optical = n_optical + 1
          }
        _ => ()
      }
    }
  }
  // 单端分桶
  let single_groups : Map[String, Array[Int32]] = Map([], capacity=0)
  for i in single_list {
    let r = records[i]
    let key = "\{r.ref_id}:\{bam_pp_five_prime(r)}:\{(r.flag & 0x10) != 0}"
    match single_groups.get(key) {
      Some(arr) => arr.push(i)
      None => single_groups.set(key, [i])
    }
  }
  let mut n_single_groups = 0
  let mut n_dup_single = 0
  for skey in single_groups.keys() {
    let group = single_groups.get(skey).unwrap()
    n_single_groups = n_single_groups + 1
    if group.length() <= 1 {
      continue
    }
    let mut best = group[0]
    let mut best_mapq = records[group[0]].mapq
    let mut best_qname = records[group[0]].qname
    for gi = 1; gi < group.length(); gi = gi + 1 {
      let cand = records[group[gi]]
      let better = cand.mapq > best_mapq ||
        (cand.mapq == best_mapq && cand.qname < best_qname)
      if better {
        best = group[gi]
        best_mapq = cand.mapq
        best_qname = cand.qname
      }
    }
    for gi = 0; gi < group.length(); gi = gi + 1 {
      if group[gi] == best {
        continue
      }
      marked[group[gi]] = true
      n_marked_reads = n_marked_reads + 1
      n_dup_single = n_dup_single + 1
    }
  }
  // 应用标记
  let out : Array[BamRecord] = Array(capacity=records.length())
  for k = 0; k < records.length(); k = k + 1 {
    if marked[k] {
      out.push(bam_pp_with_mate(records[k], flag=Some(records[k].flag | 0x400)))
    } else {
      out.push(records[k])
    }
  }
  let dup_rate : Double = if n_total_reads > 0 {
    n_marked_reads.to_double() / n_total_reads.to_double()
  } else {
    0.0
  }
  let stats = BamMarkdupStats::{
    n_pair_groups,
    n_dup_pairs,
    n_optical,
    n_single_groups,
    n_dup_single,
    dup_rate,
  }
  (BamFile::{ header: bam.header, records: out, }, stats)
}

// ---------------------------------------------------------------------------
// BQSR
// ---------------------------------------------------------------------------

///|
/// BQSR 重校正报告.
pub(all) struct BqsrReport {
  /// 参与训练的碱基观察数(不含 known-sites 豁免).
  n_observations : Int32
  /// 其中错配数.
  n_errors : Int32
  /// covariate 表条目数(quality x context).
  n_table_entries : Int32
  /// 重估前平均碱基质量.
  mean_qual_before : Double
  /// 重估后平均碱基质量.
  mean_qual_after : Double
}

///|
/// 读段某参考位置是否落在 known-sites 区间(0-based 半开).
fn bam_pp_in_known_sites(
  known : Array[(Int32, Int32, Int32)],
  ref_id : Int32,
  pos : Int32,
) -> Bool {
  for site in known {
    if site.0 == ref_id && pos >= site.1 && pos < site.2 {
      return true
    }
  }
  false
}

///|
/// 读段第 i 位的三碱基上下文(边界取自身).
fn bam_pp_context(seq : String, i : Int32) -> String {
  let left = if i > 0 {
    seq.unsafe_get(i - 1).unsafe_to_char()
  } else {
    seq.unsafe_get(i).unsafe_to_char()
  }
  let right = if i + 1 < seq.length() {
    seq.unsafe_get(i + 1).unsafe_to_char()
  } else {
    seq.unsafe_get(i).unsafe_to_char()
  }
  left.to_string() +
  seq.unsafe_get(i).unsafe_to_char().to_string() +
  right.to_string()
}

///|
/// GATK BQSR 简化(本地两遍):
/// 1) 训练遍: 逐读段 CIGAR walk, 对 M/=/X 碱基与参考比对, 按
///    (quality, 三碱基上下文) covariate 累计观察/错配; known-sites
///    区间内的差异豁免(视为真实变异);
/// 2) 应用遍: 经验质量 emp = -10 log10((err+1)/(obs+1)) 夹到
///    [min_qual, max_qual], 重写 M/=/X/S 位置的质量字符.
/// 训练仅用双端均映射读段的首读(避免 mate 重复计数).
pub fn bam_bqsr(
  bam : BamFile,
  reference : String,
  known_sites~ : Array[(Int32, Int32, Int32)] = [],
  min_qual~ : Int32 = 2,
  max_qual~ : Int32 = 40,
) -> (BamFile, BqsrReport) {
  // ---- 训练遍 ----
  let obs_table : Map[String, Int32] = Map([], capacity=0)
  let err_table : Map[String, Int32] = Map([], capacity=0)
  let mut n_obs = 0
  let mut n_err = 0
  let mut qual_sum_before = 0.0
  let mut qual_count = 0
  for rec in bam.records {
    if rec.is_unmapped() || rec.ref_id != 0 {
      continue
    }
    // 仅首读参与训练(避免 mate 双计)
    if (rec.flag & 0x1) != 0 && (rec.flag & 0x80) != 0 {
      continue
    }
    let mut ref_pos = rec.pos
    let mut qi = 0
    for ci = 0; ci < rec.cigar.length(); ci = ci + 1 {
      let elem = rec.cigar[ci]
      let op = elem.op
      let len = elem.length
      if op is Match || op is Equal || op is Diff {
        for k = 0; k < len; k = k + 1 {
          let q = if qi < rec.qual.length() {
            rec.qual.unsafe_get(qi).to_int() - 33
          } else {
            30
          }
          qual_sum_before = qual_sum_before + q.to_double()
          qual_count = qual_count + 1
          let known = bam_pp_in_known_sites(known_sites, rec.ref_id, ref_pos)
          if known == false && ref_pos >= 0 && ref_pos < reference.length() {
            let rb = reference.unsafe_get(ref_pos).unsafe_to_char()
            let qb = rec.seq.unsafe_get(qi).unsafe_to_char()
            let ctx = bam_pp_context(rec.seq, qi)
            let key = "\{q}|\{ctx}"
            let o = obs_table.get(key).unwrap_or(0)
            obs_table.set(key, o + 1)
            n_obs = n_obs + 1
            if qb != rb && qb != 'N' && rb != 'N' {
              let e = err_table.get(key).unwrap_or(0)
              err_table.set(key, e + 1)
              n_err = n_err + 1
            }
          }
          ref_pos = ref_pos + 1
          qi = qi + 1
        }
      } else if op is Insert || op is SoftClip {
        for k = 0; k < len; k = k + 1 {
          if op is SoftClip {
            let q = if qi < rec.qual.length() {
              rec.qual.unsafe_get(qi).to_int() - 33
            } else {
              30
            }
            qual_sum_before = qual_sum_before + q.to_double()
            qual_count = qual_count + 1
          }
          qi = qi + 1
        }
      } else if op is Delete || op is Skip {
        ref_pos = ref_pos + len
      }
    }
  }
  // ---- 应用遍 ----
  let ln10 = @math.ln(10.0)
  let out_records : Array[BamRecord] = Array(capacity=bam.records.length())
  let mut qual_sum_after = 0.0
  for rec in bam.records {
    if rec.is_unmapped() || rec.qual.length() == 0 || rec.cigar.length() == 0 {
      out_records.push(rec)
      continue
    }
    let qual_sb = StringBuilder()
    let mut qi = 0
    for ci = 0; ci < rec.cigar.length(); ci = ci + 1 {
      let elem = rec.cigar[ci]
      let op = elem.op
      let len = elem.length
      if op is Match ||
        op is Equal ||
        op is Diff ||
        op is Insert ||
        op is SoftClip {
        for k = 0; k < len; k = k + 1 {
          let q = if qi < rec.qual.length() {
            rec.qual.unsafe_get(qi).to_int() - 33
          } else {
            30
          }
          let ctx = if qi < rec.seq.length() {
            bam_pp_context(rec.seq, qi)
          } else {
            "AAA"
          }
          let key = "\{q}|\{ctx}"
          let new_q : Int32 = match obs_table.get(key) {
            Some(o) => {
              let e = err_table.get(key).unwrap_or(0)
              let rate = (e.to_double() + 1.0) / (o.to_double() + 1.0)
              let emp = -10.0 * @math.ln(rate) / ln10
              let clamped = if emp < min_qual.to_double() {
                min_qual.to_double()
              } else if emp > max_qual.to_double() {
                max_qual.to_double()
              } else {
                emp
              }
              clamped.round().to_int()
            }
            None => q
          }
          let qch = (new_q + 33).to_char().unwrap_or('!')
          qual_sb.write_string(qch.to_string())
          qual_sum_after = qual_sum_after + new_q.to_double()
          qi = qi + 1
        }
      } else {
        // Delete/Skip/HardClip/Pad 无质量位
        ignore(elem)
      }
    }
    out_records.push(bam_pp_with_mate(rec, qual=Some(qual_sb.to_string())))
  }
  let mean_before : Double = if qual_count > 0 {
    qual_sum_before / qual_count.to_double()
  } else {
    0.0
  }
  let mean_after : Double = if qual_count > 0 {
    qual_sum_after / qual_count.to_double()
  } else {
    0.0
  }
  let report = BqsrReport::{
    n_observations: n_obs,
    n_errors: n_err,
    n_table_entries: obs_table.length(),
    mean_qual_before: mean_before,
    mean_qual_after: mean_after,
  }
  (BamFile::{ header: bam.header, records: out_records, }, report)
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// markdup 统计摘要行.
pub fn bam_markdup_stats_line(stats : BamMarkdupStats) -> String {
  "markdup: pair_groups=\{stats.n_pair_groups}, dup_pairs=\{stats.n_dup_pairs}, optical=\{stats.n_optical}, single_groups=\{stats.n_single_groups}, dup_single=\{stats.n_dup_single}, dup_rate=\{stats.dup_rate}"
}

///|
/// BQSR 报告摘要行.
pub fn bam_bqsr_report_line(report : BqsrReport) -> String {
  "bqsr: observations=\{report.n_observations}, errors=\{report.n_errors}, table_entries=\{report.n_table_entries}, mean_qual \{report.mean_qual_before} -> \{report.mean_qual_after}"
}

///|
/// BQSR 报告 -> JSON 值.
pub fn bam_bqsr_value(report : BqsrReport) -> JsonValue {
  JsonObj(
    Map([
      ("n_observations", JsonNum(report.n_observations.to_double())),
      ("n_errors", JsonNum(report.n_errors.to_double())),
      ("n_table_entries", JsonNum(report.n_table_entries.to_double())),
      ("mean_qual_before", JsonNum(report.mean_qual_before)),
      ("mean_qual_after", JsonNum(report.mean_qual_after)),
    ]),
  )
}
