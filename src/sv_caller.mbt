///|
/// SV 长读长检测(minimap 比对水平).
///
/// 定位: 基于 minimap 比对的 SV 检出, 支持 PacBio/ONT 长读长数据.
/// 复用 minimap.mbt 的 minimizer 种子查找 + 链化 DP 能力,
/// 从 CIGAR 操作提取 SV 候选, 按支持读数过滤.
///
/// 依赖: minimap.mbt, sam.mbt
/// 数据不出本机.

// ---------------------------------------------------------------------------
// SV 调用结果类型
// ---------------------------------------------------------------------------

///|
/// 长读长 SV 调用结果.
pub struct SvCall {
  chrom : String
  sv_type : SvType
  start : Int32
  end : Int32
  support_reads : Int32
  total_reads : Int32
  score : Double
}

///|
/// 染色体.
pub fn SvCall::chrom(self : SvCall) -> String {
  self.chrom
}

///|
/// SV 类型.
pub fn SvCall::sv_type(self : SvCall) -> SvType {
  self.sv_type
}

///|
/// 起始位置(0-based 半开).
pub fn SvCall::start(self : SvCall) -> Int32 {
  self.start
}

///|
/// 终止位置(0-based 半开).
pub fn SvCall::end(self : SvCall) -> Int32 {
  self.end
}

///|
/// 支持读数.
pub fn SvCall::support_reads(self : SvCall) -> Int32 {
  self.support_reads
}

///|
/// 总覆盖读数.
pub fn SvCall::total_reads(self : SvCall) -> Int32 {
  self.total_reads
}

///|
/// 置信度得分 [0,1].
pub fn SvCall::score(self : SvCall) -> Double {
  self.score
}

///|
/// 高置信度判定: 得分 >= 0.8.
pub fn SvCall::is_high_confidence(self : SvCall) -> Bool {
  self.score >= 0.8
}

// ---------------------------------------------------------------------------
// SV 检测参数
// ---------------------------------------------------------------------------

///|
/// SV 检测参数.
pub struct SvParams {
  min_support : Int32
  min_sv_len : Int32
  max_ins_len : Int32
  max_gap : Int32
}

///|
/// 默认参数.
pub fn SvParams::default() -> SvParams {
  SvParams::{
    min_support: 3,
    min_sv_len: 50,
    max_ins_len: 1000,
    max_gap: 1000,
  }
}

///|
/// 自定义参数.
pub fn SvParams::new(
  min_support~ : Int32 = 3,
  min_sv_len~ : Int32 = 50,
  max_ins_len~ : Int32 = 1000,
  max_gap~ : Int32 = 1000,
) -> SvParams {
  SvParams::{ min_support, min_sv_len, max_ins_len, max_gap, }
}

// ---------------------------------------------------------------------------
// SV 检测核心(基于 CIGAR 路径)
// ---------------------------------------------------------------------------

///|
/// 从单条 SAM 记录的 CIGAR 提取 SV 候选.
/// 返回 SV 调用列表.
pub fn sv_candidates_from_sam(
  chrom : String,
  pos : Int32,
  cigar : Array[CigarElement],
  params : SvParams,
) -> Array[SvCall] {
  let calls : Array[SvCall] = Array(capacity=0)
  let mut ref_pos = pos
  let mut i = 0
  while i < cigar.length() {
    let elem = cigar[i]
    let op = elem.op
    let len = elem.length
    match op {
      CigarOp::Match => ref_pos = ref_pos + len
      CigarOp::Delete => {
        // 缺失: 参考序列中跳过的碱基
        if len >= params.min_sv_len {
          calls.push(SvCall::{
            chrom,
            sv_type: SvType::DEL,
            start: ref_pos,
            end: ref_pos + len,
            support_reads: 1,
            total_reads: 1,
            score: 1.0,
          })
        }
        ref_pos = ref_pos + len
      }
      CigarOp::Insert =>
        // 插入: 读段中额外碱基
        if len >= params.min_sv_len && len <= params.max_ins_len {
          calls.push(SvCall::{
            chrom,
            sv_type: SvType::INS,
            start: ref_pos,
            end: ref_pos,
            support_reads: 1,
            total_reads: 1,
            score: 1.0,
          })
        }
      CigarOp::Skip => {
        // 跳过: 内含子/大片段缺失
        if len >= params.min_sv_len {
          calls.push(SvCall::{
            chrom,
            sv_type: SvType::DEL,
            start: ref_pos,
            end: ref_pos + len,
            support_reads: 1,
            total_reads: 1,
            score: 1.0,
          })
        }
        ref_pos = ref_pos + len
      }
      CigarOp::SoftClip => ()
      CigarOp::HardClip => ()
      CigarOp::Pad => ()
      CigarOp::Equal => ref_pos = ref_pos + len
      CigarOp::Diff => ref_pos = ref_pos + len
    }
    i = i + 1
  }
  calls
}

///|
/// 合并相同位置的 SV 候选(按支持读数聚合).
///
/// depth_by_chrom 提供各染色体总覆盖读数, 用于计算置信度
/// score = support_reads / total_reads. 未提供时总读数暂按支持读数计.
pub fn merge_sv_candidates(
  calls : Array[SvCall],
  depth_by_chrom~ : Map[String, Int32] = Map([]),
) -> Array[SvCall] {
  if calls.length() == 0 {
    return calls
  }
  // 按 (chrom, start, end, sv_type) 分组
  let groups : Map[String, Array[SvCall]] = Map([], capacity=0)
  for c in calls {
    let key = c.chrom +
      ":" +
      c.start.to_string() +
      "-" +
      c.end.to_string() +
      "-" +
      sv_type_label(c.sv_type)
    match groups.get(key) {
      Some(arr) => arr.push(c)
      None => {
        let arr : Array[SvCall] = Array(capacity=0)
        arr.push(c)
        groups[key] = arr
      }
    }
  }
  let merged : Array[SvCall] = Array(capacity=0)
  for key in groups.keys() {
    let arr = groups[key]
    if arr.length() == 0 {
      continue
    }
    let first = arr[0]
    let support = arr.length()
    let depth = match depth_by_chrom.get(first.chrom) {
      Some(d) => d
      None => 0
    }
    let total = if depth > support { depth } else { support }
    merged.push(SvCall::{
      chrom: first.chrom,
      sv_type: first.sv_type,
      start: first.start,
      end: first.end,
      support_reads: support,
      total_reads: total,
      score: support.to_double() / total.to_double(),
    })
  }
  merged
}

///|
/// 过滤低支持度 SV.
pub fn filter_sv_calls(
  calls : Array[SvCall],
  min_support : Int32,
) -> Array[SvCall] {
  let out : Array[SvCall] = Array(capacity=0)
  for c in calls {
    if c.support_reads >= min_support {
      out.push(c)
    }
  }
  out
}

///|
/// 按得分降序排序.
pub fn sort_sv_calls_by_score(calls : Array[SvCall]) -> Array[SvCall] {
  let result = calls
  let mut i = 0
  while i < result.length() {
    let mut j = i + 1
    while j < result.length() {
      let vi = result[i]
      let vj = result[j]
      if vj.score > vi.score {
        result[i] = vj
        result[j] = vi
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ---------------------------------------------------------------------------
// 批量 SV 检测入口
// ---------------------------------------------------------------------------

///|
/// 批量 SAM 记录 SV 检测.
/// 输入: (chrom, pos, cigar) 列表.
pub fn sv_caller_batch(
  records : Array[(String, Int32, Array[CigarElement])],
  params : SvParams,
  depth_by_chrom~ : Map[String, Int32] = Map([]),
) -> Array[SvCall] {
  let all_candidates : Array[SvCall] = Array(capacity=0)
  for record in records {
    let (chrom, pos, cigar) = record
    let candidates = sv_candidates_from_sam(chrom, pos, cigar, params)
    for c in candidates {
      all_candidates.push(c)
    }
  }
  let merged = merge_sv_candidates(all_candidates, depth_by_chrom~)
  let filtered = filter_sv_calls(merged, params.min_support)
  sort_sv_calls_by_score(filtered)
}

// ---------------------------------------------------------------------------
// SA tag(supplementary alignment)解析与 DUP/INV/BND 检出
// ---------------------------------------------------------------------------

///|
/// 补充比对片段(来自 SAM SA tag).
pub struct SaSegment {
  chrom : String
  pos : Int32
  strand : Int32
}

///|
/// 染色体.
pub fn SaSegment::chrom(self : SaSegment) -> String {
  self.chrom
}

///|
/// 1-based 起始位置.
pub fn SaSegment::pos(self : SaSegment) -> Int32 {
  self.pos
}

///|
/// 链: 0 正链, 1 反链.
pub fn SaSegment::strand(self : SaSegment) -> Int32 {
  self.strand
}

///|
/// 解析 SAM SA tag 为补充比对片段列表.
/// 格式: rname,pos,strand,CIGAR,mapQ,NM;...
pub fn parse_sa_tag(sa : String) -> Array[SaSegment] {
  let out : Array[SaSegment] = Array(capacity=0)
  if sa.length() == 0 {
    return out
  }
  for entry in sa.split(";") {
    if entry.length() == 0 {
      continue
    }
    let fields = entry.split(",").to_array()
    if fields.length() < 3 {
      continue
    }
    let chrom = fields[0].to_owned()
    let pos = parse_int(fields[1].to_owned())
    let strand = if fields[2].to_owned() == "-" { 1 } else { 0 }
    out.push(SaSegment::{ chrom, pos, strand, })
  }
  out
}

///|
/// 计算 CIGAR 在参考序列上的跨度(消耗参考的操作长度之和).
pub fn cigar_ref_span(cigar : Array[CigarElement]) -> Int32 {
  let mut span = 0
  for e in cigar {
    match e.op {
      CigarOp::Match
      | CigarOp::Delete
      | CigarOp::Skip
      | CigarOp::Equal
      | CigarOp::Diff => span = span + e.length
      CigarOp::Insert | CigarOp::SoftClip | CigarOp::HardClip | CigarOp::Pad =>
        ()
    }
  }
  span
}

///|
/// 依据 SA tag 检出 DUP/INV/BND.
/// 判定规则:
/// - 补充片段位于不同染色体 -> BND
/// - 同染色体且链相反 -> INV
/// - 同染色体同链且补充片段在参考上更靠左(读段倒挂) -> DUP
/// 同染色体同链且补充片段在右侧视为普通分段比对, 不报 SV.
pub fn sv_candidates_from_sa(
  chrom : String,
  pos : Int32,
  strand : Int32,
  sa_tag : String,
  params : SvParams,
) -> Array[SvCall] {
  let out : Array[SvCall] = Array(capacity=0)
  for seg in parse_sa_tag(sa_tag) {
    if seg.chrom != chrom {
      out.push(SvCall::{
        chrom,
        sv_type: SvType::BND,
        start: pos,
        end: pos,
        support_reads: 1,
        total_reads: 1,
        score: 1.0,
      })
      continue
    }
    let low = if seg.pos < pos { seg.pos } else { pos }
    let high = if seg.pos < pos { pos } else { seg.pos }
    if seg.strand != strand {
      if high - low >= params.min_sv_len {
        out.push(SvCall::{
          chrom,
          sv_type: SvType::INV,
          start: low,
          end: high,
          support_reads: 1,
          total_reads: 1,
          score: 1.0,
        })
      }
      continue
    }
    if seg.pos < pos && pos - seg.pos >= params.min_sv_len {
      out.push(SvCall::{
        chrom,
        sv_type: SvType::DUP,
        start: seg.pos,
        end: pos,
        support_reads: 1,
        total_reads: 1,
        score: 1.0,
      })
    }
  }
  out
}
