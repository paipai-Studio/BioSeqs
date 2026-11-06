///|
/// phase_lite.mbt — 基于读段的 VCF 定相(read-backed phasing).
///
/// 将 VCF 中的杂合 SNP 与比对读段关联:
/// - 按 CIGAR 把读段碱基投影到参考坐标
/// - 每条跨多个杂合位点的读段给出一组"同相/反相"约束
/// - 用 BFS 在约束图上传播相位，输出 phased GT 与单倍型块
/// - 统计读段与最终相位冲突的比例

///|
/// 定相输入读段(已比对).
pub(all) struct PhaseRead {
  /// 读段名称
  name : String
  /// 参考序列名称
  chrom : String
  /// 比对起点(0-based，与 SAM 一致)
  start : Int32
  /// CIGAR 字符串
  cigar : String
  /// 读段序列
  sequence : String
}

///|
/// 构造定相读段.
pub fn phase_read_new(
  name : String,
  chrom : String,
  start : Int32,
  cigar : String,
  sequence : String,
) -> PhaseRead {
  PhaseRead::{ name, chrom, start, cigar, sequence, }
}

///|
/// 单个位点的定相结果.
pub(all) struct PhaseCall {
  /// 参考序列名称
  chrom : String
  /// VCF 位置(1-based)
  pos : Int32
  /// 参考碱基
  ref_seq : String
  /// 替代碱基
  alt : String
  /// 相位化基因型(如 `1|0`);未定相时为原始 GT
  genotype : String
  /// 单倍型块编号(未定相为 None)
  block : Int32?
}

///|
/// 定相结果.
pub(all) struct PhaseResult {
  /// 单倍型块数量
  blocks : Int32
  /// 读段与最终相位冲突的比例
  conflict_rate : Double
  /// 逐位点结果(保持 VCF 原顺序)
  calls : Array[PhaseCall]
}

///|
/// 定相用的位点内部表示.
priv struct PhaseVariant {
  index : Int32
  chrom : String
  pos : Int32
  ref_seq : String
  alt : String
  base : Char
  alt_char : Char
}

///|
/// 读段碱基 -> 参考坐标(0-based)的投影表.
fn phase_projection(read : PhaseRead) -> Map[Int32, Int32] {
  let projection : Map[Int32, Int32] = Map([], capacity=read.sequence.length())
  let cigar = parse_cigar(read.cigar)
  let mut ref_pos = read.start
  let mut query_pos = 0
  for element in cigar {
    let len = element.length
    match element.op {
      Match | Equal | Diff => {
        for offset = 0; offset < len; offset = offset + 1 {
          if query_pos + offset < read.sequence.length() {
            projection.set(ref_pos + offset, query_pos + offset)
          }
        }
        ref_pos = ref_pos + len
        query_pos = query_pos + len
      }
      Insert | SoftClip => query_pos = query_pos + len
      Delete | Skip => ref_pos = ref_pos + len
      HardClip | Pad => ()
    }
  }
  projection
}

///|
/// 读段在某位点观测到的等位基因(0=ref, 1=alt);无法判断返回 None.
fn phase_observed_allele(
  projection : Map[Int32, Int32],
  sequence : String,
  variant : PhaseVariant,
) -> Int32? {
  if variant.ref_seq.length() != 1 || variant.alt.length() != 1 {
    return None
  }
  match projection.get(variant.pos - 1) {
    Some(query_pos) =>
      if query_pos < 0 || query_pos >= sequence.length() {
        None
      } else {
        let observed = sequence.unsafe_get(query_pos).unsafe_to_char()
        if observed == variant.alt_char {
          Some(1)
        } else if observed == variant.base {
          Some(0)
        } else {
          None
        }
      }
    None => None
  }
}

///|
/// 解析 GT 字符串,返回等位基因数组(非二倍体返回 None).
fn phase_parse_gt(gt : String) -> Array[String]? {
  let normalized = gt.replace_all(old="/", new="|")
  let parts = normalized.split("|").collect()
  if parts.length() != 2 {
    return None
  }
  let alleles : Array[String] = Array(capacity=2)
  for part in parts {
    alleles.push(part.to_owned())
  }
  Some(alleles)
}

///|
/// 基于读段对 VCF 进行定相,返回逐位点结果.
///
/// `sample` 为 VCF 样本列名。只有 0/1 杂合 SNP 会被尝试定相,
/// 缺少读段支持或与其它位点无连接的杂合位点保持原 GT。
pub fn phase_vcf_records(
  vcf : VcfFile,
  reads : Array[PhaseRead],
  sample : String,
) -> PhaseResult {
  let variants : Array[PhaseVariant] = Array(capacity=vcf.records.length())
  for i = 0; i < vcf.records.length(); i = i + 1 {
    let rec = vcf.records[i]
    let mut het = false
    match rec.get_genotype(sample) {
      Some(data) =>
        if data.length() > 0 {
          match phase_parse_gt(data[0]) {
            Some(alleles) =>
              if rec.is_snp() &&
                (
                  (alleles[0] == "0" && alleles[1] == "1") ||
                  (alleles[0] == "1" && alleles[1] == "0")
                ) {
                het = true
              }
            None => ()
          }
        }
      None => ()
    }
    if het {
      variants.push(PhaseVariant::{
        index: i,
        chrom: rec.chrom,
        pos: rec.pos,
        ref_seq: rec.ref_seq,
        alt: rec.alt[0],
        base: rec.ref_seq.unsafe_get(0).unsafe_to_char(),
        alt_char: rec.alt[0].unsafe_get(0).unsafe_to_char(),
      })
    }
  }

  let variant_of_index : Map[Int32, Int32] = Map([], capacity=variants.length())
  for i = 0; i < variants.length(); i = i + 1 {
    variant_of_index.set(variants[i].index, i)
  }
  let n = variants.length()

  // 约束图: phase[a] xor phase[b] = weight.
  let adjacency : Array[Array[(Int32, Int32)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adjacency.push([])
  }
  let mut total_links = 0
  let mut conflicting_links = 0

  for read in reads {
    let projection = phase_projection(read)
    let observed : Array[(Int32, Int32)] = []
    for j = 0; j < n; j = j + 1 {
      let variant = variants[j]
      if variant.chrom == read.chrom {
        match phase_observed_allele(projection, read.sequence, variant) {
          Some(allele) => observed.push((j, allele))
          None => ()
        }
      }
    }
    if observed.length() >= 2 {
      let anchor = observed[0]
      for k = 1; k < observed.length(); k = k + 1 {
        let other = observed[k]
        let weight = anchor.1 ^ other.1
        adjacency[anchor.0].push((other.0, weight))
        adjacency[other.0].push((anchor.0, weight))
        total_links = total_links + 1
      }
    }
  }

  // 在约束图上做 BFS 传播相位.
  let phase : Array[Int32?] = Array::make(n, None)
  for start = 0; start < n; start = start + 1 {
    if phase[start] is Some(_) || adjacency[start].length() == 0 {
      continue
    }
    phase[start] = Some(0)
    let queue : Array[Int32] = [start]
    let mut head = 0
    while head < queue.length() {
      let current = queue[head]
      head = head + 1
      let current_phase = phase[current].unwrap()
      for edge in adjacency[current] {
        let target = edge.0
        let expected = current_phase ^ edge.1
        match phase[target] {
          None => {
            phase[target] = Some(expected)
            queue.push(target)
          }
          Some(existing) =>
            if existing != expected && current < target {
              conflicting_links = conflicting_links + 1
            }
        }
      }
    }
  }

  // 为每个连通分量分配块编号.
  let block_of_variant : Array[Int32?] = Array::make(n, None)
  let visited : Array[Bool] = Array::make(n, false)
  let mut next_block = 0
  for start = 0; start < n; start = start + 1 {
    if phase[start] is None || visited[start] || adjacency[start].length() == 0 {
      continue
    }
    let queue : Array[Int32] = [start]
    visited[start] = true
    let mut head = 0
    while head < queue.length() {
      let current = queue[head]
      head = head + 1
      for edge in adjacency[current] {
        if !visited[edge.0] {
          visited[edge.0] = true
          queue.push(edge.0)
        }
      }
    }
    for id in queue {
      block_of_variant[id] = Some(next_block)
    }
    next_block = next_block + 1
  }

  let calls : Array[PhaseCall] = Array(capacity=vcf.records.length())
  for i = 0; i < vcf.records.length(); i = i + 1 {
    let rec = vcf.records[i]
    let gt = match rec.get_genotype(sample) {
      Some(data) => if data.length() > 0 { data[0] } else { "." }
      None => "."
    }
    let phased = match variant_of_index.get(i) {
      Some(internal) =>
        match (phase[internal], block_of_variant[internal]) {
          (Some(polarity), Some(block)) => {
            let genotype = if polarity == 0 { "1|0" } else { "0|1" }
            (genotype, Some(block))
          }
          _ => (gt, None)
        }
      None => (gt, None)
    }
    calls.push(PhaseCall::{
      chrom: rec.chrom,
      pos: rec.pos,
      ref_seq: rec.ref_seq,
      alt: if rec.alt.length() > 0 {
        rec.alt[0]
      } else {
        "."
      },
      genotype: phased.0,
      block: phased.1,
    })
  }

  let conflict_rate = if total_links > 0 {
    conflicting_links.to_double() / total_links.to_double()
  } else {
    0.0
  }
  PhaseResult::{ blocks: next_block, conflict_rate, calls, }
}

///|
/// 定相结果 -> JSON 文本.
pub fn phase_result_json(result : PhaseResult) -> String {
  let calls : Array[JsonValue] = Array(capacity=result.calls.length())
  for call in result.calls {
    let block = match call.block {
      Some(value) => JsonNum(value.to_double())
      None => JsonNull
    }
    calls.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(call.chrom)),
          ("pos", JsonNum(call.pos.to_double())),
          ("ref", JsonStr(call.ref_seq)),
          ("alt", JsonStr(call.alt)),
          ("genotype", JsonStr(call.genotype)),
          ("block", block),
        ]),
      ),
    )
  }
  json_to_string(
    JsonObj(
      Map([
        ("blocks", JsonNum(result.blocks.to_double())),
        ("conflict_rate", JsonNum(result.conflict_rate)),
        ("calls", JsonArr(calls)),
      ]),
    ),
  )
}

///|
/// 定相结果 -> TSV 文本.
pub fn phase_result_tsv(result : PhaseResult) -> String {
  let sb = StringBuilder()
  sb.write_string("chrom\tpos\tref\talt\tgenotype\tblock\n")
  for call in result.calls {
    sb.write_string(call.chrom)
    sb.write_char('\t')
    sb.write_string(call.pos.to_string())
    sb.write_char('\t')
    sb.write_string(call.ref_seq)
    sb.write_char('\t')
    sb.write_string(call.alt)
    sb.write_char('\t')
    sb.write_string(call.genotype)
    sb.write_char('\t')
    match call.block {
      Some(value) => sb.write_string(value.to_string())
      None => sb.write_string(".")
    }
    sb.write_char('\n')
  }
  sb.to_string()
}
