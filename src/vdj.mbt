///|
/// vdj.mbt — 免疫组库分析(V(D)J 重组指认 / CDR3 提取 / 克隆型聚类 / 多样性 / TCRdist-lite).
///
/// 对齐 ImmunARCH / scirpy / Change-O 的核心分析语义:
/// - germline 基因段指认(V/D/J call): k-mer 种子对角线比对 + identity 阈值过滤
/// - CDR3 提取: IMGT 定义, 从 V 区保守 Cys104(含)到 J 区保守 Phe/Trp118(含)
/// - productive 判定: 框内(长度 %3==0) + 无终止密码子 + C 起始 F/W 结尾
/// - 克隆型聚类: (cdr3_aa, v_family, j_family) 三元组分组, 输出计数/频率
/// - 多样性: Shannon / Simpson / Gini / clonality / D50
/// - TCRdist-lite: CDR3 氨基酸全局比对 + IMGT 化学相似组打分
///
/// 内置 germline 库为教学级精简合成参考(各基因段含唯一区分核心 +
/// IMGT 保守端序), 生产使用建议经 `vdj_build_gene_db` 导入 IMGT
/// 完整 germline 序列. 数据不出本机.
pub suberror VdjError {
  VdjError(String)
}

///|
/// 基因段类别.
pub(all) enum VdjGeneKind {
  VGene
  DGene
  JGene
} derive(Eq, Debug)

///|
pub extend VdjGeneKind with Eq::{equal, not_equal}

///|
pub extend VdjGeneKind with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// 一条 germline 基因段参考序列.
pub(all) struct VdjGene {
  name : String
  locus : String
  kind : VdjGeneKind
  sequence : String
}

///|
/// germline 基因库: V/D/J 段集合与预构建的 k-mer 索引.
pub struct VdjGeneDb {
  v_genes : Array[VdjGene]
  d_genes : Array[VdjGene]
  j_genes : Array[VdjGene]
  k : Int32
  /// 每条基因的 k-mer -> 位置列表 索引, 与 (v+d+j) 拼接顺序对齐
  indices : Array[Map[String, Array[Int32]]]
  /// 与 indices 对齐的目标序列缓存
  targets : Array[String]
}

///|
/// 一段比对结果(读序对 germline 基因段).
pub(all) struct VdjAlignment {
  gene : String
  identity : Double
  read_start : Int32
  read_end : Int32
}

///|
/// 一条完成指认的免疫组库链.
pub(all) struct VdjChain {
  sequence_id : String
  v_call : String
  d_call : String
  j_call : String
  v_identity : Double
  j_identity : Double
  cdr3_nt : String
  cdr3_aa : String
  in_frame : Bool
  stop_codon : Bool
  productive : Bool
}

///|
/// 一个克隆型(CDR3 氨基酸 + V/J 基因家族分组).
pub(all) struct VdjClonotype {
  clonotype_id : Int32
  cdr3_aa : String
  v_family : String
  j_family : String
  count : Int32
  frequency : Double
  member_ids : Array[String]
}

///|
/// 组库多样性指标.
pub(all) struct VdjDiversity {
  total_chains : Int32
  total_clonotypes : Int32
  shannon : Double
  simpson : Double
  gini : Double
  clonality : Double
  d50 : Int32
}

// ---------------------------------------------------------------------------
// 合成 germline 参考序列构造(教学级)
// ---------------------------------------------------------------------------

///|
fn vdj_syn_v(core : String) -> String {
  // FR3 保守前导 + 唯一区分核心 + 以编码 Cys104 的 TGT 收尾
  "GACACAGCGGCCACCTGGAAAAGCTGCTG" + core + "CCCCAGACTGGCCACCTGT"
}

///|
fn vdj_syn_j(core : String) -> String {
  // J-REGION 5' 保守 W-G-G(TGGGG) + 唯一核心 + J 常量区
  "TGGGG" + core + "ACGAGCAGCGCCCTGGGCCAAA"
}

///|
fn vdj_gene(
  name : String,
  locus : String,
  kind : VdjGeneKind,
  seq : String,
) -> VdjGene {
  { name, locus, kind, sequence: seq, }
}

///|
/// 教学级精简人 TRB germline 库(6V/3D/6J, 合成参考序列).
pub fn vdj_default_trb_db() -> VdjGeneDb {
  let v : Array[VdjGene] = [
    vdj_gene("TRBV01", "TRB", VGene, vdj_syn_v("AAGGTCGACCTGGAA")),
    vdj_gene("TRBV02", "TRB", VGene, vdj_syn_v("TTCCAGGCCTGGAT")),
    vdj_gene("TRBV03", "TRB", VGene, vdj_syn_v("GACATTCAGGTCCT")),
    vdj_gene("TRBV04", "TRB", VGene, vdj_syn_v("CTGGGACAGGTGTCG")),
    vdj_gene("TRBV05", "TRB", VGene, vdj_syn_v("AGCTGGGACCTTCAG")),
    vdj_gene("TRBV06", "TRB", VGene, vdj_syn_v("TCAGGACCTGGTGCA")),
  ]
  let d : Array[VdjGene] = [
    vdj_gene("TRBD1", "TRB", DGene, "AGCGGGGCTCGGACCGGGAC"),
    vdj_gene("TRBD2", "TRB", DGene, "GTCCTCGGAGCAGAGGCCA"),
    vdj_gene("TRBD3", "TRB", DGene, "AGCTGACGGCCGAGT"),
  ]
  let j : Array[VdjGene] = [
    vdj_gene("TRBJ01", "TRB", JGene, vdj_syn_j("ACCCAGCTGAGT")),
    vdj_gene("TRBJ02", "TRB", JGene, vdj_syn_j("TCGAGCCCCCGGG")),
    vdj_gene("TRBJ03", "TRB", JGene, vdj_syn_j("GTGACCTCGACC")),
    vdj_gene("TRBJ04", "TRB", JGene, vdj_syn_j("CTGGGCTGAAGC")),
    vdj_gene("TRBJ05", "TRB", JGene, vdj_syn_j("GAGCCCCGACGT")),
    vdj_gene("TRBJ06", "TRB", JGene, vdj_syn_j("TGGCCGAGACCA")),
  ]
  vdj_build_gene_db(v, d, j, k=8)
}

///|
/// 教学级精简人 IGH germline 库(4V/2D/4J, 合成参考序列).
pub fn vdj_default_igh_db() -> VdjGeneDb {
  let v : Array[VdjGene] = [
    vdj_gene("IGHV01", "IGH", VGene, vdj_syn_v("GGTTCGCCTGGGAC")),
    vdj_gene("IGHV02", "IGH", VGene, vdj_syn_v("CAGGCTGGAAGTCC")),
    vdj_gene("IGHV03", "IGH", VGene, vdj_syn_v("TTGACCTGGGCCAG")),
    vdj_gene("IGHV04", "IGH", VGene, vdj_syn_v("ACCGGGTTTCAGG")),
  ]
  let d : Array[VdjGene] = [
    vdj_gene("IGHD1", "IGH", DGene, "GGAAGTAGTCGAGCCGGG"),
    vdj_gene("IGHD2", "IGH", DGene, "AGCTGGGGACCGTCG"),
  ]
  let j : Array[VdjGene] = [
    vdj_gene("IGHJ01", "IGH", JGene, vdj_syn_j("ACGACCTGGGGC")),
    vdj_gene("IGHJ02", "IGH", JGene, vdj_syn_j("GGCCCCAGCTG")),
    vdj_gene("IGHJ03", "IGH", JGene, vdj_syn_j("TTGACCCTGGAC")),
    vdj_gene("IGHJ04", "IGH", JGene, vdj_syn_j("CCTGGGACGAGT")),
  ]
  vdj_build_gene_db(v, d, j, k=8)
}

///|
/// 从外部 germline 序列构建基因库(生产推荐: 导入 IMGT 完整库).
pub fn vdj_build_gene_db(
  v : Array[VdjGene],
  d : Array[VdjGene],
  j : Array[VdjGene],
  k~ : Int32 = 8,
) -> VdjGeneDb {
  let all : Array[VdjGene] = Array(capacity=0)
  for g in v {
    all.push(g)
  }
  for g in d {
    all.push(g)
  }
  for g in j {
    all.push(g)
  }
  let indices : Array[Map[String, Array[Int32]]] = Array(capacity=0)
  let targets : Array[String] = Array(capacity=0)
  for g in all {
    let idx : Map[String, Array[Int32]] = Map([], capacity=0)
    let n = g.sequence.length()
    let kk = if k > 0 && k < 16 { k } else { 8 }
    for p = 0; p + kk <= n; p = p + 1 {
      let key = g.sequence[p:p + kk].to_owned()
      match idx.get(key) {
        Some(arr) => arr.push(p)
        None => {
          let arr : Array[Int32] = Array(capacity=0)
          arr.push(p)
          idx.set(key, arr)
        }
      }
    }
    indices.push(idx)
    targets.push(g.sequence)
  }
  { v_genes: v, d_genes: d, j_genes: j, k, indices, targets, }
}

// ---------------------------------------------------------------------------
// 种子对角线比对(指认核心)
// ---------------------------------------------------------------------------

///|
/// 读序对单条 germline 序列的 k-mer 种子对角线比对.
/// 选锚点最多的对角线并逐碱基统计 identity; 无种子命中返回 None.
fn vdj_seed_align(
  query : String,
  target : String,
  index : Map[String, Array[Int32]],
  k : Int32,
) -> (Double, Int32, Int32)? {
  let n = query.length()
  let m = target.length()
  if n < k || m < k {
    return None
  }
  let anchors : Map[Int32, Int32] = Map([], capacity=0)
  for qi = 0; qi + k <= n; qi = qi + 1 {
    let key = query[qi:qi + k].to_owned()
    match index.get(key) {
      Some(tposs) =>
        for tp in tposs {
          let d = qi - tp
          match anchors.get(d) {
            Some(c) => anchors.set(d, c + 1)
            None => anchors.set(d, 1)
          }
        }
      None => ()
    }
  }
  if anchors.length() == 0 {
    return None
  }
  let mut best_d : Int32 = 0
  let mut best_cnt : Int32 = -1
  for d in anchors.keys() {
    let c = anchors.get(d).unwrap()
    if c > best_cnt {
      best_cnt = c
      best_d = d
    }
  }
  // 对角线区间: query[qs, qe) vs target[qs - best_d, qe - best_d)
  let mut qs : Int32 = 0
  let mut qe : Int32 = n
  if best_d >= 0 {
    qs = best_d
  } else {
    qs = -best_d
  }
  let t_hi = m + best_d
  if qe > t_hi {
    qe = t_hi
  }
  if qe <= qs {
    return None
  }
  let mut matched = 0
  let mut mismatches = 0
  for i = qs; i < qe; i = i + 1 {
    let tp = i - best_d
    if query.unsafe_get(i) == target.unsafe_get(tp) {
      matched = matched + 1
    } else {
      mismatches = mismatches + 1
    }
  }
  let total = matched + mismatches
  if total == 0 {
    return None
  }
  Some((matched.to_double() / total.to_double(), qs, qe))
}

///|
/// 基因家族: 等位基因标识前的主名(TRBV05*01 -> TRBV05).
pub fn vdj_gene_family(gene_call : String) -> String {
  let s = gene_call.trim()
  if s.length() == 0 {
    return ""
  }
  let mut cut = s.length()
  let chars : Array[Char] = s.to_array()
  for i = 0; i < chars.length(); i = i + 1 {
    if chars[i] == '*' {
      cut = i
      break
    }
  }
  s[0:cut].to_owned()
}

///|
/// 在 V/D/J 各段集合中找最佳命中.
fn vdj_best_in(
  read : String,
  db : VdjGeneDb,
  genes : Array[VdjGene],
  base : Int32,
) -> VdjAlignment? {
  let mut best : VdjAlignment? = None
  let mut best_identity = 0.0
  for gi = 0; gi < genes.length(); gi = gi + 1 {
    let target = db.targets[base + gi]
    let index = db.indices[base + gi]
    match vdj_seed_align(read, target, index, db.k) {
      Some((identity, rs, re)) => {
        if identity > best_identity {
          best_identity = identity
          best = Some({
            gene: genes[gi].name,
            identity,
            read_start: rs,
            read_end: re,
          })
        }
        ignore(rs)
        ignore(re)
      }
      None => ()
    }
  }
  if best_identity < 0.55 {
    return None
  }
  best
}

// ---------------------------------------------------------------------------
// CDR3 提取(IMGT 定义)
// ---------------------------------------------------------------------------

///|
/// 在读序中寻找编码 Cys104 的保守三联体(TGT/TGC)位置.
/// 返回(位置, 之后首个 J 保守 motif 起点)的候选.
fn vdj_find_cys_positions(read : String) -> Array[Int32] {
  let n = read.length()
  let pos : Array[Int32] = Array(capacity=0)
  for i = 0; i + 3 <= n; i = i + 1 {
    let codon = read[i:i + 3].to_owned()
    if codon == "TGT" || codon == "TGC" {
      pos.push(i)
    }
  }
  pos
}

///|
/// 在位置 from 之后寻找 J-REGION 5' 保守 motif(W-G-G / F-G-G).
/// 返回 motif 中 F/W 密码子的起点(即 CDR3 终点+1 的编码位).
fn vdj_find_j_motif(read : String, from : Int32) -> Int32? {
  // J-REGION 5' 保守 6nt motif: (TGG|TTC) + GGx (W/F + Gly)
  let n = read.length()
  for i = from; i + 6 <= n; i = i + 1 {
    let fw = read[i:i + 3].to_owned()
    if fw != "TGG" && fw != "TTC" {
      continue
    }
    let g1 = read[i + 3:i + 6].to_owned()
    if g1 != "GGT" && g1 != "GGC" && g1 != "GGA" && g1 != "GGG" {
      continue
    }
    return Some(i)
  }
  None
}

///|
/// 提取 CDR3 核苷酸(IMGT: Cys104 含 -> Phe/Trp118 含).
/// 从读序中定位最后一个有效 Cys(其后存在 J motif)并返回 [cys, j+3) 片段.
pub fn vdj_extract_cdr3_nt(read : String) -> String? {
  let cys_positions = vdj_find_cys_positions(read)
  let mut i = cys_positions.length() - 1
  while i >= 0 {
    let c = cys_positions[i]
    match vdj_find_j_motif(read, c + 3) {
      Some(j) => {
        let seg = read[c:j + 3].to_owned()
        return Some(seg)
      }
      None => i = i - 1
    }
  }
  None
}

///|
/// CDR3 核苷酸翻译并判定 productive.
/// 返回 (cdr3_aa, in_frame, stop_codon, productive).
pub fn vdj_cdr3_aa(cdr3_nt : String) -> (String, Bool, Bool, Bool) {
  let aa = biostrings_translate(cdr3_nt)
  let n = cdr3_nt.length()
  let in_frame = n % 3 == 0 && n >= 6
  let stop_codon = aa.contains("*")
  let productive = if in_frame && !stop_codon && aa.length() >= 2 {
    let first = aa.unsafe_get(0).unsafe_to_char()
    let last = aa.unsafe_get(aa.length() - 1).unsafe_to_char()
    first == 'C' && (last == 'F' || last == 'W')
  } else {
    false
  }
  (aa, in_frame, stop_codon, productive)
}

// ---------------------------------------------------------------------------
// V(D)J 指认
// ---------------------------------------------------------------------------

///|
/// 对一条读序做 V/D/J 基因段指认并提取 CDR3.
/// V 命中时优先在 V 比对终点附近定位 Cys; 否则全读序扫描.
pub fn vdj_assign(
  read : String,
  db : VdjGeneDb,
  sequence_id~ : String = "",
) -> VdjChain {
  let n_v = db.v_genes.length()
  let n_d = db.d_genes.length()
  let v_align = vdj_best_in(read, db, db.v_genes, 0)
  let d_align = vdj_best_in(read, db, db.d_genes, n_v)
  let j_align = vdj_best_in(read, db, db.j_genes, n_v + n_d)
  let v_call = match v_align {
    Some(a) => a.gene
    None => ""
  }
  let d_call = match d_align {
    Some(a) => a.gene
    None => ""
  }
  let j_call = match j_align {
    Some(a) => a.gene
    None => ""
  }
  let v_identity = match v_align {
    Some(a) => a.identity
    None => 0.0
  }
  let j_identity = match j_align {
    Some(a) => a.identity
    None => 0.0
  }
  // CDR3: 优先取 V 终点附近的 Cys(V 段以 TGT 收尾时即 read 终点-3)
  let cdr3_nt : String? = match v_align {
    Some(a) => {
      let near = a.read_end - 3
      let cys_positions = vdj_find_cys_positions(read)
      let mut chosen : Int32? = None
      let mut best_dist : Int32 = 0x7FFFFFFF
      for p in cys_positions {
        let dist = if p > near { p - near } else { near - p }
        if dist < best_dist {
          best_dist = dist
          chosen = Some(p)
        }
      }
      match chosen {
        Some(c) =>
          match vdj_find_j_motif(read, c + 3) {
            Some(j) => Some(read[c:j + 3].to_owned())
            None => vdj_extract_cdr3_nt(read)
          }
        None => vdj_extract_cdr3_nt(read)
      }
    }
    None => vdj_extract_cdr3_nt(read)
  }
  let (cdr3_aa, in_frame, stop_codon, productive) = match cdr3_nt {
    Some(nt) => vdj_cdr3_aa(nt)
    None => ("", false, false, false)
  }
  {
    sequence_id,
    v_call,
    d_call,
    j_call,
    v_identity,
    j_identity,
    cdr3_nt: match cdr3_nt {
      Some(nt) => nt
      None => ""
    },
    cdr3_aa,
    in_frame,
    stop_codon,
    productive,
  }
}

///|
/// 批量指认.
pub fn vdj_assign_all(
  reads : Array[(String, String)],
  db : VdjGeneDb,
) -> Array[VdjChain] {
  let out : Array[VdjChain] = Array(capacity=0)
  for pair in reads {
    out.push(vdj_assign(pair.0, db, sequence_id=pair.1))
  }
  out
}

// ---------------------------------------------------------------------------
// 克隆型聚类
// ---------------------------------------------------------------------------

///|
/// 克隆型聚类: 按 (cdr3_aa, v_family, j_family) 分组.
/// 仅 productive 链参与; 输出按计数降序, clonotype_id 从 1 起.
pub fn vdj_clonotypes(chains : Array[VdjChain]) -> Array[VdjClonotype] {
  let total = chains.length()
  if total == 0 {
    return []
  }
  let groups : Map[String, Array[String]] = Map([], capacity=0)
  let meta : Map[String, (String, String)] = Map([], capacity=0)
  for ch in chains {
    if !ch.productive || ch.cdr3_aa.length() == 0 {
      continue
    }
    let vf = vdj_gene_family(ch.v_call)
    let jf = vdj_gene_family(ch.j_call)
    let key = ch.cdr3_aa + "|" + vf + "|" + jf
    match groups.get(key) {
      Some(ids) => ids.push(ch.sequence_id)
      None => {
        let ids : Array[String] = Array(capacity=0)
        ids.push(ch.sequence_id)
        groups.set(key, ids)
        meta.set(key, (vf, jf))
      }
    }
  }
  let entries : Array[(String, Array[String], String, String, String)] = Array(
    capacity=0,
  )
  for key in groups.keys() {
    let ids = groups.get(key).unwrap()
    let (vf, jf) = meta.get(key).unwrap()
    let cdr3 = key[0:key.length() - vf.length() - jf.length() - 2].to_owned()
    entries.push((key, ids, vf, jf, cdr3))
  }
  entries.sort_by(fn(a, b) {
    if a.1.length() != b.1.length() {
      b.1.length() - a.1.length()
    } else {
      compare(a.0, b.0)
    }
  })
  let out : Array[VdjClonotype] = Array(capacity=0)
  let mut cid = 0
  for entry in entries {
    cid = cid + 1
    out.push({
      clonotype_id: cid,
      cdr3_aa: entry.4,
      v_family: entry.2,
      j_family: entry.3,
      count: entry.1.length(),
      frequency: entry.1.length().to_double() / total.to_double(),
      member_ids: entry.1,
    })
  }
  out
}

// ---------------------------------------------------------------------------
// 多样性指标
// ---------------------------------------------------------------------------

///|
/// 计算组库多样性指标(基于克隆型计数分布).
pub fn vdj_diversity(chains : Array[VdjChain]) -> VdjDiversity {
  let clonotypes = vdj_clonotypes(chains)
  let total = chains.length()
  let k = clonotypes.length()
  if k == 0 {
    return {
      total_chains: total,
      total_clonotypes: 0,
      shannon: 0.0,
      simpson: 0.0,
      gini: 0.0,
      clonality: 0.0,
      d50: 0,
    }
  }
  let counts : Array[Double] = Array(capacity=0)
  for c in clonotypes {
    counts.push(c.count.to_double())
  }
  // Shannon 熵
  let mut h = 0.0
  for c in counts {
    let p = c / total.to_double()
    if p > 0.0 {
      h = h - p * @math.ln(p)
    }
  }
  // clonality = 1 - H / ln(k)
  let hmax = @math.ln(k.to_double())
  // k=1(单克隆型)时 ln(1)=0, clonality 惯例取 1(完全扩张)
  let clonality = if hmax > 0.0 {
    1.0 - h / hmax
  } else if k == 1 {
    1.0
  } else {
    0.0
  }
  // Simpson (1 - sum p^2)
  let mut s2 = 0.0
  for c in counts {
    let p = c / total.to_double()
    s2 = s2 + p * p
  }
  let simpson = 1.0 - s2
  // Gini(洛伦兹曲线)
  let sorted_counts : Array[Double] = Array(capacity=0)
  for c in counts {
    sorted_counts.push(c)
  }
  sorted_counts.sort()
  let cum : Array[Double] = Array(capacity=0)
  let mut acc = 0.0
  for c in sorted_counts {
    acc = acc + c
    cum.push(acc)
  }
  let grand = acc
  let mut area = 0.0
  let kk = sorted_counts.length().to_double()
  for i = 0; i < sorted_counts.length(); i = i + 1 {
    let y0 = if i == 0 { 0.0 } else { cum[i - 1] / grand }
    let y1 = cum[i] / grand
    let x0 = i.to_double() / kk
    let x1 = (i + 1).to_double() / kk
    area = area + (y0 + y1) / 2.0 * (x1 - x0)
  }
  let gini = 1.0 - 2.0 * area
  // D50: 累计频率达 50% 所需的最小克隆型数
  let mut d50 : Int32 = 0
  let mut freq_acc = 0.0
  for c in clonotypes {
    if freq_acc >= 0.5 {
      break
    }
    freq_acc = freq_acc + c.frequency
    d50 = d50 + 1
  }
  {
    total_chains: total,
    total_clonotypes: k,
    shannon: h,
    simpson,
    gini,
    clonality,
    d50,
  }
}

// ---------------------------------------------------------------------------
// TCRdist-lite(CDR3 氨基酸距离)
// ---------------------------------------------------------------------------

///|
/// IMGT 化学相似组(同组残基相互替换得高分).
fn vdj_aa_group(c : Char) -> Int32 {
  match c {
    'F' | 'W' | 'Y' => 1
    'I' | 'L' | 'M' | 'V' | 'A' => 2
    'G' | 'S' | 'T' | 'P' => 3
    'D' | 'E' => 4
    'N' | 'Q' => 5
    'K' | 'R' | 'H' => 6
    'C' => 7
    _ => 0
  }
}

///|
/// TCRdist-lite: 两条 CDR3 氨基酸序列的全局比对距离.
/// 打分: 同残基 +5, 同化学组 +2, 异组 -4, gap -6;
/// 距离 = (5 * max(la, lb) - score), 归一化到 [0, +inf), 相同序列为 0.
pub fn vdj_cdr3_distance(a : String, b : String) -> Double {
  let la = a.length()
  let lb = b.length()
  if la == 0 && lb == 0 {
    return 0.0
  }
  if la == 0 || lb == 0 {
    return 5.0 * (if la > lb { la } else { lb }).to_double()
  }
  let n = la
  let m = lb
  // DP: dp[i][j] = 最佳比对得分
  let dp : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i <= n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j <= m; j = j + 1 {
      row.push(0.0)
    }
    dp.push(row)
  }
  for i = 1; i <= n; i = i + 1 {
    dp[i][0] = dp[i - 1][0] - 6.0
  }
  for j = 1; j <= m; j = j + 1 {
    dp[0][j] = dp[0][j - 1] - 6.0
  }
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let ca = a.unsafe_get(i - 1).unsafe_to_char()
      let cb = b.unsafe_get(j - 1).unsafe_to_char()
      let sub = if ca == cb {
        5.0
      } else if vdj_aa_group(ca) == vdj_aa_group(cb) && vdj_aa_group(ca) != 0 {
        2.0
      } else {
        -4.0
      }
      let diag = dp[i - 1][j - 1] + sub
      let up = dp[i - 1][j] - 6.0
      let left = dp[i][j - 1] - 6.0
      let best = if diag >= up && diag >= left {
        diag
      } else if up >= left {
        up
      } else {
        left
      }
      dp[i][j] = best
    }
  }
  let max_len = if la > lb { la } else { lb }
  let max_score = 5.0 * max_len.to_double()
  let d = max_score - dp[n][m]
  if d < 0.0 {
    0.0
  } else {
    d
  }
}

///|
/// TCRdist-lite 矩阵: 克隆型(或链)间两两 CDR3 距离.
pub fn vdj_tcrdist_matrix(aas : Array[String]) -> Array[Array[Double]] {
  let n = aas.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else if j < i {
        row.push(out[j][i])
      } else {
        row.push(vdj_cdr3_distance(aas[i], aas[j]))
      }
    }
    out.push(row)
  }
  out
}
