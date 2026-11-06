///|
/// vg_lite.mbt — vg 风格泛基因组变异图基因分型
/// (本地确定性实现).
///
/// 数据流:
/// - VCF 变异集 + 参考序列 -> 变异图(节点 = 序列块, 边 = 邻接):
///   变异位点切分参考为块, 各等位(含 ref)成替代路径节点;
/// - 图路径枚举: source(入度 0)到 sink(出度 0)的全部路径,
///   每条路径即一条单倍型序列;
/// - 图上基因分型: 读段对各路径序列做滑窗最小 Hamming,
///   分配给最佳路径(确定性 tie-break 取路径索引小者),
///   路径读段计数前二构成二倍体 GT(ref 路径索引为 0).
///
/// 已知简化: 只处理双等位 SNV/indel 的无重叠变异;
/// 图比对用滑窗 Hamming(vg 为全动态规划图比对).
/// 数据不出本机.

///|
/// 图节点(序列块).
pub(all) struct VgNode {
  id : Int32
  seq : String
  /// "ref" 或 "alt".
  kind : String
}

///|
/// 有向边.
pub(all) struct VgEdge {
  from_id : Int32
  to_id : Int32
}

///|
/// 变异图.
pub(all) struct VgGraph {
  chrom : String
  nodes : Array[VgNode]
  edges : Array[VgEdge]
}

///|
/// 一条图路径(节点 id 序列).
pub(all) struct VgPath {
  node_ids : Array[Int32]
  /// 拼接后的单倍型序列.
  seq : String
}

///|
/// 基因分型结果.
pub(all) struct VgGenotype {
  /// GT 字符串(ref 路径为 "0", 其余按路径序).
  gt : String
  /// 两条等位路径的读段支持计数.
  counts : Array[Int32]
  /// 全部路径读段计数(与 gt 对应前二为等位).
  path_counts : Array[Int32]
  n_reads : Int32
}

// ---------------------------------------------------------------------------
// 图构造
// ---------------------------------------------------------------------------

///|
/// 从参考 + VCF 记录构造变异图(双等位, 变异位点按坐标升序且不重叠).
/// 输入未排序时内部排序; 重叠变异跳过并保持图一致性.
pub fn vg_build_graph(
  chrom : String,
  reference : String,
  records : Array[VcfRecord],
) -> VgGraph {
  // 过滤同染色体记录并按 pos 排序
  let vars : Array[VcfRecord] = Array(capacity=0)
  for rec in records {
    if rec.chrom == chrom && rec.alt.length() > 0 {
      vars.push(rec)
    }
  }
  vars.sort_by(fn(a, b) { a.pos - b.pos })
  // 去重叠(保留先出现的)
  let kept : Array[VcfRecord] = Array(capacity=0)
  let mut last_end : Int32 = -1
  for rec in vars {
    let start = rec.pos // 1-based
    let end_ = rec.pos + rec.ref_seq.length() - 1
    if start > last_end {
      kept.push(rec)
      last_end = end_
    }
  }
  let nodes : Array[VgNode] = Array(capacity=0)
  let edges : Array[VgEdge] = Array(capacity=0)
  let mut next_id : Int32 = 0
  let add_node = fn(seq : String, kind : String, n : Array[VgNode]) -> Int32 {
    let id = next_id
    next_id = next_id + 1
    n.push(VgNode::{ id, seq, kind, })
    id
  }
  let add_edge = fn(from_id : Int32, to_id : Int32, e : Array[VgEdge]) -> Unit {
    e.push(VgEdge::{ from_id, to_id, })
  }
  // 参考块: [prev_end+1, pos-1] (1-based 闭区间 -> 0-based 半开)
  let mut prev_ref_end : Int32 = 0 // 0-based 半开, 已消耗参考长度
  // 待连接到下一块的节点集合(上一变异的全部等位末端或上一参考块)
  let mut pending : Array[Int32] = []
  for rec in kept {
    let pos0 = rec.pos - 1
    // 参考前缀块 [prev_ref_end, pos0)
    if pos0 > prev_ref_end {
      let prefix = reference[prev_ref_end:pos0].to_owned()
      let bid = add_node(prefix, "ref", nodes)
      for p in pending {
        add_edge(p, bid, edges)
      }
      pending = [bid]
    }
    // ref 等位块 + alt 等位块(全部从 pending 连入)
    let allele_ends : Array[Int32] = Array(capacity=0)
    let ref_id = add_node(rec.ref_seq, "ref", nodes)
    for p in pending {
      add_edge(p, ref_id, edges)
    }
    allele_ends.push(ref_id)
    for a in rec.alt {
      if a == rec.ref_seq || a == "." {
        continue
      }
      let alt_id = add_node(a, "alt", nodes)
      for p in pending {
        add_edge(p, alt_id, edges)
      }
      allele_ends.push(alt_id)
    }
    pending = allele_ends
    prev_ref_end = pos0 + rec.ref_seq.length()
  }
  // 尾参考块(从全部等位末端连入)
  if prev_ref_end < reference.length() {
    let suffix = reference[prev_ref_end:].to_owned()
    let bid = add_node(suffix, "ref", nodes)
    for p in pending {
      add_edge(p, bid, edges)
    }
    pending = [bid]
  }
  ignore(pending)
  VgGraph::{ chrom, nodes, edges, }
}

// ---------------------------------------------------------------------------
// 路径枚举
// ---------------------------------------------------------------------------

///|
fn vg_out_edges(graph : VgGraph) -> Map[Int32, Array[Int32]] {
  let adj : Map[Int32, Array[Int32]] = Map([], capacity=0)
  for e in graph.edges {
    match adj.get(e.from_id) {
      Some(arr) => arr.push(e.to_id)
      None => adj.set(e.from_id, [e.to_id])
    }
  }
  adj
}

///|
fn vg_node_seq(graph : VgGraph, id : Int32) -> String {
  for n in graph.nodes {
    if n.id == id {
      return n.seq
    }
  }
  ""
}

///|
/// 顶层递归 DFS 枚举 source->sink 全部路径.
/// (source = 入度 0 节点; sink = 出边为空或指向已访问节点的叶子)
fn vg_dfs(
  graph : VgGraph,
  adj : Map[Int32, Array[Int32]],
  node : Int32,
  visited : Array[Bool],
  path : Array[Int32],
  out : Array[VgPath],
) -> Unit {
  visited[node] = true
  path.push(node)
  let nexts = adj.get(node)
  let is_sink = match nexts {
    Some(arr) =>
      arr.length() == 0 ||
      ({
        let mut all_visited = true
        for nxt in arr {
          if visited[nxt] == false {
            all_visited = false
            break
          }
        }
        all_visited
      })
    None => true
  }
  if is_sink {
    let seq_sb = StringBuilder()
    for id in path {
      seq_sb.write_string(vg_node_seq(graph, id))
    }
    out.push(VgPath::{ node_ids: path.copy(), seq: seq_sb.to_string(), })
  } else {
    match nexts {
      Some(arr) =>
        for nxt in arr {
          if visited[nxt] == false {
            vg_dfs(graph, adj, nxt, visited, path, out)
          }
        }
      None => ()
    }
  }
  ignore(path.pop())
  visited[node] = false
}

///|
/// 枚举全部 source -> sink 路径(每条为一条单倍型序列).
/// 变异图无环(参考坐标单调), 路径数 = 各位点等位数之积.
pub fn vg_paths(graph : VgGraph) -> Array[VgPath] {
  let adj = vg_out_edges(graph)
  // 入度表
  let in_deg : Map[Int32, Int32] = Map([], capacity=0)
  for e in graph.edges {
    let d = in_deg.get(e.to_id).unwrap_or(0)
    in_deg.set(e.to_id, d + 1)
  }
  let max_id : Int32 = {
    let mut m = -1
    for n in graph.nodes {
      if n.id > m {
        m = n.id
      }
    }
    m
  }
  if max_id < 0 {
    return []
  }
  let visited : Array[Bool] = Array::make(max_id + 1, false)
  let out : Array[VgPath] = Array(capacity=0)
  let path : Array[Int32] = Array(capacity=0)
  for n in graph.nodes {
    if in_deg.get(n.id).unwrap_or(0) == 0 {
      vg_dfs(graph, adj, n.id, visited, path, out)
    }
  }
  out
}

// ---------------------------------------------------------------------------
// 图上基因分型
// ---------------------------------------------------------------------------

///|
/// 读段对序列的滑窗最小 Hamming(窗口 = 读段长;
/// 序列短于读段时返回读段长即完全失配).
fn vg_min_hamming(read : String, seq : String) -> Int32 {
  let n = read.length()
  let m = seq.length()
  if n == 0 || m < n {
    return n
  }
  let mut best = n
  let mut start = 0
  while start + n <= m {
    let mut h = 0
    let mut i = 0
    while i < n {
      if read.unsafe_get(i) != seq.unsafe_get(start + i) {
        h = h + 1
      }
      i = i + 1
    }
    if h < best {
      best = h
    }
    start = start + 1
  }
  best
}

///|
/// 图上基因分型: 读段分派给滑窗最小 Hamming 的路径
/// (tie 取路径索引小者, 确定性), 计数前二(降序, tie 取索引小)
/// 构成二倍体等位; 路径 0(全 ref)为 "0".
/// GT: 两等位相同 -> 纯合 "x/x"; 不同 -> "x/y"(x <= y).
pub fn vg_genotype(paths : Array[VgPath], reads : Array[String]) -> VgGenotype? {
  if paths.length() == 0 {
    return None
  }
  let counts : Array[Int32] = Array::make(paths.length(), 0)
  for read in reads {
    if read.length() == 0 {
      continue
    }
    let mut best_idx = 0
    let mut best_h = vg_min_hamming(read, paths[0].seq)
    for pi = 1; pi < paths.length(); pi = pi + 1 {
      let h = vg_min_hamming(read, paths[pi].seq)
      if h < best_h {
        best_h = h
        best_idx = pi
      }
    }
    counts[best_idx] = counts[best_idx] + 1
  }
  // 前二等位(降序; tie 取索引小)
  let mut a = 0
  let mut b : Int32? = None
  for pi = 1; pi < paths.length(); pi = pi + 1 {
    if counts[pi] > counts[a] {
      b = Some(a)
      a = pi
    } else {
      match b {
        Some(bi) => if counts[pi] > counts[bi] { b = Some(pi) }
        None => b = Some(pi)
      }
    }
  }
  let second = match b {
    Some(x) => x
    None => a
  }
  // 次高等位计数为 0 时为纯合(仅一条路径有读段支持)
  let gt = if a == second || counts[second] == 0 {
    "\{a}/\{a}"
  } else {
    let lo = if a < second { a } else { second }
    let hi = if a < second { second } else { a }
    "\{lo}/\{hi}"
  }
  Some(VgGenotype::{
    gt,
    counts: [counts[a], counts[second]],
    path_counts: counts,
    n_reads: reads.length(),
  })
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 图 -> GFA 文本(S 节点 / L 边, 衔接 gfa_io 生态).
pub fn vg_to_gfa(graph : VgGraph) -> String {
  let sb = StringBuilder()
  sb.write_string("H\tVN:Z:1.0\n")
  for n in graph.nodes {
    sb.write_string("S\t\{n.id}\t\{n.seq}\n")
  }
  for e in graph.edges {
    sb.write_string("L\t\{e.from_id}\t+\t\{e.to_id}\t+\t*\n")
  }
  sb.to_string()
}

///|
/// 图 + 路径 + 分型 -> JSON 值.
pub fn vg_value(
  graph : VgGraph,
  paths : Array[VgPath],
  genotype : VgGenotype?,
) -> JsonValue {
  let nodes_arr : Array[JsonValue] = Array(capacity=graph.nodes.length())
  for n in graph.nodes {
    nodes_arr.push(
      JsonObj(
        Map([
          ("id", JsonNum(n.id.to_double())),
          ("seq", JsonStr(n.seq)),
          ("kind", JsonStr(n.kind)),
        ]),
      ),
    )
  }
  let edges_arr : Array[JsonValue] = Array(capacity=graph.edges.length())
  for e in graph.edges {
    edges_arr.push(
      JsonArr([JsonNum(e.from_id.to_double()), JsonNum(e.to_id.to_double())]),
    )
  }
  let paths_arr : Array[JsonValue] = Array(capacity=paths.length())
  for p in paths {
    paths_arr.push(JsonStr(p.seq))
  }
  let gt_obj : JsonValue = match genotype {
    Some(g) =>
      JsonObj(
        Map([
          ("gt", JsonStr(g.gt)),
          (
            "counts",
            JsonArr([
              JsonNum(g.counts[0].to_double()),
              JsonNum(g.counts[1].to_double()),
            ]),
          ),
          ("n_reads", JsonNum(g.n_reads.to_double())),
        ]),
      )
    None => JsonNull
  }
  JsonObj(
    Map([
      ("chrom", JsonStr(graph.chrom)),
      ("nodes", JsonArr(nodes_arr)),
      ("edges", JsonArr(edges_arr)),
      ("paths", JsonArr(paths_arr)),
      ("genotype", gt_obj),
    ]),
  )
}

///|
/// 摘要行.
pub fn vg_summary_line(graph : VgGraph, paths : Array[VgPath]) -> String {
  "vg-lite: nodes=\{graph.nodes.length()}, edges=\{graph.edges.length()}, paths=\{paths.length()}"
}
