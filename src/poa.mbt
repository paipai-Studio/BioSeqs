///|
/// POA(Partial Order Alignment)图多序列比对.
///
/// 将序列逐条比对并并入 DAG(Lee et al. 2002): 相同字符并入已有节点
/// (权重+1), 突变/插入产生平行路径节点. 相比 mafft/muscle 的扁平 MSA,
/// POA 天然表达重复/重组/indel 的分支结构, 适合多样本共识序列与
/// 长读长纠错场景.
///
/// 用法: `poa_add` 逐条并入(首条经 `poa_new` 建图), `poa_consensus`
/// 按拓扑列多数投票输出共识, `poa_msa` 输出各序列相对图列的对齐行.
pub struct PoaGraph {
  nodes : Array[PoaNode]
  sequences : Array[String]
  traces : Array[Array[Int32?]]
  mut last_score : Int32
}

///|
/// 图节点: 字符/支持权重/拓扑列/邻接.
pub struct PoaNode {
  letter : UInt16
  mut weight : Int32
  level : Int32
  prev : Array[Int32]
  next : Array[Int32]
}

///|
/// 默认打分: match +2, mismatch -1, gap -2.
pub fn PoaGraph::new(
  seed : String,
  match_score~ : Int32 = 2,
  mismatch_score~ : Int32 = -1,
  gap~ : Int32 = -2,
) -> PoaGraph {
  let graph = {
    nodes: Array(capacity=seed.length()),
    sequences: [seed],
    traces: [[]],
    last_score: 0,
  }
  let trace : Array[Int32?] = Array::make(seed.length(), None)
  let mut last : Int32? = None
  for j in 0..<seed.length() {
    let u = poa_new_node(
      graph,
      seed.unsafe_get(j),
      1,
      match last {
        Some(p) => [p]
        None => []
      },
      match last {
        Some(p) => graph.nodes[p].level + 1
        None => 0
      },
    )
    poa_link(graph, last, u)
    trace[j] = Some(u)
    last = Some(u)
  }
  graph.traces[0] = trace
  graph.last_score = seed.length() * match_score
  ignore(mismatch_score)
  ignore(gap)
  graph
}

///|
fn poa_new_node(
  graph : PoaGraph,
  letter : UInt16,
  weight : Int32,
  prev : Array[Int32],
  level : Int32,
) -> Int32 {
  let idx = graph.nodes.length()
  graph.nodes.push({ letter, weight, level, prev, next: [], })
  idx
}

///|
fn poa_link(graph : PoaGraph, a : Int32?, b : Int32) -> Unit {
  match a {
    Some(p) => {
      let has = {
        let mut found = false
        for x in graph.nodes[p].next {
          if x == b {
            found = true
          }
        }
        found
      }
      if !has {
        graph.nodes[p].next.push(b)
        graph.nodes[b].prev.push(p)
      }
    }
    None => ()
  }
}

///|
/// 将一条序列比对并并入图, 返回该序列与图的最佳对齐得分.
pub fn poa_add(
  graph : PoaGraph,
  seq : String,
  match_score~ : Int32 = 2,
  mismatch_score~ : Int32 = -1,
  gap~ : Int32 = -2,
) -> Int32 {
  let m = seq.length()
  let n = graph.nodes.length()
  let s_idx = graph.sequences.length()
  graph.sequences.push(seq)
  graph.traces.push([])
  if m == 0 {
    graph.last_score = 0
    return 0
  }
  // DP: dp[v+1][j], 行 0 为虚拟 START
  let dp : Array[Array[Int32]] = Array(capacity=n + 1)
  for _v in 0..<=n {
    dp.push(Array::make(m + 1, 0))
  }
  for j in 1..<=m {
    dp[0][j] = j * gap
  }
  // 节点按 id 即拓扑序(pred id 恒小于自身), 顺序填表
  for v in 0..<n {
    let row = v + 1
    let preds = graph.nodes[v].prev
    // dp[row][0]: 图首部 gap 累计
    dp[row][0] = if preds.length() == 0 {
      0
    } else {
      let mut best = preds[0]
      for p in preds {
        if dp[p + 1][0] > dp[best + 1][0] {
          best = p
        }
      }
      dp[best + 1][0] + gap
    }
    for j in 1..<=m {
      let c = seq.unsafe_get(j - 1)
      let sc = if graph.nodes[v].letter == c {
        match_score
      } else {
        mismatch_score
      }
      // diagonal: v 对齐 s[j-1]
      let diag = if preds.length() == 0 {
        dp[0][j - 1] + sc
      } else {
        let mut best = preds[0]
        for p in preds {
          if dp[p + 1][j - 1] > dp[best + 1][j - 1] {
            best = p
          }
        }
        dp[best + 1][j - 1] + sc
      }
      // 序列插入: s[j-1] 无图节点
      let ins_s = dp[row][j - 1] + gap
      // 图 gap: v 对齐 '-'
      let ins_g = if preds.length() == 0 {
        dp[0][j] + gap
      } else {
        let mut best = preds[0]
        for p in preds {
          if dp[p + 1][j] > dp[best + 1][j] {
            best = p
          }
        }
        dp[best + 1][j] + gap
      }
      let mut value = diag
      if ins_s > value {
        value = ins_s
      }
      if ins_g > value {
        value = ins_g
      }
      dp[row][j] = value
    }
  }
  // 终点: 序列全部消费, 图任意节点可为止点
  let mut best_v = -1
  let mut best_score = dp[0][m]
  for v in 0..<n {
    if dp[v + 1][m] > best_score {
      best_score = dp[v + 1][m]
      best_v = v
    }
  }
  // 回溯(不存回溯矩阵, 由 S 值反推), 路径元素: (节点?, 序列位置?)
  let path : Array[(Int32?, Int32?)] = Array(capacity=0)
  let mut cur_v = best_v
  let mut cur_j = m
  let mut steps = 0
  let steps_limit = n + m + 8
  while (cur_v >= 0 || cur_j > 0) && steps < steps_limit {
    steps = steps + 1
    if cur_v < 0 {
      // START 行: 剩余为首部插入
      path.push((None, Some(cur_j - 1)))
      cur_j = cur_j - 1
      continue
    }
    let v = cur_v
    let row = v + 1
    let preds = graph.nodes[v].prev
    if cur_j == 0 {
      // 图首部: v 对齐 gap, 沿最优 pred 回溯
      path.push((Some(v), None))
      if preds.length() == 0 {
        cur_v = -1
      } else {
        let mut best = preds[0]
        for p in preds {
          if dp[p + 1][0] > dp[best + 1][0] {
            best = p
          }
        }
        cur_v = best
      }
      continue
    }
    let c = seq.unsafe_get(cur_j - 1)
    let sc = if graph.nodes[v].letter == c {
      match_score
    } else {
      mismatch_score
    }
    let diag = if preds.length() == 0 {
      dp[0][cur_j - 1] + sc
    } else {
      let mut best = preds[0]
      for p in preds {
        if dp[p + 1][cur_j - 1] > dp[best + 1][cur_j - 1] {
          best = p
        }
      }
      dp[best + 1][cur_j - 1] + sc
    }
    let ins_s = dp[row][cur_j - 1] + gap
    let ins_g = if preds.length() == 0 {
      dp[0][cur_j] + gap
    } else {
      let mut best = preds[0]
      for p in preds {
        if dp[p + 1][cur_j] > dp[best + 1][cur_j] {
          best = p
        }
      }
      dp[best + 1][cur_j] + gap
    }
    let cur = dp[row][cur_j]
    if cur == diag {
      path.push((Some(v), Some(cur_j - 1)))
      if preds.length() == 0 {
        cur_v = -1
      } else {
        let mut best = preds[0]
        for p in preds {
          if dp[p + 1][cur_j - 1] > dp[best + 1][cur_j - 1] {
            best = p
          }
        }
        cur_v = best
      }
      cur_j = cur_j - 1
    } else if cur == ins_s {
      path.push((None, Some(cur_j - 1)))
      cur_j = cur_j - 1
    } else if cur == ins_g {
      path.push((Some(v), None))
      if preds.length() == 0 {
        cur_v = -1
      } else {
        let mut best = preds[0]
        for p in preds {
          if dp[p + 1][cur_j] > dp[best + 1][cur_j] {
            best = p
          }
        }
        cur_v = best
      }
    } else {
      // 数值不一致(不应发生), 强制终止
      path.push((Some(v), Some(cur_j - 1)))
      cur_j = cur_j - 1
    }
  }
  let path_fwd : Array[(Int32?, Int32?)] = Array(capacity=path.length())
  let mut k = path.length() - 1
  while k >= 0 {
    path_fwd.push(path[k])
    k = k - 1
  }
  // 沿路径建图: match 并入, mismatch/插入新建节点
  let trace : Array[Int32?] = Array::make(m, None)
  let mut last : Int32? = None
  for step in path_fwd {
    let (vopt, jopt) = step
    match (vopt, jopt) {
      (Some(v), Some(j)) => {
        let c = seq.unsafe_get(j)
        let u = if graph.nodes[v].letter == c {
          graph.nodes[v].weight = graph.nodes[v].weight + 1
          v
        } else {
          poa_new_node(
            graph,
            c,
            1,
            match last {
              Some(p) => [p]
              None => []
            },
            graph.nodes[v].level,
          )
        }
        poa_link(graph, last, u)
        trace[j] = Some(u)
        last = Some(u)
      }
      (None, Some(j)) => {
        let c = seq.unsafe_get(j)
        let u = poa_new_node(
          graph,
          c,
          1,
          match last {
            Some(p) => [p]
            None => []
          },
          match last {
            Some(p) => graph.nodes[p].level + 1
            None => 0
          },
        )
        poa_link(graph, last, u)
        trace[j] = Some(u)
        last = Some(u)
      }
      (_, None) => ()
    }
  }
  graph.traces[s_idx] = trace
  graph.last_score = best_score
  best_score
}

///|
/// 图节点数.
pub fn poa_node_count(graph : PoaGraph) -> Int32 {
  graph.nodes.length()
}

///|
/// 最近一次 poa_add 的对齐得分.
pub fn poa_last_score(graph : PoaGraph) -> Int32 {
  graph.last_score
}

///|
/// 共识序列: 按拓扑列取支持权重最高的节点字符(并列取先建者).
pub fn poa_consensus(graph : PoaGraph) -> String {
  if graph.nodes.length() == 0 {
    return ""
  }
  let mut max_level = 0
  for node in graph.nodes {
    if node.level > max_level {
      max_level = node.level
    }
  }
  let buckets : Array[Array[Int32]] = Array::make(max_level + 1, [])
  for i in 0..<=max_level {
    buckets[i] = []
  }
  for i, node in graph.nodes {
    buckets[node.level].push(i)
  }
  let sb = StringBuilder()
  for bucket in buckets {
    let mut best = bucket[0]
    for idx in bucket {
      if graph.nodes[idx].weight > graph.nodes[best].weight {
        best = idx
      }
    }
    sb.write_char(graph.nodes[best].letter.unsafe_to_char())
  }
  sb.to_string()
}

///|
/// 输出 MSA: 每条序列相对图拓扑列的对齐行('-' 为 gap).
pub fn poa_msa(graph : PoaGraph) -> Array[String] {
  let mut max_level = -1
  for node in graph.nodes {
    if node.level > max_level {
      max_level = node.level
    }
  }
  let ncols = max_level + 1
  let rows : Array[String] = Array(capacity=graph.sequences.length())
  for s_idx, seq in graph.sequences {
    let row : Array[UInt16] = Array::make(ncols, '-'.to_int().to_uint16())
    let trace = graph.traces[s_idx]
    for j in 0..<seq.length() {
      match trace[j] {
        Some(u) => {
          let lv = graph.nodes[u].level
          if lv >= 0 && lv < ncols {
            row[lv] = seq.unsafe_get(j)
          }
        }
        None => ()
      }
    }
    let sb = StringBuilder()
    for c in row {
      sb.write_char(c.unsafe_to_char())
    }
    rows.push(sb.to_string())
  }
  rows
}
