///|
/// assembly_graph.mbt — 组装图简化与共识抛光 (racon-lite / 图清理).
///
/// 两部分功能:
///
/// 1. **GFA 图简化** (assembly_simplify):
///    - containment 坍缩: C 记录中被包含的 segment 视为冗余并移除;
///    - 短死端裁剪 (tip removal): 入度或出度为 0 且长度低于阈值的末端
///      segment 迭代移除, 直到图稳定;
///    - 简单 bubble 检测与解析: 入口节点两条互不相交的单出边分支在
///      同一出口汇合, 按分支支持度 (DP/KC 标签之和) 与长度选择保留分支.
///
///    拓扑分析把 GFA link 视作有向边 `from -> to`, 不翻转序列方向
///    (简化级实现; 需要完整 GFA 方向语义的场景请先规范化 link).
///
/// 2. **参考导向共识抛光** (racon_polish):
///    - k-mer 锚定把每条读段定位到参考坐标 (正链锚定弱时尝试反向互补);
///    - 在参考窗口内做 Needleman-Wunsch 比对, 逐列堆叠读段碱基;
///    - 逐位点多数投票: 替换/删除需达到最小覆盖度与最小支持比例,
///      插入槽按插入串整体计数, 避免少量读段的随机插入被采纳.
///
///    抛光只修改有足够读段支持的位点, 无覆盖区域保持原参考碱基.

// ============================================================================
// 错误与配置
// ============================================================================

///|
/// 组装图简化 / 抛光参数错误.
pub suberror AssemblyGraphError {
  AssemblyGraphError(String)
}

///|
/// 图简化参数.
pub struct AssemblySimplifyConfig {
  /// 短于该长度的死端 segment 会被裁剪.
  min_tip_len : Int32
  /// bubble 两条分支总长上限 (超过视为复杂结构, 保留).
  max_bubble_len : Int32
  /// 是否移除 containment 中被包含的 segment.
  remove_contained : Bool
  /// 是否解析简单 bubble.
  resolve_bubbles : Bool
  /// 解析 bubble 时落选分支的最大支持度; 落选分支支持度不得高于该值.
  /// 设为负值表示不限制 (始终按支持度比较).
  min_bubble_support : Double
}

///|
/// 构造图简化参数并校验.
pub fn assembly_simplify_config(
  min_tip_len~ : Int32 = 100,
  max_bubble_len~ : Int32 = 1000,
  remove_contained~ : Bool = true,
  resolve_bubbles~ : Bool = true,
  min_bubble_support~ : Double = -1.0,
) -> AssemblySimplifyConfig raise AssemblyGraphError {
  if min_tip_len < 0 {
    raise AssemblyGraphError("min_tip_len must be non-negative")
  }
  if max_bubble_len <= 0 {
    raise AssemblyGraphError("max_bubble_len must be positive")
  }
  {
    min_tip_len,
    max_bubble_len,
    remove_contained,
    resolve_bubbles,
    min_bubble_support,
  }
}

///|
/// 图简化统计.
pub struct AssemblySimplifyStats {
  /// 简化前 segment 数.
  segments_before : Int32
  /// 简化后 segment 数.
  segments_after : Int32
  /// 简化前 link 数.
  links_before : Int32
  /// 简化后 link 数.
  links_after : Int32
  /// 裁剪的短死端 segment 数.
  tips_removed : Int32
  /// 坍缩的被包含 segment 数.
  contained_removed : Int32
  /// 检测到的简单 bubble 数.
  bubbles_detected : Int32
  /// 已解析 (删除落选分支) 的 bubble 数.
  bubbles_resolved : Int32
  /// 因引用被删 segment 而丢弃的 path 数.
  paths_dropped : Int32
}

///|
pub fn AssemblySimplifyStats::segments_before(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.segments_before
}

///|
pub fn AssemblySimplifyStats::segments_after(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.segments_after
}

///|
pub fn AssemblySimplifyStats::links_before(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.links_before
}

///|
pub fn AssemblySimplifyStats::links_after(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.links_after
}

///|
pub fn AssemblySimplifyStats::tips_removed(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.tips_removed
}

///|
pub fn AssemblySimplifyStats::contained_removed(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.contained_removed
}

///|
pub fn AssemblySimplifyStats::bubbles_detected(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.bubbles_detected
}

///|
pub fn AssemblySimplifyStats::bubbles_resolved(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.bubbles_resolved
}

///|
pub fn AssemblySimplifyStats::paths_dropped(
  self : AssemblySimplifyStats,
) -> Int32 {
  self.paths_dropped
}

///|
/// 简化结果: 新图 + 统计.
pub struct AssemblySimplification {
  graph : GfaGraph
  stats : AssemblySimplifyStats
}

///|
pub fn AssemblySimplification::graph(self : AssemblySimplification) -> GfaGraph {
  self.graph
}

///|
pub fn AssemblySimplification::stats(
  self : AssemblySimplification,
) -> AssemblySimplifyStats {
  self.stats
}

// ============================================================================
// 图拓扑辅助
// ============================================================================

///|
/// segment 长度: 有真实序列时取序列长度, 否则读取 LN:i 标签.
pub fn assembly_segment_length(segment : GfaSegment) -> Int32 {
  if segment.sequence != "*" {
    return segment.sequence.length()
  }
  for tag in segment.tags {
    if tag.name == "LN" {
      return assembly_parse_int(tag.value)
    }
  }
  0
}

///|
fn assembly_parse_int(s : String) -> Int32 {
  let mut value = 0
  let mut negative = false
  let mut i = 0
  if s.length() > 0 && s.unsafe_get(0).unsafe_to_char() == '-' {
    negative = true
    i = 1
  }
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c >= '0' && c <= '9' {
      value = value * 10 + (c.to_int() - '0'.to_int())
    }
    i = i + 1
  }
  if negative {
    -value
  } else {
    value
  }
}

///|
/// segment 深度: 优先读取 DP:i 标签, 其次 KC:i, 缺省为 1.0.
pub fn assembly_segment_depth(segment : GfaSegment) -> Double {
  for tag in segment.tags {
    if tag.name == "DP" {
      return assembly_parse_double(tag.value)
    }
  }
  for tag in segment.tags {
    if tag.name == "KC" {
      return assembly_parse_double(tag.value)
    }
  }
  1.0
}

///|
fn assembly_parse_double(s : String) -> Double {
  let mut value = 0.0
  let mut scale = 0.1
  let mut negative = false
  let mut in_fraction = false
  let mut i = 0
  if s.length() > 0 && s.unsafe_get(0).unsafe_to_char() == '-' {
    negative = true
    i = 1
  }
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c == '.' {
      in_fraction = true
    } else if c >= '0' && c <= '9' {
      let digit = (c.to_int() - '0'.to_int()).to_double()
      if in_fraction {
        value = value + digit * scale
        scale = scale * 0.1
      } else {
        value = value * 10.0 + digit
      }
    }
    i = i + 1
  }
  if negative {
    -value
  } else {
    value
  }
}

///|
fn assembly_live_map(graph : GfaGraph) -> Map[String, Bool] {
  let alive : Map[String, Bool] = Map([])
  for segment in graph.segments {
    alive.set(segment.name, true)
  }
  alive
}

///|
fn assembly_is_alive(alive : Map[String, Bool], name : String) -> Bool {
  match alive.get(name) {
    Some(flag) => flag
    None => false
  }
}

///|
fn assembly_out_neighbors(
  graph : GfaGraph,
  alive : Map[String, Bool],
  name : String,
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  for link in graph.links {
    if link.from_segment == name &&
      assembly_is_alive(alive, link.from_segment) &&
      assembly_is_alive(alive, link.to_segment) {
      match seen.get(link.to_segment) {
        Some(_) => ()
        None => {
          seen.set(link.to_segment, true)
          out.push(link.to_segment)
        }
      }
    }
  }
  out
}

///|
fn assembly_degrees(
  graph : GfaGraph,
  alive : Map[String, Bool],
) -> (Map[String, Int32], Map[String, Int32]) {
  let indeg : Map[String, Int32] = Map([])
  let outdeg : Map[String, Int32] = Map([])
  for segment in graph.segments {
    if assembly_is_alive(alive, segment.name) {
      indeg.set(segment.name, 0)
      outdeg.set(segment.name, 0)
    }
  }
  for link in graph.links {
    if assembly_is_alive(alive, link.from_segment) &&
      assembly_is_alive(alive, link.to_segment) {
      let current_out = match outdeg.get(link.from_segment) {
        Some(n) => n
        None => 0
      }
      outdeg.set(link.from_segment, current_out + 1)
      let current_in = match indeg.get(link.to_segment) {
        Some(n) => n
        None => 0
      }
      indeg.set(link.to_segment, current_in + 1)
    }
  }
  (indeg, outdeg)
}

///|
fn assembly_live_count(alive : Map[String, Bool]) -> Int32 {
  let mut count = 0
  for _name in alive.keys() {
    count = count + 1
  }
  count
}

///|
fn assembly_segment_by_name(graph : GfaGraph, name : String) -> GfaSegment? {
  for segment in graph.segments {
    if segment.name == name {
      return Some(segment)
    }
  }
  None
}

// ============================================================================
// 1. containment 坍缩
// ============================================================================

///|
fn assembly_remove_contained(
  graph : GfaGraph,
  alive : Map[String, Bool],
) -> Int32 {
  let mut removed = 0
  for containment in graph.containments {
    if assembly_is_alive(alive, containment.contained) {
      alive.set(containment.contained, false)
      removed = removed + 1
    }
  }
  if removed >= assembly_live_count(alive) {
    // 防止整图被坍缩为空: 回滚本步.
    for containment in graph.containments {
      alive.set(containment.contained, true)
    }
    return 0
  }
  removed
}

// ============================================================================
// 2. 短死端裁剪
// ============================================================================

///|
fn assembly_trim_tips(
  graph : GfaGraph,
  alive : Map[String, Bool],
  min_tip_len : Int32,
) -> Int32 {
  let mut removed = 0
  let mut changed = true
  while changed {
    changed = false
    let (indeg, outdeg) = assembly_degrees(graph, alive)
    let candidates : Array[String] = Array(capacity=0)
    for segment in graph.segments {
      if !assembly_is_alive(alive, segment.name) {
        continue
      }
      let d_in = match indeg.get(segment.name) {
        Some(n) => n
        None => 0
      }
      let d_out = match outdeg.get(segment.name) {
        Some(n) => n
        None => 0
      }
      let isolated = d_in == 0 && d_out == 0
      let leaf = (d_in == 0 && d_out == 1) || (d_in == 1 && d_out == 0)
      let mut is_tip = isolated
      if leaf && !isolated {
        // 仅当唯一邻居是分叉点 (总度 >= 3) 时才是旁支 tip;
        // 线性链末端视为 contig 端点, 不裁剪.
        let mut neighbor = ""
        for link in graph.links {
          if assembly_is_alive(alive, link.from_segment) &&
            assembly_is_alive(alive, link.to_segment) {
            if link.from_segment == segment.name {
              neighbor = link.to_segment
            } else if link.to_segment == segment.name {
              neighbor = link.from_segment
            }
          }
        }
        if neighbor != "" {
          let n_in = match indeg.get(neighbor) {
            Some(n) => n
            None => 0
          }
          let n_out = match outdeg.get(neighbor) {
            Some(n) => n
            None => 0
          }
          if n_in + n_out >= 3 {
            is_tip = true
          }
        }
      }
      if is_tip && assembly_segment_length(segment) < min_tip_len {
        candidates.push(segment.name)
      }
    }
    if candidates.length() == 0 ||
      candidates.length() >= assembly_live_count(alive) {
      break
    }
    for name in candidates {
      if assembly_is_alive(alive, name) {
        alive.set(name, false)
        removed = removed + 1
        changed = true
      }
    }
  }
  removed
}

// ============================================================================
// 3. 简单 bubble 检测与解析
// ============================================================================

///|
/// 一个简单 bubble: 入口 -> 两条互不相交分支 -> 同一出口.
pub struct AssemblyBubble {
  entry : String
  exit : String
  branch_a : Array[String]
  branch_b : Array[String]
  support_a : Double
  support_b : Double
}

///|
/// bubble 入口 segment 名.
pub fn AssemblyBubble::entry(self : AssemblyBubble) -> String {
  self.entry
}

///|
/// bubble 出口 segment 名.
pub fn AssemblyBubble::exit(self : AssemblyBubble) -> String {
  self.exit
}

///|
/// 分支 A 的节点序列 (含出口).
pub fn AssemblyBubble::branch_a(self : AssemblyBubble) -> Array[String] {
  self.branch_a.copy()
}

///|
/// 分支 B 的节点序列 (含出口).
pub fn AssemblyBubble::branch_b(self : AssemblyBubble) -> Array[String] {
  self.branch_b.copy()
}

///|
/// 分支 A 支持度 (内部节点 DP/KC 之和).
pub fn AssemblyBubble::support_a(self : AssemblyBubble) -> Double {
  self.support_a
}

///|
/// 分支 B 支持度 (内部节点 DP/KC 之和).
pub fn AssemblyBubble::support_b(self : AssemblyBubble) -> Double {
  self.support_b
}

///|
fn assembly_path_contains(path : Array[String], name : String) -> Bool {
  for node in path {
    if node == name {
      return true
    }
  }
  false
}

///|
fn assembly_trace_branch(
  graph : GfaGraph,
  alive : Map[String, Bool],
  start : String,
  max_steps : Int32,
) -> Array[String]? {
  let (indeg, outdeg) = assembly_degrees(graph, alive)
  let path : Array[String] = Array(capacity=0)
  path.push(start)
  let mut current = start
  let mut steps = 1
  while steps <= max_steps {
    if !assembly_is_alive(alive, current) {
      return None
    }
    let d_out = match outdeg.get(current) {
      Some(n) => n
      None => 0
    }
    if d_out != 1 {
      return None
    }
    let nexts = assembly_out_neighbors(graph, alive, current)
    if nexts.length() != 1 {
      return None
    }
    let next = nexts[0]
    if assembly_path_contains(path, next) {
      return None
    }
    let d_in_next = match indeg.get(next) {
      Some(n) => n
      None => 0
    }
    if d_in_next != 1 {
      path.push(next)
      return Some(path)
    }
    path.push(next)
    current = next
    steps = steps + 1
  }
  None
}

///|
fn assembly_path_length(graph : GfaGraph, path : Array[String]) -> Int32 {
  let mut total = 0
  for name in path {
    match assembly_segment_by_name(graph, name) {
      Some(segment) => total = total + assembly_segment_length(segment)
      None => ()
    }
  }
  total
}

///|
fn assembly_path_support(graph : GfaGraph, path : Array[String]) -> Double {
  let mut total = 0.0
  // 末端节点为出口, 不计入分支支持度.
  let n = path.length()
  for i in 0..<(n - 1) {
    match assembly_segment_by_name(graph, path[i]) {
      Some(segment) => total = total + assembly_segment_depth(segment)
      None => ()
    }
  }
  total
}

///|
fn assembly_bubble_key(
  entry : String,
  exit : String,
  a : Array[String],
  b : Array[String],
) -> String {
  // 两条分支的无序集合表示, 保证同一 bubble 只记录一次.
  let names : Array[String] = Array(capacity=0)
  for i in 0..<(a.length() - 1) {
    names.push(a[i])
  }
  for i in 0..<(b.length() - 1) {
    names.push(b[i])
  }
  names.push(entry)
  names.push(exit)
  let mut sorted = names.copy()
  sorted = assembly_sort_strings(sorted)
  let sb = StringBuilder()
  for i, name in sorted {
    if i > 0 {
      sb.write_char('|')
    }
    sb.write_string(name)
  }
  sb.to_string()
}

///|
fn assembly_sort_strings(items : Array[String]) -> Array[String] {
  let out = items.copy()
  let n = out.length()
  for i in 1..<n {
    let value = out[i]
    let mut j = i
    while j > 0 && out[j - 1] > value {
      out[j] = out[j - 1]
      j = j - 1
    }
    out[j] = value
  }
  out
}

///|
/// 检测图中的简单 bubble.
pub fn assembly_detect_bubbles(
  graph : GfaGraph,
  max_bubble_len~ : Int32 = 1000,
) -> Array[AssemblyBubble] {
  assembly_detect_bubbles_alive(graph, assembly_live_map(graph), max_bubble_len)
}

///|
fn assembly_detect_bubbles_alive(
  graph : GfaGraph,
  alive : Map[String, Bool],
  max_bubble_len : Int32,
) -> Array[AssemblyBubble] {
  let (indeg, _outdeg) = assembly_degrees(graph, alive)
  let bubbles : Array[AssemblyBubble] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  let max_steps = graph.segments.length() + 1
  for entry in graph.segments {
    if !assembly_is_alive(alive, entry.name) {
      continue
    }
    let succs = assembly_out_neighbors(graph, alive, entry.name)
    if succs.length() != 2 {
      continue
    }
    let s1 = succs[0]
    let s2 = succs[1]
    let d_in_1 = match indeg.get(s1) {
      Some(n) => n
      None => 0
    }
    let d_in_2 = match indeg.get(s2) {
      Some(n) => n
      None => 0
    }
    if d_in_1 != 1 || d_in_2 != 1 {
      continue
    }
    let traced1 = assembly_trace_branch(graph, alive, s1, max_steps)
    let traced2 = assembly_trace_branch(graph, alive, s2, max_steps)
    match (traced1, traced2) {
      (Some(p1), Some(p2)) => {
        if p1.length() < 2 || p2.length() < 2 {
          continue
        }
        let exit = p1[p1.length() - 1]
        if exit != p2[p2.length() - 1] || exit == entry.name {
          continue
        }
        // 分支内部节点不得相交.
        let mut disjoint = true
        for i in 0..<(p1.length() - 1) {
          for j in 0..<(p2.length() - 1) {
            if p1[i] == p2[j] {
              disjoint = false
            }
          }
        }
        if !disjoint {
          continue
        }
        let exit_len = match assembly_segment_by_name(graph, exit) {
          Some(segment) => assembly_segment_length(segment)
          None => 0
        }
        let total_len = assembly_path_length(graph, p1) +
          assembly_path_length(graph, p2) -
          exit_len
        if total_len > max_bubble_len {
          continue
        }
        let key = assembly_bubble_key(entry.name, exit, p1, p2)
        match seen.get(key) {
          Some(_) => continue
          None => ()
        }
        seen.set(key, true)
        bubbles.push({
          entry: entry.name,
          exit,
          branch_a: p1,
          branch_b: p2,
          support_a: assembly_path_support(graph, p1),
          support_b: assembly_path_support(graph, p2),
        })
      }
      _ => ()
    }
  }
  bubbles
}

///|
/// 解析 bubble: 删除落选分支的内部节点 (入口与出口保留).
///
/// 选择规则: 支持度高者胜出; 支持度相同时更长分支胜出; 仍相同时保留 A.
pub fn assembly_resolve_bubble(
  graph : GfaGraph,
  alive : Map[String, Bool],
  bubble : AssemblyBubble,
  min_bubble_support~ : Double = -1.0,
) -> Array[String] {
  let keep_a = if bubble.support_a > bubble.support_b {
    true
  } else if bubble.support_b > bubble.support_a {
    false
  } else {
    assembly_path_length(graph, bubble.branch_a) >=
    assembly_path_length(graph, bubble.branch_b)
  }
  let loser = if keep_a { bubble.branch_b } else { bubble.branch_a }
  let loser_support = if keep_a { bubble.support_b } else { bubble.support_a }
  let removed : Array[String] = Array(capacity=0)
  if min_bubble_support >= 0.0 && loser_support > min_bubble_support {
    return removed
  }
  for i in 0..<(loser.length() - 1) {
    let name = loser[i]
    if assembly_is_alive(alive, name) {
      alive.set(name, false)
      removed.push(name)
    }
  }
  removed
}

// ============================================================================
// 4. 图构建与主入口
// ============================================================================

///|
fn assembly_path_names(path : GfaPath) -> Array[String] {
  let names : Array[String] = Array(capacity=path.segment_list.length())
  for item in path.segment_list {
    if item.length() > 0 {
      let last = item.unsafe_get(item.length() - 1).unsafe_to_char()
      if last == '+' || last == '-' {
        names.push(item[0:item.length() - 1].to_owned())
      } else {
        names.push(item)
      }
    } else {
      names.push(item)
    }
  }
  names
}

///|
fn assembly_path_has_removed(path : GfaPath, alive : Map[String, Bool]) -> Bool {
  for name in assembly_path_names(path) {
    if !assembly_is_alive(alive, name) {
      return true
    }
  }
  false
}

///|
fn assembly_build_graph(
  graph : GfaGraph,
  alive : Map[String, Bool],
) -> (GfaGraph, Int32) {
  let out = GfaGraph::new()
  for header in graph.headers {
    out.add_header(header)
  }
  for segment in graph.segments {
    if assembly_is_alive(alive, segment.name) {
      out.add_segment(segment)
    }
  }
  for link in graph.links {
    if assembly_is_alive(alive, link.from_segment) &&
      assembly_is_alive(alive, link.to_segment) {
      out.add_link(link)
    }
  }
  for containment in graph.containments {
    if assembly_is_alive(alive, containment.container) &&
      assembly_is_alive(alive, containment.contained) {
      out.add_containment(containment)
    }
  }
  let mut dropped = 0
  for path in graph.paths {
    if assembly_path_has_removed(path, alive) {
      dropped = dropped + 1
    } else {
      out.add_path(path)
    }
  }
  (out, dropped)
}

///|
/// 简化组装图: containment 坍缩 -> 短死端裁剪 -> 简单 bubble 解析 ->
/// 再次裁剪 (bubble 删除可能产生新的短死端).
///
/// 返回新的 GfaGraph, 原图不变. 若某一步会清空整图, 该步被跳过.
pub fn assembly_simplify(
  graph : GfaGraph,
  config~ : AssemblySimplifyConfig = {
    min_tip_len: 100,
    max_bubble_len: 1000,
    remove_contained: true,
    resolve_bubbles: true,
    min_bubble_support: -1.0,
  },
) -> AssemblySimplification {
  let alive = assembly_live_map(graph)
  let mut contained_removed = 0
  if config.remove_contained {
    contained_removed = assembly_remove_contained(graph, alive)
  }
  let mut tips_removed = assembly_trim_tips(graph, alive, config.min_tip_len)
  let mut bubbles_detected = 0
  let mut bubbles_resolved = 0
  if config.resolve_bubbles {
    let bubbles = assembly_detect_bubbles_alive(
      graph,
      alive,
      config.max_bubble_len,
    )
    for bubble in bubbles {
      bubbles_detected = bubbles_detected + 1
      let removed = assembly_resolve_bubble(
        graph,
        alive,
        bubble,
        min_bubble_support=config.min_bubble_support,
      )
      if removed.length() > 0 {
        bubbles_resolved = bubbles_resolved + 1
      }
    }
    tips_removed = tips_removed +
      assembly_trim_tips(graph, alive, config.min_tip_len)
  }
  let (simplified, paths_dropped) = assembly_build_graph(graph, alive)
  let stats = {
    segments_before: graph.segments.length(),
    segments_after: simplified.segments.length(),
    links_before: graph.links.length(),
    links_after: simplified.links.length(),
    tips_removed,
    contained_removed,
    bubbles_detected,
    bubbles_resolved,
    paths_dropped,
  }
  { graph: simplified, stats, }
}

///|
/// 图简化统计 JSON 输出.
pub fn assembly_simplify_stats_to_json(stats : AssemblySimplifyStats) -> String {
  let sb = StringBuilder()
  sb.write_string("{")
  sb.write_string("\"segments_before\":")
  sb.write_string(stats.segments_before.to_string())
  sb.write_string(",\"segments_after\":")
  sb.write_string(stats.segments_after.to_string())
  sb.write_string(",\"links_before\":")
  sb.write_string(stats.links_before.to_string())
  sb.write_string(",\"links_after\":")
  sb.write_string(stats.links_after.to_string())
  sb.write_string(",\"tips_removed\":")
  sb.write_string(stats.tips_removed.to_string())
  sb.write_string(",\"contained_removed\":")
  sb.write_string(stats.contained_removed.to_string())
  sb.write_string(",\"bubbles_detected\":")
  sb.write_string(stats.bubbles_detected.to_string())
  sb.write_string(",\"bubbles_resolved\":")
  sb.write_string(stats.bubbles_resolved.to_string())
  sb.write_string(",\"paths_dropped\":")
  sb.write_string(stats.paths_dropped.to_string())
  sb.write_string("}")
  sb.to_string()
}

// ============================================================================
// 5. racon-lite 共识抛光
// ============================================================================

///|
/// 共识抛光参数.
pub struct PolishConfig {
  /// 锚定 k-mer 长度 (参考更短时自动缩小).
  k : Int32
  /// 采纳修正所需的最小位点覆盖度.
  min_coverage : Int32
  /// 采纳修正所需的最小支持比例 (0~1).
  min_fraction : Double
  /// 参考窗口两侧余量.
  flank : Int32
}

///|
/// 构造抛光参数并校验.
pub fn polish_config(
  k~ : Int32 = 15,
  min_coverage~ : Int32 = 3,
  min_fraction~ : Double = 0.6,
  flank~ : Int32 = 50,
) -> PolishConfig raise AssemblyGraphError {
  if k < 4 {
    raise AssemblyGraphError("k must be at least 4")
  }
  if min_coverage < 1 {
    raise AssemblyGraphError("min_coverage must be at least 1")
  }
  if min_fraction <= 0.0 || min_fraction > 1.0 {
    raise AssemblyGraphError("min_fraction must be in (0, 1]")
  }
  if flank < 0 {
    raise AssemblyGraphError("flank must be non-negative")
  }
  { k, min_coverage, min_fraction, flank, }
}

///|
/// 单条修正记录.
pub struct PolishEdit {
  /// 修正类型: substitution / deletion / insertion.
  kind : String
  /// 参考坐标 (0-based); 插入记录为其插入槽位置.
  position : Int32
  /// 参考碱基 (插入为空串).
  ref_base : String
  /// 修正后碱基 (删除为空串).
  new_base : String
  /// 支持该修正的读段数.
  support : Int32
  /// 该位点/插入槽的覆盖度.
  coverage : Int32
}

///|
pub fn PolishEdit::kind(self : PolishEdit) -> String {
  self.kind
}

///|
pub fn PolishEdit::position(self : PolishEdit) -> Int32 {
  self.position
}

///|
pub fn PolishEdit::ref_base(self : PolishEdit) -> String {
  self.ref_base
}

///|
pub fn PolishEdit::new_base(self : PolishEdit) -> String {
  self.new_base
}

///|
pub fn PolishEdit::support(self : PolishEdit) -> Int32 {
  self.support
}

///|
pub fn PolishEdit::coverage(self : PolishEdit) -> Int32 {
  self.coverage
}

///|
/// 抛光结果.
pub struct PolishResult {
  /// 原参考序列.
  reference : String
  /// 抛光后序列.
  consensus : String
  /// 采纳的修正列表 (按坐标升序).
  edits : Array[PolishEdit]
  /// 成功锚定并参与投票的读段数.
  reads_used : Int32
  /// 输入读段总数.
  reads_total : Int32
  /// 参考各位点平均覆盖度.
  mean_coverage : Double
}

///|
pub fn PolishResult::reference(self : PolishResult) -> String {
  self.reference
}

///|
pub fn PolishResult::consensus(self : PolishResult) -> String {
  self.consensus
}

///|
pub fn PolishResult::edits(self : PolishResult) -> Array[PolishEdit] {
  self.edits.copy()
}

///|
pub fn PolishResult::reads_used(self : PolishResult) -> Int32 {
  self.reads_used
}

///|
pub fn PolishResult::reads_total(self : PolishResult) -> Int32 {
  self.reads_total
}

///|
pub fn PolishResult::mean_coverage(self : PolishResult) -> Double {
  self.mean_coverage
}

///|
fn polish_base_index(c : Char) -> Int32 {
  if c == 'A' || c == 'a' {
    0
  } else if c == 'C' || c == 'c' {
    1
  } else if c == 'G' || c == 'g' {
    2
  } else if c == 'T' || c == 't' {
    3
  } else {
    -1
  }
}

///|
fn polish_base_char(idx : Int32) -> Char {
  if idx == 0 {
    'A'
  } else if idx == 1 {
    'C'
  } else if idx == 2 {
    'G'
  } else {
    'T'
  }
}

///|
fn polish_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let comp = if c == 'A' || c == 'a' {
      'T'
    } else if c == 'C' || c == 'c' {
      'G'
    } else if c == 'G' || c == 'g' {
      'C'
    } else if c == 'T' || c == 't' {
      'A'
    } else {
      'N'
    }
    sb.write_char(comp)
    i = i - 1
  }
  sb.to_string()
}

///|
fn polish_build_index(
  reference : String,
  k : Int32,
) -> Map[String, Array[Int32]] {
  let index : Map[String, Array[Int32]] = Map([])
  let n = reference.length()
  let mut i = 0
  while i + k <= n {
    let kmer = reference[i:i + k].to_owned()
    match index.get(kmer) {
      Some(positions) => positions.push(i)
      None => index.set(kmer, [i])
    }
    i = i + 1
  }
  index
}

///|
fn polish_count_offset(
  index : Map[String, Array[Int32]],
  seq : String,
  k : Int32,
) -> (Int32?, Int32) {
  let votes : Map[Int32, Int32] = Map([])
  let n = seq.length()
  let mut i = 0
  while i + k <= n {
    let kmer = seq[i:i + k].to_owned()
    match index.get(kmer) {
      Some(positions) =>
        for pos in positions {
          let delta = pos - i
          let current = match votes.get(delta) {
            Some(c) => c
            None => 0
          }
          votes.set(delta, current + 1)
        }
      None => ()
    }
    i = i + 1
  }
  let mut best_delta : Int32? = None
  let mut best_votes = 0
  for delta in votes.keys() {
    let count = match votes.get(delta) {
      Some(c) => c
      None => 0
    }
    let better = if count > best_votes {
      true
    } else if count == best_votes {
      match best_delta {
        Some(current) => delta < current
        None => true
      }
    } else {
      false
    }
    if better {
      best_votes = count
      best_delta = Some(delta)
    }
  }
  (best_delta, best_votes)
}

///|
fn polish_reverse_string(s : String) -> String {
  let sb = StringBuilder()
  let mut i = s.length() - 1
  while i >= 0 {
    sb.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  sb.to_string()
}

///|
/// 半全局比对: 读段整体对齐到参考窗口, 窗口两端未覆盖部分不罚分.
///
/// 与全局 Needleman-Wunsch 的差别: 第一行得分全为 0 (窗口前缀免费),
/// 回溯从最后一行的最高分位置开始 (窗口后缀免费). 这样参考窗口两侧的
/// 余量不会把读段末尾错配到随机序列上, 避免伪插入/伪缺失.
fn polish_semiglobal_align(read : String, window : String) -> (String, String) {
  let n = read.length()
  let m = window.length()
  if n == 0 {
    return ("", window)
  }
  if m == 0 {
    return (read, "")
  }
  let gap_penalty = -2
  let dp : Array[Array[Int32]] = Array(capacity=n + 1)
  for _i in 0..<=n {
    dp.push(Array::make(m + 1, 0))
  }
  for i in 1..<=n {
    dp[i][0] = i * gap_penalty
  }
  for i in 1..<=n {
    for j in 1..<=m {
      let a = read.unsafe_get(i - 1).unsafe_to_char()
      let b = window.unsafe_get(j - 1).unsafe_to_char()
      let sub = if a == b { 2 } else { -1 }
      let diag = dp[i - 1][j - 1] + sub
      let up = dp[i - 1][j] + gap_penalty
      let left = dp[i][j - 1] + gap_penalty
      dp[i][j] = if diag >= up && diag >= left {
        diag
      } else if up >= left {
        up
      } else {
        left
      }
    }
  }
  let mut best_j = 0
  let mut best_score = dp[n][0]
  for j in 1..<=m {
    if dp[n][j] > best_score {
      best_score = dp[n][j]
      best_j = j
    }
  }
  let read_mid = StringBuilder()
  let ref_mid = StringBuilder()
  let mut i = n
  let mut j = best_j
  while i > 0 {
    if j == 0 {
      read_mid.write_char(read.unsafe_get(i - 1).unsafe_to_char())
      ref_mid.write_char('-')
      i = i - 1
    } else {
      let a = read.unsafe_get(i - 1).unsafe_to_char()
      let b = window.unsafe_get(j - 1).unsafe_to_char()
      let sub = if a == b { 2 } else { -1 }
      if dp[i][j] == dp[i - 1][j - 1] + sub {
        read_mid.write_char(a)
        ref_mid.write_char(b)
        i = i - 1
        j = j - 1
      } else if dp[i][j] == dp[i - 1][j] + gap_penalty {
        read_mid.write_char(a)
        ref_mid.write_char('-')
        i = i - 1
      } else {
        read_mid.write_char('-')
        ref_mid.write_char(b)
        j = j - 1
      }
    }
  }
  let prefix = window[0:j].to_owned()
  let suffix = window[best_j:m].to_owned()
  let read_aligned = StringBuilder()
  for _k in 0..<prefix.length() {
    read_aligned.write_char('-')
  }
  read_aligned.write_string(polish_reverse_string(read_mid.to_string()))
  for _k in 0..<suffix.length() {
    read_aligned.write_char('-')
  }
  let ref_aligned = StringBuilder()
  ref_aligned.write_string(prefix)
  ref_aligned.write_string(polish_reverse_string(ref_mid.to_string()))
  ref_aligned.write_string(suffix)
  (read_aligned.to_string(), ref_aligned.to_string())
}

///|
fn polish_anchor_read(
  index : Map[String, Array[Int32]],
  read : String,
  k : Int32,
) -> (String, Int32)? {
  let (forward_delta, forward_votes) = polish_count_offset(index, read, k)
  let rc = polish_revcomp(read)
  let (reverse_delta, reverse_votes) = polish_count_offset(index, rc, k)
  if forward_votes == 0 && reverse_votes == 0 {
    return None
  }
  let use_forward = if reverse_votes > forward_votes {
    false
  } else if forward_votes > reverse_votes {
    true
  } else {
    // 平票时优先正链, 确定性行为.
    true
  }
  if use_forward {
    match forward_delta {
      Some(delta) => Some((read, delta))
      None => None
    }
  } else {
    match reverse_delta {
      Some(delta) => Some((rc, delta))
      None => None
    }
  }
}

///|
fn polish_arr_add(array : Array[Int32], idx : Int32, amount : Int32) -> Unit {
  if idx >= 0 && idx < array.length() {
    array[idx] = array[idx] + amount
  }
}

///|
fn polish_map_add(map : Map[String, Int32], key : String) -> Unit {
  let current = match map.get(key) {
    Some(c) => c
    None => 0
  }
  map.set(key, current + 1)
}

///|
/// racon-lite 参考导向共识抛光.
///
/// 读段经 k-mer 锚定定位, 在参考窗口内与参考做全局比对后逐位点投票:
/// - 替换: 覆盖度 >= min_coverage 且目标碱基支持比例 >= min_fraction;
/// - 删除: gap 支持比例 >= min_fraction 且 gap 票数高于任何碱基;
/// - 插入: 插入串整体计数, 支持数 >= min_coverage 且占相邻列最大覆盖度
///   的比例 >= min_fraction.
///
/// 读段中的非 ACGT 碱基(如 N)不参与投票(避免被误计为删除票);
/// 参考位点为非 ACGT 时保守保留该碱基, 不做替换/删除.
pub fn racon_polish(
  reference : String,
  reads : Array[String],
  config~ : PolishConfig = {
    k: 15,
    min_coverage: 3,
    min_fraction: 0.6,
    flank: 50,
  },
) -> PolishResult raise AssemblyGraphError {
  let ref_len = reference.length()
  if ref_len == 0 {
    raise AssemblyGraphError("reference sequence is empty")
  }
  let k = if config.k > ref_len / 2 {
    let shrunk = ref_len / 2
    if shrunk < 4 {
      raise AssemblyGraphError("reference too short for k-mer anchoring")
    }
    shrunk
  } else {
    config.k
  }
  let index = polish_build_index(reference, k)
  let coverage : Array[Int32] = Array::make(ref_len, 0)
  let base_counts : Array[Array[Int32]] = Array(capacity=ref_len)
  for _i in 0..<ref_len {
    base_counts.push([0, 0, 0, 0])
  }
  let insert_slots : Array[Map[String, Int32]] = Array(capacity=ref_len + 1)
  for _i in 0..<=ref_len {
    insert_slots.push(Map([]))
  }
  let mut reads_used = 0
  for read in reads {
    if read.length() < k {
      continue
    }
    match polish_anchor_read(index, read, k) {
      None => ()
      Some((oriented, delta)) => {
        let ref_start = if delta - config.flank < 0 {
          0
        } else {
          delta - config.flank
        }
        let ref_end = if delta + oriented.length() + config.flank > ref_len {
          ref_len
        } else {
          delta + oriented.length() + config.flank
        }
        if ref_start >= ref_end {
          continue
        }
        let window = reference[ref_start:ref_end].to_owned()
        let (aligned_read, aligned_ref) = polish_semiglobal_align(
          oriented, window,
        )
        let alen = aligned_read.length()
        if alen != aligned_ref.length() {
          continue
        }
        reads_used = reads_used + 1
        let mut first_non_gap = -1
        let mut last_non_gap = -1
        for jj in 0..<alen {
          if aligned_read.unsafe_get(jj).unsafe_to_char() != '-' {
            if first_non_gap < 0 {
              first_non_gap = jj
            }
            last_non_gap = jj
          }
        }
        if first_non_gap < 0 {
          continue
        }
        let mut ref_consumed = 0
        let mut j = 0
        while j < alen {
          if aligned_ref.unsafe_get(j).unsafe_to_char() == '-' {
            let sb = StringBuilder()
            let mut j2 = j
            while j2 < alen &&
                  aligned_ref.unsafe_get(j2).unsafe_to_char() == '-' {
              let base = aligned_read.unsafe_get(j2).unsafe_to_char()
              if base != '-' {
                sb.write_char(base)
              }
              j2 = j2 + 1
            }
            let inserted = sb.to_string()
            let slot = ref_start + ref_consumed
            if inserted != "" &&
              j >= first_non_gap &&
              j <= last_non_gap &&
              slot >= 0 &&
              slot <= ref_len {
              polish_map_add(insert_slots[slot], inserted)
            }
            j = j2
          } else {
            let pos = ref_start + ref_consumed
            if pos >= 0 &&
              pos < ref_len &&
              j >= first_non_gap &&
              j <= last_non_gap {
              let base = aligned_read.unsafe_get(j).unsafe_to_char()
              if base == '-' {
                // 真实比对 gap: 计为删除票.
                coverage[pos] = coverage[pos] + 1
              } else {
                let idx = polish_base_index(base)
                if idx >= 0 {
                  coverage[pos] = coverage[pos] + 1
                  polish_arr_add(base_counts[pos], idx, 1)
                }
                // 非 ACGT 碱基(如 N)不投票, 避免被误计为删除票.
              }
            }
            ref_consumed = ref_consumed + 1
            j = j + 1
          }
        }
      }
    }
  }
  let deleted : Array[Bool] = Array::make(ref_len, false)
  let corrected : Array[String] = Array(capacity=ref_len)
  for i in 0..<ref_len {
    corrected.push(reference[i:i + 1].to_owned())
  }
  let edits : Array[PolishEdit] = Array(capacity=0)
  for pos in 0..<ref_len {
    if coverage[pos] < config.min_coverage {
      continue
    }
    let ref_char = reference.unsafe_get(pos).unsafe_to_char()
    let ref_idx = polish_base_index(ref_char)
    if ref_idx < 0 {
      // 参考碱基非 ACGT(如 N): 投票语义不明, 保守保留参考碱基.
      continue
    }
    let total = coverage[pos]
    let counts = base_counts[pos]
    let mut best_idx = 0
    let mut best_count = counts[0]
    for b in 1..<4 {
      if counts[b] > best_count {
        best_count = counts[b]
        best_idx = b
      }
    }
    let mut gap_count = total
    for b in 0..<4 {
      gap_count = gap_count - counts[b]
    }
    let gap_fraction = gap_count.to_double() / total.to_double()
    let best_fraction = best_count.to_double() / total.to_double()
    if gap_count > best_count && gap_fraction >= config.min_fraction {
      deleted[pos] = true
      corrected[pos] = ""
      edits.push({
        kind: "deletion",
        position: pos,
        ref_base: reference[pos:pos + 1].to_owned(),
        new_base: "",
        support: gap_count,
        coverage: total,
      })
    } else if best_count > 0 &&
      best_fraction >= config.min_fraction &&
      best_idx != ref_idx {
      let base = polish_base_char(best_idx).to_string()
      corrected[pos] = base
      edits.push({
        kind: "substitution",
        position: pos,
        ref_base: reference[pos:pos + 1].to_owned(),
        new_base: base,
        support: best_count,
        coverage: total,
      })
    }
  }
  let insertions : Array[String?] = Array::make(ref_len + 1, None)
  let insert_support : Array[Int32] = Array::make(ref_len + 1, 0)
  let insert_coverage : Array[Int32] = Array::make(ref_len + 1, 0)
  for slot in 0..<=ref_len {
    let slot_map = insert_slots[slot]
    let left_cov = if slot > 0 { coverage[slot - 1] } else { 0 }
    let right_cov = if slot < ref_len { coverage[slot] } else { 0 }
    let denom = if left_cov > right_cov { left_cov } else { right_cov }
    if denom < config.min_coverage {
      continue
    }
    let mut best_insert = ""
    let mut best_count = 0
    for candidate in slot_map.keys() {
      let count = match slot_map.get(candidate) {
        Some(c) => c
        None => 0
      }
      if count > best_count {
        best_count = count
        best_insert = candidate
      } else if count == best_count && candidate < best_insert {
        // 平票时按字典序取较小串, 保证确定性.
        best_insert = candidate
      }
    }
    if best_count >= config.min_coverage &&
      best_count.to_double() / denom.to_double() >= config.min_fraction {
      insertions[slot] = Some(best_insert)
      insert_support[slot] = best_count
      insert_coverage[slot] = denom
    }
  }
  let sb = StringBuilder()
  match insertions[0] {
    Some(ins) => {
      sb.write_string(ins)
      edits.push({
        kind: "insertion",
        position: 0,
        ref_base: "",
        new_base: ins,
        support: insert_support[0],
        coverage: insert_coverage[0],
      })
    }
    None => ()
  }
  for pos in 0..<ref_len {
    if !deleted[pos] {
      sb.write_string(corrected[pos])
    }
    match insertions[pos + 1] {
      Some(ins) => {
        sb.write_string(ins)
        edits.push({
          kind: "insertion",
          position: pos + 1,
          ref_base: "",
          new_base: ins,
          support: insert_support[pos + 1],
          coverage: insert_coverage[pos + 1],
        })
      }
      None => ()
    }
  }
  let mut total_cov = 0
  for c in coverage {
    total_cov = total_cov + c
  }
  let mean_coverage = total_cov.to_double() / ref_len.to_double()
  {
    reference,
    consensus: sb.to_string(),
    edits,
    reads_used,
    reads_total: reads.length(),
    mean_coverage,
  }
}

///|
/// 抛光结果 JSON 输出.
pub fn polish_result_to_json(result : PolishResult) -> String {
  let sb = StringBuilder()
  sb.write_string("{")
  sb.write_string("\"reference_length\":")
  sb.write_string(result.reference.length().to_string())
  sb.write_string(",\"consensus_length\":")
  sb.write_string(result.consensus.length().to_string())
  sb.write_string(",\"reads_total\":")
  sb.write_string(result.reads_total.to_string())
  sb.write_string(",\"reads_used\":")
  sb.write_string(result.reads_used.to_string())
  sb.write_string(",\"mean_coverage\":")
  sb.write_string(result.mean_coverage.to_string())
  sb.write_string(",\"edits\":[")
  for i, edit in result.edits {
    if i > 0 {
      sb.write_string(",")
    }
    sb.write_string("{\"kind\":\"")
    sb.write_string(edit.kind)
    sb.write_string("\",\"position\":")
    sb.write_string(edit.position.to_string())
    sb.write_string(",\"ref_base\":\"")
    sb.write_string(edit.ref_base)
    sb.write_string("\",\"new_base\":\"")
    sb.write_string(edit.new_base)
    sb.write_string("\",\"support\":")
    sb.write_string(edit.support.to_string())
    sb.write_string(",\"coverage\":")
    sb.write_string(edit.coverage.to_string())
    sb.write_string("}")
  }
  sb.write_string("]}")
  sb.to_string()
}
