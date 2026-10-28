///|
/// mmseqs2 风格本地蛋白质搜索.
///
/// 流水线对齐 mmseqs2 的 prefilter→align 两段式:
/// (1) k-mer 倒排索引种子检索 + 对角线合并 + ungapped x-drop 延伸
///     (prefilter, 无空位快速评分);
/// (2) 候选靶标 Smith-Waterman 精化 + Karlin-Altschul E-value 排序.
/// 另提供 linclust 风格贪心预聚类(共享 k-mer 覆盖率聚类).
///
/// 与 mmseqs2 的差异(定位: 筛查级, 非其替身): k-mer 采用 5bit/aa
/// 精确编码(k<=12)而非 3bit reduced alphabet + 替换剖面; ungapped
/// 延伸单点对角线而非 profile; E-value 采用 BLOSUM62 ungapped
/// lambda/K(blosum62_stats), 线性 gap 近似 affine 11/1.
pub suberror MmseqsError {
  MmseqsError(String)
}

///|
/// 蛋白序列库.
pub struct MmseqsDB {
  ids : Array[String]
  seqs : Array[String]
}

///|
/// 新建空库.
pub fn MmseqsDB::new() -> MmseqsDB {
  { ids: [], seqs: [], }
}

///|
/// 添加序列(id 缺省为 "seq<i>"; 重复 id 抛 MmseqsError).
pub fn MmseqsDB::add(
  self : MmseqsDB,
  seq : String,
  id~ : String = "",
) -> MmseqsDB raise MmseqsError {
  let rid = if id.is_empty() {
    "seq" + self.seqs.length().to_string()
  } else {
    id
  }
  for existing in self.ids {
    if existing == rid {
      raise MmseqsError("duplicate sequence id: " + rid)
    }
  }
  self.ids.push(rid)
  self.seqs.push(seq)
  self
}

///|
/// 库大小(序列数).
pub fn MmseqsDB::size(self : MmseqsDB) -> Int32 {
  self.seqs.length()
}

///|
/// 序列(下标越界返回 None).
pub fn MmseqsDB::seq(self : MmseqsDB, i : Int32) -> String? {
  if i < 0 || i >= self.seqs.length() {
    None
  } else {
    Some(self.seqs[i])
  }
}

///|
/// 序列 id.
pub fn MmseqsDB::seq_id(self : MmseqsDB, i : Int32) -> String? {
  if i < 0 || i >= self.ids.length() {
    None
  } else {
    Some(self.ids[i])
  }
}

///|
/// 库总残基数(E-value 的 n 项).
pub fn MmseqsDB::total_residues(self : MmseqsDB) -> Int32 {
  let mut total = 0
  for s in self.seqs {
    total = total + s.length()
  }
  total
}

///|
/// k-mer 倒排索引(key → (序列下标, 位置) 列表).
pub struct MmseqsIndex {
  db : MmseqsDB
  k : Int32
  table : Map[Int64, Array[(Int32, Int32)]]
}

///|
/// 构建索引(k 缺省 6, mmseqs2 快速档蛋白默认; 范围 3..12 受 5bit/aa
/// 编码与 Int64 容量约束).
pub fn MmseqsIndex::build(
  db : MmseqsDB,
  k~ : Int32 = 6,
) -> MmseqsIndex raise MmseqsError {
  if k < 3 || k > 12 {
    raise MmseqsError("k must be in 3..12 (5bit/aa in Int64)")
  }
  let table : Map[Int64, Array[(Int32, Int32)]] = Map([], capacity=0)
  for si in 0..<db.seqs.length() {
    let s = db.seqs[si]
    let n = s.length()
    let mut pos = 0
    while pos + k <= n {
      match mmseqs_kmer_key(s, pos, k) {
        Some(key) =>
          match table.get(key) {
            Some(arr) => arr.push((si, pos))
            None => {
              let arr : Array[(Int32, Int32)] = Array(capacity=1)
              arr.push((si, pos))
              table[key] = arr
            }
          }
        None => ()
      }
      pos = pos + 1
    }
  }
  { db, k, table, }
}

///|
/// 索引的 k 值.
pub fn MmseqsIndex::k(self : MmseqsIndex) -> Int32 {
  self.k
}

///|
/// 索引条目数(不同 k-mer 数).
pub fn MmseqsIndex::num_kmers(self : MmseqsIndex) -> Int32 {
  self.table.size()
}

///|
/// prefilter 结果(单靶标最优对角线的 ungapped 分数).
pub struct MmseqsPrefilter {
  target : Int32
  diagonal : Int32
  ungapped_score : Double
  seeds : Int32
}

///|
/// prefilter 结果对应的库内序列下标.
pub fn MmseqsPrefilter::target(self : MmseqsPrefilter) -> Int32 {
  self.target
}

///|
/// 对角线(靶标位置 - 查询位置).
pub fn MmseqsPrefilter::diagonal(self : MmseqsPrefilter) -> Int32 {
  self.diagonal
}

///|
/// ungapped 延伸分数(BLOSUM62, x-drop 截断).
pub fn MmseqsPrefilter::ungapped_score(self : MmseqsPrefilter) -> Double {
  self.ungapped_score
}

///|
/// 该对角线命中的种子数.
pub fn MmseqsPrefilter::seeds(self : MmseqsPrefilter) -> Int32 {
  self.seeds
}

///|
/// 最终命中(SW 精化 + 统计).
pub struct MmseqsHit {
  target : Int32
  target_id : String
  raw_score : Double
  bitscore : Double
  evalue : Double
  identity : Double
  aligned_query : String
  aligned_target : String
  match_line : String
  query_range : (Int32, Int32)
  target_range : (Int32, Int32)
}

///|
/// 靶标序列下标.
pub fn MmseqsHit::target(self : MmseqsHit) -> Int32 {
  self.target
}

///|
/// 靶标序列 id.
pub fn MmseqsHit::target_id(self : MmseqsHit) -> String {
  self.target_id
}

///|
/// SW 原始分.
pub fn MmseqsHit::raw_score(self : MmseqsHit) -> Double {
  self.raw_score
}

///|
/// 位分(bit score).
pub fn MmseqsHit::bitscore(self : MmseqsHit) -> Double {
  self.bitscore
}

///|
/// E-value(Karlin-Altschul).
pub fn MmseqsHit::evalue(self : MmseqsHit) -> Double {
  self.evalue
}

///|
/// 一致度(0..1).
pub fn MmseqsHit::identity(self : MmseqsHit) -> Double {
  self.identity
}

///|
/// 比对结果行(查询/匹配线/靶标)格式化.
pub fn MmseqsHit::format(self : MmseqsHit) -> String {
  "  query  " +
  self.query_range.0.to_string() +
  "-" +
  self.query_range.1.to_string() +
  "  " +
  self.aligned_query +
  "\n          " +
  self.match_line +
  "\n  target " +
  self.target_range.0.to_string() +
  "-" +
  self.target_range.1.to_string() +
  "  " +
  self.aligned_target +
  "\n"
}

///|
/// prefilter: k-mer 种子 → 对角线分桶 → 相邻种子去重 → x-drop ungapped
/// 延伸, 返回每靶标最优对角线(分数降序, 同分按 target 升序).
/// xdrop 为延伸终止落差; min_score<=0 的对角线丢弃.
pub fn mmseqs_prefilter(
  index : MmseqsIndex,
  query : String,
  xdrop~ : Double = 25.0,
  min_score~ : Double = 30.0,
) -> Array[MmseqsPrefilter] {
  let k = index.k
  let matrix = @data.blosum62()
  // 1. 种子检索 → target → diag 两级分桶
  let by_target : Map[Int32, Map[Int32, Array[(Int32, Int32)]]] = Map(
    [],
    capacity=0,
  )
  let qn = query.length()
  let mut qpos = 0
  while qpos + k <= qn {
    match mmseqs_kmer_key(query, qpos, k) {
      Some(key) =>
        match index.table.get(key) {
          Some(posts) =>
            for p in posts {
              let dmap = match by_target.get(p.0) {
                Some(m) => m
                None => {
                  let m : Map[Int32, Array[(Int32, Int32)]] = Map(
                    [],
                    capacity=0,
                  )
                  by_target[p.0] = m
                  m
                }
              }
              let diag = p.1 - qpos
              match dmap.get(diag) {
                Some(arr) => arr.push((qpos, p.1))
                None => {
                  let arr : Array[(Int32, Int32)] = Array(capacity=1)
                  arr.push((qpos, p.1))
                  dmap[diag] = arr
                }
              }
            }
          None => ()
        }
      None => ()
    }
    qpos = qpos + 1
  }
  // 2. target 升序 → diag 升序, 保证确定性遍历
  let targets : Array[Int32] = []
  for entry in by_target {
    targets.push(entry.0)
  }
  mmseqs_sort_int(targets)
  let out : Array[MmseqsPrefilter] = []
  for tid in targets {
    let dmap = by_target[tid]
    let tseq = index.db.seqs[tid]
    let diags : Array[Int32] = []
    for entry in dmap {
      diags.push(entry.0)
    }
    mmseqs_sort_int(diags)
    // 3. 每对角线: 种子去重(相邻间隔 < k 合并) + 逐种子 x-drop 延伸取最优
    for diag in diags {
      let hits = dmap[diag]
      mmseqs_sort_pairs(hits)
      let mut best = 0.0
      let mut nseeds = 0
      let mut prev_q = -1000000
      for h in hits {
        if h.0 - prev_q < k {
          continue
        }
        prev_q = h.0
        nseeds = nseeds + 1
        let score = mmseqs_ungapped_extend(
          query,
          tseq,
          h.0,
          h.1,
          k,
          matrix,
          xdrop,
        )
        if score > best {
          best = score
        }
      }
      if nseeds > 0 && best >= min_score {
        out.push({
          target: tid,
          diagonal: diag,
          ungapped_score: best,
          seeds: nseeds,
        })
      }
    }
  }
  // 4. 分数降序, 同分 target 升序再 diagonal 升序
  mmseqs_sort_prefilters(out)
  out
}

///|
/// 完整搜索: prefilter 取前 max_targets 候选 → Smith-Waterman 精化 →
/// bit score/E-value → 阈值过滤 → E-value 升序输出.
/// E-value 采用 BLOSUM62 ungapped 统计(blosum62_stats), 线性 gap 近似.
pub fn mmseqs_search(
  index : MmseqsIndex,
  query : String,
  gap_penalty~ : Double = -10.0,
  evalue_threshold~ : Double = 1.0e-3,
  max_targets~ : Int32 = 50,
) -> Array[MmseqsHit] {
  let pre = mmseqs_prefilter(index, query)
  let stats = blosum62_stats()
  let m = query.length()
  let n = index.db.total_residues()
  let hits : Array[MmseqsHit] = []
  let limit = if max_targets < pre.length() {
    max_targets
  } else {
    pre.length()
  }
  for i in 0..<limit {
    let tid = pre[i].target
    let tseq = index.db.seqs[tid]
    let sw = smith_waterman_with_matrix(
      query,
      tseq,
      @data.blosum62(),
      gap_penalty,
    )
    let bs = bit_score(sw.score(), stats)
    let ev = evalue(sw.score(), m, n, stats)
    if ev <= evalue_threshold {
      hits.push({
        target: tid,
        target_id: index.db.ids[tid],
        raw_score: sw.score(),
        bitscore: bs,
        evalue: ev,
        identity: sw.percent_identity() / 100.0,
        aligned_query: sw.aligned1(),
        aligned_target: sw.aligned2(),
        match_line: sw.match_line(),
        query_range: sw.pos1(),
        target_range: sw.pos2(),
      })
    }
  }
  // E-value 升序, 同值按 target 升序
  mmseqs_sort_hits(hits)
  hits
}

///|
/// linclust 风格贪心预聚类: 序列按长度降序, 依次与已有簇代表比较
/// 共享 k-mer 覆盖率(共享 k-mer 数 / 新序列 k-mer 总数), 达阈值归簇,
/// 否则自成新簇. 返回每序列(输入序)的簇编号.
pub fn mmseqs_cluster(
  db : MmseqsDB,
  k~ : Int32 = 6,
  coverage~ : Double = 0.7,
) -> Array[Int32] raise MmseqsError {
  if k < 3 || k > 12 {
    raise MmseqsError("k must be in 3..12 (5bit/aa in Int64)")
  }
  if coverage <= 0.0 || coverage > 1.0 {
    raise MmseqsError("coverage must be in (0, 1]")
  }
  let n = db.seqs.length()
  // 长度降序(稳定: 同长按下标升序)
  let order : Array[Int32] = Array(capacity=n)
  for i in 0..<n {
    order.push(i)
  }
  let mut i = 0
  while i < order.length() {
    let mut best = i
    let mut j = i + 1
    while j < order.length() {
      let lj = db.seqs[order[j]].length()
      let lb = db.seqs[order[best]].length()
      if lj > lb || (lj == lb && order[j] < order[best]) {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = order[i]
      order[i] = order[best]
      order[best] = t
    }
    i = i + 1
  }
  // 贪心聚类
  let assign : Array[Int32] = Array::make(n, -1)
  let reps : Array[Int32] = []
  let rep_kmers : Array[Map[Int64, Bool]] = []
  for idx in order {
    let s = db.seqs[idx]
    let keys : Array[Int64] = []
    let mut pos = 0
    while pos + k <= s.length() {
      match mmseqs_kmer_key(s, pos, k) {
        Some(key) => keys.push(key)
        None => ()
      }
      pos = pos + 1
    }
    let mut joined = false
    for ci in 0..<reps.length() {
      let rep_set = rep_kmers[ci]
      let mut shared = 0
      for key in keys {
        if rep_set.contains(key) {
          shared = shared + 1
        }
      }
      let cov = if keys.length() > 0 {
        shared.to_double() / keys.length().to_double()
      } else {
        0.0
      }
      if cov >= coverage {
        assign[idx] = ci
        joined = true
        break
      }
    }
    if !joined {
      let set : Map[Int64, Bool] = Map([], capacity=0)
      for key in keys {
        set[key] = true
      }
      rep_kmers.push(set)
      reps.push(idx)
      assign[idx] = reps.length() - 1
    }
  }
  assign
}

///|
/// 标准氨基酸(单字母大写)→ 3bit 外编码 0..19(BLAST 顺序
/// A R N D C Q E G H I L K M F P S T W Y V); 非标准字符 None.
fn mmseqs_aa_code(code : Int32) -> Int64? {
  let c = if code >= 97 && code <= 122 { code - 32 } else { code }
  if c == 65 {
    Some(0L)
  } else if c == 82 {
    Some(1L)
  } else if c == 78 {
    Some(2L)
  } else if c == 68 {
    Some(3L)
  } else if c == 67 {
    Some(4L)
  } else if c == 81 {
    Some(5L)
  } else if c == 69 {
    Some(6L)
  } else if c == 71 {
    Some(7L)
  } else if c == 72 {
    Some(8L)
  } else if c == 73 {
    Some(9L)
  } else if c == 76 {
    Some(10L)
  } else if c == 75 {
    Some(11L)
  } else if c == 77 {
    Some(12L)
  } else if c == 70 {
    Some(13L)
  } else if c == 80 {
    Some(14L)
  } else if c == 83 {
    Some(15L)
  } else if c == 84 {
    Some(16L)
  } else if c == 87 {
    Some(17L)
  } else if c == 89 {
    Some(18L)
  } else if c == 86 {
    Some(19L)
  } else {
    None
  }
}

///|
/// 序列 [pos, pos+k) 的 5bit/aa k-mer key(含非标字符返回 None).
fn mmseqs_kmer_key(seq : String, pos : Int32, k : Int32) -> Int64? {
  let mut key = 0L
  for i in 0..<k {
    match mmseqs_aa_code(seq.unsafe_get(pos + i).to_int()) {
      None => return None
      Some(v) => key = (key << 5) | v
    }
  }
  Some(key)
}

///|
/// k-mer key 计算(测试/诊断辅助: 暴露编码语义供黑盒验证).
pub fn mmseqs_kmer_key_of(seq : String, pos : Int32, k : Int32) -> Int64? {
  mmseqs_kmer_key(seq, pos, k)
}

///|
/// ungapped x-drop 双向延伸: 从种子 (q0, t0) 起左右延伸累计 BLOSUM62
/// 分数, 各向最大值; 分数跌至该向最优 - xdrop 以下终止.
/// 返回 best_left + seed + best_right.
fn mmseqs_ungapped_extend(
  query : String,
  tseq : String,
  q0 : Int32,
  t0 : Int32,
  k : Int32,
  matrix : @data.SubstitutionMatrix,
  xdrop : Double,
) -> Double {
  // 种子段分数(全同段, 恒为对角线 k 步)
  let mut seed = 0.0
  for i in 0..<k {
    seed = seed + mmseqs_sub(query, q0 + i, tseq, t0 + i, matrix)
  }
  // 右延伸
  let mut right = 0.0
  let mut r_best = 0.0
  let mut qi = q0 + k
  let mut ti = t0 + k
  while qi < query.length() && ti < tseq.length() {
    right = right + mmseqs_sub(query, qi, tseq, ti, matrix)
    if right > r_best {
      r_best = right
    }
    if right < r_best - xdrop {
      break
    }
    qi = qi + 1
    ti = ti + 1
  }
  // 左延伸
  let mut left = 0.0
  let mut l_best = 0.0
  let mut qj = q0 - 1
  let mut tj = t0 - 1
  while qj >= 0 && tj >= 0 {
    left = left + mmseqs_sub(query, qj, tseq, tj, matrix)
    if left > l_best {
      l_best = left
    }
    if left < l_best - xdrop {
      break
    }
    qj = qj - 1
    tj = tj - 1
  }
  l_best + seed + r_best
}

///|
/// 单残基替换分.
fn mmseqs_sub(
  q : String,
  qi : Int32,
  t : String,
  ti : Int32,
  matrix : @data.SubstitutionMatrix,
) -> Double {
  let c1 = q.unsafe_get(qi).unsafe_to_char().to_string()
  let c2 = t.unsafe_get(ti).unsafe_to_char().to_string()
  matrix.get_score_case_insensitive(c1, c2).to_double()
}

///|
/// Int32 数组升序排序(选择排序, 规模小).
fn mmseqs_sort_int(arr : Array[Int32]) -> Unit {
  let mut i = 0
  while i < arr.length() {
    let mut best = i
    let mut j = i + 1
    while j < arr.length() {
      if arr[j] < arr[best] {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = arr[i]
      arr[i] = arr[best]
      arr[best] = t
    }
    i = i + 1
  }
}

///|
/// (qpos, tpos) 对按 qpos 升序排序.
fn mmseqs_sort_pairs(arr : Array[(Int32, Int32)]) -> Unit {
  let mut i = 0
  while i < arr.length() {
    let mut best = i
    let mut j = i + 1
    while j < arr.length() {
      if arr[j].0 < arr[best].0 {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = arr[i]
      arr[i] = arr[best]
      arr[best] = t
    }
    i = i + 1
  }
}

///|
/// prefilter 排序: 分数降序, 同分 target 升序再 diagonal 升序.
fn mmseqs_sort_prefilters(arr : Array[MmseqsPrefilter]) -> Unit {
  let mut i = 0
  while i < arr.length() {
    let mut best = i
    let mut j = i + 1
    while j < arr.length() {
      let a = arr[j]
      let b = arr[best]
      if a.ungapped_score > b.ungapped_score ||
        (
          a.ungapped_score == b.ungapped_score &&
          (
            a.target < b.target ||
            (a.target == b.target && a.diagonal < b.diagonal)
          )
        ) {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = arr[i]
      arr[i] = arr[best]
      arr[best] = t
    }
    i = i + 1
  }
}

///|
/// hits 排序: E-value 升序, 同值按 target 升序.
fn mmseqs_sort_hits(arr : Array[MmseqsHit]) -> Unit {
  let mut i = 0
  while i < arr.length() {
    let mut best = i
    let mut j = i + 1
    while j < arr.length() {
      let a = arr[j]
      let b = arr[best]
      if a.evalue < b.evalue || (a.evalue == b.evalue && a.target < b.target) {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = arr[i]
      arr[i] = arr[best]
      arr[best] = t
    }
    i = i + 1
  }
}
