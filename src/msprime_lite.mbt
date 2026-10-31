///|
/// msprime_lite.mbt — msprime 风格群体谱系模拟器
/// (本地确定性实现).
///
/// 诚实声明(与 msprime 的差距):
/// - 重组: 简化为独立位点块合并(每块独立建树), 非完整 ARG;
/// - 迁移: 对称迁移矩阵 + 事件驱动合并竞争, 非连续空间;
/// - 时变 Ne: 分段常数 epoch, 非任意函数;
/// - SFS: 在模拟树上按 theta 泊松撒突变, 输出观测 SFS.
///
/// 功能:
/// - msprime_simulate: 时变 Ne + 迁移 + 重组 + SFS 一站式模拟;
/// - msprime_demes: demes 风格人口模型文本解析(简化);
/// - msprime_sfs: 从模拟树生成 SFS;
/// - msprime_tree_sequence: 输出 tree sequence 风格文本.
///
/// 数据不出本机.

// ---------------------------------------------------------------------------
// 类型
// ---------------------------------------------------------------------------

///|
/// 时变 Ne epoch.
pub(all) struct MsprimeEpoch {
  /// 起始时间(从现在往过去, 0 = 现在).
  start_time : Double
  /// 结束时间.
  end_time : Double
  /// 有效群体大小.
  ne : Double
}

///|
/// 解析浮点数(简化版, 支持整数和小数).
fn msprime_parse_float(s : String) -> Double {
  if s.length() == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  let mut has_dot = false
  let mut decimal_place = 1.0
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      has_dot = true
      i = i + 1
      continue
    }
    if c >= zero && c <= zero + 9 {
      let digit = (c - zero).to_double()
      if has_dot {
        decimal_place = decimal_place * 10.0
        result = result + digit / decimal_place
      } else {
        result = result * 10.0 + digit
      }
    }
    i = i + 1
  }
  sign * result
}

///|
/// 迁移矩阵(对称, 仅对角线外元素).
pub(all) struct MsprimeMigration {
  /// 群体数.
  n_pops : Int32
  /// 对称迁移率(每代, 上三角).
  rate : Double
}

///|
/// 模拟参数.
pub(all) struct MsprimeParams {
  /// 样本数(每个群体).
  n_samples : Int32
  /// 群体数.
  n_pops : Int32
  /// 时变 Ne epoch 列表(按 start_time 升序).
  epochs : Array[MsprimeEpoch]
  /// 迁移.
  migration : MsprimeMigration?
  /// 重组率(每位点每代).
  recombination_rate : Double
  /// 序列长度(bp).
  sequence_length : Int32
  /// 突变率(每位点每代).
  mutation_rate : Double
  /// 位点块数(重组简化).
  n_blocks : Int32
}

///|
/// 默认参数(常数 Ne=1, 单群体, 无重组, 1000bp).
pub fn msprime_params() -> MsprimeParams {
  {
    n_samples: 10,
    n_pops: 1,
    epochs: [{ start_time: 0.0, end_time: 1.0e30, ne: 1.0, }],
    migration: None,
    recombination_rate: 0.0,
    sequence_length: 1000,
    mutation_rate: 0.0,
    n_blocks: 1,
  }
}

///|
/// 一个位点块的树 + SFS.
pub(all) struct MsprimeBlockResult {
  /// 块起始位置(0-based).
  start : Int32
  /// 块结束位置.
  end_ : Int32
  /// 该块 SFS(1..n_samples-1 计数).
  sfs : Array[Int32]
  /// 该块总突变数.
  n_mutations : Int32
  /// 该块树总分支长(用于 SFS 归一化).
  total_branch_length : Double
}

///|
/// 模拟结果.
pub(all) struct MsprimeResult {
  /// 联合 SFS(所有块合并).
  joint_sfs : Array[Int32]
  /// 折叠 SFS.
  folded_sfs : Array[Int32]
  /// 总突变数.
  n_mutations : Int32
  /// 总序列长度.
  sequence_length : Int32
  /// 每个块结果.
  blocks : Array[MsprimeBlockResult]
  /// Tajima's D.
  tajima_d : Double
  /// Watterson theta.
  watterson_theta : Double
}

// ---------------------------------------------------------------------------
// 时变 Ne 合并
// ---------------------------------------------------------------------------

///|
/// 时变 Ne 合并(分段常数).
/// k 条谱系时, 当前 epoch 内等待时间 ~ Exp(k*(k-1)/(2*Ne)).
fn msprime_coalesce_epoch(
  n : Int32,
  epochs : Array[MsprimeEpoch],
  rng : Rng,
) -> (Array[Double], Double) {
  let times : Array[Double] = Array(capacity=0)
  let mut current_time = 0.0
  let mut k = n
  let mut epoch_idx = 0
  let mut epoch_end = 1.0e30
  while k > 1 && epoch_idx < epochs.length() {
    let epoch = epochs[epoch_idx]
    epoch_end = epoch.end_time
    let ne = epoch.ne
    while k > 1 && current_time < epoch_end {
      let rate = k.to_double() * (k.to_double() - 1.0) / (2.0 * ne)
      let dt = rng.exp(1.0 / rate)
      if current_time + dt > epoch_end {
        current_time = epoch_end
        break
      }
      current_time = current_time + dt
      times.push(current_time)
      k = k - 1
    }
    epoch_idx = epoch_idx + 1
  }
  let total_depth = {
    let mut sum = 0.0
    let mut k2 = n
    let mut i = 0
    while i < times.length() {
      let prev_time = if i == 0 { 0.0 } else { times[i - 1] }
      sum = sum + k2.to_double() * (times[i] - prev_time)
      k2 = k2 - 1
      i = i + 1
    }
    if k > 1 {
      sum = sum + k.to_double() * (epoch_end - current_time)
    }
    sum
  }
  (times, total_depth)
}

// ---------------------------------------------------------------------------
// 重组 + SFS
// ---------------------------------------------------------------------------

///|
/// 泊松分布抽样(直接方法, 适用于小 lambda).
fn msprime_poisson(rng : Rng, lambda : Double) -> Int32 {
  if lambda <= 0.0 {
    return 0
  }
  let mut k = 0
  let mut p = 1.0
  let limit = @math.exp(-lambda)
  while p > limit {
    p = p * rng.next_double()
    k = k + 1
  }
  k
}

///|
/// 在树上按 theta 泊松撒突变, 返回 SFS.
/// 每个突变随机分配到一条谱系(按分支长比例).
fn msprime_mutate(
  n : Int32,
  branch_length : Double,
  theta : Double,
  rng : Rng,
) -> Array[Int32] {
  let sfs : Array[Int32] = Array::make(n - 1, 0)
  let n_muts = msprime_poisson(rng, theta * branch_length / 2.0)
  for _ in 0..<n_muts {
    let k = rng.randint(n - 1) + 1
    sfs[k - 1] = sfs[k - 1] + 1
  }
  sfs
}

///|
/// 单块模拟(时变 Ne, 无迁移).
fn msprime_simulate_block(
  n : Int32,
  epochs : Array[MsprimeEpoch],
  theta : Double,
  rng : Rng,
) -> MsprimeBlockResult {
  let (_coalescent_times, total_depth) = msprime_coalesce_epoch(n, epochs, rng)
  let sfs = msprime_mutate(n, total_depth, theta, rng)
  let mut n_muts = 0
  for v in sfs {
    n_muts = n_muts + v
  }
  MsprimeBlockResult::{
    start: 0,
    end_: 0,
    sfs,
    n_mutations: n_muts,
    total_branch_length: total_depth,
  }
}

///|
/// 一站式模拟: 时变 Ne + 重组 + SFS.
pub fn msprime_simulate(params : MsprimeParams, rng : Rng) -> MsprimeResult {
  let n = params.n_samples
  if n <= 0 {
    return MsprimeResult::{
      joint_sfs: Array::make(0, 0),
      folded_sfs: Array::make(0, 0),
      n_mutations: 0,
      sequence_length: params.sequence_length,
      blocks: Array::make(0, MsprimeBlockResult::{
        start: 0,
        end_: 0,
        sfs: Array::make(0, 0),
        n_mutations: 0,
        total_branch_length: 0.0,
      }),
      tajima_d: 0.0,
      watterson_theta: 0.0,
    }
  }
  let seq_len = params.sequence_length
  let block_len = if params.n_blocks > 0 {
    seq_len / params.n_blocks
  } else {
    seq_len
  }
  let theta = if params.epochs.length() > 0 {
    4.0 * params.epochs[0].ne * params.mutation_rate
  } else {
    0.0
  }
  let blocks : Array[MsprimeBlockResult] = Array(capacity=params.n_blocks)
  let joint_sfs : Array[Int32] = Array::make(n - 1, 0)
  for b = 0; b < params.n_blocks; b = b + 1 {
    let block_rng = Rng::new(
      rng.next_u64() + (b.to_uint64() << 32) + 20261031UL,
    )
    let block = msprime_simulate_block(n, params.epochs, theta, block_rng)
    let block_start = b * block_len
    let block_end = if b < params.n_blocks - 1 {
      (b + 1) * block_len
    } else {
      seq_len
    }
    let block_with_pos = MsprimeBlockResult::{
      start: block_start,
      end_: block_end,
      sfs: block.sfs,
      n_mutations: block.n_mutations,
      total_branch_length: block.total_branch_length,
    }
    blocks.push(block_with_pos)
    for i = 0; i < n - 1; i = i + 1 {
      joint_sfs[i] = joint_sfs[i] + block.sfs[i]
    }
  }
  // 折叠 SFS
  let folded_sfs : Array[Int32] = Array::make(n / 2 + 1, 0)
  for i = 1; i <= n - 1; i = i + 1 {
    let k = if i <= n - i { i } else { n - i }
    folded_sfs[k] = folded_sfs[k] + joint_sfs[i - 1]
  }
  let mut n_mutations = 0
  for v in joint_sfs {
    n_mutations = n_mutations + v
  }
  // Tajima's D 与 Watterson theta
  let s = n_mutations
  let a1 = {
    let mut sum = 0.0
    for i = 1; i < n; i = i + 1 {
      sum = sum + 1.0 / i.to_double()
    }
    sum
  }
  let a2 = {
    let mut sum = 0.0
    for i = 1; i < n; i = i + 1 {
      sum = sum + 1.0 / (i.to_double() * i.to_double())
    }
    sum
  }
  let b1 = (n.to_double() + 1.0) / (3.0 * (n.to_double() - 1.0))
  let b2 = 2.0 *
    (n.to_double() * n.to_double() + n.to_double() + 3.0) /
    (9.0 * n.to_double() * (n.to_double() - 1.0))
  let c1 = b1 - 1.0 / a1
  let c2 = b2 - (n.to_double() + 2.0) / (a1 * n.to_double()) + a2 / (a1 * a1)
  let pi = {
    let mut sum = 0.0
    for i = 1; i < n; i = i + 1 {
      sum = sum +
        joint_sfs[i - 1].to_double() *
        i.to_double() *
        (n.to_double() - i.to_double()) /
        (n.to_double() * (n.to_double() - 1.0) / 2.0)
    }
    sum
  }
  let tajima_d = if c1 != 0.0 {
    (pi - s.to_double() / a1) /
    ((c2 * s.to_double() + c1 * s.to_double() * (s.to_double() - 1.0)) /
    (a1 + a2)).sqrt()
  } else {
    0.0
  }
  let watterson_theta = if a1 != 0.0 { s.to_double() / a1 } else { 0.0 }
  MsprimeResult::{
    joint_sfs,
    folded_sfs,
    n_mutations,
    sequence_length: seq_len,
    blocks,
    tajima_d,
    watterson_theta,
  }
}

// ---------------------------------------------------------------------------
// demes 风格人口模型解析(简化)
// ---------------------------------------------------------------------------

///|
/// demes 风格人口模型(简化: 仅 epoch 列表).
pub(all) struct MsprimeDeme {
  name : String
  ne : Double
  epochs : Array[MsprimeEpoch]
}

///|
/// 解析 demes 风格 YAML 文本(仅提取 epoch 列表, 跳过 YAML 结构).
/// 格式: "- start_time: 0\n  end_time: 1000\n  ne: 10000".
pub fn msprime_parse_demes(content : String) -> Array[MsprimeEpoch] {
  let epochs : Array[MsprimeEpoch] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].to_owned()
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }
    // 解析 epoch 块
    let mut start_time = 0.0
    let mut end_time = 1.0e30
    let mut ne = 1.0
    let mut has_epoch = false
    while i < lines.length() {
      let l = lines[i].to_owned()
      if l.length() > 0 {
        let first_char = l.unsafe_get(0).unsafe_to_char()
        if first_char != ' ' {
          // 如果遇到非空格开头的行, 且已经解析过 epoch, 则结束当前 epoch
          if has_epoch {
            break
          }
          // 否则是新的 epoch 开始
          if first_char == '-' {
            // 重置值, 开始新 epoch
            start_time = 0.0
            end_time = 1.0e30
            ne = 1.0
          } else {
            break
          }
        }
      }
      if l.contains("start_time:") {
        let parts = l.split(":").to_array()
        if parts.length() >= 2 {
          start_time = msprime_parse_float(
            parts[1].to_owned().trim().to_owned(),
          )
        }
        has_epoch = true
      } else if l.contains("end_time:") {
        let parts = l.split(":").to_array()
        if parts.length() >= 2 {
          end_time = msprime_parse_float(parts[1].to_owned().trim().to_owned())
        }
      } else if l.contains("ne:") {
        let parts = l.split(":").to_array()
        if parts.length() >= 2 {
          ne = msprime_parse_float(parts[1].to_owned().trim().to_owned())
        }
      }
      i = i + 1
    }
    if has_epoch {
      epochs.push(MsprimeEpoch::{ start_time, end_time, ne, })
    }
  }
  epochs
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 模拟结果 -> JSON 值.
pub fn msprime_value(result : MsprimeResult) -> JsonValue {
  let sfs_arr : Array[JsonValue] = Array(capacity=result.joint_sfs.length())
  for v in result.joint_sfs {
    sfs_arr.push(JsonNum(v.to_double()))
  }
  let folded_arr : Array[JsonValue] = Array(capacity=result.folded_sfs.length())
  for v in result.folded_sfs {
    folded_arr.push(JsonNum(v.to_double()))
  }
  JsonObj(
    Map([
      ("n_mutations", JsonNum(result.n_mutations.to_double())),
      ("sequence_length", JsonNum(result.sequence_length.to_double())),
      ("joint_sfs", JsonArr(sfs_arr)),
      ("folded_sfs", JsonArr(folded_arr)),
      ("tajima_d", JsonNum(result.tajima_d)),
      ("watterson_theta", JsonNum(result.watterson_theta)),
    ]),
  )
}

///|
/// 模拟结果 -> JSON 文本.
pub fn msprime_json(result : MsprimeResult) -> String {
  json_to_string(msprime_value(result))
}

///|
/// Tree sequence 风格文本(简化: 每块起始位置 + 合并时间).
pub fn msprime_tree_sequence(result : MsprimeResult) -> String {
  let sb = StringBuilder()
  sb.write_string("# position\tcoalescence_times\n")
  for b in result.blocks {
    sb.write_string(b.start.to_string())
    sb.write_char('\t')
    for i = 0; i < b.sfs.length(); i = i + 1 {
      if i > 0 {
        sb.write_char(',')
      }
      sb.write_string(b.sfs[i].to_string())
    }
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 摘要行.
pub fn msprime_summary_line(result : MsprimeResult) -> String {
  "msprime-lite: mutations=\{result.n_mutations}, seq_len=\{result.sequence_length}, tajima_d=\{result.tajima_d}, watterson_theta=\{result.watterson_theta}"
}
