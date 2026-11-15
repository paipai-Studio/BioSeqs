///|
/// PERMANOVA Power 分析（Monte Carlo bootstrap 功效估计）。
///
/// 在给定 pilot 距离矩阵和分组的前提下，估计 PERMANOVA 检验在不同样本量
/// 下能以 α 水平检测到真实效应的统计功效。实现参考 scikit-bio 0.7.3
/// `skbio.stats.power.permanova_power`：
///
/// 1. PCoA 对 pilot 距离矩阵降维，得到样本在欧式 PCoA 空间的坐标；
/// 2. 估计各组在 PCoA 各轴上的均值 μ_{g,j} 以及组内平均方差 σ²；
/// 3. 对每个候选样本量 n，用独立同方差近似在 PCoA 空间采样：
///    每个样本在每个轴上独立从 N(μ, σ) 采样（不需要 Cholesky）；
/// 4. 从采样坐标重建 pairwise Euclidean 距离矩阵，跑 PERMANOVA；
/// 5. Monte Carlo 重复 B 次，Power = #(p < α) / B。
///
/// 这种做法是 PERMANOVA Power 的标准非参数近似：假设 PCoA 前 k 个轴
/// 独立且同方差（scikit-bio 默认 k = rank of distance matrix）。
///

///|
/// Power 曲线的一个点。
pub struct SkbioPowerPoint {
  sample_size : Int32
  power : Double
}

// ---------------------------------------------------------------------------
// 内部工具
// ---------------------------------------------------------------------------

///|
/// Box–Muller 采一个标准正态样本。
fn skbp_box_muller(seed : Int32) -> (Int32, Double) {
  let lcg = fn(s : Int32) -> Int32 { (s * 1103515245 + 12345) & 0x7FFFFFFF }
  let mut s = seed
  s = lcg(s)
  let u1 = (s.to_double() + 1.0) / 2147483647.0
  s = lcg(s)
  let u2 = (s.to_double() + 1.0) / 2147483647.0
  let two_pi = 6.283185307179586
  let r = (-2.0 * @math.ln(u1)).sqrt()
  let theta = two_pi * u2
  let z = r * @math.cos(theta)
  (s, z)
}

///|
/// 计算各 group 的样本 index 列表。
fn skbp_group_indices(grouping : Array[Int32]) -> Array[Array[Int32]] {
  let levels : Array[Int32] = Array(capacity=0)
  for g in grouping {
    let mut dup = false
    for lv in levels {
      if lv == g { dup = true; break }
    }
    if !dup { levels.push(g) }
  }
  let out : Array[Array[Int32]] = Array(capacity=0)
  for lv in levels {
    let members : Array[Int32] = Array(capacity=0)
    for i = 0; i < grouping.length(); i = i + 1 {
      if grouping[i] == lv { members.push(i) }
    }
    out.push(members)
  }
  out
}

///|
/// PCoA 坐标 → pairwise Euclidean distance matrix。
fn skbp_coords_to_dm(coords : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = coords.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      let k = coords[0].length()
      let mut sum_sq = 0.0
      let a = coords[i]
      let b = coords[j]
      for ax = 0; ax < k; ax = ax + 1 {
        let d = a[ax] - b[ax]
        sum_sq = sum_sq + d * d
      }
      row.push(sum_sq.sqrt())
    }
    out.push(row)
  }
  out
}

// ---------------------------------------------------------------------------
// PCoA 空间的 pilot 估计
// ---------------------------------------------------------------------------

struct PilotEst {
  group_means : Array[Array[Double]]
  sigma2 : Double
  k : Int32
  group_labels : Array[Int32]
  group_sizes : Array[Int]
}

fn skbp_estimate_pilot(
  samples : Array[Array[Double]],
  grouping : Array[Int32],
) -> PilotEst {
  let n = samples.length()
  let k = if n > 0 { samples[0].length() } else { 0 }
  let groups = skbp_group_indices(grouping)
  let ng = groups.length()
  let group_means : Array[Array[Double]] = Array(capacity=ng)
  let group_sizes : Array[Int32] = Array(capacity=ng)
  let group_labels : Array[Int32] = Array(capacity=ng)
  let mut total_ss = 0.0
  let mut total_count = 0
  for gi = 0; gi < ng; gi = gi + 1 {
    let members = groups[gi]
    let m = members.length()
    group_sizes.push(m)
    group_labels.push(grouping[members[0]])
    let gmean : Array[Double] = Array(capacity=k)
    for ax = 0; ax < k; ax = ax + 1 {
      let mut sum_a = 0.0
      for idx in members { sum_a = sum_a + samples[idx][ax] }
      gmean.push(sum_a / m.to_double())
    }
    group_means.push(gmean)
    for idx in members {
      for ax = 0; ax < k; ax = ax + 1 {
        let d = samples[idx][ax] - gmean[ax]
        total_ss = total_ss + d * d
      }
      total_count = total_count + k
    }
  }
  let sigma2 = if total_count > 0 { total_ss / total_count.to_double() } else { 0.0 }
  PilotEst::{ group_means, sigma2, k, group_labels, group_sizes }
}

// ---------------------------------------------------------------------------
// PERMANOVA Power 主入口
// ---------------------------------------------------------------------------

///|
/// PERMANOVA Power 曲线：给定 pilot 距离矩阵 + grouping，估计各候选样本量
/// 下的功效。
///
/// - `dm`: n×n pilot 距离矩阵
/// - `grouping`: n 个样本的组标签（Int32 编码）
/// - `candidate_ns`: 候选总样本量
/// - `alpha`: 显著性水平（默认 0.05）
/// - `permutations`: 每次 PERMANOVA 用多少次 permutation（默认 999）
/// - `n_bootstrap`: Monte Carlo 次数（默认 1000）
/// - `seed`: RNG 种子
pub fn skbio_permanova_power(
  dm : Array[Array[Double]],
  grouping : Array[Int32],
  candidate_ns : Array[Int32],
  alpha~ : Double = 0.05,
  permutations~ : Int32 = 999,
  n_bootstrap~ : Int32 = 1000,
  seed~ : Int32 = 42,
) -> Array[SkbioPowerPoint] {
  let out : Array[SkbioPowerPoint] = Array(capacity=candidate_ns.length())
  let n = dm.length()
  if n == 0 || grouping.length() != n || candidate_ns.length() == 0 { return out }
  let pcoa = skbio_pcoa(dm)
  let samples = pcoa.samples
  if samples.length() == 0 { return out }
  let pilot = skbp_estimate_pilot(samples, grouping)
  let k = pilot.k
  let sigma = pilot.sigma2.sqrt()
  let ng = pilot.group_labels.length()
  if ng < 2 { return out }

  let mut state = seed
  for cand_n in candidate_ns {
    let mut total_pilot_sum : Int = 0
    for g in pilot.group_sizes { total_pilot_sum = total_pilot_sum + g }
    let new_group_sizes : Array[Int] = Array(capacity=ng)
    let mut assigned : Int = 0
    for gi = 0; gi < ng; gi = gi + 1 {
      if gi == ng - 1 {
        new_group_sizes.push(cand_n.to_int64().to_int() - assigned)
      } else {
        let size_d = (cand_n.to_double() * pilot.group_sizes[gi].to_double()) / total_pilot_sum.to_double()
        new_group_sizes.push(size_d.floor().to_int64().to_int())
        assigned = assigned + new_group_sizes[gi]
      }
    }

    let mut count_reject = 0
    for _b = 0; _b < n_bootstrap; _b = _b + 1 {
      let new_coords : Array[Array[Double]] = []
      let new_grouping : Array[Int32] = []
      for gi = 0; gi < ng; gi = gi + 1 {
        let size_g = new_group_sizes[gi]
        for _s = 0; _s < size_g; _s = _s + 1 {
          let row : Array[Double] = Array(capacity=k)
          for ax = 0; ax < k; ax = ax + 1 {
            let (ns, z) = skbp_box_muller(state)
            state = ns
            row.push(pilot.group_means[gi][ax] + sigma * z)
          }
          new_coords.push(row)
          new_grouping.push(pilot.group_labels[gi])
        }
      }
      let new_dm = skbp_coords_to_dm(new_coords)
      let res = skbio_permanova(new_dm, new_grouping, permutations, state)
      state = (state * 1103515245 + 12345) & 0x7FFFFFFF
      if res.p_value < alpha {
        count_reject = count_reject + 1
      }
    }
    let power = count_reject.to_double() / n_bootstrap.to_double()
    out.push(SkbioPowerPoint::{ sample_size: cand_n, power })
  }
  out
}
