// ============================================================================
// C2: Lifetime Analyses (Kaplan-Meier Survival Analysis)
// C3: Trajectory Clustering (AtomGroup 接口包装)
// C4: SASA 时间序列 (AtomGroup wrapper)
// ============================================================================

// ============================================================================
// C2: Kaplan-Meier Survival Estimator
// ============================================================================

///|
/// 从一组"生存时间"样本估计生存函数 S(t)（Kaplan-Meier 无参数估计）。
///
/// 输入：每个条目是 (duration, is_censored)。
///   - duration > 0：事件发生时间
///   - is_censored = false：事件在此时刻被观测到
///   - is_censored = true：到此时间为止事件未发生（右删失）
///
/// 返回 Array[(time, S(t))]。
pub fn mda_kaplan_meier(
  samples : Array[(Double, Bool)],
) -> Array[(Double, Double)] {
  let n = samples.length()
  if n == 0 {
    return Array(capacity=0)
  }
  // 按时间升序排序（sort_by 返回 Int）
  let data : Array[(Double, Bool)] = samples
  data.sort_by(mda_km_cmp)
  let out : Array[(Double, Double)] = Array(capacity=n)
  let mut cum_surv = 1.0
  let mut i = 0
  while i < n {
    // 合并同一时间点的事件
    let ti = data[i].0
    let mut d_i = 0 // 事件数（not censored）
    let mut c_i = 0 // 删失数
    let mut j = i
    while j < n && (data[j].0 - ti).abs() < 1.0e-10 {
      if data[j].1 {
        c_i = c_i + 1
      } else {
        d_i = d_i + 1
      }
      j = j + 1
    }
    let n_remaining = (n - i).to_double()
    if d_i > 0 && n_remaining > 0.0 {
      cum_surv = cum_surv * (1.0 - d_i.to_double() / n_remaining)
    }
    out.push((ti, cum_surv))
    i = j // 跳到下一个时间点
  }
  out
}

///|
fn mda_km_cmp(a : (Double, Bool), b : (Double, Bool)) -> Int32 {
  if a.0 < b.0 {
    -1
  } else if a.0 > b.0 {
    1
  } else {
    0
  }
}

///|
/// 从中断事件数组（frames × N 个连续段）估算 Mean Lifetime。
/// 输入：每段的持续帧数。
pub fn mda_mean_lifetime(segments : Array[Int32]) -> Double {
  let n = segments.length()
  if n == 0 {
    return 0.0
  }
  let mut total = 0
  let mut i = 0
  while i < n {
    total = total + segments[i]
    i = i + 1
  }
  total.to_double() / n.to_double()
}

// ============================================================================
// C3: AtomGroup 轨迹聚类 (RMSD + K-means / cutoff)
// ============================================================================

///|
fn mda_traj_to_dcd_format(
  traj : Array[Array[Vector3]],
) -> Array[Array[Array[Double]]] {
  let out : Array[Array[Array[Double]]] = Array(capacity=traj.length())
  let mut f = 0
  while f < traj.length() {
    let fr : Array[Array[Double]] = Array(capacity=traj[f].length())
    let mut i = 0
    while i < traj[f].length() {
      let a : Array[Double] = Array(capacity=3)
      a.push(traj[f][i].x)
      a.push(traj[f][i].y)
      a.push(traj[f][i].z)
      fr.push(a)
      i = i + 1
    }
    out.push(fr)
    f = f + 1
  }
  out
}

///|
/// 用 AtomGroup 的 selected atoms 做 K-means 轨迹聚类（RMSD 距离）。
/// 返回 `MdaTrajCluster`（底层已有类型：assignments/centroids/sizes/n_clusters）。
pub fn AtomGroup::trajectory_cluster_kmeans(
  self : AtomGroup,
  k~ : Int32 = 3,
) -> MdaTrajCluster {
  let nf = self.universe.n_frames()
  let idx = self.indices
  let n_sel = idx.length()
  let traj_xyz : Array[Array[Vector3]] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    let frame_vec : Array[Vector3] = Array(capacity=n_sel)
    let mut i = 0
    while i < n_sel {
      let id = idx[i]
      frame_vec.push(
        if id < fr.length() {
          fr[id]
        } else {
          Vector3::new(0.0, 0.0, 0.0)
        },
      )
      i = i + 1
    }
    traj_xyz.push(frame_vec)
    f = f + 1
  }
  let dcd_fmt = mda_traj_to_dcd_format(traj_xyz)
  mda_trajectory_kmeans(dcd_fmt, k)
}

///|
/// 用 AtomGroup 做 cutoff-based 轨迹聚类（RMSD 阈值，比 K-means 更自动）。
pub fn AtomGroup::trajectory_cluster_by_cutoff(
  self : AtomGroup,
  cutoff~ : Double = 2.0,
) -> MdaTrajCluster {
  let nf = self.universe.n_frames()
  let idx = self.indices
  let n_sel = idx.length()
  let traj_xyz : Array[Array[Vector3]] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    let frame_vec : Array[Vector3] = Array(capacity=n_sel)
    let mut i = 0
    while i < n_sel {
      let id = idx[i]
      frame_vec.push(
        if id < fr.length() {
          fr[id]
        } else {
          Vector3::new(0.0, 0.0, 0.0)
        },
      )
      i = i + 1
    }
    traj_xyz.push(frame_vec)
    f = f + 1
  }
  let dcd_fmt = mda_traj_to_dcd_format(traj_xyz)
  mda_trajectory_cluster_by_cutoff(dcd_fmt, cutoff)
}

///|
/// 从已聚类的轨迹中取第 k 个簇的"代表性帧"（RMSD 到簇质心最小的帧索引）。
pub fn AtomGroup::trajectory_representative_frame(
  self : AtomGroup,
  cluster_result : MdaTrajCluster,
  k : Int32,
) -> Int32 {
  // cluster_result.centroids[k] 已经是质心帧索引
  let cs = cluster_result.centroids()
  if k >= 0 && k < cs.length() {
    cs[k]
  } else {
    0
  }
}

///|
/// 从轨迹里提取一个代表性帧（返回新的 AtomGroup 或 Universe 子集不方便，
/// 这里直接返回帧的 frame index，调用者可以 frame() 拿到坐标）。
pub fn AtomGroup::trajectory_frame_rmsd_to(
  self : AtomGroup,
  frame_a : Int32,
  frame_b : Int32,
) -> Double {
  let fr_a = self.universe.frame(frame_a)
  let fr_b = self.universe.frame(frame_b)
  let idx = self.indices
  let n = idx.length()
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    let id = idx[i]
    if id < fr_a.length() && id < fr_b.length() {
      let dx = fr_a[id].x - fr_b[id].x
      let dy = fr_a[id].y - fr_b[id].y
      let dz = fr_a[id].z - fr_b[id].z
      sum = sum + dx * dx + dy * dy + dz * dz
    }
    i = i + 1
  }
  (sum / n.to_double()).sqrt()
}

// ============================================================================
// C4: AtomGroup::sasa_trajectory() — Shrake-Rupley 单帧的时间序列版
// ============================================================================

///|
/// SASA 总表面积时间序列。
/// 每帧独立算一次 Shrake-Rupley（与 sasa(frame=) 相同底层）。
pub fn AtomGroup::sasa_trajectory(
  self : AtomGroup,
  probe~ : Double = 1.4,
) -> Array[Double] {
  let nf = self.universe.n_frames()
  let out : Array[Double] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    out.push(self.sasa(frame=f, probe~))
    f = f + 1
  }
  out
}

///|
/// SASA 波动（per-atom 每帧方差，判断哪些残基表面积不稳定）。
pub fn AtomGroup::sasa_fluctuation(
  self : AtomGroup,
  probe~ : Double = 1.4,
) -> Array[Double] {
  let n = self.indices.length()
  let nf = self.universe.n_frames()
  if n == 0 || nf < 2 {
    let empty : Array[Double] = Array(capacity=n)
    let mut i = 0
    while i < n {
      empty.push(0.0)
      i = i + 1
    }
    return empty
  }
  let per_atom_traj : Array[Array[Double]] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    per_atom_traj.push(self.sasa_per_atom(frame=f, probe~))
    f = f + 1
  }
  // 每原子方差
  let out : Array[Double] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let mut sx = 0.0
    f = 0
    while f < nf {
      sx = sx + per_atom_traj[f][i]
      f = f + 1
    }
    let mean = sx / nf.to_double()
    let mut s2 = 0.0
    f = 0
    while f < nf {
      let d = per_atom_traj[f][i] - mean
      s2 = s2 + d * d
      f = f + 1
    }
    out.push((s2 / (nf - 1).to_double()).sqrt())
    i = i + 1
  }
  out
}
