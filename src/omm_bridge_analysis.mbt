// ============================================================================
// B1 + A7: OpenMM 模拟 ↔ MDAnalysis 分析一键桥接
// ============================================================================

///|
/// 模拟轨迹采样：positions × energies × (温度/压力可选)。
pub struct OmmTrajectorySample {
  frames : Array[Array[Vector3]]
  pe : Array[Double]
  ke : Array[Double]
  temperature : Array[Double]
  pressure : Array[Double]
  density : Array[Double]
  bond : Array[Double]
  angle : Array[Double]
  dihedral : Array[Double]
  lj : Array[Double]
  coulomb : Array[Double]
  time : Array[Double]
  steps : Array[Int32]
}

fn md_boltzmann() -> Double { 8.314462618e-3 }

fn md_safe_double(d : Double) -> Double {
  if d.is_nan() { 0.0 } else { d }
}

///|
/// 从 OmmSimulation 运行 n_steps 步，每 report_interval 步收集一帧。
/// 返回采样结果（可转换为 Universe 喂 MDAnalysis 全套分析）。
pub fn omm_simulation_run_to_universe(
  sim : OmmSimulation,
  system : OmmSystem,
  n_steps~ : Int32 = 1000,
  report_interval~ : Int32 = 10,
) -> OmmTrajectorySample {
  let n_frames_expected : Int32 = n_steps / report_interval
  let s = OmmTrajectorySample::{
    frames : Array(capacity=n_frames_expected),
    pe : Array(capacity=n_frames_expected),
    ke : Array(capacity=n_frames_expected),
    temperature : Array(capacity=n_frames_expected),
    pressure : Array(capacity=n_frames_expected),
    density : Array(capacity=n_frames_expected),
    bond : Array(capacity=n_frames_expected),
    angle : Array(capacity=n_frames_expected),
    dihedral : Array(capacity=n_frames_expected),
    lj : Array(capacity=n_frames_expected),
    coulomb : Array(capacity=n_frames_expected),
    time : Array(capacity=n_frames_expected),
    steps : Array(capacity=n_frames_expected),
  }
  let mut i : Int32 = 0
  while i < n_steps {
    sim.step(n_steps=report_interval)
    let st = sim.state()
    let n_atoms = st.positions.length()
    let frame_vecs : Array[Vector3] = Array(capacity=n_atoms)
    let mut k = 0
    while k < n_atoms {
      frame_vecs.push(Vector3::new(st.positions[k][0], st.positions[k][1], st.positions[k][2]))
      k = k + 1
    }
    s.frames.push(frame_vecs)
    let (bond_e, ang_e, dih_e, lj_e, coul_e) = system.energy_decomposed(st.positions)
    let pe = st.potential_energy
    let ke = st.kinetic_energy
    s.pe.push(pe); s.ke.push(ke)
    s.bond.push(bond_e); s.angle.push(ang_e); s.dihedral.push(dih_e)
    s.lj.push(lj_e); s.coulomb.push(coul_e)
    let natoms_d = n_atoms.to_double()
    let temp = if natoms_d > 0.0 { md_safe_double((2.0 / 3.0) * ke / (natoms_d * md_boltzmann())) } else { 0.0 }
    s.temperature.push(temp)
    let box = st.box
    let vol = box[0] * box[1] * box[2]
    let pressure = if vol > 0.0 && natoms_d > 0.0 {
      md_safe_double((natoms_d * md_boltzmann() * temp / vol) * 602.214)
    } else { 0.0 }
    s.pressure.push(pressure)
    let density = if vol > 0.0 { md_safe_double(natoms_d / vol) } else { 0.0 }
    s.density.push(density)
    s.time.push(st.time)
    s.steps.push(st.step)
    i = i + report_interval
  }
  s
}

// ============================================================================
// OmmTrajectorySample 方法（每个方法单独 fn Type::method 声明）
// ============================================================================

pub fn OmmTrajectorySample::n_frames(self : OmmTrajectorySample) -> Int32 { self.frames.length() }
pub fn OmmTrajectorySample::total_potential_energy(self : OmmTrajectorySample) -> Array[Double] { self.pe }
pub fn OmmTrajectorySample::kinetic_energy(self : OmmTrajectorySample) -> Array[Double] { self.ke }
pub fn OmmTrajectorySample::bond(self : OmmTrajectorySample) -> Array[Double] { self.bond }
pub fn OmmTrajectorySample::angle(self : OmmTrajectorySample) -> Array[Double] { self.angle }
pub fn OmmTrajectorySample::dihedral(self : OmmTrajectorySample) -> Array[Double] { self.dihedral }
pub fn OmmTrajectorySample::lj(self : OmmTrajectorySample) -> Array[Double] { self.lj }
pub fn OmmTrajectorySample::coulomb(self : OmmTrajectorySample) -> Array[Double] { self.coulomb }
pub fn OmmTrajectorySample::temperature_trajectory(self : OmmTrajectorySample) -> Array[Double] { self.temperature }
pub fn OmmTrajectorySample::pressure_trajectory(self : OmmTrajectorySample) -> Array[Double] { self.pressure }
pub fn OmmTrajectorySample::density_trajectory(self : OmmTrajectorySample) -> Array[Double] { self.density }
pub fn OmmTrajectorySample::time(self : OmmTrajectorySample) -> Array[Double] { self.time }

pub fn OmmTrajectorySample::total_energy(self : OmmTrajectorySample) -> Array[Double] {
  let out : Array[Double] = Array(capacity=self.pe.length())
  let mut i = 0
  while i < self.pe.length() { out.push(self.pe[i] + self.ke[i]); i = i + 1 }
  out
}

///|
/// 能量统计：(mean, stddev, min, max)。
pub fn OmmTrajectorySample::energy_stats(self : OmmTrajectorySample, e : Array[Double]) -> (Double, Double, Double, Double) {
  let n = e.length()
  if n == 0 { return (0.0, 0.0, 0.0, 0.0) }
  let mut sum = 0.0
  let mut mn = 1.0e30
  let mut mx = -1.0e30
  let mut i = 0
  while i < n {
    sum = sum + e[i]
    if e[i] < mn { mn = e[i] }
    if e[i] > mx { mx = e[i] }
    i = i + 1
  }
  let mean = sum / n.to_double()
  let mut s2 = 0.0
  i = 0
  while i < n { let d = e[i] - mean; s2 = s2 + d * d; i = i + 1 }
  let sd = if n > 1 { (s2 / (n - 1).to_double()).sqrt() } else { 0.0 }
  (mean, sd, mn, mx)
}

///|
/// 把采样的轨迹坐标 → 可直接喂 mda AtomGroup 的 Universe。
pub fn OmmTrajectorySample::to_universe(self : OmmTrajectorySample, reference : Universe) -> Universe {
  let frames : Array[Array[Array[Double]]] = Array(capacity=self.frames.length())
  let mut f = 0
  while f < self.frames.length() {
    let n = self.frames[f].length()
    let row : Array[Array[Double]] = Array(capacity=n)
    let mut i = 0
    while i < n {
      let v = self.frames[f][i]
      row.push([v.x, v.y, v.z])
      i = i + 1
    }
    frames.push(row); f = f + 1
  }
  Universe::from_dcd_frames(reference.atoms, frames)
}
