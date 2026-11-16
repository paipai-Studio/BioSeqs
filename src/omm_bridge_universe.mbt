// ============================================================================
// B3: Universe ↔ OmmSystem 原生互转 + 模拟入口
// ============================================================================

///|
/// Universe → OmmSystem。原生方法，用户无需记忆 omm_bridge 命名。
pub fn Universe::to_omm_system(univ : Universe) -> OmmSystem {
  omm_system_from_universe(univ)
}

///|
/// 用 Universe 的当前帧初始化 OmmSimulation（跑模拟）。
/// OmmSimulation 内部使用 Leapfrog 积分器 + Langevin 恒温器。
/// 可传 dt / temperature / friction。
pub fn Universe::to_omm_simulation(
  univ : Universe,
  dt~ : Double = 0.001,
  temperature~ : Double = 300.0,
  friction~ : Double = 1.0,
) -> (OmmSimulation, OmmSystem) {
  let sys = omm_system_from_universe(univ)
  let nf = univ.n_frames()
  let start_pos_opt : Array[Array[Double]]? = if nf > 0 {
    let fr = univ.frame(0)
    let out : Array[Array[Double]] = Array(capacity=fr.length())
    let mut i = 0
    while i < fr.length() {
      out.push([fr[i].x, fr[i].y, fr[i].z])
      i = i + 1
    }
    Some(out)
  } else {
    None
  }
  let sim = OmmSimulation::new(
    system=sys,
    dt~,
    temperature~,
    friction~,
    positions=start_pos_opt,
  )
  (sim, sys)
}
