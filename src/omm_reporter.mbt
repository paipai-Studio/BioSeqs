/// Reporter framework for MD simulation output.
///
/// OpenMM equivalent: `DCDReporter`, `StateDataReporter`, `CheckpointReporter`
///
/// Provides hooks to collect trajectory, energy, and state data
/// during a simulation at regular intervals.

///|
/// Simulation state data collected per step
pub struct OmmStateData {
  step : Int32
  time : Double
  potential_energy : Double
  kinetic_energy : Double
  total_energy : Double
  temperature : Double
  volume : Double
}

///|
pub fn OmmStateData::step(self : OmmStateData) -> Int32 {
  self.step
}

///|
pub fn OmmStateData::time(self : OmmStateData) -> Double {
  self.time
}

///|
pub fn OmmStateData::total_energy(self : OmmStateData) -> Double {
  self.total_energy
}

///|
pub fn OmmStateData::temperature(self : OmmStateData) -> Double {
  self.temperature
}

///|
/// A reporter that collects state data during simulation
pub struct OmmStateReporter {
  report_interval : Int32
  data : Array[OmmStateData]
}

///|
pub fn OmmStateReporter::new(interval~ : Int32 = 100) -> OmmStateReporter {
  { report_interval: interval, data: [], }
}

///|
pub fn OmmStateReporter::should_report(
  self : OmmStateReporter,
  step : Int32,
) -> Bool {
  step % self.report_interval == 0
}

///|
pub fn OmmStateReporter::report(
  self : OmmStateReporter,
  step : Int32,
  time : Double,
  potential : Double,
  kinetic : Double,
  temperature : Double,
  volume~ : Double = 0.0,
) -> Unit {
  self.data.push({
    step,
    time,
    potential_energy: potential,
    kinetic_energy: kinetic,
    total_energy: potential + kinetic,
    temperature,
    volume,
  })
}

///|
pub fn OmmStateReporter::data(self : OmmStateReporter) -> Array[OmmStateData] {
  self.data
}

///|
pub fn OmmStateReporter::n_reports(self : OmmStateReporter) -> Int32 {
  self.data.length()
}

///|
/// DCD trajectory reporter
pub struct OmmDCDReporter {
  report_interval : Int32
  frames : Array[Array[Array[Double]]] // [frame][atom][xyz]
}

///|
pub fn OmmDCDReporter::new(interval~ : Int32 = 100) -> OmmDCDReporter {
  { report_interval: interval, frames: [], }
}

///|
pub fn OmmDCDReporter::should_report(
  self : OmmDCDReporter,
  step : Int32,
) -> Bool {
  step % self.report_interval == 0
}

///|
pub fn OmmDCDReporter::report(
  self : OmmDCDReporter,
  positions : Array[Array[Double]],
) -> Unit {
  let n = positions.length()
  let frame : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    frame[i] = [positions[i][0], positions[i][1], positions[i][2]]
  }
  self.frames.push(frame)
}

///|
pub fn OmmDCDReporter::frames(
  self : OmmDCDReporter,
) -> Array[Array[Array[Double]]] {
  self.frames
}

///|
pub fn OmmDCDReporter::n_frames(self : OmmDCDReporter) -> Int32 {
  self.frames.length()
}

///|
/// Checkpoint reporter (stores full state: positions + velocities)
pub struct OmmCheckpointReporter {
  report_interval : Int32
  checkpoints : Array[(Array[Array[Double]], Array[Array[Double]])]
}

///|
pub fn OmmCheckpointReporter::new(
  interval~ : Int32 = 1000,
) -> OmmCheckpointReporter {
  { report_interval: interval, checkpoints: [], }
}

///|
pub fn OmmCheckpointReporter::should_report(
  self : OmmCheckpointReporter,
  step : Int32,
) -> Bool {
  step % self.report_interval == 0
}

///|
pub fn OmmCheckpointReporter::report(
  self : OmmCheckpointReporter,
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
) -> Unit {
  let n = positions.length()
  let pos_copy : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let vel_copy : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos_copy[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel_copy[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  self.checkpoints.push((pos_copy, vel_copy))
}

///|
pub fn OmmCheckpointReporter::last_checkpoint(
  self : OmmCheckpointReporter,
) -> (Array[Array[Double]], Array[Array[Double]])? {
  if self.checkpoints.length() == 0 {
    None
  } else {
    Some(self.checkpoints[self.checkpoints.length() - 1])
  }
}

///|
pub fn omm_format_state_table(data : Array[OmmStateData]) -> String {
  let mut result = "Step    Time(ps)    PE(kJ/mol)    KE(kJ/mol)    Total(kJ/mol)    Temp(K)\n"
  for d in data {
    result = result +
      d.step.to_string() +
      "    " +
      d.time.to_string() +
      "    " +
      d.potential_energy.to_string() +
      "    " +
      d.kinetic_energy.to_string() +
      "    " +
      d.total_energy.to_string() +
      "    " +
      d.temperature.to_string() +
      "\n"
  }
  result
}

// ============================================================================
// O1: CheckpointReporter 完整版本（带 box/time/step + 序列化）
// ============================================================================

///|
/// 完整快照：step, time, box, positions, velocities。
pub struct OmmCheckpoint {
  step : Int32
  time : Double
  box : Array[Double]
  positions : Array[Array[Double]]
  velocities : Array[Array[Double]]
}

///|
/// 带时间/步长/box 的完整 checkpoint reporter。
/// 支持 serialize 成文本（每行一行原子：x,y,z,vx,vy,vz）方便落盘。
pub struct OmmFullCheckpointReporter {
  report_interval : Int32
  checkpoints : Array[OmmCheckpoint]
}

///|
pub fn OmmFullCheckpointReporter::new(
  interval~ : Int32 = 1000,
) -> OmmFullCheckpointReporter {
  { report_interval: interval, checkpoints: [], }
}

///|
pub fn OmmFullCheckpointReporter::should_report(
  self : OmmFullCheckpointReporter,
  step : Int32,
) -> Bool {
  step % self.report_interval == 0
}

///|
pub fn OmmFullCheckpointReporter::report(
  self : OmmFullCheckpointReporter,
  step : Int32,
  time : Double,
  box : Array[Double],
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
) -> Unit {
  let n = positions.length()
  let pos_copy : Array[Array[Double]] = Array(capacity=n)
  let vel_copy : Array[Array[Double]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    pos_copy.push([positions[i][0], positions[i][1], positions[i][2]])
    vel_copy.push([velocities[i][0], velocities[i][1], velocities[i][2]])
    i = i + 1
  }
  let box_copy : Array[Double] = [box[0], box[1], box[2]]
  self.checkpoints.push({
    step,
    time,
    box: box_copy,
    positions: pos_copy,
    velocities: vel_copy,
  })
}

///|
pub fn OmmFullCheckpointReporter::last(
  self : OmmFullCheckpointReporter,
) -> OmmCheckpoint? {
  if self.checkpoints.length() > 0 {
    Some(self.checkpoints[self.checkpoints.length() - 1])
  } else {
    None
  }
}

///|
/// 把最后一个 checkpoint 应用到 simulation 的 state（positions/velocities/time/step）。
/// box 字段不是 OmmState.mut，也不能整体替换 OmmSimulation.state，
/// 所以 box 更新需要通过 omm_simulation_from_checkpoint 重新构造。
pub fn OmmFullCheckpointReporter::restore(
  self : OmmFullCheckpointReporter,
  sim : OmmSimulation,
) -> Unit {
  match self.last() {
    Some(cp) => {
      sim.set_positions(cp.positions)
      sim.set_velocities(cp.velocities)
      // time 和 step 也不是 mut，用户需要 omm_simulation_from_checkpoint 完整恢复
      ignore(cp.time)
      ignore(cp.step)
    }
    None => ()
  }
}

///|
/// 把最后一个 checkpoint 序列化成文本（复用 omm_checkpoint_save）。
pub fn OmmFullCheckpointReporter::serialize_last(
  self : OmmFullCheckpointReporter,
) -> String {
  match self.last() {
    Some(cp) => {
      let st = OmmState::new(
        positions=cp.positions,
        velocities=cp.velocities,
        box=cp.box,
        time=cp.time,
        step=cp.step,
      )
      omm_checkpoint_save(st)
    }
    None => ""
  }
}

///|
/// 从已有 checkpoint 文本恢复（复用 omm_checkpoint_load）。
pub fn omm_checkpoint_restore(text : String) -> OmmFullCheckpointReporter {
  let rep = OmmFullCheckpointReporter::new()
  match omm_checkpoint_load(text) {
    Ok(st) => rep.report(st.step, st.time, st.box, st.positions, st.velocities)
    Err(_) => ()
  }
  rep
}
