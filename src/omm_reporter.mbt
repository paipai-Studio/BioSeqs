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
/// Format state data as a table (for logging)
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
