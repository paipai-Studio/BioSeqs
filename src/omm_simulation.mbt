/// omm_simulation.mbt — OpenMM `Simulation` / `State` application layer.
///
/// Mirrors `openmm.app.Simulation`: a facade over System + integrator +
/// reporters with the standard workflow
/// `minimize -> equilibrate -> production`, plus checkpoint serialization /
/// deserialization (`Simulation::saveState` / `loadState`).
///
/// The integrator is a Langevin (BAOAB) step built from `md_langevin_step` /
/// `md_langevin_finish`; `friction = 0` reduces it to plain velocity Verlet.
/// Units: kcal/mol, A, ps, K.

///|
/// Full simulation state snapshot (OpenMM `State`): positions, velocities,
/// box, time, step and energies.
pub(all) struct OmmState {
  mut positions : Array[Array[Double]]
  mut velocities : Array[Array[Double]]
  box : Array[Double]
  mut time : Double
  mut step : Int32
  mut potential_energy : Double
  mut kinetic_energy : Double
}

///|
pub fn OmmState::new(
  positions~ : Array[Array[Double]],
  velocities~ : Array[Array[Double]],
  box~ : Array[Double] = [0.0, 0.0, 0.0],
  time~ : Double = 0.0,
  step~ : Int32 = 0,
  potential_energy~ : Double = 0.0,
  kinetic_energy~ : Double = 0.0,
) -> OmmState {
  { positions, velocities, box, time, step, potential_energy, kinetic_energy, }
}

///|
/// Number of atoms covered by the state.
pub fn OmmState::n_atoms(self : OmmState) -> Int32 {
  self.positions.length()
}

///|
/// Total (potential + kinetic) energy.
pub fn OmmState::total_energy(self : OmmState) -> Double {
  self.potential_energy + self.kinetic_energy
}

///|
/// Deep copy of the state (all rows duplicated).
pub fn OmmState::copy(self : OmmState) -> OmmState {
  let n = self.positions.length()
  let pos : Array[Array[Double]] = md_new_forces(n)
  let vel : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    pos[i] = [self.positions[i][0], self.positions[i][1], self.positions[i][2]]
    vel[i] = [
      self.velocities[i][0],
      self.velocities[i][1],
      self.velocities[i][2],
    ]
  }
  {
    positions: pos,
    velocities: vel,
    box: [self.box[0], self.box[1], self.box[2]],
    time: self.time,
    step: self.step,
    potential_energy: self.potential_energy,
    kinetic_energy: self.kinetic_energy,
  }
}

///|
/// Kinetic energy from velocities: KE = 0.5 * sum m v^2.
pub fn omm_kinetic_energy(
  velocities : Array[Array[Double]],
  masses : Array[Double],
) -> Double {
  let mut ke = 0.0
  for i = 0; i < velocities.length(); i = i + 1 {
    let v2 = velocities[i][0] * velocities[i][0] +
      velocities[i][1] * velocities[i][1] +
      velocities[i][2] * velocities[i][2]
    ke = ke + 0.5 * masses[i] * v2
  }
  ke
}

///|
/// Temperature from kinetic energy: T = 2 KE / (N_dof * kB).
pub fn omm_temperature_from_energy(ke : Double, n_dof : Int32) -> Double {
  if n_dof <= 0 {
    return 0.0
  }
  2.0 * ke / (MD_KB * n_dof.to_double())
}

///|
/// Simulation facade (OpenMM `Simulation`): system + Langevin integrator +
/// state, with minimize / equilibrate / production workflow.
pub struct OmmSimulation {
  system : OmmSystem
  dt : Double
  temperature : Double
  friction : Double
  state : OmmState
}

///|
pub fn OmmSimulation::new(
  system~ : OmmSystem,
  dt~ : Double = 0.001,
  temperature~ : Double = 300.0,
  friction~ : Double = 1.0,
  positions~ : Array[Array[Double]]? = None,
) -> OmmSimulation {
  let n = system.num_particles()
  let pos : Array[Array[Double]] = md_new_forces(n)
  let vel : Array[Array[Double]] = md_new_forces(n)
  match positions {
    Some(p) =>
      for i = 0; i < n && i < p.length(); i = i + 1 {
        pos[i] = [p[i][0], p[i][1], p[i][2]]
      }
    None => ()
  }
  let st = OmmState::new(positions=pos, velocities=vel, box=[
    system.box[0],
    system.box[1],
    system.box[2],
  ])
  let sim = { system, dt, temperature, friction, state: st, }
  sim.refresh_energies()
  sim
}

///|
/// Access the underlying system.
pub fn OmmSimulation::system(self : OmmSimulation) -> OmmSystem {
  self.system
}

///|
/// Snapshot (deep copy) of the current state.
pub fn OmmSimulation::state(self : OmmSimulation) -> OmmState {
  self.state.copy()
}

///|
/// Replace positions (and reset cached energies).
pub fn OmmSimulation::set_positions(
  self : OmmSimulation,
  positions : Array[Array[Double]],
) -> Unit {
  let n = self.state.positions.length()
  for i = 0; i < n && i < positions.length(); i = i + 1 {
    self.state.positions[i] = [
      positions[i][0],
      positions[i][1],
      positions[i][2],
    ]
  }
  self.state.potential_energy = self.system.potential_energy(
    self.state.positions,
  )
}

///|
/// Replace velocities.
pub fn OmmSimulation::set_velocities(
  self : OmmSimulation,
  velocities : Array[Array[Double]],
) -> Unit {
  let n = self.state.velocities.length()
  for i = 0; i < n && i < velocities.length(); i = i + 1 {
    self.state.velocities[i] = [
      velocities[i][0],
      velocities[i][1],
      velocities[i][2],
    ]
  }
  self.state.kinetic_energy = omm_kinetic_energy(
    self.state.velocities,
    self.system.get_masses(),
  )
}

///|
/// Draw Maxwell-Boltzmann velocities at `temperature` (deterministic given
/// `seed`), then remove the center-of-mass momentum.
pub fn OmmSimulation::initialize_velocities(
  self : OmmSimulation,
  temperature~ : Double = 300.0,
  seed~ : Int64 = 42L,
) -> Unit {
  md_seed_rng(seed)
  let masses = self.system.get_masses()
  let n = self.state.velocities.length()
  let kT = MD_KB * temperature
  for i = 0; i < n; i = i + 1 {
    let s = (kT / masses[i]).sqrt()
    for d = 0; d < 3; d = d + 1 {
      self.state.velocities[i][d] = s * md_gaussian()
    }
  }
  // Remove COM momentum (3 translational DOF).
  let mut px = 0.0
  let mut py = 0.0
  let mut pz = 0.0
  let mut mtot = 0.0
  for i = 0; i < n; i = i + 1 {
    px = px + masses[i] * self.state.velocities[i][0]
    py = py + masses[i] * self.state.velocities[i][1]
    pz = pz + masses[i] * self.state.velocities[i][2]
    mtot = mtot + masses[i]
  }
  if mtot > 0.0 {
    for i = 0; i < n; i = i + 1 {
      self.state.velocities[i][0] = self.state.velocities[i][0] - px / mtot
      self.state.velocities[i][1] = self.state.velocities[i][1] - py / mtot
      self.state.velocities[i][2] = self.state.velocities[i][2] - pz / mtot
    }
  }
  self.state.kinetic_energy = omm_kinetic_energy(self.state.velocities, masses)
}

///|
/// Refresh cached energies (PE from the force field, KE from velocities).
pub fn OmmSimulation::refresh_energies(self : OmmSimulation) -> Unit {
  self.state.potential_energy = self.system.potential_energy(
    self.state.positions,
  )
  self.state.kinetic_energy = omm_kinetic_energy(
    self.state.velocities,
    self.system.get_masses(),
  )
}

///|
/// LocalEnergyMinimizer: minimize the potential energy.
/// `algorithm` is "lbfgs" (default) or "steepest_descent".
/// Returns the energy decrease E_before - E_after.
pub fn OmmSimulation::minimize(
  self : OmmSimulation,
  max_iter~ : Int32 = 200,
  tolerance~ : Double = 0.01,
  algorithm~ : String = "lbfgs",
) -> Double {
  let system = self.system
  fn force_fn(pos : Array[Array[Double]]) -> Array[Array[Double]] {
    system.forces(pos)
  }
  let e_before = self.system.potential_energy(self.state.positions)
  let minimized = if algorithm == "steepest_descent" {
    omm_minimize_steepest_descent(
      self.state.positions,
      force_fn,
      max_iter~,
      tolerance~,
    )
  } else {
    omm_minimize_lbfgs(self.state.positions, force_fn, max_iter~, tolerance~)
  }
  self.state.positions = minimized
  self.state.velocities = Array::make(minimized.length(), [0.0, 0.0, 0.0])
  self.refresh_energies()
  e_before - self.state.potential_energy
}

///|
/// Advance the simulation by `n_steps` Langevin (BAOAB) steps.
/// With friction = 0 this is plain velocity Verlet (NVE).
pub fn OmmSimulation::step(self : OmmSimulation, n_steps~ : Int32 = 1) -> Unit {
  let masses = self.system.get_masses()
  let system = self.system
  fn force_fn(pos : Array[Array[Double]]) -> Array[Array[Double]] {
    system.forces(pos)
  }
  for _s = 0; _s < n_steps; _s = _s + 1 {
    let forces = force_fn(self.state.positions)
    let (x_new, v_after_o) = md_langevin_step(
      self.state.positions,
      self.state.velocities,
      forces,
      masses,
      self.dt,
      self.temperature,
      self.friction,
    )
    let forces_new = force_fn(x_new)
    let v_new = md_langevin_finish(v_after_o, forces_new, masses, self.dt)
    self.state.positions = x_new
    self.state.velocities = v_new
    self.state.time = self.state.time + self.dt
    self.state.step = self.state.step + 1
    self.state.potential_energy = self.system.potential_energy(x_new)
    self.state.kinetic_energy = omm_kinetic_energy(v_new, masses)
  }
}

///|
/// Shared run loop for equilibrate/production: steps + optional state-data
/// reporting. `report_interval <= 0` disables reporting.
fn omm_simulation_run(
  sim : OmmSimulation,
  n_steps : Int32,
  report_interval : Int32,
) -> Array[OmmStateData] {
  let reporter = OmmStateReporter::new(
    interval=if report_interval > 0 { report_interval } else { 1 },
  )
  // Reporting counts steps relative to the start of this run so that
  // equilibrate/production phases report on their own cadence.
  let step0 = sim.state.step
  for _s = 0; _s < n_steps; _s = _s + 1 {
    sim.step(n_steps=1)
    if report_interval > 0 && (sim.state.step - step0) % report_interval == 0 {
      let n_dof = 3 * sim.state.positions.length() - 3
      reporter.report(
        sim.state.step,
        sim.state.time,
        sim.state.potential_energy,
        sim.state.kinetic_energy,
        omm_temperature_from_energy(sim.state.kinetic_energy, n_dof),
        volume=sim.state.box[0] * sim.state.box[1] * sim.state.box[2],
      )
    }
  }
  reporter.data()
}

///|
/// Equilibration phase: run steps with thermostat, optionally collecting
/// state data. Returns collected OmmStateData rows (empty when
/// `report_interval <= 0`).
pub fn OmmSimulation::equilibrate(
  self : OmmSimulation,
  n_steps~ : Int32,
  report_interval~ : Int32 = 0,
) -> Array[OmmStateData] {
  omm_simulation_run(self, n_steps, report_interval)
}

///|
/// Production phase: run steps collecting OmmStateData every
/// `report_interval` steps.
pub fn OmmSimulation::production(
  self : OmmSimulation,
  n_steps~ : Int32,
  report_interval~ : Int32 = 100,
) -> Array[OmmStateData] {
  omm_simulation_run(self, n_steps, report_interval)
}

///|
/// Current temperature from kinetic energy (3N-3 DOF).
pub fn OmmSimulation::temperature(self : OmmSimulation) -> Double {
  let n_dof = 3 * self.state.positions.length() - 3
  omm_temperature_from_energy(self.state.kinetic_energy, n_dof)
}

// ---------------------------------------------------------------------------
// Checkpoint serialization (Simulation::saveState / loadState)
// ---------------------------------------------------------------------------

///|
/// Serialize a state to a text checkpoint (positions + velocities + box +
/// time/step/energies), one atom per line. The text is platform independent
/// and round-trips through `omm_checkpoint_load`.
pub fn omm_checkpoint_save(state : OmmState) -> String {
  let sb : StringBuilder = StringBuilder()
  sb.write_string("OMMCHECKPOINT 1\n")
  sb.write_string("NATOMS " + state.positions.length().to_string() + "\n")
  sb.write_string(
    "BOX " +
    state.box[0].to_string() +
    " " +
    state.box[1].to_string() +
    " " +
    state.box[2].to_string() +
    "\n",
  )
  sb.write_string("TIME " + state.time.to_string() + "\n")
  sb.write_string("STEP " + state.step.to_string() + "\n")
  sb.write_string("PE " + state.potential_energy.to_string() + "\n")
  sb.write_string("KE " + state.kinetic_energy.to_string() + "\n")
  for i = 0; i < state.positions.length(); i = i + 1 {
    sb.write_string(
      state.positions[i][0].to_string() +
      " " +
      state.positions[i][1].to_string() +
      " " +
      state.positions[i][2].to_string() +
      " " +
      state.velocities[i][0].to_string() +
      " " +
      state.velocities[i][1].to_string() +
      " " +
      state.velocities[i][2].to_string() +
      "\n",
    )
  }
  sb.to_string()
}

///|
/// Parse one Double from the tail of a line ("KEY value").
fn omm_ckpt_double(rest : String) -> Double? {
  let parts = rest.split(" ").to_array()
  for p in parts {
    let t = p.to_owned().trim()
    if t.length() > 0 {
      return parse_double(t.to_owned())
    }
  }
  None
}

///|
/// Deserialize a text checkpoint produced by `omm_checkpoint_save`.
/// Returns Err with a diagnostic on malformed input.
pub fn omm_checkpoint_load(content : String) -> Result[OmmState, String] {
  let lines = content.split("\n").to_array()
  let mut n_atoms = -1
  let pos : Array[Array[Double]] = Array(capacity=0)
  let vel : Array[Array[Double]] = Array(capacity=0)
  let mut box_ = [0.0, 0.0, 0.0]
  let mut time = 0.0
  let mut step = 0
  let mut pe = 0.0
  let mut ke = 0.0
  let mut seen_header = false
  for line in lines {
    let l = line.to_owned().trim().to_owned()
    if l.length() == 0 {
      continue
    }
    if !seen_header {
      if !l.has_prefix("OMMCHECKPOINT") {
        return Err("checkpoint: missing OMMCHECKPOINT header")
      }
      seen_header = true
      continue
    }
    if l.has_prefix("NATOMS") {
      match omm_ckpt_double(l[6:].to_owned()) {
        Some(v) => n_atoms = v.to_int()
        None => return Err("checkpoint: bad NATOMS")
      }
      continue
    }
    if l.has_prefix("BOX") {
      let parts = l.split(" ").to_array()
      let vals : Array[Double] = Array(capacity=0)
      for p in parts {
        let t = p.to_owned().trim()
        if t.length() > 0 && t != "BOX" {
          match parse_double(t.to_owned()) {
            Some(v) => vals.push(v)
            None => ()
          }
        }
      }
      if vals.length() >= 3 {
        box_ = [vals[0], vals[1], vals[2]]
      }
      continue
    }
    if l.has_prefix("TIME") {
      match omm_ckpt_double(l[4:].to_owned()) {
        Some(v) => time = v
        None => ()
      }
      continue
    }
    if l.has_prefix("STEP") {
      match omm_ckpt_double(l[4:].to_owned()) {
        Some(v) => step = v.to_int()
        None => ()
      }
      continue
    }
    if l.has_prefix("PE") {
      match omm_ckpt_double(l[2:].to_owned()) {
        Some(v) => pe = v
        None => ()
      }
      continue
    }
    if l.has_prefix("KE") {
      match omm_ckpt_double(l[2:].to_owned()) {
        Some(v) => ke = v
        None => ()
      }
      continue
    }
    // Atom line: x y z vx vy vz
    let parts = l.split(" ").to_array()
    let vals : Array[Double] = Array(capacity=0)
    for p in parts {
      let t = p.to_owned().trim()
      if t.length() > 0 {
        match parse_double(t.to_owned()) {
          Some(v) => vals.push(v)
          None => ()
        }
      }
    }
    if vals.length() >= 6 {
      pos.push([vals[0], vals[1], vals[2]])
      vel.push([vals[3], vals[4], vals[5]])
    }
  }
  if n_atoms < 0 {
    return Err("checkpoint: missing NATOMS")
  }
  if pos.length() != n_atoms {
    return Err("checkpoint: atom count mismatch")
  }
  Ok(
    OmmState::new(
      positions=pos,
      velocities=vel,
      box=box_,
      time~,
      step~,
      potential_energy=pe,
      kinetic_energy=ke,
    ),
  )
}
