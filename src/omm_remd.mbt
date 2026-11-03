/// Replica Exchange Molecular Dynamics (REMD).
///
/// OpenMM equivalent: `ReplicaExchange` (via OpenMMTools)
///
/// Runs multiple replicas at different temperatures in parallel,
/// periodically swapping configurations between adjacent temperatures
/// to overcome energy barriers.

///|
/// A single replica in REMD
pub struct OmmReplica {
  mut temperature : Double
  positions : Array[Array[Double]]
  velocities : Array[Array[Double]]
  energy : Double
}

///|
pub fn OmmReplica::temperature(self : OmmReplica) -> Double {
  self.temperature
}

///|
pub fn OmmReplica::energy(self : OmmReplica) -> Double {
  self.energy
}

///|
/// Create replicas at geometrically spaced temperatures.
///
/// n_replicas: number of replicas
/// t_min: lowest temperature
/// t_max: highest temperature
/// initial_positions: starting positions (same for all replicas)
pub fn omm_create_replicas(
  n_replicas : Int32,
  t_min : Double,
  t_max : Double,
  initial_positions : Array[Array[Double]],
) -> Array[OmmReplica] {
  if n_replicas <= 0 {
    return []
  }
  let replicas : Array[OmmReplica] = Array(capacity=0)
  for i = 0; i < n_replicas; i = i + 1 {
    let temp = if n_replicas == 1 {
      t_min
    } else {
      let ratio = t_max / t_min
      t_min * @math.pow(ratio, i.to_double() / (n_replicas - 1).to_double())
    }
    let n = initial_positions.length()
    let pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
    let vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
    for j = 0; j < n; j = j + 1 {
      pos[j] = [
        initial_positions[j][0],
        initial_positions[j][1],
        initial_positions[j][2],
      ]
    }
    replicas.push({
      temperature: temp,
      positions: pos,
      velocities: vel,
      energy: 0.0,
    })
  }
  replicas
}

///|
/// Attempt a swap between two adjacent replicas.
/// Metropolis criterion: accept if exp((1/kT_i - 1/kT_j) * (E_j - E_i)) > rand
///
/// Returns true if swap was accepted.
pub fn omm_attempt_replica_swap(
  replica_i : OmmReplica,
  replica_j : OmmReplica,
) -> Bool {
  let kB = 0.001987204
  let beta_i = 1.0 / (kB * replica_i.temperature)
  let beta_j = 1.0 / (kB * replica_j.temperature)
  let dE = replica_j.energy - replica_i.energy
  let dBeta = beta_i - beta_j
  let delta = dBeta * dE
  if delta < 0.0 {
    true
  } else {
    md_rand_double() < @math.exp(-delta)
  }
}

///|
/// Perform one round of replica exchanges.
/// Attempts swaps between all adjacent pairs (odd-even decomposition).
///
/// Returns the swap pattern (array of bools, true = swapped).
pub fn omm_replica_exchange_round(
  replicas : Array[OmmReplica],
  parity~ : Int32 = 0, // 0 = even pairs, 1 = odd pairs
) -> Array[Bool] {
  let n = replicas.length()
  let swaps : Array[Bool] = Array::make(n - 1, false)
  let start = parity % 2
  let mut i = start
  while i < n - 1 {
    if omm_attempt_replica_swap(replicas[i], replicas[i + 1]) {
      swaps[i] = true
      // Swap temperatures (configurations stay, temperatures swap)
      let temp_i = replicas[i].temperature
      replicas[i].temperature = replicas[i + 1].temperature
      replicas[i + 1].temperature = temp_i
    }
    i = i + 2
  }
  swaps
}

///|
/// Run REMD for n_steps.
///
/// replicas: initial replicas
/// force_fn: force calculation
/// masses: per-atom masses
/// dt: time step
/// n_steps: total MD steps
/// swap_interval: steps between swap attempts
/// step_fn: function to advance one replica by one MD step
///
/// Returns final replicas and swap acceptance statistics.
pub fn omm_run_remd(
  replicas : Array[OmmReplica],
  n_steps : Int32,
  swap_interval~ : Int32 = 100,
  step_fn : (OmmReplica, Double) -> OmmReplica, // (replica, dt) -> new replica
) -> (Array[OmmReplica], Int32, Int32) {
  let mut total_swaps = 0
  let mut total_attempts = 0
  for step = 0; step < n_steps; step = step + 1 {
    // Advance each replica
    for i = 0; i < replicas.length(); i = i + 1 {
      replicas[i] = step_fn(replicas[i], 1.0) // dt=1.0 placeholder
    }
    // Attempt swaps
    if step > 0 && step % swap_interval == 0 {
      let parity = step / swap_interval % 2
      let swaps = omm_replica_exchange_round(replicas, parity~)
      for s in swaps {
        total_attempts = total_attempts + 1
        if s {
          total_swaps = total_swaps + 1
        }
      }
    }
  }
  (replicas, total_swaps, total_attempts)
}

///|
/// Compute swap acceptance rate.
pub fn omm_swap_acceptance_rate(n_swaps : Int32, n_attempts : Int32) -> Double {
  if n_attempts == 0 {
    return 0.0
  }
  n_swaps.to_double() / n_attempts.to_double()
}

///|
/// Optimal temperature spacing for REMD.
/// Uses the formula from Kone & Kofke (2005):
/// T_{i+1}/T_i = 1 + sqrt(1/(C * N_dof)) where C ~ 1 for acceptance ~ 30%
pub fn omm_optimal_temperatures(
  n_replicas : Int32,
  t_min : Double,
  n_dof : Int32,
) -> Array[Double] {
  if n_replicas <= 0 {
    return []
  }
  let temps : Array[Double] = Array::make(n_replicas, 0.0)
  let ratio = 1.0 + (1.0 / n_dof.to_double()).sqrt()
  temps[0] = t_min
  for i = 1; i < n_replicas; i = i + 1 {
    temps[i] = temps[i - 1] * ratio
  }
  temps
}
