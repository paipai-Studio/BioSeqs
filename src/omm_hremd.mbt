/// omm_hremd.mbt — Hamiltonian replica exchange (HREMD).
///
/// OpenMM equivalent: `ForceField`-swapped replica exchange beyond pure
/// temperature scaling (lambda/Hamiltonian replicas, mixed CV states).
///
/// Each replica k has its own inverse temperature beta_k and Hamiltonian.
/// Swapping the configurations x_i, x_j between replicas i and j follows
/// the Metropolis criterion
///
///   Delta = beta_i [U_i(x_j) - U_i(x_i)] + beta_j [U_j(x_i) - U_j(x_j)]
///   P(accept) = min(1, exp(-Delta))
///
/// which reduces to the temperature-REMD criterion when all replicas share
/// one Hamiltonian. Units: kcal/mol (energies), 1/(kcal/mol) (betas).

///|
/// One attempted replica-pair swap.
pub struct OmmHREMDSwap {
  i : Int32
  j : Int32
  accepted : Bool
  log_prob : Double
}

///|
/// Log Metropolis acceptance probability for swapping configurations
/// between replicas i and j.
/// `u_ii` = energy of config i on replica i, `u_ij` = config j on replica
/// i, `u_ji` = config i on replica j, `u_jj` = config j on replica j.
pub fn omm_hremd_log_accept(
  u_ii : Double,
  u_ij : Double,
  u_ji : Double,
  u_jj : Double,
  beta_i : Double,
  beta_j : Double,
) -> Double {
  let delta = beta_i * (u_ij - u_ii) + beta_j * (u_ji - u_jj)
  -delta
}

///|
/// Attempt one swap given a uniform random number in [0, 1).
pub fn omm_hremd_attempt(
  u_ii : Double,
  u_ij : Double,
  u_ji : Double,
  u_jj : Double,
  beta_i : Double,
  beta_j : Double,
  uniform~ : Double,
) -> Bool {
  let log_p = omm_hremd_log_accept(u_ii, u_ij, u_ji, u_jj, beta_i, beta_j)
  if uniform <= 0.0 {
    return log_p >= 0.0
  }
  @math.ln(uniform) < log_p
}

///|
/// One even-pair swap cycle over `n` replicas.
///
/// `energies[k][c]` is the energy of configuration c on replica k's
/// Hamiltonian. Configurations start assigned to the replica of the same
/// index (occupancy[c] = c). Pairs (0,1), (2,3), ... are attempted in
/// order; `uniforms` supplies one uniform deviate per pair (extra entries
/// are ignored; missing entries default to 0.5 = deterministic accept at
/// log_p >= 0).
///
/// Returns the swap records in attempt order.
pub fn omm_hremd_swap_cycle(
  energies : Array[Array[Double]],
  betas : Array[Double],
  uniforms~ : Array[Double] = [],
) -> Array[OmmHREMDSwap] {
  let n = energies.length()
  let occupancy : Array[Int32] = Array::make(n, 0)
  for k = 0; k < n; k = k + 1 {
    occupancy[k] = k
  }
  let swaps : Array[OmmHREMDSwap] = Array(capacity=0)
  let mut pair = 0
  let mut i = 0
  while i + 1 < n {
    let j = i + 1
    let ci = occupancy[i]
    let cj = occupancy[j]
    let u_ii = energies[i][ci]
    let u_ij = energies[i][cj]
    let u_ji = energies[j][ci]
    let u_jj = energies[j][cj]
    let log_p = omm_hremd_log_accept(u_ii, u_ij, u_ji, u_jj, betas[i], betas[j])
    let u = if pair < uniforms.length() { uniforms[pair] } else { 0.5 }
    let accepted = omm_hremd_attempt(
      u_ii,
      u_ij,
      u_ji,
      u_jj,
      betas[i],
      betas[j],
      uniform=u,
    )
    if accepted {
      occupancy[i] = cj
      occupancy[j] = ci
    }
    swaps.push({ i, j, accepted, log_prob: log_p, })
    i = i + 2
    pair = pair + 1
  }
  swaps
}

///|
/// Final configuration assignment after the last `swap_cycle` call on the
/// same occupancy bookkeeping. Returns occupancy[replica] = config index.
pub fn omm_hremd_final_occupancy(
  n_replicas : Int32,
  swaps : Array[OmmHREMDSwap],
) -> Array[Int32] {
  let occupancy : Array[Int32] = Array::make(n_replicas, 0)
  for k = 0; k < n_replicas; k = k + 1 {
    occupancy[k] = k
  }
  for s in swaps {
    if s.accepted && s.i < n_replicas && s.j < n_replicas {
      let tmp = occupancy[s.i]
      occupancy[s.i] = occupancy[s.j]
      occupancy[s.j] = tmp
    }
  }
  occupancy
}
