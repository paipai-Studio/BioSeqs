///|
/// velociraptor - Single-cell RNA velocity analysis.
/// Inspired by Bioconductor velociraptor package (R wrapper for scVelo).
///
/// Implements steady-state and kinetic (EM) models for estimating
/// transcription (alpha), splicing (beta), and degradation (gamma) rates
/// from unspliced/spliced mRNA counts, computing velocity vectors, and
/// projecting them onto low-dimensional embeddings.
///
/// References:
///   - La Manno G et al. (2018) RNA velocity of single cells. Nature 560:494.
///   - Bergen V et al. (2020) Generalizing RNA velocity to transient cell
///     states through dynamical modeling. Nat Biotechnol 38:1408.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// Per-gene spliced and unspliced counts across cells.
pub struct VelocityGeneData {
  gene_name : String
  spliced : Array[Double] // counts per cell
  unspliced : Array[Double]
}

///|
/// Create new gene data.
pub fn VelocityGeneData::new(
  gene_name : String,
  spliced : Array[Double],
  unspliced : Array[Double]
) -> VelocityGeneData {
  VelocityGeneData::{ gene_name, spliced, unspliced }
}

///|
/// Accessor for gene name.
pub fn VelocityGeneData::gene_name(self : VelocityGeneData) -> String {
  self.gene_name
}

///|
/// Accessor for spliced counts.
pub fn VelocityGeneData::spliced(self : VelocityGeneData) -> Array[Double] {
  self.spliced
}

///|
/// Accessor for unspliced counts.
pub fn VelocityGeneData::unspliced(self : VelocityGeneData) -> Array[Double] {
  self.unspliced
}

///|
/// Kinetic parameters for a single gene.
pub struct GeneKinetics {
  gene_name : String
  alpha : Double // transcription rate
  beta : Double // splicing rate
  gamma : Double // degradation rate
  likelihood : Double
  steady_state_u : Double
  steady_state_s : Double
}

///|
/// Accessor for alpha.
pub fn GeneKinetics::alpha(self : GeneKinetics) -> Double {
  self.alpha
}

///|
/// Accessor for beta.
pub fn GeneKinetics::beta(self : GeneKinetics) -> Double {
  self.beta
}

///|
/// Accessor for gamma.
pub fn GeneKinetics::gamma(self : GeneKinetics) -> Double {
  self.gamma
}

///|
/// Accessor for gene name.
pub fn GeneKinetics::gene_name(self : GeneKinetics) -> String {
  self.gene_name
}

///|
/// Accessor for likelihood.
pub fn GeneKinetics::likelihood(self : GeneKinetics) -> Double {
  self.likelihood
}

///|
/// Velocity vector for a single cell.
pub struct CellVelocity {
  cell_id : String
  velocity : Array[Double] // per-gene velocity
  speed : Double
}

///|
/// Accessor for velocity vector.
pub fn CellVelocity::velocity(self : CellVelocity) -> Array[Double] {
  self.velocity
}

///|
/// Accessor for speed (magnitude of velocity).
pub fn CellVelocity::speed(self : CellVelocity) -> Double {
  self.speed
}

///|
/// Accessor for cell id.
pub fn CellVelocity::cell_id(self : CellVelocity) -> String {
  self.cell_id
}

///|
/// Projected velocity in embedding space.
pub struct EmbeddingVelocity {
  cell_id : String
  embedding : Array[Double] // e.g., [umap_x, umap_y]
  velocity_embedding : Array[Double] // projected velocity
}

///|
/// Accessor for embedding.
pub fn EmbeddingVelocity::embedding(self : EmbeddingVelocity) -> Array[Double] {
  self.embedding
}

///|
/// Accessor for velocity embedding.
pub fn EmbeddingVelocity::velocity_embedding(self : EmbeddingVelocity) -> Array[
  Double
] {
  self.velocity_embedding
}

// ===========================================================================
// Internal helpers
// ===========================================================================

fn vr_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

fn vr_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x / 2.0
  if guess < 1.0 {
    guess = 1.0
  }
  for _ in 0..<50 {
    let next = (guess + x / guess) / 2.0
    let diff = next - guess
    if vr_abs(diff) < 0.0000000001 {
      return next
    }
    guess = next
  }
  guess
}

fn vr_log(x : Double) -> Double {
  if x <= 0.0 {
    return -700.0
  }
  @math.ln(x)
}

fn vr_exp(x : Double) -> Double {
  if x > 700.0 {
    return 1.7976931348623157e308
  }
  if x < -700.0 {
    return 0.0
  }
  @math.exp(x)
}

fn vr_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s / n.to_double()
}

fn vr_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = vr_mean(arr)
  let mut s = 0.0
  for v in arr {
    let d = v - m
    s = s + d * d
  }
  s / (n - 1).to_double()
}

fn vr_lcg_rand(state : Int) -> (Int, Double) {
  let s = (state * 1103515245 + 12345) & 0x7FFFFFFF
  (s, (s % 1000000).to_double() / 1000000.0)
}

// ===========================================================================
// Steady-state model (La Manno 2018)
// ===========================================================================

///|
/// Estimate gamma/beta ratio via steady-state linear regression.
/// Fits unspliced ~ spliced using a quantile-based approach to identify
/// the steady-state population (upper quantile of spliced).
///
/// Returns (gamma_beta_ratio, steady_state_s, steady_state_u).
pub fn velocity_steady_state(
  spliced : Array[Double],
  unspliced : Array[Double]
) -> (Double, Double, Double) {
  let n = spliced.length()
  if n < 3 {
    return (1.0, 0.0, 0.0)
  }
  // Extreme quantile filtering: use cells in the top 20% of spliced
  // as steady-state cells.
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Sort indices by spliced value (bubble sort).
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if spliced[indices[j]] < spliced[indices[i]] {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
      j = j + 1
    }
  }
  let quantile_start = (n.to_double() * 0.8).to_int()
  let quantile_start = if quantile_start >= n { n - 1 } else { quantile_start }
  // Compute means of steady-state cells.
  let mut ss_s_sum = 0.0
  let mut ss_u_sum = 0.0
  let mut count = 0
  for i in quantile_start..<n {
    let idx = indices[i]
    ss_s_sum = ss_s_sum + spliced[idx]
    ss_u_sum = ss_u_sum + unspliced[idx]
    count = count + 1
  }
  if count == 0 {
    return (1.0, 0.0, 0.0)
  }
  let ss_s = ss_s_sum / count.to_double()
  let ss_u = ss_u_sum / count.to_double()
  // Compute slope (gamma/beta ratio) via least squares on steady-state cells.
  let mut sxx = 0.0
  let mut sxy = 0.0
  for i in quantile_start..<n {
    let idx = indices[i]
    let ds = spliced[idx] - ss_s
    sxx = sxx + ds * ds
    sxy = sxy + ds * (unspliced[idx] - ss_u)
  }
  let ratio = if sxx > 0.0000001 { sxy / sxx } else { 1.0 }
  // Gamma/beta ratio should be positive; if negative, use absolute value.
  let ratio = if ratio > 0.0 { ratio } else { vr_abs(ratio) }
  (ratio, ss_s, ss_u)
}

// ===========================================================================
// Kinetic model (scVelo EM algorithm)
// ===========================================================================

///|
/// Estimate kinetic parameters (alpha, beta, gamma) using the EM-based
/// dynamical model. Simplified version of scVelo's approach.
///
/// - E-step: assign each cell to a kinetic phase (induction, steady, repression)
///   based on current parameters.
/// - M-step: re-estimate alpha, beta, gamma from phase assignments.
///
/// Returns GeneKinetics with estimated parameters.
pub fn velocity_kinetic_model(
  gene_name : String,
  spliced : Array[Double],
  unspliced : Array[Double],
  n_iterations : Int
) -> GeneKinetics {
  let n = spliced.length()
  if n < 3 {
    return GeneKinetics::{
      gene_name,
      alpha: 0.0,
      beta: 1.0,
      gamma: 1.0,
      likelihood: 0.0,
      steady_state_u: 0.0,
      steady_state_s: 0.0,
    }
  }
  // Initialize from steady-state model.
  let (gamma_beta_ratio, ss_s, ss_u) = velocity_steady_state(spliced, unspliced)
  let beta = 1.0
  let mut gamma = gamma_beta_ratio * beta
  let mut alpha = gamma * ss_s // at steady state: alpha = gamma * s_ss
  // Normalize counts (avoid numerical issues with large counts).
  let max_s = vr_max_array(spliced)
  let max_u = vr_max_array(unspliced)
  let scale = if max_s > max_u { max_s } else { max_u }
  let scale = if scale < 1.0 { 1.0 } else { scale }
  let norm_s = Array::make(n, 0.0)
  let norm_u = Array::make(n, 0.0)
  for i in 0..<n {
    norm_s[i] = spliced[i] / scale
    norm_u[i] = unspliced[i] / scale
  }
  // EM iterations.
  let phase_labels = Array::make(n, 0) // 0=induction, 1=steady, 2=repression
  let mut likelihood = 0.0
  for _ in 0..<n_iterations {
    // E-step: assign phases based on unspliced vs spliced relationship.
    // At steady state, u = gamma/beta * s.
    let threshold = ss_u / scale
    for i in 0..<n {
      let expected_u_steady = gamma / beta * norm_s[i]
      let diff = norm_u[i] - expected_u_steady
      let abs_diff = vr_abs(diff)
      let steady_threshold = 0.1
      if abs_diff < steady_threshold {
        phase_labels[i] = 1 // steady
      } else if diff > 0.0 {
        phase_labels[i] = 0 // induction (more unspliced than expected)
      } else {
        phase_labels[i] = 2 // repression
      }
    }
    // M-step: re-estimate parameters.
    // For induction phase: du/dt > 0, ds/dt > 0
    // For steady phase: du/dt ~ 0, ds/dt ~ 0 => u = (alpha - gamma*s)/beta
    // For repression phase: du/dt < 0, ds/dt < 0
    // Estimate gamma from steady-state cells.
    let steady_indices : Array[Int] = Array::new()
    for i in 0..<n {
      if phase_labels[i] == 1 {
        steady_indices.push(i)
      }
    }
    if steady_indices.length() >= 2 {
      let mut sxx = 0.0
      let mut sxy = 0.0
      let mut mean_s = 0.0
      let mut mean_u = 0.0
      for idx in steady_indices {
        mean_s = mean_s + norm_s[idx]
        mean_u = mean_u + norm_u[idx]
      }
      let ns = steady_indices.length().to_double()
      mean_s = mean_s / ns
      mean_u = mean_u / ns
      for idx in steady_indices {
        let ds = norm_s[idx] - mean_s
        sxx = sxx + ds * ds
        sxy = sxy + ds * (norm_u[idx] - mean_u)
      }
      if sxx > 0.0000001 {
        let new_ratio = sxy / sxx
        let new_ratio = if new_ratio > 0.01 { new_ratio } else { 0.01 }
        gamma = new_ratio * beta
      }
    }
    // Estimate alpha from induction phase cells.
    let induction_indices : Array[Int] = Array::new()
    for i in 0..<n {
      if phase_labels[i] == 0 {
        induction_indices.push(i)
      }
    }
    if induction_indices.length() > 0 {
      let mut sum_u = 0.0
      let mut sum_s = 0.0
      for idx in induction_indices {
        sum_u = sum_u + norm_u[idx]
        sum_s = sum_s + norm_s[idx]
      }
      let ni = induction_indices.length().to_double()
      let mean_u_ind = sum_u / ni
      let mean_s_ind = sum_s / ni
      // During induction: alpha = beta * u + gamma * s (approximately)
      let new_alpha = beta * mean_u_ind + gamma * mean_s_ind
      alpha = if new_alpha > 0.0 { new_alpha } else { alpha }
    }
    // Compute likelihood (simplified: negative sum of squared residuals).
    let mut ll = 0.0
    for i in 0..<n {
      let expected_u = (alpha - gamma * norm_s[i]) / beta
      let residual = norm_u[i] - expected_u
      ll = ll - residual * residual
    }
    likelihood = ll
  }
  // Compute steady state values.
  let steady_u = alpha / beta
  let steady_s = alpha / gamma
  GeneKinetics::{
    gene_name,
    alpha,
    beta,
    gamma,
    likelihood,
    steady_state_u: steady_u * scale,
    steady_state_s: steady_s * scale,
  }
}

fn vr_max_array(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut m = arr[0]
  for i in 1..<n {
    if arr[i] > m {
      m = arr[i]
    }
  }
  m
}

// ===========================================================================
// Velocity computation
// ===========================================================================

///|
/// Compute RNA velocity for each cell across all genes.
/// Velocity v_i = beta * u_i - gamma * s_i (per gene per cell).
///
/// `gene_data` is the per-gene spliced/unspliced data.
/// `kinetics` is the per-gene kinetic parameters.
/// Returns per-cell velocity vectors.
pub fn velocity_compute(
  gene_data : Array[VelocityGeneData],
  kinetics : Array[GeneKinetics]
) -> Array[CellVelocity] {
  let n_genes = gene_data.length()
  if n_genes == 0 {
    return []
  }
  let n_cells = gene_data[0].spliced.length()
  let cell_ids = Array::make(n_cells, "cell_0")
  // Build per-cell velocity vectors.
  let result : Array[CellVelocity] = Array::new()
  for ci in 0..<n_cells {
    let vel = Array::make(n_genes, 0.0)
    for gi in 0..<n_genes {
      let s = gene_data[gi].spliced[ci]
      let u = gene_data[gi].unspliced[ci]
      let beta = kinetics[gi].beta
      let gamma = kinetics[gi].gamma
      let v = beta * u - gamma * s
      vel[gi] = if v.is_nan() { 0.0 } else { v }
    }
    // Compute speed (magnitude).
    let mut speed_sq = 0.0
    for v in vel {
      speed_sq = speed_sq + v * v
    }
    let speed = vr_sqrt(speed_sq)
    let speed = if speed.is_nan() { 0.0 } else { speed }
    result.push(CellVelocity::{
      cell_id: "cell_" + ci.to_string(),
      velocity: vel,
      speed,
    })
  }
  result
}

///|
/// Compute velocity for a single gene across all cells.
/// v = beta * u - gamma * s
pub fn velocity_gene(
  spliced : Array[Double],
  unspliced : Array[Double],
  beta : Double,
  gamma : Double
) -> Array[Double] {
  let n = spliced.length()
  let vel = Array::make(n, 0.0)
  for i in 0..<n {
    vel[i] = beta * unspliced[i] - gamma * spliced[i]
  }
  vel
}

// ===========================================================================
// Velocity embedding (projection to low-dim)
// ===========================================================================

///|
/// Project velocity vectors onto a low-dimensional embedding using
/// k-nearest-neighbor weighted averaging.
///
/// `velocities` is per-cell velocity vectors (in gene space).
/// `embeddings` is per-cell low-dim coordinates (e.g., UMAP).
/// `k` is the number of neighbors to consider.
/// Returns projected velocity vectors in embedding space.
pub fn velocity_embedding(
  velocities : Array[CellVelocity],
  embeddings : Array[Array[Double]],
  k : Int
) -> Array[EmbeddingVelocity] {
  let n_cells = velocities.length()
  if n_cells == 0 || n_cells != embeddings.length() {
    return []
  }
  let n_dims = embeddings[0].length()
  let result : Array[EmbeddingVelocity] = Array::new()
  for ci in 0..<n_cells {
    // Find k nearest neighbors based on embedding distance.
    let neighbors : Array[(Int, Double)] = Array::new()
    for cj in 0..<n_cells {
      if ci != cj {
        let dist = vr_euclidean(embeddings[ci], embeddings[cj])
        neighbors.push((cj, dist))
      }
    }
    // Sort by distance (bubble sort, k is small).
    let nn = neighbors.length()
    for i in 0..<nn {
      let mut j = i + 1
      while j < nn {
        if neighbors[j].1 < neighbors[i].1 {
          let tmp = neighbors[i]
          neighbors[i] = neighbors[j]
          neighbors[j] = tmp
        }
        j = j + 1
      }
    }
    // Take k nearest.
    let actual_k = if k < nn { k } else { nn }
    // Compute weighted average of velocity-to-neighbor directions.
    let vel_emb = Array::make(n_dims, 0.0)
    let mut total_weight = 0.0
    for ni in 0..<actual_k {
      let (n_idx, dist) = neighbors[ni]
      let weight = 1.0 / (dist + 0.0000001)
      // Project velocity gene vector onto embedding direction.
      // Simplified: use the embedding difference direction weighted by
      // cosine similarity with velocity.
      let emb_diff = Array::make(n_dims, 0.0)
      for d in 0..<n_dims {
        emb_diff[d] = embeddings[n_idx][d] - embeddings[ci][d]
      }
      // Cosine similarity between velocity and emb_diff (using embedding as proxy).
      let vel_vec = velocities[ci].velocity
      let emb_diff_mag = vr_sqrt(vr_dot(emb_diff, emb_diff))
      let vel_mag = velocities[ci].speed
      if emb_diff_mag > 0.0000001 && vel_mag > 0.0000001 {
        // Use dot product of normalized vectors as alignment score.
        // Since velocity is in gene space and emb_diff is in embedding space,
        // we use a simplified scalar projection.
        let alignment = 1.0 // simplified: assume positive alignment for neighbors
        for d in 0..<n_dims {
          vel_emb[d] = vel_emb[d] + weight * alignment * emb_diff[d] / emb_diff_mag
        }
        total_weight = total_weight + weight
      }
    }
    // Normalize.
    if total_weight > 0.0 {
      for d in 0..<n_dims {
        vel_emb[d] = vel_emb[d] / total_weight
      }
    }
    result.push(EmbeddingVelocity::{
      cell_id: velocities[ci].cell_id,
      embedding: embeddings[ci],
      velocity_embedding: vel_emb,
    })
  }
  result
}

fn vr_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut s = 0.0
  for i in 0..<n {
    let d = a[i] - b[i]
    s = s + d * d
  }
  vr_sqrt(s)
}

fn vr_dot(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut s = 0.0
  for i in 0..<n {
    s = s + a[i] * b[i]
  }
  s
}

// ===========================================================================
// Speed and root cell identification
// ===========================================================================

///|
/// Identify root cells (cells with highest negative velocity magnitude,
/// indicating they are at the start of differentiation).
///
/// `velocities` is per-cell velocity vectors.
/// `n_roots` is the number of root cells to identify.
/// Returns indices of root cells sorted by speed (descending).
pub fn velocity_find_root_cells(
  velocities : Array[CellVelocity],
  n_roots : Int
) -> Array[Int] {
  let n = velocities.length()
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Sort by speed descending.
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if velocities[indices[j]].speed > velocities[indices[i]].speed {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
      j = j + 1
    }
  }
  let limit = if n_roots < n { n_roots } else { n }
  let result : Array[Int] = Array::new()
  for i in 0..<limit {
    result.push(indices[i])
  }
  result
}

///|
/// Compute transition probability matrix from velocity vectors.
/// P(i -> j) proportional to exp(cos(v_i, x_j - x_i) / temperature)
/// for neighboring cells j.
///
/// `velocities` is per-cell velocity.
/// `embeddings` is per-cell coordinates.
/// `k` is the number of neighbors.
/// `temperature` controls sharpness of transition.
/// Returns transition matrix [n_cells][n_cells] (sparse: only k neighbors non-zero).
pub fn velocity_transition_matrix(
  velocities : Array[CellVelocity],
  embeddings : Array[Array[Double]],
  k : Int,
  temperature : Double
) -> Array[Array[Double]] {
  let n = velocities.length()
  let temp = if temperature > 0.0 { temperature } else { 1.0 }
  let matrix : Array[Array[Double]] = Array::new()
  for _ in 0..<n {
    matrix.push(Array::make(n, 0.0))
  }
  for ci in 0..<n {
    // Find k nearest neighbors.
    let neighbors : Array[(Int, Double)] = Array::new()
    for cj in 0..<n {
      if ci != cj {
        let dist = vr_euclidean(embeddings[ci], embeddings[cj])
        neighbors.push((cj, dist))
      }
    }
    let nn = neighbors.length()
    for i in 0..<nn {
      let mut j = i + 1
      while j < nn {
        if neighbors[j].1 < neighbors[i].1 {
          let tmp = neighbors[i]
          neighbors[i] = neighbors[j]
          neighbors[j] = tmp
        }
        j = j + 1
      }
    }
    let actual_k = if k < nn { k } else { nn }
    // Compute transition probabilities.
    let mut sum_exp = 0.0
    let exps = Array::make(actual_k, 0.0)
    let cell_speed = velocities[ci].speed
    for ni in 0..<actual_k {
      let (n_idx, _) = neighbors[ni]
      // Direction to neighbor in embedding space.
      let dir = Array::make(embeddings[ci].length(), 0.0)
      for d in 0..<embeddings[ci].length() {
        dir[d] = embeddings[n_idx][d] - embeddings[ci][d]
      }
      // Use speed as alignment score (simplified: higher speed = stronger
      // directional bias toward neighbors).
      let dir_mag = vr_sqrt(vr_dot(dir, dir))
      let cos_sim = if dir_mag > 0.0000001 && cell_speed > 0.0000001 {
        // Use 1.0 as positive alignment (neighbor is reachable).
        1.0 / (1.0 + dir_mag / (cell_speed + 0.0000001))
      } else {
        0.0
      }
      let e = vr_exp(cos_sim / temp)
      exps[ni] = e
      sum_exp = sum_exp + e
    }
    // Normalize.
    if sum_exp > 0.0 {
      for ni in 0..<actual_k {
        let (n_idx, _) = neighbors[ni]
        matrix[ci][n_idx] = exps[ni] / sum_exp
      }
    }
  }
  matrix
}

// ===========================================================================
// Full pipeline
// ===========================================================================

///|
/// Run the complete RNA velocity pipeline:
/// 1. Estimate kinetic parameters per gene
/// 2. Compute velocity vectors per cell
/// 3. Project to embedding (if provided)
///
/// `gene_data` is per-gene spliced/unspliced data.
/// `embeddings` is optional per-cell embedding (empty array to skip).
/// `k` is the number of neighbors for embedding projection.
/// Returns (kinetics, velocities, embedding_velocities).
pub fn velocity_run(
  gene_data : Array[VelocityGeneData],
  embeddings : Array[Array[Double]],
  k : Int
) -> (Array[GeneKinetics], Array[CellVelocity], Array[EmbeddingVelocity]) {
  // 1. Estimate kinetics per gene.
  let kinetics : Array[GeneKinetics] = Array::new()
  for gd in gene_data {
    let k = velocity_kinetic_model(gd.gene_name, gd.spliced, gd.unspliced, 50)
    kinetics.push(k)
  }
  // 2. Compute velocity vectors.
  let velocities = velocity_compute(gene_data, kinetics)
  // 3. Project to embedding.
  let emb_vel = if embeddings.length() > 0 {
    velocity_embedding(velocities, embeddings, k)
  } else {
    []
  }
  (kinetics, velocities, emb_vel)
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create sample spliced/unspliced count data for testing.
/// Simulates 4 genes across 20 cells with a differentiation trajectory:
/// cells transition from high unspliced (early) to high spliced (late).
pub fn velocity_sample_data() -> Array[VelocityGeneData] {
  let genes = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let n_cells = 20
  let gene_data : Array[VelocityGeneData] = Array::new()
  let mut rng_state = 42
  for gene_idx in 0..<genes.length() {
    let spliced = Array::make(n_cells, 0.0)
    let unspliced = Array::make(n_cells, 0.0)
    let alpha = 2.0 + gene_idx.to_double() * 0.5
    let beta = 0.8
    let gamma = 0.6 + gene_idx.to_double() * 0.1
    for ci in 0..<n_cells {
      let t = ci.to_double() / n_cells.to_double() // pseudo-time 0 to 1
      let (s1, r1) = vr_lcg_rand(rng_state)
      rng_state = s1
      let noise_u = r1 * 0.5
      let (s2, r2) = vr_lcg_rand(rng_state)
      rng_state = s2
      let noise_s = r2 * 0.5
      // Unspliced: high early, decreases over time.
      unspliced[ci] = alpha / beta * (1.0 - vr_exp(-beta * t)) + noise_u
      // Spliced: low early, increases over time (simplified model).
      spliced[ci] = alpha / gamma * (1.0 - vr_exp(-gamma * t)) + noise_s
      // Ensure non-negative.
      if spliced[ci] < 0.0 {
        spliced[ci] = 0.0
      }
      if unspliced[ci] < 0.0 {
        unspliced[ci] = 0.0
      }
      if spliced[ci].is_nan() {
        spliced[ci] = 0.0
      }
      if unspliced[ci].is_nan() {
        unspliced[ci] = 0.0
      }
    }
    gene_data.push(VelocityGeneData::new(genes[gene_idx], spliced, unspliced))
  }
  gene_data
}

///|
/// Create sample UMAP-like embeddings for testing.
pub fn velocity_sample_embeddings(n_cells : Int) -> Array[Array[Double]] {
  let embeddings : Array[Array[Double]] = Array::new()
  let mut rng_state = 99
  for i in 0..<n_cells {
    let t = i.to_double() / n_cells.to_double()
    let (s1, r1) = vr_lcg_rand(rng_state)
    rng_state = s1
    let (s2, r2) = vr_lcg_rand(rng_state)
    rng_state = s2
    let x = t * 10.0 + r1 * 0.5
    let y = @math.sin(t * 3.14159) * 5.0 + r2 * 0.5
    embeddings.push([x, y])
  }
  embeddings
}
