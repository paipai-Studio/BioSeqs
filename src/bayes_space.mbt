///|
/// BayesSpace module - Spatial transcriptomics clustering with spatial prior.
///
/// References:
///   - Zhao E et al. (2021) Spatial transcriptomics at subspot resolution with BayesSpace.
///     Nature Biotechnology 39:1375-1384.
///   - https://www.bioconductor.org/packages/release/bioc/html/BayesSpace.html
///
/// Core algorithm:
///   BayesSpace clusters spatial transcriptomics spots using a t-distributed mixture model
///   with a Markov random field (MRF) prior that encourages neighboring spots to share
///   cluster labels. The model is fit via Expectation-Maximization (EM):
///     - E-step: compute posterior responsibilities for each spot's cluster assignment
///       accounting for both expression likelihood (t-distribution) and the spatial prior
///       (Potts model) derived from spot neighbors.
///     - M-step: update cluster means, scales, and degrees of freedom.
///   This implementation provides a faithful, compact, deterministic (seeded) version
///   suitable for small datasets and unit tests.

///|
/// Spatial coordinates for a single spot.
pub struct SpotCoord {
  spot_id : String
  row : Int        // array row index
  col : Int        // array column index
  x : Double       // continuous x coordinate (e.g. pixel)
  y : Double       // continuous y coordinate
}

///|
/// Create a SpotCoord.
pub fn spot_coord(
  spot_id : String,
  row : Int,
  col : Int,
  x : Double,
  y : Double,
) -> SpotCoord {
  SpotCoord::{ spot_id, row, col, x, y }
}

///|
/// BayesSpace clustering result.
pub struct BayesSpaceResult {
  spot_ids : Array[String]
  clusters : Array[Int]            // cluster assignment per spot (0..q-1)
  q : Int                          // number of clusters
  cluster_centers : Array[Array[Double]]   // q × n_features matrix
  responsibilities : Array[Array[Double]]  // n_spots × q soft assignments
  n_iterations : Int
  log_likelihood : Double
  converged : Bool
}

///|
/// Build a hexagonal-grid neighbor list from spot row/col coordinates.
/// Two spots are neighbors if they are adjacent on a hex grid (6-connectivity):
/// for spot at (r, c), neighbors are (r±1, c), (r, c±1), and (r±1, c∓1) (offset varies).
/// For square-grid compatibility we also include (r±1, c±1) as neighbors.
///
/// Returns neighbor_indices[i] = list of spot indices that are neighbors of spot i.
pub fn bayes_space_hex_neighbors(
  spots : Array[SpotCoord],
) -> Array[Array[Int]] {
  let n = spots.length()
  // Build a map from (row, col) to spot index
  let coord_map : Map[String, Int] = Map::new()
  let mut i = 0
  while i < n {
    let key = spots[i].row.to_string() + "," + spots[i].col.to_string()
    coord_map.set(key, i)
    i = i + 1
  }
  let neighbors : Array[Array[Int]] = Array::new()
  let mut k = 0
  while k < n {
    neighbors.push(Array::new())
    k = k + 1
  }
  i = 0
  while i < n {
    let r = spots[i].row
    let c = spots[i].col
    // Hex 6-connectivity: (r, c-1), (r, c+1), (r-1, c), (r-1, c+1), (r+1, c-1), (r+1, c)
    // For square grid (4-connectivity), include only (r±1, c), (r, c±1).
    // We use 6-connectivity as BayesSpace default for Visium hex grids.
    let candidates = [
      (r, c - 1),
      (r, c + 1),
      (r - 1, c),
      (r - 1, c + 1),
      (r + 1, c - 1),
      (r + 1, c),
    ]
    let mut j = 0
    while j < candidates.length() {
      let (nr, nc) = candidates[j]
      let key = nr.to_string() + "," + nc.to_string()
      match coord_map.get(key) {
        Some(idx) => neighbors[i].push(idx)
        None => ()
      }
      j = j + 1
    }
    i = i + 1
  }
  neighbors
}

///|
/// Build a square-grid neighbor list (4-connectivity).
pub fn bayes_space_square_neighbors(
  spots : Array[SpotCoord],
) -> Array[Array[Int]] {
  let n = spots.length()
  let coord_map : Map[String, Int] = Map::new()
  let mut i = 0
  while i < n {
    let key = spots[i].row.to_string() + "," + spots[i].col.to_string()
    coord_map.set(key, i)
    i = i + 1
  }
  let neighbors : Array[Array[Int]] = Array::new()
  let mut k = 0
  while k < n {
    neighbors.push(Array::new())
    k = k + 1
  }
  i = 0
  while i < n {
    let r = spots[i].row
    let c = spots[i].col
    let candidates = [(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)]
    let mut j = 0
    while j < candidates.length() {
      let (nr, nc) = candidates[j]
      let key = nr.to_string() + "," + nc.to_string()
      match coord_map.get(key) {
        Some(idx) => neighbors[i].push(idx)
        None => ()
      }
      j = j + 1
    }
    i = i + 1
  }
  neighbors
}

///|
/// Simple deterministic pseudo-random number generator (Linear Congruential).
struct BayesRng {
  mut state : Int
}

///|
fn bayes_rng_next(rng : BayesRng) -> Int {
  rng.state = (rng.state * 1103515245 + 12345) & 0x7fffffff
  rng.state
}

///|
/// Return a uniform pseudo-random double in [0, 1).
fn bayes_rng_uniform(rng : BayesRng) -> Double {
  let v = bayes_rng_next(rng)
  v.to_double() / 2147483648.0
}

///|
/// Log of the Student's t-distribution density at point x given location mu,
/// scale sigma (so variance = sigma^2), and degrees of freedom df.
fn bayes_space_log_t_density(
  x : Array[Double],
  mu : Array[Double],
  sigma : Double,
  df : Double,
) -> Double {
  let d = x.length()
  let mut sq = 0.0
  let mut i = 0
  while i < d {
    let diff = x[i] - mu[i]
    sq = sq + diff * diff
    i = i + 1
  }
  // Multivariate t with diagonal covariance sigma^2 * I and df degrees of freedom.
  // log p(x) = -((df + d) / 2) * log(1 + sq / (df * sigma^2)) + const
  let sigma2 = sigma * sigma
  let arg = 1.0 + sq / (df * sigma2)
  let log_p = -((df + d.to_double()) / 2.0) * @math.ln(arg)
  log_p
}

///|
/// Initialize cluster centers using k-means++ seeding for well-separated starts.
fn bayes_space_init_centers(
  expression : Array[Array[Double]],
  q : Int,
  rng : BayesRng,
) -> Array[Array[Double]] {
  let n = expression.length()
  let d = if n > 0 { expression[0].length() } else { 0 }
  let centers : Array[Array[Double]] = Array::new()
  let mut ci = 0
  while ci < q {
    centers.push(Array::make(d, 0.0))
    ci = ci + 1
  }
  if n == 0 || q == 0 {
    return centers
  }
  // Pick first center randomly
  let first_idx = bayes_rng_next(rng) % n
  let mut j = 0
  while j < d {
    centers[0][j] = expression[first_idx][j]
    j = j + 1
  }
  // Pick subsequent centers using k-means++ (probability ∝ squared distance)
  let mut k = 1
  while k < q {
    // Compute min squared distance to any existing center for each spot
    let dists : Array[Double] = Array::make(n, 0.0)
    let mut total = 0.0
    let mut i = 0
    while i < n {
      let mut min_dist = bayes_space_euclidean(expression[i], centers[0])
      let mut c = 1
      while c < k {
        let dist = bayes_space_euclidean(expression[i], centers[c])
        if dist < min_dist {
          min_dist = dist
        }
        c = c + 1
      }
      let sq = min_dist * min_dist
      dists[i] = sq
      total = total + sq
      i = i + 1
    }
    // Pick the spot with the largest distance (deterministic version of k-means++)
    let mut best_idx = 0
    let mut best_dist = dists[0]
    i = 1
    while i < n {
      if dists[i] > best_dist {
        best_dist = dists[i]
        best_idx = i
      }
      i = i + 1
    }
    let mut j = 0
    while j < d {
      centers[k][j] = expression[best_idx][j]
      j = j + 1
    }
    k = k + 1
  }
  centers
}

///|
/// Compute Euclidean distance between two vectors.
fn bayes_space_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  let mut i = 0
  while i < a.length() {
    let d = a[i] - b[i]
    s = s + d * d
    i = i + 1
  }
  s.sqrt()
}

///|
/// Run BayesSpace clustering.
///
/// Parameters:
///   expression : n_spots × n_features matrix (typically PCA-reduced expression)
///   spots      : spot coordinates
///   q          : number of clusters
///   neighbors  : neighbor list per spot (e.g. from bayes_space_hex_neighbors)
///   max_iters  : maximum EM iterations
///   gamma      : spatial smoothing parameter (Potts coupling strength)
///   df         : degrees of freedom for the t-distribution
///   seed       : random seed for initialization
///   tol        : convergence tolerance on log-likelihood
pub fn bayes_space_run(
  expression : Array[Array[Double]],
  spots : Array[SpotCoord],
  q : Int,
  neighbors : Array[Array[Int]],
  max_iters? : Int = 50,
  gamma? : Double = 1.0,
  df? : Double = 6.0,
  seed? : Int = 42,
  tol? : Double = 1.0e-4,
) -> BayesSpaceResult {
  let n_spots = expression.length()
  let n_features = if n_spots > 0 { expression[0].length() } else { 0 }
  let max_it = max_iters
  let g = gamma
  let d_df = df
  let tolerance = tol

  let rng : BayesRng = { state: seed }

  // Initialize centers using k-means++-like seeding
  let centers = bayes_space_init_centers(expression, q, rng)
  // Initial scales: 1.0 for each cluster
  let scales : Array[Double] = Array::make(q, 1.0)
  // Initial mixing proportions: uniform
  let pis : Array[Double] = Array::make(q, 1.0 / q.to_double())

  // Initialize responsibilities by nearest-center assignment
  let resp : Array[Array[Double]] = Array::new()
  let mut ri = 0
  while ri < n_spots {
    resp.push(Array::make(q, 0.0))
    ri = ri + 1
  }
  let mut i = 0
  while i < n_spots {
    let mut best = 0
    let mut best_dist = bayes_space_euclidean(expression[i], centers[0])
    let mut k = 1
    while k < q {
      let dist = bayes_space_euclidean(expression[i], centers[k])
      if dist < best_dist {
        best_dist = dist
        best = k
      }
      k = k + 1
    }
    resp[i][best] = 1.0
    i = i + 1
  }

  let mut log_lik = -1.0e18
  let mut converged = false
  let mut iter = 0

  while iter < max_it {
    // ---- E-step: compute responsibilities with spatial prior ----
    // log_prob[i][k] = log(pi_k) + log_t_density(x_i | mu_k, sigma_k, df)
    //                 + gamma * sum over neighbors j of resp[j][k]  (Potts field)
    let log_prob : Array[Array[Double]] = Array::new()
    let mut lpi = 0
    while lpi < n_spots {
      log_prob.push(Array::make(q, 0.0))
      lpi = lpi + 1
    }
    i = 0
    while i < n_spots {
      let mut k = 0
      while k < q {
        let log_t = bayes_space_log_t_density(expression[i], centers[k], scales[k], d_df)
        let mut spatial = 0.0
        let mut ni = 0
        while ni < neighbors[i].length() {
          let j = neighbors[i][ni]
          spatial = spatial + resp[j][k]
          ni = ni + 1
        }
        log_prob[i][k] = @math.ln(pis[k]) + log_t + g * spatial
        k = k + 1
      }
      // Normalize via log-sum-exp
      let mut max_lp = log_prob[i][0]
      let mut k = 1
      while k < q {
        if log_prob[i][k] > max_lp {
          max_lp = log_prob[i][k]
        }
        k = k + 1
      }
      let mut sum_exp = 0.0
      k = 0
      while k < q {
        resp[i][k] = @math.exp(log_prob[i][k] - max_lp)
        sum_exp = sum_exp + resp[i][k]
        k = k + 1
      }
      if sum_exp > 1.0e-12 {
        k = 0
        while k < q {
          resp[i][k] = resp[i][k] / sum_exp
          k = k + 1
        }
      }
      i = i + 1
    }

    // ---- M-step: update parameters ----
    let mut k = 0
    while k < q {
      let mut nk = 0.0
      let new_center = Array::make(n_features, 0.0)
      i = 0
      while i < n_spots {
        let r = resp[i][k]
        nk = nk + r
        let mut j = 0
        while j < n_features {
          new_center[j] = new_center[j] + r * expression[i][j]
          j = j + 1
        }
        i = i + 1
      }
      if nk > 1.0e-12 {
        let mut j = 0
        while j < n_features {
          centers[k][j] = new_center[j] / nk
          j = j + 1
        }
        pis[k] = nk / n_spots.to_double()
        // Update scale: weighted average squared distance to center
        let mut sq_sum = 0.0
        i = 0
        while i < n_spots {
          let r = resp[i][k]
          let mut sq = 0.0
          let mut j = 0
          while j < n_features {
            let diff = expression[i][j] - centers[k][j]
            sq = sq + diff * diff
            j = j + 1
          }
          sq_sum = sq_sum + r * sq
          i = i + 1
        }
        // Scale = sqrt(weighted sq_sum / (nk * n_features))
        scales[k] = (sq_sum / (nk * n_features.to_double())).sqrt()
        if scales[k] < 1.0e-6 {
          scales[k] = 1.0e-6
        }
      } else {
        // Reinitialize empty cluster to a random spot
        let idx = bayes_rng_next(rng) % n_spots
        let mut j = 0
        while j < n_features {
          centers[k][j] = expression[idx][j]
          j = j + 1
        }
        pis[k] = 1.0 / n_spots.to_double()
        scales[k] = 1.0
      }
      k = k + 1
    }

    // ---- Compute log-likelihood ----
    let mut new_log_lik = 0.0
    i = 0
    while i < n_spots {
      let mut max_lp = -1.0e18
      let mut k = 0
      while k < q {
        let log_t = bayes_space_log_t_density(expression[i], centers[k], scales[k], d_df)
        let lp = @math.ln(pis[k]) + log_t
        if lp > max_lp {
          max_lp = lp
        }
        k = k + 1
      }
      let mut sum_exp = 0.0
      k = 0
      while k < q {
        let log_t = bayes_space_log_t_density(expression[i], centers[k], scales[k], d_df)
        sum_exp = sum_exp + @math.exp(@math.ln(pis[k]) + log_t - max_lp)
        k = k + 1
      }
      new_log_lik = new_log_lik + max_lp + @math.ln(sum_exp)
      i = i + 1
    }

    if (new_log_lik - log_lik).abs() < tolerance {
      log_lik = new_log_lik
      converged = true
      iter = iter + 1
      break
    }
    log_lik = new_log_lik
    iter = iter + 1
  }

  // Final hard assignments
  let clusters : Array[Int] = Array::make(n_spots, 0)
  let spot_ids : Array[String] = Array::make(n_spots, "")
  i = 0
  while i < n_spots {
    let mut best = 0
    let mut best_resp = resp[i][0]
    let mut k = 1
    while k < q {
      if resp[i][k] > best_resp {
        best_resp = resp[i][k]
        best = k
      }
      k = k + 1
    }
    clusters[i] = best
    spot_ids[i] = spots[i].spot_id
    i = i + 1
  }

  BayesSpaceResult::{
    spot_ids,
    clusters,
    q,
    cluster_centers: centers,
    responsibilities: resp,
    n_iterations: iter,
    log_likelihood: log_lik,
    converged,
  }
}

///|
/// Get cluster label for a spot by ID.
pub fn bayes_space_get_cluster(
  result : BayesSpaceResult,
  spot_id : String,
) -> Int {
  let mut i = 0
  while i < result.spot_ids.length() {
    if result.spot_ids[i] == spot_id {
      return result.clusters[i]
    }
    i = i + 1
  }
  -1
}

///|
/// Get all spot IDs assigned to a given cluster.
pub fn bayes_space_spots_in_cluster(
  result : BayesSpaceResult,
  cluster : Int,
) -> Array[String] {
  let out : Array[String] = []
  let mut i = 0
  while i < result.spot_ids.length() {
    if result.clusters[i] == cluster {
      out.push(result.spot_ids[i])
    }
    i = i + 1
  }
  out
}

///|
/// Count spots per cluster.
pub fn bayes_space_cluster_counts(result : BayesSpaceResult) -> Array[Int] {
  let counts = Array::make(result.q, 0)
  let mut i = 0
  while i < result.clusters.length() {
    counts[result.clusters[i]] = counts[result.clusters[i]] + 1
    i = i + 1
  }
  counts
}

///|
/// Render a 2D cluster map as ASCII art. Spots are placed by their (row, col)
/// coordinates; cells are labeled by their cluster (0-9 then A-Z).
pub fn bayes_space_render_clusters(
  result : BayesSpaceResult,
  spots : Array[SpotCoord],
) -> String {
  if spots.length() == 0 {
    return ""
  }
  let mut min_r = spots[0].row
  let mut max_r = spots[0].row
  let mut min_c = spots[0].col
  let mut max_c = spots[0].col
  let mut i = 1
  while i < spots.length() {
    if spots[i].row < min_r { min_r = spots[i].row }
    if spots[i].row > max_r { max_r = spots[i].row }
    if spots[i].col < min_c { min_c = spots[i].col }
    if spots[i].col > max_c { max_c = spots[i].col }
    i = i + 1
  }
  let label_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  let mut out = ""
  let mut r = min_r
  while r <= max_r {
    let mut c = min_c
    while c <= max_c {
      let mut ch = "."
      let mut i = 0
      while i < spots.length() {
        if spots[i].row == r && spots[i].col == c {
          let cl = result.clusters[i]
          if cl >= 0 && cl < label_chars.length() {
            ch = label_chars[cl].to_string()
          } else {
            ch = "?"
          }
          break
        }
        i = i + 1
      }
      out = out + ch + " "
      c = c + 1
    }
    out = out + "\n"
    r = r + 1
  }
  out
}

///|
/// Generate a small synthetic spatial dataset for testing and demos.
/// Returns (expression, spots) where spots are arranged on a 4×4 grid with
/// two spatially coherent clusters (top-left vs bottom-right).
pub fn bayes_space_sample_data() -> (Array[Array[Double]], Array[SpotCoord]) {
  let expression : Array[Array[Double]] = []
  let spots : Array[SpotCoord] = []
  let mut r = 0
  // Deterministic noise source
  let mut state = 7
  while r < 4 {
    let mut c = 0
    while c < 4 {
      // Top-left quadrant (r < 2, c < 2): cluster 0 (low values)
      // Bottom-right quadrant (r >= 2, c >= 2): cluster 1 (high values)
      // Boundaries are mixed
      let is_top_left = r < 2 && c < 2
      let is_bottom_right = r >= 2 && c >= 2
      state = (state * 1103515245 + 12345) & 0x7fffffff
      let noise = (state % 100).to_double() / 100.0 - 0.5
      let base = if is_top_left {
        0.5
      } else if is_bottom_right {
        3.5
      } else {
        2.0
      }
      // 3-feature expression vector per spot
      expression.push([base + noise, base * 0.8 + noise * 0.5, base * 1.2 + noise * 0.3])
      spots.push(spot_coord("spot_" + r.to_string() + "_" + c.to_string(), r, c, c.to_double(), r.to_double()))
      c = c + 1
    }
    r = r + 1
  }
  (expression, spots)
}
