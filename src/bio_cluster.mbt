///|
/// Bio.Cluster - Numerical clustering library, ported from Biopython's
/// Bio.Cluster module (the C Clustering Library, de Hoon et al. 2004).
///
/// Implements:
///   - 8 distance metrics (Euclidean, City Block, (absolute/uncentered)
///     Pearson correlation, Spearman rank correlation, Kendall's tau)
///   - Hierarchical clustering (single/complete/average/centroid linkage)
///   - k-means / k-median clustering (kcluster)
///   - k-medoids clustering
///   - Self-organizing maps (somcluster)
///   - Principal component analysis (pca, via Golub-Reinsch SVD)
///   - Cluster centroids and cluster distances
///   - Cluster/TreeView (Michael Eisen) file format parsing (Record)

// ---------------------------------------------------------------------------
// Constants and small utilities
// ---------------------------------------------------------------------------

///|
let bc_DBL_MAX : Double = 1.7976931348623157e308

///|
fn bc_max(a : Double, b : Double) -> Double {
  if a > b {
    a
  } else {
    b
  }
}

///|
/// Parse a floating point number from a string (mirrors Python float()).
fn bc_parse_double(s : String) -> Double {
  let sv = s.trim().to_string()
  if sv.length() == 0 {
    return 0.0
  }
  let n = sv.length()
  let mut result = 0.0
  let mut decimal = 0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  if sv.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  } else if sv.unsafe_get(0).to_int() == '+'.to_int() {
    i = 1
  }
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      break
    }
    result = result * 10.0 + (c - zero).to_double()
    i = i + 1
  }
  if i < n && sv.unsafe_get(i) == '.' {
    i = i + 1
    while i < n {
      let c = sv.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 {
        break
      }
      result = result * 10.0 + (c - zero).to_double()
      decimal = decimal + 1
      i = i + 1
    }
  }
  while decimal > 0 {
    result = result / 10.0
    decimal = decimal - 1
  }
  // Handle scientific notation (e.g., 1.5e-30)
  let mut exponent = 0
  let mut exp_sign = 1
  if i < n && (sv.unsafe_get(i) == 'e' || sv.unsafe_get(i) == 'E') {
    i = i + 1
    if i < n && sv.unsafe_get(i) == '-' {
      exp_sign = -1
      i = i + 1
    } else if i < n && sv.unsafe_get(i) == '+' {
      i = i + 1
    }
    while i < n {
      let c = sv.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 {
        break
      }
      exponent = exponent * 10 + (c - zero)
      i = i + 1
    }
    exponent = exponent * exp_sign
  }
  result * sign * @math.pow(10.0, exponent.to_double())
}

// ---------------------------------------------------------------------------
// Deterministic random number generator
//
// The C Clustering Library uses random initial assignments for k-means,
// k-medoids and SOM. Here a Park-Miller style generator (Schrage's method)
// is used so that results are reproducible; use cluster_seed() to control it.
// ---------------------------------------------------------------------------

///|
let bc_rng_state : Ref[Int] = Ref::new(987654321)

///|
/// Seed the internal random number generator used by kcluster (npass > 1),
/// kmedoids (npass > 1) and somcluster. The generator is deterministic.
pub fn cluster_seed(seed : Int) -> Unit {
  bc_rng_state.val = if seed <= 0 { 987654321 } else { seed % 2147483647 }
}

///|
/// Uniform random number in (0, 1).
fn bc_uniform() -> Double {
  let a = 48271
  let m = 2147483647
  let q = 44488
  let r = 3399
  let s = bc_rng_state.val
  let hi = s / q
  let lo = s % q
  let mut t = a * lo - r * hi
  if t <= 0 {
    t = t + m
  }
  bc_rng_state.val = t
  t.to_double() / m.to_double()
}

///|
/// Binomially distributed random number (counting method).
fn bc_binomial(n : Int, p : Double) -> Int {
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    if bc_uniform() < p {
      count = count + 1
    }
  }
  count
}

///|
/// Randomly assign elements to clusters, making sure each cluster receives
/// at least one element (port of randomassign in cluster.c).
fn bc_randomassign(
  nclusters : Int,
  nelements : Int,
  clusterid : Array[Int],
) -> Unit {
  let mut k = 0
  let mut n = nelements - nclusters
  let mut i = 0
  // Draw the number of elements in each cluster from a multinomial
  // distribution, reserving nclusters elements to set independently.
  while i < nclusters - 1 {
    let p = 1.0 / (nclusters - i).to_double()
    let mut j = bc_binomial(n, p)
    n = n - j
    j = j + k + 1 // Assign at least one element to cluster i
    while k < j {
      clusterid[k] = i
      k = k + 1
    }
    i = i + 1
  }
  // Assign the remaining elements to the last cluster
  while k < nelements {
    clusterid[k] = i
    k = k + 1
  }
  // Create a random permutation of the cluster assignments
  for i2 = 0; i2 < nelements; i2 = i2 + 1 {
    let j = i2 + ((nelements - i2).to_double() * bc_uniform()).to_int()
    let t = clusterid[j]
    clusterid[j] = clusterid[i2]
    clusterid[i2] = t
  }
}

// ---------------------------------------------------------------------------
// Basic statistics: median, sort_index, getrank
// ---------------------------------------------------------------------------

///|
/// Median of the first n values in x (result identical to median() in
/// cluster.c; implemented via sorting for clarity).
fn bc_median(n : Int, x : Array[Double]) -> Double {
  if n < 1 {
    return 0.0
  }
  if n == 1 {
    return x[0]
  }
  if n == 2 {
    return 0.5 * (x[0] + x[1])
  }
  let arr : Array[Double] = Array::new()
  for i = 0; i < n; i = i + 1 {
    arr.push(x[i])
  }
  arr.sort()
  let nr = n / 2
  if n % 2 == 1 {
    arr[nr]
  } else {
    0.5 * (arr[nr - 1] + arr[nr])
  }
}

///|
/// Fill index with the positions of the data sorted in increasing order.
fn bc_sort_index(n : Int, data : Array[Double], index : Array[Int]) -> Unit {
  for i = 0; i < n; i = i + 1 {
    index[i] = i
  }
  index.sort_by(fn(a : Int, b : Int) -> Int {
    if data[a] < data[b] {
      -1
    } else if data[a] > data[b] {
      1
    } else {
      0
    }
  })
}

///|
/// Calculate the ranks of the elements in the data array. Elements with the
/// same value get the same rank, equal to the average of the ranks had the
/// elements been different (port of getrank in cluster.c).
fn bc_getrank(
  n : Int,
  data : Array[Double],
  weight : Array[Double],
) -> Array[Double] {
  let rank : Array[Double] = Array::make(n, 0.0)
  if n < 1 {
    return rank
  }
  let index : Array[Int] = Array::make(n, 0)
  bc_sort_index(n, data, index)
  let mut total = 0.0
  let mut k = 0
  let mut current = data[index[0]]
  let mut subtotal = weight[index[0]]
  for i = 1; i < n; i = i + 1 {
    let j = index[i]
    let value = data[j]
    if value != current {
      current = value
      let r = total + (subtotal + 1.0) / 2.0
      for l = k; l < i; l = l + 1 {
        rank[index[l]] = r
      }
      k = i
      total = total + subtotal
      subtotal = 0.0
    }
    subtotal = subtotal + weight[j]
  }
  let r = total + (subtotal + 1.0) / 2.0
  for l = k; l < n; l = l + 1 {
    rank[index[l]] = r
  }
  rank
}

// ---------------------------------------------------------------------------
// Distance metrics
// ---------------------------------------------------------------------------

///|
/// Weighted Euclidean distance between two rows/columns (port of euclid).
fn bc_euclid(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  let mut result = 0.0
  let mut tweight = 0.0
  if transpose == 0 {
    for i = 0; i < n; i = i + 1 {
      if mask1[index1][i] != 0 && mask2[index2][i] != 0 {
        let term = data1[index1][i] - data2[index2][i]
        result = result + weight[i] * term * term
        tweight = tweight + weight[i]
      }
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      if mask1[i][index1] != 0 && mask2[i][index2] != 0 {
        let term = data1[i][index1] - data2[i][index2]
        result = result + weight[i] * term * term
        tweight = tweight + weight[i]
      }
    }
  }
  if tweight == 0.0 {
    return 0.0 // usually due to empty clusters
  }
  result / tweight
}

///|
/// Weighted "City Block" distance between two rows/columns (port of cityblock).
fn bc_cityblock(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  let mut result = 0.0
  let mut tweight = 0.0
  if transpose == 0 {
    for i = 0; i < n; i = i + 1 {
      if mask1[index1][i] != 0 && mask2[index2][i] != 0 {
        let term = data1[index1][i] - data2[index2][i]
        result = result + weight[i] * term.abs()
        tweight = tweight + weight[i]
      }
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      if mask1[i][index1] != 0 && mask2[i][index2] != 0 {
        let term = data1[i][index1] - data2[i][index2]
        result = result + weight[i] * term.abs()
        tweight = tweight + weight[i]
      }
    }
  }
  if tweight == 0.0 {
    return 0.0
  }
  result / tweight
}

///|
/// Compute the weighted Pearson-correlation-based distance; `absolute`
/// selects |r| and `uncentered` uses zero means (port of correlation,
/// acorrelation, ucorrelation and uacorrelation).
fn bc_correlation_core(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
  absolute~ : Bool,
  uncentered~ : Bool,
) -> Double {
  let mut result = 0.0
  let mut sum1 = 0.0
  let mut sum2 = 0.0
  let mut denom1 = 0.0
  let mut denom2 = 0.0
  let mut tweight = 0.0
  let mut flag = false
  if transpose == 0 {
    for i = 0; i < n; i = i + 1 {
      if mask1[index1][i] != 0 && mask2[index2][i] != 0 {
        let term1 = data1[index1][i]
        let term2 = data2[index2][i]
        let w = weight[i]
        if !uncentered {
          sum1 = sum1 + w * term1
          sum2 = sum2 + w * term2
        }
        result = result + w * term1 * term2
        denom1 = denom1 + w * term1 * term1
        denom2 = denom2 + w * term2 * term2
        tweight = tweight + w
        flag = true
      }
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      if mask1[i][index1] != 0 && mask2[i][index2] != 0 {
        let term1 = data1[i][index1]
        let term2 = data2[i][index2]
        let w = weight[i]
        if !uncentered {
          sum1 = sum1 + w * term1
          sum2 = sum2 + w * term2
        }
        result = result + w * term1 * term2
        denom1 = denom1 + w * term1 * term1
        denom2 = denom2 + w * term2 * term2
        tweight = tweight + w
        flag = true
      }
    }
  }
  if uncentered {
    // ucorrelation / uacorrelation
    if !flag {
      return 0.0
    }
    if denom1 == 0.0 {
      return 1.0
    }
    if denom2 == 0.0 {
      return 1.0
    }
    result = result / (denom1 * denom2).sqrt()
    if absolute {
      result = result.abs()
    }
    return 1.0 - result
  }
  if tweight == 0.0 {
    return 0.0
  }
  result = result - sum1 * sum2 / tweight
  denom1 = denom1 - sum1 * sum1 / tweight
  denom2 = denom2 - sum2 * sum2 / tweight
  if denom1 <= 0.0 {
    return 1.0 // deal with roundoff errors
  }
  if denom2 <= 0.0 {
    return 1.0
  }
  result = result / (denom1 * denom2).sqrt()
  if absolute {
    result = result.abs()
  }
  1.0 - result
}

///|
fn bc_correlation(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  bc_correlation_core(
    n,
    data1,
    data2,
    mask1,
    mask2,
    weight,
    index1,
    index2,
    transpose,
    absolute=false,
    uncentered=false,
  )
}

///|
fn bc_acorrelation(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  bc_correlation_core(
    n,
    data1,
    data2,
    mask1,
    mask2,
    weight,
    index1,
    index2,
    transpose,
    absolute=true,
    uncentered=false,
  )
}

///|
fn bc_ucorrelation(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  bc_correlation_core(
    n,
    data1,
    data2,
    mask1,
    mask2,
    weight,
    index1,
    index2,
    transpose,
    absolute=false,
    uncentered=true,
  )
}

///|
fn bc_uacorrelation(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  bc_correlation_core(
    n,
    data1,
    data2,
    mask1,
    mask2,
    weight,
    index1,
    index2,
    transpose,
    absolute=true,
    uncentered=true,
  )
}

///|
/// Spearman rank correlation distance, 1 - r_s (port of spearman).
fn bc_spearman(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  let mut m = 0
  let tdata1 : Array[Double] = Array::new()
  let tdata2 : Array[Double] = Array::new()
  if transpose == 0 {
    for i = 0; i < n; i = i + 1 {
      if mask1[index1][i] != 0 && mask2[index2][i] != 0 {
        tdata1.push(data1[index1][i])
        tdata2.push(data2[index2][i])
        m = m + 1
      }
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      if mask1[i][index1] != 0 && mask2[i][index2] != 0 {
        tdata1.push(data1[i][index1])
        tdata2.push(data2[i][index2])
        m = m + 1
      }
    }
  }
  if m == 0 {
    return 0.0
  }
  let rank1 = bc_getrank(m, tdata1, weight)
  let rank2 = bc_getrank(m, tdata2, weight)
  let mut result = 0.0
  let mut denom1 = 0.0
  let mut denom2 = 0.0
  let mut sum1 = 0.0
  let mut sum2 = 0.0
  let mut totalweight = 0.0
  for i = 0; i < m; i = i + 1 {
    let term1 = rank1[i]
    let term2 = rank2[i]
    let w = weight[i]
    sum1 = sum1 + term1 * w
    sum2 = sum2 + term2 * w
    result = result + term1 * term2 * w
    denom1 = denom1 + term1 * term1 * w
    denom2 = denom2 + term2 * term2 * w
    totalweight = totalweight + w
  }
  if totalweight == 0.0 {
    return 0.0
  }
  result = result - sum1 * sum2 / totalweight
  denom1 = denom1 - sum1 * sum1 / totalweight
  denom2 = denom2 - sum2 * sum2 / totalweight
  if denom1 <= 0.0 {
    return 1.0
  }
  if denom2 <= 0.0 {
    return 1.0
  }
  result = result / (denom1 * denom2).sqrt()
  1.0 - result
}

///|
/// Kendall's tau distance, 1 - tau (port of kendall).
fn bc_kendall(
  n : Int,
  data1 : Array[Array[Double]],
  data2 : Array[Array[Double]],
  mask1 : Array[Array[Int]],
  mask2 : Array[Array[Int]],
  weight : Array[Double],
  index1 : Int,
  index2 : Int,
  transpose : Int,
) -> Double {
  let mut con = 0.0
  let mut dis = 0.0
  let mut exx = 0.0
  let mut exy = 0.0
  let mut flag = false
  if transpose == 0 {
    for i = 0; i < n; i = i + 1 {
      if mask1[index1][i] != 0 && mask2[index2][i] != 0 {
        for j = 0; j < i; j = j + 1 {
          if mask1[index1][j] != 0 && mask2[index2][j] != 0 {
            let x1 = data1[index1][i]
            let x2 = data1[index1][j]
            let y1 = data2[index2][i]
            let y2 = data2[index2][j]
            let w = weight[i] * weight[j]
            if x1 < x2 && y1 < y2 {
              con = con + w
            } else if x1 > x2 && y1 > y2 {
              con = con + w
            } else if x1 < x2 && y1 > y2 {
              dis = dis + w
            } else if x1 > x2 && y1 < y2 {
              dis = dis + w
            } else if x1 == x2 && y1 != y2 {
              exx = exx + w
            } else if x1 != x2 && y1 == y2 {
              exy = exy + w
            }
            flag = true
          }
        }
      }
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      if mask1[i][index1] != 0 && mask2[i][index2] != 0 {
        for j = 0; j < i; j = j + 1 {
          if mask1[j][index1] != 0 && mask2[j][index2] != 0 {
            let x1 = data1[i][index1]
            let x2 = data1[j][index1]
            let y1 = data2[i][index2]
            let y2 = data2[j][index2]
            let w = weight[i] * weight[j]
            if x1 < x2 && y1 < y2 {
              con = con + w
            } else if x1 > x2 && y1 > y2 {
              con = con + w
            } else if x1 < x2 && y1 > y2 {
              dis = dis + w
            } else if x1 > x2 && y1 < y2 {
              dis = dis + w
            } else if x1 == x2 && y1 != y2 {
              exx = exx + w
            } else if x1 != x2 && y1 == y2 {
              exy = exy + w
            }
            flag = true
          }
        }
      }
    }
  }
  if !flag {
    return 0.0
  }
  let denomx = con + dis + exx
  let denomy = con + dis + exy
  if denomx == 0.0 {
    return 1.0
  }
  if denomy == 0.0 {
    return 1.0
  }
  let tau = (con - dis) / (denomx * denomy).sqrt()
  1.0 - tau
}

///|
priv type BCMetric = (
  Int,
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Array[Int]],
  Array[Array[Int]],
  Array[Double],
  Int,
  Int,
  Int,
) -> Double

///|
/// Select the distance function indicated by dist (port of setmetric):
/// 'e' Euclidean, 'b' City Block, 'c' Pearson, 'a' absolute Pearson,
/// 'u' uncentered Pearson, 'x' absolute uncentered Pearson,
/// 's' Spearman, 'k' Kendall. Other values default to Euclidean.
fn bc_setmetric(dist : String) -> BCMetric {
  match dist {
    "b" => bc_cityblock
    "c" => bc_correlation
    "a" => bc_acorrelation
    "u" => bc_ucorrelation
    "x" => bc_uacorrelation
    "s" => bc_spearman
    "k" => bc_kendall
    _ => bc_euclid
  }
}

// ---------------------------------------------------------------------------
// Distance matrix
// ---------------------------------------------------------------------------

///|
/// Calculate the distance matrix between the rows (transpose=false) or
/// columns (transpose=true) of the data matrix. As the matrix is symmetric
/// with zeros on the diagonal, only the lower triangular half is stored:
/// element i of the result contains i values (matrix[i][j], j < i).
pub fn distancematrix(
  data : Array[Array[Double]],
  mask? : Array[Array[Int]],
  weight? : Array[Double],
  transpose? : Bool,
  dist? : String,
) -> Array[Array[Double]] {
  let tr = transpose.unwrap_or(false)
  let d = dist.unwrap_or("e")
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let n = if tr { ncolumns } else { nrows }
  let ndata = if tr { nrows } else { ncolumns }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let w = weight.unwrap_or(bc_ones_weight(ndata))
  let metric = bc_setmetric(d)
  let matrix : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    matrix.push(row)
  }
  for i = 1; i < n; i = i + 1 {
    for j = 0; j < i; j = j + 1 {
      matrix[i].push(
        metric(ndata, data, data, m, m, w, i, j, if tr { 1 } else { 0 }),
      )
    }
  }
  matrix
}

///|
/// Ragged lower-triangular matrix filled with zeros.
fn bc_zero_distmatrix(n : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    matrix.push(row)
  }
  for i = 1; i < n; i = i + 1 {
    for _j = 0; _j < i; _j = _j + 1 {
      matrix[i].push(0.0)
    }
  }
  matrix
}

///|
fn bc_ones_mask(nrows : Int, ncolumns : Int) -> Array[Array[Int]] {
  let mask : Array[Array[Int]] = Array::new()
  for _i = 0; _i < nrows; _i = _i + 1 {
    let row : Array[Int] = Array::new()
    for _j = 0; _j < ncolumns; _j = _j + 1 {
      row.push(1)
    }
    mask.push(row)
  }
  mask
}

///|
fn bc_ones_weight(n : Int) -> Array[Double] {
  Array::make(n, 1.0)
}

///|
/// Deep copy of a ragged lower-triangular distance matrix.
fn bc_copy_distmatrix(
  distmatrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let copy : Array[Array[Double]] = Array::new()
  for i = 0; i < distmatrix.length(); i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < distmatrix[i].length(); j = j + 1 {
      row.push(distmatrix[i][j])
    }
    copy.push(row)
  }
  copy
}

///|
/// Find the shortest distance between pairs of elements
/// (port of find_closest_pair); returns (distance, ip, jp).
fn bc_find_closest_pair(
  n : Int,
  distmatrix : Array[Array[Double]],
) -> (Double, Int, Int) {
  let mut distance = distmatrix[1][0]
  let mut ip = 1
  let mut jp = 0
  for i = 1; i < n; i = i + 1 {
    for j = 0; j < i; j = j + 1 {
      let temp = distmatrix[i][j]
      if temp < distance {
        distance = temp
        ip = i
        jp = j
      }
    }
  }
  (distance, ip, jp)
}

// ---------------------------------------------------------------------------
// Hierarchical clustering
// ---------------------------------------------------------------------------

///|
/// One node in a hierarchical clustering tree. Elements are numbered
/// 0..nelements-1; nodes are numbered -1..-(nelements-1).
pub struct ClusterNode {
  mut left : Int
  mut right : Int
  mut distance : Double
}

///|
/// Hierarchical clustering tree consisting of ClusterNode elements.
pub struct ClusterTree {
  nodes : Array[ClusterNode]
  nelements : Int
}

///|
/// Number of nodes in the tree (nelements - 1 for a complete tree).
pub fn ClusterTree::length(self : ClusterTree) -> Int {
  self.nodes.length()
}

///|
/// Access node i of the tree.
pub fn ClusterTree::node(self : ClusterTree, i : Int) -> ClusterNode {
  self.nodes[i]
}

///|
/// All nodes of the tree.
pub fn ClusterTree::all_nodes(self : ClusterTree) -> Array[ClusterNode] {
  self.nodes
}

///|
/// Divide the elements in the hierarchical clustering tree into clusters
/// and return the cluster number to which each element was assigned
/// (port of cuttree). Clusters are numbered in the left-to-right order in
/// which they appear in the tree.
pub fn ClusterTree::cut(self : ClusterTree, nclusters? : Int) -> Array[Int] {
  let n = nclusters.unwrap_or(self.nelements)
  let k = if n > self.nelements { self.nelements } else { n }
  bc_cuttree(self.nelements, self.nodes, k)
}

///|
/// Sort the tree by switching the left and right subnodes such that the
/// elements in left-to-right order tend to have increasing order values.
/// Returns the indices of the elements in left-to-right order in the
/// dendrogram (port of sorttree). This mutates the tree.
pub fn ClusterTree::sort(
  self : ClusterTree,
  order? : Array[Double],
) -> Array[Int] {
  let nnodes = self.nodes.length()
  let nelements = nnodes + 1
  let ord = order.unwrap_or(Array::make(nelements, 1.0))
  let indices : Array[Int] = Array::make(nelements, 0)
  bc_sorttree(nnodes, self.nodes, ord, indices)
  indices
}

///|
/// Port of cuttree in cluster.c.
fn bc_cuttree(
  nelements : Int,
  tree : Array[ClusterNode],
  nclusters : Int,
) -> Array[Int] {
  let clusterid : Array[Int] = Array::make(nelements, 0)
  if nclusters == 1 {
    for i = 0; i < nelements; i = i + 1 {
      clusterid[i] = 0
    }
    return clusterid
  }
  if nclusters >= nelements {
    for i = 0; i < nelements; i = i + 1 {
      clusterid[i] = i
    }
    return clusterid
  }
  let parents : Array[Int] = Array::make(nelements - 1, 0)
  let mut i = -nelements + 1 // top node
  let mut k = -1
  let mut previous = nelements
  let n = nelements - nclusters // number of nodes to join
  while true {
    if i >= 0 {
      clusterid[i] = k
      let j = i
      i = previous
      previous = j
    } else {
      let j = -i - 1
      if previous == tree[j].left {
        previous = i
        i = tree[j].right
        if j >= n && (i >= 0 || -i - 1 < n) {
          k = k + 1
        }
      } else if previous == tree[j].right {
        previous = i
        i = parents[j]
        if i == nelements {
          break
        }
      } else {
        parents[j] = previous
        previous = i
        i = tree[j].left
        if j >= n && (i >= 0 || -i - 1 < n) {
          k = k + 1
        }
      }
    }
  }
  clusterid
}

///|
/// Port of sorttree in cluster.c.
fn bc_sorttree(
  nnodes : Int,
  tree : Array[ClusterNode],
  order : Array[Double],
  indices : Array[Int],
) -> Unit {
  let nodecounts : Array[Int] = Array::make(nnodes, 0)
  let nodeorder : Array[Double] = Array::make(nnodes, 0.0)
  for i = 0; i < nnodes; i = i + 1 {
    let i1 = tree[i].left
    let i2 = tree[i].right
    // i1 and i2 are the elements that are to be joined
    let order1 = if i1 < 0 {
      let index = -i1 - 1
      nodeorder[index]
    } else {
      order[i1]
    }
    let counts1 = if i1 < 0 { nodecounts[-i1 - 1] } else { 1 }
    let order2 = if i2 < 0 {
      let index = -i2 - 1
      nodeorder[index]
    } else {
      order[i2]
    }
    let counts2 = if i2 < 0 { nodecounts[-i2 - 1] } else { 1 }
    if order1 > order2 {
      tree[i].left = i2
      tree[i].right = i1
    }
    nodecounts[i] = counts1 + counts2
    nodeorder[i] = (counts1.to_double() * order1 + counts2.to_double() * order2) /
      (counts1 + counts2).to_double()
  }
  let mut i = nnodes - 1
  nodecounts[i] = 0
  while i >= 0 {
    let i1 = tree[i].left
    let i2 = tree[i].right
    let counts1 = if i1 < 0 { nodecounts[-i1 - 1] } else { 1 }
    let mut index = nodecounts[i]
    if i1 >= 0 {
      indices[index] = i1
    } else {
      nodecounts[-i1 - 1] = index
    }
    index = index + counts1
    if i2 >= 0 {
      indices[index] = i2
    } else {
      nodecounts[-i2 - 1] = index
    }
    i = i - 1
  }
}

///|
/// Pairwise single-linkage hierarchical clustering (SLINK algorithm,
/// port of pslcluster in cluster.c, distance-matrix path).
fn bc_pslcluster(
  nelements : Int,
  distmatrix : Array[Array[Double]],
) -> Array[ClusterNode] {
  let nnodes = nelements - 1
  let temp : Array[Double] = Array::make(nnodes, 0.0)
  let vector : Array[Int] = Array::make(nnodes, 0)
  let result : Array[ClusterNode] = Array::new()
  for _i = 0; _i < nelements; _i = _i + 1 {
    result.push({ left: 0, right: 0, distance: bc_DBL_MAX, })
  }
  for i = 0; i < nnodes; i = i + 1 {
    vector[i] = i
  }
  for i = 0; i < nelements; i = i + 1 {
    result[i].distance = bc_DBL_MAX
    for j = 0; j < i; j = j + 1 {
      temp[j] = distmatrix[i][j]
    }
    for j = 0; j < i; j = j + 1 {
      let k = vector[j]
      if result[j].distance >= temp[j] {
        if result[j].distance < temp[k] {
          temp[k] = result[j].distance
        }
        result[j].distance = temp[j]
        vector[j] = i
      } else if temp[j] < temp[k] {
        temp[k] = temp[j]
      }
    }
    for j = 0; j < i; j = j + 1 {
      if result[j].distance >= result[vector[j]].distance {
        vector[j] = i
      }
    }
  }
  // Sort the nodes by distance and link them into a tree
  for i = 0; i < nnodes; i = i + 1 {
    result[i].left = i
  }
  let sorted : Array[ClusterNode] = Array::new()
  for i = 0; i < nnodes; i = i + 1 {
    sorted.push(result[i])
  }
  sorted.sort_by(fn(a : ClusterNode, b : ClusterNode) -> Int {
    if a.distance < b.distance {
      -1
    } else if a.distance > b.distance {
      1
    } else if a.left < b.left {
      -1
    } else if a.left > b.left {
      1
    } else {
      0
    }
  })
  let index : Array[Int] = Array::make(nelements, 0)
  for i = 0; i < nelements; i = i + 1 {
    index[i] = i
  }
  for i = 0; i < nnodes; i = i + 1 {
    let j = sorted[i].left
    let k = vector[j]
    sorted[i].left = index[j]
    sorted[i].right = index[k]
    index[k] = -i - 1
  }
  // Trim to nnodes entries
  let trimmed : Array[ClusterNode] = Array::new()
  for i = 0; i < nnodes; i = i + 1 {
    trimmed.push(sorted[i])
  }
  trimmed
}

///|
/// Pairwise maximum- (complete-) linkage clustering (port of pmlcluster).
/// Note: this function modifies the contents of distmatrix.
fn bc_pmlcluster(
  nelements : Int,
  distmatrix : Array[Array[Double]],
) -> Array[ClusterNode] {
  let clusterid : Array[Int] = Array::make(nelements, 0)
  let result : Array[ClusterNode] = Array::new()
  for _i = 0; _i < nelements - 1; _i = _i + 1 {
    result.push({ left: 0, right: 0, distance: 0.0, })
  }
  // Setup a list specifying to which cluster an element belongs
  for j = 0; j < nelements; j = j + 1 {
    clusterid[j] = j
  }
  let mut n = nelements
  while n > 1 {
    let (dist, is_idx0, js0) = bc_find_closest_pair(n, distmatrix)
    let is_idx = is_idx0
    let js = js0
    result[nelements - n].distance = dist
    // Fix the distances
    for j = 0; j < js; j = j + 1 {
      distmatrix[js][j] = bc_max(distmatrix[is_idx][j], distmatrix[js][j])
    }
    for j = js + 1; j < is_idx; j = j + 1 {
      distmatrix[j][js] = bc_max(distmatrix[is_idx][j], distmatrix[j][js])
    }
    for j = is_idx + 1; j < n; j = j + 1 {
      distmatrix[j][js] = bc_max(distmatrix[j][is_idx], distmatrix[j][js])
    }
    for j = 0; j < is_idx; j = j + 1 {
      distmatrix[is_idx][j] = distmatrix[n - 1][j]
    }
    for j = is_idx + 1; j < n - 1; j = j + 1 {
      distmatrix[j][is_idx] = distmatrix[n - 1][j]
    }
    // Update clusterids
    result[nelements - n].left = clusterid[is_idx]
    result[nelements - n].right = clusterid[js]
    clusterid[js] = n - nelements - 1
    clusterid[is_idx] = clusterid[n - 1]
    n = n - 1
  }
  result
}

///|
/// Pairwise average-linkage clustering (port of palcluster).
/// Note: this function modifies the contents of distmatrix.
fn bc_palcluster(
  nelements : Int,
  distmatrix : Array[Array[Double]],
) -> Array[ClusterNode] {
  let clusterid : Array[Int] = Array::make(nelements, 0)
  let number : Array[Int] = Array::make(nelements, 0)
  let result : Array[ClusterNode] = Array::new()
  for _i = 0; _i < nelements - 1; _i = _i + 1 {
    result.push({ left: 0, right: 0, distance: 0.0, })
  }
  // Setup a list specifying to which cluster an element belongs, and keep
  // track of the number of elements in each cluster.
  for j = 0; j < nelements; j = j + 1 {
    number[j] = 1
    clusterid[j] = j
  }
  let mut n = nelements
  while n > 1 {
    let (dist, is_idx0, js0) = bc_find_closest_pair(n, distmatrix)
    let is_idx = is_idx0
    let js = js0
    result[nelements - n].distance = dist
    // Save result
    result[nelements - n].left = clusterid[is_idx]
    result[nelements - n].right = clusterid[js]
    // Fix the distances
    let sum = number[is_idx] + number[js]
    for j = 0; j < js; j = j + 1 {
      distmatrix[js][j] = (
          distmatrix[is_idx][j] * number[is_idx].to_double() +
          distmatrix[js][j] * number[js].to_double()
        ) /
        sum.to_double()
    }
    for j = js + 1; j < is_idx; j = j + 1 {
      distmatrix[j][js] = (
          distmatrix[is_idx][j] * number[is_idx].to_double() +
          distmatrix[j][js] * number[js].to_double()
        ) /
        sum.to_double()
    }
    for j = is_idx + 1; j < n; j = j + 1 {
      distmatrix[j][js] = (
          distmatrix[j][is_idx] * number[is_idx].to_double() +
          distmatrix[j][js] * number[js].to_double()
        ) /
        sum.to_double()
    }
    for j = 0; j < is_idx; j = j + 1 {
      distmatrix[is_idx][j] = distmatrix[n - 1][j]
    }
    for j = is_idx + 1; j < n - 1; j = j + 1 {
      distmatrix[j][is_idx] = distmatrix[n - 1][j]
    }
    // Update number of elements in the clusters
    number[js] = sum
    number[is_idx] = number[n - 1]
    // Update clusterids
    clusterid[js] = n - nelements - 1
    clusterid[is_idx] = clusterid[n - 1]
    n = n - 1
  }
  result
}

///|
/// Pairwise centroid-linkage clustering (port of pclcluster).
/// Note: this function modifies the contents of distmatrix.
fn bc_pclcluster(
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  weight : Array[Double],
  distmatrix : Array[Array[Double]],
  dist : String,
  transpose : Int,
) -> Array[ClusterNode] {
  let nelements = if transpose == 0 { nrows } else { ncolumns }
  let ndata = if transpose == 0 { ncolumns } else { nrows }
  let nnodes = nelements - 1
  let metric = bc_setmetric(dist)
  let distid : Array[Int] = Array::make(nelements, 0)
  let result : Array[ClusterNode] = Array::new()
  for _i = 0; _i < nnodes; _i = _i + 1 {
    result.push({ left: 0, right: 0, distance: 0.0, })
  }
  // Storage for node data (a private copy of data and mask)
  let newdata : Array[Array[Double]] = Array::new()
  let newmask : Array[Array[Int]] = Array::new()
  for i = 0; i < nelements; i = i + 1 {
    let rowd : Array[Double] = Array::new()
    let rowm : Array[Int] = Array::new()
    for j = 0; j < ndata; j = j + 1 {
      if transpose == 0 {
        rowd.push(data[i][j])
        rowm.push(mask[i][j])
      } else {
        rowd.push(data[j][i])
        rowm.push(mask[j][i])
      }
    }
    newdata.push(rowd)
    newmask.push(rowm)
  }
  for i = 0; i < nelements; i = i + 1 {
    distid[i] = i
  }
  for inode = 0; inode < nnodes; inode = inode + 1 {
    // Find the pair with the shortest distance
    let (d, is_idx, js) = bc_find_closest_pair(nelements - inode, distmatrix)
    result[inode].distance = d
    result[inode].left = distid[js]
    result[inode].right = distid[is_idx]
    // Make node js the new node
    for i = 0; i < ndata; i = i + 1 {
      newdata[js][i] = newdata[js][i] * newmask[js][i].to_double() +
        newdata[is_idx][i] * newmask[is_idx][i].to_double()
      newmask[js][i] = newmask[js][i] + newmask[is_idx][i]
      if newmask[js][i] > 0 {
        newdata[js][i] = newdata[js][i] / newmask[js][i].to_double()
      }
    }
    // Move the last active row into position is_idx
    newdata[is_idx] = newdata[nnodes - inode]
    newmask[is_idx] = newmask[nnodes - inode]
    // Fix the distances
    distid[is_idx] = distid[nnodes - inode]
    for i = 0; i < is_idx; i = i + 1 {
      distmatrix[is_idx][i] = distmatrix[nnodes - inode][i]
    }
    for i = is_idx + 1; i < nnodes - inode; i = i + 1 {
      distmatrix[i][is_idx] = distmatrix[nnodes - inode][i]
    }
    distid[js] = -inode - 1
    for i = 0; i < js; i = i + 1 {
      distmatrix[js][i] = metric(
        ndata, newdata, newdata, newmask, newmask, weight, js, i, 0,
      )
    }
    for i = js + 1; i < nnodes - inode; i = i + 1 {
      distmatrix[i][js] = metric(
        ndata, newdata, newdata, newmask, newmask, weight, js, i, 0,
      )
    }
  }
  result
}

///|
/// Perform hierarchical clustering on the rows (transpose=false) or the
/// columns (transpose=true) of the data matrix, and return a ClusterTree.
///
/// Methods:
///   - method = "s": single pairwise linkage
///   - method = "m": complete (maximum) pairwise linkage (default)
///   - method = "a": average pairwise linkage
///   - method = "c": centroid linkage
///
/// Distance measures (dist):
///   - "e" Euclidean, "b" City Block, "c" Pearson, "a" absolute Pearson,
///     "u" uncentered Pearson, "x" absolute uncentered Pearson,
///     "s" Spearman, "k" Kendall.
pub fn treecluster(
  data : Array[Array[Double]],
  mask? : Array[Array[Int]],
  weight? : Array[Double],
  transpose? : Bool,
  method? : String,
  dist? : String,
) -> ClusterTree {
  let tr = transpose.unwrap_or(false)
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let nelements = if tr { ncolumns } else { nrows }
  if nelements < 2 {
    return { nodes: Array::new(), nelements, }
  }
  let ndata = if tr { nrows } else { ncolumns }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let w = weight.unwrap_or(bc_ones_weight(ndata))
  let d = dist.unwrap_or("e")
  let meth = method.unwrap_or("m")
  // Calculate the distance matrix (ragged lower triangle)
  let distmatrix : Array[Array[Double]] = Array::new()
  for _i = 0; _i < nelements; _i = _i + 1 {
    distmatrix.push(Array::new())
  }
  let metric = bc_setmetric(d)
  for i = 1; i < nelements; i = i + 1 {
    for j = 0; j < i; j = j + 1 {
      distmatrix[i].push(
        metric(ndata, data, data, m, m, w, i, j, if tr { 1 } else { 0 }),
      )
    }
  }
  let nodes = match meth {
    "s" => bc_pslcluster(nelements, distmatrix)
    "a" => bc_palcluster(nelements, distmatrix)
    "c" =>
      bc_pclcluster(
        nrows,
        ncolumns,
        data,
        m,
        w,
        distmatrix,
        d,
        if tr {
          1
        } else {
          0
        },
      )
    "m" => bc_pmlcluster(nelements, distmatrix)
    _ => Array::new()
  }
  { nodes, nelements, }
}

///|
/// Perform hierarchical clustering from a pre-computed distance matrix
/// (a ragged lower-triangular matrix, as returned by distancematrix).
/// Only single ("s"), complete ("m", default) and average ("a") linkage
/// are available; centroid linkage requires the original data.
pub fn treecluster_from_distance(
  distance : Array[Array[Double]],
  method? : String,
) -> ClusterTree {
  let nelements = distance.length()
  if nelements < 2 {
    return { nodes: Array::new(), nelements, }
  }
  let meth = method.unwrap_or("m")
  let distmatrix = bc_copy_distmatrix(distance)
  let nodes = match meth {
    "s" => bc_pslcluster(nelements, distmatrix)
    "a" => bc_palcluster(nelements, distmatrix)
    "m" => bc_pmlcluster(nelements, distmatrix)
    _ => Array::new()
  }
  { nodes, nelements, }
}

// ---------------------------------------------------------------------------
// Cluster centroids (means / medians / medoids)
// ---------------------------------------------------------------------------

///|
/// Result of clustercentroids: the cluster centroids and their mask.
pub struct CentroidsResult {
  cdata : Array[Array[Double]]
  cmask : Array[Array[Int]]
}

///|
pub fn CentroidsResult::cdata(self : CentroidsResult) -> Array[Array[Double]] {
  self.cdata
}

///|
pub fn CentroidsResult::cmask(self : CentroidsResult) -> Array[Array[Int]] {
  self.cmask
}

///|
/// Calculate the cluster centroids as the mean over all elements belonging
/// to each cluster, for each dimension (port of getclustermeans).
fn bc_getclustermeans(
  nclusters : Int,
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  clusterid : Array[Int],
  cdata : Array[Array[Double]],
  cmask : Array[Array[Int]],
  transpose : Int,
) -> Unit {
  if transpose == 0 {
    for i = 0; i < nclusters; i = i + 1 {
      for j = 0; j < ncolumns; j = j + 1 {
        cmask[i][j] = 0
        cdata[i][j] = 0.0
      }
    }
    for k = 0; k < nrows; k = k + 1 {
      let i = clusterid[k]
      for j = 0; j < ncolumns; j = j + 1 {
        if mask[k][j] != 0 {
          cdata[i][j] = cdata[i][j] + data[k][j]
          cmask[i][j] = cmask[i][j] + 1
        }
      }
    }
    for i = 0; i < nclusters; i = i + 1 {
      for j = 0; j < ncolumns; j = j + 1 {
        if cmask[i][j] > 0 {
          cdata[i][j] = cdata[i][j] / cmask[i][j].to_double()
          cmask[i][j] = 1
        }
      }
    }
  } else {
    for i = 0; i < nrows; i = i + 1 {
      for j = 0; j < nclusters; j = j + 1 {
        cdata[i][j] = 0.0
        cmask[i][j] = 0
      }
    }
    for k = 0; k < ncolumns; k = k + 1 {
      let i = clusterid[k]
      for j = 0; j < nrows; j = j + 1 {
        if mask[j][k] != 0 {
          cdata[j][i] = cdata[j][i] + data[j][k]
          cmask[j][i] = cmask[j][i] + 1
        }
      }
    }
    for i = 0; i < nrows; i = i + 1 {
      for j = 0; j < nclusters; j = j + 1 {
        if cmask[i][j] > 0 {
          cdata[i][j] = cdata[i][j] / cmask[i][j].to_double()
          cmask[i][j] = 1
        }
      }
    }
  }
}

///|
/// Calculate the cluster centroids as the median over all elements
/// belonging to each cluster, for each dimension (port of getclustermedians).
fn bc_getclustermedians(
  nclusters : Int,
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  clusterid : Array[Int],
  cdata : Array[Array[Double]],
  cmask : Array[Array[Int]],
  transpose : Int,
  cache : Array[Double],
) -> Unit {
  if transpose == 0 {
    for i = 0; i < nclusters; i = i + 1 {
      for j = 0; j < ncolumns; j = j + 1 {
        let mut count = 0
        for k = 0; k < nrows; k = k + 1 {
          if i == clusterid[k] && mask[k][j] != 0 {
            cache[count] = data[k][j]
            count = count + 1
          }
        }
        if count > 0 {
          cdata[i][j] = bc_median(count, cache)
          cmask[i][j] = 1
        } else {
          cdata[i][j] = 0.0
          cmask[i][j] = 0
        }
      }
    }
  } else {
    for i = 0; i < nclusters; i = i + 1 {
      for j = 0; j < nrows; j = j + 1 {
        let mut count = 0
        for k = 0; k < ncolumns; k = k + 1 {
          if i == clusterid[k] && mask[j][k] != 0 {
            cache[count] = data[j][k]
            count = count + 1
          }
        }
        if count > 0 {
          cdata[j][i] = bc_median(count, cache)
          cmask[j][i] = 1
        } else {
          cdata[j][i] = 0.0
          cmask[j][i] = 0
        }
      }
    }
  }
}

///|
/// Calculate the medoids: the element of each cluster with the smallest sum
/// of distances to the other elements (port of getclustermedoids).
fn bc_getclustermedoids(
  nclusters : Int,
  nelements : Int,
  distmatrix : Array[Array[Double]],
  clusterid : Array[Int],
  centroids : Array[Int],
  errors : Array[Double],
) -> Unit {
  for j = 0; j < nclusters; j = j + 1 {
    errors[j] = bc_DBL_MAX
  }
  for i = 0; i < nelements; i = i + 1 {
    let mut d = 0.0
    let j = clusterid[i]
    for k = 0; k < nelements; k = k + 1 {
      if i == k || clusterid[k] != j {
        continue
      }
      d = d + (if i < k { distmatrix[k][i] } else { distmatrix[i][k] })
      if d > errors[j] {
        break
      }
    }
    if d < errors[j] {
      errors[j] = d
      centroids[j] = i
    }
  }
}

// ---------------------------------------------------------------------------
// k-means / k-median clustering
// ---------------------------------------------------------------------------

///|
/// EM iteration of k-means or k-median clustering (port of kmeans and
/// kmedians in cluster.c). Returns (ifound, error).
fn bc_kmeans(
  nclusters : Int,
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  weight : Array[Double],
  transpose : Int,
  npass : Int,
  dist : String,
  cdata : Array[Array[Double]],
  cmask : Array[Array[Int]],
  clusterid : Array[Int],
  use_medians : Bool,
) -> (Int, Double) {
  let nelements = if transpose == 0 { nrows } else { ncolumns }
  let ndata = if transpose == 0 { ncolumns } else { nrows }
  let metric = bc_setmetric(dist)
  let mut ifound = 1
  let mut error = bc_DBL_MAX
  let tclusterid : Array[Int] = Array::make(nelements, 0)
  let counts : Array[Int] = Array::make(nclusters, 0)
  let mapping : Array[Int] = Array::make(nclusters, -1)
  let saved : Array[Int] = Array::make(nelements, 0)
  let cache : Array[Double] = Array::make(nelements, 0.0)
  let mut ipass = 0
  while true {
    let mut total = bc_DBL_MAX
    let mut counter = 0
    let mut period = 10
    // Perform the EM algorithm.
    // First, randomly assign elements to clusters.
    if npass != 0 {
      bc_randomassign(nclusters, nelements, tclusterid)
    } else {
      for i = 0; i < nelements; i = i + 1 {
        tclusterid[i] = clusterid[i]
      }
    }
    for i = 0; i < nclusters; i = i + 1 {
      counts[i] = 0
    }
    for i = 0; i < nelements; i = i + 1 {
      counts[tclusterid[i]] = counts[tclusterid[i]] + 1
    }
    // Start the loop
    while true {
      let previous = total
      total = 0.0
      if counter % period == 0 {
        // Save the current cluster assignments
        for i = 0; i < nelements; i = i + 1 {
          saved[i] = tclusterid[i]
        }
        if period < 1073741823 {
          period = period * 2
        }
      }
      counter = counter + 1
      // Find the center
      if use_medians {
        bc_getclustermedians(
          nclusters, nrows, ncolumns, data, mask, tclusterid, cdata, cmask, transpose,
          cache,
        )
      } else {
        bc_getclustermeans(
          nclusters, nrows, ncolumns, data, mask, tclusterid, cdata, cmask, transpose,
        )
      }
      for i = 0; i < nelements; i = i + 1 {
        // Calculate the distances
        let k = tclusterid[i]
        if counts[k] == 1 {
          continue
        }
        // No reassignment if that would lead to an empty cluster.
        // Treat the present cluster as a special case.
        let mut distance = metric(
          ndata, data, cdata, mask, cmask, weight, i, k, transpose,
        )
        for j = 0; j < nclusters; j = j + 1 {
          if j == k {
            continue
          }
          let tdistance = metric(
            ndata, data, cdata, mask, cmask, weight, i, j, transpose,
          )
          if tdistance < distance {
            distance = tdistance
            counts[tclusterid[i]] = counts[tclusterid[i]] - 1
            tclusterid[i] = j
            counts[j] = counts[j] + 1
          }
        }
        total = total + distance
      }
      if total >= previous {
        break
      }
      // Check if an identical solution was found
      let mut same = true
      for i = 0; i < nelements; i = i + 1 {
        if saved[i] != tclusterid[i] {
          same = false
          break
        }
      }
      if same {
        break
      }
    }
    if npass <= 1 {
      error = total
      // The clustering solution is the current assignment
      for i = 0; i < nelements; i = i + 1 {
        clusterid[i] = tclusterid[i]
      }
      break
    }
    for i = 0; i < nclusters; i = i + 1 {
      mapping[i] = -1
    }
    let mut broke = false
    for i = 0; i < nelements; i = i + 1 {
      let j = tclusterid[i]
      let k = clusterid[i]
      if mapping[k] == -1 {
        mapping[k] = j
      } else if mapping[k] != j {
        if total < error {
          ifound = 1
          error = total
          for j2 = 0; j2 < nelements; j2 = j2 + 1 {
            clusterid[j2] = tclusterid[j2]
          }
        }
        broke = true
        break
      }
    }
    if !broke {
      ifound = ifound + 1 // identical solution found
    }
    ipass = ipass + 1
    if ipass >= npass {
      break
    }
  }
  (ifound, error)
}

///|
/// Result of kcluster: cluster assignments, within-cluster sum of distances
/// and the number of times the optimal solution was found.
pub struct KClusterResult {
  clusterid : Array[Int]
  error : Double
  ifound : Int
}

///|
pub fn KClusterResult::clusterid(self : KClusterResult) -> Array[Int] {
  self.clusterid
}

///|
pub fn KClusterResult::error(self : KClusterResult) -> Double {
  self.error
}

///|
pub fn KClusterResult::ifound(self : KClusterResult) -> Int {
  self.ifound
}

///|
/// Perform k-means ("a") or k-median ("m") clustering on the rows
/// (transpose=false) or the columns (transpose=true) of the data matrix
/// (port of kcluster in cluster.c).
///
/// If initialid is provided, the EM algorithm is performed once starting
/// from the given initial clustering (fully deterministic). Otherwise the
/// algorithm is run npass times with random initial assignments, and the
/// solution with the lowest within-cluster sum of distances is kept.
pub fn kcluster(
  data : Array[Array[Double]],
  nclusters? : Int,
  mask? : Array[Array[Int]],
  weight? : Array[Double],
  transpose? : Bool,
  npass? : Int,
  method? : String,
  dist? : String,
  initialid? : Array[Int],
) -> KClusterResult {
  let tr = transpose.unwrap_or(false)
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let nelements = if tr { ncolumns } else { nrows }
  let ndata = if tr { nrows } else { ncolumns }
  let k = nclusters.unwrap_or(2)
  let mut passes = npass.unwrap_or(1)
  let meth = method.unwrap_or("a")
  let d = dist.unwrap_or("e")
  let clusterid : Array[Int] = match initialid {
    Some(init) => {
      passes = 0 // deterministic single run from the given clustering
      let arr : Array[Int] = Array::new()
      for i = 0; i < nelements; i = i + 1 {
        arr.push(init[i])
      }
      arr
    }
    None => Array::make(nelements, 0)
  }
  if nelements < k {
    // More clusters asked for than elements available
    return { clusterid, error: 0.0, ifound: 0, }
  }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let w = weight.unwrap_or(bc_ones_weight(ndata))
  // Allocate space to store the centroid data
  let (crow, ccol) = if tr { (ndata, k) } else { (k, ndata) }
  let cdata : Array[Array[Double]] = Array::new()
  let cmask : Array[Array[Int]] = Array::new()
  for _i = 0; _i < crow; _i = _i + 1 {
    cdata.push(Array::make(ccol, 0.0))
    cmask.push(Array::make(ccol, 0))
  }
  let use_medians = meth == "m"
  let (ifound, error) = bc_kmeans(
    k,
    nrows,
    ncolumns,
    data,
    m,
    w,
    if tr {
      1
    } else {
      0
    },
    passes,
    d,
    cdata,
    cmask,
    clusterid,
    use_medians,
  )
  { clusterid, error, ifound, }
}

// ---------------------------------------------------------------------------
// k-medoids clustering
// ---------------------------------------------------------------------------

///|
/// Result of kmedoids: cluster assignments (cluster numbers are the item
/// numbers of the cluster medoids), within-cluster sum of distances, and
/// the number of times the optimal solution was found.
pub struct KMedoidsResult {
  clusterid : Array[Int]
  error : Double
  ifound : Int
}

///|
pub fn KMedoidsResult::clusterid(self : KMedoidsResult) -> Array[Int] {
  self.clusterid
}

///|
pub fn KMedoidsResult::error(self : KMedoidsResult) -> Double {
  self.error
}

///|
pub fn KMedoidsResult::ifound(self : KMedoidsResult) -> Int {
  self.ifound
}

///|
/// Perform k-medoids clustering on the given (ragged lower-triangular)
/// distance matrix (port of kmedoids in cluster.c).
pub fn kmedoids(
  distance : Array[Array[Double]],
  nclusters? : Int,
  npass? : Int,
  initialid? : Array[Int],
) -> KMedoidsResult {
  let nelements = distance.length()
  let k = nclusters.unwrap_or(2)
  let mut passes = npass.unwrap_or(1)
  let clusterid : Array[Int] = match initialid {
    Some(init) => {
      passes = 0
      let arr : Array[Int] = Array::new()
      for i = 0; i < nelements; i = i + 1 {
        arr.push(init[i])
      }
      arr
    }
    None => Array::make(nelements, -1)
  }
  if nelements < k {
    // More clusters asked for than elements available
    return { clusterid, error: 0.0, ifound: 0, }
  }
  let tclusterid : Array[Int] = Array::make(nelements, 0)
  let saved : Array[Int] = Array::make(nelements, 0)
  let centroids : Array[Int] = Array::make(k, 0)
  let errors : Array[Double] = Array::make(k, 0.0)
  let mut ifound = 1
  let mut error = bc_DBL_MAX
  let mut ipass = 0
  while true {
    let mut total = bc_DBL_MAX
    let mut counter = 0
    let mut period = 10
    if passes != 0 {
      bc_randomassign(k, nelements, tclusterid)
    } else {
      for i = 0; i < nelements; i = i + 1 {
        tclusterid[i] = clusterid[i]
      }
    }
    while true {
      let previous = total
      total = 0.0
      if counter % period == 0 {
        for i = 0; i < nelements; i = i + 1 {
          saved[i] = tclusterid[i]
        }
        if period < 1073741823 {
          period = period * 2
        }
      }
      counter = counter + 1
      // Find the center
      bc_getclustermedoids(
        k, nelements, distance, tclusterid, centroids, errors,
      )
      for i = 0; i < nelements; i = i + 1 {
        // Find the closest cluster
        let mut dist = bc_DBL_MAX
        for icluster = 0; icluster < k; icluster = icluster + 1 {
          let j = centroids[icluster]
          if i == j {
            dist = 0.0
            tclusterid[i] = icluster
            break
          }
          let tdistance = if i > j { distance[i][j] } else { distance[j][i] }
          if tdistance < dist {
            dist = tdistance
            tclusterid[i] = icluster
          }
        }
        total = total + dist
      }
      if total >= previous {
        break
      }
      let mut same = true
      for i = 0; i < nelements; i = i + 1 {
        if saved[i] != tclusterid[i] {
          same = false
          break
        }
      }
      if same {
        break
      }
    }
    if passes <= 1 {
      ifound = 1
      error = total
      // Replace by the centroid in each cluster.
      for j = 0; j < nelements; j = j + 1 {
        clusterid[j] = centroids[tclusterid[j]]
      }
      break
    }
    let mut broke = false
    for i = 0; i < nelements; i = i + 1 {
      if clusterid[i] != centroids[tclusterid[i]] {
        if total < error {
          ifound = 1
          error = total
          for j = 0; j < nelements; j = j + 1 {
            clusterid[j] = centroids[tclusterid[j]]
          }
        }
        broke = true
        break
      }
    }
    if !broke {
      ifound = ifound + 1
    }
    ipass = ipass + 1
    if ipass >= passes {
      break
    }
  }
  { clusterid, error, ifound, }
}

// ---------------------------------------------------------------------------
// clustercentroids and clusterdistance
// ---------------------------------------------------------------------------

///|
/// Calculate the centroid of each cluster, given to which cluster each item
/// belongs. The centroid is defined as the mean (method "a", default) or the
/// median (method "m") over all elements belonging to each cluster, for each
/// dimension (port of clustercentroids / getclustercentroids).
///
/// If clusterid is omitted, all items are assumed to belong to cluster 0.
/// If transpose is false, cdata has dimensions nclusters x ncolumns;
/// otherwise nrows x nclusters.
pub fn clustercentroids(
  data : Array[Array[Double]],
  mask? : Array[Array[Int]],
  clusterid? : Array[Int],
  method? : String,
  transpose? : Bool,
) -> CentroidsResult {
  let tr = transpose.unwrap_or(false)
  let transpose_i = if tr { 1 } else { 0 }
  let meth = method.unwrap_or("a")
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let nelements = if tr { ncolumns } else { nrows }
  let ids : Array[Int] = match clusterid {
    Some(ids) => ids
    None => Array::make(nelements, 0)
  }
  let mut nclusters = 0
  for i = 0; i < nelements; i = i + 1 {
    if ids[i] + 1 > nclusters {
      nclusters = ids[i] + 1
    }
  }
  let cdata : Array[Array[Double]] = Array::new()
  let cmask : Array[Array[Int]] = Array::new()
  if transpose_i == 0 {
    for _i = 0; _i < nclusters; _i = _i + 1 {
      cdata.push(Array::make(ncolumns, 0.0))
      cmask.push(Array::make(ncolumns, 0))
    }
  } else {
    for _i = 0; _i < nrows; _i = _i + 1 {
      cdata.push(Array::make(nclusters, 0.0))
      cmask.push(Array::make(nclusters, 0))
    }
  }
  match meth {
    "m" => {
      let cache : Array[Double] = Array::make(nelements, 0.0)
      bc_getclustermedians(
        nclusters, nrows, ncolumns, data, m, ids, cdata, cmask, transpose_i, cache,
      )
    }
    "a" =>
      bc_getclustermeans(
        nclusters, nrows, ncolumns, data, m, ids, cdata, cmask, transpose_i,
      )
    _ => () // unknown method: return zero-filled centroids
  }
  { cdata, cmask, }
}

///|
/// Calculate the distance between two clusters of genes (transpose=false) or
/// samples (transpose=true), using the measured expression vectors
/// (port of clusterdistance).
///
/// method:
///   "a" (default) : distance between the arithmetic means of the clusters
///   "m"           : distance between the medians of the clusters
///   "s"           : smallest pairwise distance between members
///   "x"           : largest pairwise distance between members
///   "v"           : average of the pairwise distances between members
///
/// Returns -1.0 if a cluster is empty or an index is out of bounds, and
/// -2.0 if the method is unknown (mirroring the C Clustering Library).
pub fn clusterdistance(
  data : Array[Array[Double]],
  mask? : Array[Array[Int]],
  weight? : Array[Double],
  index1? : Array[Int],
  index2? : Array[Int],
  method? : String,
  dist? : String,
  transpose? : Bool,
) -> Double {
  let tr = transpose.unwrap_or(false)
  let transpose_i = if tr { 1 } else { 0 }
  let meth = method.unwrap_or("a")
  let d = dist.unwrap_or("e")
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let ndata = if tr { nrows } else { ncolumns }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let w = weight.unwrap_or(bc_ones_weight(ndata))
  let metric = bc_setmetric(d)
  let idx1 = index1.unwrap_or([0])
  let idx2 = index2.unwrap_or([0])
  let n1 = idx1.length()
  let n2 = idx2.length()
  // If one or both clusters are empty, return -1
  if n1 < 1 || n2 < 1 {
    return -1.0
  }
  // Check the indices
  let limit = if tr { ncolumns } else { nrows }
  for i = 0; i < n1; i = i + 1 {
    if idx1[i] < 0 || idx1[i] >= limit {
      return -1.0
    }
  }
  for i = 0; i < n2; i = i + 1 {
    if idx2[i] < 0 || idx2[i] >= limit {
      return -1.0
    }
  }
  match meth {
    "a" => {
      // Distance between the arithmetic means of the two clusters
      let cdata : Array[Array[Double]] = Array::new()
      let cmask : Array[Array[Int]] = Array::new()
      let count : Array[Array[Int]] = Array::new()
      if transpose_i == 0 {
        for _i = 0; _i < 2; _i = _i + 1 {
          cdata.push(Array::make(ncolumns, 0.0))
          cmask.push(Array::make(ncolumns, 0))
          count.push(Array::make(ncolumns, 0))
        }
        for i = 0; i < n1; i = i + 1 {
          let k = idx1[i]
          for j = 0; j < ncolumns; j = j + 1 {
            if m[k][j] != 0 {
              cdata[0][j] = cdata[0][j] + data[k][j]
              count[0][j] = count[0][j] + 1
            }
          }
        }
        for i = 0; i < n2; i = i + 1 {
          let k = idx2[i]
          for j = 0; j < ncolumns; j = j + 1 {
            if m[k][j] != 0 {
              cdata[1][j] = cdata[1][j] + data[k][j]
              count[1][j] = count[1][j] + 1
            }
          }
        }
        for i = 0; i < 2; i = i + 1 {
          for j = 0; j < ncolumns; j = j + 1 {
            if count[i][j] > 0 {
              cdata[i][j] = cdata[i][j] / count[i][j].to_double()
              cmask[i][j] = 1
            } else {
              cmask[i][j] = 0
            }
          }
        }
        metric(ncolumns, cdata, cdata, cmask, cmask, w, 0, 1, 0)
      } else {
        for _i = 0; _i < nrows; _i = _i + 1 {
          cdata.push(Array::make(2, 0.0))
          cmask.push(Array::make(2, 0))
          count.push(Array::make(2, 0))
        }
        for i = 0; i < n1; i = i + 1 {
          let k = idx1[i]
          for j = 0; j < nrows; j = j + 1 {
            if m[j][k] != 0 {
              cdata[j][0] = cdata[j][0] + data[j][k]
              count[j][0] = count[j][0] + 1
            }
          }
        }
        for i = 0; i < n2; i = i + 1 {
          let k = idx2[i]
          for j = 0; j < nrows; j = j + 1 {
            if m[j][k] != 0 {
              cdata[j][1] = cdata[j][1] + data[j][k]
              count[j][1] = count[j][1] + 1
            }
          }
        }
        for i = 0; i < nrows; i = i + 1 {
          for j = 0; j < 2; j = j + 1 {
            if count[i][j] > 0 {
              cdata[i][j] = cdata[i][j] / count[i][j].to_double()
              cmask[i][j] = 1
            } else {
              cmask[i][j] = 0
            }
          }
        }
        metric(nrows, cdata, cdata, cmask, cmask, w, 0, 1, 1)
      }
    }
    "m" => {
      // Distance between the medians of the two clusters
      let cdata : Array[Array[Double]] = Array::new()
      let cmask : Array[Array[Int]] = Array::new()
      if transpose_i == 0 {
        for _i = 0; _i < 2; _i = _i + 1 {
          cdata.push(Array::make(ncolumns, 0.0))
          cmask.push(Array::make(ncolumns, 0))
        }
        let temp : Array[Double] = Array::make(nrows, 0.0)
        for j = 0; j < ncolumns; j = j + 1 {
          let mut count = 0
          for k = 0; k < n1; k = k + 1 {
            let i = idx1[k]
            if m[i][j] != 0 {
              temp[count] = data[i][j]
              count = count + 1
            }
          }
          if count > 0 {
            cdata[0][j] = bc_median(count, temp)
            cmask[0][j] = 1
          } else {
            cdata[0][j] = 0.0
            cmask[0][j] = 0
          }
        }
        for j = 0; j < ncolumns; j = j + 1 {
          let mut count = 0
          for k = 0; k < n2; k = k + 1 {
            let i = idx2[k]
            if m[i][j] != 0 {
              temp[count] = data[i][j]
              count = count + 1
            }
          }
          if count > 0 {
            cdata[1][j] = bc_median(count, temp)
            cmask[1][j] = 1
          } else {
            cdata[1][j] = 0.0
            cmask[1][j] = 0
          }
        }
        metric(ncolumns, cdata, cdata, cmask, cmask, w, 0, 1, 0)
      } else {
        for _i = 0; _i < nrows; _i = _i + 1 {
          cdata.push(Array::make(2, 0.0))
          cmask.push(Array::make(2, 0))
        }
        let temp : Array[Double] = Array::make(ncolumns, 0.0)
        for j = 0; j < nrows; j = j + 1 {
          let mut count = 0
          for k = 0; k < n1; k = k + 1 {
            let i = idx1[k]
            if m[j][i] != 0 {
              temp[count] = data[j][i]
              count = count + 1
            }
          }
          if count > 0 {
            cdata[j][0] = bc_median(count, temp)
            cmask[j][0] = 1
          } else {
            cdata[j][0] = 0.0
            cmask[j][0] = 0
          }
        }
        for j = 0; j < nrows; j = j + 1 {
          let mut count = 0
          for k = 0; k < n2; k = k + 1 {
            let i = idx2[k]
            if m[j][i] != 0 {
              temp[count] = data[j][i]
              count = count + 1
            }
          }
          if count > 0 {
            cdata[j][1] = bc_median(count, temp)
            cmask[j][1] = 1
          } else {
            cdata[j][1] = 0.0
            cmask[j][1] = 0
          }
        }
        metric(nrows, cdata, cdata, cmask, cmask, w, 0, 1, 1)
      }
    }
    "s" => {
      // Smallest pairwise distance between members of the two clusters
      let mut mindistance = bc_DBL_MAX
      for i1 = 0; i1 < n1; i1 = i1 + 1 {
        for i2 = 0; i2 < n2; i2 = i2 + 1 {
          let j1 = idx1[i1]
          let j2 = idx2[i2]
          let distance = metric(ndata, data, data, m, m, w, j1, j2, transpose_i)
          if distance < mindistance {
            mindistance = distance
          }
        }
      }
      mindistance
    }
    "x" => {
      // Largest pairwise distance between members of the two clusters
      let mut maxdistance = 0.0
      for i1 = 0; i1 < n1; i1 = i1 + 1 {
        for i2 = 0; i2 < n2; i2 = i2 + 1 {
          let j1 = idx1[i1]
          let j2 = idx2[i2]
          let distance = metric(ndata, data, data, m, m, w, j1, j2, transpose_i)
          if distance > maxdistance {
            maxdistance = distance
          }
        }
      }
      maxdistance
    }
    "v" => {
      // Average of the pairwise distances between members of the two clusters
      let mut distance = 0.0
      for i1 = 0; i1 < n1; i1 = i1 + 1 {
        for i2 = 0; i2 < n2; i2 = i2 + 1 {
          let j1 = idx1[i1]
          let j2 = idx2[i2]
          distance = distance +
            metric(ndata, data, data, m, m, w, j1, j2, transpose_i)
        }
      }
      distance / (n1 * n2).to_double()
    }
    _ => -2.0
  }
}

// ---------------------------------------------------------------------------
// Self-organizing maps (somcluster)
// ---------------------------------------------------------------------------

///|
/// Result of somcluster: the grid coordinates of the cell to which each item
/// was assigned, and the gene expression data vector for each node (cell) in
/// the 2D grid.
pub struct SomResult {
  clusterid : Array[(Int, Int)]
  celldata : Array[Array[Array[Double]]]
}

///|
pub fn SomResult::clusterid(self : SomResult) -> Array[(Int, Int)] {
  self.clusterid
}

///|
pub fn SomResult::celldata(self : SomResult) -> Array[Array[Array[Double]]] {
  self.celldata
}

///|
/// Train the self-organizing map (port of somworker in cluster.c).
fn bc_somworker(
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  weights : Array[Double],
  transpose : Int,
  nxgrid : Int,
  nygrid : Int,
  inittau : Double,
  celldata : Array[Array[Array[Double]]],
  niter : Int,
  dist : String,
) -> Unit {
  let nelements = if transpose == 0 { nrows } else { ncolumns }
  let ndata = if transpose == 0 { ncolumns } else { nrows }
  let metric = bc_setmetric(dist)
  // Maximum radius in which nodes are adjusted
  let maxradius = (nxgrid * nxgrid + nygrid * nygrid).to_double().sqrt()
  // Calculate the standard deviation for each row or column
  let stddata : Array[Double] = Array::make(nelements, 0.0)
  if transpose == 0 {
    for i = 0; i < nelements; i = i + 1 {
      let mut n = 0
      for j = 0; j < ndata; j = j + 1 {
        if mask[i][j] != 0 {
          stddata[i] = stddata[i] + data[i][j] * data[i][j]
          n = n + 1
        }
      }
      stddata[i] = if stddata[i] > 0.0 {
        (stddata[i] / n.to_double()).sqrt()
      } else {
        1.0
      }
    }
  } else {
    for i = 0; i < nelements; i = i + 1 {
      let mut n = 0
      for j = 0; j < ndata; j = j + 1 {
        if mask[j][i] != 0 {
          stddata[i] = stddata[i] + data[j][i] * data[j][i]
          n = n + 1
        }
      }
      stddata[i] = if stddata[i] > 0.0 {
        (stddata[i] / n.to_double()).sqrt()
      } else {
        1.0
      }
    }
  }
  // Dummy mask (all ones) used when comparing data with cell vectors
  let dummymask : Array[Array[Int]] = Array::new()
  if transpose == 0 {
    for _i = 0; _i < nygrid; _i = _i + 1 {
      dummymask.push(Array::make(ndata, 1))
    }
  } else {
    for _i = 0; _i < ndata; _i = _i + 1 {
      dummymask.push(Array::make(1, 1))
    }
  }
  // Randomly initialize the nodes
  for ix = 0; ix < nxgrid; ix = ix + 1 {
    for iy = 0; iy < nygrid; iy = iy + 1 {
      let mut sum = 0.0
      for i = 0; i < ndata; i = i + 1 {
        let term = -1.0 + 2.0 * bc_uniform()
        celldata[ix][iy][i] = term
        sum = sum + term * term
      }
      sum = (sum / ndata.to_double()).sqrt()
      for i = 0; i < ndata; i = i + 1 {
        celldata[ix][iy][i] = celldata[ix][iy][i] / sum
      }
    }
  }
  // Randomize the order in which genes or arrays will be used
  let index : Array[Int] = Array::make(nelements, 0)
  for i = 0; i < nelements; i = i + 1 {
    index[i] = i
  }
  for i = 0; i < nelements; i = i + 1 {
    let j = i + ((nelements - i).to_double() * bc_uniform()).to_int()
    let t = index[j]
    index[j] = index[i]
    index[i] = t
  }
  // Start the iteration
  for iter = 0; iter < niter; iter = iter + 1 {
    let mut ixbest = 0
    let mut iybest = 0
    let iobject = index[iter % nelements]
    let radius = maxradius * (1.0 - iter.to_double() / niter.to_double())
    let tau = inittau * (1.0 - iter.to_double() / niter.to_double())
    if transpose == 0 {
      let mut closest = metric(
        ndata,
        data,
        celldata[ixbest],
        mask,
        dummymask,
        weights,
        iobject,
        iybest,
        transpose,
      )
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          let distance = metric(
            ndata,
            data,
            celldata[ix],
            mask,
            dummymask,
            weights,
            iobject,
            iy,
            transpose,
          )
          if distance < closest {
            ixbest = ix
            iybest = iy
            closest = distance
          }
        }
      }
      // Update the cells within the neighbourhood radius
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          let dx = ix - ixbest
          let dy = iy - iybest
          if (dx * dx + dy * dy).to_double().sqrt() < radius {
            let mut sum = 0.0
            for i = 0; i < ndata; i = i + 1 {
              if mask[iobject][i] != 0 {
                celldata[ix][iy][i] = celldata[ix][iy][i] +
                  tau *
                  (data[iobject][i] / stddata[iobject] - celldata[ix][iy][i])
              }
            }
            for i = 0; i < ndata; i = i + 1 {
              let term = celldata[ix][iy][i]
              sum = sum + term * term
            }
            if sum > 0.0 {
              sum = (sum / ndata.to_double()).sqrt()
              for i = 0; i < ndata; i = i + 1 {
                celldata[ix][iy][i] = celldata[ix][iy][i] / sum
              }
            }
          }
        }
      }
    } else {
      // The cell vector is presented to the metric as a matrix with one
      // column; row i contains celldata[ix][iy][i].
      let cellvector : Array[Array[Double]] = Array::new()
      for _i = 0; _i < ndata; _i = _i + 1 {
        cellvector.push([0.0])
      }
      for i = 0; i < ndata; i = i + 1 {
        cellvector[i][0] = celldata[ixbest][iybest][i]
      }
      let mut closest = metric(
        ndata, data, cellvector, mask, dummymask, weights, iobject, 0, transpose,
      )
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          // Note: the C Clustering Library refills cellvector from
          // celldata[ixbest][iybest] inside the search loop (instead of
          // celldata[ix][iy]); this behaviour is mirrored exactly.
          for i = 0; i < ndata; i = i + 1 {
            cellvector[i][0] = celldata[ixbest][iybest][i]
          }
          let distance = metric(
            ndata, data, cellvector, mask, dummymask, weights, iobject, 0, transpose,
          )
          if distance < closest {
            ixbest = ix
            iybest = iy
            closest = distance
          }
        }
      }
      // Update the cells within the neighbourhood radius
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          let dx = ix - ixbest
          let dy = iy - iybest
          if (dx * dx + dy * dy).to_double().sqrt() < radius {
            let mut sum = 0.0
            for i = 0; i < ndata; i = i + 1 {
              if mask[i][iobject] != 0 {
                celldata[ix][iy][i] = celldata[ix][iy][i] +
                  tau *
                  (data[i][iobject] / stddata[iobject] - celldata[ix][iy][i])
              }
            }
            for i = 0; i < ndata; i = i + 1 {
              let term = celldata[ix][iy][i]
              sum = sum + term * term
            }
            if sum > 0.0 {
              sum = (sum / ndata.to_double()).sqrt()
              for i = 0; i < ndata; i = i + 1 {
                celldata[ix][iy][i] = celldata[ix][iy][i] / sum
              }
            }
          }
        }
      }
    }
  }
}

///|
/// Assign each item to the closest cell of the trained SOM
/// (port of somassign in cluster.c).
fn bc_somassign(
  nrows : Int,
  ncolumns : Int,
  data : Array[Array[Double]],
  mask : Array[Array[Int]],
  weights : Array[Double],
  transpose : Int,
  nxgrid : Int,
  nygrid : Int,
  celldata : Array[Array[Array[Double]]],
  dist : String,
) -> Array[(Int, Int)] {
  let ndata = if transpose == 0 { ncolumns } else { nrows }
  let metric = bc_setmetric(dist)
  let clusterid : Array[(Int, Int)] = Array::new()
  if transpose == 0 {
    let dummymask : Array[Array[Int]] = Array::new()
    for _i = 0; _i < nygrid; _i = _i + 1 {
      dummymask.push(Array::make(ncolumns, 1))
    }
    for i = 0; i < nrows; i = i + 1 {
      let mut ixbest = 0
      let mut iybest = 0
      let mut closest = metric(
        ndata,
        data,
        celldata[ixbest],
        mask,
        dummymask,
        weights,
        i,
        iybest,
        transpose,
      )
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          let distance = metric(
            ndata,
            data,
            celldata[ix],
            mask,
            dummymask,
            weights,
            i,
            iy,
            transpose,
          )
          if distance < closest {
            ixbest = ix
            iybest = iy
            closest = distance
          }
        }
      }
      clusterid.push((ixbest, iybest))
    }
  } else {
    let dummymask : Array[Array[Int]] = Array::new()
    for _i = 0; _i < nrows; _i = _i + 1 {
      dummymask.push(Array::make(1, 1))
    }
    for i = 0; i < ncolumns; i = i + 1 {
      let mut ixbest = 0
      let mut iybest = 0
      let cellvector : Array[Array[Double]] = Array::new()
      for _j = 0; _j < ndata; _j = _j + 1 {
        cellvector.push([0.0])
      }
      for j = 0; j < ndata; j = j + 1 {
        cellvector[j][0] = celldata[ixbest][iybest][j]
      }
      let mut closest = metric(
        ndata, data, cellvector, mask, dummymask, weights, i, 0, transpose,
      )
      for ix = 0; ix < nxgrid; ix = ix + 1 {
        for iy = 0; iy < nygrid; iy = iy + 1 {
          for j = 0; j < ndata; j = j + 1 {
            cellvector[j][0] = celldata[ix][iy][j]
          }
          let distance = metric(
            ndata, data, cellvector, mask, dummymask, weights, i, 0, transpose,
          )
          if distance < closest {
            ixbest = ix
            iybest = iy
            closest = distance
          }
        }
      }
      clusterid.push((ixbest, iybest))
    }
  }
  clusterid
}

///|
/// Implement a self-organizing map (Kohonen) on a rectangular grid
/// (port of somcluster). Rows (transpose=false, default) or columns
/// (transpose=true) of the data matrix are mapped to cells in an
/// nxgrid x nygrid grid.
///
/// The default values follow the Python API: nxgrid = 2, nygrid = 1,
/// inittau = 0.02, niter = 1, dist = "e".
///
/// Returns the cell coordinates for each item and the centroid vector of
/// each cell. The result depends on the internal random generator; use
/// cluster_seed() to make it reproducible.
pub fn somcluster(
  data : Array[Array[Double]],
  mask? : Array[Array[Int]],
  weight? : Array[Double],
  transpose? : Bool,
  nxgrid? : Int,
  nygrid? : Int,
  inittau? : Double,
  niter? : Int,
  dist? : String,
) -> SomResult {
  let tr = transpose.unwrap_or(false)
  let transpose_i = if tr { 1 } else { 0 }
  let nx = nxgrid.unwrap_or(2)
  let ny = nygrid.unwrap_or(1)
  let tau = inittau.unwrap_or(0.02)
  let iters = niter.unwrap_or(1)
  let d = dist.unwrap_or("e")
  let nrows = data.length()
  let ncolumns = if nrows > 0 { data[0].length() } else { 0 }
  let nobjects = if transpose_i == 0 { nrows } else { ncolumns }
  let ndata = if transpose_i == 0 { ncolumns } else { nrows }
  let m = mask.unwrap_or(bc_ones_mask(nrows, ncolumns))
  let w = weight.unwrap_or(bc_ones_weight(ndata))
  let celldata : Array[Array[Array[Double]]] = Array::new()
  for _ix = 0; _ix < nx; _ix = _ix + 1 {
    let plane : Array[Array[Double]] = Array::new()
    for _iy = 0; _iy < ny; _iy = _iy + 1 {
      plane.push(Array::make(ndata, 0.0))
    }
    celldata.push(plane)
  }
  let clusterid : Array[(Int, Int)] = Array::new()
  if nobjects >= 2 {
    bc_somworker(
      nrows, ncolumns, data, m, w, transpose_i, nx, ny, tau, celldata, iters, d,
    )
    let ids = bc_somassign(
      nrows, ncolumns, data, m, w, transpose_i, nx, ny, celldata, d,
    )
    for i = 0; i < ids.length(); i = i + 1 {
      clusterid.push(ids[i])
    }
  }
  { clusterid, celldata, }
}

// ---------------------------------------------------------------------------
// Principal component analysis (pca, via Golub-Reinsch SVD)
// ---------------------------------------------------------------------------

///|
/// Result of pca: the mean of each column, the coordinates of the data along
/// the principal components, the principal components themselves, and the
/// associated eigenvalues (sorted by magnitude, largest first).
pub struct ClusterPcaResult {
  columnmean : Array[Double]
  coordinates : Array[Array[Double]]
  pc : Array[Array[Double]]
  eigenvalues : Array[Double]
}

///|
pub fn ClusterPcaResult::columnmean(self : ClusterPcaResult) -> Array[Double] {
  self.columnmean
}

///|
pub fn ClusterPcaResult::coordinates(
  self : ClusterPcaResult,
) -> Array[Array[Double]] {
  self.coordinates
}

///|
pub fn ClusterPcaResult::pc(self : ClusterPcaResult) -> Array[Array[Double]] {
  self.pc
}

///|
pub fn ClusterPcaResult::eigenvalues(self : ClusterPcaResult) -> Array[Double] {
  self.eigenvalues
}

///|
/// Singular value decomposition A = U S V of a real m x n rectangular matrix,
/// following the Algol procedure svd of Golub and Reinsch (1970) as ported to
/// C by M. de Hoon (svd in cluster.c). Householder bidiagonalization and a
/// variant of the QR algorithm are used.
///
/// On input, u contains the m x n matrix A. On output, u contains U,
/// vt contains V, and w contains the n (non-negative) singular values.
///
/// Returns 0 on success, and k (> 0) if the k-th singular value was not
/// determined after 30 iterations.
fn bc_svd(
  m : Int,
  n : Int,
  u : Array[Array[Double]],
  w : Array[Double],
  vt : Array[Array[Double]],
) -> Int {
  let rv1 : Array[Double] = Array::make(n, 0.0)
  let mut ierr = 0
  let mut g = 0.0
  let mut scale = 0.0
  let mut anorm = 0.0
  let mut f = 0.0
  let mut h = 0.0
  let mut s = 0.0
  let mut x = 0.0
  let mut y = 0.0
  let mut z = 0.0
  let mut c = 0.0
  let mut l = 0
  if m >= n {
    // Householder reduction to bidiagonal form
    for i = 0; i < n; i = i + 1 {
      l = i + 1
      rv1[i] = scale * g
      g = 0.0
      s = 0.0
      scale = 0.0
      for k = i; k < m; k = k + 1 {
        scale = scale + u[k][i].abs()
      }
      if scale != 0.0 {
        for k = i; k < m; k = k + 1 {
          u[k][i] = u[k][i] / scale
          s = s + u[k][i] * u[k][i]
        }
        f = u[i][i]
        g = if f >= 0.0 { -s.sqrt() } else { s.sqrt() }
        h = f * g - s
        u[i][i] = f - g
        if i < n - 1 {
          for j = l; j < n; j = j + 1 {
            s = 0.0
            for k = i; k < m; k = k + 1 {
              s = s + u[k][i] * u[k][j]
            }
            f = s / h
            for k = i; k < m; k = k + 1 {
              u[k][j] = u[k][j] + f * u[k][i]
            }
          }
        }
        for k = i; k < m; k = k + 1 {
          u[k][i] = u[k][i] * scale
        }
      }
      w[i] = scale * g
      g = 0.0
      s = 0.0
      scale = 0.0
      if i < n - 1 {
        for k = l; k < n; k = k + 1 {
          scale = scale + u[i][k].abs()
        }
        if scale != 0.0 {
          for k = l; k < n; k = k + 1 {
            u[i][k] = u[i][k] / scale
            s = s + u[i][k] * u[i][k]
          }
          f = u[i][l]
          g = if f >= 0.0 { -s.sqrt() } else { s.sqrt() }
          h = f * g - s
          u[i][l] = f - g
          for k = l; k < n; k = k + 1 {
            rv1[k] = u[i][k] / h
          }
          for j = l; j < m; j = j + 1 {
            s = 0.0
            for k = l; k < n; k = k + 1 {
              s = s + u[j][k] * u[i][k]
            }
            for k = l; k < n; k = k + 1 {
              u[j][k] = u[j][k] + s * rv1[k]
            }
          }
          for k = l; k < n; k = k + 1 {
            u[i][k] = u[i][k] * scale
          }
        }
      }
      let t = w[i].abs() + rv1[i].abs()
      if t > anorm {
        anorm = t
      }
    }
    // Accumulation of right-hand transformations
    for i = n - 1; i >= 0; i = i - 1 {
      if i < n - 1 {
        if g != 0.0 {
          for j = l; j < n; j = j + 1 {
            vt[i][j] = u[i][j] / u[i][l] / g
          }
          // Double division avoids possible underflow
          for j = l; j < n; j = j + 1 {
            s = 0.0
            for k = l; k < n; k = k + 1 {
              s = s + u[i][k] * vt[j][k]
            }
            for k = l; k < n; k = k + 1 {
              vt[j][k] = vt[j][k] + s * vt[i][k]
            }
          }
        }
      }
      for j = l; j < n; j = j + 1 {
        vt[j][i] = 0.0
        vt[i][j] = 0.0
      }
      vt[i][i] = 1.0
      g = rv1[i]
      l = i
    }
    // Accumulation of left-hand transformations
    for i = n - 1; i >= 0; i = i - 1 {
      l = i + 1
      g = w[i]
      if i != n - 1 {
        for j = l; j < n; j = j + 1 {
          u[i][j] = 0.0
        }
      }
      if g != 0.0 {
        if i != n - 1 {
          for j = l; j < n; j = j + 1 {
            s = 0.0
            for k = l; k < m; k = k + 1 {
              s = s + u[k][i] * u[k][j]
            }
            // Double division avoids possible underflow
            f = s / u[i][i] / g
            for k = i; k < m; k = k + 1 {
              u[k][j] = u[k][j] + f * u[k][i]
            }
          }
        }
        for j = i; j < m; j = j + 1 {
          u[j][i] = u[j][i] / g
        }
      } else {
        for j = i; j < m; j = j + 1 {
          u[j][i] = 0.0
        }
      }
      u[i][i] = u[i][i] + 1.0
    }
    // Diagonalization of the bidiagonal form
    for k = n - 1; k >= 0; k = k - 1 {
      let k1 = k - 1
      let mut its = 0
      while true {
        // Test for splitting
        l = k
        while l >= 0 {
          let l1 = l - 1
          if rv1[l].abs() + anorm == anorm {
            break
          }
          // rv1[0] is always zero, so there is no exit
          // through the bottom of the loop
          if w[l1].abs() + anorm == anorm {
            // Cancellation of rv1[l] if l greater than 0
            c = 0.0
            s = 1.0
            for i = l; i <= k; i = i + 1 {
              f = s * rv1[i]
              rv1[i] = rv1[i] * c
              if f.abs() + anorm == anorm {
                break
              }
              g = w[i]
              h = (f * f + g * g).sqrt()
              w[i] = h
              c = g / h
              s = -f / h
              for j = 0; j < m; j = j + 1 {
                y = u[j][l1]
                z = u[j][i]
                u[j][l1] = y * c + z * s
                u[j][i] = -y * s + z * c
              }
            }
            break
          }
          l = l - 1
        }
        // Test for convergence
        z = w[k]
        if l == k {
          // Convergence
          if z < 0.0 {
            // w[k] is made non-negative
            w[k] = -z
            for j = 0; j < n; j = j + 1 {
              vt[k][j] = -vt[k][j]
            }
          }
          break
        } else if its == 30 {
          ierr = k
          break
        } else {
          // Shift from bottom 2 by 2 minor
          its = its + 1
          x = w[l]
          y = w[k1]
          g = rv1[k1]
          h = rv1[k]
          f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y)
          g = (f * f + 1.0).sqrt()
          f = (
              (x - z) * (x + z) +
              h * (y / (f + (if f >= 0.0 { g } else { -g })) - h)
            ) /
            x
          // Next QR transformation
          c = 1.0
          s = 1.0
          for i1 = l; i1 <= k1; i1 = i1 + 1 {
            let i = i1 + 1
            g = rv1[i]
            y = w[i]
            h = s * g
            g = c * g
            z = (f * f + h * h).sqrt()
            rv1[i1] = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            for j = 0; j < n; j = j + 1 {
              x = vt[i1][j]
              z = vt[i][j]
              vt[i1][j] = x * c + z * s
              vt[i][j] = -x * s + z * c
            }
            z = (f * f + h * h).sqrt()
            w[i1] = z
            // Rotation can be arbitrary if z is zero
            if z != 0.0 {
              c = f / z
              s = h / z
            }
            f = c * g + s * y
            x = -s * g + c * y
            for j = 0; j < m; j = j + 1 {
              y = u[j][i1]
              z = u[j][i]
              u[j][i1] = y * c + z * s
              u[j][i] = -y * s + z * c
            }
          }
          rv1[l] = 0.0
          rv1[k] = f
          w[k] = x
        }
      }
    }
  } else {
    // m < n
    // Householder reduction to bidiagonal form
    for i = 0; i < m; i = i + 1 {
      l = i + 1
      rv1[i] = scale * g
      g = 0.0
      s = 0.0
      scale = 0.0
      for k = i; k < n; k = k + 1 {
        scale = scale + u[i][k].abs()
      }
      if scale != 0.0 {
        for k = i; k < n; k = k + 1 {
          u[i][k] = u[i][k] / scale
          s = s + u[i][k] * u[i][k]
        }
        f = u[i][i]
        g = if f >= 0.0 { -s.sqrt() } else { s.sqrt() }
        h = f * g - s
        u[i][i] = f - g
        if i < m - 1 {
          for j = l; j < m; j = j + 1 {
            s = 0.0
            for k = i; k < n; k = k + 1 {
              s = s + u[i][k] * u[j][k]
            }
            f = s / h
            for k = i; k < n; k = k + 1 {
              u[j][k] = u[j][k] + f * u[i][k]
            }
          }
        }
        for k = i; k < n; k = k + 1 {
          u[i][k] = u[i][k] * scale
        }
      }
      w[i] = scale * g
      g = 0.0
      s = 0.0
      scale = 0.0
      if i < m - 1 {
        for k = l; k < m; k = k + 1 {
          scale = scale + u[k][i].abs()
        }
        if scale != 0.0 {
          for k = l; k < m; k = k + 1 {
            u[k][i] = u[k][i] / scale
            s = s + u[k][i] * u[k][i]
          }
          f = u[l][i]
          g = if f >= 0.0 { -s.sqrt() } else { s.sqrt() }
          h = f * g - s
          u[l][i] = f - g
          for k = l; k < m; k = k + 1 {
            rv1[k] = u[k][i] / h
          }
          for j = l; j < n; j = j + 1 {
            s = 0.0
            for k = l; k < m; k = k + 1 {
              s = s + u[k][j] * u[k][i]
            }
            for k = l; k < m; k = k + 1 {
              u[k][j] = u[k][j] + s * rv1[k]
            }
          }
          for k = l; k < m; k = k + 1 {
            u[k][i] = u[k][i] * scale
          }
        }
      }
      let t = w[i].abs() + rv1[i].abs()
      if t > anorm {
        anorm = t
      }
    }
    // Accumulation of right-hand transformations
    for i = m - 1; i >= 0; i = i - 1 {
      if i < m - 1 {
        if g != 0.0 {
          for j = l; j < m; j = j + 1 {
            vt[j][i] = u[j][i] / u[l][i] / g
          }
          // Double division avoids possible underflow
          for j = l; j < m; j = j + 1 {
            s = 0.0
            for k = l; k < m; k = k + 1 {
              s = s + u[k][i] * vt[k][j]
            }
            for k = l; k < m; k = k + 1 {
              vt[k][j] = vt[k][j] + s * vt[k][i]
            }
          }
        }
      }
      for j = l; j < m; j = j + 1 {
        vt[i][j] = 0.0
        vt[j][i] = 0.0
      }
      vt[i][i] = 1.0
      g = rv1[i]
      l = i
    }
    // Accumulation of left-hand transformations
    for i = m - 1; i >= 0; i = i - 1 {
      l = i + 1
      g = w[i]
      if i != m - 1 {
        for j = l; j < m; j = j + 1 {
          u[j][i] = 0.0
        }
      }
      if g != 0.0 {
        if i != m - 1 {
          for j = l; j < m; j = j + 1 {
            s = 0.0
            for k = l; k < n; k = k + 1 {
              s = s + u[i][k] * u[j][k]
            }
            // Double division avoids possible underflow
            f = s / u[i][i] / g
            for k = i; k < n; k = k + 1 {
              u[j][k] = u[j][k] + f * u[i][k]
            }
          }
        }
        for j = i; j < n; j = j + 1 {
          u[i][j] = u[i][j] / g
        }
      } else {
        for j = i; j < n; j = j + 1 {
          u[i][j] = 0.0
        }
      }
      u[i][i] = u[i][i] + 1.0
    }
    // Diagonalization of the bidiagonal form
    for k = m - 1; k >= 0; k = k - 1 {
      let k1 = k - 1
      let mut its = 0
      while true {
        // Test for splitting
        l = k
        while l >= 0 {
          let l1 = l - 1
          if rv1[l].abs() + anorm == anorm {
            break
          }
          // rv1[0] is always zero, so there is no exit
          // through the bottom of the loop
          if w[l1].abs() + anorm == anorm {
            // Cancellation of rv1[l] if l greater than 0
            c = 0.0
            s = 1.0
            for i = l; i <= k; i = i + 1 {
              f = s * rv1[i]
              rv1[i] = rv1[i] * c
              if f.abs() + anorm == anorm {
                break
              }
              g = w[i]
              h = (f * f + g * g).sqrt()
              w[i] = h
              c = g / h
              s = -f / h
              for j = 0; j < n; j = j + 1 {
                y = u[l1][j]
                z = u[i][j]
                u[l1][j] = y * c + z * s
                u[i][j] = -y * s + z * c
              }
            }
            break
          }
          l = l - 1
        }
        // Test for convergence
        z = w[k]
        if l == k {
          // Convergence
          if z < 0.0 {
            // w[k] is made non-negative
            w[k] = -z
            for j = 0; j < m; j = j + 1 {
              vt[j][k] = -vt[j][k]
            }
          }
          break
        } else if its == 30 {
          ierr = k
          break
        } else {
          // Shift from bottom 2 by 2 minor
          its = its + 1
          x = w[l]
          y = w[k1]
          g = rv1[k1]
          h = rv1[k]
          f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y)
          g = (f * f + 1.0).sqrt()
          f = (
              (x - z) * (x + z) +
              h * (y / (f + (if f >= 0.0 { g } else { -g })) - h)
            ) /
            x
          // Next QR transformation
          c = 1.0
          s = 1.0
          for i1 = l; i1 <= k1; i1 = i1 + 1 {
            let i = i1 + 1
            g = rv1[i]
            y = w[i]
            h = s * g
            g = c * g
            z = (f * f + h * h).sqrt()
            rv1[i1] = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            for j = 0; j < m; j = j + 1 {
              x = vt[j][i1]
              z = vt[j][i]
              vt[j][i1] = x * c + z * s
              vt[j][i] = -x * s + z * c
            }
            z = (f * f + h * h).sqrt()
            w[i1] = z
            // Rotation can be arbitrary if z is zero
            if z != 0.0 {
              c = f / z
              s = h / z
            }
            f = c * g + s * y
            x = -s * g + c * y
            for j = 0; j < n; j = j + 1 {
              y = u[i1][j]
              z = u[i][j]
              u[i1][j] = y * c + z * s
              u[i][j] = -y * s + z * c
            }
          }
          rv1[l] = 0.0
          rv1[k] = f
          w[k] = x
        }
      }
    }
  }
  ierr
}

///|
/// Principal component analysis of the mean-centered data matrix
/// (port of pca in cluster.c). On return, u contains the coordinates
/// (nrows >= ncolumns) or the principal component vectors (nrows <
/// ncolumns); see bc_pca for the exact layout. Returns 0 on success.
fn bc_pca_core(
  nrows : Int,
  ncolumns : Int,
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  w : Array[Double],
) -> Int {
  let error = bc_svd(nrows, ncolumns, u, w, v)
  if error == 0 {
    if nrows >= ncolumns {
      for j = 0; j < ncolumns; j = j + 1 {
        let s = w[j]
        for i = 0; i < nrows; i = i + 1 {
          u[i][j] = u[i][j] * s
        }
      }
      let index : Array[Int] = Array::make(ncolumns, 0)
      bc_sort_index(ncolumns, w, index)
      for i = 0; i < ncolumns / 2; i = i + 1 {
        let j = index[i]
        index[i] = index[ncolumns - 1 - i]
        index[ncolumns - 1 - i] = j
      }
      let temp : Array[Double] = Array::make(ncolumns, 0.0)
      for i = 0; i < nrows; i = i + 1 {
        for j = 0; j < ncolumns; j = j + 1 {
          temp[j] = u[i][index[j]]
        }
        for j = 0; j < ncolumns; j = j + 1 {
          u[i][j] = temp[j]
        }
      }
      for i = 0; i < ncolumns; i = i + 1 {
        for j = 0; j < ncolumns; j = j + 1 {
          temp[j] = v[index[j]][i]
        }
        for j = 0; j < ncolumns; j = j + 1 {
          v[j][i] = temp[j]
        }
      }
      for i = 0; i < ncolumns; i = i + 1 {
        temp[i] = w[index[i]]
      }
      for i = 0; i < ncolumns; i = i + 1 {
        w[i] = temp[i]
      }
    } else {
      // nrows < ncolumns
      for j = 0; j < nrows; j = j + 1 {
        let s = w[j]
        for i = 0; i < nrows; i = i + 1 {
          v[i][j] = v[i][j] * s
        }
      }
      let index : Array[Int] = Array::make(nrows, 0)
      bc_sort_index(nrows, w, index)
      for i = 0; i < nrows / 2; i = i + 1 {
        let j = index[i]
        index[i] = index[nrows - 1 - i]
        index[nrows - 1 - i] = j
      }
      let temp : Array[Double] = Array::make(nrows, 0.0)
      for j = 0; j < ncolumns; j = j + 1 {
        for i = 0; i < nrows; i = i + 1 {
          temp[i] = u[index[i]][j]
        }
        for i = 0; i < nrows; i = i + 1 {
          u[i][j] = temp[i]
        }
      }
      for j = 0; j < nrows; j = j + 1 {
        for i = 0; i < nrows; i = i + 1 {
          temp[i] = v[j][index[i]]
        }
        for i = 0; i < nrows; i = i + 1 {
          v[j][i] = temp[i]
        }
      }
      for i = 0; i < nrows; i = i + 1 {
        temp[i] = w[index[i]]
      }
      for i = 0; i < nrows; i = i + 1 {
        w[i] = temp[i]
      }
    }
  }
  error
}

///|
/// Error raised when the singular value decomposition in pca fails to
/// converge, or when an Eisen Cluster/TreeView file is malformed.
pub suberror ClusterError {
  ClusterError(String)
}

///|
/// Perform principal component analysis (port of the Python Bio.Cluster.pca).
///
/// The function returns the mean of each column, the principal components as
/// an nmin x ncolumns array (nmin = min(nrows, ncolumns)), the coordinates of
/// the data along the principal components as an nrows x nmin array, and the
/// associated eigenvalues. The principal components, the coordinates, and the
/// eigenvalues are sorted by the magnitude of the eigenvalue, with the largest
/// eigenvalues appearing first. Adding the column means to the dot product of
/// the coordinates and the principal components recreates the data matrix:
///   data ~= columnmean + dot(coordinates, pc)
pub fn cluster_pca(
  data : Array[Array[Double]],
) -> ClusterPcaResult raise ClusterError {
  let nrows = data.length()
  let ncols = if nrows > 0 { data[0].length() } else { 0 }
  if nrows == 0 || ncols == 0 {
    raise ClusterError("data matrix is empty")
  }
  let nmin = if nrows < ncols { nrows } else { ncols }
  // Calculate the mean of each column
  let columnmean : Array[Double] = Array::make(ncols, 0.0)
  for j = 0; j < ncols; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < nrows; i = i + 1 {
      s = s + data[i][j]
    }
    columnmean[j] = s / nrows.to_double()
  }
  // Subtract the mean of each column. Depending on the shape of the matrix,
  // the mean-centered data is loaded into the coordinates buffer (nrows >=
  // ncolumns) or into the pc buffer (nrows < ncolumns).
  let coordinates : Array[Array[Double]] = Array::new()
  let pc : Array[Array[Double]] = Array::new()
  if nrows >= ncols {
    for i = 0; i < nrows; i = i + 1 {
      let row : Array[Double] = Array::new()
      for j = 0; j < ncols; j = j + 1 {
        row.push(data[i][j] - columnmean[j])
      }
      coordinates.push(row)
    }
    for _i = 0; _i < ncols; _i = _i + 1 {
      pc.push(Array::make(ncols, 0.0))
    }
  } else {
    for i = 0; i < nrows; i = i + 1 {
      let row : Array[Double] = Array::new()
      for j = 0; j < ncols; j = j + 1 {
        row.push(data[i][j] - columnmean[j])
      }
      pc.push(row)
    }
    for _i = 0; _i < nrows; _i = _i + 1 {
      coordinates.push(Array::make(nrows, 0.0))
    }
  }
  let eigenvalues : Array[Double] = Array::make(nmin, 0.0)
  let error = bc_pca_core(
    nrows,
    ncols,
    if nrows >= ncols {
      coordinates
    } else {
      pc
    },
    if nrows >= ncols {
      pc
    } else {
      coordinates
    },
    eigenvalues,
  )
  if error > 0 {
    raise ClusterError("Singular value decomposition failed to converge")
  }
  { columnmean, coordinates, pc, eigenvalues, }
}

// ---------------------------------------------------------------------------
// Record: Cluster/TreeView (Michael Eisen) file format
// ---------------------------------------------------------------------------

///|
/// Store gene expression data and related information contained in a data
/// file following the file format defined for Michael Eisen's
/// Cluster/TreeView program (port of Bio.Cluster.Record).
///
/// Fields (via getters):
///   - data      : the gene expression data matrix
///   - mask      : 1's and 0's denoting which values are present (1) or
///                 missing (0); None if no data are missing
///   - geneid    : a unique identifier for each gene (e.g., ORF name)
///   - genename  : an additional description for each gene (e.g., gene name)
///   - gweight   : the weight to be used for each gene when calculating
///                 the distance
///   - gorder    : the preferred order of the genes in the output file
///   - expid     : a unique identifier for each sample
///   - eweight   : the weight to be used for each sample when calculating
///                 the distance
///   - eorder    : the preferred order of the samples in the output file
///   - uniqid    : the string that was used instead of UNIQID in the
///                 input file
pub struct Record {
  data : Array[Array[Double]]
  mask : Array[Array[Int]]?
  geneid : Array[String]
  genename : Array[String]?
  gweight : Array[Double]?
  gorder : Array[Double]?
  expid : Array[String]
  eweight : Array[Double]?
  eorder : Array[Double]?
  uniqid : String
}

///|
pub fn Record::data(self : Record) -> Array[Array[Double]] {
  self.data
}

///|
pub fn Record::mask(self : Record) -> Array[Array[Int]]? {
  self.mask
}

///|
pub fn Record::geneid(self : Record) -> Array[String] {
  self.geneid
}

///|
pub fn Record::genename(self : Record) -> Array[String]? {
  self.genename
}

///|
pub fn Record::gweight(self : Record) -> Array[Double]? {
  self.gweight
}

///|
pub fn Record::gorder(self : Record) -> Array[Double]? {
  self.gorder
}

///|
pub fn Record::expid(self : Record) -> Array[String] {
  self.expid
}

///|
pub fn Record::eweight(self : Record) -> Array[Double]? {
  self.eweight
}

///|
pub fn Record::eorder(self : Record) -> Array[Double]? {
  self.eorder
}

///|
pub fn Record::uniqid(self : Record) -> String {
  self.uniqid
}

///|
/// Remove '\r' and '\n' characters from both ends of a line
/// (mirrors Python's line.strip("\r\n")).
fn bc_strip_eol(line : String) -> String {
  let n = line.length()
  let mut end = n
  while end > 0 {
    let c = line.unsafe_get(end - 1)
    if c == '\r' || c == '\n' {
      end = end - 1
    } else {
      break
    }
  }
  let mut start = 0
  while start < end {
    let c = line.unsafe_get(start)
    if c == '\r' || c == '\n' {
      start = start + 1
    } else {
      break
    }
  }
  if start == 0 && end == n {
    line
  } else {
    line.substring(start~, end~).to_string()
  }
}

///|
/// Split a string on tab characters into an array of Strings.
fn bc_split_tabs(line : String) -> Array[String] {
  let parts : Array[String] = Array::new()
  for part in line.split("\t") {
    parts.push(part.to_string())
  }
  parts
}

///|
/// Read gene expression data from a file in the format defined for Michael
/// Eisen's Cluster/TreeView program, and return a Record (port of
/// Bio.Cluster.read / Record.__init__). The content should be the complete
/// text file; lines are separated by newlines and columns by tabs.
pub fn bio_cluster_read(content : String) -> Record raise ClusterError {
  let lines : Array[String] = Array::new()
  for part in content.split("\n") {
    lines.push(bc_strip_eol(part.to_string()))
  }
  // Drop trailing empty lines
  while lines.length() > 0 && lines[lines.length() - 1].length() == 0 {
    let _ = lines.pop()
  }
  if lines.length() == 0 {
    raise ClusterError("Empty input file")
  }
  let header = bc_split_tabs(lines[0])
  let n = header.length()
  let uniqid = header[0]
  let expid : Array[String] = Array::new()
  // Column layout: 0 is always the GENEID column; NAME, GWEIGHT and GORDER
  // columns are optional and identified by their header word.
  let mut name_col = -1
  let mut gweight_col = -1
  let mut gorder_col = -1
  for i = 1; i < n; i = i + 1 {
    let word = header[i]
    if word == "NAME" {
      name_col = i
    } else if word == "GWEIGHT" {
      gweight_col = i
    } else if word == "GORDER" {
      gorder_col = i
    } else {
      expid.push(word)
    }
  }
  let geneid : Array[String] = Array::new()
  let genename : Array[String] = Array::new()
  let gweight : Array[Double] = Array::new()
  let gorder : Array[Double] = Array::new()
  let data : Array[Array[Double]] = Array::new()
  let mask : Array[Array[Int]] = Array::new()
  let mut eweight : Array[Double]? = None
  let mut eorder : Array[Double]? = None
  let mut needmask = false
  for li = 1; li < lines.length(); li = li + 1 {
    let line = bc_split_tabs(lines[li])
    if line.length() != n {
      raise ClusterError(
        "Line with \{line.length()} columns found (expected \{n})",
      )
    }
    if line[0] == "EWEIGHT" {
      let row : Array[Double] = Array::new()
      for i = n - expid.length(); i < n; i = i + 1 {
        row.push(bc_parse_double(line[i]))
      }
      eweight = Some(row)
      continue
    }
    if line[0] == "EORDER" {
      let row : Array[Double] = Array::new()
      for i = n - expid.length(); i < n; i = i + 1 {
        row.push(bc_parse_double(line[i]))
      }
      eorder = Some(row)
      continue
    }
    let rowdata : Array[Double] = Array::new()
    let rowmask : Array[Int] = Array::new()
    for i = 0; i < n; i = i + 1 {
      let word = line[i]
      if i == 0 {
        geneid.push(word)
      } else if i == name_col {
        genename.push(word)
      } else if i == gweight_col {
        gweight.push(bc_parse_double(word))
      } else if i == gorder_col {
        gorder.push(bc_parse_double(word))
      } else if word.length() == 0 {
        rowdata.push(0.0)
        rowmask.push(0)
        needmask = true
      } else {
        rowdata.push(bc_parse_double(word))
        rowmask.push(1)
      }
    }
    data.push(rowdata)
    mask.push(rowmask)
  }
  let mask_opt : Array[Array[Int]]? = if needmask { Some(mask) } else { None }
  let gweight_opt : Array[Double]? = if gweight.length() > 0 {
    Some(gweight)
  } else {
    None
  }
  let gorder_opt : Array[Double]? = if gorder.length() > 0 {
    Some(gorder)
  } else {
    None
  }
  let genename_opt : Array[String]? = if genename.length() > 0 {
    Some(genename)
  } else {
    None
  }
  {
    data,
    mask: mask_opt,
    geneid,
    genename: genename_opt,
    gweight: gweight_opt,
    gorder: gorder_opt,
    expid,
    eweight,
    eorder,
    uniqid,
  }
}

///|
/// Apply hierarchical clustering to the expression data in this Record and
/// return a ClusterTree. The pairwise single ("s"), complete ("m", default),
/// centroid ("c"), and average ("a") linkage methods are available.
pub fn Record::treecluster(
  self : Record,
  transpose? : Bool,
  method? : String,
  dist? : String,
) -> ClusterTree {
  let tr = transpose.unwrap_or(false)
  let weight = if tr { self.gweight } else { self.eweight }
  treecluster(self.data, mask?=self.mask, weight?, transpose=tr, method?, dist?)
}

///|
/// Apply k-means or k-median clustering to the data in this Record
/// (see kcluster).
pub fn Record::kcluster(
  self : Record,
  nclusters? : Int,
  transpose? : Bool,
  npass? : Int,
  method? : String,
  dist? : String,
  initialid? : Array[Int],
) -> KClusterResult {
  let tr = transpose.unwrap_or(false)
  let weight = if tr { self.gweight } else { self.eweight }
  kcluster(
    self.data,
    nclusters?,
    mask?=self.mask,
    weight?,
    transpose=tr,
    npass?,
    method?,
    dist?,
    initialid?,
  )
}

///|
/// Calculate a self-organizing map on the data in this Record
/// (see somcluster).
pub fn Record::somcluster(
  self : Record,
  transpose? : Bool,
  nxgrid? : Int,
  nygrid? : Int,
  inittau? : Double,
  niter? : Int,
  dist? : String,
) -> SomResult {
  let tr = transpose.unwrap_or(false)
  let weight = if tr { self.gweight } else { self.eweight }
  somcluster(
    self.data,
    mask?=self.mask,
    weight?,
    transpose=tr,
    nxgrid?,
    nygrid?,
    inittau?,
    niter?,
    dist?,
  )
}

///|
/// Calculate the centroid of each cluster in this Record
/// (see clustercentroids).
pub fn Record::clustercentroids(
  self : Record,
  clusterid? : Array[Int],
  method? : String,
  transpose? : Bool,
) -> CentroidsResult {
  clustercentroids(self.data, mask?=self.mask, clusterid?, method?, transpose?)
}

///|
/// Calculate the distance between two clusters in this Record
/// (see clusterdistance).
pub fn Record::clusterdistance(
  self : Record,
  index1? : Array[Int],
  index2? : Array[Int],
  method? : String,
  dist? : String,
  transpose? : Bool,
) -> Double {
  let tr = transpose.unwrap_or(false)
  let weight = if tr { self.gweight } else { self.eweight }
  clusterdistance(
    self.data,
    mask?=self.mask,
    weight?,
    index1?,
    index2?,
    method?,
    dist?,
    transpose=tr,
  )
}

///|
/// Calculate the distance matrix of the data in this Record
/// (see distancematrix).
pub fn Record::distancematrix(
  self : Record,
  transpose? : Bool,
  dist? : String,
) -> Array[Array[Double]] {
  let tr = transpose.unwrap_or(false)
  let weight = if tr { self.gweight } else { self.eweight }
  distancematrix(self.data, mask?=self.mask, weight?, transpose=tr, dist?)
}
