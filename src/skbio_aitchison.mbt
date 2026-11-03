/// Aitchison distance and compositional distance metrics.
///
/// scikit-bio equivalent: `composition.aitchison_distance()`
///
/// The Aitchison distance is the Euclidean distance in CLR space:
/// d(x, y) = ||clr(x) - clr(y)||

///|
/// Compute CLR (centered log-ratio) transform.
fn aitchison_clr(x : Array[Double]) -> Array[Double] {
  let n = x.length()
  if n == 0 {
    return []
  }
  // Geometric mean
  let mut log_sum = 0.0
  for xi in x {
    log_sum = log_sum + @math.ln(xi)
  }
  let gm = log_sum / n.to_double()
  let result : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    result[i] = @math.ln(x[i]) - gm
  }
  result
}

///|
/// Compute Aitchison distance between two compositions.
/// d(x, y) = sqrt(sum((clr(x_i) - clr(y_i))^2))
pub fn skbio_aitchison_distance(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() || x.length() == 0 {
    return 0.0
  }
  let clr_x = aitchison_clr(x)
  let clr_y = aitchison_clr(y)
  let mut sum = 0.0
  for i = 0; i < clr_x.length(); i = i + 1 {
    let d = clr_x[i] - clr_y[i]
    sum = sum + d * d
  }
  sum.sqrt()
}

///|
/// Compute Aitchison inner product.
/// <x, y>_a = sum(clr(x_i) * clr(y_i))
pub fn skbio_aitchison_inner(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() || x.length() == 0 {
    return 0.0
  }
  let clr_x = aitchison_clr(x)
  let clr_y = aitchison_clr(y)
  let mut sum = 0.0
  for i = 0; i < clr_x.length(); i = i + 1 {
    sum = sum + clr_x[i] * clr_y[i]
  }
  sum
}

///|
/// Compute Aitchison norm.
/// ||x||_a = sqrt(<x, x>_a)
pub fn skbio_aitchison_norm(x : Array[Double]) -> Double {
  skbio_aitchison_inner(x, x).sqrt()
}

///|
/// Compute Aitchison angle (cosine of angle in Aitchison geometry).
pub fn skbio_aitchison_cosine(x : Array[Double], y : Array[Double]) -> Double {
  let inner = skbio_aitchison_inner(x, y)
  let norm_x = skbio_aitchison_norm(x)
  let norm_y = skbio_aitchison_norm(y)
  if norm_x == 0.0 || norm_y == 0.0 {
    return 0.0
  }
  inner / (norm_x * norm_y)
}

///|
/// Compute pairwise Aitchison distance matrix.
pub fn skbio_aitchison_distance_matrix(
  compositions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = compositions.length()
  let result : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    result[i] = Array::make(n, 0.0)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let d = skbio_aitchison_distance(compositions[i], compositions[j])
      result[i][j] = d
      result[j][i] = d
    }
  }
  result
}

///|
/// Compute compositional distance using a specified method.
/// Methods: "aitchison", "bray_curtis", "jaccard", "euclidean"
pub fn skbio_composition_distance(
  x : Array[Double],
  y : Array[Double],
  method~ : String = "aitchison",
) -> Double {
  match method {
    "aitchison" => skbio_aitchison_distance(x, y)
    "euclidean" => {
      let mut sum = 0.0
      for i = 0; i < x.length(); i = i + 1 {
        if i < y.length() {
          let d = x[i] - y[i]
          sum = sum + d * d
        }
      }
      sum.sqrt()
    }
    "bray_curtis" => {
      let mut numerator = 0.0
      let mut denominator = 0.0
      for i = 0; i < x.length(); i = i + 1 {
        if i < y.length() {
          numerator = numerator + (x[i] - y[i]).abs()
          denominator = denominator + x[i] + y[i]
        }
      }
      if denominator > 0.0 {
        numerator / denominator
      } else {
        0.0
      }
    }
    _ => skbio_aitchison_distance(x, y)
  }
}
