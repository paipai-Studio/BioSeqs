///|
/// Bio.LogisticRegression - Binary logistic regression classifier using
/// Newton-Raphson optimization, inspired by Biopython Bio.LogisticRegression.
/// Reference: Chang J (2002) Biopython LogisticRegression module.

// ===========================================================================
// Model struct
// ===========================================================================

///|
/// Logistic regression model holding the weight vector beta.
/// beta[0] is the bias (intercept) term; beta[1..] are feature weights.
pub struct LogisticRegressionModel {
  beta : Array[Double]
}

///|
/// Construct a new LogisticRegressionModel with the given beta weights.
pub fn LogisticRegressionModel::new(
  beta : Array[Double]
) -> LogisticRegressionModel {
  LogisticRegressionModel::{ beta }
}

///|
/// Accessor for the beta weight vector.
pub fn LogisticRegressionModel::beta(self : LogisticRegressionModel) -> Array[
  Double
] {
  self.beta
}

///|
/// Number of weights (including bias).
pub fn LogisticRegressionModel::n_weights(
  self : LogisticRegressionModel
) -> Int {
  self.beta.length()
}

// ===========================================================================
// Internal helpers
// ===========================================================================

///|
/// Sigmoid function: 1 / (1 + exp(-x)), with overflow protection.
fn lr_sigmoid(x : Double) -> Double {
  if x >= 40.0 {
    return 1.0
  }
  if x <= -40.0 {
    return 0.0000000001 // ~0 but non-zero to avoid log(0)
  }
  1.0 / (1.0 + @math.exp(-x))
}

///|
/// Safe natural log, clamping the argument to avoid -inf.
fn lr_safe_log(x : Double) -> Double {
  if x <= 0.0 {
    return -40.0
  }
  @math.ln(x)
}

///|
/// Dot product of two double arrays.
fn lr_dot(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut s = 0.0
  for i in 0..<n {
    s = s + a[i] * b[i]
  }
  s
}

///|
/// Solve a linear system A*x = b using Gaussian elimination with partial
/// pivoting. A is n×n (row-major), b is length n. Returns the solution x,
/// or None if A is singular.
fn lr_solve(matrix : Array[Array[Double]], b : Array[Double]) -> Array[Double]? {
  let n = matrix.length()
  if n == 0 {
    return Some([])
  }
  // Augmented matrix [matrix | b]
  let m : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(n + 1, 0.0)
    for j in 0..<n {
      row[j] = matrix[i][j]
    }
    row[n] = b[i]
    m.push(row)
  }
  // Forward elimination with partial pivoting.
  for col in 0..<n {
    // Find pivot row.
    let mut pivot = col
    let mut max_val = (m[col][col])
    if max_val < 0.0 {
      max_val = -max_val
    }
    let mut r = col + 1
    while r < n {
      let val = (m[r][col])
      let abs_val = if val < 0.0 { -val } else { val }
      if abs_val > max_val {
        max_val = abs_val
        pivot = r
      }
      r = r + 1
    }
    if max_val < 0.0000000001 {
      return None // singular
    }
    if pivot != col {
      let tmp = m[col]
      m[col] = m[pivot]
      m[pivot] = tmp
    }
    // Eliminate below.
    let mut r2 = col + 1
    while r2 < n {
      let factor = m[r2][col] / m[col][col]
      for c in col..<(n + 1) {
        m[r2][c] = m[r2][c] - factor * m[col][c]
      }
      r2 = r2 + 1
    }
  }
  // Back substitution.
  let x = Array::make(n, 0.0)
  let mut i = n - 1
  while i >= 0 {
    let mut s = m[i][n]
    let mut j = i + 1
    while j < n {
      s = s - m[i][j] * x[j]
      j = j + 1
    }
    x[i] = s / m[i][i]
    i = i - 1
  }
  Some(x)
}

// ===========================================================================
// Public API
// ===========================================================================

///|
/// Train a logistic regression classifier on a training set.
///
/// `xs` is an array of observations (each observation is an array of Double
/// feature values). `ys` is an array of class assignments (0 or 1).
/// Uses Newton-Raphson optimization with step-size halving.
///
/// Returns a trained `LogisticRegressionModel`.
pub fn logistic_regression_train(
  xs : Array[Array[Double]],
  ys : Array[Int]
) -> LogisticRegressionModel {
  let n = xs.length()
  // Each observation gets a leading constant 1.0 for the bias term.
  let ndims = if n > 0 { xs[0].length() + 1 } else { 1 }
  // Build design matrix x_design (n × ndims) with leading 1.0.
  let x_design : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(ndims, 1.0)
    for j in 1..<ndims {
      row[j] = xs[i][j - 1]
    }
    x_design.push(row)
  }
  // Initialize beta to zeros.
  let beta = Array::make(ndims, 0.0)
  let max_iterations = 500
  let converge_threshold = 0.01
  let mut stepsize = 1.0
  let mut old_llik : Double? = None
  let old_beta : Array[Double] = Array::make(ndims, 0.0)
  for _ in 0..<max_iterations {
    // Calculate probabilities p_i = sigmoid(beta · x_design_i).
    let p = Array::make(n, 0.5)
    for i in 0..<n {
      p[i] = lr_sigmoid(lr_dot(beta, x_design[i]))
    }
    // Log-likelihood.
    let mut llik = 0.0
    for i in 0..<n {
      let yi = ys[i].to_double()
      llik = llik + yi * lr_safe_log(p[i]) + (1.0 - yi) * lr_safe_log(1.0 - p[i])
    }
    // Check convergence / adjust step size.
    match old_llik {
      Some(old) => {
        if llik < old {
          stepsize = stepsize / 2.0
          // Restore old beta.
          for k in 0..<ndims {
            beta[k] = old_beta[k]
          }
        }
        let diff = llik - old
        let abs_diff = if diff < 0.0 { -diff } else { diff }
        if abs_diff <= converge_threshold {
          break
        }
      }
      None => ()
    }
    old_llik = Some(llik)
    for k in 0..<ndims {
      old_beta[k] = beta[k]
    }
    // Gradient: x_design^T (y - p).
    let grad = Array::make(ndims, 0.0)
    for j in 0..<ndims {
      let mut s = 0.0
      for i in 0..<n {
        s = s + x_design[i][j] * (ys[i].to_double() - p[i])
      }
      grad[j] = s
    }
    // Hessian: x_design^T W x_design, where W = diag(p*(1-p)).
    let h : Array[Array[Double]] = Array::new()
    for _ in 0..<ndims {
      h.push(Array::make(ndims, 0.0))
    }
    for i in 0..<n {
      let w = p[i] * (1.0 - p[i])
      for j in 0..<ndims {
        for k in 0..<ndims {
          h[j][k] = h[j][k] + x_design[i][j] * w * x_design[i][k]
        }
      }
    }
    // Solve h * delta = grad.
    match lr_solve(h, grad) {
      Some(delta) => {
        for k in 0..<ndims {
          let d = if stepsize != 1.0 { delta[k] * stepsize } else { delta[k] }
          beta[k] = beta[k] + d
        }
      }
      None => {
        break
      }
    }
  }
  LogisticRegressionModel::{ beta }
}

///|
/// Calculate the probability for each class given an observation.
/// Returns [P(class=0), P(class=1)].
pub fn logistic_regression_calculate(
  model : LogisticRegressionModel,
  x : Array[Double]
) -> Array[Double] {
  // Prepend constant 1.0 for bias.
  let xfull = Array::make(x.length() + 1, 1.0)
  for i in 0..<x.length() {
    xfull[i + 1] = x[i]
  }
  let score = lr_dot(model.beta, xfull)
  let p = lr_sigmoid(score)
  [1.0 - p, p]
}

///|
/// Classify an observation into class 0 or 1.
pub fn logistic_regression_classify(
  model : LogisticRegressionModel,
  x : Array[Double]
) -> Int {
  let probs = logistic_regression_calculate(model, x)
  if probs[0] > probs[1] {
    0
  } else {
    1
  }
}

// ===========================================================================
// Sample data: operon prediction (Bacillus subtilis gene pairs)
// ===========================================================================

///|
/// Sample training data for operon prediction (Biopython tutorial example).
/// Each observation: [intergenic_distance, expression_score].
/// Label 1 = same operon (OP), 0 = different operon (NOP).
pub fn logistic_regression_sample_data() -> (Array[
  Array[Double]
], Array[Int]) {
  let xs = [
    [-53.0, -200.78],
    [117.0, -267.14],
    [57.0, -163.47],
    [16.0, -190.30],
    [11.0, -220.94],
    [85.0, -193.94],
    [16.0, -182.71],
    [15.0, -180.41],
    [-26.0, -181.73],
    [58.0, -259.87],
    [126.0, -414.53],
    [191.0, -249.57],
    [113.0, -265.28],
    [145.0, -312.99],
    [154.0, -213.83],
    [147.0, -380.85],
    [93.0, -291.13],
  ]
  let ys = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
  (xs, ys)
}
