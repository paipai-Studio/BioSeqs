// scvi_autograd.mbt - Minimal reverse-mode automatic differentiation for scVI models
//
// The engine is deliberately small: it supports vectors (Array[Double]) and
// scalar nodes (length-one vectors), dense affine layers, element-wise
// activations and the likelihood terms needed by scVI-style variational
// models. Parameters are stored in the tape as matrices / vectors so that
// gradients can be consumed by the optimisers in scvi_optim.mbt.

// ===========================================================================
// Operations
// ===========================================================================

///|
/// Computation-graph operation recorded for each node.
pub enum ScviAgOp {
  AgLeaf
  AgAdd(Int32, Int32)
  AgSub(Int32, Int32)
  AgMul(Int32, Int32)
  AgScale(Int32, Double)
  AgDense(Int32, Int32, Int32) // input node, weight matrix, bias vector
  AgRelu(Int32)
  AgSigmoid(Int32)
  AgSoftmax(Int32)
  AgLog(Int32)
  AgExp(Int32)
  AgSquare(Int32)
  AgNeg(Int32)
  AgSum(Int32)
  AgMean(Int32)
  AgSelect(Int32, Int32)
  AgKl(Int32, Int32)
  AgNbNll(Int32, Array[Double], Array[Double]) // mu node, counts, theta
  AgPoissonNll(Int32, Array[Double]) // mu node, counts
}

// ===========================================================================
// Parameters and tape
// ===========================================================================

///|
/// Trainable matrix parameter with its accumulated gradient.
pub struct ScviAgMatrix {
  value : Array[Array[Double]]
  grad : Array[Array[Double]]
}

///|
/// Trainable vector parameter with its accumulated gradient.
pub struct ScviAgVector {
  value : Array[Double]
  grad : Array[Double]
}

///|
/// A reverse-mode tape holding forward values, gradients and parameter slots.
pub struct ScviAgTape {
  vals : Array[Array[Double]]
  grads : Array[Array[Double]]
  ops : Array[ScviAgOp]
  mats : Array[ScviAgMatrix]
  vecs : Array[ScviAgVector]
}

///|
fn scvi_ag_zero_matrix(rows : Int32, cols : Int32) -> Array[Array[Double]] {
  let out = Array(capacity=rows)
  for _i in 0..<rows {
    let row = Array::make(cols, 0.0)
    out.push(row)
  }
  out
}

///|
fn scvi_ag_push(t : ScviAgTape, value : Array[Double], op : ScviAgOp) -> Int32 {
  let id = t.vals.length()
  t.vals.push(value)
  t.grads.push(Array::make(value.length(), 0.0))
  t.ops.push(op)
  id
}

///|
fn scvi_ag_accum(t : ScviAgTape, id : Int32, g : Array[Double]) -> Unit {
  let grad = t.grads[id]
  let n = if grad.length() < g.length() { grad.length() } else { g.length() }
  for i in 0..<n {
    grad[i] += g[i]
  }
}

///|
fn scvi_ag_zero_grads(t : ScviAgTape) -> Unit {
  for i in 0..<t.grads.length() {
    let g = t.grads[i]
    for j in 0..<g.length() {
      g[j] = 0.0
    }
  }
  for i in 0..<t.mats.length() {
    let g = t.mats[i].grad
    for r in 0..<g.length() {
      for c in 0..<g[r].length() {
        g[r][c] = 0.0
      }
    }
  }
  for i in 0..<t.vecs.length() {
    let g = t.vecs[i].grad
    for j in 0..<g.length() {
      g[j] = 0.0
    }
  }
}

///|
fn scvi_ag_elementwise_binary(
  t : ScviAgTape,
  a_id : Int32,
  b_id : Int32,
  sub : Bool,
  op : ScviAgOp,
) -> Int32 {
  let a = t.vals[a_id]
  let b = t.vals[b_id]
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let out = Array(capacity=n)
  for i in 0..<n {
    if sub {
      out.push(a[i] - b[i])
    } else {
      out.push(a[i] + b[i])
    }
  }
  scvi_ag_push(t, out, op)
}

///|
fn scvi_ag_unary(
  t : ScviAgTape,
  a : Int32,
  op : ScviAgOp,
  f : (Double) -> Double,
) -> Int32 {
  let av = t.vals[a]
  let out = Array(capacity=av.length())
  for i in 0..<av.length() {
    out.push(f(av[i]))
  }
  scvi_ag_push(t, out, op)
}

// ===========================================================================
// Tape construction and parameter registration
// ===========================================================================

///|
pub fn ScviAgTape::new() -> ScviAgTape {
  {
    vals: Array(capacity=0),
    grads: Array(capacity=0),
    ops: Array(capacity=0),
    mats: Array(capacity=0),
    vecs: Array(capacity=0),
  }
}

///|
/// Register a non-trainable leaf (input data, constant vector).
pub fn ScviAgTape::leaf(self : ScviAgTape, value : Array[Double]) -> Int32 {
  scvi_ag_push(self, value, AgLeaf)
}

///|
/// Register a trainable matrix parameter. The tape keeps a reference to the
/// supplied array, so in-place updates of `value` are visible to the model.
pub fn ScviAgTape::matrix_param(
  self : ScviAgTape,
  value : Array[Array[Double]],
) -> Int32 {
  let rows = value.length()
  let cols = if rows == 0 { 0 } else { value[0].length() }
  let id = self.mats.length()
  self.mats.push({ value, grad: scvi_ag_zero_matrix(rows, cols), })
  id
}

///|
/// Register a trainable vector parameter.
pub fn ScviAgTape::vector_param(
  self : ScviAgTape,
  value : Array[Double],
) -> Int32 {
  let id = self.vecs.length()
  self.vecs.push({ value, grad: Array::make(value.length(), 0.0), })
  id
}

///|
/// Number of registered matrix parameters.
pub fn ScviAgTape::n_matrices(self : ScviAgTape) -> Int32 {
  self.mats.length()
}

///|
/// Number of registered vector parameters.
pub fn ScviAgTape::n_vectors(self : ScviAgTape) -> Int32 {
  self.vecs.length()
}

///|
/// Read a node value.
pub fn ScviAgTape::value(self : ScviAgTape, id : Int32) -> Array[Double] {
  self.vals[id]
}

///|
/// Read a node gradient.
pub fn ScviAgTape::gradient(self : ScviAgTape, id : Int32) -> Array[Double] {
  self.grads[id]
}

///|
/// Read a matrix parameter gradient.
pub fn ScviAgTape::matrix_grad(
  self : ScviAgTape,
  id : Int32,
) -> Array[Array[Double]] {
  self.mats[id].grad
}

///|
/// Read a vector parameter gradient.
pub fn ScviAgTape::vector_grad(self : ScviAgTape, id : Int32) -> Array[Double] {
  self.vecs[id].grad
}

// ===========================================================================
// Forward operators
// ===========================================================================

///|
pub fn ScviAgTape::add(self : ScviAgTape, a : Int32, b : Int32) -> Int32 {
  scvi_ag_elementwise_binary(self, a, b, false, AgAdd(a, b))
}

///|
pub fn ScviAgTape::sub(self : ScviAgTape, a : Int32, b : Int32) -> Int32 {
  scvi_ag_elementwise_binary(self, a, b, true, AgSub(a, b))
}

///|
pub fn ScviAgTape::mul(self : ScviAgTape, a : Int32, b : Int32) -> Int32 {
  let av = self.vals[a]
  let bv = self.vals[b]
  let n = if av.length() < bv.length() { av.length() } else { bv.length() }
  let out = Array(capacity=n)
  for i in 0..<n {
    out.push(av[i] * bv[i])
  }
  scvi_ag_push(self, out, AgMul(a, b))
}

///|
pub fn ScviAgTape::scale(self : ScviAgTape, a : Int32, s : Double) -> Int32 {
  let av = self.vals[a]
  let out = Array(capacity=av.length())
  for i in 0..<av.length() {
    out.push(av[i] * s)
  }
  scvi_ag_push(self, out, AgScale(a, s))
}

///|
pub fn ScviAgTape::dense(
  self : ScviAgTape,
  x_id : Int32,
  w_id : Int32,
  b_id : Int32,
) -> Int32 {
  let x = self.vals[x_id]
  let w = self.mats[w_id].value
  let b = self.vecs[b_id].value
  let out = Array(capacity=w.length())
  for o in 0..<w.length() {
    let mut s = b[o]
    let row = w[o]
    let n_in = if row.length() < x.length() { row.length() } else { x.length() }
    for i in 0..<n_in {
      s += row[i] * x[i]
    }
    out.push(s)
  }
  scvi_ag_push(self, out, AgDense(x_id, w_id, b_id))
}

///|
pub fn ScviAgTape::relu(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgRelu(a), fn(x) { if x > 0.0 { x } else { 0.0 } })
}

///|
pub fn ScviAgTape::sigmoid(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgSigmoid(a), fn(x) {
    if x >= 0.0 {
      1.0 / (1.0 + @math.exp(-x))
    } else {
      let e = @math.exp(x)
      e / (1.0 + e)
    }
  })
}

///|
pub fn ScviAgTape::log(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgLog(a), fn(x) { @math.ln(scvi_dist_clamp(x)) })
}

///|
pub fn ScviAgTape::exp(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgExp(a), fn(x) { @math.exp(x) })
}

///|
pub fn ScviAgTape::square(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgSquare(a), fn(x) { x * x })
}

///|
pub fn ScviAgTape::neg(self : ScviAgTape, a : Int32) -> Int32 {
  scvi_ag_unary(self, a, AgNeg(a), fn(x) { -x })
}

///|
pub fn ScviAgTape::softmax(self : ScviAgTape, a : Int32) -> Int32 {
  let av = self.vals[a]
  let n = av.length()
  let mut m = if n == 0 { 0.0 } else { av[0] }
  for i in 1..<n {
    if av[i] > m {
      m = av[i]
    }
  }
  let out = Array(capacity=n)
  let mut denom = 0.0
  for i in 0..<n {
    let e = @math.exp(av[i] - m)
    out.push(e)
    denom += e
  }
  for i in 0..<n {
    out[i] = out[i] / denom
  }
  scvi_ag_push(self, out, AgSoftmax(a))
}

///|
pub fn ScviAgTape::sum(self : ScviAgTape, a : Int32) -> Int32 {
  let av = self.vals[a]
  let mut s = 0.0
  for i in 0..<av.length() {
    s += av[i]
  }
  scvi_ag_push(self, [s], AgSum(a))
}

///|
pub fn ScviAgTape::mean(self : ScviAgTape, a : Int32) -> Int32 {
  let av = self.vals[a]
  let mut s = 0.0
  for i in 0..<av.length() {
    s += av[i]
  }
  let n = if av.length() == 0 { 1.0 } else { av.length().to_double() }
  scvi_ag_push(self, [s / n], AgMean(a))
}

///|
pub fn ScviAgTape::select(self : ScviAgTape, a : Int32, k : Int32) -> Int32 {
  let av = self.vals[a]
  let idx = if k < 0 {
    0
  } else if k >= av.length() {
    av.length() - 1
  } else {
    k
  }
  scvi_ag_push(self, [av[idx]], AgSelect(a, idx))
}

///|
/// KL(q(z|x) || N(0, I)) for diagonal Gaussian q with mean `mu` and
/// log-variance `log_var`.
pub fn ScviAgTape::kl_divergence(
  self : ScviAgTape,
  mu : Int32,
  log_var : Int32,
) -> Int32 {
  let m = self.vals[mu]
  let lv = self.vals[log_var]
  let n = if m.length() < lv.length() { m.length() } else { lv.length() }
  let mut kl = 0.0
  for i in 0..<n {
    kl += 0.5 * (lv[i] - @math.exp(lv[i]) - m[i] * m[i] + 1.0)
  }
  scvi_ag_push(self, [kl], AgKl(mu, log_var))
}

///|
/// Negative binomial negative log-likelihood (up to parameter-independent
/// constants). `theta` is a fixed per-gene inverse dispersion vector.
pub fn ScviAgTape::nb_nll(
  self : ScviAgTape,
  mu : Int32,
  x : Array[Double],
  theta : Array[Double],
) -> Int32 {
  let m = self.vals[mu]
  let n = if m.length() < x.length() { m.length() } else { x.length() }
  let mut loss = 0.0
  for i in 0..<n {
    let mu_i = scvi_dist_clamp(m[i])
    let theta_i = scvi_dist_clamp(theta[i])
    loss += -(x[i] * @math.ln(mu_i) -
      mu_i -
      theta_i * @math.ln(theta_i + mu_i) +
      x[i] * @math.ln(theta_i + mu_i))
  }
  scvi_ag_push(self, [loss], AgNbNll(mu, x, theta))
}

///|
/// Poisson negative log-likelihood (up to a constant).
pub fn ScviAgTape::poisson_nll(
  self : ScviAgTape,
  mu : Int32,
  x : Array[Double],
) -> Int32 {
  let m = self.vals[mu]
  let n = if m.length() < x.length() { m.length() } else { x.length() }
  let mut loss = 0.0
  for i in 0..<n {
    let mu_i = scvi_dist_clamp(m[i])
    loss += mu_i - x[i] * @math.ln(mu_i)
  }
  scvi_ag_push(self, [loss], AgPoissonNll(mu, x))
}

// ===========================================================================
// Backward pass
// ===========================================================================

///|
/// Run reverse-mode differentiation starting from a scalar node.
pub fn ScviAgTape::backward(self : ScviAgTape, out_id : Int32) -> Unit {
  scvi_ag_zero_grads(self)
  let out_grad = self.grads[out_id]
  if out_grad.length() > 0 {
    out_grad[0] = 1.0
  }
  let n_nodes = self.ops.length()
  for step in 0..<n_nodes {
    let node = n_nodes - 1 - step
    let g = self.grads[node]
    match self.ops[node] {
      AgLeaf => ()
      AgAdd(a, b) => {
        scvi_ag_accum(self, a, g)
        scvi_ag_accum(self, b, g)
      }
      AgSub(a, b) => {
        scvi_ag_accum(self, a, g)
        let ng = Array(capacity=g.length())
        for i in 0..<g.length() {
          ng.push(-g[i])
        }
        scvi_ag_accum(self, b, ng)
      }
      AgMul(a, b) => {
        let av = self.vals[a]
        let bv = self.vals[b]
        let ga = Array(capacity=g.length())
        let gb = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(g[i] * bv[i])
          gb.push(g[i] * av[i])
        }
        scvi_ag_accum(self, a, ga)
        scvi_ag_accum(self, b, gb)
      }
      AgScale(a, s) => {
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(g[i] * s)
        }
        scvi_ag_accum(self, a, ga)
      }
      AgDense(x_id, w_id, b_id) => {
        let x = self.vals[x_id]
        let w = self.mats[w_id].value
        let wg = self.mats[w_id].grad
        let bg = self.vecs[b_id].grad
        let dx = Array::make(x.length(), 0.0)
        for o in 0..<g.length() {
          let go = g[o]
          let row = w[o]
          let n_in = if row.length() < x.length() {
            row.length()
          } else {
            x.length()
          }
          for i in 0..<n_in {
            dx[i] += row[i] * go
            wg[o][i] += go * x[i]
          }
          bg[o] += go
        }
        scvi_ag_accum(self, x_id, dx)
      }
      AgRelu(a) => {
        let av = self.vals[a]
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          if av[i] > 0.0 {
            ga.push(g[i])
          } else {
            ga.push(0.0)
          }
        }
        scvi_ag_accum(self, a, ga)
      }
      AgSigmoid(a) => {
        let av = self.vals[node]
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(g[i] * av[i] * (1.0 - av[i]))
        }
        scvi_ag_accum(self, a, ga)
      }
      AgSoftmax(a) => {
        let av = self.vals[node]
        let mut dot = 0.0
        for i in 0..<g.length() {
          dot += g[i] * av[i]
        }
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(av[i] * (g[i] - dot))
        }
        scvi_ag_accum(self, a, ga)
      }
      AgLog(a) => {
        let av = self.vals[a]
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(g[i] / scvi_dist_clamp(av[i]))
        }
        scvi_ag_accum(self, a, ga)
      }
      AgExp(a) => {
        let av = self.vals[node]
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(g[i] * av[i])
        }
        scvi_ag_accum(self, a, ga)
      }
      AgSquare(a) => {
        let av = self.vals[a]
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(2.0 * av[i] * g[i])
        }
        scvi_ag_accum(self, a, ga)
      }
      AgNeg(a) => {
        let ga = Array(capacity=g.length())
        for i in 0..<g.length() {
          ga.push(-g[i])
        }
        scvi_ag_accum(self, a, ga)
      }
      AgSum(a) => {
        let av = self.vals[a]
        let ga = Array::make(av.length(), g[0])
        scvi_ag_accum(self, a, ga)
      }
      AgMean(a) => {
        let av = self.vals[a]
        let n = if av.length() == 0 { 1.0 } else { av.length().to_double() }
        let ga = Array::make(av.length(), g[0] / n)
        scvi_ag_accum(self, a, ga)
      }
      AgSelect(a, k) => {
        let av = self.vals[a]
        let ga = Array::make(av.length(), 0.0)
        ga[k] = g[0]
        scvi_ag_accum(self, a, ga)
      }
      AgKl(mu, log_var) => {
        let m = self.vals[mu]
        let lv = self.vals[log_var]
        let n = if m.length() < lv.length() { m.length() } else { lv.length() }
        let gm = Array(capacity=n)
        let gl = Array(capacity=n)
        for i in 0..<n {
          gm.push(-m[i] * g[0])
          gl.push(0.5 * (1.0 - @math.exp(lv[i])) * g[0])
        }
        scvi_ag_accum(self, mu, gm)
        scvi_ag_accum(self, log_var, gl)
      }
      AgNbNll(mu, x, theta) => {
        let m = self.vals[mu]
        let n = if m.length() < x.length() { m.length() } else { x.length() }
        let gm = Array(capacity=n)
        for i in 0..<n {
          let mu_i = scvi_dist_clamp(m[i])
          let theta_i = scvi_dist_clamp(theta[i])
          let d = -(x[i] / mu_i - 1.0 + (x[i] - theta_i) / (theta_i + mu_i))
          gm.push(d * g[0])
        }
        scvi_ag_accum(self, mu, gm)
      }
      AgPoissonNll(mu, x) => {
        let m = self.vals[mu]
        let n = if m.length() < x.length() { m.length() } else { x.length() }
        let gm = Array(capacity=n)
        for i in 0..<n {
          let mu_i = scvi_dist_clamp(m[i])
          gm.push((1.0 - x[i] / mu_i) * g[0])
        }
        scvi_ag_accum(self, mu, gm)
      }
    }
  }
}

///|
/// Global gradient-norm clipping for all registered parameters.
pub fn ScviAgTape::clip_grad_norm(
  self : ScviAgTape,
  max_norm : Double,
) -> Double {
  let mut sq = 0.0
  for i in 0..<self.mats.length() {
    let g = self.mats[i].grad
    for r in 0..<g.length() {
      for c in 0..<g[r].length() {
        sq += g[r][c] * g[r][c]
      }
    }
  }
  for i in 0..<self.vecs.length() {
    let g = self.vecs[i].grad
    for j in 0..<g.length() {
      sq += g[j] * g[j]
    }
  }
  let norm = @math.pow(sq, 0.5)
  if norm > max_norm && norm > 0.0 {
    let factor = max_norm / norm
    for i in 0..<self.mats.length() {
      let g = self.mats[i].grad
      for r in 0..<g.length() {
        for c in 0..<g[r].length() {
          g[r][c] = g[r][c] * factor
        }
      }
    }
    for i in 0..<self.vecs.length() {
      let g = self.vecs[i].grad
      for j in 0..<g.length() {
        g[j] = g[j] * factor
      }
    }
  }
  norm
}
