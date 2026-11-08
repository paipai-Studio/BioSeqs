// scvi_optim.mbt - Optimisers for the scVI autograd engine
// Mirrors the role of torch.optim in scvi-tools: SGD, momentum SGD and Adam.

///|
fn scvi_optim_zero_matrix(rows : Int32, cols : Int32) -> Array[Array[Double]] {
  let out = Array(capacity=rows)
  for _i in 0..<rows {
    out.push(Array::make(cols, 0.0))
  }
  out
}

///|
pub struct ScviSgdState {
  lr : Double
  weight_decay : Double
}

///|
pub fn scvi_sgd_init(lr : Double) -> ScviSgdState {
  { lr, weight_decay: 0.0, }
}

///|
/// Plain stochastic gradient descent step (in place).
pub fn scvi_sgd_step(tape : ScviAgTape, state : ScviSgdState) -> Unit {
  for i in 0..<tape.mats.length() {
    let param = tape.mats[i]
    let value = param.value
    let grad = param.grad
    for r in 0..<value.length() {
      for c in 0..<value[r].length() {
        let decay = state.weight_decay * value[r][c]
        value[r][c] = value[r][c] - state.lr * (grad[r][c] + decay)
      }
    }
  }
  for i in 0..<tape.vecs.length() {
    let param = tape.vecs[i]
    let value = param.value
    let grad = param.grad
    for j in 0..<value.length() {
      let decay = state.weight_decay * value[j]
      value[j] = value[j] - state.lr * (grad[j] + decay)
    }
  }
}

///|
pub struct ScviMomentumState {
  lr : Double
  momentum : Double
  mats_velocity : Array[Array[Array[Double]]]
  vecs_velocity : Array[Array[Double]]
}

///|
pub fn scvi_momentum_init(
  tape : ScviAgTape,
  lr : Double,
  momentum : Double,
) -> ScviMomentumState {
  let mats = Array(capacity=tape.mats.length())
  for i in 0..<tape.mats.length() {
    let shape = tape.mats[i].value
    let rows = shape.length()
    let cols = if rows == 0 { 0 } else { shape[0].length() }
    mats.push(scvi_optim_zero_matrix(rows, cols))
  }
  let vecs = Array(capacity=tape.vecs.length())
  for i in 0..<tape.vecs.length() {
    vecs.push(Array::make(tape.vecs[i].value.length(), 0.0))
  }
  { lr, momentum, mats_velocity: mats, vecs_velocity: vecs, }
}

///|
pub fn scvi_momentum_step(tape : ScviAgTape, state : ScviMomentumState) -> Unit {
  for i in 0..<tape.mats.length() {
    let param = tape.mats[i]
    let value = param.value
    let grad = param.grad
    let v = state.mats_velocity[i]
    for r in 0..<value.length() {
      for c in 0..<value[r].length() {
        v[r][c] = state.momentum * v[r][c] + grad[r][c]
        value[r][c] = value[r][c] - state.lr * v[r][c]
      }
    }
  }
  for i in 0..<tape.vecs.length() {
    let param = tape.vecs[i]
    let value = param.value
    let grad = param.grad
    let v = state.vecs_velocity[i]
    for j in 0..<value.length() {
      v[j] = state.momentum * v[j] + grad[j]
      value[j] = value[j] - state.lr * v[j]
    }
  }
}

///|
pub struct ScviAdamState {
  lr : Double
  beta1 : Double
  beta2 : Double
  eps : Double
  weight_decay : Double
  max_grad_norm : Double
  mut step_count : Int32
  mats_m : Array[Array[Array[Double]]]
  mats_v : Array[Array[Array[Double]]]
  vecs_m : Array[Array[Double]]
  vecs_v : Array[Array[Double]]
}

///|
pub fn scvi_adam_init(tape : ScviAgTape, lr : Double) -> ScviAdamState {
  let mats_m = Array(capacity=tape.mats.length())
  let mats_v = Array(capacity=tape.mats.length())
  for i in 0..<tape.mats.length() {
    let shape = tape.mats[i].value
    let rows = shape.length()
    let cols = if rows == 0 { 0 } else { shape[0].length() }
    mats_m.push(scvi_optim_zero_matrix(rows, cols))
    mats_v.push(scvi_optim_zero_matrix(rows, cols))
  }
  let vecs_m = Array(capacity=tape.vecs.length())
  let vecs_v = Array(capacity=tape.vecs.length())
  for i in 0..<tape.vecs.length() {
    let n = tape.vecs[i].value.length()
    vecs_m.push(Array::make(n, 0.0))
    vecs_v.push(Array::make(n, 0.0))
  }
  {
    lr,
    beta1: 0.9,
    beta2: 0.999,
    eps: 1.0e-8,
    weight_decay: 0.0,
    max_grad_norm: 0.0,
    step_count: 0,
    mats_m,
    mats_v,
    vecs_m,
    vecs_v,
  }
}

///|
/// Adam update (in place). Call after `tape.backward(loss_id)`.
pub fn scvi_adam_step(tape : ScviAgTape, state : ScviAdamState) -> Unit {
  state.step_count = state.step_count + 1
  let bias1 = 1.0 - @math.pow(state.beta1, state.step_count.to_double())
  let bias2 = 1.0 - @math.pow(state.beta2, state.step_count.to_double())
  for i in 0..<tape.mats.length() {
    let param = tape.mats[i]
    let value = param.value
    let grad = param.grad
    let m = state.mats_m[i]
    let v = state.mats_v[i]
    for r in 0..<value.length() {
      for c in 0..<value[r].length() {
        let g = grad[r][c]
        m[r][c] = state.beta1 * m[r][c] + (1.0 - state.beta1) * g
        v[r][c] = state.beta2 * v[r][c] + (1.0 - state.beta2) * g * g
        let m_hat = m[r][c] / bias1
        let v_hat = v[r][c] / bias2
        let update = state.lr * m_hat / (@math.pow(v_hat, 0.5) + state.eps)
        let decay = state.weight_decay * value[r][c]
        value[r][c] = value[r][c] - update - state.lr * decay
      }
    }
  }
  for i in 0..<tape.vecs.length() {
    let param = tape.vecs[i]
    let value = param.value
    let grad = param.grad
    let m = state.vecs_m[i]
    let v = state.vecs_v[i]
    for j in 0..<value.length() {
      let g = grad[j]
      m[j] = state.beta1 * m[j] + (1.0 - state.beta1) * g
      v[j] = state.beta2 * v[j] + (1.0 - state.beta2) * g * g
      let m_hat = m[j] / bias1
      let v_hat = v[j] / bias2
      let update = state.lr * m_hat / (@math.pow(v_hat, 0.5) + state.eps)
      let decay = state.weight_decay * value[j]
      value[j] = value[j] - update - state.lr * decay
    }
  }
}

///|
/// Project a scalar loss value out of a tape.
pub fn scvi_ag_loss(tape : ScviAgTape, loss_id : Int32) -> Double {
  let v = tape.vals[loss_id]
  if v.length() > 0 {
    v[0]
  } else {
    0.0
  }
}
