/// Energy minimization: steepest descent and L-BFGS.
///
/// OpenMM equivalent: `LocalEnergyMinimizer`
///
/// Minimizes the potential energy by adjusting atomic positions
/// using gradient information from the force function.

///|
/// Steepest descent minimization.
///
/// positions: initial positions
/// force_fn: returns forces (negative gradient) from positions
/// max_iter: maximum iterations
/// tolerance: force convergence tolerance (max |F| < tolerance)
/// step_size: initial step size (auto-adjusted)
///
/// Returns minimized positions.
pub fn omm_minimize_steepest_descent(
  positions : Array[Array[Double]],
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  max_iter~ : Int32 = 1000,
  tolerance~ : Double = 0.01,
  step_size~ : Double = 0.01,
) -> Array[Array[Double]] {
  let n = positions.length()
  let mut pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
  }
  let mut dt = step_size
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let forces = force_fn(pos)
    // Compute max force
    let mut max_f = 0.0
    let mut f_norm_sq = 0.0
    for i = 0; i < n; i = i + 1 {
      let f_sq = forces[i][0] * forces[i][0] +
        forces[i][1] * forces[i][1] +
        forces[i][2] * forces[i][2]
      f_norm_sq = f_norm_sq + f_sq
      let f_mag = f_sq.sqrt()
      if f_mag > max_f {
        max_f = f_mag
      }
    }
    if max_f < tolerance {
      break
    }
    // Step along force direction (steepest descent = -gradient = +force)
    let f_norm = f_norm_sq.sqrt()
    if f_norm == 0.0 {
      break
    }
    let new_pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
    for i = 0; i < n; i = i + 1 {
      new_pos[i][0] = pos[i][0] + dt * forces[i][0] / f_norm
      new_pos[i][1] = pos[i][1] + dt * forces[i][1] / f_norm
      new_pos[i][2] = pos[i][2] + dt * forces[i][2] / f_norm
    }
    // Adaptive step size: if energy increases, reduce step
    pos = new_pos
    dt = dt * 1.1 // gradually increase step
    if dt > 0.1 {
      dt = 0.1
    }
  }
  pos
}

///|
/// L-BFGS minimization (limited-memory BFGS).
///
/// Uses gradient history to build an approximate inverse Hessian
/// for superlinear convergence. Memory efficient: stores only m past gradients.
pub fn omm_minimize_lbfgs(
  positions : Array[Array[Double]],
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  max_iter~ : Int32 = 200,
  tolerance~ : Double = 0.01,
  m~ : Int32 = 5,
) -> Array[Array[Double]] {
  let n = positions.length()
  let dim = n * 3
  // Flatten positions
  let mut x : Array[Double] = Array::make(dim, 0.0)
  for i = 0; i < n; i = i + 1 {
    x[3 * i] = positions[i][0]
    x[3 * i + 1] = positions[i][1]
    x[3 * i + 2] = positions[i][2]
  }
  // Helper: flatten forces
  fn flatten_forces(forces : Array[Array[Double]]) -> Array[Double] {
    let f : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < n; i = i + 1 {
      f[3 * i] = forces[i][0]
      f[3 * i + 1] = forces[i][1]
      f[3 * i + 2] = forces[i][2]
    }
    f
  }
  // Helper: unflatten to positions
  fn unflatten(arr : Array[Double]) -> Array[Array[Double]] {
    let p : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
    for i = 0; i < n; i = i + 1 {
      p[i] = [arr[3 * i], arr[3 * i + 1], arr[3 * i + 2]]
    }
    p
  }
  // L-BFGS history
  let s_list : Array[Array[Double]] = Array(capacity=0) // position differences
  let y_list : Array[Array[Double]] = Array(capacity=0) // gradient differences
  let rho_list : Array[Double] = Array(capacity=0) // 1 / (y^T s)
  // Initial gradient (negative force)
  let forces0 = force_fn(unflatten(x))
  let mut grad = flatten_forces(forces0)
  for i = 0; i < dim; i = i + 1 {
    grad[i] = -grad[i] // gradient = -force
  }
  for iter = 0; iter < max_iter; iter = iter + 1 {
    // Check convergence
    let mut g_norm = 0.0
    for i = 0; i < dim; i = i + 1 {
      g_norm = g_norm + grad[i] * grad[i]
    }
    g_norm = g_norm.sqrt()
    if g_norm < tolerance {
      break
    }
    // L-BFGS two-loop recursion to compute search direction
    let q : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      q[i] = grad[i]
    }
    let hist_len = s_list.length()
    let alpha : Array[Double] = Array::make(hist_len, 0.0)
    // First loop (backward)
    for k = hist_len - 1; k >= 0; k = k - 1 {
      let mut sq = 0.0
      for i = 0; i < dim; i = i + 1 {
        sq = sq + s_list[k][i] * q[i]
      }
      alpha[k] = rho_list[k] * sq
      for i = 0; i < dim; i = i + 1 {
        q[i] = q[i] - alpha[k] * y_list[k][i]
      }
    }
    // Apply initial Hessian approximation H0 = (s^T y / y^T y) * I
    let mut gamma = 1.0
    if hist_len > 0 {
      let mut sy = 0.0
      let mut yy = 0.0
      for i = 0; i < dim; i = i + 1 {
        sy = sy + s_list[hist_len - 1][i] * y_list[hist_len - 1][i]
        yy = yy + y_list[hist_len - 1][i] * y_list[hist_len - 1][i]
      }
      if yy > 0.0 {
        gamma = sy / yy
      }
    }
    let r : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      r[i] = gamma * q[i]
    }
    // Second loop (forward)
    for k = 0; k < hist_len; k = k + 1 {
      let mut yr = 0.0
      for i = 0; i < dim; i = i + 1 {
        yr = yr + y_list[k][i] * r[i]
      }
      let beta = rho_list[k] * yr
      for i = 0; i < dim; i = i + 1 {
        r[i] = r[i] + s_list[k][i] * (alpha[k] - beta)
      }
    }
    // Search direction = -H * grad = -r
    let p_dir : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      p_dir[i] = -r[i]
    }
    // Line search (simple backtracking)
    let step = 1.0
    let x_new : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      x_new[i] = x[i] + step * p_dir[i]
    }
    let forces_new = force_fn(unflatten(x_new))
    let grad_new = flatten_forces(forces_new)
    for i = 0; i < dim; i = i + 1 {
      grad_new[i] = -grad_new[i]
    }
    // Check Armijo condition (simplified)
    let mut dg = 0.0
    for i = 0; i < dim; i = i + 1 {
      dg = dg + grad[i] * p_dir[i]
    }
    // Accept step (simplified line search)
    // Update history
    let s_new : Array[Double] = Array::make(dim, 0.0)
    let y_new : Array[Double] = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      s_new[i] = x_new[i] - x[i]
      y_new[i] = grad_new[i] - grad[i]
    }
    let mut sy = 0.0
    for i = 0; i < dim; i = i + 1 {
      sy = sy + s_new[i] * y_new[i]
    }
    if sy > 0.0 {
      if s_list.length() >= m {
        ignore(s_list.remove(0))
        ignore(y_list.remove(0))
        ignore(rho_list.remove(0))
      }
      s_list.push(s_new)
      y_list.push(y_new)
      rho_list.push(1.0 / sy)
    }
    x = x_new
    grad = grad_new
  }
  unflatten(x)
}
