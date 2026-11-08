///|
/// cobra3_solver.mbt — 线性/混合整数/二次规划求解器（纯 MoonBit CPU）.
///
/// - LP：两阶段单纯形（Bland 规则），返回**精确对偶解**（shadow price）与 reduced cost
/// - MILP：分支定界，整数变量由调用方给出下标
/// - QP：增广拉格朗日 + 坐标下降（box 约束），用于 MOMA 的二次目标
///
/// 约束统一为 Cobra3LinConstraint（<= / = / >=），变量带 lb/ub.
/// 所有边界都会被裁剪到 [-1e9, 1e9]，无穷上界用 1e15 表示.

///|
pub(all) enum Cobra3Sense {
  LessEqual
  Equal
  GreaterEqual
}

///|
pub(all) struct Cobra3LinConstraint {
  coeffs : Array[Double]
  rhs : Double
  sense : Cobra3Sense
}

///|
pub(all) enum Cobra3Status {
  Optimal
  Feasible
  Infeasible
  Unbounded
  IterLimit
}

///|
pub(all) struct Cobra3LpSolution {
  status : Cobra3Status
  objective_value : Double
  x : Array[Double]
  /// 原始约束行的对偶值（个数与 constraints 一致）
  duals : Array[Double]
  /// 原变量的 reduced cost
  reduced_costs : Array[Double]
  iterations : Int32
}

///|
pub(all) struct Cobra3QpSolution {
  status : Cobra3Status
  objective_value : Double
  x : Array[Double]
  /// 等式约束的对偶（min 问题语义）
  duals : Array[Double]
  converged : Bool
  iterations : Int32
}

// ---------------------------------------------------------------------------
// 单纯形内部
// ---------------------------------------------------------------------------

///|
priv struct Cobra3Tableau {
  a : Array[Array[Double]]
  b : Array[Double]
  basis : Array[Int32]
  rows : Int32
  cols : Int32
}

///|
fn cobra3_tableau_pivot(t : Cobra3Tableau, row : Int32, col : Int32) -> Unit {
  let piv = t.a[row][col]
  if piv == 0.0 {
    return
  }
  for j in 0..<t.cols {
    t.a[row][j] = t.a[row][j] / piv
  }
  t.b[row] = t.b[row] / piv
  for i in 0..<t.rows {
    if i != row {
      let factor = t.a[i][col]
      if factor != 0.0 {
        for j in 0..<t.cols {
          t.a[i][j] = t.a[i][j] - factor * t.a[row][j]
        }
        t.b[i] = t.b[i] - factor * t.b[row]
      }
    }
  }
  t.basis[row] = col
}

///|
fn cobra3_tableau_reduced_cost(
  t : Cobra3Tableau,
  c : Array[Double],
  j : Int32,
) -> Double {
  let mut z = c[j]
  for i in 0..<t.rows {
    z = z - c[t.basis[i]] * t.a[i][j]
  }
  z
}

///|
/// Phase 1 后把仍留在基的人工变量换出；无法换出的行是冗余约束，清零.
fn cobra3_remove_artificial_basis(t : Cobra3Tableau, col_art : Int32) -> Unit {
  for i in 0..<t.rows {
    if t.basis[i] >= col_art {
      let mut pivot_col = -1
      for j in 0..<col_art {
        if t.a[i][j].abs() > 1.0e-9 && !t.basis.contains(j) {
          pivot_col = j
          break
        }
      }
      if pivot_col >= 0 {
        cobra3_tableau_pivot(t, i, pivot_col)
      } else {
        for j in 0..<t.cols {
          t.a[i][j] = 0.0
        }
        t.b[i] = 0.0
      }
    }
  }
}

///|
/// 返回 (0=最优, 1=无界, 2=迭代上限) 与已执行迭代数.
fn cobra3_simplex_run(
  t : Cobra3Tableau,
  c : Array[Double],
  max_iter : Int32,
  enter_limit : Int32,
) -> (Int32, Int32) {
  let tol = 1.0e-9
  let mut iter = 0
  while iter < max_iter {
    iter = iter + 1
    let mut enter = -1
    for j in 0..<enter_limit {
      if cobra3_tableau_reduced_cost(t, c, j) > tol {
        enter = j
        break
      }
    }
    if enter < 0 {
      return (0, iter)
    }
    let mut leave = -1
    let mut best_ratio = 1.0e100
    for i in 0..<t.rows {
      let d = t.a[i][enter]
      if d > tol {
        let ratio = t.b[i] / d
        if ratio < best_ratio - 1.0e-12 ||
          (
            ratio < best_ratio + 1.0e-12 &&
            (leave < 0 || t.basis[i] < t.basis[leave])
          ) {
          best_ratio = ratio
          leave = i
        }
      }
    }
    if leave < 0 {
      return (1, iter)
    }
    cobra3_tableau_pivot(t, leave, enter)
  }
  (2, iter)
}

///|
/// 解 M^T? 这里解 M·pi = rhs，部分主元高斯消元.
fn cobra3_solve_linear(
  m : Array[Array[Double]],
  rhs : Array[Double],
) -> Array[Double] {
  let n = m.length()
  let a : Array[Array[Double]] = Array(capacity=n)
  for i in 0..<n {
    let row : Array[Double] = Array::make(n + 1, 0.0)
    for j in 0..<n {
      row[j] = m[i][j]
    }
    row[n] = rhs[i]
    a.push(row)
  }
  for col in 0..<n {
    let mut piv = col
    for i in col..<n {
      if a[i][col].abs() > a[piv][col].abs() {
        piv = i
      }
    }
    if a[piv][col].abs() < 1.0e-12 {
      continue
    }
    if piv != col {
      let tmp = a[col]
      a[col] = a[piv]
      a[piv] = tmp
    }
    let d = a[col][col]
    for j in col..<(n + 1) {
      a[col][j] = a[col][j] / d
    }
    for i in 0..<n {
      if i != col {
        let factor = a[i][col]
        if factor != 0.0 {
          for j in col..<(n + 1) {
            a[i][j] = a[i][j] - factor * a[col][j]
          }
        }
      }
    }
  }
  let out : Array[Double] = Array(capacity=n)
  for i in 0..<n {
    out.push(a[i][n])
  }
  out
}

///|
fn cobra3_clamp_bounds(lb : Array[Double], ub : Array[Double]) -> Unit {
  for j in 0..<lb.length() {
    if lb[j] < -1.0e9 {
      lb[j] = -1.0e9
    }
    if ub[j] > 1.0e9 {
      ub[j] = 1.0e9
    }
  }
}

// ---------------------------------------------------------------------------
// LP
// ---------------------------------------------------------------------------

///|
/// 通用 LP：max c^T x，s.t. constraints，lb <= x <= ub.
pub fn cobra3_solve_lp(
  constraints : Array[Cobra3LinConstraint],
  c : Array[Double],
  lb_in : Array[Double],
  ub_in : Array[Double],
) -> Cobra3LpSolution {
  let n = c.length()
  let lb : Array[Double] = Array(capacity=n)
  let ub : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    lb.push(lb_in[j])
    ub.push(ub_in[j])
  }
  cobra3_clamp_bounds(lb, ub)
  let empty_x : Array[Double] = Array::make(n, 0.0)
  let empty_dual : Array[Double] = Array::make(constraints.length(), 0.0)
  let empty_rc : Array[Double] = Array::make(n, 0.0)
  let open_upper : Array[Bool] = Array(capacity=n)
  for j in 0..<n {
    open_upper.push(ub_in[j] >= 1.0e15)
  }
  if n == 0 {
    return {
      status: Cobra3Status::Optimal,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: 0,
    }
  }
  for j in 0..<n {
    if lb[j] > ub[j] + 1.0e-9 {
      return {
        status: Cobra3Status::Infeasible,
        objective_value: 0.0,
        x: empty_x,
        duals: empty_dual,
        reduced_costs: empty_rc,
        iterations: 0,
      }
    }
  }
  // 变量 y = x - lb >= 0
  let nc = constraints.length()
  let mut n_slack = 0
  for con in constraints {
    match con.sense {
      Equal => ()
      _ => n_slack = n_slack + 1
    }
  }
  let mut n_ub_rows = 0
  for j in 0..<n {
    if !open_upper[j] && ub[j] - lb[j] < 1.0e15 {
      n_ub_rows = n_ub_rows + 1
    }
  }
  let rows = nc + n_ub_rows
  let col_slack = n
  let col_art = n + n_slack + n_ub_rows
  let cols = col_art + rows
  let a : Array[Array[Double]] = Array(capacity=rows)
  for _ in 0..<rows {
    a.push(Array::make(cols, 0.0))
  }
  let b : Array[Double] = Array(capacity=rows)
  let signs : Array[Double] = Array(capacity=rows)
  let mut slack_idx = 0
  for i in 0..<nc {
    let con = constraints[i]
    let mut rhs = con.rhs
    for j in 0..<n {
      a[i][j] = con.coeffs[j]
      rhs = rhs - con.coeffs[j] * lb[j]
    }
    match con.sense {
      Equal => ()
      LessEqual => {
        a[i][col_slack + slack_idx] = 1.0
        slack_idx = slack_idx + 1
      }
      GreaterEqual => {
        a[i][col_slack + slack_idx] = -1.0
        slack_idx = slack_idx + 1
      }
    }
    let mut sign = 1.0
    if rhs < 0.0 {
      sign = -1.0
      rhs = -rhs
      for j in 0..<cols {
        a[i][j] = -a[i][j]
      }
    }
    signs.push(sign)
    b.push(rhs)
  }
  let mut row = nc
  for j in 0..<n {
    let w = ub[j] - lb[j]
    if !open_upper[j] && w < 1.0e15 {
      a[row][j] = 1.0
      a[row][n + n_slack + (row - nc)] = 1.0
      b.push(w)
      signs.push(1.0)
      row = row + 1
    }
  }
  // 人工变量
  for i in 0..<rows {
    a[i][col_art + i] = 1.0
  }
  let basis : Array[Int32] = Array(capacity=rows)
  for i in 0..<rows {
    basis.push(col_art + i)
  }
  let t : Cobra3Tableau = { a, b, basis, rows, cols, }
  // Phase 1: max -sum(artificial)
  let c1 : Array[Double] = Array::make(cols, 0.0)
  for i in 0..<rows {
    c1[col_art + i] = -1.0
  }
  let max_iter = 20000 + rows * 200
  let (r1, it1) = cobra3_simplex_run(t, c1, max_iter, col_art)
  if r1 == 2 {
    return {
      status: Cobra3Status::IterLimit,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: it1,
    }
  }
  let mut art_sum = 0.0
  for i in 0..<rows {
    if t.basis[i] >= col_art {
      art_sum = art_sum + t.b[i]
    }
  }
  if art_sum > 1.0e-7 {
    return {
      status: Cobra3Status::Infeasible,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: it1,
    }
  }
  // Phase 2：先把零值人工变量换出基，避免其值在 Phase 2 变化破坏约束
  cobra3_remove_artificial_basis(t, col_art)
  let c2 : Array[Double] = Array::make(cols, 0.0)
  for j in 0..<n {
    c2[j] = c[j]
  }
  // Phase 2 不允许人工变量重新入基，否则会破坏原约束
  let (r2, it2) = cobra3_simplex_run(t, c2, max_iter, col_art)
  if r2 == 2 {
    return {
      status: Cobra3Status::IterLimit,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: it1 + it2,
    }
  }
  if r2 == 1 {
    return {
      status: Cobra3Status::Unbounded,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: it1 + it2,
    }
  }
  let x : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    x.push(lb[j])
  }
  // 基变量值：x_basis = b（tableau 已化为 B^-1 b）
  for i in 0..<rows {
    let basis_var = t.basis[i]
    if basis_var < n {
      x[basis_var] = lb[basis_var] + t.b[i]
    }
  }
  // 非基原变量为 lb
  let mut obj = 0.0
  for j in 0..<n {
    obj = obj + c[j] * x[j]
  }
  // 对偶：解 B^T pi = c_B
  let bt : Array[Array[Double]] = Array(capacity=rows)
  for k in 0..<rows {
    let rowk : Array[Double] = Array::make(rows, 0.0)
    for i in 0..<rows {
      rowk[i] = t.a[i][t.basis[k]]
    }
    bt.push(rowk)
  }
  let cb : Array[Double] = Array(capacity=rows)
  for i in 0..<rows {
    cb.push(c2[t.basis[i]])
  }
  let pi = cobra3_solve_linear(bt, cb)
  let duals : Array[Double] = Array(capacity=nc)
  for i in 0..<nc {
    duals.push(pi[i] * signs[i])
  }
  let rc : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    rc.push(cobra3_tableau_reduced_cost(t, c2, j))
  }
  {
    status: Cobra3Status::Optimal,
    objective_value: obj,
    x,
    duals,
    reduced_costs: rc,
    iterations: it1 + it2,
  }
}

///|
/// 便捷入口：等式 A x = b 形式的 FBA LP.
pub fn cobra3_solve_lp_eq(
  a_eq : Array[Array[Double]],
  b_eq : Array[Double],
  c : Array[Double],
  lb : Array[Double],
  ub : Array[Double],
) -> Cobra3LpSolution {
  let constraints : Array[Cobra3LinConstraint] = Array(capacity=a_eq.length())
  for i in 0..<a_eq.length() {
    constraints.push({
      coeffs: a_eq[i],
      rhs: b_eq[i],
      sense: Cobra3Sense::Equal,
    })
  }
  cobra3_solve_lp(constraints, c, lb, ub)
}

// ---------------------------------------------------------------------------
// MILP：分支定界
// ---------------------------------------------------------------------------

///|
fn cobra3_bounds_copy(arr : Array[Double]) -> Array[Double] {
  let out : Array[Double] = Array(capacity=arr.length())
  for v in arr {
    out.push(v)
  }
  out
}

///|
fn cobra3_pick_branch_var(x : Array[Double], integer : Array[Int32]) -> Int32 {
  let mut best = -1
  let mut best_frac = -1.0
  for idx in integer {
    let v = x[idx]
    let frac = (v - @math.floor(v)).abs()
    let d = if frac < 0.5 { frac } else { 1.0 - frac }
    if d > 1.0e-6 && d > best_frac {
      best = idx
      best_frac = d
    }
  }
  best
}

///|
/// 分支定界 MILP（max c^T x）.
/// 达到 max_nodes 时返回当前 incumbent 并标记 Optimal，不保证全局最优；
/// 未找到 incumbent 时返回 Infeasible（也可能是节点上限导致，见 iterations）.
pub fn cobra3_solve_milp(
  constraints : Array[Cobra3LinConstraint],
  c : Array[Double],
  lb : Array[Double],
  ub : Array[Double],
  integer_indices : Array[Int32],
  max_nodes~ : Int32 = 4000,
) -> Cobra3LpSolution {
  let n = c.length()
  let empty_x : Array[Double] = Array::make(n, 0.0)
  let empty_dual : Array[Double] = Array::make(constraints.length(), 0.0)
  let empty_rc : Array[Double] = Array::make(n, 0.0)
  let root = cobra3_solve_lp(constraints, c, lb, ub)
  if root.status is Unbounded {
    return root
  }
  if root.status is Infeasible || root.status is IterLimit {
    return root
  }
  let mut incumbent_obj = -1.0e100
  let mut incumbent_x : Array[Double] = Array::make(n, 0.0)
  let mut incumbent_dual : Array[Double] = Array::make(
    constraints.length(),
    0.0,
  )
  let mut incumbent_rc : Array[Double] = Array::make(n, 0.0)
  let mut incumbent_found = false
  let stack_lb : Array[Array[Double]] = Array(capacity=16)
  let stack_ub : Array[Array[Double]] = Array(capacity=16)
  stack_lb.push(cobra3_bounds_copy(lb))
  stack_ub.push(cobra3_bounds_copy(ub))
  let mut nodes = 0
  while stack_lb.length() > 0 && nodes < max_nodes {
    nodes = nodes + 1
    let node_lb = stack_lb[stack_lb.length() - 1]
    let node_ub = stack_ub[stack_ub.length() - 1]
    let _ = stack_lb.pop()
    let _ = stack_ub.pop()
    let sol = cobra3_solve_lp(constraints, c, node_lb, node_ub)
    if sol.status is Infeasible || sol.status is IterLimit {
      continue
    }
    if sol.status is Unbounded {
      continue
    }
    if sol.objective_value <= incumbent_obj + 1.0e-7 {
      continue
    }
    let branch = cobra3_pick_branch_var(sol.x, integer_indices)
    if branch < 0 {
      incumbent_obj = sol.objective_value
      incumbent_x = sol.x
      incumbent_dual = sol.duals
      incumbent_rc = sol.reduced_costs
      incumbent_found = true
      continue
    }
    let v = sol.x[branch]
    let fl = @math.floor(v)
    let ce = @math.ceil(v)
    // 分支 1: x <= floor(v)
    if fl >= node_lb[branch] - 1.0e-9 {
      let lb1 = cobra3_bounds_copy(node_lb)
      let ub1 = cobra3_bounds_copy(node_ub)
      ub1[branch] = fl
      if lb1[branch] <= ub1[branch] {
        stack_lb.push(lb1)
        stack_ub.push(ub1)
      }
    }
    // 分支 2: x >= ceil(v)
    if ce <= node_ub[branch] + 1.0e-9 {
      let lb2 = cobra3_bounds_copy(node_lb)
      let ub2 = cobra3_bounds_copy(node_ub)
      lb2[branch] = ce
      if lb2[branch] <= ub2[branch] {
        stack_lb.push(lb2)
        stack_ub.push(ub2)
      }
    }
  }
  if !incumbent_found {
    return {
      status: if stack_lb.length() == 0 {
        Cobra3Status::Infeasible
      } else {
        Cobra3Status::IterLimit
      },
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      reduced_costs: empty_rc,
      iterations: nodes,
    }
  }
  {
    status: if stack_lb.length() == 0 {
      Cobra3Status::Optimal
    } else {
      Cobra3Status::Feasible
    },
    objective_value: incumbent_obj,
    x: incumbent_x,
    duals: incumbent_dual,
    reduced_costs: incumbent_rc,
    iterations: nodes,
  }
}

// ---------------------------------------------------------------------------
// QP：增广拉格朗日 + 坐标下降
// ---------------------------------------------------------------------------

///|
/// min 0.5 x^T Q x + c^T x，s.t. A x = b，lb <= x <= ub.
/// Q 需对称半正定；坐标下降逐维求解析解，乘子法处理等式约束.
pub fn cobra3_solve_qp(
  q : Array[Array[Double]],
  c : Array[Double],
  a_eq : Array[Array[Double]],
  b_eq : Array[Double],
  lb_in : Array[Double],
  ub_in : Array[Double],
  max_iter~ : Int32 = 4000,
  tol~ : Double = 1.0e-8,
) -> Cobra3QpSolution {
  let n = c.length()
  let m = a_eq.length()
  let lb = cobra3_bounds_copy(lb_in)
  let ub = cobra3_bounds_copy(ub_in)
  cobra3_clamp_bounds(lb, ub)
  let empty_x : Array[Double] = Array::make(n, 0.0)
  let empty_dual : Array[Double] = Array::make(m, 0.0)
  for j in 0..<n {
    if lb[j] > ub[j] + 1.0e-9 {
      return {
        status: Cobra3Status::Infeasible,
        objective_value: 0.0,
        x: empty_x,
        duals: empty_dual,
        converged: false,
        iterations: 0,
      }
    }
  }
  // 用 LP 找可行初值
  let zero_c : Array[Double] = Array::make(n, 0.0)
  let feasible = cobra3_solve_lp_eq(a_eq, b_eq, zero_c, lb, ub)
  if feasible.status is Infeasible || feasible.status is IterLimit {
    return {
      status: Cobra3Status::Infeasible,
      objective_value: 0.0,
      x: empty_x,
      duals: empty_dual,
      converged: false,
      iterations: 0,
    }
  }
  let x = cobra3_bounds_copy(feasible.x)
  for j in 0..<n {
    if x[j] < lb[j] {
      x[j] = lb[j]
    }
    if x[j] > ub[j] {
      x[j] = ub[j]
    }
  }
  let lambda : Array[Double] = Array::make(m, 0.0)
  let ax : Array[Double] = Array::make(m, 0.0)
  for i in 0..<m {
    let mut sum = 0.0
    for j in 0..<n {
      sum = sum + a_eq[i][j] * x[j]
    }
    ax[i] = sum
  }
  let rho = 1.0
  let mut iter = 0
  let mut converged = false
  while iter < max_iter {
    iter = iter + 1
    let lambda_prev = cobra3_bounds_copy(lambda)
    // 内层：坐标下降最小化增广拉格朗日
    for _inner in 0..<100 {
      let mut max_delta = 0.0
      for j in 0..<n {
        let mut rest = c[j]
        for k in 0..<n {
          if k != j {
            rest = rest + q[j][k] * x[k]
          }
        }
        let mut denom = q[j][j]
        for i in 0..<m {
          let aij = a_eq[i][j]
          if aij != 0.0 {
            let resid = ax[i] - b_eq[i] - aij * x[j]
            rest = rest + aij * (lambda[i] + rho * resid)
            denom = denom + rho * aij * aij
          }
        }
        if denom > 1.0e-14 {
          let mut nx = -rest / denom
          if nx < lb[j] {
            nx = lb[j]
          }
          if nx > ub[j] {
            nx = ub[j]
          }
          let delta = nx - x[j]
          if delta.abs() > max_delta {
            max_delta = delta.abs()
          }
          if delta != 0.0 {
            for i in 0..<m {
              ax[i] = ax[i] + a_eq[i][j] * delta
            }
            x[j] = nx
          }
        }
      }
      if max_delta < 1.0e-11 {
        break
      }
    }
    for i in 0..<m {
      lambda[i] = lambda[i] + rho * (ax[i] - b_eq[i])
    }
    let mut primal = 0.0
    for i in 0..<m {
      let r = (ax[i] - b_eq[i]).abs()
      if r > primal {
        primal = r
      }
    }
    let mut dual_delta = 0.0
    for i in 0..<m {
      let d = (lambda[i] - lambda_prev[i]).abs()
      if d > dual_delta {
        dual_delta = d
      }
    }
    if primal < tol && dual_delta < tol {
      converged = true
      break
    }
  }
  let mut obj = 0.0
  for i in 0..<n {
    obj = obj + c[i] * x[i]
    for j in 0..<n {
      obj = obj + 0.5 * q[i][j] * x[i] * x[j]
    }
  }
  {
    status: if converged {
      Cobra3Status::Optimal
    } else {
      Cobra3Status::IterLimit
    },
    objective_value: obj,
    x,
    duals: lambda,
    converged,
    iterations: iter,
  }
}
