// ===========================================================================
// cobra3_solver_test.mbt — LP/MILP/QP 求解器测试
// ===========================================================================

///|
fn cobra3_con_le(
  coeffs : Array[Double],
  rhs : Double,
) -> @src.Cobra3LinConstraint {
  { coeffs, rhs, sense: @src.Cobra3Sense::LessEqual, }
}

///|
fn cobra3_con_eq(
  coeffs : Array[Double],
  rhs : Double,
) -> @src.Cobra3LinConstraint {
  { coeffs, rhs, sense: @src.Cobra3Sense::Equal, }
}

///|
fn cobra3_con_ge(
  coeffs : Array[Double],
  rhs : Double,
) -> @src.Cobra3LinConstraint {
  { coeffs, rhs, sense: @src.Cobra3Sense::GreaterEqual, }
}

///|
test "cobra3 lp classic max 3x+5y" {
  let cons = [
    cobra3_con_le([1.0, 0.0], 4.0),
    cobra3_con_le([0.0, 2.0], 12.0),
    cobra3_con_le([3.0, 2.0], 18.0),
  ]
  let sol = @src.cobra3_solve_lp(cons, [3.0, 5.0], [0.0, 0.0], [1.0e18, 1.0e18])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 36.0).abs() < 1.0e-6)
  assert_true((sol.x[0] - 2.0).abs() < 1.0e-6)
  assert_true((sol.x[1] - 6.0).abs() < 1.0e-6)
}

///|
test "cobra3 lp dual of x<=2" {
  let cons = [cobra3_con_le([1.0], 2.0)]
  let sol = @src.cobra3_solve_lp(cons, [1.0], [0.0], [1.0e18])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 2.0).abs() < 1.0e-9)
  assert_true((sol.duals[0] - 1.0).abs() < 1.0e-9)
}

///|
test "cobra3 lp equality dual" {
  // max x s.t. x = 3, 0<=x<=10 -> obj 3, dual 1
  let cons = [cobra3_con_eq([1.0], 3.0)]
  let sol = @src.cobra3_solve_lp(cons, [1.0], [0.0], [10.0])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 3.0).abs() < 1.0e-9)
  assert_true((sol.duals[0] - 1.0).abs() < 1.0e-9)
}

///|
test "cobra3 lp infeasible" {
  let cons = [cobra3_con_eq([1.0], 5.0)]
  let sol = @src.cobra3_solve_lp(cons, [1.0], [0.0], [2.0])
  assert_true(sol.status is @src.Cobra3Status::Infeasible)
}

///|
test "cobra3 lp unbounded" {
  let cons = [cobra3_con_ge([1.0], 0.0)]
  let sol = @src.cobra3_solve_lp(cons, [1.0], [0.0], [1.0e18])
  assert_true(sol.status is @src.Cobra3Status::Unbounded)
}

///|
test "cobra3 lp negative lower bound shift" {
  // max -x, x >= -4, x <= 1 -> x=-4, obj 4
  let cons : Array[@src.Cobra3LinConstraint] = []
  let sol = @src.cobra3_solve_lp(cons, [-1.0], [-4.0], [1.0])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.x[0] + 4.0).abs() < 1.0e-9)
  assert_true((sol.objective_value - 4.0).abs() < 1.0e-9)
}

///|
test "cobra3 milp integer knapsack" {
  // max 3a+2b, 2a+1b<=4, a,b binary
  let cons = [cobra3_con_le([2.0, 1.0], 4.0)]
  let sol = @src.cobra3_solve_milp(cons, [3.0, 2.0], [0.0, 0.0], [1.0, 1.0], [
    0, 1,
  ])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 5.0).abs() < 1.0e-6)
  assert_true((sol.x[0] - 1.0).abs() < 1.0e-6)
  assert_true((sol.x[1] - 1.0).abs() < 1.0e-6)
}

///|
test "cobra3 milp chooses fractional rounding" {
  // max 5a+4b, 6a+4b<=8, binary -> only b (obj 4) 优于 a (obj 5)? 5a 单独可行
  // 6a<=8 -> a=1 可行 obj 5; b=1 -> 4; a=b -> 1.33 不可行
  let cons = [cobra3_con_le([6.0, 4.0], 8.0)]
  let sol = @src.cobra3_solve_milp(cons, [5.0, 4.0], [0.0, 0.0], [1.0, 1.0], [
    0, 1,
  ])
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 5.0).abs() < 1.0e-6)
}

///|
test "cobra3 qp box equality" {
  // min (x-3)^2+(y-2)^2 s.t. x+y=4, 0<=x,y<=3
  // 最优 x=2.5,y=1.5
  let q = [[2.0, 0.0], [0.0, 2.0]]
  let sol = @src.cobra3_solve_qp(
    q,
    [-6.0, -4.0],
    [[1.0, 1.0]],
    [4.0],
    [0.0, 0.0],
    [3.0, 3.0],
  )
  assert_true(sol.converged)
  assert_true((sol.x[0] - 2.5).abs() < 1.0e-4)
  assert_true((sol.x[1] - 1.5).abs() < 1.0e-4)
}

///|
test "cobra3 qp respects box" {
  // min (x-10)^2 s.t. no equality, 0<=x<=4 -> x=4
  let sol = @src.cobra3_solve_qp([[2.0]], [-20.0], [], [], [0.0], [4.0])
  assert_true(sol.converged)
  assert_true((sol.x[0] - 4.0).abs() < 1.0e-6)
}
