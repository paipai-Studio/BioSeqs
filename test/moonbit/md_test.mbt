///|
/// Tests for the minimal MD engine (md.mbt).

/// 数值梯度（中心差分）用于验证解析力 = -∇E。
fn num_gradient(
  positions : Array[Array[Double]],
  bonds : Array[(Int, Int, Double, Double)],
  i : Int,
  d : Int,
  h : Double,
) -> Double {
  let pos_plus = positions.copy()
  let pos_minus = positions.copy()
  pos_plus[i] = pos_plus[i].copy()
  pos_minus[i] = pos_minus[i].copy()
  pos_plus[i][d] = pos_plus[i][d] + h
  pos_minus[i][d] = pos_minus[i][d] - h
  (@src.md_bond_energy(pos_plus, bonds) - @src.md_bond_energy(pos_minus, bonds)) / (2.0 * h)
}

test "md_bond_energy" {
  let pos = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]
  let bonds = [(0, 1, 1.0, 10.0)]
  // 0.5 * 10 * (1.1-1)^2 = 0.05
  assert_true(assert_close(@src.md_bond_energy(pos, bonds), 0.05, 1.0e-9))
  // 平衡位置 → 0
  assert_true(assert_close(@src.md_bond_energy([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], bonds), 0.0, 1.0e-9))
}

test "md_bond_forces_match_gradient" {
  let pos = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]
  let bonds = [(0, 1, 1.0, 10.0)]
  let f = @src.md_bond_forces(pos, bonds)
  // 解析力 = -梯度
  assert_true(assert_close(f[0][0], -num_gradient(pos, bonds, 0, 0, 1.0e-6), 1.0e-5))
  assert_true(assert_close(f[1][0], -num_gradient(pos, bonds, 1, 0, 1.0e-6), 1.0e-5))
  // 作用力=反作用力
  assert_true(assert_close(f[0][0], -f[1][0], 1.0e-12))
}

test "md_lj_energy_minimum" {
  // LJ 最小值在 r = 2^(1/6) σ，E = -ε
  let r_min = @math.pow(2.0, 1.0 / 6.0)
  let pos = [[0.0, 0.0, 0.0], [r_min, 0.0, 0.0]]
  let pairs = [(0, 1, 1.0, 1.0)]
  assert_true(assert_close(@src.md_lj_energy(pos, pairs), -1.0, 1.0e-9))
}

test "md_coulomb_energy" {
  // ke * q1 q2 / r = 332.0637 / 2
  let pos = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let pairs = [(0, 1, 1.0, 1.0)]
  assert_true(assert_close(@src.md_coulomb_energy(pos, pairs), 332.0637 / 2.0, 1.0e-6))
}

test "md_energy_conservation" {
  // 双原子谐振子：NVE 下总能量应近似守恒
  let pos = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]
  let vel = [[0.3, 0.0, 0.0], [-0.3, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let bonds = [(0, 1, 1.0, 10.0)]
  let energies = @src.md_simulate_harmonic(pos, vel, masses, bonds, 0.0005, 1000)
  let e0 = energies[0]
  let mut max_drift = 0.0
  for e in energies {
    let drift = (e - e0).abs()
    if drift > max_drift {
      max_drift = drift
    }
  }
  // Verlet 能量漂移应远小于能量本身
  assert_true(max_drift < (0.01 * e0.abs() + 1.0e-6))
}
