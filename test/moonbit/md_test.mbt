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


// ============================================================
// angle + dihedral
// ============================================================

// 数值梯度 for angle energy
fn angle_num_gradient(
  positions : Array[Array[Double]],
  angles : Array[(Int, Int, Int, Double, Double)],
  i : Int,
  d : Int,
  h : Double,
) -> Double {
  let plus = positions.copy()
  let minus = positions.copy()
  plus[i] = plus[i].copy()
  minus[i] = minus[i].copy()
  plus[i][d] = plus[i][d] + h
  minus[i][d] = minus[i][d] - h
  (@src.md_angle_energy(plus, angles) - @src.md_angle_energy(minus, angles)) / (2.0 * h)
}

test "md_angle_energy" {
  // i=(1,0,0) j=(0,0,0) k=(0,1,0)：∠i-j-k = 90° = π/2
  let pos = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  // θ0 = π/2 → 平衡，能量 0
  assert_true(assert_close(@src.md_angle_energy(pos, [(0, 1, 2, 1.5707963267948966, 1.0)]), 0.0, 1.0e-9))
  // θ0 = π/3（60°）→ θ 偏离 π/6，E = 0.5 * (π/6)^2
  let e = @src.md_angle_energy(pos, [(0, 1, 2, 1.0471975511965976, 1.0)])
  assert_true(assert_close(e, 0.5 * 0.5235987755982988 * 0.5235987755982988, 1.0e-9))
}

test "md_angle_forces_match_gradient" {
  let pos = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let angles = [(0, 1, 2, 1.0471975511965976, 2.5)]
  let f = @src.md_angle_forces(pos, angles)
  assert_true(assert_close(f[0][0], -angle_num_gradient(pos, angles, 0, 0, 1.0e-6), 1.0e-5))
  assert_true(assert_close(f[2][1], -angle_num_gradient(pos, angles, 2, 1, 1.0e-6), 1.0e-5))
  // 内力守恒：三个原子合力 ≈ 0
  let fsum = [
    f[0][0] + f[1][0] + f[2][0],
    f[0][1] + f[1][1] + f[2][1],
    f[0][2] + f[1][2] + f[2][2],
  ]
  assert_true(assert_close(fsum[0], 0.0, 1.0e-10))
  assert_true(assert_close(fsum[1], 0.0, 1.0e-10))
}

test "md_dihedral_energy_trans_cis" {
  // trans：i=(1,0,0) j=(0,0,0) k=(0,1,0) l=(1,1,0) → φ=π
  let trans = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
  let cis = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]
  // n=3, δ=0, k=1：trans φ=π → E=1+cos(3π)=0；cis φ=0 → E=1+cos(0)=2
  let d_list = [(0, 1, 2, 3, 1.0, 3, 0.0)]
  assert_true(assert_close(@src.md_dihedral_energy(trans, d_list), 0.0, 1.0e-9))
  assert_true(assert_close(@src.md_dihedral_energy(cis, d_list), 2.0, 1.0e-9))
}

test "md_dihedral_forces_conserved" {
  // l=(0,1,1) → φ=π/2（非极值点，dE/dφ ≠ 0，力非零）
  let mid = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
  let d_list = [(0, 1, 2, 3, 1.0, 3, 0.0)]
  let f = @src.md_dihedral_forces(mid, d_list)
  // 四个原子合力 ≈ 0（内力守恒；数值梯度精度 ~1e-6）
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  for v in f {
    sx = sx + v[0]
    sy = sy + v[1]
    sz = sz + v[2]
  }
  assert_true(assert_close(sx, 0.0, 1.0e-6))
  assert_true(assert_close(sy, 0.0, 1.0e-6))
  assert_true(assert_close(sz, 0.0, 1.0e-6))
  // 力非零
  let mut mag = 0.0
  for v in f {
    mag = mag + v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
  }
  assert_true(mag > 1.0e-3)
}
