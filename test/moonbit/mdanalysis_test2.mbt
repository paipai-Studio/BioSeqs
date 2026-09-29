///|
/// Tests for mdanalysis batch 2: hydrogen bonds and phi/psi dihedrals.

test "mdanalysis_hydrogen_bonds" {
  // D=atom0(0,0,0) H=atom1(1,0,0)
  // A=atom2(2,0,0): D-A=2.0≤3.0, D-H-A=180°≥150° → 通过
  // A=atom3(4,0,0): D-A=4.0>3.0 → 不通过
  // A=atom4(1,1,0): D-A=√2≤3.0, 但 D-H-A=90°<150° → 不通过
  let coords = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [4.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
  ]
  let pairs = [(0, 1)]
  let acceptors = [2, 3, 4]
  let hbonds = @src.mdanalysis_hydrogen_bonds(coords, pairs, acceptors, 3.0, 150.0)
  assert_eq(hbonds.length(), 1)
  let (d, h, a, dist, ang) = hbonds[0]
  assert_eq(d, 0)
  assert_eq(h, 1)
  assert_eq(a, 2)
  assert_true(assert_close(dist, 2.0, 1.0e-9))
  assert_true(assert_close(ang, 180.0, 1.0e-6))
}

test "mdanalysis_phi_psi_planar" {
  // 3 残基全在 z=0 平面 → 所有二面角为 0° 或 ±180°（共面）
  let bb = [
    ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
    ([0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [-1.0, 2.0, 0.0]),
    ([-1.0, 1.0, 0.0], [-2.0, 1.0, 0.0], [-2.0, 0.0, 0.0]),
  ]
  let result = @src.mdanalysis_phi_psi(bb)
  assert_eq(result.length(), 3)
  // 首残基 phi = None
  assert_true(result[0].0 is None)
  // 尾残基 psi = None
  assert_true(result[2].1 is None)
  // 中间残基 phi/psi 都是 Some 且有限
  let phi1 = result[1].0
  assert_true(phi1 is Some(_))
  match phi1 {
    Some(phi) => assert_true(phi.abs() <= 180.0)
    None => ()
  }
  let psi0 = result[0].1
  assert_true(psi0 is Some(_))
  match psi0 {
    Some(psi) => assert_true(psi.abs() <= 180.0)
    None => ()
  }
  // phi2 和 psi1 也应该是 Some
  assert_true(result[2].0 is Some(_))
  assert_true(result[1].1 is Some(_))
}
