///|
test "mdanalysis tors alpha" {
  // Simple 4-atom chain along z-axis
  // O3' at (0,0,0), P at (1,0,0), O5' at (1,1,0), C5' at (1,1,1)
  let a = @src.mdanalysis_tors_alpha([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0])
  // Should be a valid angle in [0, 360)
  assert_true(a >= 0.0)
  assert_true(a < 360.0)
}

///|
test "mdanalysis tors all 7" {
  // Build a simple nucleic acid residue with atoms
  let prev : Map[String, Array[Double]] = Map([])
  prev["O3'"] = [0.0, 0.0, 0.0]
  let atoms : Map[String, Array[Double]] = Map([])
  atoms["P"] = [1.5, 0.0, 0.0]
  atoms["O5'"] = [2.0, 1.0, 0.0]
  atoms["C5'"] = [2.0, 1.0, 1.5]
  atoms["C4'"] = [3.0, 1.0, 1.5]
  atoms["C3'"] = [3.0, 2.0, 1.5]
  atoms["O3'"] = [3.0, 2.0, 3.0]
  atoms["O4'"] = [2.5, 0.5, 1.5]
  atoms["C1'"] = [2.5, 0.5, 2.5]
  atoms["N9"] = [3.5, 0.5, 2.5]
  atoms["C4"] = [3.5, 1.5, 2.5]
  let next : Map[String, Array[Double]] = Map([])
  next["P"] = [4.0, 2.0, 3.0]
  next["O5'"] = [4.5, 3.0, 3.0]
  let (alpha, beta, gamma, delta, eps, zeta, chi) = @src.mdanalysis_tors_all(prev, atoms, next)
  // All should be in [0, 360)
  for a in [alpha, beta, gamma, delta, eps, zeta, chi] {
    assert_true(a >= 0.0)
    assert_true(a < 360.0)
  }
}

///|
test "mdanalysis phase cp planar" {
  // 5 atoms in a plane — phase should be well-defined
  // Regular pentagon in xy-plane
  let atoms : Array[Array[Double]] = []
  for i = 0; i < 5; i = i + 1 {
    let angle = 2.0 * @math.pi * i.to_double() / 5.0
    atoms.push([@math.cos(angle), @math.sin(angle), 0.0])
  }
  let phase = @src.mdanalysis_phase_cp(atoms)
  assert_true(phase >= 0.0)
  assert_true(phase < 360.0)
}

///|
test "mdanalysis phase cp c3_endo" {
  // C3'-endo pucker (North, typical in A-DNA): phase ~18°
  // Construct a ribose ring with C3' puckered up
  // Base positions roughly in a plane, C3' displaced
  let o4p = [0.0, 0.0, 0.0]
  let c1p = [1.5, 0.0, 0.0]
  let c2p = [2.0, 1.2, 0.0]
  let c3p = [1.0, 2.0, 0.5]  // C3' puckered up (endo)
  let c4p = [0.0, 1.2, 0.0]
  let phase = @src.mdanalysis_phase_cp([o4p, c1p, c2p, c3p, c4p])
  assert_true(phase >= 0.0)
  assert_true(phase < 360.0)
}

///|
test "mdanalysis phase as" {
  // Same ribose ring as phase_cp test
  let o4p = [0.0, 0.0, 0.0]
  let c1p = [1.5, 0.0, 0.0]
  let c2p = [2.0, 1.2, 0.0]
  let c3p = [1.0, 2.0, 0.5]
  let c4p = [0.0, 1.2, 0.0]
  let phase = @src.mdanalysis_phase_as([o4p, c1p, c2p, c3p, c4p])
  assert_true(phase >= 0.0)
  assert_true(phase < 360.0)
}

///|

///|
test "mdanalysis leaflet two layers" {
  // 6 lipid heads: 3 in upper leaflet (z=20), 3 in lower (z=-20)
  let heads = [
    [0.0, 0.0, 20.0],
    [5.0, 0.0, 20.0],
    [0.0, 5.0, 20.0],
    [0.0, 0.0, -20.0],
    [5.0, 0.0, -20.0],
    [0.0, 5.0, -20.0],
  ]
  let labels = @src.mdanalysis_leaflet_finder(heads, 15.0)
  assert_eq(labels.length(), 6)
  // Upper layer atoms (0,1,2) should have same label
  assert_eq(labels[0], labels[1])
  assert_eq(labels[1], labels[2])
  // Lower layer atoms (3,4,5) should have same label
  assert_eq(labels[3], labels[4])
  assert_eq(labels[4], labels[5])
  // Upper and lower should be different
  assert_true(labels[0] != labels[3])
}

///|
test "mdanalysis leaflet single layer" {
  // 3 atoms all close together → 1 leaflet
  let heads = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]]
  let labels = @src.mdanalysis_leaflet_finder(heads, 5.0)
  assert_eq(labels[0], labels[1])
  assert_eq(labels[1], labels[2])
}

///|
test "mdanalysis leaflet empty" {
  let labels = @src.mdanalysis_leaflet_finder([], 15.0)
  assert_eq(labels.length(), 0)
}

///|
test "mdanalysis bat root" {
  // 3 root atoms
  let p0 = [0.0, 0.0, 0.0]
  let p1 = [1.0, 0.0, 0.0]
  let p2 = [1.0, 1.0, 0.0]
  let (phi, theta, omega, r01, r12, a012) = @src.mdanalysis_bat_root(p0, p1, p2)
  // r01 = dist(p0,p1) = 1.0
  assert_true((r01 - 1.0).abs() < 1.0e-9)
  // r12 = dist(p1,p2) = 1.0
  assert_true((r12 - 1.0).abs() < 1.0e-9)
  // a012 = angle p0-p1-p2 = 90° = π/2
  assert_true((a012 - @math.pi / 2.0).abs() < 1.0e-9)
  // phi = atan2(0, 1) = 0 (v01 points along x)
  assert_true(phi.abs() < 1.0e-9)
  // theta = acos(0) = π/2 (v01 has z=0)
  assert_true((theta - @math.pi / 2.0).abs() < 1.0e-9)
}

///|
test "mdanalysis bat internal" {
  // 2 sets of 4 consecutive atoms
  let ag1 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let ag2 = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let ag3 = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
  let ag4 = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
  let (bonds, angles, torsions) = @src.mdanalysis_bat_internal(ag1, ag2, ag3, ag4)
  assert_eq(bonds.length(), 2)
  assert_eq(angles.length(), 2)
  assert_eq(torsions.length(), 2)
  // Bond = 1.0
  assert_true((bonds[0] - 1.0).abs() < 1.0e-9)
  // Angle = 90° = π/2
  assert_true((angles[0] - @math.pi / 2.0).abs() < 1.0e-9)
}
