///|
test "mdanalysis contacts basic" {
  // Two groups of atoms, some within cutoff
  let a = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
  let b = [[1.0, 0.0, 0.0], [6.0, 0.0, 0.0], [20.0, 0.0, 0.0]]
  let pairs = @src.mdanalysis_contacts(a, b, 3.0)
  // a[0]-b[0]: dist=1 ≤ 3 ✓
  // a[1]-b[1]: dist=1 ≤ 3 ✓
  // a[2]-b[2]: dist=10 > 3 ✗
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0], (0, 0))
  assert_eq(pairs[1], (1, 1))
}

///|
test "mdanalysis contacts empty" {
  let a = [[0.0, 0.0, 0.0]]
  let b = [[100.0, 0.0, 0.0]]
  let pairs = @src.mdanalysis_contacts(a, b, 5.0)
  assert_eq(pairs.length(), 0)
}

///|
test "mdanalysis contact fraction hard cut" {
  // Reference: 2 contacts at dist 1.0 and 1.0
  let ref_a = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
  let ref_b = [[1.0, 0.0, 0.0], [6.0, 0.0, 0.0]]
  // Frame 1: both contacts maintained
  let f1_a = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
  let f1_b = [[1.0, 0.0, 0.0], [6.0, 0.0, 0.0]]
  // Frame 2: first contact broken (moved far)
  let f2_a = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
  let f2_b = [[50.0, 0.0, 0.0], [6.0, 0.0, 0.0]]
  let frames_a = [f1_a, f2_a]
  let frames_b = [f1_b, f2_b]
  let q = @src.mdanalysis_contact_fraction(ref_a, ref_b, frames_a, frames_b, 3.0, true)
  assert_eq(q.length(), 2)
  // Frame 1: both contacts present → Q = 1.0
  assert_true((q[0] - 1.0).abs() < 1.0e-9)
  // Frame 2: only second contact → Q = 0.5
  assert_true((q[1] - 0.5).abs() < 1.0e-9)
}

///|
test "mdanalysis ramachandran region" {
  // α-helix: φ=-57°, ψ=-47° → core
  assert_eq(@src.mdanalysis_ramachandran_region(-57.0, -47.0), "core")
  // β-sheet: φ=-120°, ψ=120° → core
  assert_eq(@src.mdanalysis_ramachandran_region(-120.0, 120.0), "core")
  // Left-handed helix: φ=60°, ψ=45° → core
  assert_eq(@src.mdanalysis_ramachandran_region(60.0, 45.0), "core")
  // Disallowed: φ=0°, ψ=180° → outlier
  assert_eq(@src.mdanalysis_ramachandran_region(0.0, 180.0), "outlier")
  // Borderline allowed: φ=-80°, ψ=20° → allowed
  assert_eq(@src.mdanalysis_ramachandran_region(-80.0, 20.0), "allowed")
}

///|
test "mdanalysis ramachandran multi frame" {
  // Build a simple 3-residue backbone
  // Residue i: N, CA, C positions (roughly linear for simplicity)
  let n1 = [0.0, 0.0, 0.0]
  let ca1 = [1.5, 0.0, 0.0]
  let c1 = [2.5, 1.0, 0.0]
  let n2 = [3.5, 0.5, 0.0]
  let ca2 = [4.5, 0.5, 0.0]
  let c2 = [5.5, 1.5, 0.0]
  let n3 = [6.5, 1.0, 0.0]
  let ca3 = [7.5, 1.0, 0.0]
  let c3 = [8.5, 2.0, 0.0]
  let frame1 = [(n1, ca1, c1), (n2, ca2, c2), (n3, ca3, c3)]
  // Slightly perturbed frame 2
  let n2b = [3.5, 0.6, 0.1]
  let ca2b = [4.5, 0.4, 0.1]
  let c2b = [5.5, 1.6, 0.1]
  let frame2 = [(n1, ca1, c1), (n2b, ca2b, c2b), (n3, ca3, c3)]
  let (angles, regions) = @src.mdanalysis_ramachandran([frame1, frame2])
  assert_eq(angles.length(), 2)
  assert_eq(regions.length(), 2)
  // Each frame has 3 residues
  assert_eq(angles[0].length(), 3)
  assert_eq(regions[0].length(), 3)
  // First residue has no phi (None), last has no psi (None)
  let (p0, _) = angles[0][0]
  let (_, s2) = angles[0][2]
  assert_true(p0.is_empty())
  assert_true(s2.is_empty())
  // Middle residue has both phi and psi
  let (p1, _) = angles[0][1]
  let (_, s1) = angles[0][1]
  assert_true(p1.is_some())
  assert_true(s1.is_some())
  // Middle residue region should be Some
  assert_true(regions[0][1].is_some())
  // First/last residue region is None (missing phi or psi)
  assert_true(regions[0][0].is_empty())
  assert_true(regions[0][2].is_empty())
}

///|
test "mdanalysis linear density 1d" {
  // 10 atoms evenly spaced from 0 to 10 along x
  let coords : Array[Array[Double]] = []
  for i = 0; i < 10; i = i + 1 {
    coords.push([i.to_double(), 0.0, 0.0])
  }
  // 5 bins over [0, 10) → each bin width = 2
  let hist = @src.mdanalysis_linear_density_1d(coords, 0, 5, 0.0, 10.0, None)
  assert_eq(hist.length(), 5)
  // bin 0: x=0,1 → 2 atoms
  // bin 1: x=2,3 → 2 atoms
  // bin 2: x=4,5 → 2 atoms
  // bin 3: x=6,7 → 2 atoms
  // bin 4: x=8,9 → 2 atoms
  for i = 0; i < 5; i = i + 1 {
    assert_eq(hist[i], 2.0)
  }
}

///|
test "mdanalysis linear density weighted" {
  let coords = [[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let weights = [2.0, 3.0]
  let hist = @src.mdanalysis_linear_density_1d(coords, 0, 2, 0.0, 2.0, Some(weights))
  // bin 0: x=0.5, weight=2.0
  // bin 1: x=1.5, weight=3.0
  assert_eq(hist[0], 2.0)
  assert_eq(hist[1], 3.0)
}

///|
test "mdanalysis linear density 3d" {
  let coords = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let box = [10.0, 10.0, 10.0]
  let (hx, hy, hz) = @src.mdanalysis_linear_density(coords, 5, box, None)
  assert_eq(hx.length(), 5)
  assert_eq(hy.length(), 5)
  assert_eq(hz.length(), 5)
  // x=1 → bin 0, x=4 → bin 2
  assert_eq(hx[0], 1.0)
  assert_eq(hx[2], 1.0)
  // y=2 → bin 1, y=5 → bin 2
  assert_eq(hy[1], 1.0)
  assert_eq(hy[2], 1.0)
  // z=3 → bin 1, z=6 → bin 3
  assert_eq(hz[1], 1.0)
  assert_eq(hz[3], 1.0)
}

///|
test "mdanalysis bond autocorrelation" {
  // Straight line along x: bond vectors all (1,0,0)
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [4.0, 0.0, 0.0],
  ]
  let (ac, lb) = @src.mdanalysis_bond_autocorrelation(positions)
  // 4 bonds, lb = 1.0
  assert_true((lb - 1.0).abs() < 1.0e-9)
  assert_eq(ac.length(), 4)
  // All bonds parallel → C(n) = 1.0 for all n
  for i = 0; i < 4; i = i + 1 {
    assert_true((ac[i] - 1.0).abs() < 1.0e-9)
  }
}

///|
test "mdanalysis bond autocorrelation bent" {
  // L-shape: first bond along x, second along y → 90° bend
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
  let (ac, lb) = @src.mdanalysis_bond_autocorrelation(positions)
  assert_eq(ac.length(), 2)
  // C(0) = 1.0 (self-correlation)
  assert_true((ac[0] - 1.0).abs() < 1.0e-9)
  // C(1) = cos(90°) = 0.0
  assert_true(ac[1].abs() < 1.0e-9)
  assert_true((lb - 1.0).abs() < 1.0e-9)
}

///|
test "mdanalysis persistence length straight" {
  // Straight line: lp should be very large (no decay)
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [4.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
  ]
  let (lp, lb, ac, fit) = @src.mdanalysis_persistence_length(positions)
  assert_eq(ac.length(), 5)
  assert_true((lb - 1.0).abs() < 1.0e-9)
  // For a straight line, C(n) = 1 for all n, ln(C) = 0, slope ≈ 0
  // lp = -lb/slope → very large or 0 (slope ≥ 0 edge case)
  // Just verify lp is non-negative and fit has correct length
  assert_true(lp >= 0.0)
  assert_eq(fit.length(), 5)
}

///|
test "mdanalysis persistence length decay" {
  // Gradually rotating chain: each bond rotates 30° from previous
  // This gives exponential decay in C(n)
  let positions : Array[Array[Double]] = [[0.0, 0.0, 0.0]]
  let mut angle = 0.0
  let bond_len = 1.0
  let mut x = 0.0
  let mut y = 0.0
  for i = 0; i < 10; i = i + 1 {
    x = x + bond_len * @math.cos(angle)
    y = y + bond_len * @math.sin(angle)
    positions.push([x, y, 0.0])
    angle = angle + 0.5 // 0.5 rad per bond
  }
  let (lp, lb, ac, fit) = @src.mdanalysis_persistence_length(positions)
  assert_eq(ac.length(), 10)
  assert_true((lb - 1.0).abs() < 1.0e-9)
  // C(0) = 1.0
  assert_true((ac[0] - 1.0).abs() < 1.0e-9)
  // C(n) should decrease monotonically
  assert_true(ac[1] < ac[0])
  // lp should be positive and finite
  assert_true(lp > 0.0)
  assert_true(lp < 100.0)
  // Fit curve should match at n=0
  assert_true((fit[0] - 1.0).abs() < 1.0e-6)
}
