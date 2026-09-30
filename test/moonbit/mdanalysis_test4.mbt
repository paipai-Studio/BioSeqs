///|
test "mdanalysis pca 2 atoms" {
  // 2 atoms, 3 frames with different fluctuations
  // Atom 0 moves along x, atom 1 moves along y
  let f1 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let f2 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let f3 = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
  let (evals, evects, var_frac) = @src.mdanalysis_pca([f1, f2, f3])
  assert_eq(evals.length(), 2)
  assert_eq(evects.length(), 2)
  assert_eq(var_frac.length(), 2)
  // Both atoms have same variance, so eigenvalues should be equal
  // var per atom = (1+1)/2 = 1.0 for each atom
  assert_true(evals[0] > 0.0)
  assert_true(evals[1] > 0.0)
  // Cumulative variance fraction should end at 1.0
  assert_true((var_frac[1] - 1.0).abs() < 1.0e-9)
}

///|
test "mdanalysis pca single frame" {
  let f1 = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
  let (evals, _, _) = @src.mdanalysis_pca([f1])
  // Single frame → no variance
  assert_eq(evals.length(), 0)
}

///|
test "mdanalysis pca dominant motion" {
  // 3 atoms, 5 frames
  // Atom 0 has large motion, atoms 1 and 2 are static
  let frames : Array[Array[Array[Double]]] = []
  for i = 0; i < 5; i = i + 1 {
    let x = (i - 2).to_double() * 2.0 // large motion for atom 0
    frames.push([[x, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
  }
  let (evals, _, var_frac) = @src.mdanalysis_pca(frames)
  assert_eq(evals.length(), 3)
  // First eigenvalue should be much larger than the rest
  assert_true(evals[0] > evals[1])
  assert_true(evals[0] > evals[2])
  // First PC should capture most variance
  assert_true(var_frac[0] > 0.99)
}

///|
test "mdanalysis gnm kirchhoff" {
  // 3 atoms in a line: 0-1-2, cutoff connects neighbors
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let k = @src.mdanalysis_gnm_kirchhoff(positions, 1.5)
  // 0-1 connected, 1-2 connected, 0-2 not connected (dist=2 > 1.5)
  // Γ = [[1, -1, 0], [-1, 2, -1], [0, -1, 1]]
  assert_eq(k[0][0], 1.0)
  assert_eq(k[0][1], -1.0)
  assert_eq(k[0][2], 0.0)
  assert_eq(k[1][0], -1.0)
  assert_eq(k[1][1], 2.0)
  assert_eq(k[1][2], -1.0)
  assert_eq(k[2][0], 0.0)
  assert_eq(k[2][1], -1.0)
  assert_eq(k[2][2], 1.0)
}

///|
test "mdanalysis gnm eigenvalues" {
  // 3 atoms in a line
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let (evals, evects) = @src.mdanalysis_gnm(positions, 1.5)
  assert_eq(evals.length(), 3)
  assert_eq(evects.length(), 3)
  // Smallest eigenvalue should be ~0 (rigid-body mode)
  assert_true(evals[0].abs() < 1.0e-6)
  // Other eigenvalues should be positive
  assert_true(evals[1] > 0.0)
  assert_true(evals[2] > 0.0)
  // For this simple chain, eigenvalues are 0, 1, 3
  assert_true((evals[1] - 1.0).abs() < 1.0e-6)
  assert_true((evals[2] - 3.0).abs() < 1.0e-6)
}

///|
test "mdanalysis gnm bfactors" {
  // 4 atoms in a line
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
  ]
  let (evals, evects) = @src.mdanalysis_gnm(positions, 1.5)
  let bf = @src.mdanalysis_gnm_bfactors(evals, evects)
  assert_eq(bf.length(), 4)
  // B-factors should be positive (more flexibility = higher B-factor)
  for i = 0; i < 4; i = i + 1 {
    assert_true(bf[i] > 0.0)
  }
  // Terminal atoms should have higher B-factors than middle atoms
  assert_true(bf[0] > bf[1])
  assert_true(bf[3] > bf[2])
}

///|
test "mdanalysis helix analysis ideal alpha" {
  // Generate an ideal α-helix: 3.6 residues/turn, rise=1.5Å, radius=2.3Å
  // φ = -57°, ψ = -47° → twist ≈ 100° per residue
  let n_res = 10
  let positions : Array[Array[Double]] = []
  let rise = 1.5
  let radius = 2.3
  let twist = 100.0 * @math.pi / 180.0 // 100° per residue
  for i = 0; i < n_res; i = i + 1 {
    let angle = i.to_double() * twist
    let x = radius * @math.cos(angle)
    let y = radius * @math.sin(angle)
    let z = i.to_double() * rise
    positions.push([x, y, z])
  }
  let result = @src.mdanalysis_helix_analysis(positions)
  // N-3 = 7 windows
  assert_eq(result.local_twists.length(), 7)
  assert_eq(result.local_nres_per_turn.length(), 7)
  assert_eq(result.local_axes.length(), 7)
  assert_eq(result.local_heights.length(), 7)
  assert_eq(result.local_radii.length(), 7)
  // N-6 = 4 bends
  assert_eq(result.local_bends.length(), 4)
  // Ideal helix: twist ≈ 100°, nres/turn ≈ 3.6, rise ≈ 1.5, radius ≈ 2.3
  for i = 0; i < 7; i = i + 1 {
    assert_true((result.local_twists[i] - 100.0).abs() < 5.0)
    assert_true((result.local_nres_per_turn[i] - 3.6).abs() < 0.5)
    assert_true((result.local_heights[i] - 1.5).abs() < 0.3)
    assert_true((result.local_radii[i] - 2.3).abs() < 0.5)
  }
  // Ideal helix: bends should be near 0
  for i = 0; i < 4; i = i + 1 {
    assert_true(result.local_bends[i] < 10.0)
  }
}

///|
test "mdanalysis helix analysis too short" {
  // 3 atoms — not enough for a 4-atom window
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let result = @src.mdanalysis_helix_analysis(positions)
  assert_eq(result.local_twists.length(), 0)
}

///|
test "mdanalysis helix analysis linear" {
  // Linear chain (not a helix): all atoms on z-axis
  let positions : Array[Array[Double]] = []
  for i = 0; i < 8; i = i + 1 {
    positions.push([0.0, 0.0, i.to_double()])
  }
  let result = @src.mdanalysis_helix_analysis(positions)
  // N-3 = 5 windows
  assert_eq(result.local_twists.length(), 5)
  // For a linear chain, bisectors are all parallel → twist ≈ 0
  // This is a degenerate case; just verify it doesn't crash
  // and produces finite values
  for i = 0; i < 5; i = i + 1 {
    assert_true(result.local_twists[i] >= 0.0)
    assert_true(result.local_twists[i] <= 180.0)
  }
}
