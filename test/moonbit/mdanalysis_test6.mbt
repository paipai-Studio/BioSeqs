///|
test "mdanalysis density basic" {
  // 2 frames, 2 atoms each, in a 2×2×2 grid
  let f1 = [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]
  let f2 = [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]
  let result = @src.mdanalysis_density([f1, f2], [0.0, 0.0, 0.0], [2.0, 1.0, 1.0], (2, 1, 1))
  assert_eq(result.bins.0, 2)
  assert_eq(result.bins.1, 1)
  assert_eq(result.bins.2, 1)
  // Voxel volume = (2/2) * (1/1) * (1/1) = 1.0
  assert_true((result.voxel_volume - 1.0).abs() < 1.0e-9)
  // Grid has 2 bins: atom 0 in bin 0, atom 1 in bin 1
  // Each bin gets 1 atom per frame, 2 frames → 2 counts
  // density = 2 / (2 frames * 1.0 voxel) = 1.0
  assert_eq(result.grid.length(), 2)
  assert_true((result.grid[0] - 1.0).abs() < 1.0e-9)
  assert_true((result.grid[1] - 1.0).abs() < 1.0e-9)
}

///|
test "mdanalysis density out of bounds" {
  // Atom outside grid → not counted
  let f1 = [[5.0, 5.0, 5.0]]
  let result = @src.mdanalysis_density([f1], [0.0, 0.0, 0.0], [2.0, 2.0, 2.0], (2, 2, 2))
  // All bins should be 0
  for i = 0; i < 8; i = i + 1 {
    assert_true(result.grid[i] < 1.0e-15)
  }
}

///|
test "mdanalysis density empty" {
  let result = @src.mdanalysis_density([], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], (1, 1, 1))
  assert_eq(result.grid.length(), 0)
}

///|
test "mdanalysis density voxel volume" {
  // 3×3×3 grid in a 3×3×3 box → voxel = 1.0
  let result = @src.mdanalysis_density([[[0.5, 0.5, 0.5]]], [0.0, 0.0, 0.0], [3.0, 3.0, 3.0], (3, 3, 3))
  assert_true((result.voxel_volume - 1.0).abs() < 1.0e-9)
  // 1 atom in bin (0,0,0), 1 frame → density = 1/(1*1.0) = 1.0
  assert_true((result.grid[0] - 1.0).abs() < 1.0e-9)
}

///|
test "mdanalysis diffusionmap basic" {
  // 3 frames with known distances
  // Frame 0-1: dist=1, Frame 1-2: dist=1, Frame 0-2: dist=2
  let dist = [
    [0.0, 1.0, 2.0],
    [1.0, 0.0, 1.0],
    [2.0, 1.0, 0.0],
  ]
  let result = @src.mdanalysis_diffusionmap(dist, 1.0)
  assert_eq(result.eigenvalues.length(), 3)
  assert_eq(result.eigenvectors.length(), 3)
  // Largest eigenvalue should be ~1 (stochastic matrix)
  assert_true(result.eigenvalues[0] > result.eigenvalues[1])
  assert_true(result.eigenvalues[1] > result.eigenvalues[2])
  // Largest eigenvalue of a diffusion matrix is ~1
  assert_true((result.eigenvalues[0] - 1.0).abs() < 0.1)
}

///|
test "mdanalysis diffusionmap identical frames" {
  // All frames identical → all distances 0 → kernel all 1
  let dist = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  let result = @src.mdanalysis_diffusionmap(dist, 1.0)
  assert_eq(result.eigenvalues.length(), 3)
  // All eigenvalues should be ~1 (uniform kernel → identity-like)
  // Actually for uniform kernel, diffusion matrix is rank-1
  // Largest eigenvalue = 1, rest = 0
  assert_true((result.eigenvalues[0] - 1.0).abs() < 0.1)
}

///|
test "mdanalysis diffusionmap transform" {
  let dist = [
    [0.0, 1.0, 2.0],
    [1.0, 0.0, 1.0],
    [2.0, 1.0, 0.0],
  ]
  let result = @src.mdanalysis_diffusionmap(dist, 1.0)
  let embedding = @src.mdanalysis_diffusionmap_transform(result, 2, 1.0)
  // 3 frames × 2 eigenvectors
  assert_eq(embedding.length(), 3)
  assert_eq(embedding[0].length(), 2)
}

///|
test "mdanalysis diffusionmap too small" {
  let result = @src.mdanalysis_diffusionmap([[0.0]], 1.0)
  assert_eq(result.eigenvalues.length(), 0)
}

///|
test "mdanalysis dielectric basic" {
  // 2 atoms with charges +1 and -1
  // Frame 1: positions (0,0,0) and (1,0,0) → M = (0+(-1), 0, 0) = (-1, 0, 0)
  // Frame 2: positions (0,0,0) and (-1,0,0) → M = (0-(-1), 0, 0) = (1, 0, 0)
  let charges = [1.0, -1.0]
  let f1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let f2 = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
  let result = @src.mdanalysis_dielectric(charges, [f1, f2], 300.0, 1000.0)
  // Mean dipole: ((-1+1)/2, 0, 0) = (0, 0, 0)
  assert_true(result.mean_dipole[0].abs() < 1.0e-9)
  assert_true(result.mean_dipole[1].abs() < 1.0e-9)
  assert_true(result.mean_dipole[2].abs() < 1.0e-9)
  // Fluctuation x: <Mx^2> - <Mx>^2 = (1+1)/2 - 0 = 1.0
  assert_true((result.fluctuation[0] - 1.0).abs() < 1.0e-9)
  // y and z fluctuations should be 0
  assert_true(result.fluctuation[1].abs() < 1.0e-9)
  assert_true(result.fluctuation[2].abs() < 1.0e-9)
  // eps_x should be > 1 (fluctuation > 0)
  assert_true(result.eps[0] > 1.0)
  // eps_y and eps_z should be ~1 (no fluctuation)
  assert_true((result.eps[1] - 1.0).abs() < 1.0e-6)
  assert_true((result.eps[2] - 1.0).abs() < 1.0e-6)
}

///|
test "mdanalysis dielectric no fluctuation" {
  // Same dipole every frame → no fluctuation → eps = 1
  let charges = [1.0, -1.0]
  let f1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let f2 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let result = @src.mdanalysis_dielectric(charges, [f1, f2], 300.0, 1000.0)
  // All eps should be ~1
  for i = 0; i < 3; i = i + 1 {
    assert_true((result.eps[i] - 1.0).abs() < 1.0e-6)
  }
  assert_true((result.eps_mean - 1.0).abs() < 1.0e-6)
}

///|
test "mdanalysis dielectric empty" {
  let result = @src.mdanalysis_dielectric([], [], 300.0, 1000.0)
  assert_true((result.eps_mean - 1.0).abs() < 1.0e-9)
}
