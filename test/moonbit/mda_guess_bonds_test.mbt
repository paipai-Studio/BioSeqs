///|
/// Tests for covalent bond guessing (mda_guess_bonds.mbt).
test "mda_guess_bonds_water" {
  // O 在原点，两个 H 在 0.96 A 处成 90° 角：
  // O-H 截止 0.66+0.31+0.45=1.42 → 成键；H-H 距离 1.358 > 0.31+0.31+0.45=1.07 → 不成键
  let coords = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [0.0, 0.96, 0.0]]
  let elements = ["O", "H", "H"]
  let bonds = @src.mda_guess_bonds(coords, elements)
  assert_eq(bonds.length(), 2)
  let (a0, b0) = bonds[0]
  let (a1, b1) = bonds[1]
  assert_eq(a0, 0)
  assert_eq(b0, 1)
  assert_eq(a1, 0)
  assert_eq(b1, 2)
}

///|
test "mda_guess_bonds_pbc" {
  // C-C 跨 box 边界：0.5 与 9.8 相距 9.3，最小镜像后 0.7 <= 0.76+0.76+0.45=1.97 → 成键
  let coords = [[0.5, 0.0, 0.0], [9.8, 0.0, 0.0]]
  let elements = ["C", "C"]
  let with_pbc = @src.mda_guess_bonds(
    coords,
    elements,
    box=Some([10.0, 10.0, 10.0]),
  )
  assert_eq(with_pbc.length(), 1)
  // 不开 PBC 时 9.3 远超截止 → 无键
  let no_pbc = @src.mda_guess_bonds(coords, elements)
  assert_eq(no_pbc.length(), 0)
}

///|
test "mda_cov_radius_values" {
  assert_true(@moonbit.assert_num_close(@src.mda_cov_radius("H"), 0.31, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(@src.mda_cov_radius("Zn"), 1.22, 1.0e-9),
  )
  assert_true(@moonbit.assert_num_close(@src.mda_cov_radius("I"), 1.39, 1.0e-9))
  // 未知元素使用默认半径 0.77
  assert_true(
    @moonbit.assert_num_close(@src.mda_cov_radius("Xx"), 0.77, 1.0e-9),
  )
}

///|
test "mda_bond_list_adjacency" {
  // 链式键 0-1-2-3 的邻接表
  let adj = @src.mda_bond_list_adjacency(4, [(0, 1), (1, 2), (2, 3)])
  assert_eq(adj.length(), 4)
  assert_eq(adj[0].length(), 1)
  assert_eq(adj[0][0], 1)
  assert_eq(adj[1].length(), 2)
  assert_eq(adj[1][0], 0)
  assert_eq(adj[1][1], 2)
  assert_eq(adj[2].length(), 2)
  assert_eq(adj[3].length(), 1)
  assert_eq(adj[3][0], 2)
}
