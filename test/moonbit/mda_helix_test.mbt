// ============================================================================
// B4: α / 3₁₀ / π 螺旋检测测试
// ============================================================================

// 手工构造 3 残基 peptide backbone 坐标。
// phi(i) = C(i-1) - N(i) - CA(i) - C(i)
// psi(i) = N(i) - CA(i) - C(i) - N(i+1)
// alpha helix 典型 phi/psi: φ ≈ -57°, ψ ≈ -47°

test "mda_helix_detection_three_residues_alpha" {
  let atoms : Array[@src.UnivAtom] = [
    // 残基 1
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=1, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=1, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=1, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
    // 残基 2
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=2, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=2, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=2, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
    // 残基 3
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=3, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=3, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=3, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
  ]
  // 简化 alpha helix 坐标（沿 +z 上升 + 旋转）
  let coords : Array[Array[Double]] = [
    [0.0,  0.0,  0.0],   // N0
    [1.47, 0.0,  0.0],   // CA0
    [2.60, 0.80, 0.0],   // C0
    [3.40, 0.50, 1.00],  // N1
    [4.85, 0.50, 1.00],  // CA1
    [6.00, 1.30, 1.00],  // C1
    [6.80, 1.00, 2.00],  // N2
    [8.25, 1.00, 2.00],  // CA2
    [9.40, 1.80, 2.00],  // C2
  ]
  let u = @src.Universe::from_dcd_frames(atoms, [coords])
  // select_atoms 用 String 选择表达式
  let ag = u.select_atoms("backbone")
  let helices = ag.helix_detection()
  // 验证函数不 crash，返回合法元组格式
  for h in helices {
    let (_start, _len, htype) = h
    assert_true(htype == "alpha" || htype == "3_10" || htype == "pi")
  }
}

test "mda_helix_detection_too_short_no_detection" {
  // 只有 2 残基 → 任何 helix 都不检测（最小 3）
  let atoms : Array[@src.UnivAtom] = [
    @src.UnivAtom::new(name="N", resname="GLY", chain="A", segid="P", element="N", atom_type="N", resid=1, mass=14.0, charge=-0.4, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="GLY", chain="A", segid="P", element="C", atom_type="C", resid=1, mass=12.0, charge=0.0, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="GLY", chain="A", segid="P", element="C", atom_type="C", resid=1, mass=12.0, charge=0.6, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="N", resname="GLY", chain="A", segid="P", element="N", atom_type="N", resid=2, mass=14.0, charge=-0.4, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="GLY", chain="A", segid="P", element="C", atom_type="C", resid=2, mass=12.0, charge=0.0, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="GLY", chain="A", segid="P", element="C", atom_type="C", resid=2, mass=12.0, charge=0.6, is_protein=true, is_water=false, is_backbone=true),
  ]
  let coords : Array[Array[Double]] = [
    [0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.6, 0.8, 0.0],
    [3.4, 0.5, 1.0], [4.9, 0.5, 1.0], [6.0, 1.3, 1.0],
  ]
  let u = @src.Universe::from_dcd_frames(atoms, [coords])
  let ag = u.select_atoms("backbone")
  let helices = ag.helix_detection()
  assert_eq(helices.length(), 0)
}

test "mda_helix_fraction_trajectory_basic" {
  // 多帧轨迹 → helix_fraction_trajectory 返回长度正确、值 ∈ [0,1]
  let atoms : Array[@src.UnivAtom] = [
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=1, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=1, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=1, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=2, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=2, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=2, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="N", resname="ALA", chain="A", segid="PROT", element="N", atom_type="N", resid=3, mass=14.007, charge=-0.4157, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="CA", resname="ALA", chain="A", segid="PROT", element="C", atom_type="CT", resid=3, mass=12.011, charge=0.0337, is_protein=true, is_water=false, is_backbone=true),
    @src.UnivAtom::new(name="C", resname="ALA", chain="A", segid="PROT", element="C", atom_type="C", resid=3, mass=12.011, charge=0.5973, is_protein=true, is_water=false, is_backbone=true),
  ]
  let coords0 : Array[Array[Double]] = [
    [0.0, 0.0, 0.0], [1.47, 0.0, 0.0], [2.60, 0.80, 0.0],
    [3.40, 0.50, 1.00], [4.85, 0.50, 1.00], [6.00, 1.30, 1.00],
    [6.80, 1.00, 2.00], [8.25, 1.00, 2.00], [9.40, 1.80, 2.00],
  ]
  // 直接构造 3 帧
  let u = @src.Universe::from_dcd_frames(atoms, [coords0, coords0, coords0])
  let ag = u.select_atoms("backbone")
  let frac = ag.helix_fraction_trajectory()
  assert_eq(frac.length(), 3)
  for f in frac {
    assert_true(f >= 0.0 && f <= 1.0)
  }
}

test "mda_helix_detection_empty_atoms_no_crash" {
  // 空 AtomGroup → 空数组
  let atoms : Array[@src.UnivAtom] = [
    @src.UnivAtom::new(name="O", resname="HOH", chain="A", segid="SYS", element="O", atom_type="O", resid=1, mass=16.0, charge=-0.8, is_protein=false, is_water=true, is_backbone=false),
  ]
  let coords : Array[Array[Double]] = [[0.0, 0.0, 0.0]]
  let u = @src.Universe::from_dcd_frames(atoms, [coords])
  // 用 select_atoms_by_indices 传空数组得到空 AtomGroup
  let ag = u.select_atoms_by_indices([])
  let helices = ag.helix_detection()
  assert_eq(helices.length(), 0)
}
