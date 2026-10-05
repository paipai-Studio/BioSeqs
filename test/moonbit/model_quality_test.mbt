///|
/// Helper: generate non-collinear CA-like coordinates (ideal alpha-helix).
fn make_helix_coords(n : Int) -> Array[@src.Vector3] {
  let coords = Array(capacity=n)
  let radius = 2.3
  let rise = 1.5
  let mut i = 0
  while i < n {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    let x = radius * @math.cos(angle)
    let y = radius * @math.sin(angle)
    let z = i.to_double() * rise
    coords.push(@src.Vector3::new(x, y, z))
    i = i + 1
  }
  coords
}

///|
test "model_quality: GDT_TS identity (same structure)" {
  let coords = make_helix_coords(10)
  let result = @src.calc_gdt_ts(coords, coords)
  assert_true(result.gdt_ts > 99.0)
  assert_true(result.gdt_ha > 99.0)
}

///|
test "model_quality: GDT_TS with rigid translation" {
  let ref = make_helix_coords(8)
  // Translate by (10, 20, 30)
  let model = Array(capacity=8)
  for v in ref {
    model.push(@src.Vector3::new(v.x + 10.0, v.y + 20.0, v.z + 30.0))
  }
  let result = @src.calc_gdt_ts(ref, model)
  assert_true(result.gdt_ts > 99.0)
}

///|
test "model_quality: GDT_TS with small perturbation" {
  let ref = make_helix_coords(8)
  // Perturb by < 0.5 Å each
  let model = Array(capacity=8)
  let mut i = 0
  while i < 8 {
    let v = ref[i]
    let dx = if i % 3 == 0 { 0.3 } else { 0.0 }
    let dy = if i % 3 == 1 { 0.3 } else { 0.0 }
    let dz = if i % 3 == 2 { 0.3 } else { 0.0 }
    model.push(@src.Vector3::new(v.x + dx, v.y + dy, v.z + dz))
    i = i + 1
  }
  let result = @src.calc_gdt_ts(ref, model)
  // Should have high GDT_TS since perturbation < 1 Å
  assert_true(result.gdt_ts > 50.0)
}

///|
test "model_quality: GDT_TS empty arrays" {
  let result = @src.calc_gdt_ts([], [])
  assert_eq(result.gdt_ts, 0.0)
}

///|
test "model_quality: GDT_TS mismatched lengths" {
  let ref = [@src.Vector3::new(0.0, 0.0, 0.0)]
  let model = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let result = @src.calc_gdt_ts(ref, model)
  assert_eq(result.gdt_ts, 0.0)
}

///|
test "model_quality: GDT_HA identity" {
  let coords = make_helix_coords(6)
  let result = @src.calc_gdt_ts(coords, coords)
  assert_true(result.gdt_ha > 99.0)
}

///|
test "model_quality: lDDT identity" {
  let coords = make_helix_coords(8)
  let result = @src.calc_lddt(coords, coords)
  assert_true((result.lddt - 1.0).abs() < 0.001)
}

///|
test "model_quality: lDDT with local distortion" {
  let ref = make_helix_coords(6)
  // Distort one atom by 3 Å
  let model = Array(capacity=6)
  let mut i = 0
  while i < 6 {
    let v = ref[i]
    if i == 2 {
      model.push(@src.Vector3::new(v.x, v.y + 3.0, v.z))
    } else {
      model.push(v)
    }
    i = i + 1
  }
  let result = @src.calc_lddt(ref, model)
  assert_true(result.lddt < 1.0)
  assert_true(result.lddt > 0.0)
}

///|
test "model_quality: lDDT per-residue" {
  let ref = make_helix_coords(5)
  let result = @src.calc_lddt(ref, ref)
  assert_eq(result.per_residue.length(), 5)
  for v in result.per_residue {
    assert_true((v - 1.0).abs() < 0.001)
  }
}

///|
test "model_quality: lDDT empty" {
  let result = @src.calc_lddt([], [])
  assert_eq(result.lddt, 0.0)
  assert_eq(result.per_residue.length(), 0)
}

///|
test "model_quality: DRMSD identity" {
  let coords = make_helix_coords(6)
  let result = @src.calc_drmsd(coords, coords)
  assert_true(result < 0.001)
}

///|
test "model_quality: DRMSD with known displacement" {
  let ref = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(10.0, 0.0, 0.0),
  ]
  // Move one atom by 1 Å perpendicular — changes one distance
  let model = [
    @src.Vector3::new(0.0, 1.0, 0.0),
    @src.Vector3::new(10.0, 0.0, 0.0),
  ]
  let result = @src.calc_drmsd(ref, model)
  // d_ref = 10, d_model = sqrt(100+1) ≈ 10.05
  // DRMSD = |10.05 - 10| ≈ 0.05
  assert_true(result > 0.01)
  assert_true(result < 0.2)
}

///|
test "model_quality: DRMSD empty" {
  assert_eq(@src.calc_drmsd([], []), 0.0)
}

///|
test "model_quality: DRMSD single point" {
  let coords = [@src.Vector3::new(1.0, 2.0, 3.0)]
  assert_eq(@src.calc_drmsd(coords, coords), 0.0)
}

///|
test "model_quality: DockQ perfect docking" {
  // Non-collinear receptor
  let receptor = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 3.0, 0.0),
  ]
  let ligand = [
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 3.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 3.0),
  ]
  let result = @src.calc_dockq(receptor, ligand, receptor, ligand)
  assert_true(result.fnat > 0.99)
  assert_true(result.dockq > 0.9)
}

///|
test "model_quality: DockQ poor docking" {
  let receptor_ref = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 3.0, 0.0),
  ]
  let ligand_ref = [
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 3.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 3.0),
  ]
  // Move ligand far away
  let ligand_model = [
    @src.Vector3::new(50.0, 0.0, 0.0),
    @src.Vector3::new(50.0, 3.0, 0.0),
    @src.Vector3::new(50.0, 0.0, 3.0),
  ]
  let result = @src.calc_dockq(receptor_ref, ligand_ref, receptor_ref, ligand_model)
  assert_true(result.fnat < 0.1)
  assert_true(result.dockq < 0.5)
}

///|
test "model_quality: DockQ empty" {
  let result = @src.calc_dockq([], [], [], [])
  assert_eq(result.dockq, 0.0)
}

///|
test "model_quality: DockQ mismatched lengths" {
  let rec = [@src.Vector3::new(0.0, 0.0, 0.0)]
  let lig = [@src.Vector3::new(5.0, 0.0, 0.0)]
  let rec2 = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 0.0, 0.0),
  ]
  let result = @src.calc_dockq(rec, lig, rec2, lig)
  assert_eq(result.dockq, 0.0)
}

///|
test "model_quality: GDT_TS with rotation" {
  let ref = make_helix_coords(8)
  // Rotate 90° around z: (x,y,z) → (-y,x,z)
  let model = Array(capacity=8)
  for v in ref {
    model.push(@src.Vector3::new(-v.y, v.x, v.z))
  }
  let result = @src.calc_gdt_ts(ref, model)
  assert_true(result.gdt_ts > 99.0)
}
