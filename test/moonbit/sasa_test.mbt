///|
/// Tests for Bio.PDB.SASA module.

///|
/// Helper: build an atom with explicit element on chain A.
fn make_atom(
  name : String,
  x : Double,
  y : Double,
  z : Double,
  element : String,
  resname : String,
  resseq : Int,
) -> @src.Atom {
  @src.Atom::new(
    name=name,
    coord=@src.Vector3::new(x, y, z),
    resname=resname,
    chainid='A',
    resseq=resseq,
    element=element,
  )
}

///|
/// Helper: place all given atoms into a single residue on chain A.
fn make_structure(atoms : Array[@src.Atom]) -> @src.Structure {
  let res = @src.Residue::new(resname="UNK", chainid='A', resseq=1, atoms=atoms)
  let chain = @src.Chain::new(id='A', residues=[res])
  let model = @src.Model::new(id=1, chains=[chain])
  @src.Structure::new(id="test", models=[model])
}

///|
/// Pi constant matching the implementation.
let sasa_pi : Double = 3.14159265358979323846

///|
test "sasa_vdw_radius_h" {
  assert_eq(@src.sasa_vdw_radius("H"), 1.20)
}

///|
test "sasa_vdw_radius_c" {
  assert_eq(@src.sasa_vdw_radius("C"), 1.70)
}

///|
test "sasa_vdw_radius_n" {
  assert_eq(@src.sasa_vdw_radius("N"), 1.55)
}

///|
test "sasa_vdw_radius_o" {
  assert_eq(@src.sasa_vdw_radius("O"), 1.52)
}

///|
test "sasa_vdw_radius_s_and_p" {
  assert_eq(@src.sasa_vdw_radius("S"), 1.80)
  assert_eq(@src.sasa_vdw_radius("P"), 1.80)
}

///|
test "sasa_vdw_radius_default" {
  assert_eq(@src.sasa_vdw_radius("X"), 1.70)
  assert_eq(@src.sasa_vdw_radius(""), 1.70)
}

///|
test "sasa_vdw_radius_case_insensitive" {
  assert_eq(@src.sasa_vdw_radius("h"), 1.20)
  assert_eq(@src.sasa_vdw_radius("c"), 1.70)
  assert_eq(@src.sasa_vdw_radius("cl"), 1.75)
}

///|
test "sasa_vdw_radius_halogen" {
  assert_eq(@src.sasa_vdw_radius("F"), 1.47)
  assert_eq(@src.sasa_vdw_radius("Cl"), 1.75)
  assert_eq(@src.sasa_vdw_radius("Br"), 1.85)
  assert_eq(@src.sasa_vdw_radius("I"), 1.98)
}

///|
test "sasa_generate_sphere_points_count" {
  let pts = @src.sasa_generate_sphere_points(100)
  assert_eq(pts.length(), 100)
}

///|
test "sasa_generate_sphere_points_unit_length" {
  let pts = @src.sasa_generate_sphere_points(100)
  for p in pts {
    let norm_sq = p.x * p.x + p.y * p.y + p.z * p.z
    assert_true((norm_sq - 1.0).abs() < 1.0e-9)
  }
}

///|
test "sasa_generate_sphere_points_n1" {
  let pts = @src.sasa_generate_sphere_points(1)
  assert_eq(pts.length(), 1)
  let p = pts[0]
  let norm_sq = p.x * p.x + p.y * p.y + p.z * p.z
  assert_true((norm_sq - 1.0).abs() < 1.0e-9)
}

///|
test "sasa_generate_sphere_points_n0_clamped" {
  // n <= 0 is clamped to 1.
  let pts = @src.sasa_generate_sphere_points(0)
  assert_eq(pts.length(), 1)
}

///|
test "sasa_single_atom_sasa" {
  let s = make_structure([make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1)])
  let res = @src.sasa_calc(s, 200, 1.4)
  assert_eq(res.length(), 1)
  let r = 1.70 + 1.4
  let expected = 4.0 * sasa_pi * r * r
  // Isolated atom: all points accessible -> exact 4*pi*r^2.
  assert_true((res[0].sasa - expected).abs() < 1.0e-6)
}

///|
test "sasa_single_atom_probe_zero" {
  let s = make_structure([make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1)])
  let res = @src.sasa_calc(s, 200, 0.0)
  let r = 1.70
  let expected = 4.0 * sasa_pi * r * r
  assert_true((res[0].sasa - expected).abs() < 1.0e-6)
}

///|
test "sasa_two_atoms_far_apart" {
  let s = make_structure([
    make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1),
    make_atom("CB", 50.0, 0.0, 0.0, "C", "ALA", 1),
  ])
  let res = @src.sasa_calc(s, 200, 1.4)
  let r = 1.70 + 1.4
  let single = 4.0 * sasa_pi * r * r
  let total = res[0].sasa + res[1].sasa
  // No overlap: each atom fully accessible.
  assert_true((total - 2.0 * single).abs() < 1.0e-6)
}

///|
test "sasa_two_atoms_close_reduced" {
  let close = make_structure([
    make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1),
    make_atom("CB", 1.0, 0.0, 0.0, "C", "ALA", 1),
  ])
  let far = make_structure([
    make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1),
    make_atom("CB", 50.0, 0.0, 0.0, "C", "ALA", 1),
  ])
  let close_res = @src.sasa_calc(close, 200, 1.4)
  let far_res = @src.sasa_calc(far, 200, 1.4)
  let close_total = close_res[0].sasa + close_res[1].sasa
  let far_total = far_res[0].sasa + far_res[1].sasa
  assert_true(close_total < far_total)
  assert_true(close_total > 0.0)
}

///|
test "sasa_empty_structure" {
  let s = @src.Structure::new(id="empty", models=[@src.Model::new(id=1)])
  let res = @src.sasa_calc(s, 100, 1.4)
  assert_eq(res.length(), 0)
  let summary = @src.sasa_calc_total(res)
  assert_eq(summary.get_num_atoms(), 0)
  assert_eq(summary.get_num_residues(), 0)
  assert_true(summary.get_total_sasa() == 0.0)
}

///|
test "sasa_zero_points" {
  let s = make_structure([
    make_atom("CA", 0.0, 0.0, 0.0, "C", "ALA", 1),
    make_atom("CB", 1.0, 0.0, 0.0, "C", "ALA", 1),
  ])
  let res = @src.sasa_calc(s, 0, 1.4)
  assert_eq(res.length(), 2)
  for r in res {
    assert_true(r.sasa == 0.0)
  }
}

///|
test "sasa_residue_sasa_basic" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let res_res = @src.sasa_calc_residue(atom_res, s)
  assert_eq(res_res.length(), 3)
  assert_eq(res_res[0].get_resname(), "ALA")
  assert_eq(res_res[1].get_resname(), "GLY")
  assert_eq(res_res[2].get_resname(), "SER")
  assert_eq(res_res[0].get_resseq(), 1)
  assert_eq(res_res[0].get_chainid(), 'A')
}

///|
test "sasa_residue_backbone_sidechain_split" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let res_res = @src.sasa_calc_residue(atom_res, s)
  for r in res_res {
    assert_true(
      (r.get_backbone_sasa() + r.get_sidechain_sasa() - r.get_total_sasa()).abs() <
      1.0e-9,
    )
  }
  // GLY has only backbone atoms (N, CA, C, O): sidechain must be zero.
  let gly = res_res[1]
  assert_true(gly.get_sidechain_sasa() == 0.0)
  assert_true((gly.get_backbone_sasa() - gly.get_total_sasa()).abs() < 1.0e-9)
  // ALA and SER have sidechain atoms.
  assert_true(res_res[0].get_sidechain_sasa() > 0.0)
  assert_true(res_res[2].get_sidechain_sasa() > 0.0)
}

///|
test "sasa_total_sasa_sum" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let summary = @src.sasa_calc_total(atom_res)
  let mut sum = 0.0
  for r in atom_res {
    sum = sum + r.sasa
  }
  assert_true((summary.get_total_sasa() - sum).abs() < 1.0e-9)
}

///|
test "sasa_total_backbone_sidechain" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let summary = @src.sasa_calc_total(atom_res)
  assert_true(
    (summary.get_backbone_sasa() + summary.get_sidechain_sasa() -
      summary.get_total_sasa()).abs() <
    1.0e-9,
  )
  assert_true(summary.get_backbone_sasa() > 0.0)
  assert_true(summary.get_sidechain_sasa() > 0.0)
}

///|
test "sasa_total_num_residues" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let summary = @src.sasa_calc_total(atom_res)
  assert_eq(summary.get_num_atoms(), s.get_num_atoms())
  assert_eq(summary.get_num_residues(), 3)
}

///|
test "sasa_sample_structure_shape" {
  let s = @src.sasa_sample_structure()
  assert_eq(s.get_id(), "sample")
  assert_eq(s.get_num_models(), 1)
  assert_eq(s.get_num_chains(), 1)
  assert_eq(s.get_num_residues(), 3)
  // ALA(5) + GLY(4) + SER(6) = 15 atoms.
  assert_eq(s.get_num_atoms(), 15)
}

///|
test "sasa_sample_structure_sasa_positive" {
  let s = @src.sasa_sample_structure()
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  assert_eq(atom_res.length(), s.get_num_atoms())
  for r in atom_res {
    assert_true(r.sasa >= 0.0)
  }
  let summary = @src.sasa_calc_total(atom_res)
  assert_true(summary.get_total_sasa() > 0.0)
}

///|
test "sasa_single_residue_structure" {
  let s = make_structure([
    make_atom("N", 0.0, 0.0, 0.0, "N", "ALA", 1),
    make_atom("CA", 1.46, 0.0, 0.0, "C", "ALA", 1),
    make_atom("C", 2.01, 1.42, 0.0, "C", "ALA", 1),
    make_atom("O", 1.23, 2.36, 0.0, "O", "ALA", 1),
    make_atom("CB", 1.05, -0.79, 1.21, "C", "ALA", 1),
  ])
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let res_res = @src.sasa_calc_residue(atom_res, s)
  assert_eq(res_res.length(), 1)
  assert_eq(res_res[0].get_num_atoms(), 5)
  let summary = @src.sasa_calc_total(atom_res)
  assert_eq(summary.get_num_residues(), 1)
}

///|
test "sasa_all_backbone_atoms_sidechain_zero" {
  // Residue with only N, CA, C, O: sidechain SASA must be zero.
  let s = make_structure([
    make_atom("N", 0.0, 0.0, 0.0, "N", "GLY", 1),
    make_atom("CA", 1.46, 0.0, 0.0, "C", "GLY", 1),
    make_atom("C", 2.01, 1.42, 0.0, "C", "GLY", 1),
    make_atom("O", 1.23, 2.36, 0.0, "O", "GLY", 1),
  ])
  let atom_res = @src.sasa_calc(s, 100, 1.4)
  let res_res = @src.sasa_calc_residue(atom_res, s)
  assert_eq(res_res.length(), 1)
  assert_true(res_res[0].get_sidechain_sasa() == 0.0)
  assert_true((res_res[0].get_backbone_sasa() - res_res[0].get_total_sasa()).abs() < 1.0e-9)
}

///|
test "sasa_result_accessors" {
  let r = @src.SasaResult::new(
    atom_index=3,
    atom_name="CA",
    resname="ALA",
    chainid='A',
    resseq=1,
    sasa=42.5,
  )
  assert_eq(r.get_atom_index(), 3)
  assert_eq(r.get_atom_name(), "CA")
  assert_eq(r.get_resname(), "ALA")
  assert_eq(r.get_chainid(), 'A')
  assert_eq(r.get_resseq(), 1)
  assert_true((r.get_sasa() - 42.5).abs() < 1.0e-9)
}

///|
test "sasa_residue_sasa_accessors" {
  let r = @src.ResidueSasa::new(
    resname="SER",
    chainid='B',
    resseq=7,
    total_sasa=100.0,
    backbone_sasa=40.0,
    sidechain_sasa=60.0,
    num_atoms=6,
  )
  assert_eq(r.get_resname(), "SER")
  assert_eq(r.get_chainid(), 'B')
  assert_eq(r.get_resseq(), 7)
  assert_true((r.get_total_sasa() - 100.0).abs() < 1.0e-9)
  assert_true((r.get_backbone_sasa() - 40.0).abs() < 1.0e-9)
  assert_true((r.get_sidechain_sasa() - 60.0).abs() < 1.0e-9)
  assert_eq(r.get_num_atoms(), 6)
}

///|
test "sasa_summary_accessors" {
  let s = @src.SasaSummary::new(
    total_sasa=250.0,
    backbone_sasa=100.0,
    sidechain_sasa=150.0,
    num_atoms=20,
    num_residues=4,
  )
  assert_true((s.get_total_sasa() - 250.0).abs() < 1.0e-9)
  assert_true((s.get_backbone_sasa() - 100.0).abs() < 1.0e-9)
  assert_true((s.get_sidechain_sasa() - 150.0).abs() < 1.0e-9)
  assert_eq(s.get_num_atoms(), 20)
  assert_eq(s.get_num_residues(), 4)
}

///|
test "sasa_atom_result_indexing_matches_atoms" {
  let s = @src.sasa_sample_structure()
  let atoms = s.get_atoms()
  let res = @src.sasa_calc(s, 50, 1.4)
  assert_eq(res.length(), atoms.length())
  for i in 0..<res.length() {
    assert_eq(res[i].atom_index, i)
    assert_eq(res[i].atom_name, atoms[i].name)
    assert_eq(res[i].resname, atoms[i].resname)
    assert_eq(res[i].resseq, atoms[i].resseq)
  }
}
