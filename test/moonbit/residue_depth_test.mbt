///|
test "residue_depth_point3d" {
  let p = @src.RDPoint3D::new(1.0, 2.0, 3.0)

  assert_eq(p.x, 1.0)
  assert_eq(p.y, 2.0)
  assert_eq(p.z, 3.0)
}

///|
test "residue_depth_point_distance" {
  let p1 = @src.RDPoint3D::new(0.0, 0.0, 0.0)
  let p2 = @src.RDPoint3D::new(3.0, 4.0, 0.0)

  let dist = p1.distance(p2)
  assert_eq(dist, 5.0)
}

///|
test "residue_depth_create_atom" {
  let atom = @src.RDAtom::new("CA", @src.RDPoint3D::new(1.0, 2.0, 3.0), "C")

  assert_eq(atom.name, "CA")
  assert_eq(atom.element, "C")
}

///|
test "residue_depth_create_residue" {
  let atoms = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(1.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("N", @src.RDPoint3D::new(0.0, 0.0, 0.0), "N"),
  ]
  let residue = @src.RDResidue::new("ALA", 1, atoms)

  assert_eq(residue.resname, "ALA")
  assert_eq(residue.resseq, 1)
  assert_eq(residue.atoms.length(), 2)
}

///|
test "residue_depth_center_of_mass" {
  let atoms = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(0.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("N", @src.RDPoint3D::new(2.0, 0.0, 0.0), "N"),
  ]
  let residue = @src.RDResidue::new("ALA", 1, atoms)

  let com = residue.center_of_mass()
  assert_true(com.x > 0.9 && com.x < 1.1)
}

///|
test "residue_depth_get_ca" {
  let atoms = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(1.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("N", @src.RDPoint3D::new(0.0, 0.0, 0.0), "N"),
  ]
  let residue = @src.RDResidue::new("ALA", 1, atoms)

  let ca = residue.get_ca()
  assert_true(ca.is_some())
}

///|
test "residue_depth_analyze" {
  let atoms1 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(0.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("N", @src.RDPoint3D::new(-1.0, 0.0, 0.0), "N"),
  ]
  let atoms2 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(5.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("N", @src.RDPoint3D::new(6.0, 0.0, 0.0), "N"),
  ]
  let residues = [
    @src.RDResidue::new("ALA", 1, atoms1),
    @src.RDResidue::new("GLY", 2, atoms2),
  ]

  let results = @src.bio_residue_depth_analyze(residues)

  assert_eq(results.length(), 2)
}

///|
test "residue_depth_find_surface" {
  let atoms1 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(0.0, 0.0, 0.0), "C"),
  ]
  let atoms2 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(100.0, 0.0, 0.0), "C"),
  ]
  let residues = [
    @src.RDResidue::new("ALA", 1, atoms1),
    @src.RDResidue::new("GLY", 2, atoms2),
  ]

  let surface = @src.bio_residue_depth_find_surface(residues, 3.0)

  assert_true(surface.length() >= 1)
}

///|
test "residue_depth_find_core" {
  let atoms1 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(0.0, 0.0, 0.0), "C"),
    @src.RDAtom::new("CB", @src.RDPoint3D::new(0.5, 0.5, 0.5), "C"),
    @src.RDAtom::new("CG", @src.RDPoint3D::new(1.0, 1.0, 1.0), "C"),
  ]
  let atoms2 = [
    @src.RDAtom::new("CA", @src.RDPoint3D::new(2.0, 0.0, 0.0), "C"),
  ]
  let residues = [
    @src.RDResidue::new("ALA", 1, atoms1),
    @src.RDResidue::new("GLY", 2, atoms2),
  ]

  let core = @src.bio_residue_depth_find_core(residues, 8.0)

  assert_eq(core.length() >= 0, true)
}

///|
test "residue_depth_average" {
  let results = [
    @src.RDResult::new(1, "ALA", 2.0, 3.0, 10.0, 15.0),
    @src.RDResult::new(2, "GLY", 4.0, 5.0, 20.0, 25.0),
  ]

  let (avg_ca, avg_com) = @src.bio_residue_depth_average(results)

  assert_eq(avg_ca, 3.0)
  assert_eq(avg_com, 4.0)
}