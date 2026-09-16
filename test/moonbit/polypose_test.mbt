///|
/// Tests for polypose module (Bio.PDB.Polypose - multi-conformer structures).

fn ppy_example_conformer_1() -> @src.PpyConformer {
  let a1 = @src.PpyAtom::new("CA", 1.0, 0.0, 0.0)
  let a2 = @src.PpyAtom::new("CA", 2.0, 0.0, 0.0)
  let a3 = @src.PpyAtom::new("CA", 3.0, 0.0, 0.0)
  let r1 = @src.PpyResidue::new("ALA", 1, [a1])
  let r2 = @src.PpyResidue::new("GLY", 2, [a2])
  let r3 = @src.PpyResidue::new("LEU", 3, [a3])
  @src.PpyConformer::new("conf1", [r1, r2, r3])
}

fn ppy_example_conformer_2() -> @src.PpyConformer {
  let a1 = @src.PpyAtom::new("CA", 1.0, 0.1, 0.0)
  let a2 = @src.PpyAtom::new("CA", 2.1, 0.0, 0.0)
  let a3 = @src.PpyAtom::new("CA", 2.9, -0.1, 0.0)
  let r1 = @src.PpyResidue::new("ALA", 1, [a1])
  let r2 = @src.PpyResidue::new("GLY", 2, [a2])
  let r3 = @src.PpyResidue::new("LEU", 3, [a3])
  @src.PpyConformer::new("conf2", [r1, r2, r3])
}

test "ppy_conformer_residues" {
  let c = ppy_example_conformer_1()
  let rs = @src.ppy_conformer_residues(c)
  assert_eq(rs.length(), 3)
  assert_eq(rs[0].resname, "ALA")
}

test "ppy_find_residue" {
  let c = ppy_example_conformer_1()
  match @src.ppy_find_residue(c, 2) {
    Some(r) => assert_eq(r.resname, "GLY")
    None => assert_true(false)
  }
}

test "ppy_atom" {
  let c = ppy_example_conformer_1()
  match @src.ppy_atom(@src.ppy_find_residue(c, 1).unwrap(), "CA") {
    Some(a) => assert_true((a.x - 1.0).abs() < 1.0e-9)
    None => assert_true(false)
  }
}

test "ppy_coords" {
  let c = ppy_example_conformer_1()
  let cs = @src.ppy_coords(c, "CA")
  assert_eq(cs.length(), 3)
  assert_true((cs[0].x - 1.0).abs() < 1.0e-9)
  assert_true((cs[2].x - 3.0).abs() < 1.0e-9)
}

test "ppy_kabsch_rmsd" {
  let a = @src.ppy_coords(ppy_example_conformer_1(), "CA")
  let b = @src.ppy_coords(ppy_example_conformer_2(), "CA")
  let d = @src.ppy_kabsch_rmsd(a, b)
  assert_true(d >= 0.0)
  assert_true(d < 0.5)
}

test "ppy_pairwise_rmsds" {
  let st = @src.ppy_new("nmr")
  let st1 = @src.ppy_add_conformer(st, ppy_example_conformer_1())
  let st2 = @src.ppy_add_conformer(st1, ppy_example_conformer_2())
  let mat_opt = @src.ppy_pairwise_rmsds(st2)
  match mat_opt {
    Some(mat) => {
      assert_eq(mat.length(), 2)
      assert_eq(mat[0].length(), 2)
      assert_true((mat[0][1] - mat[1][0]).abs() < 1.0e-9)
    }
    None => assert_true(false)
  }
}

test "ppy_convergence" {
  let st = @src.ppy_new("nmr")
  let st1 = @src.ppy_add_conformer(st, ppy_example_conformer_1())
  let st2 = @src.ppy_add_conformer(st1, ppy_example_conformer_2())
  match @src.ppy_convergence(st2) {
    Some(c) => {
      assert_eq(c.n_conformers, 2)
      assert_eq(c.n_calpha, 3)
      assert_true(c.mean_rmsd >= 0.0)
    }
    None => assert_true(false)
  }
}

test "ppy_altloc_residue_count" {
  let a1 = @src.PpyAtom::new("CA", 1.0, 0.0, 0.0, altloc="A")
  let a2 = @src.PpyAtom::new("CA", 1.0, 0.1, 0.0, altloc="B")
  let r = @src.PpyResidue::new("ALA", 1, [a1, a2])
  let n = @src.ppy_altloc_residue_count(@src.PpyConformer::new("c", [r]))
  assert_eq(n, 1)
}

test "ppy_structure_coords_consistent" {
  let st = @src.ppy_new("nmr")
  let st1 = @src.ppy_add_conformer(st, ppy_example_conformer_1())
  let st2 = @src.ppy_add_conformer(st1, ppy_example_conformer_2())
  let opt = @src.ppy_structure_coords(st2, "CA")
  match opt {
    Some(coords) => {
      assert_eq(coords.length(), 2)
      assert_eq(coords[0].length(), 3)
    }
    None => assert_true(false)
  }
}
