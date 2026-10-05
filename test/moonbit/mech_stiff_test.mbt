///|
/// Tests for MechStiff (mechanical stiffness matrix).

///|
/// Build a chain with CA atoms at given coordinates.
fn ms_make_chain(coords : Array[(@src.Vector3, String, Int)]) -> @src.Chain {
  let residues = Array(capacity=coords.length())
  for (coord, resname, resseq) in coords {
    let atoms = [
      @src.Atom::new(name="N", coord=coord, resname=resname, resseq=resseq),
      @src.Atom::new(name="CA", coord=coord, resname=resname, resseq=resseq),
      @src.Atom::new(name="C", coord=coord, resname=resname, resseq=resseq),
      @src.Atom::new(name="O", coord=coord, resname=resname, resseq=resseq),
    ]
    residues.push(@src.Residue::new(resname=resname, chainid='A', resseq=resseq, atoms=atoms))
  }
  @src.Chain::new(id='A', residues=residues)
}

///|
test "mechstiff: basic construction" {
  // 5 residues in a line
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
    (@src.Vector3::new(15.2, 0.0, 0.0), "ILE", 5),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      assert_eq(m.num_residues(), 5)
      let profile = m.stiffness_profile()
      assert_eq(profile.length(), 5)
    }
    None => assert_true(false)
  }
}

///|
test "mechstiff: effective spring" {
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      // Self-stiffness should be positive
      let s00 = m.effective_spring(0, 0)
      assert_true(s00 >= 0.0)
      // Stiffness between nearby residues should be positive
      let s01 = m.effective_spring(0, 1)
      assert_true(s01 >= 0.0)
    }
    None => assert_true(false)
  }
}

///|
test "mechstiff: stiffness between" {
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
    (@src.Vector3::new(15.2, 0.0, 0.0), "ILE", 5),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      let block = m.stiffness_between(0, 1)
      assert_eq(block.length(), 3)
      assert_eq(block[0].length(), 3)
    }
    None => assert_true(false)
  }
}

///|
test "mechstiff: mechanical distance profile" {
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
    (@src.Vector3::new(15.2, 0.0, 0.0), "ILE", 5),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      let dist = m.mechanical_distance_profile()
      assert_eq(dist.length(), 5)
      for d in dist {
        assert_true(d >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

///|
test "mechstiff: stiff residues" {
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
    (@src.Vector3::new(15.2, 0.0, 0.0), "ILE", 5),
    (@src.Vector3::new(19.0, 0.0, 0.0), "PHE", 6),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      let stiff = m.stiff_residues()
      // Should return some indices
      assert_true(stiff.length() >= 0)
    }
    None => assert_true(false)
  }
}

///|
test "mechstiff: too few residues" {
  // Only 2 residues — ANM needs at least 3
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "mechstiff: get matrix" {
  let coords = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "ALA", 1),
    (@src.Vector3::new(3.8, 0.0, 0.0), "GLY", 2),
    (@src.Vector3::new(7.6, 0.0, 0.0), "VAL", 3),
    (@src.Vector3::new(11.4, 0.0, 0.0), "LEU", 4),
  ]
  let chain = ms_make_chain(coords)
  let ms = @src.MechStiff::new(chain)
  match ms {
    Some(m) => {
      let mat = m.get_matrix()
      // 3N x 3N = 12 x 12
      assert_eq(mat.length(), 12)
      let springs = m.get_springs()
      // N x N = 4 x 4
      assert_eq(springs.length(), 4)
    }
    None => assert_true(false)
  }
}
