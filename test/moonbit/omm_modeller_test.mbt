///|
test "modeller_add_solvent_fills_box" {
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(
        name="CA",
        resname="CAL",
        element="Ca",
        mass=40.08,
        charge=2.0,
      ),
      position=[2.0, 2.0, 2.0],
    ),
  )
  m.set_box(8.0, 8.0, 8.0)
  let n_water = m.add_solvent(spacing=3.1, min_dist=2.0)
  // Lattice 1.55/4.65 x 1.55/4.65 x 1.55/4.65 = 8 sites; the site nearest
  // the solute (1.55,1.55,1.55) is rejected by the min_dist filter.
  assert_eq(n_water, 7)
  assert_eq(m.num_atoms(), 1 + 3 * n_water)
  // Every water oxygen is inside the box and at least min_dist from the
  // solute under the minimum-image convention.
  let bx = m.box_lengths()
  let solute = [2.0, 2.0, 2.0]
  let atoms = m.atom_list()
  let pos = m.positions()
  for i = 0; i < m.num_atoms(); i = i + 1 {
    if atoms[i].resname == "HOH" && atoms[i].name == "O" {
      let p = pos[i]
      assert_true(p[0] >= 0.0 && p[0] < bx[0])
      assert_true(p[1] >= 0.0 && p[1] < bx[1])
      assert_true(p[2] >= 0.0 && p[2] < bx[2])
      let d = @src.md_pbc_distance(p, solute, bx)
      assert_true(d >= 2.0 - 1.0e-9)
    }
  }
  // Net charge of the solvated (non-neutralized) model is conserved.
  let q = m.net_charge()
  assert_true(@moonbit.assert_num_close(q, 2.0 + 7.0 * 0.0, 1.0e-9))
}

///|
test "modeller_add_solvent_derives_box" {
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(name="X", element="C", mass=12.011),
      position=[0.0, 0.0, 0.0],
    ),
  )
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(name="Y", element="C", mass=12.011),
      position=[4.0, 2.0, 1.0],
    ),
  )
  let n = m.add_solvent(spacing=2.0, min_dist=1.5, padding=1.0)
  assert_true(n > 0)
  let bx = m.box_lengths()
  // Bounding box (0..4, 0..2, 0..1) + 2*padding = (6, 4, 3)
  assert_true(@moonbit.assert_num_close(bx[0], 6.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(bx[1], 4.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(bx[2], 3.0, 1.0e-9))
}

///|
test "modeller_center_in_box" {
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(atom=@src.OmmModelAtom::new(name="A", element="C", mass=12.0), position=[
      8.0, 8.0, 8.0,
    ]),
  )
  ignore(
    m.add_atom(atom=@src.OmmModelAtom::new(name="B", element="C", mass=12.0), position=[
      10.0, 10.0, 10.0,
    ]),
  )
  m.set_box(10.0, 10.0, 10.0)
  m.center_in_box()
  let p = m.positions()
  assert_true(@moonbit.assert_num_close((p[0][0] + p[1][0]) / 2.0, 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close((p[0][1] + p[1][1]) / 2.0, 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close((p[0][2] + p[1][2]) / 2.0, 5.0, 1.0e-9))
  // No-op without a box.
  let m2 = @src.OmmModeller::new()
  ignore(
    m2.add_atom(atom=@src.OmmModelAtom::new(name="A", element="C", mass=12.0), position=[
      1.0, 2.0, 3.0,
    ]),
  )
  m2.center_in_box()
  assert_true(@moonbit.assert_num_close(m2.positions()[0][0], 1.0, 1.0e-9))
}

///|
test "modeller_delete_water" {
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(
        name="CL",
        resname="CL",
        element="Cl",
        mass=35.453,
        charge=-1.0,
      ),
      position=[0.0, 0.0, 0.0],
    ),
  )
  // Two explicit waters.
  let ws : Array[Array[Double]] = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]]
  for k = 0; k < 2; k = k + 1 {
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="O",
          resname="HOH",
          element="O",
          mass=15.9994,
          charge=-0.834,
        ),
        position=ws[k],
      ),
    )
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="H1",
          resname="HOH",
          element="H",
          mass=1.008,
          charge=0.417,
        ),
        position=[ws[k][0], ws[k][1] + 0.757, ws[k][2] + 0.587],
      ),
    )
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="H2",
          resname="HOH",
          element="H",
          mass=1.008,
          charge=0.417,
        ),
        position=[ws[k][0], ws[k][1] - 0.757, ws[k][2] + 0.587],
      ),
    )
  }
  assert_eq(m.num_atoms(), 7)
  assert_eq(m.delete_water(), 6)
  assert_eq(m.num_atoms(), 1)
}

///|
test "modeller_neutralize_replaces_water_with_ions" {
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(
        name="CL",
        resname="CL",
        element="Cl",
        mass=35.453,
        charge=-1.0,
      ),
      position=[0.0, 0.0, 0.0],
    ),
  )
  // One water near the ion site to be replaced.
  let wpos = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0]]
  for k = 0; k < 2; k = k + 1 {
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="O",
          resname="HOH",
          element="O",
          mass=15.9994,
          charge=-0.834,
        ),
        position=wpos[k],
      ),
    )
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="H1",
          resname="HOH",
          element="H",
          mass=1.008,
          charge=0.417,
        ),
        position=[wpos[k][0], wpos[k][1] + 0.757, wpos[k][2] + 0.587],
      ),
    )
    ignore(
      m.add_atom(
        atom=@src.OmmModelAtom::new(
          name="H2",
          resname="HOH",
          element="H",
          mass=1.008,
          charge=0.417,
        ),
        position=[wpos[k][0], wpos[k][1] - 0.757, wpos[k][2] + 0.587],
      ),
    )
  }
  let ions = m.neutralize()
  assert_eq(ions, 1)
  assert_true(@moonbit.assert_num_close(m.net_charge(), 0.0, 1.0e-9))
  // One water lost its two hydrogens; the ion is named NA.
  let atoms = m.atom_list()
  let mut n_na = 0
  let mut n_h = 0
  for a in atoms {
    if a.name == "NA" {
      n_na = n_na + 1
    }
    if a.element == "H" {
      n_h = n_h + 1
    }
  }
  assert_eq(n_na, 1)
  assert_eq(n_h, 2)
  // Neutralizing an already neutral model is a no-op.
  assert_eq(m.neutralize(), 0)
}

///|
test "modeller_add_hydrogens_valence" {
  // Two bonded carbons: each needs 3 hydrogens.
  let m = @src.OmmModeller::new()
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(name="C1", element="C", mass=12.011),
      position=[0.0, 0.0, 0.0],
    ),
  )
  ignore(
    m.add_atom(
      atom=@src.OmmModelAtom::new(name="C2", element="C", mass=12.011),
      position=[1.54, 0.0, 0.0],
    ),
  )
  let added = m.add_hydrogens()
  assert_eq(added, 6)
  assert_eq(m.num_atoms(), 8)
  // All C-H distances ~1.09 A; no hydrogen collisions.
  let pos = m.positions()
  let atoms = m.atom_list()
  for i = 0; i < m.num_atoms(); i = i + 1 {
    if atoms[i].element == "H" {
      let d0 = @src.md_dist(pos[i], pos[0])
      let d1 = @src.md_dist(pos[i], pos[1])
      let near = if d0 < d1 { d0 } else { d1 }
      assert_true(@moonbit.assert_num_close(near, 1.09, 0.05))
    }
  }
  for i = 2; i < m.num_atoms(); i = i + 1 {
    for j = i + 1; j < m.num_atoms(); j = j + 1 {
      let d = @src.md_dist(pos[i], pos[j])
      assert_true(d > 0.3)
    }
  }
  // Isolated oxygen: valence 2 -> two hydrogens.
  let m2 = @src.OmmModeller::new()
  ignore(
    m2.add_atom(
      atom=@src.OmmModelAtom::new(name="O", element="O", mass=15.999),
      position=[0.0, 0.0, 0.0],
    ),
  )
  assert_eq(m2.add_hydrogens(), 2)
}

///|
test "modeller_add_atom_resid_tracking" {
  let m = @src.OmmModeller::new()
  let i0 = m.add_atom(
    atom=@src.OmmModelAtom::new(
      name="A",
      resname="ALA",
      resid=3,
      element="C",
      mass=12.0,
    ),
    position=[0.0, 0.0, 0.0],
  )
  assert_eq(i0, 0)
  assert_eq(m.num_atoms(), 1)
  // box_lengths starts as zeros (non-periodic).
  let bx = m.box_lengths()
  assert_true(@moonbit.assert_num_close(bx[0], 0.0, 1.0e-12))
}
