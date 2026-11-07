///|
/// Tests for CHARMM/XPLOR/NAMD PSF topology parsing (structure_psf.mbt).

///|
test "psf_parse_complete" {
  let text = "PSF EXT\n\n       2 !NTITLE\n REMARKS example topology\n REMARKS second line\n\n       4 !NATOM\n   1 A     1  ALA  N    NH1   -0.415700       14.0100           0\n   2 A     1  ALA  H    H      0.271900        1.0080           0\n   3 A     1  ALA  CA   CT1    0.021300       12.0110           0\n   4 B     2  HOH  O    OH2   -0.834000       16.0000           0\n\n       3 !NBOND: bonds\n   1     2     2     3     3     4\n\n       2 !NTHETA: angles\n   1     2     3     2     3     4\n\n       2 !NPHI: dihedrals\n   1     2     3     4     2     3     4     1\n\n       1 !NIMPHI: impropers\n   1     2     3     4\n\n       1 !NDON: donors\n   1     2\n\n       1 !NACC: acceptors\n   3     4\n"
  let topo = @src.psf_parse(text)
  assert_eq(topo.natoms(), 4)
  assert_eq(topo.titles.length(), 2)
  assert_eq(topo.titles[0], "REMARKS example topology")
  assert_eq(topo.titles[1], "REMARKS second line")
  assert_eq(topo.atoms[0].atomid, 1)
  assert_eq(topo.atoms[0].segid, "A")
  assert_eq(topo.atoms[0].resid, 1)
  assert_eq(topo.atoms[0].resname, "ALA")
  assert_eq(topo.atoms[0].atom_name, "N")
  assert_eq(topo.atoms[0].atom_type, "NH1")
  assert_eq(topo.atoms[3].segid, "B")
  assert_eq(topo.atoms[3].resname, "HOH")
  assert_eq(topo.atoms[3].atom_type, "OH2")
  // 3 bond pairs streamed from a single physical line
  assert_eq(topo.bonds.length(), 3)
  assert_eq(topo.bonds[0], (1, 2))
  assert_eq(topo.bonds[1], (2, 3))
  assert_eq(topo.bonds[2], (3, 4))
  // 2 angle triples streamed from a single physical line
  assert_eq(topo.angles.length(), 2)
  assert_eq(topo.angles[0], (1, 2, 3))
  assert_eq(topo.angles[1], (2, 3, 4))
  assert_eq(topo.dihedrals.length(), 2)
  assert_eq(topo.dihedrals[0], (1, 2, 3, 4))
  assert_eq(topo.dihedrals[1], (2, 3, 4, 1))
  assert_eq(topo.impropers.length(), 1)
  assert_eq(topo.impropers[0], (1, 2, 3, 4))
  assert_eq(topo.donors.length(), 1)
  assert_eq(topo.donors[0], (1, 2))
  assert_eq(topo.acceptors.length(), 1)
  assert_eq(topo.acceptors[0], (3, 4))
  assert_true(@moonbit.assert_num_close(topo.atoms[0].charge, -0.4157, 1.0e-9))
  assert_true(@moonbit.assert_num_close(topo.atoms[1].charge, 0.2719, 1.0e-9))
  assert_true(@moonbit.assert_num_close(topo.atoms[2].mass, 12.011, 1.0e-9))
  assert_true(@moonbit.assert_num_close(topo.atoms[3].mass, 16.0, 1.0e-9))
}

///|
test "psf_parse_missing_sections" {
  // Only a NATOM section: all connectivity tables and titles stay empty.
  let text = "PSF\n\n       1 !NATOM\n   1 SYS     1  GLY  CA   CT2    0.100000       12.0110           0\n"
  let topo = @src.psf_parse(text)
  assert_eq(topo.natoms(), 1)
  assert_eq(topo.atoms[0].segid, "SYS")
  assert_eq(topo.atoms[0].resname, "GLY")
  assert_eq(topo.n_bonds(), 0)
  assert_eq(topo.n_angles(), 0)
  assert_eq(topo.n_dihedrals(), 0)
  assert_eq(topo.n_impropers(), 0)
  assert_eq(topo.donors.length(), 0)
  assert_eq(topo.acceptors.length(), 0)
  assert_eq(topo.titles.length(), 0)
}

///|
test "psf_write_roundtrip" {
  let atoms : Array[@src.PsfAtom] = [
    @src.PsfAtom::{
      atomid: 1,
      segid: "A",
      resid: 1,
      resname: "ALA",
      atom_name: "N",
      atom_type: "NH1",
      charge: -0.4157,
      mass: 14.01,
    },
    @src.PsfAtom::{
      atomid: 2,
      segid: "A",
      resid: 1,
      resname: "ALA",
      atom_name: "CA",
      atom_type: "CT1",
      charge: 0.0213,
      mass: 12.011,
    },
    @src.PsfAtom::{
      atomid: 3,
      segid: "B",
      resid: 2,
      resname: "HOH",
      atom_name: "O",
      atom_type: "OH2",
      charge: -0.834,
      mass: 16.0,
    },
  ]
  let topo = @src.PsfTopology::new(
    atoms,
    [(1, 2)],
    [(1, 2, 3)],
    [],
    [],
    [(1, 2)],
    [],
    ["REMARKS round trip"],
  )
  let text = @src.psf_write(topo)
  let back = @src.psf_parse(text)
  assert_eq(back.natoms(), 3)
  assert_eq(back.n_bonds(), 1)
  assert_eq(back.n_angles(), 1)
  assert_eq(back.n_dihedrals(), 0)
  assert_eq(back.n_impropers(), 0)
  assert_eq(back.bonds[0], (1, 2))
  assert_eq(back.angles[0], (1, 2, 3))
  assert_eq(back.donors[0], (1, 2))
  assert_eq(back.acceptors.length(), 0)
  assert_eq(back.titles[0], "REMARKS round trip")
  assert_eq(back.atoms[1].segid, "A")
  assert_eq(back.atoms[1].resid, 1)
  assert_eq(back.atoms[1].atom_name, "CA")
  assert_eq(back.atoms[1].atom_type, "CT1")
  assert_eq(back.atoms[2].resname, "HOH")
  assert_true(@moonbit.assert_num_close(back.atoms[0].charge, -0.4157, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back.atoms[0].mass, 14.01, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back.atoms[1].charge, 0.0213, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back.atoms[2].charge, -0.834, 1.0e-9))
}

///|
test "psf_natoms_total_charge" {
  let text = "PSF EXT\n\n       4 !NATOM\n   1 A     1  ALA  N    NH1   -0.415700       14.0100           0\n   2 A     1  ALA  H    H      0.271900        1.0080           0\n   3 A     1  ALA  CA   CT1    0.021300       12.0110           0\n   4 B     2  HOH  O    OH2   -0.834000       16.0000           0\n"
  let topo = @src.psf_parse(text)
  assert_eq(topo.natoms(), 4)
  // -0.4157 + 0.2719 + 0.0213 - 0.834 = -0.9565
  assert_true(@moonbit.assert_num_close(topo.total_charge(), -0.9565, 1.0e-9))
  assert_eq(@src.psf_is_protein_resname("ALA"), true)
  assert_eq(@src.psf_is_protein_resname("TRP"), true)
  assert_eq(@src.psf_is_protein_resname("HOH"), false)
  assert_eq(@src.psf_is_protein_resname("DA"), false)
}

///|
test "psf_parse_bad_atom_row_raises" {
  // A NATOM entry with only 4 of the 9 columns must raise.
  let text = "PSF EXT\n\n       1 !NATOM\n   1 A     1  ALA\n"
  let raised = try {
    let _ = @src.psf_parse(text)
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "psf_parse_unknown_section_skipped" {
  // An unknown keyword section skips exactly its declared entry count.
  let text = "PSF EXT\n\n       1 !NATOM\n   1 SYS     1  GLY  CA   CT2    0.100000       12.0110           0\n\n       2 !NUNKNOWN: extra\n   1     2\n   3     4\n\n       1 !NBOND: bonds\n   1     1\n"
  let topo = @src.psf_parse(text)
  assert_eq(topo.natoms(), 1)
  assert_eq(topo.n_bonds(), 1)
  assert_eq(topo.bonds[0], (1, 1))
}
