///|
fn xml_ff_water_string() -> String {
  "<ForceField name=\"water.xml\">\n" +
  "  <AtomTypes>\n" +
  "    <Type name=\"OW\" element=\"O\" mass=\"15.9994\"/>\n" +
  "    <Type name=\"HW\" element=\"H\" mass=\"1.008\"/>\n" +
  "  </AtomTypes>\n" +
  "  <Residues>\n" +
  "    <Residue name=\"HOH\">\n" +
  "      <Atom name=\"O\" type=\"OW\"/>\n" +
  "      <Atom name=\"H1\" type=\"HW\"/>\n" +
  "      <Atom name=\"H2\" type=\"HW\"/>\n" +
  "      <Bond from=\"0\" to=\"1\"/>\n" +
  "      <Bond from=\"0\" to=\"2\"/>\n" +
  "    </Residue>\n" +
  "  </Residues>\n" +
  "  <HarmonicBondForce>\n" +
  "    <Bond type1=\"OW\" type2=\"HW\" k=\"553.0\" r0=\"0.9572\"/>\n" +
  "  </HarmonicBondForce>\n" +
  "  <NonbondedForce>\n" +
  "    <Atom type=\"OW\" charge=\"-0.834\" sigma=\"0.316557\" epsilon=\"0.650194\"/>\n" +
  "  </NonbondedForce>\n" +
  "</ForceField>\n"
}

///|
test "omm_parse_xml_ff_water" {
  let ff = @src.omm_parse_xml_ff(xml_ff_water_string())
  assert_eq(ff.name, "water.xml")
  assert_eq(ff.atom_types.length(), 2)
  assert_eq(ff.atom_types[0].name, "OW")
  assert_eq(ff.atom_types[0].element, "O")
  assert_true(
    @moonbit.assert_num_close(ff.atom_types[0].mass, 15.9994, 1.0e-12),
  )
  assert_eq(ff.atom_types[1].name, "HW")
  assert_eq(ff.atom_types[1].element, "H")
  assert_true(@moonbit.assert_num_close(ff.atom_types[1].mass, 1.008, 1.0e-12))
  assert_eq(ff.residues.length(), 1)
  let res = ff.residues[0]
  assert_eq(res.name, "HOH")
  assert_eq(res.atoms.length(), 3)
  assert_eq(res.bonds.length(), 2)
  let (an, at) = res.atoms[1]
  assert_eq(an, "H1")
  assert_eq(at, "HW")
  let (b1, b2) = res.bonds[1]
  assert_eq(b1, "0")
  assert_eq(b2, "2")
  assert_eq(ff.bond_forces.length(), 1)
  let (t1, t2, k, r0) = ff.bond_forces[0]
  assert_eq(t1, "OW")
  assert_eq(t2, "HW")
  assert_true(@moonbit.assert_num_close(k, 553.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r0, 0.9572, 1.0e-12))
  assert_eq(ff.nonbonded.length(), 1)
  let (nt, nq, ns, ne) = ff.nonbonded[0]
  assert_eq(nt, "OW")
  assert_true(@moonbit.assert_num_close(nq, -0.834, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ns, 0.316557, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ne, 0.650194, 1.0e-12))
}

///|
test "omm_match_residue_template_assignment" {
  let ff = @src.omm_parse_xml_ff(xml_ff_water_string())
  let assigned = @src.omm_match_residue_template(ff, "hoh", [
    "O", "H1", "H2", "EP",
  ])
  assert_eq(assigned.length(), 4)
  let (a0n, a0t) = assigned[0]
  assert_eq(a0n, "O")
  assert_eq(a0t, "OW")
  let (a1n, a1t) = assigned[1]
  assert_eq(a1n, "H1")
  assert_eq(a1t, "HW")
  let (a3n, a3t) = assigned[3]
  assert_eq(a3n, "EP")
  assert_eq(a3t, "X")
  // unknown residue yields all-unknown types
  let miss = @src.omm_match_residue_template(ff, "XYZ", ["A", "B"])
  assert_eq(miss.length(), 2)
  let (_, mt0) = miss[0]
  assert_eq(mt0, "X")
  let (_, mt1) = miss[1]
  assert_eq(mt1, "X")
}

///|
test "omm_lookup_nonbonded" {
  let ff = @src.omm_parse_xml_ff(xml_ff_water_string())
  let (q, s, e) = @src.omm_lookup_nonbonded(ff, "OW")
  assert_true(@moonbit.assert_num_close(q, -0.834, 1.0e-12))
  assert_true(@moonbit.assert_num_close(s, 0.316557, 1.0e-12))
  assert_true(@moonbit.assert_num_close(e, 0.650194, 1.0e-12))
  // unknown type falls back to zeros
  let (q0, s0, e0) = @src.omm_lookup_nonbonded(ff, "ZZ")
  assert_true(@moonbit.assert_num_close(q0, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(s0, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(e0, 0.0, 1.0e-12))
}
