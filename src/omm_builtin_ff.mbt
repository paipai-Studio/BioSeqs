/// omm_builtin_ff.mbt — built-in amber99sb-inspired force field subset.
///
/// Item 6: a small teaching-grade protein force field bundled with the
/// library so that `omm_xml_ff_build_system` can parameterize the 20
/// standard amino acids without external files.
///
/// Scope: 20 residue templates (neutral backbone form; no terminal
/// variants), a shared atom-type table, and amber99sb-style parameters in
/// OpenMM XML units (kJ/mol, nm, degrees). Partial charges are all zero —
/// per-residue charge sets must come from a custom force field. Parameter
/// values are representative amber99sb numbers for the covered type
/// combinations; uncovered combinations fall back to the generic values
/// documented in `omm_xml_ff_build_system`.

///|
/// One `<Atom name=".." type=".."/>` line.
fn omm_ff_xatom(name : String, ty : String) -> String {
  "<Atom name=\"" + name + "\" type=\"" + ty + "\"/>"
}

///|
/// One `<Bond from="i" to="j"/>` line (template atom indices).
fn omm_ff_xbond(a : Int32, b : Int32) -> String {
  "<Bond from=\"" + a.to_string() + "\" to=\"" + b.to_string() + "\"/>"
}

///|
/// Assemble one residue block from (name, type) atom pairs and index pairs.
fn omm_ff_xres(
  name : String,
  atoms : Array[(String, String)],
  bonds : Array[(Int32, Int32)],
) -> String {
  let sb : StringBuilder = StringBuilder()
  sb.write_string("<Residue name=\"" + name + "\">\n")
  for (an, ty) in atoms {
    sb.write_string(omm_ff_xatom(an, ty) + "\n")
  }
  for (a, b) in bonds {
    sb.write_string(omm_ff_xbond(a, b) + "\n")
  }
  sb.write_string("</Residue>\n")
  sb.to_string()
}

///|
/// Standard-residue backbone: N, H, CA, HA, then the side chain, then C, O.
/// Returns (atoms, backbone bonds, index of the first side-chain atom).
fn omm_ff_std_backbone(
  side : Array[(String, String)],
) -> (Array[(String, String)], Array[(Int32, Int32)], Int32) {
  let atoms : Array[(String, String)] = Array(capacity=0)
  atoms.push(("N", "N"))
  atoms.push(("H", "H"))
  atoms.push(("CA", "CT"))
  atoms.push(("HA", "HC"))
  let side0 = atoms.length()
  for a in side {
    atoms.push(a)
  }
  let c_idx = atoms.length()
  atoms.push(("C", "C"))
  atoms.push(("O", "O"))
  let bonds : Array[(Int32, Int32)] = Array(capacity=0)
  bonds.push((0, 1)) // N-H
  bonds.push((0, 2)) // N-CA
  bonds.push((2, 3)) // CA-HA
  bonds.push((2, c_idx)) // CA-C
  bonds.push((c_idx, c_idx + 1)) // C-O
  bonds.push((2, side0)) // CA-CB (first side-chain atom)
  (atoms, bonds, side0)
}

///|
/// The built-in XML force field text (amber99sb teaching subset).
pub fn omm_builtin_amber99sb_xml() -> String {
  let sb : StringBuilder = StringBuilder()
  sb.write_string("<ForceField name=\"amber99sb-builtin\">\n")
  sb.write_string("<AtomTypes>\n")
  let types : Array[(String, String, Double)] = [
    ("C", "C", 12.011),
    ("CT", "C", 12.011),
    ("CA", "C", 12.011),
    ("N", "N", 14.007),
    ("N3", "N", 14.007),
    ("O", "O", 15.999),
    ("O2", "O", 15.999),
    ("OG", "O", 15.999),
    ("S", "S", 32.06),
    ("SH", "S", 32.06),
    ("H", "H", 1.008),
    ("HC", "H", 1.008),
    ("HO", "H", 1.008),
    ("HS", "H", 1.008),
  ]
  for (ty, el, mass) in types {
    sb.write_string(
      "<Type name=\"" +
      ty +
      "\" element=\"" +
      el +
      "\" mass=\"" +
      mass.to_string() +
      "\"/>\n",
    )
  }
  sb.write_string("</AtomTypes>\n")
  sb.write_string("<Residues>\n")
  // --- residue templates ---
  // GLY (two HA instead of a side chain)
  let gly_atoms : Array[(String, String)] = [
    ("N", "N"),
    ("H", "H"),
    ("CA", "CT"),
    ("HA1", "HC"),
    ("HA2", "HC"),
    ("C", "C"),
    ("O", "O"),
  ]
  let gly_bonds : Array[(Int32, Int32)] = [
    (0, 1),
    (0, 2),
    (2, 3),
    (2, 4),
    (2, 5),
    (5, 6),
  ]
  sb.write_string(omm_ff_xres("GLY", gly_atoms, gly_bonds))
  // ALA
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("HB3", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  sb.write_string(omm_ff_xres("ALA", a, b))
  // SER
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("OG", "OG"),
    ("HG", "HO"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  sb.write_string(omm_ff_xres("SER", a, b))
  // CYS
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("SG", "SH"),
    ("HG", "HS"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  sb.write_string(omm_ff_xres("CYS", a, b))
  // THR
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB", "HC"),
    ("OG1", "OG"),
    ("HG1", "HO"),
    ("CG2", "CT"),
    ("HG21", "HC"),
    ("HG22", "HC"),
    ("HG23", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((6, 7))
  b.push((4, 8))
  b.push((8, 9))
  b.push((8, 10))
  b.push((8, 11))
  sb.write_string(omm_ff_xres("THR", a, b))
  // VAL
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB", "HC"),
    ("CG1", "CT"),
    ("HG11", "HC"),
    ("HG12", "HC"),
    ("HG13", "HC"),
    ("CG2", "CT"),
    ("HG21", "HC"),
    ("HG22", "HC"),
    ("HG23", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((6, 7))
  b.push((6, 8))
  b.push((6, 9))
  b.push((4, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((10, 13))
  sb.write_string(omm_ff_xres("VAL", a, b))
  // LEU
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG", "HC"),
    ("CD1", "CT"),
    ("HD11", "HC"),
    ("HD12", "HC"),
    ("HD13", "HC"),
    ("CD2", "CT"),
    ("HD21", "HC"),
    ("HD22", "HC"),
    ("HD23", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((9, 10))
  b.push((9, 11))
  b.push((9, 12))
  b.push((7, 13))
  b.push((13, 14))
  b.push((13, 15))
  b.push((13, 16))
  sb.write_string(omm_ff_xres("LEU", a, b))
  // ILE
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB", "HC"),
    ("CG1", "CT"),
    ("HG11", "HC"),
    ("HG12", "HC"),
    ("HG13", "HC"),
    ("CG2", "CT"),
    ("HG21", "HC"),
    ("HG22", "HC"),
    ("HG23", "HC"),
    ("CD1", "CT"),
    ("HD11", "HC"),
    ("HD12", "HC"),
    ("HD13", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((6, 7))
  b.push((6, 8))
  b.push((6, 9))
  b.push((4, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((10, 13))
  b.push((6, 14))
  b.push((14, 15))
  b.push((14, 16))
  b.push((14, 17))
  sb.write_string(omm_ff_xres("ILE", a, b))
  // PRO (no H on N; ring closure CD-N)
  let pro_atoms : Array[(String, String)] = [
    ("N", "N"),
    ("CA", "CT"),
    ("HA", "HC"),
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("CD", "CT"),
    ("HD1", "HC"),
    ("HD2", "HC"),
    ("C", "C"),
    ("O", "O"),
  ]
  let pro_bonds : Array[(Int32, Int32)] = [
    (0, 1),
    (1, 2),
    (1, 3),
    (3, 4),
    (3, 5),
    (3, 6),
    (6, 7),
    (6, 8),
    (6, 9),
    (9, 10),
    (9, 11),
    (9, 0),
    (1, 12),
    (12, 13),
  ]
  sb.write_string(omm_ff_xres("PRO", pro_atoms, pro_bonds))
  // PHE
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CA"),
    ("CD1", "CA"),
    ("HD1", "HC"),
    ("CD2", "CA"),
    ("HD2", "HC"),
    ("CE1", "CA"),
    ("HE1", "HC"),
    ("CE2", "CA"),
    ("HE2", "HC"),
    ("CZ", "CA"),
    ("HZ", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((8, 9))
  b.push((8, 12))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 14))
  b.push((12, 13))
  b.push((12, 16))
  b.push((14, 15))
  b.push((14, 16))
  b.push((16, 17))
  sb.write_string(omm_ff_xres("PHE", a, b))
  // TYR (PHE ring + OH/HH)
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CA"),
    ("CD1", "CA"),
    ("HD1", "HC"),
    ("CD2", "CA"),
    ("HD2", "HC"),
    ("CE1", "CA"),
    ("HE1", "HC"),
    ("CE2", "CA"),
    ("HE2", "HC"),
    ("CZ", "CA"),
    ("OH", "OG"),
    ("HH", "HO"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((8, 9))
  b.push((8, 12))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 14))
  b.push((12, 13))
  b.push((12, 16))
  b.push((14, 15))
  b.push((14, 16))
  b.push((16, 17))
  b.push((17, 18))
  sb.write_string(omm_ff_xres("TYR", a, b))
  // TRP
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CA"),
    ("CD1", "CA"),
    ("HD1", "HC"),
    ("NE1", "N"),
    ("HE1", "H"),
    ("CE2", "CA"),
    ("CD2", "CA"),
    ("CE3", "CA"),
    ("HE3", "HC"),
    ("CZ2", "CA"),
    ("HZ2", "HC"),
    ("CZ3", "CA"),
    ("HZ3", "HC"),
    ("CH2", "CA"),
    ("HH2", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((8, 9))
  b.push((8, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((12, 13))
  b.push((13, 7))
  b.push((13, 14))
  b.push((14, 15))
  b.push((14, 18))
  b.push((18, 19))
  b.push((18, 20))
  b.push((20, 21))
  b.push((20, 16))
  b.push((16, 17))
  b.push((16, 12))
  sb.write_string(omm_ff_xres("TRP", a, b))
  // ASP
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "C"),
    ("OD1", "O2"),
    ("OD2", "O2"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  sb.write_string(omm_ff_xres("ASP", a, b))
  // GLU
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("CD", "C"),
    ("OE1", "O2"),
    ("OE2", "O2"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 12))
  sb.write_string(omm_ff_xres("GLU", a, b))
  // ASN
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "C"),
    ("OD1", "O2"),
    ("ND2", "N"),
    ("HD21", "H"),
    ("HD22", "H"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((9, 10))
  b.push((9, 11))
  sb.write_string(omm_ff_xres("ASN", a, b))
  // GLN
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("CD", "C"),
    ("OE1", "O2"),
    ("NE2", "N"),
    ("HE21", "H"),
    ("HE22", "H"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((12, 13))
  b.push((12, 14))
  sb.write_string(omm_ff_xres("GLN", a, b))
  // HIS
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CA"),
    ("ND1", "N"),
    ("HD1", "H"),
    ("CE1", "CA"),
    ("HE1", "HC"),
    ("NE2", "N"),
    ("CD2", "CA"),
    ("HD2", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((8, 9))
  b.push((8, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((12, 13))
  b.push((13, 14))
  b.push((13, 7))
  sb.write_string(omm_ff_xres("HIS", a, b))
  // LYS
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("CD", "CT"),
    ("HD1", "HC"),
    ("HD2", "HC"),
    ("CE", "CT"),
    ("HE1", "HC"),
    ("HE2", "HC"),
    ("NZ", "N3"),
    ("HZ1", "HC"),
    ("HZ2", "HC"),
    ("HZ3", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((10, 13))
  b.push((13, 14))
  b.push((13, 15))
  b.push((13, 16))
  b.push((16, 17))
  b.push((16, 18))
  b.push((16, 19))
  sb.write_string(omm_ff_xres("LYS", a, b))
  // MET
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("SD", "S"),
    ("CE", "CT"),
    ("HE1", "HC"),
    ("HE2", "HC"),
    ("HE3", "HC"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((7, 10))
  b.push((10, 11))
  b.push((11, 12))
  b.push((11, 13))
  b.push((11, 14))
  sb.write_string(omm_ff_xres("MET", a, b))
  // ARG
  let (a, b, _) = omm_ff_std_backbone([
    ("CB", "CT"),
    ("HB1", "HC"),
    ("HB2", "HC"),
    ("CG", "CT"),
    ("HG1", "HC"),
    ("HG2", "HC"),
    ("CD", "CT"),
    ("HD1", "HC"),
    ("HD2", "HC"),
    ("NE", "N"),
    ("HE", "H"),
    ("CZ", "C"),
    ("NH1", "N"),
    ("HH11", "H"),
    ("HH12", "H"),
    ("NH2", "N"),
    ("HH21", "H"),
    ("HH22", "H"),
  ])
  b.push((4, 5))
  b.push((4, 6))
  b.push((4, 7))
  b.push((7, 8))
  b.push((7, 9))
  b.push((7, 10))
  b.push((10, 11))
  b.push((10, 12))
  b.push((10, 13))
  b.push((13, 14))
  b.push((13, 15))
  b.push((15, 16))
  b.push((16, 17))
  b.push((16, 18))
  b.push((15, 19))
  b.push((19, 20))
  b.push((19, 21))
  sb.write_string(omm_ff_xres("ARG", a, b))
  sb.write_string("</Residues>\n")
  // --- bond parameters (amber99sb; k kJ/mol/nm^2, r0 nm) ---
  sb.write_string("<HarmonicBondForce>\n")
  let bond_params : Array[(String, String, Double, Double)] = [
    ("C", "N", 1966.49, 0.1335),
    ("C", "O", 2124.96, 0.1229),
    ("C", "O2", 2124.96, 0.125),
    ("C", "CT", 1834.23, 0.1522),
    ("C", "CA", 1834.23, 0.1499),
    ("CA", "CA", 1924.56, 0.1396),
    ("CA", "CT", 1730.58, 0.151),
    ("CA", "HC", 1678.40, 0.108),
    ("CT", "HC", 1678.40, 0.109),
    ("CT", "CT", 1730.58, 0.1526),
    ("N", "CT", 1413.60, 0.147),
    ("N", "H", 3631.75, 0.101),
    ("CT", "OG", 1388.11, 0.141),
    ("OG", "HO", 1458.55, 0.096),
    ("CT", "S", 1183.94, 0.181),
    ("CT", "SH", 1183.94, 0.181),
    ("SH", "HS", 1183.94, 0.134),
    ("CT", "N3", 1413.60, 0.147),
    ("N3", "HC", 1678.40, 0.109),
  ]
  for (t1, t2, k, r0) in bond_params {
    sb.write_string(
      "<Bond type1=\"" +
      t1 +
      "\" type2=\"" +
      t2 +
      "\" k=\"" +
      k.to_string() +
      "\" r0=\"" +
      r0.to_string() +
      "\"/>\n",
    )
  }
  sb.write_string("</HarmonicBondForce>\n")
  // --- angle parameters (k kJ/mol/rad^2, theta0 deg) ---
  sb.write_string("<HarmonicAngleForce>\n")
  let angle_params : Array[(String, String, String, Double, Double)] = [
    ("C", "N", "CT", 288.70, 121.9),
    ("N", "CT", "C", 288.70, 110.1),
    ("CT", "C", "N", 422.58, 116.6),
    ("CT", "C", "O", 422.58, 120.4),
    ("O", "C", "N", 422.58, 122.9),
    ("CT", "CT", "HC", 131.80, 110.1),
    ("CT", "CT", "CT", 131.80, 111.0),
    ("CT", "CT", "OG", 131.80, 111.0),
    ("HC", "CT", "OG", 131.80, 109.5),
    ("CT", "CT", "S", 131.80, 114.7),
    ("CT", "CT", "SH", 131.80, 114.7),
    ("CA", "CA", "CA", 131.80, 120.0),
    ("CA", "CA", "HC", 131.80, 120.0),
    ("N", "CT", "HC", 131.80, 109.5),
    ("CT", "C", "O2", 422.58, 118.0),
  ]
  for (t1, t2, t3, k, th) in angle_params {
    sb.write_string(
      "<Angle type1=\"" +
      t1 +
      "\" type2=\"" +
      t2 +
      "\" type3=\"" +
      t3 +
      "\" k=\"" +
      k.to_string() +
      "\" theta0=\"" +
      th.to_string() +
      "\"/>\n",
    )
  }
  sb.write_string("</HarmonicAngleForce>\n")
  // --- torsion parameters (k kJ/mol, phase deg, periodicity) ---
  sb.write_string("<PeriodicTorsionForce>\n")
  let torsion_params : Array[
    (String, String, String, String, Double, Double, Int32),
  ] = [
    ("C", "N", "CT", "C", 0.4602, 180.0, 2),
    ("N", "CT", "C", "N", 0.4602, 0.0, 3),
    ("CT", "CT", "CT", "HC", 0.8033, 0.0, 3),
    ("CT", "CT", "CT", "CT", 1.2552, 0.0, 3),
    ("C", "CT", "CT", "CT", 0.7531, 0.0, 3),
  ]
  for (t1, t2, t3, t4, k, phi, n) in torsion_params {
    sb.write_string(
      "<Proper type1=\"" +
      t1 +
      "\" type2=\"" +
      t2 +
      "\" type3=\"" +
      t3 +
      "\" type4=\"" +
      t4 +
      "\" periodicity1=\"" +
      n.to_string() +
      "\" phase1=\"" +
      phi.to_string() +
      "\" k1=\"" +
      k.to_string() +
      "\"/>\n",
    )
  }
  sb.write_string("</PeriodicTorsionForce>\n")
  // --- nonbonded parameters (charge, sigma nm, epsilon kJ/mol) ---
  sb.write_string("<NonbondedForce>\n")
  let nb_params : Array[(String, Double, Double, Double)] = [
    ("C", 0.0, 0.33997, 0.359824),
    ("CT", 0.0, 0.33997, 0.457729),
    ("CA", 0.0, 0.33997, 0.359824),
    ("N", 0.0, 0.32500, 0.711280),
    ("N3", 0.0, 0.32500, 0.711280),
    ("O", 0.0, 0.29600, 0.878640),
    ("O2", 0.0, 0.29600, 0.878640),
    ("OG", 0.0, 0.30665, 0.711280),
    ("S", 0.0, 0.35636, 1.046000),
    ("SH", 0.0, 0.35636, 1.046000),
    ("H", 0.0, 0.10691, 0.0656888),
    ("HC", 0.0, 0.26495, 0.0656888),
    ("HO", 0.0, 0.15319, 0.0656888),
    ("HS", 0.0, 0.10691, 0.0656888),
  ]
  for (ty, q, s, e) in nb_params {
    sb.write_string(
      "<Atom type=\"" +
      ty +
      "\" charge=\"" +
      q.to_string() +
      "\" sigma=\"" +
      s.to_string() +
      "\" epsilon=\"" +
      e.to_string() +
      "\"/>\n",
    )
  }
  sb.write_string("</NonbondedForce>\n")
  sb.write_string("</ForceField>\n")
  sb.to_string()
}

///|
/// Parse the built-in amber99sb subset into an OmmXmlForceField.
pub fn omm_builtin_amber99sb() -> OmmXmlForceField {
  omm_parse_xml_ff(omm_builtin_amber99sb_xml())
}
