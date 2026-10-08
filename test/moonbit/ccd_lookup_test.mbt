test "ccd_lookup_amino_acid" {
  let ala = @src.ccd_lookup("ALA")
  assert_true(ala is Some(_))
  match ala {
    Some(e) => {
      assert_eq(e.id, "ALA")
      assert_eq(e.name, "ALANINE")
      assert_eq(e.one_letter_code, "A")
      assert_true(@src.ccd_is_amino_acid("ALA"))
    }
    None => fail("expected ALA entry")
  }
}

test "ccd_one_letter_code" {
  assert_eq(@src.ccd_one_letter("ALA"), "A")
  assert_eq(@src.ccd_one_letter("ARG"), "R")
  assert_eq(@src.ccd_one_letter("ASN"), "N")
  assert_eq(@src.ccd_one_letter("ASP"), "D")
  assert_eq(@src.ccd_one_letter("CYS"), "C")
  assert_eq(@src.ccd_one_letter("GLN"), "Q")
  assert_eq(@src.ccd_one_letter("GLU"), "E")
  assert_eq(@src.ccd_one_letter("GLY"), "G")
  assert_eq(@src.ccd_one_letter("HIS"), "H")
  assert_eq(@src.ccd_one_letter("ILE"), "I")
  assert_eq(@src.ccd_one_letter("LEU"), "L")
  assert_eq(@src.ccd_one_letter("LYS"), "K")
  assert_eq(@src.ccd_one_letter("MET"), "M")
  assert_eq(@src.ccd_one_letter("PHE"), "F")
  assert_eq(@src.ccd_one_letter("PRO"), "P")
  assert_eq(@src.ccd_one_letter("SER"), "S")
  assert_eq(@src.ccd_one_letter("THR"), "T")
  assert_eq(@src.ccd_one_letter("TRP"), "W")
  assert_eq(@src.ccd_one_letter("TYR"), "Y")
  assert_eq(@src.ccd_one_letter("VAL"), "V")
  assert_eq(@src.ccd_one_letter("XXX"), "X")
}

test "ccd_is_amino_acid" {
  assert_true(@src.ccd_is_amino_acid("ALA"))
  assert_true(@src.ccd_is_amino_acid("GLY"))
  assert_true(@src.ccd_is_amino_acid("TRP"))
  assert_true(not(@src.ccd_is_amino_acid("HOH")))
  assert_true(not(@src.ccd_is_amino_acid("ATP")))
  assert_true(not(@src.ccd_is_amino_acid("XYZ")))
}

test "ccd_is_nucleotide" {
  assert_true(@src.ccd_is_nucleotide("DA"))
  assert_true(@src.ccd_is_nucleotide("DT"))
  assert_true(@src.ccd_is_nucleotide("A"))
  assert_true(@src.ccd_is_nucleotide("U"))
  assert_true(not(@src.ccd_is_nucleotide("ALA")))
  assert_true(not(@src.ccd_is_nucleotide("HOH")))
}

test "ccd_is_water" {
  assert_true(@src.ccd_is_water("HOH"))
  assert_true(@src.ccd_is_water("WAT"))
  assert_true(@src.ccd_is_water("H2O"))
  assert_true(not(@src.ccd_is_water("ALA")))
}

test "ccd_atom_element" {
  // Standard backbone atoms
  assert_eq(@src.ccd_atom_element("ALA", "N"), "N")
  assert_eq(@src.ccd_atom_element("ALA", "CA"), "C")
  assert_eq(@src.ccd_atom_element("ALA", "C"), "C")
  assert_eq(@src.ccd_atom_element("ALA", "O"), "O")
  assert_eq(@src.ccd_atom_element("ALA", "CB"), "C")
  assert_eq(@src.ccd_atom_element("CYS", "SG"), "S")
}

test "ccd_guess_element" {
  // Direct lookup for known residue
  assert_eq(@src.ccd_atom_element("ALA", "N"), "N")
  // Unknown residue → guess from atom name
  assert_eq(@src.ccd_atom_element("XXX", "CA"), "C")
  assert_eq(@src.ccd_atom_element("XXX", "FE"), "FE")
  assert_eq(@src.ccd_atom_element("XXX", "ZN"), "ZN")
  assert_eq(@src.ccd_atom_element("XXX", "CL"), "CL")
}

test "ccd_lookup_unknown" {
  assert_true(@src.ccd_lookup("XYZ") is None)
}

test "ccd_lookup_nucleotide" {
  let da = @src.ccd_lookup("DA")
  assert_true(da is Some(_))
  match da {
    Some(e) => assert_eq(e.id, "DA")
    None => fail("expected DA entry")
  }
}

test "ccd_lookup_ion" {
  let ca = @src.ccd_lookup("CA")
  assert_true(ca is Some(_))
  match ca {
    Some(e) => {
      assert_eq(e.name, "CALCIUM ION")
      assert_true(not(@src.ccd_is_amino_acid("CA")))
    }
    None => fail("expected CA entry")
  }
}
