///|
/// Tests for Proteomics module - proteolytic cleavage, peptide mass, and isotope patterns.

///|
test "proteomics_trypsin_cleavage_sites" {
  let sites = @src.proteomics_cleavage_sites("PEPTIDEKTEST", "trypsin")
  assert_eq(sites.length(), 1)
  assert_eq(sites[0], 8)
}

///|
test "proteomics_trypsin_cleavage_sites_with_proline_exception" {
  let sites = @src.proteomics_cleavage_sites("PEPTKPKTEST", "trypsin")
  assert_eq(sites.length(), 1)
  assert_eq(sites[0], 7)
}

///|
test "proteomics_trypsin_cleavage_multiple_sites" {
  let sites = @src.proteomics_cleavage_sites("KRKR", "trypsin")
  assert_eq(sites.length(), 3)
}

///|
test "proteomics_chymotrypsin_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPTFYWLTEST", "chymotrypsin")
  assert_eq(sites.length(), 4)
}

///|
test "proteomics_chymotrypsin_proline_exception" {
  let sites = @src.proteomics_cleavage_sites("PEPTFPTEST", "chymotrypsin")
  assert_eq(sites.length(), 0)
}

///|
test "proteomics_pepsin_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPTFLTEST", "pepsin_ph2")
  assert_eq(sites.length(), 2)
}

///|
test "proteomics_lysc_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPTKTEST", "lysc")
  assert_eq(sites.length(), 1)
}

///|
test "proteomics_lysc_proline_exception" {
  let sites = @src.proteomics_cleavage_sites("PEPTKPTEST", "lysc")
  assert_eq(sites.length(), 0)
}

///|
test "proteomics_argc_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPRTEST", "argc")
  assert_eq(sites.length(), 1)
}

///|
test "proteomics_cnbr_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPMTEST", "cnbr")
  assert_eq(sites.length(), 1)
}

///|
test "proteomics_gluc_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPDETEST", "gluc")
  assert_eq(sites.length(), 4)
}

///|
test "proteomics_aspn_cleavage" {
  let sites = @src.proteomics_cleavage_sites("PEPDTEST", "aspn")
  assert_eq(sites.length(), 1)
  assert_eq(sites[0], 3)
}

///|
test "proteomics_digest_trypsin" {
  let result = @src.proteomics_digest("PEPTIDEKTEST", "trypsin")
  assert_eq(result.num_peptides, 2)
  assert_eq(result.peptides[0].peptide_sequence, "PEPTIDEK")
  assert_eq(result.peptides[1].peptide_sequence, "TEST")
}

///|
test "proteomics_digest_chymotrypsin" {
  let result = @src.proteomics_digest("PEPTFWLTEST", "chymotrypsin")
  assert_eq(result.num_peptides, 4)
}

///|
test "proteomics_calculate_mass_monoisotopic" {
  let mass = @src.proteomics_calculate_mass("A", true)
  assert_true(mass > 0.0)
  let expected = 71.03711 + 18.010565
  assert_true(abs(mass - expected) < 0.01)
}

///|
test "proteomics_calculate_mass_average" {
  let mass = @src.proteomics_calculate_mass("A", false)
  assert_true(mass > 0.0)
  let expected = 89.0935 + 18.01524
  assert_true(abs(mass - expected) < 0.01)
}

///|
test "proteomics_molecular_weight" {
  let mw = @src.proteomics_molecular_weight("A")
  let mass = @src.proteomics_calculate_mass("A", false)
  assert_true(abs(mw - mass) < 0.001)
}

///|
test "proteomics_mass_empty_peptide" {
  let mass = @src.proteomics_calculate_mass("", true)
  assert_eq(mass, 0.0)
  let mass2 = @src.proteomics_calculate_mass("", false)
  assert_eq(mass2, 0.0)
}

///|
test "proteomics_isotope_pattern" {
  let pattern = @src.proteomics_isotope_pattern("A", 1)
  assert_true(pattern.masses.length() >= 3)
  assert_true(pattern.abundances.length() >= 3)
  assert_eq(pattern.charge, 1)
  assert_true(pattern.abundances[0] == 1.0)
}

///|
test "proteomics_isotope_pattern_charge" {
  let pattern = @src.proteomics_isotope_pattern("GGR", 2)
  assert_eq(pattern.charge, 2)
  assert_true(pattern.masses.length() >= 3)
}

///|
test "proteomics_fragment_ions" {
  let frags = @src.proteomics_fragment_ions("PEPT", 1)
  assert_true(frags.length() >= 6)
}

///|
test "proteomics_fragment_ions_charge" {
  let frags = @src.proteomics_fragment_ions("PEPT", 2)
  assert_true(frags.length() >= 6)
}

///|
test "proteomics_select_enzyme" {
  let enzyme = @src.proteomics_select_enzyme("trypsin")
  assert_true(enzyme is Some(_))
  assert_eq(enzyme.unwrap().name, "trypsin")
}

///|
test "proteomics_select_enzyme_invalid" {
  let enzyme = @src.proteomics_select_enzyme("invalid_enzyme")
  assert_true(enzyme is None)
}

///|
test "proteomics_available_enzymes" {
  let enzymes = @src.proteomics_available_enzymes()
  assert_true(enzymes.length() >= 8)
}

///|
test "proteomics_trypsin_digest_no_missed" {
  let result = @src.proteomics_trypsin_digest("PEPTIDEKTESTKR", 0)
  assert_true(result.num_peptides >= 3)
}

///|
test "proteomics_trypsin_digest_missed_cleavages" {
  let result = @src.proteomics_trypsin_digest("PEPTIDEKTESTKR", 1)
  assert_true(result.num_peptides > 3)
}

///|
test "proteomics_digest_empty_sequence" {
  let result = @src.proteomics_digest("", "trypsin")
  assert_eq(result.num_peptides, 0)
  assert_eq(result.protein_mass, 0.0)
}

///|
test "proteomics_digest_single_residue" {
  let result = @src.proteomics_digest("A", "trypsin")
  assert_eq(result.num_peptides, 1)
  assert_eq(result.peptides[0].peptide_sequence, "A")
}

///|
test "proteomics_aspn_cleavage_before_d" {
  let sites = @src.proteomics_cleavage_sites("PEPDTEST", "aspn")
  assert_eq(sites.length(), 1)
  assert_eq(sites[0], 3)
}

///|
test "proteomics_fragment_ions_empty" {
  let frags = @src.proteomics_fragment_ions("", 1)
  assert_eq(frags.length(), 0)
}

///|
test "proteomics_fragment_ions_single_residue" {
  let frags = @src.proteomics_fragment_ions("A", 1)
  assert_eq(frags.length(), 0)
}

///|
fn abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}