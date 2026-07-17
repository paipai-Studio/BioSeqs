///|
test "restriction_builtin_enzymes" {
  let enzymes = @src.builtin_enzymes()

  assert_true(enzymes.keys().collect().length() > 0)
}

///|
test "restriction_get_enzyme_by_name" {
  let enzyme = @src.get_enzyme_by_name("EcoRI")

  assert_eq(enzyme.is_some(), true)
}

///|
test "restriction_find_cut_sites" {
  let enzyme = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let sites = @src.find_cut_sites("GAATTCGAATTC", enzyme)

  assert_eq(sites.length(), 2)
}

///|
test "restriction_digest" {
  let enzyme = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let fragments = @src.digest("GAATTCGAATTC", enzyme)

  assert_eq(fragments.length(), 3)
}

///|
test "restriction_digest_with_multiple" {
  let enzyme1 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let enzyme2 = @src.restriction_enzyme("HindIII", "AAGCTT", 1, "AGCT", "")
  let fragments = @src.digest_with_multiple("GAATTCAAGCTT", [enzyme1, enzyme2])

  assert_eq(fragments.length(), 3)
}

///|
test "restriction_fragment_sizes" {
  let fragments = ["ACGT", "ACAT"]
  let sizes = @src.fragment_sizes(fragments)

  assert_eq(sizes[0], 4)
}

///|
test "restriction_is_compatible_overhang" {
  let enzyme1 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let enzyme2 = @src.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")

  assert_true(@src.is_compatible_overhang(enzyme1, enzyme2))
}

///|
test "restriction_list_enzymes" {
  let list = @src.list_enzymes()

  assert_true(list.length() > 0)
}

///|
test "restriction_enzymes_by_overhang" {
  let enzymes = @src.enzymes_by_overhang("AATT")

  assert_true(enzymes.length() > 0)
}

///|
test "restriction_enzymes_by_site_length" {
  let enzymes = @src.enzymes_by_site_length(6, 6)

  assert_true(enzymes.length() > 0)
}
