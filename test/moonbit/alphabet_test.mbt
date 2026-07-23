// Tests for Bio.Alphabet module

test "IUPAC unambiguous DNA alphabet" {
  let dna = @src.iupac_unambiguous_dna()
  assert_true(dna.name == "IUPACUnambiguousDNA")
  assert_true(dna.is_gapped == false)
  assert_true(dna.letters.length() == 4)
  assert_true(dna.is_valid("A"))
  assert_true(dna.is_valid("C"))
  assert_true(dna.is_valid("G"))
  assert_true(dna.is_valid("T"))
  assert_false(dna.is_valid("U"))
  assert_false(dna.is_valid("X"))
}

test "IUPAC unambiguous RNA alphabet" {
  let rna = @src.iupac_unambiguous_rna()
  assert_true(rna.name == "IUPACUnambiguousRNA")
  assert_true(rna.is_gapped == false)
  assert_true(rna.letters.length() == 4)
  assert_true(rna.is_valid("A"))
  assert_true(rna.is_valid("C"))
  assert_true(rna.is_valid("G"))
  assert_true(rna.is_valid("U"))
  assert_false(rna.is_valid("T"))
}

test "IUPAC ambiguous DNA alphabet" {
  let dna = @src.iupac_ambiguous_dna()
  assert_true(dna.name == "IUPACAmbiguousDNA")
  assert_true(dna.letters.length() == 15)
  assert_true(dna.is_valid("R"))
  assert_true(dna.is_valid("Y"))
  assert_true(dna.is_valid("N"))
}

test "IUPAC protein alphabet" {
  let protein = @src.iupac_protein()
  assert_true(protein.name == "IUPACProtein")
  assert_true(protein.letters.length() == 23)
  assert_true(protein.is_valid("A"))
  assert_true(protein.is_valid("W"))
  assert_true(protein.is_valid("X"))
}

test "Gapped DNA alphabet" {
  let dna = @src.gapped_dna()
  assert_true(dna.is_gapped)
  assert_true(dna.is_valid("-"))
  assert_true(dna.is_valid("A"))
}

test "Reduced protein alphabet" {
  let protein = @src.reduced_protein()
  assert_true(protein.letters.length() == 20)
  assert_true(protein.is_valid("A"))
  assert_false(protein.is_valid("X"))
}

test "Three to one letter code" {
  let mapping = @src.three_to_one()
  let ala = mapping.get("Ala")
  assert_true(ala is Some("A"))
  let gly = mapping.get("Gly")
  assert_true(gly is Some("G"))
  let ter = mapping.get("Ter")
  assert_true(ter is Some("*"))
}

test "One to three letter code" {
  let mapping = @src.one_to_three()
  let a = mapping.get("A")
  assert_true(a is Some("Ala"))
  let ter = mapping.get("*")
  assert_true(ter is Some("Ter"))
}
