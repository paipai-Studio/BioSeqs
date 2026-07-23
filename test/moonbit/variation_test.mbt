///|
test "variation_analyze_snp" {
  let snp = @src.analyze_snp(100, "A", "T", 50.0, "PASS")
  assert_eq(snp.position, 100)
  assert_eq(snp.ref_allele, "A")
  assert_eq(snp.alt_allele, "T")
}

///|
test "variation_detect_mutations" {
  let dna = "ATGTCTGGCTTA"
  let mutations = @src.detect_mutations("BRCA1", "BRCA1_HUMAN", dna, 4, "A")
  assert_eq(mutations.length(), 1)
}

///|
test "variation_analyze_amino_acid_substitution" {
  let sub = @src.analyze_amino_acid_substitution("G", "A")
  assert_eq(sub.original, "G")
  assert_eq(sub.mutated, "A")
}

///|
test "variation_predict_damaging" {
  let sub = @src.analyze_amino_acid_substitution("G", "W")
  assert_true(@src.predict_damaging(sub))
}

///|
test "variation_parse_vcf_line" {
  let line = "17\t7577538\t.\tG\tA\t.\tPASS\t."
  let result = @src.parse_vcf_line(line)
  assert_true(result.is_some())
}

///|
test "variation_create_example_data" {
  let records = @src.create_example_variation_data()
  assert_eq(records.length(), 2)
}