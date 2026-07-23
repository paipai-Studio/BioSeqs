///|

test "variant_filtering_create_variant" {
  let info = @hashmap.HashMap([("DP", "20"), ("AF", "0.5")])
  let genotypes = [("sample1", "0/1:30"), ("sample2", "1/1:25")]
  
  let variant = @src.Variant::new(
    "chr1", 1000, "rs123", "A", "T", 50.0, "PASS", info, genotypes
  )
  
  assert_eq(variant.chr, "chr1")
  assert_eq(variant.pos, 1000)
  assert_eq(variant.id, "rs123")
  assert_eq(variant.ref_allele, "A")
  assert_eq(variant.alt_allele, "T")
  assert_eq(variant.qual, 50.0)
  assert_eq(variant.filter, "PASS")
  assert_eq(variant.genotypes.length(), 2)
}

test "variant_filtering_filter" {
  let info_pass = @hashmap.HashMap([("DP", "20"), ("AF", "0.5")])
  let info_fail = @hashmap.HashMap([("DP", "5"), ("AF", "0.001")])
  
  let variants = [
    @src.Variant::new("chr1", 1000, "rs1", "A", "T", 50.0, "PASS", info_pass, [("s1", "0/1:30")]),
    @src.Variant::new("chr1", 2000, "rs2", "G", "C", 20.0, "FAIL", info_fail, [("s1", "0/1:10")]),
    @src.Variant::new("chr1", 3000, "rs3", "T", "A", 60.0, "PASS", info_pass, [("s1", "0/1:35")]),
  ]
  
  let params = @src.VariantFilteringParam::new()
  let result = @src.bio_vfilter_filter(variants, params)
  
  assert_eq(result.passed_variants.length(), 2)
  assert_eq(result.filtered_variants.length(), 1)
}

test "variant_filtering_check_autosomal_dominant" {
  let genotypes = [
    ("proband", "1/1"),
    ("father", "0/1"),
    ("mother", "0/0"),
    ("sibling", "0/0"),
  ]
  
  let info = @hashmap.HashMap([])
  let variant = @src.Variant::new("chr1", 1000, "rs1", "A", "T", 50.0, "PASS", info, genotypes)
  
  let phenotypes = @hashmap.HashMap([("proband", true), ("father", false), ("mother", false), ("sibling", false)])
  let model = @src.GeneticModel::new("autosomal_dominant")
  
  let passes = @src.bio_vfilter_check_inheritance(variant, model, phenotypes)
  
  assert_true(passes)
}

test "variant_filtering_check_autosomal_recessive" {
  let genotypes = [
    ("proband", "1/1"),
    ("father", "0/1"),
    ("mother", "0/1"),
    ("sibling", "0/0"),
  ]
  
  let info = @hashmap.HashMap([])
  let variant = @src.Variant::new("chr1", 1000, "rs1", "A", "T", 50.0, "PASS", info, genotypes)
  
  let phenotypes = @hashmap.HashMap([("proband", true), ("father", false), ("mother", false), ("sibling", false)])
  let model = @src.GeneticModel::new("autosomal_recessive")
  
  let passes = @src.bio_vfilter_check_inheritance(variant, model, phenotypes)
  
  assert_true(passes)
}

test "variant_filtering_summary" {
  let info_pass = @hashmap.HashMap([("DP", "20"), ("AF", "0.5")])
  let info_fail = @hashmap.HashMap([("DP", "5"), ("AF", "0.001")])
  
  let variants = [
    @src.Variant::new("chr1", 1000, "rs1", "A", "T", 50.0, "PASS", info_pass, [("s1", "0/1:30")]),
    @src.Variant::new("chr1", 2000, "rs2", "G", "C", 20.0, "FAIL", info_fail, [("s1", "0/1:10")]),
  ]
  
  let params = @src.VariantFilteringParam::new()
  let result = @src.bio_vfilter_filter(variants, params)
  
  let summary = @src.bio_vfilter_summary(result)
  
  assert_true(summary.contains("Total variants"))
  assert_true(summary.contains("Passed"))
  assert_true(summary.contains("Filtered"))
}

test "variant_filtering_predict_consequence" {
  let info = @hashmap.HashMap([])
  let genotypes = [("s1", "0/1")]
  
  let snp = @src.Variant::new("chr1", 1000, "rs1", "A", "T", 50.0, "PASS", info, genotypes)
  let insertion = @src.Variant::new("chr1", 2000, "rs2", "A", "AT", 50.0, "PASS", info, genotypes)
  let deletion = @src.Variant::new("chr1", 3000, "rs3", "AT", "A", 50.0, "PASS", info, genotypes)
  
  assert_eq(@src.bio_vfilter_predict_consequence(snp), "missense_variant")
  assert_eq(@src.bio_vfilter_predict_consequence(insertion), "insertion")
  assert_eq(@src.bio_vfilter_predict_consequence(deletion), "deletion")
}