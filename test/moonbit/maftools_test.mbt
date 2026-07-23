///|
test "maftools_create_mutation" {
  let mutation = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  assert_eq(mutation.hugo_symbol, "TP53")
  assert_eq(mutation.chromosome, "chr17")
  assert_eq(mutation.start_position, 7577121)
  assert_eq(mutation.end_position, 7577121)
  assert_eq(mutation.reference_allele, "C")
  assert_eq(mutation.tumor_seq_allele1, "T")
  assert_eq(mutation.tumor_sample_barcode, "Sample1")
  assert_eq(mutation.variant_classification, "Missense_Mutation")
}

///|
test "maftools_snv_detection" {
  let snv = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  assert_true(snv.is_snv())
  assert_false(snv.is_indel())
  assert_false(snv.is_complex())
}

///|
test "maftools_indel_detection" {
  let insertion = @src.MAFMutation::new(
    "BRCA1",
    "chr17",
    43091134,
    43091135,
    "CT",
    "CTT",
    "Sample1",
    "Frame_Shift_Ins"
  )
  
  assert_false(insertion.is_snv())
  assert_true(insertion.is_indel())
  assert_false(insertion.is_complex())
}

///|
test "maftools_transition_detection" {
  let transition = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  assert_true(transition.is_transition())
  assert_false(transition.is_transversion())
}

///|
test "maftools_transversion_detection" {
  let transversion = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "A",
    "Sample1",
    "Missense_Mutation"
  )
  
  assert_false(transversion.is_transition())
  assert_true(transversion.is_transversion())
}

///|
test "maftools_maf_data_operations" {
  let mut maf = @src.MAFData::new()
  
  let mutation1 = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  let mutation2 = @src.MAFMutation::new(
    "BRCA1",
    "chr17",
    43091134,
    43091134,
    "A",
    "G",
    "Sample2",
    "Nonsense_Mutation"
  )
  
  maf = maf.add_mutation(mutation1)
  maf = maf.add_mutation(mutation2)
  
  assert_eq(maf.count_mutations(), 2)
  assert_eq(maf.count_unique_genes(), 2)
  assert_eq(maf.count_unique_samples(), 2)
}

///|
test "maftools_filter_genes" {
  let mut maf = @src.MAFData::new()
  
  let mutation1 = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  let mutation2 = @src.MAFMutation::new(
    "BRCA1",
    "chr17",
    43091134,
    43091134,
    "A",
    "G",
    "Sample2",
    "Nonsense_Mutation"
  )
  
  maf = maf.add_mutation(mutation1)
  maf = maf.add_mutation(mutation2)
  
  let filtered = maf.filter_genes(["TP53"])
  assert_eq(filtered.count_mutations(), 1)
  assert_eq(filtered.count_unique_genes(), 1)
}

///|
test "maftools_filter_samples" {
  let mut maf = @src.MAFData::new()
  
  let mutation1 = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  let mutation2 = @src.MAFMutation::new(
    "BRCA1",
    "chr17",
    43091134,
    43091134,
    "A",
    "G",
    "Sample2",
    "Nonsense_Mutation"
  )
  
  maf = maf.add_mutation(mutation1)
  maf = maf.add_mutation(mutation2)
  
  let filtered = maf.filter_samples(["Sample1"])
  assert_eq(filtered.count_mutations(), 1)
}

///|
test "maftools_mutation_spectrum" {
  let mut maf = @src.MAFData::new()
  
  let snv1 = @src.MAFMutation::new(
    "TP53",
    "chr17",
    7577121,
    7577121,
    "C",
    "T",
    "Sample1",
    "Missense_Mutation"
  )
  
  let snv2 = @src.MAFMutation::new(
    "BRCA1",
    "chr17",
    43091134,
    43091134,
    "A",
    "G",
    "Sample1",
    "Missense_Mutation"
  )
  
  let indel = @src.MAFMutation::new(
    "EGFR",
    "chr7",
    55086714,
    55086715,
    "CT",
    "CTT",
    "Sample2",
    "Frame_Shift_Ins"
  )
  
  maf = maf.add_mutation(snv1)
  maf = maf.add_mutation(snv2)
  maf = maf.add_mutation(indel)
  
  let spectrum = @src.calculate_mutation_spectrum(maf)
  
  assert_eq(spectrum.snv_count, 2)
  assert_eq(spectrum.indel_count, 1)
  assert_eq(spectrum.complex_count, 0)
}

///|
test "maftools_tmb_calculation" {
  let mut maf = @src.MAFData::new()
  
  let mut i = 0
  while i < 100 {
    let mutation = @src.MAFMutation::new(
      "TP53",
      "chr17",
      7577121 + i,
      7577121 + i,
      "C",
      "T",
      "Sample1",
      "Missense_Mutation"
    )
    maf = maf.add_mutation(mutation)
    i = i + 1
  }
  
  let tmb_result = @src.calculate_tmb(maf, 3.0e7)
  
  assert_eq(tmb_result.total_mutations, 100)
  assert_eq(tmb_result.coding_region_size, 3.0e7)
  
  // TMB should be 100 / 30000000 = 0.00000333...
  let expected_tmb = 100.0 / 3.0e7
  let diff = tmb_result.tmb - expected_tmb
  assert_true(diff.abs() < 1.0e-12)
}

///|
test "maftools_co_occurrence_analysis" {
  let mut maf = @src.MAFData::new()
  
  // Sample1 has both TP53 and BRCA1 mutations
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 100, 100, "C", "T", "Sample1", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("BRCA1", "chr17", 200, 200, "A", "G", "Sample1", "Missense_Mutation"))
  
  // Sample2 has only TP53 mutation
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 300, 300, "C", "T", "Sample2", "Missense_Mutation"))
  
  // Sample3 has only BRCA1 mutation
  maf = maf.add_mutation(@src.MAFMutation::new("BRCA1", "chr17", 400, 400, "A", "G", "Sample3", "Missense_Mutation"))
  
  let result = @src.analyze_co_occurrence(maf, "TP53", "BRCA1")
  
  assert_eq(result.co_occurrence, 1)
  assert_eq(result.mutual_exclusivity, 2)
}

///|
test "maftools_create_example" {
  let maf = @src.create_example_maf()
  
  assert_true(maf.count_mutations() > 0)
  assert_true(maf.count_unique_genes() > 0)
  assert_true(maf.count_unique_samples() > 0)
}

///|
test "maftools_summarize" {
  let maf = @src.create_example_maf()
  let summary = @src.summarize_maf(maf)
  
  assert_true(summary.contains("MAF Summary"))
  assert_true(summary.contains("Total mutations:"))
  assert_true(summary.contains("TMB:"))
}

///|
test "maftools_oncoplot_data" {
  let maf = @src.create_example_maf()
  let oncoplot = @src.generate_oncoplot_data(maf, 5)
  
  assert_true(oncoplot.genes.length() <= 5)
  assert_true(oncoplot.samples.length() > 0)
  assert_true(oncoplot.mutation_matrix.length() > 0)
}

///|
test "maftools_get_gene_mutation_counts" {
  let mut maf = @src.MAFData::new()
  
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 100, 100, "C", "T", "Sample1", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 200, 200, "A", "G", "Sample2", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("BRCA1", "chr17", 300, 300, "C", "T", "Sample1", "Missense_Mutation"))
  
  let counts = @src.get_gene_mutation_counts(maf)
  
  assert_eq(counts.get("TP53").unwrap_or(0), 2)
  assert_eq(counts.get("BRCA1").unwrap_or(0), 1)
}

///|
test "maftools_get_sample_mutation_counts" {
  let mut maf = @src.MAFData::new()
  
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 100, 100, "C", "T", "Sample1", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("BRCA1", "chr17", 200, 200, "A", "G", "Sample1", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("EGFR", "chr7", 300, 300, "C", "T", "Sample2", "Missense_Mutation"))
  
  let counts = @src.get_sample_mutation_counts(maf)
  
  assert_eq(counts.get("Sample1").unwrap_or(0), 2)
  assert_eq(counts.get("Sample2").unwrap_or(0), 1)
}

///|
test "maftools_filter_variant_type" {
  let mut maf = @src.MAFData::new()
  
  maf = maf.add_mutation(@src.MAFMutation::new("TP53", "chr17", 100, 100, "C", "T", "Sample1", "Missense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("BRCA1", "chr17", 200, 200, "A", "G", "Sample2", "Nonsense_Mutation"))
  maf = maf.add_mutation(@src.MAFMutation::new("EGFR", "chr7", 300, 300, "C", "T", "Sample3", "Missense_Mutation"))
  
  let filtered = maf.filter_variant_type("Nonsense")
  
  assert_eq(filtered.count_mutations(), 1)
  assert_eq(filtered.get_mutation_gene(0), "BRCA1")
}

///|
test "maftools_parse_maf_content" {
  let content = "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Sample_Barcode\tVariant_Classification\nTP53\tchr17\t7577121\t7577121\tC\tT\tSample1\tMissense_Mutation\nBRCA1\tchr17\t43091134\t43091134\tA\tG\tSample2\tNonsense_Mutation"
  
  let maf = @src.parse_maf_content(content)
  
  assert_eq(maf.count_mutations(), 2)
  assert_eq(maf.get_mutation_gene(0), "TP53")
  assert_eq(maf.get_mutation_gene(1), "BRCA1")
}
