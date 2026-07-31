///|
/// Test file for RaggedExperiment module.

test "mutation_record_creation" {
  let mr = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  )
  assert_eq(mr.sample_id, "S1")
  assert_eq(mr.gene_symbol, "TP53")
  assert_eq(mr.chrom, "chr17")
  assert_eq(mr.pos, 7577121)
  assert_eq(mr.mutation_type, @src.mutation_missense())
}

test "mutation_record_mutation_id" {
  let mr = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  )
  let id = mr.mutation_id()
  assert_true(id.contains("TP53"))
  assert_true(id.contains("chr17"))
}

test "mutation_record_is_los" {
  let missense = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
    mutation_type=@src.mutation_missense(),
  )
  assert_false(missense.is_loss_of_function())

  let nonsense = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
    mutation_type=@src.mutation_nonsense(),
  )
  assert_true(nonsense.is_loss_of_function())

  let fs = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
    mutation_type=@src.mutation_fs_del(),
  )
  assert_true(fs.is_loss_of_function())
}

test "mutation_type_to_string" {
  assert_eq(@src.mutation_missense().to_string(), "Missense_Mutation")
  assert_eq(@src.mutation_nonsense().to_string(), "Nonsense_Mutation")
  assert_eq(@src.mutation_fs_del().to_string(), "Frame_Shift_Del")
  assert_eq(@src.mutation_fs_ins().to_string(), "Frame_Shift_Ins")
  assert_eq(@src.mutation_splice_site().to_string(), "Splice_Site")
  assert_eq(@src.mutation_silent().to_string(), "Silent_Mutation")
}

test "ragged_experiment_empty" {
  let exp = @src.RaggedExperiment::new()
  assert_eq(exp.n_rows(), 0)
  assert_eq(exp.n_cols(), 0)
}

test "ragged_experiment_add_record" {
  let exp = @src.RaggedExperiment::new()
  let mr = @src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  )
  exp.add_record(mr)
  assert_eq(exp.n_rows(), 1)
  assert_eq(exp.n_cols(), 1)
}

test "ragged_experiment_multiple_genes" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="KRAS",
    chrom="chr12", pos=25398284,
    ref_allele="G", alt_allele="T",
  ))
  assert_eq(exp.n_rows(), 2)
  assert_eq(exp.n_cols(), 1)
}

test "ragged_experiment_multiple_samples" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S2", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  assert_eq(exp.n_rows(), 1)
  assert_eq(exp.n_cols(), 2)
}

test "ragged_experiment_get_records" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7578000,
    ref_allele="G", alt_allele="A",
  ))

  let recs = exp.get_records(gene="TP53", sample="S1")
  assert_eq(recs.length(), 2)
}

test "ragged_experiment_get_gene_records" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S2", gene_symbol="TP53",
    chrom="chr17", pos=7578000,
    ref_allele="G", alt_allele="A",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="KRAS",
    chrom="chr12", pos=25398284,
    ref_allele="G", alt_allele="T",
  ))

  let tp53_recs = exp.get_gene_records(gene="TP53")
  assert_eq(tp53_recs.length(), 2)
}

test "ragged_experiment_get_sample_records" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="KRAS",
    chrom="chr12", pos=25398284,
    ref_allele="G", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S2", gene_symbol="MYC",
    chrom="chr8", pos=128748315,
    ref_allele="C", alt_allele="A",
  ))

  let s1_recs = exp.get_sample_records(sample="S1")
  assert_eq(s1_recs.length(), 2)

  let s2_recs = exp.get_sample_records(sample="S2")
  assert_eq(s2_recs.length(), 1)
}

test "ragged_experiment_tmb" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="KRAS",
    chrom="chr12", pos=25398284,
    ref_allele="G", alt_allele="T",
  ))
  exp.add_record(@src.MutationRecord::new(
    sample_id="S2", gene_symbol="TP53",
    chrom="chr17", pos=7578000,
    ref_allele="G", alt_allele="A",
  ))

  let tmb = exp.get_tmb()
  assert_eq(tmb.length(), 2)
  assert_eq(tmb[0], 2.0)
  assert_eq(tmb[1], 1.0)
}

test "ragged_experiment_tmb_per_mb" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))

  let tmb_mb = exp.get_tmb_per_mb(genome_size_mb=3000.0)
  assert_true(tmb_mb.length() > 0)
  assert_true(tmb_mb[0] > 0.0)
}

test "ragged_experiment_summary" {
  let exp = @src.RaggedExperiment::new()
  exp.add_record(@src.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
  ))
  let summary = exp.summary()
  assert_true(summary.length() > 0)
}

test "ragged_experiment_filter_by_type" {
  let exp = @src.ragged_sample_data()
  let types : Array[@src.MutationType] = Array::new()
  types.push(@src.mutation_nonsense())
  let filtered = exp.filter_by_type(types=types)
  assert_true(filtered.n_rows() > 0)
}

test "ragged_experiment_filter_by_genes" {
  let exp = @src.ragged_sample_data()
  let genes : Array[String] = Array::new()
  genes.push("TP53")
  let filtered = exp.filter_by_genes(genes=genes)
  assert_eq(filtered.n_rows(), 1)
}

test "ragged_experiment_filter_by_samples" {
  let exp = @src.ragged_sample_data()
  let samples : Array[String] = Array::new()
  samples.push("Sample1")
  let filtered = exp.filter_by_samples(samples=samples)
  assert_eq(filtered.n_cols(), 1)
}

test "ragged_experiment_genes_mutated_per_sample" {
  let exp = @src.ragged_sample_data()
  let mutated = exp.genes_mutated_per_sample()
  assert_true(mutated.length() > 0)
}

test "ragged_experiment_count_matrix" {
  let exp = @src.ragged_sample_data()
  let matrix = exp.get_count_matrix()
  assert_true(matrix.length() > 0)
}

test "ragged_sample_data" {
  let exp = @src.ragged_sample_data()
  assert_true(exp.n_rows() > 0)
  assert_true(exp.n_cols() > 0)
  let tmb = exp.get_tmb()
  assert_true(tmb.length() > 0)
}
