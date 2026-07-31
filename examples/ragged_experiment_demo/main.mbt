///|
/// RaggedExperiment module example - Ragged mutation data arrays.
///
/// Run: moon run examples/ragged_experiment_demo
fn main {
  println("=== RaggedExperiment Module Demo ===")
  println("Ragged Mutation Data Arrays for Multi-Sample Analysis")
  println("")

  // Create experiment from records
  println("=== Creating RaggedExperiment ===")
  let exp = @bio.ragged_sample_data()
  println("Dimensions: " + exp.n_rows().to_string() + " genes x " + exp.n_cols().to_string() + " samples")
  println("Row names: " + exp.get_rownames().to_string())
  println("Col names: " + exp.get_colnames().to_string())
  println("")

  // Summary
  println("=== Summary ===")
  println(exp.summary())
  println("")

  // Gene-level queries
  println("=== Gene-Level Queries ===")
  let tp53_recs = exp.get_gene_records(gene="TP53")
  println("TP53 mutations: " + tp53_recs.length().to_string())
  let mut i = 0
  while i < tp53_recs.length() {
    let rec = tp53_recs[i]
    println("  " + rec.sample_id + ": " + rec.mutation_type.to_string() + " (VAF=" + rec.variant_allele_frequency.to_string() + ")")
    i = i + 1
  }
  println("")

  // Sample-level queries
  println("=== Sample-Level Queries ===")
  let s1_recs = exp.get_sample_records(sample="Sample1")
  println("Sample1 mutations: " + s1_recs.length().to_string())
  let mut j = 0
  while j < s1_recs.length() {
    let rec = s1_recs[j]
    println("  " + rec.gene_symbol + ": " + rec.mutation_type.to_string())
    j = j + 1
  }
  println("")

  // TMB calculation
  println("=== Tumor Mutational Burden (TMB) ===")
  let tmb = exp.get_tmb()
  let mut k = 0
  while k < tmb.length() {
    println("  " + exp.get_colnames()[k] + ": TMB=" + tmb[k].to_string())
    k = k + 1
  }
  println("")

  let tmb_mb = exp.get_tmb_per_mb(genome_size_mb=3000.0)
  println("TMB per Mb (genome_size=3000Mb):")
  let mut k2 = 0
  while k2 < tmb_mb.length() {
    println("  " + exp.get_colnames()[k2] + ": " + tmb_mb[k2].to_string() + " mutations/Mb")
    k2 = k2 + 1
  }
  println("")

  // Genes mutated per sample
  println("=== Genes Mutated Per Sample ===")
  let mutated = exp.genes_mutated_per_sample()
  let mut m = 0
  while m < mutated.length() {
    println("  " + exp.get_colnames()[m] + ": " + mutated[m].to_string() + " genes mutated")
    m = m + 1
  }
  println("")

  // Filtering
  println("=== Filtering by Mutation Type ===")
  let types : Array[@bio.MutationType] = Array::new()
  types.push(@bio.mutation_missense())
  types.push(@bio.mutation_nonsense())
  let filtered_type = exp.filter_by_type(types=types)
  println("After filtering (missense+nonsense): " + filtered_type.n_rows().to_string() + " genes, " + filtered_type.n_cols().to_string() + " samples")
  println("")

  println("=== Filtering by Genes ===")
  let genes : Array[String] = Array::new()
  genes.push("TP53")
  genes.push("KRAS")
  let filtered_genes = exp.filter_by_genes(genes=genes)
  println("After filtering (TP53, KRAS): " + filtered_genes.n_rows().to_string() + " genes, " + filtered_genes.n_cols().to_string() + " samples")
  println("")

  println("=== Filtering by Samples ===")
  let samples : Array[String] = Array::new()
  samples.push("Sample1")
  samples.push("Sample2")
  let filtered_samples = exp.filter_by_samples(samples=samples)
  println("After filtering (Sample1, Sample2): " + filtered_samples.n_rows().to_string() + " genes, " + filtered_samples.n_cols().to_string() + " samples")
  println("")

  // Mutation type helpers
  println("=== Mutation Types ===")
  let missense = @bio.mutation_missense()
  let nonsense = @bio.mutation_nonsense()
  let fs_ins = @bio.mutation_fs_ins()
  let fs_del = @bio.mutation_fs_del()
  let splice = @bio.mutation_splice_site()
  let silent = @bio.mutation_silent()
  println("Missense: " + missense.to_string())
  println("Nonsense: " + nonsense.to_string())
  println("Frameshift Ins: " + fs_ins.to_string())
  println("Frameshift Del: " + fs_del.to_string())
  println("Splice Site: " + splice.to_string())
  println("Silent: " + silent.to_string())
  println("")

  // Loss-of-function detection
  println("=== Loss-of-Function Mutations ===")
  let records : Array[@bio.MutationRecord] = Array::new()
  records.push(@bio.MutationRecord::new(
    sample_id="S1", gene_symbol="TP53",
    chrom="chr17", pos=7577121,
    ref_allele="C", alt_allele="T",
    mutation_type=@bio.mutation_missense(),
  ))
  records.push(@bio.MutationRecord::new(
    sample_id="S1", gene_symbol="BRCA1",
    chrom="chr17", pos=43045629,
    ref_allele="C", alt_allele="T",
    mutation_type=@bio.mutation_nonsense(),
  ))
  records.push(@bio.MutationRecord::new(
    sample_id="S1", gene_symbol="KRAS",
    chrom="chr12", pos=25398284,
    ref_allele="G", alt_allele="T",
    mutation_type=@bio.mutation_fs_del(),
  ))

  let mut lo = 0
  while lo < records.length() {
    let r = records[lo]
    println("  " + r.gene_symbol + " (" + r.mutation_type.to_string() + "): LoF=" + r.is_loss_of_function().to_string())
    lo = lo + 1
  }
  println("")

  // Count matrix
  println("=== Count Matrix ===")
  let matrix = exp.get_count_matrix()
  let mut row = 0
  while row < matrix.length() {
    let mut line = "  " + exp.get_rownames()[row] + ":"
    let mut col = 0
    while col < matrix[row].length() {
      line = line + " " + matrix[row][col].to_string()
      col = col + 1
    }
    println(line)
    row = row + 1
  }
}
