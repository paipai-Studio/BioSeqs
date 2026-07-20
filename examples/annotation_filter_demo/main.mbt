///|
fn main {
  let annotations = @bio.af_create_example_annotations()

  println("Total annotations: " + annotations.length().to_string())

  let chr17_genes = @bio.af_filter_chromosome(annotations, "chr17")
  println("Genes on chr17: " + chr17_genes.length().to_string())

  let mirna_genes = @bio.af_filter_biotype(annotations, "miRNA")
  println("miRNA genes: " + mirna_genes.length().to_string())

  let minus_strand = @bio.af_filter_strand(annotations, "-")
  println("Genes on - strand: " + minus_strand.length().to_string())

  let cancer_genes = @bio.af_filter_description_contains(annotations, "cancer")
  println("Cancer-related genes: " + cancer_genes.length().to_string())

  let brca_genes = @bio.af_filter_symbol_pattern(annotations, "BRCA")
  println("BRCA genes: " + brca_genes.length().to_string())

  let summary = @bio.af_filter_summary(annotations)
  println(summary)
}
