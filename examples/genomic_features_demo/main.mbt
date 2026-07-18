///|
fn main {
  println("=== GenomicFeatures Demo ===")
  
  let annotation = @src.create_example_annotation()
  
  println("\n1. Annotation Summary:")
  println(@src.genome_annotation_summarize(annotation))
  
  println("\n2. Get gene by ID (ENSG00000130203):")
  let tp53 = @src.genome_annotation_get_gene_by_id(annotation, "ENSG00000130203")
  if tp53.is_some() {
    println(tp53.unwrap().to_string())
  }
  
  println("\n3. Get genes on chromosome chr17:")
  let chr17_genes = @src.genome_annotation_genes_on_chromosome(annotation, "chr17")
  let mut i = 0
  while i < chr17_genes.length() {
    println(chr17_genes[i].gene_name + " (" + chr17_genes[i].gene_id + ")")
    i = i + 1
  }
  
  println("\n4. Get genes in region chr17:7,000,000-8,000,000:")
  let region_genes = @src.genome_annotation_genes_in_region(annotation, "chr17", 7000000, 8000000)
  i = 0
  while i < region_genes.length() {
    println(region_genes[i].gene_name + ": " + region_genes[i].chromosome + ":" + region_genes[i].start.to_string() + "-" + region_genes[i].end.to_string())
    i = i + 1
  }
  
  println("\n5. Get transcripts for TP53:")
  let tp53_transcripts = @src.genome_annotation_get_transcripts_by_gene(annotation, "ENSG00000130203")
  i = 0
  while i < tp53_transcripts.length() {
    let tx = tp53_transcripts[i]
    println("  " + tx.transcript_id + " (" + tx.transcript_name + "), exons: " + tx.exons.length().to_string())
    i = i + 1
  }
  
  println("\n6. Parse GTF content:")
  let gtf_content = 
    "chr17\tENSEMBL\tgene\t7661779\t7687550\t.\t-\t.\tgene_id \"ENSG00000130203\"; gene_name \"TP53\"; gene_biotype \"protein_coding\";\n" + "chr17\tENSEMBL\ttranscript\t7661779\t7687550\t.\t-\t.\ttranscript_id \"ENST00000357654\"; gene_id \"ENSG00000130203\";\n"
  
  let parsed = @src.gf_parse_gtf(gtf_content)
  println("  Genes parsed: " + @src.genome_annotation_count_genes(parsed).to_string())
  println("  Transcripts parsed: " + @src.genome_annotation_count_transcripts(parsed).to_string())
  
  println("\nGenomicFeatures demo completed!")
}