///|
/// ensembldb demo.

pub fn main() {
  let db = create_example_ensdb()
  println("EnsDb created with")
  println("- version: " + db.version)
  println("- species: " + db.species)
  println("- genes: " + edb_num_genes(db).to_string())
  println("- transcripts: " + edb_num_transcripts(db).to_string())
  println("- exons: " + edb_num_exons(db).to_string())
  
  let gene = edb_get_gene_by_id(db, "ENSG000001")
  println("\nGene ENSG000001:")
  if gene.is_some() {
    let g = gene.unwrap()
    println("- name: " + g.gene_name)
    println("- biotype: " + g.biotype)
    println("- chromosome: " + g.chromosome)
  }
  
  let genes_by_name = edb_get_gene_by_name(db, "ACTB")
  println("\nGenes named ACTB: " + genes_by_name.length().to_string())
  
  let txs = edb_get_transcripts_by_gene(db, "ENSG000001")
  println("Transcripts for ENSG000001: " + txs.length().to_string())
  
  let chr17 = edb_filter_by_chromosome(db, "17")
  println("\nGenes on chromosome 17: " + edb_num_genes(chr17).to_string())
  
  let gene_length = edb_get_gene_length(db, "ENSG000001")
  println("Gene length: " + gene_length.to_string())
}