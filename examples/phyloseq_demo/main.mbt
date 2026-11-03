///|
/// phyloseq demo.

pub fn main() {
  let ps = create_example_phyloseq()
  println("phyloseq created with")
  println("- OTUs: " + ps_num_otus(ps).to_string())
  println("- samples: " + ps_num_samples(ps).to_string())
  
  let total = ps_total_abundance(ps)
  println("\nTotal abundance: " + total.to_string())
  
  let filtered = ps_filter_by_abundance(ps, 300.0)
  println("Filtered OTUs (>=300): " + ps_num_otus(filtered).to_string())
  
  let proteo = ps_filter_by_taxonomy(ps, "phylum", "Proteobacteria")
  println("Proteobacteria OTUs: " + ps_num_otus(proteo).to_string())
  
  let summary = ps_get_taxa_summary(ps, "phylum")
  println("\nPhylum summary:")
  if summary.contains("Proteobacteria") {
    println("- Proteobacteria: " + summary.get("Proteobacteria").unwrap().to_string())
  }
}