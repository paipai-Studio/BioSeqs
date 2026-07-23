///|
/// maftools Demo - Cancer Genomics MAF Analysis.
/// Inspired by Bioconductor maftools package.

///|
fn main {
  println("=== maftools Cancer Genomics Demo ===")
  println("Inspired by Bioconductor maftools package")
  println("")
  
  println("1. Creating MAF data with cancer mutations...")
  let mut maf_data = @src.MAFData::new()
  
  // Add sample mutations
  let mut i = 0
  while i < 50 {
    let gene = if i < 20 { "TP53" } else if i < 35 { "BRCA1" } else { "KRAS" }
    let chromosome = if gene == "TP53" { "chr17" } else if gene == "BRCA1" { "chr17" } else { "chr12" }
    let start = 1000000 + i * 1000
    let mut_type = if i % 3 == 0 { "Missense_Mutation" } else if i % 3 == 1 { "Nonsense_Mutation" } else { "Silent" }
    
    let mutation = @src.MAFMutation::new(
      gene,
      chromosome,
      start,
      start,
      "C",
      "T",
      "Sample_" + ((i % 5) + 1).to_string(),
      mut_type
    )
    
    maf_data = maf_data.add_mutation(mutation)
    i = i + 1
  }
  
  println("   Total mutations: \{maf_data.count_mutations().to_string()}")
  println("   Samples: \{maf_data.count_unique_samples().to_string()}")
  println("   Genes: \{maf_data.count_unique_genes().to_string()}")
  println("")
  
  println("2. Gene mutation counts:")
  let gene_counts = @src.get_gene_mutation_counts(maf_data)
  println("   TP53: \{gene_counts.get("TP53").unwrap_or(0).to_string()} mutations")
  println("   BRCA1: \{gene_counts.get("BRCA1").unwrap_or(0).to_string()} mutations")
  println("   KRAS: \{gene_counts.get("KRAS").unwrap_or(0).to_string()} mutations")
  println("")
  
  println("3. Tumor Mutational Burden (TMB) calculation:")
  let tmb_result = @src.calculate_tmb(maf_data, 3.0e7)
  println("   Total TMB: \{tmb_result.tmb.to_string()} mutations/Mb")
  println("   Total mutations: \{tmb_result.total_mutations.to_string()}")
  println("")
  
  println("4. Mutation spectrum analysis:")
  let spectrum = @src.calculate_mutation_spectrum(maf_data)
  println("   SNVs: \{spectrum.snv_count.to_string()}")
  println("   Indels: \{spectrum.indel_count.to_string()}")
  println("   Transitions: \{spectrum.transition_count.to_string()}")
  println("   Transversions: \{spectrum.transversion_count.to_string()}")
  println("")
  
  println("5. Sample mutation counts:")
  let sample_counts = @src.get_sample_mutation_counts(maf_data)
  let samples = maf_data.sample_barcodes
  let mut k = 0
  while k < samples.length() {
    let sample = samples[k]
    let count = sample_counts.get(sample).unwrap_or(0)
    println("   \{sample}: \{count.to_string()} mutations")
    k = k + 1
  }
  println("")
  
  println("6. Mutation classification:")
  let mut snv_count = 0
  let mut indel_count = 0
  let mut j = 0
  while j < maf_data.mutations.length() {
    let mutation = maf_data.mutations[j]
    if mutation.is_snv() {
      snv_count = snv_count + 1
    } else if mutation.is_indel() {
      indel_count = indel_count + 1
    }
    j = j + 1
  }
  println("   SNVs: \{snv_count.to_string()}")
  println("   Indels: \{indel_count.to_string()}")
  println("")
  
  println("=== maftools demo completed ===")
}
