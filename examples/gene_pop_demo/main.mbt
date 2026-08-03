///|
/// gene_pop module demo - GenePop file format parsing and analysis.
///
/// Demonstrates parsing of GenePop population genetics files, including
/// genotype extraction, allele frequency calculation, and heterozygosity
/// statistics.
///
/// Reference: Biopython Bio.PopGen.GenePop
/// Run: moon run examples/gene_pop_demo
fn main {
  println("=== gene_pop Module Demo ===")
  println("GenePop file format parsing and analysis")
  println("Reference: Biopython Bio.PopGen.GenePop")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample GenePop file
  // ------------------------------------------------------------------
  println("=== 1. Parse Sample GenePop File ===")
  let sample = @bio.genepop_sample()
  let record = @bio.parse_genepop(sample)
  println("Title: " + record.title)
  println("Loci digit size: " + record.locus_digit_size.to_string())
  println("Loci names: " + record.loci_names.length().to_string())
  let mut li = 0
  while li < record.loci_names.length() {
    println("  Locus " + (li + 1).to_string() + ": " + record.loci_names[li])
    li = li + 1
  }
  println("Populations: " + @bio.genepop_num_populations(record).to_string())
  println("Total individuals: " + @bio.genepop_total_individuals(record).to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Population and individual details
  // ------------------------------------------------------------------
  println("=== 2. Population and Individual Details ===")
  let mut pi = 0
  while pi < record.populations.length() {
    let pop = record.populations[pi]
    println("Population: " + pop.name + "  individuals: " + pop.individuals.length().to_string())
    let mut ii = 0
    while ii < pop.individuals.length() {
      let ind = pop.individuals[ii]
      let mut geno_str = ""
      let mut gi = 0
      while gi < ind.genotypes.length() {
        let g = ind.genotypes[gi]
        if gi > 0 {
          geno_str = geno_str + " "
        }
        match g.allele2 {
          Some(a2) => geno_str = geno_str + g.allele1.to_string() + "/" + a2.to_string()
          None => geno_str = geno_str + g.allele1.to_string()
        }
        gi = gi + 1
      }
      println("  " + ind.name + ": " + geno_str)
      ii = ii + 1
    }
    pi = pi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Allele frequencies
  // ------------------------------------------------------------------
  println("=== 3. Allele Frequencies (Locus 1) ===")
  let mut ppi = 0
  while ppi < record.populations.length() {
    let pop = record.populations[ppi]
    let freqs = @bio.genepop_allele_frequencies(pop, 0)
    println("Population " + pop.name + ":")
    let keys = freqs.keys()
    let key_arr : Array[Int] = []
    let mut ki = 0
    for k in keys {
      key_arr.push(k)
      ki = ki + 1
    }
    // Simple bubble sort for deterministic ordering.
    let mut a = 0
    while a < key_arr.length() {
      let mut b = a + 1
      while b < key_arr.length() {
        if key_arr[b] < key_arr[a] {
          let tmp = key_arr[a]
          key_arr[a] = key_arr[b]
          key_arr[b] = tmp
        }
        b = b + 1
      }
      a = a + 1
    }
    let mut si = 0
    while si < key_arr.length() {
      let allele = key_arr[si]
      let freq = freqs.get_or_default(allele, 0.0)
      println("  Allele " + allele.to_string() + ": " + freq.to_string())
      si = si + 1
    }
    ppi = ppi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Heterozygosity statistics
  // ------------------------------------------------------------------
  println("=== 4. Heterozygosity Statistics (Locus 1) ===")
  let mut hi = 0
  while hi < record.populations.length() {
    let pop = record.populations[hi]
    let obs = @bio.genepop_observed_heterozygosity(pop, 0)
    let exp = @bio.genepop_expected_heterozygosity(pop, 0)
    println(
      "Population " + pop.name +
      ": observed=" + obs.to_string() +
      "  expected=" + exp.to_string(),
    )
    hi = hi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Round-trip serialization
  // ------------------------------------------------------------------
  println("=== 5. Round-trip Serialization ===")
  let serialized = @bio.genepop_to_string(record)
  let reparsed = @bio.parse_genepop(serialized)
  println("Re-parsed populations: " + @bio.genepop_num_populations(reparsed).to_string())
  println("Re-parsed individuals: " + @bio.genepop_total_individuals(reparsed).to_string())
  println("Digit size preserved: " + reparsed.locus_digit_size.to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. Summary report
  // ------------------------------------------------------------------
  println("=== 6. Summary Report ===")
  println(@bio.genepop_summary(record))
  println("")

  println("=== Demo Complete ===")
}
