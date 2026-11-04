///|
/// BSGenome Demo
fn main {
  println("=== BSGenome Demo ===")

  println("\n1. Create Genome Object:")
  println("========================")
  let genome = @bio.create_demo_genome()
  println("Name: " + genome.name)
  println("Organism: " + genome.organism)
  println("Provider: " + genome.provider)
  println("Release: " + genome.release)
  println("Chromosome Count: " + genome.chromosome_count().to_string())
  println("Total Length: " + genome.total_length().to_string())

  println("\n2. Chromosome Names:")
  println("====================")
  let chr_names = genome.chromosome_names()
  let mut i = 0
  while i < chr_names.length() {
    let len = genome.chromosome_length(chr_names[i])
    println(
      "  - " +
      chr_names[i] +
      " (" +
      (if len is Some(_) { len.unwrap().to_string() } else { "N/A" }) +
      " bp)",
    )
    i = i + 1
  }

  println("\n3. Get Chromosome Sequence:")
  println("===========================")
  let chr1_seq = genome.chromosome_sequence("chr1")
  if chr1_seq is Some(_) {
    println("chr1: " + chr1_seq.unwrap()[0:20].to_owned() + "...")
  }

  println("\n4. Get Subsequence:")
  println("===================")
  let subseq = genome.get_sequence("chr1", 0, 10)
  if subseq is Some(_) {
    println("chr1:0-10: " + subseq.unwrap())
  }

  let subseq2 = genome.get_subseq("chr1", 10, 8)
  if subseq2 is Some(_) {
    println("chr1:10+8: " + subseq2.unwrap())
  }

  println("\n5. Extract Gene (Strand-aware):")
  println("===============================")
  let gene_plus = genome.extract_gene("chr1", 0, 12, "+")
  if gene_plus is Some(_) {
    println("Gene (+ strand): " + gene_plus.unwrap())
  }

  let gene_minus = genome.extract_gene("chr1", 0, 12, "-")
  if gene_minus is Some(_) {
    println("Gene (- strand): " + gene_minus.unwrap())
  }

  println("\n6. Check Chromosome Existence:")
  println("==============================")
  println("Has chr1: " + genome.has_chromosome("chr1").to_string())
  println("Has chr4: " + genome.has_chromosome("chr4").to_string())
  println("Has chrX: " + genome.has_chromosome("chrX").to_string())

  println("\n=== Demo Complete ===")
}
