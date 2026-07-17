///|
/// GenomeInfoDb Demo - Genome information management.
fn main {
  println("=== GenomeInfoDb Demo ===\n")

  println("1. Loading pre-built hg38 genome...")
  let hg38 = @src.genome_hg38()
  println("   \{@src.GenomeInfoDb::summary(hg38)}")

  println("\n2. Getting chromosome information...")
  let names = @src.GenomeInfoDb::seqnames(hg38)
  println("   Total chromosomes: \{names.length()}")
  println("   First 5: \{names[0:5].to_array().join(", ")}")
  match @src.GenomeInfoDb::seqLength(hg38, "chr1") {
    Some(len) => println("   chr1 length: \{len} bp")
    None => println("   chr1 not found")
  }

  println("\n3. Checking circular chromosomes...")
  match @src.GenomeInfoDb::is_circular(hg38, "chrM") {
    Some(c) => println("   chrM is circular: \{c}")
    None => println("   chrM not found")
  }
  match @src.GenomeInfoDb::is_circular(hg38, "chr1") {
    Some(c) => println("   chr1 is circular: \{c}")
    None => println("   chr1 not found")
  }

  println("\n4. Getting centromere positions...")
  match @src.GenomeInfoDb::centromere(hg38, "chr1") {
    Some(c) => println("   chr1 centromere: \{c.start}-\{c.end}")
    None => println("   chr1 centromere not found")
  }

  println("\n5. Getting chromosomal arms...")
  let arm1 = @src.GenomeInfoDb::get_arm(hg38, "chr1", 1000000)
  let arm2 = @src.GenomeInfoDb::get_arm(hg38, "chr1", 200000000)
  let arm3 = @src.GenomeInfoDb::get_arm(hg38, "chr1", 123000000)
  println("   chr1:1,000,000 -> arm \{arm1}")
  println("   chr1:200,000,000 -> arm \{arm2}")
  println("   chr1:123,000,000 -> arm \{arm3}")

  println("\n6. Comparing genome builds...")
  let hg19 = @src.genome_hg19()
  let mm10 = @src.genome_mm10()
  match @src.GenomeInfoDb::seqLength(hg38, "chr1") {
    Some(l38) =>
      match @src.GenomeInfoDb::seqLength(hg19, "chr1") {
        Some(l19) => println("   hg38 chr1: \{l38} bp, hg19 chr1: \{l19} bp")
        None => println("   hg19 chr1 not found")
      }
    None => println("   hg38 chr1 not found")
  }
  match @src.GenomeInfoDb::seqLength(mm10, "chr1") {
    Some(l) => println("   mm10 chr1: \{l} bp")
    None => println("   mm10 chr1 not found")
  }

  println("\n7. Looking up organism from build...")
  println("   hg38 -> \{@src.organism_from_build("hg38")}")
  println("   mm10 -> \{@src.organism_from_build("mm10")}")
  println("   dm6 -> \{@src.organism_from_build("dm6")}")

  println("\n8. Filtering chromosomes...")
  let chr1_like = @src.GenomeInfoDb::filter_chromosomes(hg38, "chr1")
  println("   Chromosomes matching 'chr1': \{chr1_like.join(", ")}")

  println("\n9. Building custom genome...")
  let custom = @src.GenomeInfoDb::new(
      genome_build="custom",
      organism="Test organism",
    )
    .add_chrom("chr1", 5000000)
    .add_chrom("chr2", 4000000)
    .add_chrom("chrM", 16000, is_circular=true)
  println("   \{@src.GenomeInfoDb::summary(custom)}")

  println("\n10. Supported genome builds...")
  let builds = @src.supported_genome_builds()
  println("   \{builds.join(", ")}")

  println("\n=== GenomeInfoDb Demo Complete ===")
}
