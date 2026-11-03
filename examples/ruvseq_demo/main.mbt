///|
/// RUVSeq Demo
fn main {
  println("=== RUVSeq Demo ===")

  println("\n1. Create Demo Data:")
  println("====================")
  let data = @bio.create_demo_data()
  println("Genes: " + data.genes.length().to_string())
  println("Samples: " + data.samples.length().to_string())
  println("Batches: " + data.batch.length().to_string())

  let mut i = 0
  while i < data.samples.length() {
    println("  " + data.samples[i] + " (" + data.batch[i] + ")")
    i = i + 1
  }

  println("\n2. Raw Counts (GeneA):")
  println("======================")
  let gene_a_idx = data.get_gene_index("GeneA")
  if gene_a_idx is Some(_) {
    let counts = data.counts[gene_a_idx.unwrap()]
    let mut line = ""
    i = 0
    while i < counts.length() {
      line = line + counts[i].to_string() + "\t"
      i = i + 1
    }
    println(line)
  }

  println("\n3. Normalize Data:")
  println("==================")
  let normalized_data = data.normalize()
  println("Normalized: " + normalized_data.normalized.to_string())

  if gene_a_idx is Some(_) {
    let counts = normalized_data.counts[gene_a_idx.unwrap()]
    let mut line = ""
    i = 0
    while i < counts.length() {
      let s = counts[i].to_string()
      let len = if s.length() < 6 { s.length() } else { 6 }
      line = line + s[0:len].to_owned() + "\t"
      i = i + 1
    }
    println(line)
  }

  println("\n4. Log2 Transform:")
  println("==================")
  let log_data = data.log2_transform()

  if gene_a_idx is Some(_) {
    let counts = log_data.counts[gene_a_idx.unwrap()]
    let mut line = ""
    i = 0
    while i < counts.length() {
      let s = counts[i].to_string()
      let len = if s.length() < 5 { s.length() } else { 5 }
      line = line + s[0:len].to_owned() + "\t"
      i = i + 1
    }
    println(line)
  }

  println("\n5. RUVg (Using Control Genes):")
  println("===============================")
  let controls = ["Control1", "Control2", "Control3"]
  let ruvg_factors = @bio.ruvg(data, controls, 1)
  println("Number of factors: " + ruvg_factors.k.to_string())
  println("Factor values:")
  let factors = ruvg_factors.factors[0]
  let mut line = ""
  i = 0
  while i < factors.length() {
    let s = factors[i].to_string()
    let len = if s.length() < 8 { s.length() } else { 8 }
    line = line + s[0:len].to_owned() + "\t"
    i = i + 1
  }
  println(line)

  println("\n6. RUVs (Using All Genes):")
  println("==========================")
  let ruvs_factors = @bio.ruvs(data, 1)
  println("Number of factors: " + ruvs_factors.k.to_string())
  println("Factor values:")
  let factors2 = ruvs_factors.factors[0]
  let mut line2 = ""
  i = 0
  while i < factors2.length() {
    let s = factors2[i].to_string()
    let len = if s.length() < 8 { s.length() } else { 8 }
    line2 = line2 + s[0:len].to_owned() + "\t"
    i = i + 1
  }
  println(line2)

  println("\n7. Remove Batch Effect:")
  println("=======================")
  let corrected = @bio.remove_batch_effect(data, ruvg_factors)
  println("Corrected data genes: " + corrected.genes.length().to_string())
  println("Corrected data samples: " + corrected.samples.length().to_string())

  let corr_gene_a_idx = corrected.get_gene_index("GeneA")
  if corr_gene_a_idx is Some(_) {
    let corr_counts = corrected.counts[corr_gene_a_idx.unwrap()]
    let mut line3 = ""
    i = 0
    while i < corr_counts.length() {
      let s = corr_counts[i].to_string()
      let len = if s.length() < 6 { s.length() } else { 6 }
      line3 = line3 + s[0:len].to_owned() + "\t"
      i = i + 1
    }
    println(line3)
  }

  println("\n=== Demo Complete ===")
}
