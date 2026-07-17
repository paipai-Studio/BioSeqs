///|
fn main {
  println("=== Ballgown (Transcript-level DE Analysis) Demo ===")

  let bg = @bio.sample_ballgown()

  println("")
  println("Ballgown Object Summary:")
  println("  Samples: " + bg.sample_names.length().to_string())
  println("  Transcripts: " + bg.transcripts.length().to_string())
  println("  Exons: " + bg.exons.length().to_string())
  println("  Introns: " + bg.introns.length().to_string())

  let gene_ids = bg.get_gene_ids()
  println("  Genes: " + gene_ids.length().to_string())

  println("")
  println("Sample names:")
  let mut i = 0
  while i < bg.sample_names.length() {
    println("  " + (i + 1).to_string() + ". " + bg.sample_names[i])
    i = i + 1
  }

  println("")
  println("Gene list:")
  let gene_names = bg.get_gene_names()
  i = 0
  while i < gene_names.length() {
    let gene_id = gene_ids[i]
    let txs = bg.get_gene_transcripts(gene_id)
    println(
      "  " +
      gene_names[i] +
      " (" +
      gene_id +
      "): " +
      txs.length().to_string() +
      " transcripts",
    )
    i = i + 1
  }

  println("")
  println("FPKM matrix (first 5 transcripts):")
  let fpkm_mat = bg.get_fpkm_matrix()
  i = 0
  while i < 5 && i < bg.transcripts.length() {
    let t = bg.transcripts[i]
    let mut j = 0
    let mut line = t.id + ": "
    while j < t.fpkm.length() {
      line = line + t.fpkm[j].to_string() + " "
      j = j + 1
    }
    println("  " + line)
    i = i + 1
  }

  println("")
  println("=== Differential Expression Analysis ===")
  println("Group 0: samples 0-2 (control)")
  println("Group 1: samples 3-5 (treatment)")
  let group = [0, 0, 0, 1, 1, 1]

  println("")
  println("Transcript-level DE results:")
  let de_tx = @bio.stattest(bg, group, feature="transcript")
  i = 0
  while i < de_tx.length() {
    let r = de_tx[i]
    println(
      "  " +
      r.feature_id +
      " | log2FC=" +
      r.log2_fc.to_string() +
      " | pval=" +
      r.p_value.to_string() +
      " | qval=" +
      r.q_value.to_string(),
    )
    i = i + 1
  }

  println("")
  println("Gene-level DE results:")
  let de_gene = @bio.stattest(bg, group, feature="gene")
  i = 0
  while i < de_gene.length() {
    let r = de_gene[i]
    println(
      "  " +
      r.gene_name +
      " | log2FC=" +
      r.log2_fc.to_string() +
      " | pval=" +
      r.p_value.to_string() +
      " | qval=" +
      r.q_value.to_string(),
    )
    i = i + 1
  }

  println("")
  println("Gene expression for GeneA:")
  let expr_geneA = @bio.gene_expression(bg, gene_ids[0])
  i = 0
  while i < expr_geneA.length() {
    println("  " + bg.sample_names[i] + ": " + expr_geneA[i].to_string())
    i = i + 1
  }

  println("")
  println("=== Demo Complete ===")
}
