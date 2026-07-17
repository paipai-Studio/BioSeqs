///|
fn main {
  println("=== TxDb (Transcript Database) Demo ===")
  println("Inspired by Bioconductor GenomicFeatures package")
  println("")

  let gtf_content = "chr1\tensembl\tgene\t100\t500\t.\t+\t.\tgene_id \"GeneA\"; gene_name \"GeneA\";\n" +
    "chr1\tensembl\ttranscript\t100\t500\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\";\n" +
    "chr1\tensembl\texon\t100\t150\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"1\";\n" +
    "chr1\tensembl\texon\t200\t300\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"2\";\n" +
    "chr1\tensembl\texon\t350\t500\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"3\";\n" +
    "chr1\tensembl\tCDS\t120\t150\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"1\";\n" +
    "chr1\tensembl\tCDS\t200\t300\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"2\";\n" +
    "chr1\tensembl\tCDS\t350\t420\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"3\";\n" +
    "chr1\tensembl\tgene\t1000\t1500\t.\t-\t.\tgene_id \"GeneB\"; gene_name \"GeneB\";\n" +
    "chr1\tensembl\ttranscript\t1000\t1500\t.\t-\t.\tgene_id \"GeneB\"; transcript_id \"Tx2\";\n" +
    "chr1\tensembl\texon\t1000\t1100\t.\t-\t.\tgene_id \"GeneB\"; transcript_id \"Tx2\"; exon_number \"1\";\n" +
    "chr1\tensembl\texon\t1200\t1350\t.\t-\t.\tgene_id \"GeneB\"; transcript_id \"Tx2\"; exon_number \"2\";\n" +
    "chr1\tensembl\tCDS\t1050\t1100\t.\t-\t0\tgene_id \"GeneB\"; transcript_id \"Tx2\"; exon_number \"1\";\n" +
    "chr1\tensembl\tCDS\t1200\t1320\t.\t-\t0\tgene_id \"GeneB\"; transcript_id \"Tx2\"; exon_number \"2\";\n"

  println("1. Parse GTF File")
  let txdb = @bio.parse_gtf(gtf_content)
  println("GTF parsed successfully")
  println("")

  println("2. TxDb Summary")
  println(@bio.txdb_summary(txdb))
  println("")

  println("3. Extract Genes")
  let gene_gr = @bio.genes(txdb)
  let mut i = 0
  while i < @bio.granges_length(gene_gr) {
    let strand_char = @bio.strand_to_char(gene_gr.strands[i])
    println(
      "  Gene \{i + 1}: \{gene_gr.seqnames[i]}:\{gene_gr.starts[i]}-\{gene_gr.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("4. Extract Transcripts")
  let tx_gr = @bio.transcripts(txdb)
  i = 0
  while i < @bio.granges_length(tx_gr) {
    let strand_char = @bio.strand_to_char(tx_gr.strands[i])
    println(
      "  Tx \{i + 1}: \{tx_gr.seqnames[i]}:\{tx_gr.starts[i]}-\{tx_gr.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("5. Extract Exons")
  let exon_gr = @bio.exons(txdb)
  i = 0
  while i < @bio.granges_length(exon_gr) {
    let strand_char = @bio.strand_to_char(exon_gr.strands[i])
    println(
      "  Exon \{i + 1}: \{exon_gr.seqnames[i]}:\{exon_gr.starts[i]}-\{exon_gr.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("6. Extract CDS Regions")
  let cds_gr = @bio.cds(txdb)
  i = 0
  while i < @bio.granges_length(cds_gr) {
    let strand_char = @bio.strand_to_char(cds_gr.strands[i])
    println(
      "  CDS \{i + 1}: \{cds_gr.seqnames[i]}:\{cds_gr.starts[i]}-\{cds_gr.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("7. Extract Promoters (200bp upstream, 50bp downstream)")
  let prom_gr = @bio.promoters(txdb, 200, 50)
  i = 0
  while i < @bio.granges_length(prom_gr) {
    let strand_char = @bio.strand_to_char(prom_gr.strands[i])
    println(
      "  Promoter \{i + 1}: \{prom_gr.seqnames[i]}:\{prom_gr.starts[i]}-\{prom_gr.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("8. Introns by Transcript")
  let introns_map = @bio.introns_by_transcript(txdb)
  let tx_ids = @bio.transcript_ids(txdb)
  i = 0
  while i < tx_ids.length() {
    let tx_id = tx_ids[i]
    let introns_opt = introns_map.get(tx_id)
    if introns_opt.is_some() {
      let introns = introns_opt.unwrap()
      println("  \{tx_id}:")
      let mut j = 0
      while j < introns.length() {
        println("    Intron \{j + 1}: \{introns[j].0}-\{introns[j].1}")
        j = j + 1
      }
    }
    i = i + 1
  }
  println("")

  println("9. 5' UTRs")
  let five_prime_utr = @bio.five_prime_utrs(txdb)
  i = 0
  while i < @bio.granges_length(five_prime_utr) {
    let strand_char = @bio.strand_to_char(five_prime_utr.strands[i])
    println(
      "  5'UTR \{i + 1}: \{five_prime_utr.seqnames[i]}:\{five_prime_utr.starts[i]}-\{five_prime_utr.ends[i]}(\{strand_char}) w=\{five_prime_utr.widths[i]}",
    )
    i = i + 1
  }
  println("")

  println("10. 3' UTRs")
  let three_prime_utr = @bio.three_prime_utrs(txdb)
  i = 0
  while i < @bio.granges_length(three_prime_utr) {
    let strand_char = @bio.strand_to_char(three_prime_utr.strands[i])
    println(
      "  3'UTR \{i + 1}: \{three_prime_utr.seqnames[i]}:\{three_prime_utr.starts[i]}-\{three_prime_utr.ends[i]}(\{strand_char}) w=\{three_prime_utr.widths[i]}",
    )
    i = i + 1
  }
  println("")

  println("11. ID Lists")
  println("Gene IDs: \{@bio.gene_ids(txdb)}")
  println("Transcript IDs: \{@bio.transcript_ids(txdb)}")
  println("Exon IDs: \{@bio.exon_ids(txdb)}")
}
