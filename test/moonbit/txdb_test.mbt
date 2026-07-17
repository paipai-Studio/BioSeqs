///|
/// TxDb tests

///|
test "txdb_genes" {
  let genes = [
    @src.Gene::new("GeneA", "chr1", 100, 500, @src.strand_plus()),
    @src.Gene::new("GeneB", "chr2", 200, 600, @src.strand_minus()),
  ]
  let txdb = @src.TxDb::new(genes, [], [], [], "hg19", "Homo sapiens")

  let gr = @src.genes(txdb)
  assert_eq(@src.granges_length(gr), 2)
  assert_eq(gr.seqnames[0], "chr1")
  assert_eq(gr.starts[0], 100)
  assert_eq(gr.ends[0], 500)
}

///|
test "txdb_transcripts" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 500, @src.strand_plus()),
    @src.Transcript::new("Tx2", "GeneA", "chr1", 110, 490, @src.strand_plus()),
  ]
  let txdb = @src.TxDb::new([], txs, [], [], "hg19", "Homo sapiens")

  let gr = @src.transcripts(txdb)
  assert_eq(@src.granges_length(gr), 2)
  assert_eq(gr.seqnames[1], "chr1")
}

///|
test "txdb_exons" {
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.Exon::new(
      "Ex2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      300,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], [], exons, [], "hg19", "Homo sapiens")

  let gr = @src.exons(txdb)
  assert_eq(@src.granges_length(gr), 2)
}

///|
test "txdb_cds" {
  let cds_regions = [
    @src.CDS::new(
      "CDS1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.CDS::new(
      "CDS2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      280,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], [], [], cds_regions, "hg19", "Homo sapiens")

  let gr = @src.cds(txdb)
  assert_eq(@src.granges_length(gr), 2)
}

///|
test "txdb_transcripts_by_gene" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 500, @src.strand_plus()),
    @src.Transcript::new("Tx2", "GeneA", "chr1", 110, 490, @src.strand_plus()),
    @src.Transcript::new("Tx3", "GeneB", "chr2", 200, 600, @src.strand_minus()),
  ]
  let txdb = @src.TxDb::new([], txs, [], [], "hg19", "Homo sapiens")

  let by_gene = @src.transcripts_by_gene(txdb)
  let geneA_txs = by_gene.get("GeneA")
  assert_true(geneA_txs.is_some())
  assert_eq(geneA_txs.unwrap().length(), 2)
  let geneB_txs = by_gene.get("GeneB")
  assert_true(geneB_txs.is_some())
  assert_eq(geneB_txs.unwrap().length(), 1)
}

///|
test "txdb_promoters" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 1000, 2000, @src.strand_plus()),
    @src.Transcript::new(
      "Tx2",
      "GeneB",
      "chr1",
      3000,
      4000,
      @src.strand_minus(),
    ),
  ]
  let txdb = @src.TxDb::new([], txs, [], [], "hg19", "Homo sapiens")

  let prom = @src.promoters(txdb, 200, 50)
  assert_eq(@src.granges_length(prom), 2)
  assert_eq(prom.starts[0], 800)
  assert_eq(prom.ends[0], 1049)
}

///|
test "txdb_gene_ids" {
  let genes = [
    @src.Gene::new("GeneA", "chr1", 100, 500, @src.strand_plus()),
    @src.Gene::new("GeneB", "chr2", 200, 600, @src.strand_minus()),
  ]
  let txdb = @src.TxDb::new(genes, [], [], [], "hg19", "Homo sapiens")

  let ids = @src.gene_ids(txdb)
  assert_eq(ids.length(), 2)
}

///|
test "txdb_transcript_ids" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 500, @src.strand_plus()),
    @src.Transcript::new("Tx2", "GeneA", "chr1", 110, 490, @src.strand_plus()),
  ]
  let txdb = @src.TxDb::new([], txs, [], [], "hg19", "Homo sapiens")

  let ids = @src.transcript_ids(txdb)
  assert_eq(ids.length(), 2)
}

///|
test "txdb_exon_ids" {
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.Exon::new(
      "Ex2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      300,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], [], exons, [], "hg19", "Homo sapiens")

  let ids = @src.exon_ids(txdb)
  assert_eq(ids.length(), 2)
}

///|
test "txdb_summary" {
  let genes = [@src.Gene::new("GeneA", "chr1", 100, 500, @src.strand_plus())]
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 500, @src.strand_plus()),
  ]
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
  ]
  let cds_regions = [
    @src.CDS::new(
      "CDS1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
  ]

  let txdb = @src.TxDb::new(
    genes, txs, exons, cds_regions, "hg19", "Homo sapiens",
  )

  let summary = @src.txdb_summary(txdb)
  assert_true(summary.contains("Genes: 1"))
  assert_true(summary.contains("Transcripts: 1"))
}

///|
test "txdb_parse_gtf_basic" {
  let gtf = "chr1\tensembl\tgene\t100\t500\t.\t+\t.\tgene_id \"GeneA\"; gene_name \"GeneA\";\n" +
    "chr1\tensembl\ttranscript\t100\t500\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\";\n" +
    "chr1\tensembl\texon\t100\t150\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"1\";\n" +
    "chr1\tensembl\texon\t200\t300\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"2\";\n" +
    "chr1\tensembl\tCDS\t100\t150\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"1\";\n" +
    "chr1\tensembl\tCDS\t200\t280\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"Tx1\"; exon_number \"2\";\n"

  let txdb = @src.parse_gtf(gtf)
  assert_eq(txdb.genes.length(), 1)
  assert_eq(txdb.transcripts.length(), 1)
  assert_eq(txdb.exons.length(), 2)
  assert_eq(txdb.cds.length(), 2)
}

///|
test "txdb_introns_by_transcript" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 300, @src.strand_plus()),
  ]
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.Exon::new(
      "Ex2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      300,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], txs, exons, [], "hg19", "Homo sapiens")

  let introns = @src.introns_by_transcript(txdb)
  let tx1_introns = introns.get("Tx1")
  assert_true(tx1_introns.is_some())
  assert_eq(tx1_introns.unwrap().length(), 1)
}

///|
test "txdb_five_prime_utrs" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 300, @src.strand_plus()),
  ]
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.Exon::new(
      "Ex2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      300,
      @src.strand_plus(),
      2,
    ),
  ]
  let cds_regions = [
    @src.CDS::new(
      "CDS1",
      "Tx1",
      "GeneA",
      "chr1",
      120,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.CDS::new(
      "CDS2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      280,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], txs, exons, cds_regions, "hg19", "Homo sapiens")

  let utrs = @src.five_prime_utrs(txdb)
  assert_eq(@src.granges_length(utrs), 1)
  assert_eq(utrs.starts[0], 100)
  assert_eq(utrs.ends[0], 119)
}

///|
test "txdb_three_prime_utrs" {
  let txs = [
    @src.Transcript::new("Tx1", "GeneA", "chr1", 100, 300, @src.strand_plus()),
  ]
  let exons = [
    @src.Exon::new(
      "Ex1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.Exon::new(
      "Ex2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      300,
      @src.strand_plus(),
      2,
    ),
  ]
  let cds_regions = [
    @src.CDS::new(
      "CDS1",
      "Tx1",
      "GeneA",
      "chr1",
      100,
      150,
      @src.strand_plus(),
      1,
    ),
    @src.CDS::new(
      "CDS2",
      "Tx1",
      "GeneA",
      "chr1",
      200,
      280,
      @src.strand_plus(),
      2,
    ),
  ]
  let txdb = @src.TxDb::new([], txs, exons, cds_regions, "hg19", "Homo sapiens")

  let utrs = @src.three_prime_utrs(txdb)
  assert_eq(@src.granges_length(utrs), 1)
  assert_eq(utrs.starts[0], 281)
  assert_eq(utrs.ends[0], 300)
}
