///|
test "ballgown_bgtranscript_new" {
  let fpkm = [1.0, 2.0, 3.0]
  let cov = [10.0, 20.0, 30.0]
  let tx = @src.BgTranscript::new(
    "tx1",
    "gene1",
    "GeneA",
    "chr1",
    "+",
    3,
    fpkm,
    cov,
  )

  assert_eq(tx.id, "tx1")
  assert_eq(tx.gene_id, "gene1")
  assert_eq(tx.gene_name, "GeneA")
  assert_eq(tx.chr, "chr1")
  assert_eq(tx.strand, "+")
  assert_eq(tx.num_exons, 3)
  assert_eq(tx.fpkm.length(), 3)
  assert_eq(tx.cov.length(), 3)
}

///|
test "ballgown_bgexon_new" {
  let exon = @src.BgExon::new(
    "exon1",
    "chr1",
    100,
    200,
    "+",
    ["tx1", "tx2"],
  )

  assert_eq(exon.id, "exon1")
  assert_eq(exon.chr, "chr1")
  assert_eq(exon.start, 100)
  assert_eq(exon.end, 200)
  assert_eq(exon.strand, "+")
  assert_eq(exon.transcript_ids.length(), 2)
}

///|
test "ballgown_bgintron_new" {
  let intron = @src.BgIntron::new(
    "intron1",
    "chr1",
    201,
    300,
    "+",
    ["tx1"],
  )

  assert_eq(intron.id, "intron1")
  assert_eq(intron.chr, "chr1")
  assert_eq(intron.start, 201)
  assert_eq(intron.end, 300)
  assert_eq(intron.strand, "+")
  assert_eq(intron.transcript_ids.length(), 1)
}

///|
test "ballgown_deresult_new" {
  let de = @src.DeResult::new(
    "tx1",
    "gene1",
    "GeneA",
    1.5,
    0.01,
    0.02,
    10.0,
  )

  assert_eq(de.feature_id, "tx1")
  assert_eq(de.gene_id, "gene1")
  assert_eq(de.gene_name, "GeneA")
  assert_eq(de.log2_fc, 1.5)
  assert_eq(de.p_value, 0.01)
  assert_eq(de.q_value, 0.02)
  assert_eq(de.mean_fpkm, 10.0)
}

///|
test "ballgown_sample_ballgown" {
  let bg = @src.sample_ballgown()

  assert_eq(bg.sample_names.length(), 6)
  assert_true(bg.transcripts.length() > 0)
  assert_true(bg.exons.length() > 0)
  assert_true(bg.introns.length() > 0)
}

///|
test "ballgown_get_fpkm_matrix" {
  let bg = @src.sample_ballgown()
  let fpkm_mat = bg.get_fpkm_matrix()

  assert_eq(fpkm_mat.length(), bg.transcripts.length())
  assert_eq(fpkm_mat[0].length(), 6)
}

///|
test "ballgown_get_cov_matrix" {
  let bg = @src.sample_ballgown()
  let cov_mat = bg.get_cov_matrix()

  assert_eq(cov_mat.length(), bg.transcripts.length())
  assert_eq(cov_mat[0].length(), 6)
}

///|
test "ballgown_get_transcript" {
  let bg = @src.sample_ballgown()
  let tx_id = bg.transcripts[0].id
  let tx = bg.get_transcript(tx_id)

  match tx {
    Some(t) => assert_eq(t.id, tx_id)
    None => assert_true(false)
  }
}

///|
test "ballgown_get_gene_ids" {
  let bg = @src.sample_ballgown()
  let gene_ids = bg.get_gene_ids()

  assert_eq(gene_ids.length(), 5)
}

///|
test "ballgown_get_gene_names" {
  let bg = @src.sample_ballgown()
  let gene_names = bg.get_gene_names()

  assert_eq(gene_names.length(), 5)
}

///|
test "ballgown_gene_expression" {
  let bg = @src.sample_ballgown()
  let gene_ids = bg.get_gene_ids()
  let expr = @src.gene_expression(bg, gene_ids[0])

  assert_eq(expr.length(), 6)
}

///|
test "ballgown_stattest_transcript" {
  let bg = @src.sample_ballgown()
  let group = [0, 0, 0, 1, 1, 1]

  let results = @src.stattest(bg, group, feature = "transcript")

  assert_eq(results.length(), bg.transcripts.length())
  for r in results {
    assert_true(r.p_value >= 0.0 && r.p_value <= 1.0)
    assert_true(r.q_value >= 0.0 && r.q_value <= 1.0)
  }
}

///|
test "ballgown_stattest_gene" {
  let bg = @src.sample_ballgown()
  let group = [0, 0, 0, 1, 1, 1]

  let results = @src.stattest(bg, group, feature = "gene")

  assert_eq(results.length(), 5)
  for r in results {
    assert_true(r.p_value >= 0.0 && r.p_value <= 1.0)
    assert_true(r.q_value >= 0.0 && r.q_value <= 1.0)
  }
}

///|
test "ballgown_filter_low_expression" {
  let bg = @src.sample_ballgown()
  let filtered = @src.filter_low_expression(bg, min_fpkm = 0.5)

  assert_true(filtered.transcripts.length() <= bg.transcripts.length())
}
