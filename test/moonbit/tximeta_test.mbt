///|
/// Tests for tximeta — metadata-aware transcript import.

///|
test "tximeta sample data is linked" {
  let data = @src.tximeta_sample_data()
  assert_true(data.linked)
  assert_true(@src.tximeta_validate_linked(data))
  assert_eq(data.n_samples(), 4)
  assert_eq(data.n_transcripts(), 3)
}

///|
test "tximeta unlinked validation fails" {
  let data = @src.TximetaData::new(
    sample_names=["s1"],
    transcript_names=["ENST1"],
    counts=[[5]],
    abundance=[[1.0]],
  )
  assert_false(data.linked)
  assert_false(@src.tximeta_validate_linked(data))
}

///|
test "tximeta link_txome injects provenance" {
  let data = @src.TximetaData::new(
    sample_names=["s1"],
    transcript_names=["ENST1"],
    counts=[[5]],
    abundance=[[1.0]],
  )
  let linked = @src.link_txome(
    data,
    @src.TximetaTxome::new(
      sha="deadbeef",
      genome="GRCm39",
      organism="Mus musculus",
      release="M25",
      source="GENCODE",
    ),
  )
  assert_true(@src.tximeta_validate_linked(linked))
  match linked.txome {
    Some(t) => {
      assert_eq(t.genome, "GRCm39")
      assert_eq(t.sha, "deadbeef")
    }
    None => assert_true(false)
  }
}

///|
test "tximeta metadata merge" {
  let data = @src.tximeta_sample_data()
  let meta = @src.tximeta_metadata(data)
  match meta.get("txome.genome") {
    Some(g) => assert_eq(g, "GRCh38")
    None => assert_true(false)
  }
  match meta.get("ctrl1.condition") {
    Some(c) => assert_eq(c, "control")
    None => assert_true(false)
  }
  match meta.get("treated2.condition") {
    Some(c) => assert_eq(c, "treated")
    None => assert_true(false)
  }
}

///|
test "tximeta summarize to gene" {
  let data = @src.tximeta_sample_data()
  let tx2gene = @src.tximeta_sample_tx2gene()
  let (genes, counts, abundance) = @src.tximeta_summarize_to_gene(data, tx2gene~)
  assert_eq(genes.length(), 2)
  let mut gene_a_row = -1
  for i in 0..<genes.length() {
    if genes[i] == "GENE_A" {
      gene_a_row = i
    }
  }
  assert_true(gene_a_row >= 0)
  assert_eq(counts[gene_a_row][0], 15)
  assert_eq(counts[gene_a_row][3], 67)
  assert_true((abundance[gene_a_row][0] - 1.5).abs() < 1.0e-9)
}

///|
test "tximeta summarize unknown transcript dropped" {
  let data = @src.TximetaData::new(
    sample_names=["s1"],
    transcript_names=["UNKNOWN_TX", "ENST1"],
    counts=[[3], [7]],
    abundance=[[0.3], [0.7]],
  )
  let (genes, counts, _) = @src.tximeta_summarize_to_gene(
    data,
    tx2gene=@src.tximeta_sample_tx2gene(),
  )
  assert_eq(genes.length(), 1)
  assert_eq(genes[0], "GENE_A")
  assert_eq(counts[0][0], 7)
}

///|
test "tximeta data summary string" {
  let data = @src.tximeta_sample_data()
  let s = @src.tximeta_data_summary(data)
  assert_true(s.contains("4 samples"))
  assert_true(s.contains("Homo sapiens"))
  assert_true(s.contains("GRCh38"))
}

///|
test "tximeta txome tags carried" {
  let data = @src.tximeta_sample_data()
  match data.txome {
    Some(t) => {
      assert_eq(t.versioned_tags.length(), 1)
      assert_eq(t.versioned_tags[0], "v38")
    }
    None => assert_true(false)
  }
}