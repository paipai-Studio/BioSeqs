///|
test "swissprot_sample_record" {
  let record = @src.sample_swissprot_record()

  assert_eq(record.entry_name, "GLNA_HUMAN")
  assert_true(record.accession_numbers.length() > 0)
  assert_true(record.sequence.length() > 0)
  assert_true(record.seq_length > 0)
  assert_true(record.keywords.length() > 0)
  assert_true(record.features.length() > 0)
  assert_true(record.references.length() > 0)
}

///|
test "swissprot_parse" {
  let text = "ID   TEST_HUMAN             Standard;       200 AA.\nAC   P12345; Q67890;\nDE   RecName: Full=Test protein; EC=1.2.3.4;\nGN   Name=TEST1;\nOS   Homo sapiens (Human).\nOX   NCBI_TaxID=9606;\nOC   Eukaryota; Metazoa; Chordata.\nKW   Test; Example.\nFT   DOMAIN          1     100      Test domain.\nFT   ACT_SITE       50      50      Catalytic residue.\nRN   [1]\nRX   PUBMED=1234567; DOI=10.1000/test;\nRA   Smith A., Jones B.;\nRT   \"Test paper title.\";\nRL   J. Biol. Chem. 200:123-130(1995).\nCC   -!- FUNCTION: Catalyzes a test reaction that is important for testing.\nCC   -!- SUBCELLULAR LOCATION: Cytoplasm. Mitochondrion.\nCC   -!- TISSUE SPECIFICITY: Expressed in all tissues tested.\nDT   01-JAN-1990, integrated into UniProtKB/Swiss-Prot.\nDT   15-JUL-2020, sequence version 3.\nDT   01-APR-2024, entry version 100.\nSQ   SEQUENCE   26 AA;  2843 MW;  ABC123DEF456 CRC64;\n     MALWMRLLPL LALLALWGPD PAAAFVN\n//\n"

  let record = @src.parse_swissprot(text)

  assert_true(record is Some(_))
  match record {
    Some(rec) => {
      assert_eq(rec.entry_name, "TEST_HUMAN")
      assert_eq(rec.accession_numbers.length(), 2)
      assert_eq(rec.organism, "Homo sapiens")
      assert_eq(rec.organism_id, "9606")
      assert_eq(rec.seq_length, 26)
      assert_eq(rec.keywords.length(), 2)
      assert_eq(rec.ec_numbers.length(), 1)
    }
    None => assert_true(false)
  }
}

///|
test "swissprot_primary_accession" {
  let record = @src.sample_swissprot_record()

  let primary = record.primary_accession()
  assert_true(primary is Some(_))
  match primary {
    Some(p) => assert_eq(p, "P15104")
    None => assert_true(false)
  }
}

///|
test "swissprot_primary_gene" {
  let record = @src.sample_swissprot_record()

  let gene = record.primary_gene()
  assert_true(gene is Some(_))
  match gene {
    Some(g) => assert_eq(g, "GLUL")
    None => assert_true(false)
  }
}

///|
test "swissprot_get_features" {
  let record = @src.sample_swissprot_record()

  let domains = record.get_features("DOMAIN")
  assert_eq(domains.length(), 2)

  let act_sites = record.get_features("ACT_SITE")
  assert_eq(act_sites.length(), 1)
}

///|
test "swissprot_reference_fields" {
  let record = @src.sample_swissprot_record()

  assert_true(record.references.length() > 0)
  let first_ref = record.references[0]
  assert_eq(first_ref.number, 1)
  assert_true(first_ref.title.length() > 0)
  assert_true(first_ref.authors.length() > 0)
  assert_true(first_ref.pubmed_id.length() > 0)
}

///|
test "swissprot_feature_fields" {
  let record = @src.sample_swissprot_record()

  assert_true(record.features.length() > 0)
  let first_feat = record.features[0]
  assert_true(first_feat.feature_type.length() > 0)
  assert_true(first_feat.start > 0)
  assert_true(first_feat.end > 0)
}

///|
test "swissprot_parse_empty" {
  let text = "//\n"
  let record = @src.parse_swissprot(text)
  assert_true(record is None)
}
