// Tests for the seqio_unify text/container adapters: the unified
// seqio_parse entry must return SeqRecord results aligned with the
// standalone parsers for every newly integrated format (spec 5.1.1.1a),
// raise on unknown formats (spec 5.1.3.1) and return empty arrays for
// empty content (spec 5.1.3.4).

///|
const SWISS_SAMPLE : String = "ID   TEST_HUMAN             Standard;       200 AA.\nAC   P12345; Q67890;\nDE   RecName: Full=Test protein; EC=1.2.3.4;\nGN   Name=TEST1;\nOS   Homo sapiens (Human).\nOX   NCBI_TaxID=9606;\nOC   Eukaryota; Metazoa; Chordata.\nKW   Test; Example.\nFT   DOMAIN          Test domain.        1..100\nFT   ACT_SITE       Catalytic.         50..50\nRN   [1]\nRX   PUBMED=1234567; DOI=10.1000/test;\nRA   Smith A., Jones B.;\nRT   \"Test paper title.\";\nRL   J. Biol. Chem. 200:123-130(1995).\nCC   -!- FUNCTION: Catalyzes a test reaction that is important for testing.\nCC   -!- SUBCELLULAR LOCATION: Cytoplasm. Mitochondrion.\nCC   -!- TISSUE SPECIFICITY: Expressed in all tissues tested.\nDT   01-JAN-1990, integrated into UniProtKB/Swiss-Prot.\nDT   15-JUL-2020, sequence version 3.\nDT   01-APR-2024, entry version 100.\nSQ   SEQUENCE   26 AA;  2843 MW;  ABC123DEF456 CRC64;\n     MALWMRLLPL LALLALWGPD PAAAFV\n//\n"

///|
const IMGT_SAMPLE : String = ">HLA00001|A*01:01:01:01|A*01:01:01:01|ORF|2015-03-31|1996-01-19|human|HLA-A|I|MHC|2703 bp| genomic-DNA|\nATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC\n>HLA00003|B*07:02:01|B*07:02:01|ORF|2017-07-18|1996-01-19|human|HLA-B|I|MHC|2711 bp| genomic-DNA|\nATGGCGCCCCGAACCCTCCTCCTGCTGCTCTC\n"

///|
const GB_FEATURE_SAMPLE : String = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nDEFINITION  Test gene for unit testing.\nACCESSION   NM_001\nVERSION     NM_001.1\nFEATURES             Location/Qualifiers\n     CDS             1..100\n                     /gene=\"testGene\"\n                     /translation=\"MKVL\"\nORIGIN      \n        1 atgaaattg cgtc\n//\n"

///|
test "unify_swiss_parse" {
  let records = @src.parse_swiss_records(SWISS_SAMPLE)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.id, "TEST_HUMAN")
  assert_eq(rec.seq.to_string().length(), 26)
  assert_eq(rec.annotations["organism"], "Homo sapiens")
  assert_eq(rec.annotations["keywords"], "Test, Example")
  assert_eq(rec.annotations["accessions"], "P12345, Q67890")
  // features must be mapped
  assert_eq(rec.features.length(), 2)
  assert_eq(rec.features[0].feature_type, "DOMAIN")
  assert_eq(rec.features[1].feature_type, "ACT_SITE")
  // unified entry agrees with the adapter
  let via_entry = @src.seqio_parse(SWISS_SAMPLE, "swiss")
  assert_eq(via_entry.length(), 1)
  assert_eq(via_entry[0].id, rec.id)
  assert_eq(via_entry[0].seq.to_string(), rec.seq.to_string())
}

///|
test "unify_swiss_malformed_and_empty" {
  // empty content returns an empty array, not an error
  assert_eq(@src.parse_swiss_records("").length(), 0)
  // non-empty garbage raises the documented message
  let raised = try {
    ignore(@src.parse_swiss_records("this is not swiss"))
    false
  } catch {
    SeqIOError(msg) => msg == "Malformed swiss entry"
  }
  assert_true(raised)
}

///|
test "unify_uniprot_xml_parse" {
  let xml = @src.uniprot_sample_xml()
  let records = @src.parse_uniprot_xml_records(xml)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.id, "P12345")
  assert_true(rec.seq.to_string().length() > 0)
  assert_true(rec.annotations.contains("organism"))
  // unified entry agrees
  let via_entry = @src.seqio_parse(xml, "uniprot-xml")
  assert_eq(via_entry[0].id, rec.id)
  assert_eq(via_entry[0].seq.to_string(), rec.seq.to_string())
}

///|
test "unify_seqxml_parse" {
  let content = @src.seqxml_to_xml(@src.seqxml_sample_document())
  let records = @src.parse_seqxml_records(content)
  assert_true(records.length() > 0)
  let direct = @src.seqxml_to_seq_records(@src.seqxml_parse(content))
  assert_eq(records.length(), direct.length())
  assert_eq(records[0].id, direct[0].id)
  assert_eq(records[0].seq.to_string(), direct[0].seq.to_string())
  let via_entry = @src.seqio_parse(content, "seqxml")
  assert_eq(via_entry[0].id, direct[0].id)
}

///|
test "unify_imgt_parse" {
  let records = @src.parse_imgt_records(IMGT_SAMPLE)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "A*01:01:01:01")
  assert_eq(records[0].seq.to_string(), "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC")
  // the 13 header fields are mirrored into annotations
  assert_eq(records[0].annotations["species"], "human")
  assert_eq(records[0].annotations["gene"], "HLA-A")
  assert_eq(records[1].annotations["gene"], "HLA-B")
  let via_entry = @src.seqio_parse(IMGT_SAMPLE, "imgt")
  assert_eq(via_entry[0].id, records[0].id)
}

///|
test "unify_ig_parse" {
  let content = @src.ig_sample_text()
  let records = @src.parse_ig_records(content)
  let direct = @src.ig_to_seq_records(@src.ig_parse(content))
  assert_eq(records.length(), direct.length())
  assert_eq(records[0].id, "A_U455")
  assert_eq(records[0].seq.to_string(), direct[0].seq.to_string())
  let via_entry = @src.seqio_parse(content, "ig")
  assert_eq(via_entry[0].id, "A_U455")
}

///|
test "unify_insdc_genbank_features" {
  let records = @src.parse_insdc_records(GB_FEATURE_SAMPLE)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.id, "NM_001")
  assert_eq(rec.name, "NM_001")
  assert_eq(rec.description, "Test gene for unit testing")
  // feature table preserved through SeqRecord.features (spec 5.1.1 rule 4a)
  assert_eq(rec.features.length(), 1)
  assert_eq(rec.features[0].feature_type, "CDS")
  assert_eq(rec.features[0].id, "1..100")
  let genes = rec.features[0].qualifiers.get("gene")
  assert_true(genes is Some(_))
  match genes {
    Some(vals) => assert_eq(vals[0], "testGene")
    None => assert_true(false)
  }
  // annotations carry the header metadata
  assert_eq(rec.annotations["molecule_type"], "DNA")
  assert_eq(rec.annotations["version"], "NM_001.1")
  // unified entry agrees
  let via_entry = @src.seqio_parse(GB_FEATURE_SAMPLE, "insdc")
  assert_eq(via_entry[0].id, rec.id)
  assert_eq(via_entry[0].features.length(), 1)
}

///|
test "unify_insdc_multi_record_and_embl" {
  // two GenBank records in one document
  let two = GB_FEATURE_SAMPLE + GB_FEATURE_SAMPLE
  assert_eq(@src.parse_insdc_records(two).length(), 2)
  // an EMBL record is recognized by its ID header
  let embl = "ID   X56734; SV 1; linear; mRNA; STD; PRO; 1859 BP.\nSQ   Sequence 1859 BP;\n     atgc\n//\n"
  let records = @src.parse_insdc_records(embl)
  assert_eq(records.length(), 1)
  assert_eq(records[0].name, "X56734")
  // unrecognized header raises
  let raised = try {
    ignore(@src.parse_insdc_records("WHAT   ever\n//\n"))
    false
  } catch {
    SeqIOError(msg) => msg.contains("Unrecognized INSDC record")
  }
  assert_true(raised)
}

///|
test "unify_phd_parse" {
  let content = @src.phd_sample_text()
  let records = @src.parse_phd_records(content)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.id, "sample_read")
  assert_eq(rec.seq.to_string(), "ACGTACGTAC")
  match rec.phred_quality() {
    Some(quals) => {
      assert_eq(quals.length(), 10)
      assert_eq(quals[0], 35)
      assert_eq(quals[9], 44)
    }
    None => assert_true(false)
  }
  let via_entry = @src.seqio_parse(content, "phd")
  assert_eq(via_entry[0].id, rec.id)
  assert_eq(via_entry[0].seq.to_string(), rec.seq.to_string())
}

///|
test "unify_ace_parse" {
  let content = @src.sample_ace_content()
  let records = @src.parse_ace_records(content)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.id, "contig1")
  assert_eq(rec.seq.to_string().length(), 60)
  assert_eq(rec.annotations["reads_count"], "2")
  assert_eq(rec.annotations["length"], "60")
  // base qualities match the consensus length -> phred_quality attached
  match rec.phred_quality() {
    Some(quals) => assert_eq(quals.length(), 60)
    None => assert_true(false)
  }
  let via_entry = @src.seqio_parse(content, "ace")
  assert_eq(via_entry[0].id, rec.id)
}

///|
test "unify_gfa_parse" {
  let content = "H\tVN:Z:1.0\nS\tseq1\tACGTACGT\tLN:i:8\nS\tseq2\tTTTTGGGG\tLN:i:8\nL\tseq1\t+\tseq2\t-\t8M\nP\tpath1\tseq1+,seq2-\t8M,4M\n"
  let records = @src.parse_gfa_records(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGTACGT")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "TTTTGGGG")
  let via_entry = @src.seqio_parse(content, "gfa")
  assert_eq(via_entry.length(), 2)
  assert_eq(via_entry[0].id, "seq1")
}

///|
test "unify_pdb_seqres_and_atom" {
  let seqres_text = "SEQRES   1 A    4  ALA GLY VAL LEU\n//\n"
  let records = @src.parse_pdb_seqres_records(seqres_text)
  assert_true(records.length() > 0)
  let direct = @src.pdb_seqres_parser(seqres_text)
  assert_eq(records[0].id, direct[0].id)
  assert_eq(records[0].seq.to_string(), direct[0].seq.to_string())
  let via_entry = @src.seqio_parse(seqres_text, "pdb-seqres")
  assert_eq(via_entry[0].id, direct[0].id)
  // empty atom content yields no records
  assert_eq(@src.parse_pdb_atom_records("").length(), 0)
}

///|
test "unify_new_format_branches" {
  // fasta-2line round-trips through the unified entry
  let f2l = ">s1 desc1\nACGTACGT\n"
  let recs = @src.seqio_parse(f2l, "fasta-2line")
  assert_eq(recs.length(), 1)
  assert_eq(recs[0].id, "s1")
  assert_eq(recs[0].seq.to_string(), "ACGTACGT")
  // fasta-2line output can be read back through the unified entry
  let written = @src.write_fasta_2line(recs)
  let reparsed = @src.seqio_parse(written, "fasta-2line")
  assert_eq(reparsed[0].id, "s1")
  assert_eq(reparsed[0].seq.to_string(), "ACGTACGT")
}

///|
test "unify_empty_content_returns_empty" {
  assert_eq(@src.parse_swiss_records("").length(), 0)
  assert_eq(@src.parse_uniprot_xml_records("").length(), 0)
  assert_eq(@src.parse_seqxml_records("").length(), 0)
  assert_eq(@src.parse_imgt_records("").length(), 0)
  assert_eq(@src.parse_ig_records("").length(), 0)
  assert_eq(@src.parse_insdc_records("").length(), 0)
  assert_eq(@src.parse_phd_records("").length(), 0)
  assert_eq(@src.parse_ace_records("").length(), 0)
  assert_eq(@src.parse_gfa_records("").length(), 0)
  assert_eq(@src.parse_pdb_seqres_records("").length(), 0)
  assert_eq(@src.parse_pdb_atom_records("").length(), 0)
}

///|
test "unify_unknown_format_raises" {
  let raised = try {
    ignore(@src.seqio_parse(">s1\nACGT\n", "definitely_not_a_format"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: definitely_not_a_format"
  }
  assert_true(raised)
}
