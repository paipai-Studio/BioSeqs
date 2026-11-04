///|
/// Tests for the Biopython Bio.UniGene-compatible parser.

///|
test "unigene parses core record fields" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.identifier, "Hs.12345")
  assert_eq(record.species, "Hs")
  assert_eq(record.title, "Example kinase family member")
  assert_eq(record.symbol, "EXK1")
  assert_eq(record.cytoband, "7q31.2")
  assert_eq(record.gene_id, "12345")
  assert_eq(record.locuslink, "12345")
  assert_eq(record.chromosome, "7")
  assert_eq(record.declared_sequence_count, 3)
}

///|
test "unigene parses expression and homology fields" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.expression, ["brain", "heart", "adult"])
  assert_eq(record.restricted_expression, ["brain", "adult"])
  assert_eq(record.homology, Some(true))
  assert_eq(record.genomic_terminus, "T")
}

///|
test "unigene parses protein similarity fields" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.protein_similarities.length(), 1)
  let similarity = record.protein_similarities[0]
  assert_eq(similarity.organism, "9606")
  assert_eq(similarity.protein_gi, "123456")
  assert_eq(similarity.protein_id, "NP_000001.1")
  assert_eq(similarity.percent, "99.50")
  assert_eq(similarity.alignment_length, "401")
  assert_eq(similarity.field("PCT"), Some("99.50"))
}

///|
test "unigene parses sequence entries and IMAGE clones" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.sequences.length(), 3)
  let transcript = record.sequences[0]
  assert_eq(transcript.accession, "NM_000001.2")
  assert_eq(transcript.nucleotide_id, "g123455")
  assert_eq(transcript.protein_id, "g123456")
  assert_eq(transcript.sequence_type, "mRNA")
  let image = record.sequences[1]
  assert_true(image.is_image)
  assert_eq(image.image_id, "123456")
  assert_eq(image.read_end, "5'")
  assert_eq(image.library_id, "100")
  assert_eq(image.trace, "900001")
}

///|
test "unigene parses STS and transcript map entries" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.sts.length(), 1)
  assert_eq(record.sts[0].accession, "G12345")
  assert_eq(record.sts[0].unists, "76543")
  assert_eq(record.transcript_maps.length(), 1)
  assert_eq(record.transcript_maps[0].marker, "D7S1234")
  assert_eq(record.transcript_maps[0].radiation_hybrid_panel, "GB4")
}

///|
test "unigene sequence type and accession queries" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.sequences_of_type("EST").length(), 2)
  assert_eq(record.sequences_of_type("mRNA").length(), 1)
  assert_eq(record.sequences_of_type("HTC").length(), 0)
  match record.find_sequence("AA000002.1") {
    Some(sequence) => assert_eq(sequence.mgc, "600")
    None => abort("expected sequence accession")
  }
  assert_true(record.find_sequence("missing") is None)
}

///|
test "unigene similarity and IMAGE queries" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.protein_similarities_for("9606").length(), 1)
  assert_eq(record.protein_similarities_for("10090").length(), 0)
  let image_sequences = record.image_sequences()
  assert_eq(image_sequences.length(), 1)
  assert_eq(image_sequences[0].accession, "AA000001.1")
}

///|
test "unigene preserves unknown child fields" {
  let text = "ID          Mm.7\n" +
    "PROTSIM     ORG=10090; PROTID=NP_1; EVALUE=1e-30\n" +
    "TXMAP       MARKER=D1Mit1; METHOD=RH\n" +
    "SCOUNT      1\n" +
    "SEQUENCE    ACC=NM_1; SEQTYPE=mRNA; STATUS=reviewed\n" +
    "STS         UNISTS=10 PANEL=RH\n" +
    "//\n"
  let record = @src.unigene_read(text)
  assert_eq(record.protein_similarities[0].field("EVALUE"), Some("1e-30"))
  assert_eq(record.transcript_maps[0].field("METHOD"), Some("RH"))
  assert_eq(record.sequences[0].field("STATUS"), Some("reviewed"))
  assert_eq(record.sts[0].field("PANEL"), Some("RH"))
  let reparsed = @src.unigene_read(record.to_string())
  assert_eq(reparsed.protein_similarities[0].field("EVALUE"), Some("1e-30"))
  assert_eq(reparsed.sequences[0].field("STATUS"), Some("reviewed"))
}

///|
test "unigene parses multiple records" {
  let second = "ID          Mm.2\n" +
    "TITLE       Second cluster\n" +
    "GENE        Gene2\n" +
    "HOMOL       NO\n" +
    "SCOUNT      0\n" +
    "//\n"
  let records = @src.unigene_parse(@src.unigene_sample_text() + second)
  assert_eq(records.length(), 2)
  assert_eq(records[0].identifier, "Hs.12345")
  assert_eq(records[1].identifier, "Mm.2")
  assert_eq(records[1].species, "Mm")
  assert_eq(records[1].homology, Some(false))
}

///|
test "unigene read accepts exactly one record" {
  let record = @src.unigene_read(@src.unigene_sample_text())
  assert_eq(record.identifier, "Hs.12345")
  assert_true(record.summary().contains("sequences=3"))
}

///|
test "unigene read rejects empty input" {
  let raised = try {
    ignore(@src.unigene_read(""))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene read rejects multiple records" {
  let raised = try {
    ignore(
      @src.unigene_read(@src.unigene_sample_text() + @src.unigene_sample_text()),
    )
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene validates declared sequence count" {
  let text = "ID          Hs.1\n" +
    "SCOUNT      2\n" +
    "SEQUENCE    ACC=NM_1; SEQTYPE=mRNA\n" +
    "//\n"
  let raised = try {
    ignore(@src.unigene_parse(text))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene requires sequence count" {
  let raised = try {
    ignore(@src.unigene_parse("ID          Hs.1\n//\n"))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene rejects invalid homology value" {
  let text = "ID          Hs.1\nHOMOL       MAYBE\nSCOUNT      0\n//\n"
  let raised = try {
    ignore(@src.unigene_parse(text))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene rejects unknown top-level tags" {
  let text = "ID          Hs.1\nUNKNOWN     value\nSCOUNT      0\n//\n"
  let raised = try {
    ignore(@src.unigene_parse(text))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene rejects malformed child fields" {
  let text = "ID          Hs.1\n" +
    "SCOUNT      1\n" +
    "SEQUENCE    ACC=NM_1; INVALID\n" +
    "//\n"
  let raised = try {
    ignore(@src.unigene_parse(text))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene rejects unterminated records" {
  let raised = try {
    ignore(@src.unigene_parse("ID          Hs.1\nSCOUNT      0\n"))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene rejects malformed fixed-width lines" {
  let raised = try {
    ignore(@src.unigene_parse("ID Hs.1\n"))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}

///|
test "unigene supports CRLF input" {
  let crlf = @src.unigene_sample_text().replace_all(old="\n", new="\r\n")
  let record = @src.unigene_read(crlf)
  assert_eq(record.identifier, "Hs.12345")
  assert_eq(record.sequences.length(), 3)
}

///|
test "unigene record serialization round trip" {
  let original = @src.unigene_read(@src.unigene_sample_text())
  let serialized = original.to_string()
  assert_true(serialized.contains("ID          Hs.12345"))
  assert_true(serialized.contains("SCOUNT      3"))
  let reparsed = @src.unigene_read(serialized)
  assert_true(reparsed == original)
}

///|
test "unigene constructors serialize canonical records" {
  let record = @src.UniGeneRecord::new(
    identifier="Rn.9",
    title="Constructed record",
    symbol="Gene9",
    homology=Some(false),
    sequences=[
      @src.UniGeneSequence::new(accession="NM_9", sequence_type="mRNA"),
    ],
    sts=[@src.UniGeneSTS::new(unists="999")],
  )
  let reparsed = @src.unigene_read(record.to_string())
  assert_eq(reparsed.species, "Rn")
  assert_eq(reparsed.homology, Some(false))
  assert_eq(reparsed.sequences[0].accession, "NM_9")
  assert_eq(reparsed.sts[0].unists, "999")
}

///|
test "unigene rejects duplicate sequence counts" {
  let text = "ID          Hs.1\nSCOUNT      0\nSCOUNT      0\n//\n"
  let raised = try {
    ignore(@src.unigene_parse(text))
    false
  } catch {
    UniGeneError(_) => true
  }
  assert_true(raised)
}
