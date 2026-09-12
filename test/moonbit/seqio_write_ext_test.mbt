// Tests for the seqio_write_ext adapters: every new write format must
// round-trip through its unified parse counterpart with id and sequence
// preserved (spec 5.3.1.2), phd must reject records without phred quality
// scores (spec 5.3.1.4), non-writable formats must list the writable set
// (spec 5.3.1.3), and empty record arrays produce empty output
// (spec 5.3.3.3).

///|
fn writeext_sample_record() -> @src.SeqRecord {
  let rec = @src.SeqRecord::new(
    @src.Seq::new("ACGTACGTAC"),
    id="rec1",
    name="rec1",
    description="sample record",
  )
  rec
}

///|
test "write_swiss_records_roundtrip" {
  let rec = writeext_sample_record()
  rec.annotations["organism"] = "Homo sapiens"
  rec.annotations["keywords"] = "Test, Example"
  rec.annotations["accessions"] = "P12345, Q67890"
  let text = @src.write_swiss_records([rec])
  let parsed = @src.parse_swiss_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  assert_eq(parsed[0].description, "sample record")
  assert_eq(parsed[0].annotations["organism"], "Homo sapiens")
  assert_eq(parsed[0].annotations["keywords"], "Test, Example")
  assert_eq(parsed[0].annotations["accessions"], "P12345, Q67890")
}

///|
test "write_pir_records_roundtrip" {
  let rec = writeext_sample_record()
  let text = @src.write_pir_records([rec])
  let parsed = @src.parse_pir(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  assert_eq(parsed[0].description, "sample record")
}

///|
test "write_embl_records_roundtrip" {
  let rec = writeext_sample_record()
  rec.annotations["molecule_type"] = "DNA"
  let location = @src.FeatureLocation::new_exact(1, 6)
  let quals : Map[String, Array[String]] = Map([("gene", ["testGene"])])
  rec.features.push(
    @src.SeqFeature::new(
      location=Some(location),
      feature_type="CDS",
      id="1..6",
      qualifiers=quals,
    ),
  )
  let text = @src.write_embl_records([rec])
  // features survive through the insdc adapter
  let parsed = @src.parse_insdc_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  assert_eq(parsed[0].features.length(), 1)
  assert_eq(parsed[0].features[0].feature_type, "CDS")
  // unified entry agrees
  let via_entry = @src.parse_embl(text)
  assert_eq(via_entry.length(), 1)
  assert_eq(via_entry[0].seq.to_string(), "ACGTACGTAC")
}

///|
test "write_imgt_records_roundtrip" {
  let rec = writeext_sample_record()
  rec.annotations["species"] = "human"
  rec.annotations["gene"] = "HLA-A"
  let text = @src.write_imgt_records([rec])
  let parsed = @src.parse_imgt_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  assert_eq(parsed[0].annotations["species"], "human")
  assert_eq(parsed[0].annotations["gene"], "HLA-A")
}

///|
test "write_ig_records_roundtrip" {
  let rec = writeext_sample_record()
  let text = @src.write_ig_records([rec])
  let parsed = @src.parse_ig_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
}

///|
test "write_seqxml_records_roundtrip" {
  let rec = writeext_sample_record()
  let text = @src.write_seqxml_records([rec])
  let parsed = @src.parse_seqxml_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
}

///|
test "write_phd_records_roundtrip" {
  let rec = writeext_sample_record()
  rec.letter_annotations["phred_quality"] = [
    10, 20, 30, 40, 15, 25, 35, 45, 12, 22,
  ]
  let text = @src.write_phd_records([rec])
  let parsed = @src.parse_phd_records(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  match parsed[0].phred_quality() {
    Some(quals) => {
      assert_eq(quals.length(), 10)
      assert_eq(quals[0], 10)
      assert_eq(quals[9], 22)
    }
    None => assert_true(false)
  }
}

///|
test "write_phd_records_requires_quality" {
  // missing quality track raises before producing output
  let rec = writeext_sample_record()
  let raised = try {
    ignore(@src.write_phd_records([rec]))
    false
  } catch {
    SeqIOError(msg) => msg == "record 'rec1' lacks phred quality scores"
  }
  assert_true(raised)
  // mismatched quality length raises as well
  rec.letter_annotations["phred_quality"] = [10, 20, 30]
  let mismatch = try {
    ignore(@src.write_phd_records([rec]))
    false
  } catch {
    SeqIOError(msg) => msg.contains("phred quality length does not match")
  }
  assert_true(mismatch)
}

///|
test "seqio_write_new_format_branches" {
  let rec = writeext_sample_record()
  // fasta-2line round-trips through the unified entry
  let text = @src.seqio_write([rec], "fasta-2line")
  let parsed = @src.seqio_parse(text, "fasta-2line")
  assert_eq(parsed[0].id, "rec1")
  assert_eq(parsed[0].seq.to_string(), "ACGTACGTAC")
  // the remaining branches produce parseable output
  for format in ["embl", "pir", "seqxml", "swiss", "ig", "imgt"] {
    let out = @src.seqio_write([rec], format)
    let back = @src.seqio_parse(out, format)
    assert_eq(back.length(), 1)
    assert_eq(back[0].id, "rec1")
  }
  // phd requires quality scores
  rec.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  let phd_text = @src.seqio_write([rec], "phd")
  let phd_back = @src.seqio_parse(phd_text, "phd")
  assert_eq(phd_back[0].seq.to_string(), "ACGTACGTAC")
}

///|
test "seqio_write_non_writable_format_lists_supported" {
  let rec = writeext_sample_record()
  let raised = try {
    ignore(@src.seqio_write([rec], "2bit"))
    false
  } catch {
    SeqIOError(msg) =>
      msg.contains("format does not support writing") &&
      msg.contains("2bit") &&
      msg.contains("fasta")
  }
  assert_true(raised)
  // uniprot-xml is parse-only too
  let raised2 = try {
    ignore(@src.seqio_write([rec], "uniprot-xml"))
    false
  } catch {
    SeqIOError(msg) => msg.contains("format does not support writing")
  }
  assert_true(raised2)
  // unknown formats keep the legacy message
  let raised3 = try {
    ignore(@src.seqio_write([rec], "bogus"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: bogus"
  }
  assert_true(raised3)
}

///|
test "write_ext_empty_records" {
  assert_eq(@src.write_swiss_records([]), "")
  assert_eq(@src.write_pir_records([]), "")
  assert_eq(@src.write_embl_records([]), "")
  assert_eq(@src.write_imgt_records([]), "")
  assert_eq(@src.write_ig_records([]), "")
  assert_eq(@src.write_phd_records([]), "")
}
