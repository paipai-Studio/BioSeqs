// Tests for the seqio index extension (spec 5.4): qual, swiss, and seqxml
// become randomly accessible through the unified seqio_index entry point.
// The index record count must equal the number of records in the content
// (spec 5.4.1.1a), random access must agree with a full parse of the same
// record (spec 5.4.1.1b), get_raw must return the exact byte slice from
// the original content (spec 5.4.1.3), duplicate identifiers resolve to
// the last occurrence, non-indexable formats raise, missing identifiers
// echo the key back (spec 5.4.3.1 / 5.4.3.2), and an unterminated seqxml
// entry is a boundary error pinpointing the offset (spec 5.4.3.3).

///|
test "seqio_index_qual_records_and_access" {
  let qual = ">qual1 first\n10 20 30\n>qual2 second\n40 50 60 70\n"
  let index = @src.seqio_index(qual, "qual")
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["qual1", "qual2"])
  // get_raw returns the exact slice from the original content
  assert_eq(index.get_raw(qual, "qual2"), ">qual2 second\n40 50 60 70\n")
  // random access agrees with a full parse of the record
  let rec = index.get(qual, "qual1")
  let full = @src.parse_qual(qual)
  assert_eq(rec.id, full[0].id)
  assert_eq(rec.seq.to_string(), full[0].seq.to_string())
  match (rec.phred_quality(), full[0].phred_quality()) {
    (Some(a), Some(b)) => assert_eq(a, b)
    _ => assert_true(false)
  }
}

///|
test "seqio_index_swiss_records_and_access" {
  let rec1 = @src.SeqRecord::new(
    @src.Seq::new("ACGTACGTAC"),
    id="rec1",
    name="rec1",
    description="swiss index demo",
  )
  let rec2 = @src.SeqRecord::new(
    @src.Seq::new("TTTGTTTGTT"),
    id="rec2",
    name="rec2",
    description="swiss index demo two",
  )
  let swiss = @src.write_swiss_records([rec1, rec2])
  let index = @src.seqio_index(swiss, "swiss")
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["rec1", "rec2"])
  // the raw slice is the exact entry text starting at the ID line
  let raw = index.get_raw(swiss, "rec1")
  assert_true(raw.has_prefix("ID   rec1"))
  let offset = swiss.find("ID   rec1").unwrap()
  assert_eq(raw, swiss[offset:offset + raw.length()].to_owned())
  // random access agrees with the full parse
  let rec = index.get(swiss, "rec2")
  let full = @src.parse_swiss_records(swiss)
  assert_eq(rec.id, full[1].id)
  assert_eq(rec.seq.to_string(), full[1].seq.to_string())
}

///|
test "seqio_index_seqxml_records_and_access" {
  let seqxml = "<?xml version=\"1.0\"?>\n<seqXML>\n  <entry id=\"s1\" desc=\"one\">\n    <dna>ACGTAC</dna>\n  </entry>\n  <entry id=\"s2\" desc=\"two\">\n    <dna>TTTGCA</dna>\n  </entry>\n</seqXML>\n"
  let index = @src.seqio_index(seqxml, "seqxml")
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["s1", "s2"])
  // raw slice starts at the <entry line and ends at </entry>
  let raw = index.get_raw(seqxml, "s2")
  assert_true(raw.trim().has_prefix("<entry"))
  assert_true(raw.trim().has_suffix("</entry>"))
  // random access agrees with the full parse of that entry
  let rec = index.get(seqxml, "s1")
  let full = @src.parse_seqxml_records(seqxml)
  assert_eq(rec.id, full[0].id)
  assert_eq(rec.seq.to_string(), full[0].seq.to_string())
  assert_eq(rec.seq.to_string(), "ACGTAC")
}

///|
test "seqio_index_duplicate_ids_last_wins" {
  // seqxml with two entries sharing an id: the last one wins, matching
  // the fasta/genbank index behaviour (spec 5.4.1.4)
  let dup = "<?xml version=\"1.0\"?>\n<seqXML>\n  <entry id=\"dup\">\n    <dna>AAA</dna>\n  </entry>\n  <entry id=\"dup\">\n    <dna>CCC</dna>\n  </entry>\n</seqXML>\n"
  let index = @src.seqio_index(dup, "seqxml")
  assert_eq(index.length(), 2)
  let rec = index.get(dup, "dup")
  assert_eq(rec.seq.to_string(), "CCC")
  // fasta duplicate ids behave the same way
  let fasta = ">d\nAAA\n>d\nCCC\n"
  let findex = @src.seqio_index(fasta, "fasta")
  let frec = findex.get(fasta, "d")
  assert_eq(frec.seq.to_string(), "CCC")
}

///|
test "seqio_index_non_indexable_format_raises" {
  // binary formats are not indexable (spec 5.4.3.1)
  let raised = try {
    ignore(@src.seqio_index("anything", "2bit"))
    false
  } catch {
    SeqIOError(msg) => msg == "Cannot index format: 2bit"
  }
  assert_true(raised)
  // unknown formats raise the same way
  let raised2 = try {
    ignore(@src.seqio_index("anything", "bogus"))
    false
  } catch {
    SeqIOError(msg) => msg == "Cannot index format: bogus"
  }
  assert_true(raised2)
}

///|
test "seqio_index_missing_key_echoes_identifier" {
  let qual = ">q1\n10 20\n"
  let index = @src.seqio_index(qual, "qual")
  let raised = try {
    ignore(index.get_raw(qual, "absent"))
    false
  } catch {
    SeqIOError(msg) => msg == "Key absent not found in index"
  }
  assert_true(raised)
}

///|
test "seqio_index_seqxml_unterminated_entry_raises_with_offset" {
  // <entry without </entry> is a broken record boundary; the error
  // pinpoints the offset of the unterminated record (spec 5.4.3.3)
  let broken = "<?xml version=\"1.0\"?>\n<seqXML>\n  <entry id=\"ok\">\n    <sequence>AAA</sequence>\n  </entry>\n  <entry id=\"broken\">\n    <sequence>CCC</sequence>\n"
  let raised = try {
    ignore(@src.seqio_index(broken, "seqxml"))
    false
  } catch {
    SeqIOError(msg) =>
      msg.contains("unterminated <entry>") && msg.contains("offset")
  }
  assert_true(raised)
}

///|
test "seqio_index_seqxml_self_closing_entry" {
  // a self-closing <entry ... /> record ends on its own start line
  let text = "<?xml version=\"1.0\"?>\n<seqXML>\n  <entry id=\"only\"/>\n</seqXML>\n"
  let index = @src.seqio_index(text, "seqxml")
  assert_eq(index.length(), 1)
  assert_eq(index.keys(), ["only"])
}
