// Whitebox tests for SeqIO functionality.
// These tests run within the package scope, so they can access internal
// functions (@src.parse_fasta, @src.parse_fastq, @src.write_fasta, @src.write_fastq, etc.)
// directly without a package qualifier.
// Run `moon test` to execute.

///|
test "parse_fasta_basic" {
  let records = @src.parse_fasta(">seq1 Hello\nACGT\n>seq2\nTTTT")
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].description, "seq1 Hello")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "TTTT")
}

///|
test "parse_fasta_multiline" {
  let records = @src.parse_fasta(">m1 Multi\nACGT\nTTTT\nGGGG")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "m1")
  assert_eq(records[0].description, "m1 Multi")
  assert_eq(records[0].seq.to_string(), "ACGTTTTTGGGG")
}

///|
test "parse_fasta_no_description" {
  let records = @src.parse_fasta(">seq3\nACGT")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq3")
  assert_eq(records[0].description, "seq3")
  assert_eq(records[0].seq.to_string(), "ACGT")
}

///|
test "parse_fasta_lowercase" {
  let records = @src.parse_fasta(">l1 Lower\nacgtacgt")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "l1")
  assert_eq(records[0].seq.to_string(), "acgtacgt")
}

///|
test "parse_fastq_basic" {
  let records = @src.parse_fastq("@s1 Desc\nACGT\n+\nIIII")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].description, "s1 Desc")
  assert_eq(records[0].seq.to_string(), "ACGT")
  match records[0].phred_quality() {
    Some(q) => assert_true(q == [40, 40, 40, 40])
    None => fail("expected phred_quality to be present")
  }
}

///|
test "parse_fastq_multiple" {
  let records = @src.parse_fastq(
    "@s1 Desc\nACGT\n+\nIIII\n@s2 Other\nTTTT\n+\nJJJJ",
  )
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  match records[0].phred_quality() {
    Some(q) => assert_true(q == [40, 40, 40, 40])
    None => fail("expected phred_quality for s1")
  }
  assert_eq(records[1].id, "s2")
  assert_eq(records[1].seq.to_string(), "TTTT")
  match records[1].phred_quality() {
    Some(q) => assert_true(q == [41, 41, 41, 41])
    None => fail("expected phred_quality for s2")
  }
}

///|
test "write_fasta_basic" {
  let records = [
    @src.SeqRecord::new(
      @src.Seq::new("ACGT"),
      id="r1",
      description="r1 First record",
    ),
    @src.SeqRecord::new(@src.Seq::new("TTTT"), id="r2", description="r2"),
  ]
  let out = @src.write_fasta(records)
  assert_true(out.has_prefix(">"))
  assert_true(out.contains("ACGT"))
  assert_true(out.contains("TTTT"))
  assert_eq(out, ">r1 First record\nACGT\n>r2\nTTTT\n")
}

///|
test "write_fasta_wrapping" {
  // A 70-character sequence should wrap at 60 chars per line, producing
  // two sequence lines of 60 and 10 characters.
  let long_seq = "A".repeat(70)
  let rec = @src.SeqRecord::new(
    @src.Seq::new(long_seq),
    id="long1",
    description="long1",
  )
  let out = @src.write_fasta([rec])
  let expected = ">long1\n" + "A".repeat(60) + "\n" + "A".repeat(10) + "\n"
  assert_eq(out, expected)
  assert_true(out.has_prefix(">long1\n"))
}

///|
test "write_fastq_basic" {
  let rec = @src.SeqRecord::new(
    @src.Seq::new("ACGT"),
    id="fq1",
    description="fq1 record",
  ).set_phred_quality([40, 30, 20, 10])
  let out = @src.write_fastq([rec])
  // Q40 -> 'I', Q30 -> '?', Q20 -> '5', Q10 -> '+'
  assert_eq(out, "@fq1 record\nACGT\n+\nI?5+\n")
}

///|
test "seqio_parse_fasta" {
  let records = @src.seqio_parse(">s1 Hello\nACGT", "fasta")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].description, "s1 Hello")
  assert_eq(records[0].seq.to_string(), "ACGT")
}

///|
test "seqio_parse_fastq" {
  let records = @src.seqio_parse("@s1 Desc\nACGT\n+\nIIII", "fastq")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  match records[0].phred_quality() {
    Some(q) => assert_true(q == [40, 40, 40, 40])
    None => fail("expected phred_quality")
  }
}

///|
test "seqio_parse_unknown_format" {
  let raised = try {
    ignore(@src.seqio_parse(">s1\nACGT", "xyz_unknown"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "seqio_read_single" {
  let rec = @src.seqio_read(">only One\nACGT", "fasta")
  assert_eq(rec.id, "only")
  assert_eq(rec.description, "only One")
  assert_eq(rec.seq.to_string(), "ACGT")
}

///|
test "seqio_read_multiple_error" {
  let raised = try {
    ignore(@src.seqio_read(">s1\nACGT\n>s2\nTTTT", "fasta"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "seqio_read_empty_error" {
  let raised = try {
    ignore(@src.seqio_read("", "fasta"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "seqio_to_dict" {
  let records = @src.parse_fasta(">a First\nACGT\n>b Second\nTTTT")
  let d = @src.seqio_to_dict(records)
  match d.get("a") {
    Some(r) => {
      assert_eq(r.id, "a")
      assert_eq(r.seq.to_string(), "ACGT")
    }
    None => fail("expected key 'a' in dict")
  }
  match d.get("b") {
    Some(r) => assert_eq(r.seq.to_string(), "TTTT")
    None => fail("expected key 'b' in dict")
  }
}

///|
test "seqrecord_new" {
  let rec = @src.SeqRecord::new(
    @src.Seq::new("ACGT"),
    id="myid",
    name="myname",
    description="mydesc",
  )
  assert_eq(rec.id, "myid")
  assert_eq(rec.name, "myname")
  assert_eq(rec.description, "mydesc")
  assert_eq(rec.seq.to_string(), "ACGT")
  assert_eq(rec.length(), 4)
}

///|
test "seqrecord_set_phred_quality" {
  let rec = @src.SeqRecord::from_seq(@src.Seq::new("ACGT"))
  // Initially no quality scores.
  match rec.phred_quality() {
    Some(_) => fail("expected no phred_quality initially")
    None => ()
  }
  // set_phred_quality returns a new record with the annotation.
  let rec2 = rec.set_phred_quality([40, 30, 20, 10])
  match rec2.phred_quality() {
    Some(q) => assert_true(q == [40, 30, 20, 10])
    None => fail("expected phred_quality after set")
  }
}

///|
test "roundtrip_fasta" {
  let original = ">seq1 Hello\nACGT\n>seq2\nTTTT"
  let records = @src.parse_fasta(original)
  let rewritten = @src.write_fasta(records)
  let records2 = @src.parse_fasta(rewritten)
  assert_eq(records2.length(), 2)
  assert_eq(records2[0].id, "seq1")
  assert_eq(records2[0].description, "seq1 Hello")
  assert_eq(records2[0].seq.to_string(), "ACGT")
  assert_eq(records2[1].id, "seq2")
  assert_eq(records2[1].seq.to_string(), "TTTT")
}

///|
test "roundtrip_fastq" {
  let original = "@s1 Desc\nACGT\n+\nIIII"
  let records = @src.parse_fastq(original)
  let rewritten = @src.write_fastq(records)
  let records2 = @src.parse_fastq(rewritten)
  assert_eq(records2.length(), 1)
  assert_eq(records2[0].id, "s1")
  assert_eq(records2[0].description, "s1 Desc")
  assert_eq(records2[0].seq.to_string(), "ACGT")
  match records2[0].phred_quality() {
    Some(q) => assert_true(q == [40, 40, 40, 40])
    None => fail("expected phred_quality after roundtrip")
  }
}
