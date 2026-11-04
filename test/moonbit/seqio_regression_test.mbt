// Regression tests locking down the pre-existing unified SeqIO/AlignIO
// entry-point behaviors (spec 5.7.1 rule 4 / design 2.4.1).
// These tests pin representative outputs of every existing branch of
// seqio_parse / seqio_write / alignio_parse / alignio_write and the exact
// error messages of seqio_read, so that the format-integration work in
// seqio_unify can be verified against a stable behavioral baseline.

// ============================================================================
// seqio_parse: existing 9 format branches
// ============================================================================

///|
test "regression_seqio_parse_fasta" {
  let records = @src.seqio_parse(">s1 Hello\nACGT\n>s2\nTTTT", "fasta")
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].description, "s1 Hello")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].id, "s2")
  assert_eq(records[1].seq.to_string(), "TTTT")
}

///|
test "regression_seqio_parse_fastq_sanger_alias" {
  let fastq = "@s1 Desc\nACGT\n+\nIIII"
  for fmt in ["fastq", "fastq-sanger"] {
    let records = @src.seqio_parse(fastq, fmt)
    assert_eq(records.length(), 1)
    assert_eq(records[0].id, "s1")
    assert_eq(records[0].seq.to_string(), "ACGT")
    assert_eq(records[0].letter_annotations["phred_quality"], [40, 40, 40, 40])
  }
}

///|
test "regression_seqio_parse_fastq_illumina" {
  let records = @src.seqio_parse("@seq1\nACGT\n+\n@ABC\n", "fastq-illumina")
  assert_eq(records.length(), 1)
  assert_eq(records[0].letter_annotations["phred_quality"], [0, 1, 2, 3])
}

///|
test "regression_seqio_parse_fastq_solexa" {
  let records = @src.seqio_parse("@seq1\nACGT\n+\n;>?@\n", "fastq-solexa")
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "ACGT")
}

///|
test "regression_seqio_parse_qual" {
  let records = @src.seqio_parse(">seq1\n10 20 30\n", "qual")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].letter_annotations["phred_quality"], [10, 20, 30])
}

///|
test "regression_seqio_parse_genbank_alias" {
  let gb = "LOCUS       TEST1     4 bp    DNA     linear   UNK 01-JAN-1980\n" +
    "DEFINITION  Test record.\n" +
    "ACCESSION   TEST1\n" +
    "VERSION     TEST1\n" +
    "ORIGIN      \n" +
    "        1 acgt\n" +
    "//\n"
  for fmt in ["genbank", "gb"] {
    let records = @src.seqio_parse(gb, fmt)
    assert_eq(records.length(), 1)
    assert_eq(records[0].seq.to_string(), "ACGT")
  }
}

///|
test "regression_seqio_parse_embl" {
  let embl = "ID   S1; SV 1; linear; DNA; STD; 4 BP.\n" +
    "XX\n" +
    "AC   A001;\n" +
    "XX\n" +
    "DE   A test.\n" +
    "XX\n" +
    "SQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n" +
    "     acgt\n" +
    "//\n"
  let records = @src.seqio_parse(embl, "embl")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "A001")
  assert_eq(records[0].seq.to_string(), "ACGT")
}

///|
test "regression_seqio_parse_pir_alias" {
  let pir = ">P1;P1\nA protein.\nMACE*"
  for fmt in ["pir", "nbrf"] {
    let records = @src.seqio_parse(pir, fmt)
    assert_eq(records.length(), 1)
    assert_eq(records[0].id, "P1")
    assert_eq(records[0].seq.to_string(), "MACE")
  }
}

///|
test "regression_seqio_parse_tab" {
  let records = @src.seqio_parse("s1\tACGT\ns2\tTGCA", "tab")
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "s1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].id, "s2")
  assert_eq(records[1].seq.to_string(), "TGCA")
}

// ============================================================================
// seqio_write: existing 7 format branches
// ============================================================================

///|
test "regression_seqio_write_fasta" {
  let records = @src.seqio_parse(">r1 First\nACGT", "fasta")
  assert_eq(@src.seqio_write(records, "fasta"), ">r1 First\nACGT\n")
}

///|
test "regression_seqio_write_fastq_variants" {
  let rec = @src.SeqRecord::new(
    @src.Seq::new("ACGT"),
    id="fq1",
    description="fq1 record",
  ).set_phred_quality([40, 30, 20, 10])
  assert_eq(@src.seqio_write([rec], "fastq"), "@fq1 record\nACGT\n+\nI?5+\n")
  assert_eq(
    @src.seqio_write([rec], "fastq-sanger"),
    "@fq1 record\nACGT\n+\nI?5+\n",
  )
}

///|
test "regression_seqio_write_qual" {
  let records = @src.seqio_parse(">seq1\n10 20 30\n", "qual")
  let written = @src.seqio_write(records, "qual")
  assert_true(written.contains(">seq1"))
  assert_true(written.contains("10 20 30"))
}

///|
test "regression_seqio_write_genbank_alias" {
  let records = [@src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1")]
  let gb_text = @src.seqio_write(records, "genbank")
  assert_true(gb_text.contains("LOCUS"))
  assert_eq(gb_text, @src.seqio_write(records, "gb"))
}

///|
test "regression_seqio_write_tab" {
  let records = @src.seqio_parse(">r1\nACGT\n>r2\nTTTT", "fasta")
  assert_eq(@src.seqio_write(records, "tab"), "r1\tACGT\nr2\tTTTT\n")
}

// ============================================================================
// alignio_parse / alignio_write: existing 4 format branches
// ============================================================================

///|
test "regression_alignio_parse_fasta" {
  let aligns = @src.alignio_parse(">a\nACGT\n>b\nACGT\n", "fasta")
  assert_eq(aligns.length(), 1)
  assert_eq(aligns[0].records.length(), 2)
  assert_eq(aligns[0].records[0].id, "a")
  assert_eq(aligns[0].records[0].seq.to_string(), "ACGT")
}

///|
test "regression_alignio_parse_clustal" {
  let clustal = "CLUSTAL W (2.1) multiple sequence alignment\n\n" +
    "Tax1        MKT--AILGH\n" +
    "Tax2        M-TQQAILG-\n" +
    "            **   ***  \n\n\n" +
    "Tax1        AA\n" +
    "Tax2        AA\n" +
    "            **\n\n"
  let aligns = @src.alignio_parse(clustal, "clustal")
  assert_eq(aligns.length(), 1)
  assert_eq(aligns[0].records.length(), 2)
}

///|
test "regression_alignio_parse_phylip" {
  let phylip = " 3 4\n" + "Tax1 ACGT\n" + "Tax2 ACGT\n" + "Tax3 ACGT\n"
  let aligns = @src.alignio_parse(phylip, "phylip")
  assert_eq(aligns.length(), 1)
  assert_eq(aligns[0].records.length(), 3)
}

///|
test "regression_alignio_write_fasta" {
  let aligns = @src.alignio_parse(">a\nACGT\n>b\nACGT\n", "fasta")
  let out = @src.alignio_write(aligns, "fasta")
  assert_true(out.contains(">a"))
  assert_true(out.contains("ACGT"))
}

///|
test "regression_alignio_write_clustal" {
  let aligns = @src.alignio_parse(">a\nACGT\n>b\nACGT\n", "fasta")
  let out = @src.alignio_write(aligns, "clustal")
  assert_true(out.contains("a"))
  assert_true(out.contains("ACGT"))
}

///|
test "regression_alignio_write_phylip" {
  let aligns = @src.alignio_parse(">a\nACGT\n>b\nACGT\n", "fasta")
  let out = @src.alignio_write(aligns, "phylip")
  assert_true(out.contains("ACGT"))
  let out2 = @src.alignio_write(aligns, "phylip-sequential")
  assert_true(out2.contains("ACGT"))
}

// ============================================================================
// seqio_read / unknown format: exact error messages
// ============================================================================

///|
test "regression_seqio_read_single_record" {
  let rec = @src.seqio_read(">only One\nACGT", "fasta")
  assert_eq(rec.id, "only")
  assert_eq(rec.seq.to_string(), "ACGT")
}

///|
test "regression_seqio_read_empty_message" {
  let raised = try {
    ignore(@src.seqio_read("", "fasta"))
    false
  } catch {
    SeqIOError(msg) => msg == "No records found in handle"
  }
  assert_true(raised)
}

///|
test "regression_seqio_read_multiple_message" {
  let raised = try {
    ignore(@src.seqio_read(">s1\nACGT\n>s2\nTTTT", "fasta"))
    false
  } catch {
    SeqIOError(msg) => msg == "More than one record found in handle"
  }
  assert_true(raised)
}

///|
test "regression_seqio_parse_unknown_format_message" {
  let raised = try {
    ignore(@src.seqio_parse(">s1\nACGT", "xyz_unknown"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: xyz_unknown"
  }
  assert_true(raised)
}

///|
test "regression_seqio_write_unknown_format_message" {
  let raised = try {
    ignore(@src.seqio_write([], "xyz_unknown"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: xyz_unknown"
  }
  assert_true(raised)
}

///|
test "regression_alignio_unknown_format_message" {
  let raised = try {
    ignore(@src.alignio_parse(">a\nACGT\n", "xyz_unknown"))
    false
  } catch {
    AlignIOError(msg) => msg == "Unknown format: xyz_unknown"
  }
  assert_true(raised)
}

///|
test "regression_seqio_convert_roundtrip" {
  let fasta = ">r1 Desc\nACGT\n"
  let tab = @src.seqio_convert(fasta, "fasta", "tab")
  assert_eq(tab, "r1\tACGT\n")
  // tab format carries no description, so the description is lost on the way back
  let back = @src.seqio_convert(tab, "tab", "fasta")
  assert_eq(back, ">r1\nACGT\n")
}
