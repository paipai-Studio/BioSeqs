///|
test "align_io_parse_clustalw" {
  let clustalw = "CLUSTAL W (1.83) multiple sequence alignment\n\nSeq1        ACGTACGT\nSeq2        ACATACAT\n"
  let align = @src.parse_clustalw(clustalw)

  assert_eq(@src.num_sequences(align), 2)
  assert_eq(@src.alignment_length(align), 8)
}

///|
test "align_io_parse_fasta_alignment" {
  let fasta = ">Seq1\nACGTACGT\n>Seq2\nACATACAT\n"
  let align = @src.parse_fasta_alignment(fasta)

  assert_eq(@src.num_sequences(align), 2)
}

///|
test "align_io_parse_stockholm" {
  let stockholm = "# STOCKHOLM 1.0\nSeq1 ACGTACGT\nSeq2 ACATACAT\n//\n"
  let align = @src.parse_stockholm(stockholm)

  assert_eq(@src.num_sequences(align), 2)
}

///|
test "align_io_write_clustalw" {
  let align = @src.alignment(
    [("Seq1", "ACGT"), ("Seq2", "ACAT")],
    "DNA",
    "Test",
  )
  let output = @src.write_clustalw(align)

  assert_true(output.contains("CLUSTAL"))
}

///|
test "align_io_write_fasta_alignment" {
  let align = @src.alignment(
    [("Seq1", "ACGT"), ("Seq2", "ACAT")],
    "DNA",
    "Test",
  )
  let output = @src.write_fasta_alignment(align)

  assert_true(output.contains(">Seq1"))
}

///|
test "align_io_get_sequence" {
  let align = @src.alignment(
    [("Seq1", "ACGT"), ("Seq2", "ACAT")],
    "DNA",
    "Test",
  )

  let seq = @src.get_sequence(align, "Seq1")
  assert_eq(seq.is_some(), true)
}

///|
test "align_io_consensus_sequence" {
  let align = @src.alignment(
    [("Seq1", "ACGT"), ("Seq2", "ACAT")],
    "DNA",
    "Test",
  )
  let consensus = @src.consensus_sequence(align)

  assert_eq(consensus.length(), 4)
}

///|
test "align_io_column_consensus" {
  let align = @src.alignment(
    [("Seq1", "ACGT"), ("Seq2", "ACAT")],
    "DNA",
    "Test",
  )
  let cons = @src.column_consensus(align, 0)

  assert_eq(cons, "A")
}
