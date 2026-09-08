///|
/// Tests for Bio.SeqIO.index: byte-offset indexes for random access to
/// individual records without re-parsing the whole file.

///|
test "seqio_index_fasta_keys_and_raw" {
  let fasta =
    ">seq1 first record\nACGTACGT\nTTTT\n>seq2 second\nNNNNACGT\n"
  let index = try {
    @src.seqio_index(fasta, "fasta")
  } catch {
    SeqIOError(msg) => {
      println(msg)
      @src.seqio_index("", "fasta")
    }
  }
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["seq1", "seq2"])
  assert_true(index.contains("seq1"))
  assert_true(!index.contains("seq3"))

  // get_raw returns the exact record slice.
  let raw2 = try {
    index.get_raw(fasta, "seq2")
  } catch {
    SeqIOError(_) => ""
  }
  assert_eq(raw2, ">seq2 second\nNNNNACGT\n")

  // get parses the slice into a record.
  let record = try {
    index.get(fasta, "seq1")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.id, "seq1")
  assert_eq(record.seq.to_string(), "ACGTACGTTTTT")

  // Missing key raises.
  let missing = try {
    let _ = index.get_raw(fasta, "nope")
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(missing)
}

///|
test "seqio_index_fasta_no_trailing_newline" {
  let fasta = ">a\nACGT\n>b\nTTTT"
  let index = try {
    @src.seqio_index(fasta, "fasta")
  } catch {
    SeqIOError(_) => @src.seqio_index("", "fasta")
  }
  assert_eq(index.length(), 2)
  let raw_b = try {
    index.get_raw(fasta, "b")
  } catch {
    SeqIOError(_) => ""
  }
  assert_eq(raw_b, ">b\nTTTT")
  let record = try {
    index.get(fasta, "b")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.seq.to_string(), "TTTT")
}

///|
test "seqio_index_fastq" {
  let fastq =
    "@read1 comment\nACGT\n+\nHHHH\n@read2\nTTTT\n+\n!!!!\n"
  let index = try {
    @src.seqio_index(fastq, "fastq")
  } catch {
    SeqIOError(msg) => {
      println(msg)
      @src.seqio_index("", "fasta")
    }
  }
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["read1", "read2"])

  let raw1 = try {
    index.get_raw(fastq, "read1")
  } catch {
    SeqIOError(_) => ""
  }
  assert_eq(raw1, "@read1 comment\nACGT\n+\nHHHH\n")

  let record = try {
    index.get(fastq, "read2")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.id, "read2")
  assert_eq(record.seq.to_string(), "TTTT")
}

///|
test "seqio_index_fastq_malformed" {
  // Missing the '+' separator line.
  let bad = "@read1\nACGT\nHHHH\n@read2\nTTTT\n"
  let raised = try {
    let _ = @src.seqio_index(bad, "fastq")
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "seqio_index_tab" {
  let tab = "seq1\tACGTACGT\nseq2\tNNNN\n"
  let index = try {
    @src.seqio_index(tab, "tab")
  } catch {
    SeqIOError(_) => @src.seqio_index("", "fasta")
  }
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["seq1", "seq2"])
  let raw1 = try {
    index.get_raw(tab, "seq1")
  } catch {
    SeqIOError(_) => ""
  }
  assert_eq(raw1, "seq1\tACGTACGT\n")
  let record = try {
    index.get(tab, "seq1")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.seq.to_string(), "ACGTACGT")
}

///|
test "seqio_index_embl" {
  let embl =
    "ID   SEQ1; SV 1; linear; DNA; STD; 10 BP.\n" +
    "XX\n" +
    "AC   ACC001;\n" +
    "XX\n" +
    "DE   First sequence.\n" +
    "XX\n" +
    "SQ   Sequence 10 BP; 2 A; 3 C; 2 G; 3 T; 0 other;\n" +
    "     aaccggtttt\n" +
    "//\n" +
    "ID   SEQ2; SV 1; linear; DNA; STD; 8 BP.\n" +
    "XX\n" +
    "AC   ACC002;\n" +
    "XX\n" +
    "DE   Second sequence description.\n" +
    "XX\n" +
    "SQ   Sequence 8 BP; 2 A; 2 C; 2 G; 2 T; 0 other;\n" +
    "     acgtacgt\n" +
    "//\n"
  let index = try {
    @src.seqio_index(embl, "embl")
  } catch {
    SeqIOError(msg) => {
      println(msg)
      @src.seqio_index("", "fasta")
    }
  }
  assert_eq(index.length(), 2)
  // Keys come from the ID line entry name.
  assert_eq(index.keys(), ["SEQ1", "SEQ2"])
  let raw1 = try {
    index.get_raw(embl, "SEQ1")
  } catch {
    SeqIOError(_) => ""
  }
  assert_true(raw1.has_prefix("ID   SEQ1"))
  assert_true(raw1.contains("//\n"))
  // The last record runs to EOF.
  let raw2 = try {
    index.get_raw(embl, "SEQ2")
  } catch {
    SeqIOError(_) => ""
  }
  assert_true(raw2.has_prefix("ID   SEQ2"))
  let record = try {
    index.get(embl, "SEQ2")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.seq.to_string(), "ACGTACGT")
}

///|
test "seqio_index_genbank" {
  let gb =
    "LOCUS       NC_001     10 bp    DNA     linear\n" +
    "FEATURES             Location/Qualifiers\n" +
    "     source          1..10\n" +
    "ORIGIN\n" +
    "       1 aaccggtttt\n" +
    "//\n" +
    "LOCUS       NC_002      8 bp    DNA     linear\n" +
    "FEATURES             Location/Qualifiers\n" +
    "     source          1..8\n" +
    "ORIGIN\n" +
    "       1 acgtacgt\n" +
    "//\n"
  let index = try {
    @src.seqio_index(gb, "genbank")
  } catch {
    SeqIOError(msg) => {
      println(msg)
      @src.seqio_index("", "fasta")
    }
  }
  assert_eq(index.length(), 2)
  assert_eq(index.keys(), ["NC_001", "NC_002"])
  let record = try {
    index.get(gb, "NC_001")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record.seq.to_string(), "AACCGGTTTT")
  let record2 = try {
    index.get(gb, "NC_002")
  } catch {
    SeqIOError(_) => @src.seqio_read(">x\nA", "fasta")
  }
  assert_eq(record2.seq.to_string(), "ACGTACGT")
}

///|
test "seqio_index_unsupported_format" {
  let raised = try {
    let _ = @src.seqio_index(">a\nACGT\n", "clustal")
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}
