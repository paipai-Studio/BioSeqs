///|
fn json_escape(s : String) -> String {
  let buf = StringBuilder::new(size_hint=s.length() + 2)
  for c in s.iter() {
    match c {
      '\\' => buf.write_string("\\\\")
      '"' => buf.write_string("\\\"")
      '\n' => buf.write_string("\\n")
      '\t' => buf.write_string("\\t")
      _ => buf.write_char(c)
    }
  }
  buf.to_string()
}

///|
fn json_string(s : String) -> String {
  "\"" + json_escape(s) + "\""
}

///|
fn record_to_json(rec : @lib.SeqRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string("{\"id\": ")
  buf.write_string(json_string(rec.id))
  buf.write_string(", \"name\": ")
  buf.write_string(json_string(rec.name))
  buf.write_string(", \"description\": ")
  buf.write_string(json_string(rec.description))
  buf.write_string(", \"seq\": ")
  buf.write_string(json_string(rec.seq.to_string()))
  buf.write_string(", \"length\": ")
  buf.write_string(rec.length().to_string())
  buf.write_char('}')
  buf.to_string()
}

///|
fn records_to_json(records : Array[@lib.SeqRecord]) -> String {
  let buf = StringBuilder::new()
  buf.write_char('[')
  for i, rec in records {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(record_to_json(rec))
  }
  buf.write_char(']')
  buf.to_string()
}

///|
fn alignment_to_json(align : @lib.MultipleSeqAlignment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("{\"num_records\": ")
  buf.write_string(align.num_records().to_string())
  buf.write_string(", \"alignment_length\": ")
  buf.write_string(align.get_alignment_length().to_string())
  buf.write_string(", \"records\": ")
  buf.write_string(records_to_json(align.iter()))
  buf.write_char('}')
  buf.to_string()
}

///|
fn emit(name : String, payload : String) -> Unit {
  println("=== CASE: " + name + " ===")
  println(payload)
  println("=== END ===")
}

///|
fn run_alignio_test() -> Unit raise @lib.AlignIOError {
  let clustal_sample = "CLUSTAL X (1.81) multiple sequence alignment\n\nAF191659.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG\nAF191658.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG\nAF191661.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG\nAF191660.1      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG\n\nAF191659.1      AGA\nAF191658.1      AGA\nAF191661.1      AGA\nAF191660.1      AGA\n"

  let phylip_sample = "3 20\nseq1        ATGCCGATCGATCGATCGAT\nseq2        ATGCCGATCGATCGATCGAT\nseq3        ATGCCGATCGATCGATCGAT\n"

  let phylip_seq_sample = "3 20\nseq1        ATGCCGATCGATCGATCGAT\nseq2        ATGCCGATCGATCGATCGAT\nseq3        ATGCCGATCGATCGATCGAT\n"

  let fasta_align_sample = ">seq1 Description 1\nATGCCGATCGATCGATCGATCGATCGATCGATCGATC\n>seq2 Description 2\nATGCCGATCGATCGATCGATCGATCGATCGATCGATC\n>seq3 Description 3\nATGCCGATCGATCGATCGATCGATCGATCGATCGATC\n"

  let clustal_align = @lib.alignio_read(clustal_sample, "clustal")
  emit("parse_clustal", alignment_to_json(clustal_align))

  let phylip_align = @lib.alignio_read(phylip_sample, "phylip")
  emit("parse_phylip", alignment_to_json(phylip_align))

  let phylip_seq_align = @lib.alignio_read(
    phylip_seq_sample, "phylip-sequential",
  )
  emit("parse_phylip_sequential", alignment_to_json(phylip_seq_align))

  let fasta_align = @lib.alignio_read(fasta_align_sample, "fasta")
  emit("parse_fasta_align", alignment_to_json(fasta_align))

  let write_clustal_recs = [
    @lib.SeqRecord::new(
      @lib.Seq::new("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
      id="AF191659.1",
      name="AF191659.1",
      description="AF191659.1",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
      id="AF191658.1",
      name="AF191658.1",
      description="AF191658.1",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
      id="AF191660.1",
      name="AF191660.1",
      description="AF191660.1",
    ),
  ]
  let write_clustal_align = @lib.MultipleSeqAlignment::new(write_clustal_recs)
  let clustal_out = @lib.alignio_write([write_clustal_align], "clustal")
  emit("write_clustal", json_string(clustal_out))

  let write_phylip_recs = [
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGCCGATCGATCGATCGAT"),
      id="seq1",
      name="seq1",
      description="seq1",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGCCGATCGATCGATCGAT"),
      id="seq2",
      name="seq2",
      description="seq2",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGCCGATCGATCGATCGAT"),
      id="seq3",
      name="seq3",
      description="seq3",
    ),
  ]
  let write_phylip_align = @lib.MultipleSeqAlignment::new(write_phylip_recs)
  let phylip_out = @lib.alignio_write([write_phylip_align], "phylip")
  emit("write_phylip", json_string(phylip_out))

  let write_fasta_recs = [
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGCCGATCGATCGATCGATCGATCGATCGATCGATC"),
      id="seq1",
      name="seq1",
      description="seq1 desc",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGCCGATCGATCGATCGATCGATCGATCGATCGATC"),
      id="seq2",
      name="seq2",
      description="seq2 desc",
    ),
  ]
  let write_fasta_align = @lib.MultipleSeqAlignment::new(write_fasta_recs)
  let fasta_out = @lib.alignio_write([write_fasta_align], "fasta")
  emit("write_fasta_align", json_string(fasta_out))

  emit(
    "alignment_length",
    "{\"num_records\": " +
    clustal_align.num_records().to_string() +
    ", \"alignment_length\": " +
    clustal_align.get_alignment_length().to_string() +
    "}",
  )

  let col0 = match clustal_align.column(0) {
    Some(c) => c
    None => ""
  }
  let col5 = match clustal_align.column(5) {
    Some(c) => c
    None => ""
  }
  emit(
    "alignment_column",
    "{\"column_0\": " +
    json_string(col0) +
    ", \"column_5\": " +
    json_string(col5) +
    "}",
  )

  let err_payload : String = try {
    ignore(@lib.alignio_read(clustal_sample + "\n" + clustal_sample, "clustal"))
    "{\"error\": null}"
  } catch {
    _ =>
      "{\"error\": \"ValueError\", \"message\": \"More than one alignment found in handle\"}"
  }
  emit("read_multiple_error", err_payload)
}

///|
fn alignio_test() -> Unit {
  let _ = try {
    run_alignio_test()
    ()
  } catch {
    _ => println("Test error")
  }
}
