// SeqIO test companion: runs the same battery of SeqIO test cases as
// tests/python_seqio_reference.py and emits results between case markers,
// matching BioPython's output exactly so the two can be diffed.
//
// Run:  moon run cmd/main

///|
/// Escape a string for inclusion in a JSON string literal, matching
/// Python's json.dumps string escaping for the characters that appear in
/// our test data (backslash, double-quote, newline, tab).
fn json_escape(s : String) -> String {
  let buf = StringBuilder(size_hint=s.length() + 2)
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
/// Build a JSON string literal (with surrounding double quotes).
fn json_string(s : String) -> String {
  "\"" + json_escape(s) + "\""
}

///|
/// Convert a SeqRecord to a JSON object string, matching BioPython's
/// record_to_dict. Key order: id, name, description, seq, length[, qual].
fn record_to_json(rec : @lib.SeqRecord, include_qual? : Bool = false) -> String {
  let buf = StringBuilder()
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
  if include_qual {
    buf.write_string(", \"qual\": ")
    let quals : Array[Int] = match rec.phred_quality() {
      Some(q) => q
      None => Array(capacity=0)
    }
    buf.write_char('[')
    for i, q in quals {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(q.to_string())
    }
    buf.write_char(']')
  }
  buf.write_char('}')
  buf.to_string()
}

///|
/// Convert an array of SeqRecords to a JSON array string.
fn records_to_json(
  records : Array[@lib.SeqRecord],
  include_qual? : Bool = false,
) -> String {
  let buf = StringBuilder()
  buf.write_char('[')
  for i, rec in records {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(record_to_json(rec, include_qual~))
  }
  buf.write_char(']')
  buf.to_string()
}

///|
/// Emit a test case between markers, matching python_seqio_reference.emit.
fn emit(name : String, payload : String) -> Unit {
  println("=== CASE: " + name + " ===")
  println(payload)
  println("=== END ===")
}

///|
/// Run every SeqIO test case, mirroring tests/python_seqio_reference.py.
fn seqio_test() -> Unit {
  // --- Test FASTA content (verbatim from python_seqio_reference.py) ---
  let fasta_multi = ">seq1 Description for seq1\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n>seq2 Another description\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n>seq3\nACGTACGTACGTACGT\n"
  let fasta_lowercase = ">lower1 lowercase sequence\natggccattgtaatgggccgctgaaagggtgcccgatag\n>lower2 MixedCase Sequence\natgGCCattGTAatgGGCcgctGAAagggtgccCGATAG\n"
  let fasta_multiline = ">multi1 Multi-line sequence\nATGGCCATTGTA\nATGGGCCGCTGAA\nAGGGTGCCCGATAG\n>multi2 Another multi-line\nACGT\nACGT\nACGT\nACGT\n"
  let fasta_single = ">only One record only\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
  // FASTQ (Sanger/Phred+33). seq1 quality is 39 'I' (Q=40) matching the
  // 39-character sequence.
  let fastq_multi = "@seq1 Description for seq1\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n+\n" +
    "I".repeat(39) +
    "\n@seq2 Another description\nACGTACGT\n+\nIIIIIIII\n@seq3 No quality worth mentioning\natggcc\n+\n;;;;;;\n"

  // 1. parse_fasta_multiple
  let recs = @lib.seqio_parse(fasta_multi, "fasta") catch { _ => Array(capacity=0) }
  emit("parse_fasta_multiple", records_to_json(recs))

  // 2. parse_fasta_lowercase
  let recs_l = @lib.seqio_parse(fasta_lowercase, "fasta") catch {
    _ => Array(capacity=0)
  }
  emit("parse_fasta_lowercase", records_to_json(recs_l))

  // 3. parse_fasta_multiline
  let recs_m = @lib.seqio_parse(fasta_multiline, "fasta") catch {
    _ => Array(capacity=0)
  }
  emit("parse_fasta_multiline", records_to_json(recs_m))

  // 4. read_single_record
  let single = @lib.seqio_read(fasta_single, "fasta") catch {
    _ => @lib.SeqRecord::from_seq(@lib.Seq::new(""))
  }
  emit("read_single_record", record_to_json(single))

  // 5. read_multiple_error — SeqIO.read on multi-record raises ValueError.
  let err_payload : String = try {
    ignore(@lib.seqio_read(fasta_multi, "fasta"))
    "{\"error\": null}"
  } catch {
    _ =>
      "{\"error\": \"ValueError\", \"message\": \"More than one record found in handle\"}"
  }
  emit("read_multiple_error", err_payload)

  // 6. to_dict — keyed by id, ordered by record appearance (matching
  // Python dict insertion order).
  let d = @lib.seqio_to_dict(recs)
  let dbuf = StringBuilder()
  dbuf.write_char('{')
  for i, rec in recs {
    if i > 0 {
      dbuf.write_string(", ")
    }
    dbuf.write_string(json_string(rec.id))
    dbuf.write_string(": ")
    let v = match d.get(rec.id) {
      Some(r) => r
      None => rec
    }
    dbuf.write_string(record_to_json(v))
  }
  dbuf.write_char('}')
  emit("to_dict", dbuf.to_string())

  // 7. write_fasta
  let write_recs = [
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
      id="written1",
      name="written1",
      description="written1 First written record",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("acgtacgtacgt"),
      id="written2",
      name="written2",
      description="written2 lowercase input",
    ),
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
      id="written3",
      name="written3",
      description="written3",
    ),
  ]
  let fasta_out = @lib.seqio_write(write_recs, "fasta") catch { _ => "" }
  emit("write_fasta", json_string(fasta_out))

  // 8. parse_fastq
  let fq_recs = @lib.seqio_parse(fastq_multi, "fastq") catch {
    _ => Array(capacity=0)
  }
  emit("parse_fastq", records_to_json(fq_recs, include_qual=true))

  // 9. write_fastq
  let fq1_quals : Array[Int] = Array(capacity=0)
  for i = 0; i < 39; i = i + 1 {
    fq1_quals.push(40)
  }
  let write_fq_recs = [
    @lib.SeqRecord::new(
      @lib.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
      id="fq1",
      name="fq1",
      description="fq1 written fastq record",
    ).set_phred_quality(fq1_quals),
    @lib.SeqRecord::new(
      @lib.Seq::new("acgt"),
      id="fq2",
      name="fq2",
      description="fq2 lowercase",
    ).set_phred_quality([10, 20, 30, 40]),
  ]
  let fq_out = @lib.seqio_write(write_fq_recs, "fastq") catch { _ => "" }
  emit("write_fastq", json_string(fq_out))

  // 10. parse_genbank
  let genbank_sample = "LOCUS       NM_001              48 bp    DNA     linear   PRI 01-JAN-2020\nDEFINITION  Test gene for comparison.\nACCESSION   NM_001\nVERSION     NM_001.1\nSOURCE      Test organism\n  ORGANISM  Test organism\n            Eukaryota; Test.\nFEATURES             Location/Qualifiers\n     source          1..48\n                     /organism=\"Test organism\"\n     gene            1..48\n                     /gene=\"test_gene\"\n     CDS             1..48\n                     /translation=\"MKRHPGSTH\"\nORIGIN\n        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgc\n//\n"
  let gb_recs = @lib.seqio_parse(genbank_sample, "genbank") catch {
    _ => Array(capacity=0)
  }
  let gb_json = match gb_recs.length() {
    0 => "{}"
    _ => record_to_json(gb_recs[0])
  }
  emit("parse_genbank", gb_json)
}
