// Tests for the seqio_binary_unify adapters: the unified seqio_parse
// entry must decode base64-carried binary payloads with strict validation
// (spec 4.3.1 / 5.1.3.2), reject bad magic numbers (spec 5.1.3.3), and
// agree with the standalone decoders for every supported binary format
// (spec 5.1.1.2a).

///|
test "base64_decode_strict_basic" {
  let bytes = @src.base64_decode_strict("QUJD")
  assert_eq(bytes.length(), 3)
  assert_eq(bytes[0], 65)
  assert_eq(bytes[1], 66)
  assert_eq(bytes[2], 67)
  // padded and unpadded forms decode identically
  let padded = @src.base64_decode_strict("QQ==")
  assert_eq(padded.length(), 1)
  assert_eq(padded[0], 65)
  let unpadded = @src.base64_decode_strict("QQ")
  assert_eq(unpadded[0], 65)
  // whitespace between groups is skipped
  let spaced = @src.base64_decode_strict("QUJD\n")
  assert_eq(spaced.length(), 3)
}

///|
test "base64_decode_strict_rejects" {
  // lone trailing character
  let lone = try {
    ignore(@src.base64_decode_strict("Q"))
    false
  } catch {
    SeqIOError(msg) => msg == "invalid base64 content"
  }
  assert_true(lone)
  // misplaced padding
  let misplaced = try {
    ignore(@src.base64_decode_strict("QU=D"))
    false
  } catch {
    SeqIOError(msg) => msg == "invalid base64 padding"
  }
  assert_true(misplaced)
  // invalid characters
  let invalid = try {
    ignore(@src.base64_decode_strict("@@@@"))
    false
  } catch {
    SeqIOError(msg) => msg == "invalid base64 content"
  }
  assert_true(invalid)
}

///|
test "parse_twobit_records_roundtrip" {
  let bytes = @src.twobit_to_bytes(@src.twobit_sample_file())
  let content = @src.bytes_to_base64(bytes)
  let records = @src.parse_twobit_records(content)
  let direct = @src.twobit_to_seq_records(@src.twobit_sample_file())
  assert_eq(records.length(), direct.length())
  assert_true(records.length() > 0)
  for i in 0..<records.length() {
    assert_eq(records[i].id, direct[i].id)
    assert_eq(records[i].seq.to_string(), direct[i].seq.to_string())
  }
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "2bit")
  assert_eq(via_entry.length(), records.length())
  assert_eq(via_entry[0].id, records[0].id)
}

///|
test "parse_twobit_records_bad_magic" {
  let raised = try {
    ignore(@src.parse_twobit_records(@src.bytes_to_base64([1, 2, 3, 4])))
    false
  } catch {
    SeqIOError(msg) => msg == "2bit: invalid magic number (not a 2bit file)"
  }
  assert_true(raised)
}

///|
test "parse_sff_records_roundtrip" {
  let content = @src.bytes_to_base64(@src.sff_encode(@src.sff_create_example_file()))
  let records = @src.parse_sff_records(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "READ0001")
  assert_eq(records[0].seq.to_string(), "ACGTACGTAC")
  assert_eq(records[1].id, "READ0002")
  assert_eq(records[1].seq.to_string(), "TGCATGCATG")
  match records[0].phred_quality() {
    Some(quals) => {
      assert_eq(quals.length(), 10)
      assert_eq(quals[0], 10)
    }
    None => assert_true(false)
  }
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "sff")
  assert_eq(via_entry.length(), 2)
  assert_eq(via_entry[0].id, records[0].id)
}

///|
test "parse_sff_trim_records_clip_semantics" {
  // clip quality interval [3, 8] on "ACGTACGTAC" keeps bases[2:8]
  let qualities = [10, 15, 20, 25, 30, 28, 26, 24, 22, 20]
  let read = @src.SffRead::new(
    "READTRIM",
    "ACGTACGTAC",
    qualities,
    [1.0, 0.0, 1.0, 0.0],
    [1, 2, 3, 4, 1, 2, 3, 4, 1, 2],
    3,
    8,
    0,
    0,
  )
  let header = @src.SffHeader::new(1, 4, "TACG", "TCAG")
  let content = @src.bytes_to_base64(@src.sff_encode(@src.SffFile::new(header, [read])))
  // untrimmed output keeps the full read
  let full = @src.parse_sff_records(content)
  assert_eq(full[0].seq.to_string(), "ACGTACGTAC")
  // trimmed output restricts sequence and qualities to the clip interval
  let trimmed = @src.parse_sff_trim_records(content)
  assert_eq(trimmed.length(), 1)
  assert_eq(trimmed[0].id, "READTRIM")
  assert_eq(trimmed[0].seq.to_string(), "GTACGT")
  match trimmed[0].phred_quality() {
    Some(quals) => {
      assert_eq(quals.length(), 6)
      assert_eq(quals[0], 20)
      assert_eq(quals[5], 24)
    }
    None => assert_true(false)
  }
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "sff-trim")
  assert_eq(via_entry[0].seq.to_string(), "GTACGT")
}

///|
test "parse_sff_records_malformed" {
  // wrong magic number makes the underlying parser return None
  let bytes = [
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 40,
    0, 4, 0, 4, 0, 1,
  ]
  let raised = try {
    ignore(@src.parse_sff_records(@src.bytes_to_base64(bytes)))
    false
  } catch {
    SeqIOError(msg) => msg == "Malformed SFF data"
  }
  assert_true(raised)
}

///|
test "parse_snapgene_records_roundtrip" {
  let content = @src.bytes_to_base64(@src.snapgene_to_bytes(@src.snapgene_sample()))
  let records = @src.parse_snapgene_records(content)
  assert_eq(records.length(), 1)
  let direct = @src.snapgene_to_seqrecord(@src.snapgene_sample())
  assert_eq(records[0].id, direct.id)
  assert_eq(records[0].seq.to_string(), direct.seq.to_string())
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "snapgene")
  assert_eq(via_entry.length(), 1)
  assert_eq(via_entry[0].seq.to_string(), records[0].seq.to_string())
}

///|
test "parse_gck_records_roundtrip_and_magic" {
  let content = @src.bytes_to_base64(@src.gck_to_bytes(@src.gck_sample()))
  let records = @src.parse_gck_records(content)
  assert_eq(records.length(), 1)
  let direct = @src.gck_to_seqrecord(@src.gck_sample())
  assert_eq(records[0].id, direct.id)
  assert_eq(records[0].seq.to_string(), direct.seq.to_string())
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "gck")
  assert_eq(via_entry[0].seq.to_string(), records[0].seq.to_string())
  // wrong magic number raises (spec 5.1.3.3)
  let bad = [
    71, 67, 75, 95, 70, 73, 76, 88, 0, 0, 0, 0, 0, 0,
  ]
  let raised = try {
    ignore(@src.parse_gck_records(@src.bytes_to_base64(bad)))
    false
  } catch {
    SeqIOError(msg) => msg == "gck: invalid magic number (not a GCK file)"
  }
  assert_true(raised)
}

///|
test "parse_xdna_records_roundtrip" {
  let content = @src.bytes_to_base64(@src.xdna_to_bytes(@src.xdna_sample_file()))
  let records = @src.parse_xdna_records(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ATGCATGC")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "GATTACA")
  assert_eq(records[0].annotations["sequence_type"], "DNA")
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "xdna")
  assert_eq(via_entry.length(), 2)
  assert_eq(via_entry[0].seq.to_string(), records[0].seq.to_string())
}

///|
test "parse_nib_records_single_record" {
  let content = @src.bytes_to_base64(@src.nib_encode("ATGCATGCATGC"))
  let records = @src.parse_nib_records(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "nib")
  assert_eq(records[0].seq.to_string(), "ATGCATGCATGC")
  // unified entry agrees
  let via_entry = @src.seqio_parse(content, "nib")
  assert_eq(via_entry.length(), 1)
  assert_eq(via_entry[0].seq.to_string(), "ATGCATGCATGC")
}

///|
test "parse_binary_empty_content_returns_empty" {
  assert_eq(@src.parse_twobit_records("").length(), 0)
  assert_eq(@src.parse_sff_records("").length(), 0)
  assert_eq(@src.parse_sff_trim_records("").length(), 0)
  assert_eq(@src.parse_snapgene_records("").length(), 0)
  assert_eq(@src.parse_gck_records("").length(), 0)
  assert_eq(@src.parse_xdna_records("").length(), 0)
  assert_eq(@src.parse_nib_records("").length(), 0)
}

///|
test "parse_binary_invalid_base64_names_format" {
  // decode failures terminate before entering the format parser and
  // name the format in the message (spec 5.1.3.2)
  let sff = try {
    ignore(@src.parse_sff_records("@@@"))
    false
  } catch {
    SeqIOError(msg) => msg == "invalid base64 content for sff"
  }
  assert_true(sff)
  let nib = try {
    ignore(@src.parse_nib_records("not base64!!"))
    false
  } catch {
    SeqIOError(msg) => msg.contains("for nib")
  }
  assert_true(nib)
  let gck = try {
    ignore(@src.parse_gck_records("a"))
    false
  } catch {
    SeqIOError(msg) => msg.contains("for gck")
  }
  assert_true(gck)
}
