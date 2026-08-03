///|
/// MEME motif format parsing module tests.

///|
test "meme_motif_new" {
  let pspm = [
    [0.25, 0.25, 0.25, 0.25],
    [0.5, 0.0, 0.5, 0.0],
  ]
  let motif = @src.MemeMotif::new("motif1", pspm)
  assert_eq(motif.name, "motif1")
  assert_eq(motif.alt_name, "")
  assert_eq(motif.alphabet, "ACGT")
  assert_eq(motif.width, 2)
  assert_eq(motif.nsites, 0)
  assert_eq(motif.evalue, 0.0)
  assert_eq(motif.url, "")
  assert_eq(motif.strands, "+ -")
}

///|
test "meme_motif_new_custom_alphabet" {
  let pspm = [[0.5, 0.5]]
  let motif = @src.MemeMotif::new("m1", pspm, alphabet="AC")
  assert_eq(motif.alphabet, "AC")
}

///|
test "meme_motif_consensus" {
  let pspm = [
    [0.9, 0.0, 0.1, 0.0],
    [0.0, 0.8, 0.0, 0.2],
    [0.1, 0.0, 0.9, 0.0],
  ]
  let motif = @src.MemeMotif::new("m1", pspm, alphabet="ACGT")
  let cons = motif.consensus()
  assert_eq(cons, "ACG")
}

///|
test "meme_motif_probability" {
  let pspm = [
    [0.9, 0.0, 0.1, 0.0],
    [0.0, 0.8, 0.0, 0.2],
  ]
  let motif = @src.MemeMotif::new("m1", pspm, alphabet="ACGT")
  assert_eq(motif.probability("A", 0), 0.9)
  assert_eq(motif.probability("C", 0), 0.0)
  assert_eq(motif.probability("G", 0), 0.1)
  assert_eq(motif.probability("C", 1), 0.8)
  assert_eq(motif.probability("T", 1), 0.2)
  // Out of range
  assert_eq(motif.probability("A", 5), 0.0)
}

///|
test "meme_motif_information_content" {
  let pspm = [
    [1.0, 0.0, 0.0, 0.0],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let motif = @src.MemeMotif::new("m1", pspm, alphabet="ACGT")
  let ic = motif.information_content()
  // First position has max info (2 bits), second has 0
  assert_true(ic > 1.9)
  assert_true(ic < 2.1)
}

///|
test "meme_record_new" {
  let record = @src.MemeRecord::new()
  assert_eq(record.version, "4")
  assert_eq(record.alphabet, "ACGT")
  assert_eq(record.strands, "+ -")
  assert_eq(record.background_frequencies.length(), 0)
  assert_eq(record.motifs.length(), 0)
}

///|
test "meme_record_num_motifs" {
  let record = @src.MemeRecord::new()
  assert_eq(record.num_motifs(), 0)
  let pspm = [[0.25, 0.25, 0.25, 0.25]]
  let m = @src.MemeMotif::new("m1", pspm)
  record.motifs.push(m)
  assert_eq(record.num_motifs(), 1)
}

///|
test "meme_record_find" {
  let record = @src.meme_parse(@src.meme_sample_text())
  // Find by name
  match record.find("motif1") {
    Some(found) => assert_eq(found.name, "motif1")
    None => assert_true(false)
  }
  // Find by alt_name
  match record.find("FirstMotif") {
    Some(found) => assert_eq(found.alt_name, "FirstMotif")
    None => assert_true(false)
  }
  // Not found
  match record.find("nonexistent") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "meme_parse_version" {
  let text = "MEME version 5\n"
  let record = @src.meme_parse(text)
  assert_eq(record.version, "5")
}

///|
test "meme_parse_alphabet" {
  let text = "MEME version 4\nALPHABET= ACGT\n"
  let record = @src.meme_parse(text)
  assert_eq(record.alphabet, "ACGT")
}

///|
test "meme_parse_strands" {
  let text = "MEME version 4\nstrands: + -\n"
  let record = @src.meme_parse(text)
  assert_eq(record.strands, "+ -")
}

///|
test "meme_parse_background" {
  let text = "MEME version 4\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"
  let record = @src.meme_parse(text)
  assert_eq(record.background_frequencies.length(), 4)
  assert_eq(record.background_frequencies[0].0, "A")
  assert_eq(record.background_frequencies[0].1, 0.25)
}

///|
test "meme_parse_single_motif" {
  let text = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF motif1 FirstMotif\nletter-probability matrix: alength= 4 w= 3 nsites= 20 E= 1.1e-03\n  0.750000  0.000000  0.000000  0.250000\n  0.000000  0.000000  0.750000  0.250000\n  0.800000  0.000000  0.200000  0.000000\n"
  let record = @src.meme_parse(text)
  assert_eq(record.num_motifs(), 1)
  let m = record.motifs[0]
  assert_eq(m.name, "motif1")
  assert_eq(m.alt_name, "FirstMotif")
  assert_eq(m.width, 3)
  assert_eq(m.nsites, 20)
  assert_true(m.evalue > 0.0)
  assert_eq(m.pspm.length(), 3)
  assert_eq(m.pspm[0].length(), 4)
}

///|
test "meme_parse_multiple_motifs" {
  let record = @src.meme_parse(@src.meme_sample_text())
  assert_eq(record.num_motifs(), 2)
  assert_eq(record.motifs[0].name, "motif1")
  assert_eq(record.motifs[0].alt_name, "FirstMotif")
  assert_eq(record.motifs[1].name, "motif2")
  assert_eq(record.motifs[1].alt_name, "SecondMotif")
}

///|
test "meme_parse_sample_text_motif1" {
  let record = @src.meme_parse(@src.meme_sample_text())
  let m = record.motifs[0]
  assert_eq(m.width, 6)
  assert_eq(m.nsites, 20)
  assert_eq(m.pspm.length(), 6)
  // First row should have 4 values
  assert_eq(m.pspm[0].length(), 4)
  assert_eq(m.pspm[0][0], 0.75)
}

///|
test "meme_parse_sample_text_motif2" {
  let record = @src.meme_parse(@src.meme_sample_text())
  let m = record.motifs[1]
  assert_eq(m.width, 5)
  assert_eq(m.nsites, 15)
  assert_eq(m.pspm.length(), 5)
}

///|
test "meme_sample_motif" {
  let m = @src.meme_sample_motif()
  assert_eq(m.name, "motif1")
  assert_eq(m.alt_name, "FirstMotif")
  assert_eq(m.nsites, 20)
  assert_true(m.evalue > 0.0)
  assert_eq(m.pspm.length(), 6)
}

///|
test "meme_motif_to_string" {
  let m = @src.meme_sample_motif()
  let s = @src.meme_motif_to_string(m)
  assert_true(s.contains("MOTIF motif1 FirstMotif"))
  assert_true(s.contains("letter-probability matrix"))
  assert_true(s.contains("nsites= 20"))
}

///|
test "meme_record_to_string" {
  let record = @src.meme_parse(@src.meme_sample_text())
  let s = @src.meme_record_to_string(record)
  assert_true(s.contains("MEME version 4"))
  assert_true(s.contains("MOTIF motif1"))
}

///|
test "meme_parse_empty" {
  let record = @src.meme_parse("")
  assert_eq(record.num_motifs(), 0)
}

///|
test "meme_parse_url" {
  let text = "MEME version 4\nMOTIF motif1\nletter-probability matrix: alength= 4 w= 2 nsites= 10 E= 0.01\n  0.5 0.0 0.5 0.0\n  0.0 0.5 0.0 0.5\nURL http://example.com/motif1\n"
  let record = @src.meme_parse(text)
  assert_eq(record.num_motifs(), 1)
  assert_eq(record.motifs[0].url, "http://example.com/motif1")
}

///|
test "meme_motif_alphabet_to_array" {
  let pspm = [[0.25, 0.25, 0.25, 0.25]]
  let m = @src.MemeMotif::new("m1", pspm, alphabet="ACGT")
  let arr = m.alphabet_to_array()
  assert_eq(arr.length(), 4)
  assert_eq(arr[0], "A")
  assert_eq(arr[1], "C")
  assert_eq(arr[2], "G")
  assert_eq(arr[3], "T")
}
