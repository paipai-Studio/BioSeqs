///|
/// Tests for @src.seq_standardize.mbt — U+T standardization pipeline.

///|
fn seq_std_record(s : String) -> @src.SeqRecord {
  @src.SeqRecord::new(@src.Seq::new(s), id="test")
}

///|
test "seq_standardize_uppercase_and_strip" {
  let opts = @src.seq_standardize_defaults()
  let rec = seq_std_record("acgtn-acg*")
  match @src.seq_standardize(rec, opts) {
    Ok(pair) => {
      let (clean, report) = pair
      assert_eq(clean.seq.to_string(), "ACGTNACG")
      assert_eq(report.removed_count, 2)
      let mut has_upper = false
      let mut has_strip = false
      for a in report.actions {
        if a == "uppercase" {
          has_upper = true
        }
        if a.contains("stripped") {
          has_strip = true
        }
      }
      assert_true(has_upper)
      assert_true(has_strip)
    }
    Err(_) => fail("simple cleanup should succeed")
  }
}

///|
test "seq_standardize_alphabet_detection" {
  // RNA: U present, no T
  assert_eq(@src.seq_std_detect_alphabet("ACGU"), "rna")
  // protein-specific residues
  assert_eq(@src.seq_std_detect_alphabet("MKTAYIAK"), "protein")
  // default DNA
  assert_eq(@src.seq_std_detect_alphabet("ACGT"), "dna")
}

///|
test "seq_standardize_keep_ambiguity" {
  let opts = @src.seq_standardize_defaults() // Keep policy
  let rec = seq_std_record("acgtN")
  match @src.seq_standardize(rec, opts) {
    Ok(pair) => {
      let (clean, _) = pair
      assert_eq(clean.seq.to_string(), "ACGTN")
    }
    Err(_) => fail("Keep policy must accept ambiguity codes")
  }
}

///|
test "seq_standardize_reject_ambiguity" {
  let opts = @src.SeqStandardizeOptions::new(
    ambiguity_policy=@src.Reject,
    target_alphabet="dna",
  )
  let rec = seq_std_record("ACGTN")
  match @src.seq_standardize(rec, opts) {
    Ok(_) => fail("Reject policy must refuse ambiguity codes")
    Err(report) => {
      assert_eq(report.invalid_chars.length(), 1)
      assert_eq(report.invalid_chars[0].1, "N")
    }
  }
}

///|
test "seq_standardize_invalid_chars_fail" {
  let opts = @src.seq_standardize_defaults()
  let rec = seq_std_record("ACGT1AC")
  match @src.seq_standardize(rec, opts) {
    Ok(_) => fail("digits are not a valid nucleotide")
    Err(report) => assert_eq(report.invalid_chars[0].1, "1")
  }
}

///|
test "seq_standardize_batch_splits" {
  let opts = @src.seq_standardize_defaults()
  let ok1 = seq_std_record("ACGT")
  let ok2 = seq_std_record("acgt")
  let bad = seq_std_record("ACGT!!!")
  let (accepted, rejected) = @src.seq_standardize_batch([ok1, ok2, bad], opts)
  assert_eq(accepted.length(), 2)
  assert_eq(rejected.length(), 1)
  assert_eq(rejected[0].0.seq.to_string(), "ACGT!!!")
}
