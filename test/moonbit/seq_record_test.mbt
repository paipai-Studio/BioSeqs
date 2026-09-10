// Whitebox tests for the SeqRecord structure, including the annotations and
// features fields added for the seqio_unify integration (spec 5.1.1 rule 4).

///|
test "seq_record_new_defaults" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1")
  assert_eq(rec.annotations.length(), 0)
  assert_eq(rec.features.length(), 0)
}

///|
test "seq_record_annotations_roundtrip" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1")
  rec.annotations["organism"] = "Homo sapiens"
  rec.annotations["keywords"] = "test; seqio"
  assert_eq(rec.annotations.get("organism"), Some("Homo sapiens"))
  assert_eq(rec.annotations.get("keywords"), Some("test; seqio"))
  assert_eq(rec.annotations.get("missing"), None)
}

///|
test "seq_record_features_roundtrip" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1")
  let f = @src.SeqFeature::new(
    feature_type="CDS",
    id="feat1",
  )
  rec.features.push(f)
  assert_eq(rec.features.length(), 1)
  assert_eq(rec.features[0].feature_type, "CDS")
  assert_eq(rec.features[0].id, "feat1")
}

///|
test "seq_record_from_seq_defaults" {
  let rec = @src.SeqRecord::from_seq(@src.Seq::new("TTTT"))
  assert_eq(rec.annotations.length(), 0)
  assert_eq(rec.features.length(), 0)
}