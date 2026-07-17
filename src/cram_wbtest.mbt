///|
/// CRAM format parsing white-box tests.

///|
test "cram_magic" {
  assert_true(is_cram_magic("CRAMtest"))
  assert_false(is_cram_magic("BAMtest"))
}

///|
test "cram_record_is_paired" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let record = CramRecord::{
    seq_id: 0,
    alignment_start: 100,
    read_length: 100,
    read_group: 0,
    flags: 1,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }
  assert_true(record.is_paired())
}

///|
test "cram_record_is_mapped" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let record = CramRecord::{
    seq_id: 0,
    alignment_start: 100,
    read_length: 100,
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }
  assert_true(record.is_mapped())
}

///|
test "cram_record_is_unmapped" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let record = CramRecord::{
    seq_id: -1,
    alignment_start: 0,
    read_length: 100,
    read_group: 0,
    flags: 4,
    mapping_quality: 0,
    mate_seq_id: 0,
    mate_alignment_start: 0,
    template_size: 0,
    num_reads: 2,
    read_number: 1,
    fragment_length: 0,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }
  assert_false(record.is_mapped())
}

///|
test "cram_record_is_reverse" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let record = CramRecord::{
    seq_id: 0,
    alignment_start: 100,
    read_length: 100,
    read_group: 0,
    flags: 16,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }
  assert_true(record.is_reverse())
}

///|
test "cram_file_num_references" {
  let header = BamHeader::{
    text: "",
    references: FixedArray::make(3, BamReference::{ name: "", length: 0 }),
  }
  let file = CramFile::{ header, records: Array::new() }
  assert_eq(file.num_references(), 3)
}

///|
test "cram_record_get_reference_name" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let references : FixedArray[BamReference] = FixedArray::make(2, BamReference::{
    name: "",
    length: 0,
  })
  references[0] = BamReference::{ name: "chr1", length: 1000 }
  references[1] = BamReference::{ name: "chr2", length: 2000 }
  let header = BamHeader::{ text: "", references }

  let record = CramRecord::{
    seq_id: 0,
    alignment_start: 100,
    read_length: 100,
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }

  assert_eq(record.get_reference_name(header), "chr1")
}

///|
test "cram_record_get_reference_name_unknown" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(0, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let references : FixedArray[BamReference] = FixedArray::make(2, BamReference::{
    name: "",
    length: 0,
  })
  let header = BamHeader::{ text: "", references }

  let record = CramRecord::{
    seq_id: -1,
    alignment_start: 100,
    read_length: 100,
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }

  assert_eq(record.get_reference_name(header), "*")
}

///|
test "cram_record_get_end" {
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(2, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  cigar[0] = BamCigarElement::{ op: BamCigarOp::Match, length: 50 }
  cigar[1] = BamCigarElement::{ op: BamCigarOp::Insert, length: 10 }
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    0,
    CramReadFeature::UnknownFeature(0),
  )

  let record = CramRecord::{
    seq_id: 0,
    alignment_start: 100,
    read_length: 60,
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    mate_seq_id: 0,
    mate_alignment_start: 200,
    template_size: 100,
    num_reads: 2,
    read_number: 1,
    fragment_length: 100,
    qname: "read1",
    cigar,
    seq: "ACGT",
    qual: "!!!!",
    tags: Array::new(),
    read_features: features,
  }

  assert_eq(record.get_end(), 149)
}

///|
test "read_itf8_positive" {
  let data = "\x02"
  let (val, pos) = read_itf8(data, 0)
  assert_eq(val, 1)
  assert_eq(pos, 1)
}

///|
test "read_itf8_negative" {
  let data = "\x01"
  let (val, pos) = read_itf8(data, 0)
  assert_eq(val, -1)
  assert_eq(pos, 1)
}

///|
test "read_ltf8" {
  let data = "\x02"
  let (val, pos) = read_ltf8(data, 0)
  assert_eq(val, 1)
  assert_eq(pos, 1)
}
