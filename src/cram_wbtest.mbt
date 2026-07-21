///|
/// CRAM format parsing white-box tests.

///|
test "cram_magic" {
  assert_true(is_cram_magic("CRAMtest"))
  assert_false(is_cram_magic("BAMtest"))
}

///|
test "cram_record_is_paired" {
  let cigar : Array[BamCigarElement] = Array::new()
  cigar.push(BamCigarElement::{ op: BamCigarOp::Match, length: 100 })
  
  let read_features = CramReadFeatures::{
    read_group: 0,
    flags: 1,
    mapping_quality: 60,
    read_name: "read1",
    cigar,
    sequence: "ACGT",
    quality: "!!!!",
    tags: Array::new(),
  }
  
  let record = CramRecord::{ ref_id: 0, pos: 100, read_features }
  assert_true((record.read_features.flags & 1) != 0)
}

///|
test "cram_record_is_mapped" {
  let cigar : Array[BamCigarElement] = Array::new()
  
  let read_features = CramReadFeatures::{
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    read_name: "read1",
    cigar,
    sequence: "ACGT",
    quality: "!!!!",
    tags: Array::new(),
  }
  
  let record = CramRecord::{ ref_id: 0, pos: 100, read_features }
  assert_true((record.read_features.flags & 4) == 0)
}

///|
test "cram_record_is_unmapped" {
  let cigar : Array[BamCigarElement] = Array::new()
  
  let read_features = CramReadFeatures::{
    read_group: 0,
    flags: 4,
    mapping_quality: 0,
    read_name: "read1",
    cigar,
    sequence: "ACGT",
    quality: "!!!!",
    tags: Array::new(),
  }
  
  let record = CramRecord::{ ref_id: -1, pos: 0, read_features }
  assert_true((record.read_features.flags & 4) != 0)
}

///|
test "cram_record_is_reverse" {
  let cigar : Array[BamCigarElement] = Array::new()
  
  let read_features = CramReadFeatures::{
    read_group: 0,
    flags: 16,
    mapping_quality: 60,
    read_name: "read1",
    cigar,
    sequence: "ACGT",
    quality: "!!!!",
    tags: Array::new(),
  }
  
  let record = CramRecord::{ ref_id: 0, pos: 100, read_features }
  assert_true((record.read_features.flags & 16) != 0)
}

///|
test "cram_file_num_references" {
  let references : Array[CramReference] = Array::new()
  references.push(CramReference::{ name: "chr1", length: 1000, md5: "", uri: "" })
  references.push(CramReference::{ name: "chr2", length: 2000, md5: "", uri: "" })
  references.push(CramReference::{ name: "chr3", length: 3000, md5: "", uri: "" })
  
  let header = CramHeader::{
    version: CramVersion::{ major: 3, minor: 0 },
    reference_count: 3,
    references,
    compression: CramCompression::Gzip,
    stream_count: 0,
    stream_info: Array::new(),
    sam_header: "",
  }
  let file = CramFile::{ header, records: Array::new() }
  assert_eq(file.header.references.length(), 3)
}

///|
test "cram_record_to_bam_record" {
  let cigar : Array[BamCigarElement] = Array::new()
  cigar.push(BamCigarElement::{ op: BamCigarOp::Match, length: 100 })
  
  let read_features = CramReadFeatures::{
    read_group: 0,
    flags: 0,
    mapping_quality: 60,
    read_name: "read1",
    cigar,
    sequence: "ACGT",
    quality: "!!!!",
    tags: Array::new(),
  }
  
  let record = CramRecord::{ ref_id: 0, pos: 100, read_features }
  let bam_record = record.to_bam_record()
  
  assert_eq(bam_record.ref_id, 0)
  assert_eq(bam_record.pos, 100)
  assert_eq(bam_record.mapq, 60)
  assert_eq(bam_record.qname, "read1")
}

///|
fn read_itf8(data : String, pos : Int) -> (Int, Int) {
  let byte = data[pos].to_int() & 0xFF
  let value = byte >> 1
  if (byte & 0x01) == 0 {
    (value, pos + 1)
  } else {
    (-value, pos + 1)
  }
}

///|
fn read_ltf8(data : String, pos : Int) -> (Int, Int) {
  let byte = data[pos].to_int() & 0xFF
  ((byte >> 1), pos + 1)
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
  let data = "\x03"
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