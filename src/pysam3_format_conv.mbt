// Pysam third-round: Format conversion and streaming (items 32-34)
// BAM→CRAM conversion, compression control, streaming I/O

///|
pub struct Pysam3CramHeader {
  version : String // CRAM version (e.g. "3.0")
  ref_id : String // reference sequence ID
  ref_md5 : String // MD5 of reference sequence
  ref_length : Int32 // reference length
}

///|
pub fn pysam3_cram_header_default() -> Pysam3CramHeader {
  { version: "3.0", ref_id: "", ref_md5: "", ref_length: 0, }
}

///|
pub fn pysam3_cram_header_to_bytes(header : Pysam3CramHeader) -> Array[Int32] {
  // CRAM file header encoding (simplified)
  let bytes : Array[Int32] = []
  // Magic: CRAM
  bytes.push(67)
  bytes.push(82)
  bytes.push(65)
  bytes.push(77)
  // Version
  for i in 0..<header.version.length() {
    bytes.push(header.version[i].to_int())
  }
  bytes.push(0) // null terminator
  // Ref ID
  for i in 0..<header.ref_id.length() {
    bytes.push(header.ref_id[i].to_int())
  }
  bytes.push(0)
  // Ref length (4 bytes, big-endian)
  bytes.push(header.ref_length / 16777216 % 256)
  bytes.push(header.ref_length / 65536 % 256)
  bytes.push(header.ref_length / 256 % 256)
  bytes.push(header.ref_length % 256)
  bytes
}

///|
pub fn pysam3_bam_to_cram_record(
  _chrom : String,
  pos : Int32,
  cigar : Array[(Int32, Int32)],
  seq : String,
  _qual : String,
  ref_seq : String,
) -> Array[Int32] {
  // Simplified CRAM record encoding
  // CRAM uses reference-based encoding: only differences from reference are stored
  let bytes : Array[Int32] = []
  // Position (varint-like, simplified to 4 bytes)
  bytes.push(pos % 256)
  bytes.push(pos / 256 % 256)
  bytes.push(pos / 65536 % 256)
  bytes.push(pos / 16777216 % 256)
  // CIGAR
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    bytes.push(op)
    bytes.push(length % 256)
    bytes.push(length / 256 % 256)
  }
  bytes.push(255) // CIGAR terminator
  // Sequence differences (only mismatches stored in CRAM)
  let mut ref_idx = pos
  let mut seq_idx = 0
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 0 {
      // M: compare with reference
      for k in 0..<length {
        let r = if ref_idx + k < ref_seq.length() {
          ref_seq[ref_idx + k].to_int()
        } else {
          78
        }
        let s = if seq_idx + k < seq.length() {
          seq[seq_idx + k].to_int()
        } else {
          78
        }
        if r != s {
          bytes.push(seq_idx + k) // position in read
          bytes.push(s) // base
        }
      }
      ref_idx = ref_idx + length
      seq_idx = seq_idx + length
    } else if op == 1 {
      // I: insertion — store all bases
      for k in 0..<length {
        if seq_idx + k < seq.length() {
          bytes.push(seq[seq_idx + k].to_int())
        }
      }
      seq_idx = seq_idx + length
    } else if op == 2 || op == 3 {
      ref_idx = ref_idx + length
    } else if op == 4 {
      seq_idx = seq_idx + length
    }
  }
  bytes.push(254) // record terminator
  bytes
}

///|
pub fn pysam3_cram_record_to_text(bytes : Array[Int32]) -> String {
  // Decode CRAM record back to text representation (for debugging)
  let sb = StringBuilder()
  let _ = sb.write("CRAM_record(len=")
  let _ = sb.write(bytes.length().to_string())
  let _ = sb.write(")")
  sb.to_string()
}

// --- Compression level control ---

///|
pub enum Pysam3CompressionLevel {
  Pysam3NoCompression // level 0
  Pysam3FastCompression // level 1
  Pysam3DefaultCompression // level 6
  Pysam3BestCompression // level 9
}

///|
pub fn pysam3_compression_level_value(level : Pysam3CompressionLevel) -> Int32 {
  match level {
    Pysam3NoCompression => 0
    Pysam3FastCompression => 1
    Pysam3DefaultCompression => 6
    Pysam3BestCompression => 9
  }
}

///|
pub fn pysam3_compression_level_from_int(
  level : Int32,
) -> Pysam3CompressionLevel {
  if level <= 0 {
    Pysam3NoCompression
  } else if level <= 1 {
    Pysam3FastCompression
  } else if level <= 6 {
    Pysam3DefaultCompression
  } else {
    Pysam3BestCompression
  }
}

///|
pub fn pysam3_write_bam_header(
  header_text : String,
  compression : Pysam3CompressionLevel,
) -> Array[Int32] {
  // BAM header with compression marker
  let bytes : Array[Int32] = []
  // BAM magic: BAM\1
  bytes.push(66)
  bytes.push(65)
  bytes.push(77)
  bytes.push(1)
  // Compression level
  bytes.push(pysam3_compression_level_value(compression))
  // Header text length (4 bytes, little-endian)
  let hlen = header_text.length()
  bytes.push(hlen % 256)
  bytes.push(hlen / 256 % 256)
  bytes.push(hlen / 65536 % 256)
  bytes.push(hlen / 16777216 % 256)
  // Header text
  for i in 0..<header_text.length() {
    bytes.push(header_text[i].to_int())
  }
  bytes
}

// --- Streaming I/O ---

///|
pub struct Pysam3StreamReader {
  buffer : Array[String]
  pos : Int32
}

///|
pub fn pysam3_stream_reader_new() -> Pysam3StreamReader {
  { buffer: [], pos: 0, }
}

///|
pub fn pysam3_stream_reader_push(
  reader : Pysam3StreamReader,
  line : String,
) -> Pysam3StreamReader {
  let new_buffer = []
  for i in 0..<reader.buffer.length() {
    new_buffer.push(reader.buffer[i])
  }
  new_buffer.push(line)
  { ..reader, buffer: new_buffer, }
}

///|
pub fn pysam3_stream_reader_next(
  reader : Pysam3StreamReader,
) -> (Pysam3StreamReader, String?) {
  if reader.pos >= reader.buffer.length() {
    return (reader, None)
  }
  let line = reader.buffer[reader.pos]
  ({ ..reader, pos: reader.pos + 1, }, Some(line))
}

///|
pub fn pysam3_stream_reader_has_next(reader : Pysam3StreamReader) -> Bool {
  reader.pos < reader.buffer.length()
}

///|
pub fn pysam3_stream_reader_reset(
  reader : Pysam3StreamReader,
) -> Pysam3StreamReader {
  { ..reader, pos: 0, }
}

///|
pub struct Pysam3StreamWriter {
  lines : Array[String]
}

///|
pub fn pysam3_stream_writer_new() -> Pysam3StreamWriter {
  { lines: [], }
}

///|
pub fn pysam3_stream_writer_write(
  writer : Pysam3StreamWriter,
  line : String,
) -> Pysam3StreamWriter {
  let new_lines = []
  for i in 0..<writer.lines.length() {
    new_lines.push(writer.lines[i])
  }
  new_lines.push(line)
  { lines: new_lines, }
}

///|
pub fn pysam3_stream_writer_to_string(writer : Pysam3StreamWriter) -> String {
  let sb = StringBuilder()
  for i in 0..<writer.lines.length() {
    let _ = sb.write(writer.lines[i])
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
pub fn pysam3_stream_process_sam(
  reader : Pysam3StreamReader,
  filter_fn : (String) -> Bool,
) -> Pysam3StreamWriter {
  // Stream SAM records through a filter function
  let mut writer = pysam3_stream_writer_new()
  let mut r = reader
  while pysam3_stream_reader_has_next(r) {
    let (next_r, line_opt) = pysam3_stream_reader_next(r)
    r = next_r
    match line_opt {
      Some(line) =>
        if line.has_prefix("@") || filter_fn(line) {
          writer = pysam3_stream_writer_write(writer, line)
        }
      None => break
    }
  }
  writer
}

///|
pub fn pysam3_stream_pipe(
  lines : Array[String],
  stages : Array[(String) -> String],
) -> Array[String] {
  // Apply a pipeline of transformations to each line
  let result : Array[String] = []
  for i in 0..<lines.length() {
    let mut line = lines[i]
    for j in 0..<stages.length() {
      line = stages[j](line)
    }
    result.push(line)
  }
  result
}
