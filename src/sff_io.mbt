///|
/// Bio.SeqIO.SffIO module - Standard Flowgram Format (SFF) binary parser.
/// Inspired by Biopython's Bio.SeqIO.SffIO module.
/// Parses 454/Roche SFF binary files containing flowgram data, sequences,
/// and per-base quality scores.

// ===== SFF binary format constants =====

///|
/// SFF magic number: 0x2E736666 (".sff" in ASCII).
let sff_magic : Int = 0x2E736666

///|
/// SFF version 1.0 (major=0, minor=1 encoded as 0x0001).
let sff_version : Int = 0x0001

///|
/// Flowgram format code 1 = unsigned 16-bit integer per flow value.
let sff_flowgram_format_u16 : Int = 1

// ===== Data structures =====

///|
/// SFF common header (file-level metadata).
pub struct SffHeader {
  /// Magic number (should be 0x2E736666)
  magic : Int
  /// Version (should be 0x0001)
  version : Int
  /// Index offset in the file
  index_offset : Int
  /// Index length
  index_length : Int
  /// Number of reads in the file
  n_reads : Int
  /// Header length (bytes, must be multiple of 8)
  header_length : Int
  /// Key length
  key_length : Int
  /// Number of flows
  n_flows : Int
  /// Flowgram format code (1 = u16)
  flowgram_format : Int
  /// Flow chars (e.g., "TACG" repeated)
  flow_chars : String
  /// Key sequence (e.g., "TCAG")
  key_sequence : String
}

///|
/// A single SFF read record.
pub struct SffRead {
  /// Read header length (bytes, must be multiple of 8)
  read_header_length : Int
  /// Name length
  name_length : Int
  /// Number of bases
  n_bases : Int
  /// Clip quality left (1-based position)
  clip_qual_left : Int
  /// Clip quality right
  clip_qual_right : Int
  /// Clip adapter left
  clip_adapter_left : Int
  /// Clip adapter right
  clip_adapter_right : Int
  /// Read name
  name : String
  /// Flowgram values (one per flow, typically scaled u16 / 100.0)
  flowgram : Array[Double]
  /// Flow index per base (which flow produced each base)
  flow_index : Array[Int]
  /// Called bases (ACGT, possibly N)
  bases : String
  /// Per-base quality scores (Phred)
  qualities : Array[Int]
}

///|
/// Complete SFF file representation.
pub struct SffFile {
  header : SffHeader
  reads : Array[SffRead]
}

// ===== Header methods =====

///|
/// Create a new SFF header.
pub fn SffHeader::new(n_reads : Int, n_flows : Int, flow_chars : String, key_sequence : String) -> SffHeader {
  let key_length = key_sequence.length()
  // Fixed header is 32 bytes + n_flows + key_length, padded to multiple of 8
  let base_length = 32 + n_flows + key_length
  let header_length = sff_pad_to_8(base_length)
  SffHeader::{
    magic: sff_magic,
    version: sff_version,
    index_offset: 0,
    index_length: 0,
    n_reads,
    header_length,
    key_length,
    n_flows,
    flowgram_format: sff_flowgram_format_u16,
    flow_chars,
    key_sequence
  }
}

///|
/// Get number of reads declared in the header.
pub fn SffHeader::sff_n_reads(self : SffHeader) -> Int {
  self.n_reads
}

///|
/// Get number of flows.
pub fn SffHeader::sff_n_flows(self : SffHeader) -> Int {
  self.n_flows
}

///|
/// Get flow chars.
pub fn SffHeader::sff_flow_chars(self : SffHeader) -> String {
  self.flow_chars
}

///|
/// Get key sequence.
pub fn SffHeader::sff_key_sequence(self : SffHeader) -> String {
  self.key_sequence
}

///|
/// Validate the header fields.
pub fn SffHeader::sff_validate(self : SffHeader) -> Array[String] {
  let errors : Array[String] = Array::new()
  if self.magic != sff_magic {
    errors.push("Invalid magic number: expected " + sff_magic.to_string() + ", got " + self.magic.to_string())
  }
  if self.version != sff_version {
    errors.push("Unsupported version: " + self.version.to_string())
  }
  if self.header_length % 8 != 0 {
    errors.push("Header length must be a multiple of 8, got " + self.header_length.to_string())
  }
  if self.flowgram_format != sff_flowgram_format_u16 {
    errors.push("Unsupported flowgram format: " + self.flowgram_format.to_string() + " (only format 1/u16 is supported)")
  }
  if self.flow_chars.length() != self.n_flows {
    errors.push("Flow chars length (" + self.flow_chars.length().to_string() + ") does not match n_flows (" + self.n_flows.to_string() + ")")
  }
  if self.key_sequence.length() != self.key_length {
    errors.push("Key sequence length does not match key_length")
  }
  errors
}

// ===== Read methods =====

///|
/// Create a new SFF read record.
pub fn SffRead::new(name : String, bases : String, qualities : Array[Int], flowgram : Array[Double], flow_index : Array[Int], clip_qual_left : Int, clip_qual_right : Int, clip_adapter_left : Int, clip_adapter_right : Int) -> SffRead {
  let name_length = name.length()
  let n_bases = bases.length()
  // read_header_length = 16 (fixed) + name_length, padded to multiple of 8
  let base_rh_length = 16 + name_length
  let read_header_length = sff_pad_to_8(base_rh_length)
  SffRead::{
    read_header_length,
    name_length,
    n_bases,
    clip_qual_left,
    clip_qual_right,
    clip_adapter_left,
    clip_adapter_right,
    name,
    flowgram,
    flow_index,
    bases,
    qualities
  }
}

///|
/// Get read name.
pub fn SffRead::sff_name(self : SffRead) -> String {
  self.name
}

///|
/// Get called bases.
pub fn SffRead::sff_bases(self : SffRead) -> String {
  self.bases
}

///|
/// Get number of bases.
pub fn SffRead::sff_n_bases(self : SffRead) -> Int {
  self.n_bases
}

///|
/// Get quality scores.
pub fn SffRead::sff_qualities(self : SffRead) -> Array[Int] {
  self.qualities.copy()
}

///|
/// Get flowgram values.
pub fn SffRead::sff_flowgram(self : SffRead) -> Array[Double] {
  self.flowgram.copy()
}

///|
/// Get flow index per base.
pub fn SffRead::sff_flow_index(self : SffRead) -> Array[Int] {
  self.flow_index.copy()
}

///|
/// Get the clipped sequence (after applying quality clipping).
pub fn SffRead::sff_clipped_bases(self : SffRead) -> String {
  let left = if self.clip_qual_left > 0 { self.clip_qual_left - 1 } else { 0 }
  let right = if self.clip_qual_right > 0 { self.clip_qual_right } else { self.n_bases }
  if left >= right || left >= self.n_bases {
    return ""
  }
  let end = if right > self.n_bases { self.n_bases } else { right }
  self.bases[left:end].to_owned()
}

///|
/// Get the clipped quality scores.
pub fn SffRead::sff_clipped_qualities(self : SffRead) -> Array[Int] {
  let left = if self.clip_qual_left > 0 { self.clip_qual_left - 1 } else { 0 }
  let right = if self.clip_qual_right > 0 { self.clip_qual_right } else { self.n_bases }
  let result : Array[Int] = Array::new()
  if left >= right || left >= self.n_bases {
    return result
  }
  let end = if right > self.n_bases { self.n_bases } else { right }
  let mut i = left
  while i < end {
    result.push(self.qualities[i])
    i = i + 1
  }
  result
}

///|
/// Mean quality of the clipped region.
pub fn SffRead::sff_mean_quality(self : SffRead) -> Double {
  let q = self.sff_clipped_qualities()
  if q.length() == 0 {
    return 0.0
  }
  let mut sum = 0
  let mut i = 0
  while i < q.length() {
    sum = sum + q[i]
    i = i + 1
  }
  sum.to_double() / q.length().to_double()
}

// ===== SffFile methods =====

///|
/// Create a new SFF file from header and reads.
pub fn SffFile::new(header : SffHeader, reads : Array[SffRead]) -> SffFile {
  SffFile::{ header, reads }
}

///|
/// Get the header.
pub fn SffFile::sff_header(self : SffFile) -> SffHeader {
  self.header
}

///|
/// Get all reads.
pub fn SffFile::sff_reads(self : SffFile) -> Array[SffRead] {
  self.reads.copy()
}

///|
/// Get number of reads.
pub fn SffFile::sff_n_reads(self : SffFile) -> Int {
  self.reads.length()
}

///|
/// Get a read by index.
pub fn SffFile::sff_get_read(self : SffFile, idx : Int) -> SffRead {
  self.reads[idx]
}

///|
/// Get a read by name (returns Option).
pub fn SffFile::sff_find_read(self : SffFile, name : String) -> SffRead? {
  let mut found : SffRead? = None
  let mut i = 0
  while i < self.reads.length() {
    if self.reads[i].name == name {
      found = Some(self.reads[i])
      i = self.reads.length()
    }
    i = i + 1
  }
  found
}

///|
/// Total number of bases across all reads.
pub fn SffFile::sff_total_bases(self : SffFile) -> Int {
  let mut total = 0
  let mut i = 0
  while i < self.reads.length() {
    total = total + self.reads[i].n_bases
    i = i + 1
  }
  total
}

// ===== Binary encoding helpers =====

///|
/// Pad a length up to the next multiple of 8.
fn sff_pad_to_8(n : Int) -> Int {
  let rem = n % 8
  if rem == 0 { n } else { n + (8 - rem) }
}

///|
/// Encode a 32-bit unsigned integer as 4 big-endian bytes.
fn sff_encode_u32(v : Int) -> Array[Int] {
  [
    (v / 16777216) % 256,
    (v / 65536) % 256,
    (v / 256) % 256,
    v % 256
  ]
}

///|
/// Encode a 16-bit unsigned integer as 2 big-endian bytes.
fn sff_encode_u16(v : Int) -> Array[Int] {
  [(v / 256) % 256, v % 256]
}

///|
/// Encode a 64-bit unsigned integer as 8 big-endian bytes.
/// Note: MoonBit Int is 32-bit, so we only encode the low 32 bits with 4 leading zero bytes.
fn sff_encode_u64(v : Int) -> Array[Int] {
  // High 32 bits = 0 (our values fit in 32 bits), low 32 bits = v
  let low_bytes = sff_encode_u32(v)
  [0, 0, 0, 0, low_bytes[0], low_bytes[1], low_bytes[2], low_bytes[3]]
}

///|
/// Read a 32-bit big-endian unsigned integer from a byte array at a position.
fn sff_read_u32(bytes : Array[Int], pos : Int) -> Int {
  bytes[pos] * 16777216 + bytes[pos + 1] * 65536 + bytes[pos + 2] * 256 + bytes[pos + 3]
}

///|
/// Read a 16-bit big-endian unsigned integer from a byte array at a position.
fn sff_read_u16(bytes : Array[Int], pos : Int) -> Int {
  bytes[pos] * 256 + bytes[pos + 1]
}

///|
/// Read a 64-bit big-endian unsigned integer from a byte array at a position.
/// Note: Only returns the low 32 bits (MoonBit Int is 32-bit).
fn sff_read_u64(bytes : Array[Int], pos : Int) -> Int {
  // High 32 bits at pos..pos+3 (ignored, assumed 0), low 32 bits at pos+4..pos+7
  sff_read_u32(bytes, pos + 4)
}

///|
/// Encode a string as ASCII bytes.
fn sff_encode_string(s : String) -> Array[Int] {
  let bytes : Array[Int] = Array::new()
  let mut i = 0
  while i < s.length() {
    bytes.push(s[i].to_int())
    i = i + 1
  }
  bytes
}

// ===== SFF binary encoder =====

///|
/// Encode an SFF file as a byte array (big-endian binary format).
pub fn sff_encode(file : SffFile) -> Array[Int] {
  let bytes : Array[Int] = Array::new()
  // ===== Common Header (31 bytes fixed + flow_chars + key_sequence, padded to 8) =====
  // magic (4 bytes)
  let magic_bytes = sff_encode_u32(file.header.magic)
  let mut mi = 0
  while mi < magic_bytes.length() {
    bytes.push(magic_bytes[mi])
    mi = mi + 1
  }
  // version (4 bytes; sff version is stored as major=0x00, minor=0x01 followed by 2 padding bytes)
  bytes.push(0)
  bytes.push(1)
  bytes.push(0)
  bytes.push(0)
  // index offset (8 bytes)
  let idx_off = sff_encode_u64(file.header.index_offset)
  let mut ii = 0
  while ii < idx_off.length() {
    bytes.push(idx_off[ii])
    ii = ii + 1
  }
  // index length (4 bytes)
  let idx_len = sff_encode_u32(file.header.index_length)
  let mut il = 0
  while il < idx_len.length() {
    bytes.push(idx_len[il])
    il = il + 1
  }
  // number of reads (4 bytes)
  let n_reads_b = sff_encode_u32(file.header.n_reads)
  let mut ni = 0
  while ni < n_reads_b.length() {
    bytes.push(n_reads_b[ni])
    ni = ni + 1
  }
  // header length (2 bytes)
  let hl = sff_encode_u16(file.header.header_length)
  bytes.push(hl[0])
  bytes.push(hl[1])
  // key length (2 bytes)
  let kl = sff_encode_u16(file.header.key_length)
  bytes.push(kl[0])
  bytes.push(kl[1])
  // number of flows (2 bytes)
  let nf = sff_encode_u16(file.header.n_flows)
  bytes.push(nf[0])
  bytes.push(nf[1])
  // flowgram format code (2 bytes)
  let ff = sff_encode_u16(file.header.flowgram_format)
  bytes.push(ff[0])
  bytes.push(ff[1])
  // flow chars (n_flows bytes)
  let flow_bytes = sff_encode_string(file.header.flow_chars)
  let mut fi = 0
  while fi < flow_bytes.length() {
    bytes.push(flow_bytes[fi])
    fi = fi + 1
  }
  // key sequence (key_length bytes)
  let key_bytes = sff_encode_string(file.header.key_sequence)
  let mut ki = 0
  while ki < key_bytes.length() {
    bytes.push(key_bytes[ki])
    ki = ki + 1
  }
  // Padding to multiple of 8 (header_length is already padded)
  let header_so_far = 32 + file.header.n_flows + file.header.key_length
  let pad_needed = file.header.header_length - header_so_far
  let mut p = 0
  while p < pad_needed {
    bytes.push(0)
    p = p + 1
  }
  // ===== Per-read data =====
  let mut r = 0
  while r < file.reads.length() {
    let read = file.reads[r]
    // read header length (2 bytes)
    let rhl = sff_encode_u16(read.read_header_length)
    bytes.push(rhl[0])
    bytes.push(rhl[1])
    // name length (2 bytes)
    let nl = sff_encode_u16(read.name_length)
    bytes.push(nl[0])
    bytes.push(nl[1])
    // number of bases (4 bytes)
    let nb = sff_encode_u32(read.n_bases)
    let mut nbi = 0
    while nbi < nb.length() {
      bytes.push(nb[nbi])
      nbi = nbi + 1
    }
    // clip_qual_left (2 bytes)
    let cql = sff_encode_u16(read.clip_qual_left)
    bytes.push(cql[0])
    bytes.push(cql[1])
    // clip_qual_right (2 bytes)
    let cqr = sff_encode_u16(read.clip_qual_right)
    bytes.push(cqr[0])
    bytes.push(cqr[1])
    // clip_adapter_left (2 bytes)
    let cal = sff_encode_u16(read.clip_adapter_left)
    bytes.push(cal[0])
    bytes.push(cal[1])
    // clip_adapter_right (2 bytes)
    let car = sff_encode_u16(read.clip_adapter_right)
    bytes.push(car[0])
    bytes.push(car[1])
    // name (name_length bytes)
    let name_b = sff_encode_string(read.name)
    let mut nmi = 0
    while nmi < name_b.length() {
      bytes.push(name_b[nmi])
      nmi = nmi + 1
    }
    // Padding to multiple of 8 for read header
    let rh_so_far = 16 + read.name_length
    let rh_pad = read.read_header_length - rh_so_far
    let mut rp = 0
    while rp < rh_pad {
      bytes.push(0)
      rp = rp + 1
    }
    // flowgram values (n_flows * 2 bytes, u16, scaled by 100)
    let mut fi2 = 0
    while fi2 < read.flowgram.length() {
      let fv = (read.flowgram[fi2] * 100.0 + 0.5).to_int()
      let fvb = sff_encode_u16(fv)
      bytes.push(fvb[0])
      bytes.push(fvb[1])
      fi2 = fi2 + 1
    }
    // flow index per base (n_bases bytes)
    let mut fii = 0
    while fii < read.flow_index.length() {
      bytes.push(read.flow_index[fii])
      fii = fii + 1
    }
    // bases (n_bases bytes)
    let bases_b = sff_encode_string(read.bases)
    let mut bi = 0
    while bi < bases_b.length() {
      bytes.push(bases_b[bi])
      bi = bi + 1
    }
    // quality scores (n_bases bytes)
    let mut qi = 0
    while qi < read.qualities.length() {
      bytes.push(read.qualities[qi])
      qi = qi + 1
    }
    // Padding to multiple of 8 for read data section
    let data_len = read.flowgram.length() * 2 + read.flow_index.length() + read.bases.length() + read.qualities.length()
    let data_pad = sff_pad_to_8(data_len) - data_len
    let mut dp = 0
    while dp < data_pad {
      bytes.push(0)
      dp = dp + 1
    }
    r = r + 1
  }
  bytes
}

// ===== SFF binary parser =====

///|
/// Parse an SFF byte array into an SffFile.
/// Returns Option to handle malformed input gracefully.
pub fn sff_parse(bytes : Array[Int]) -> SffFile? {
  if bytes.length() < 31 {
    return None
  }
  // ===== Common Header =====
  let magic = sff_read_u32(bytes, 0)
  if magic != sff_magic {
    return None
  }
  // version at offset 4: bytes[4]=0 (major), bytes[5]=1 (minor), bytes[6..7]=padding
  let version = bytes[4] * 256 + bytes[5]
  let index_offset = sff_read_u64(bytes, 8)
  let index_length = sff_read_u32(bytes, 16)
  let n_reads = sff_read_u32(bytes, 20)
  let header_length = sff_read_u16(bytes, 24)
  let key_length = sff_read_u16(bytes, 26)
  let n_flows = sff_read_u16(bytes, 28)
  let flowgram_format = sff_read_u16(bytes, 30)
  // flow chars at offset 32, length n_flows
  let flow_chars = sff_read_string(bytes, 32, n_flows)
  // key sequence at offset 32 + n_flows, length key_length
  let key_sequence = sff_read_string(bytes, 32 + n_flows, key_length)
  let header = SffHeader::{
    magic,
    version,
    index_offset,
    index_length,
    n_reads,
    header_length,
    key_length,
    n_flows,
    flowgram_format,
    flow_chars,
    key_sequence
  }
  // ===== Per-read data =====
  let reads : Array[SffRead] = Array::new()
  let mut pos = header_length
  let mut r = 0
  while r < n_reads {
    if pos + 16 > bytes.length() {
      break
    }
    let read_header_length = sff_read_u16(bytes, pos)
    let name_length = sff_read_u16(bytes, pos + 2)
    let n_bases = sff_read_u32(bytes, pos + 4)
    let clip_qual_left = sff_read_u16(bytes, pos + 8)
    let clip_qual_right = sff_read_u16(bytes, pos + 10)
    let clip_adapter_left = sff_read_u16(bytes, pos + 12)
    let clip_adapter_right = sff_read_u16(bytes, pos + 14)
    let name = sff_read_string(bytes, pos + 16, name_length)
    // Data section starts after read header
    let data_start = pos + read_header_length
    // flowgram values: n_flows * 2 bytes
    let flowgram : Array[Double] = Array::new()
    let mut fi = 0
    while fi < n_flows {
      let fv = sff_read_u16(bytes, data_start + fi * 2)
      flowgram.push(fv.to_double() / 100.0)
      fi = fi + 1
    }
    let flow_index_start = data_start + n_flows * 2
    // flow index per base: n_bases bytes
    let flow_index : Array[Int] = Array::new()
    let mut fii = 0
    while fii < n_bases {
      flow_index.push(bytes[flow_index_start + fii])
      fii = fii + 1
    }
    let bases_start = flow_index_start + n_bases
    // bases: n_bases bytes
    let bases = sff_read_string(bytes, bases_start, n_bases)
    let qual_start = bases_start + n_bases
    // quality scores: n_bases bytes
    let qualities : Array[Int] = Array::new()
    let mut qi = 0
    while qi < n_bases {
      qualities.push(bytes[qual_start + qi])
      qi = qi + 1
    }
    reads.push(SffRead::{
      read_header_length,
      name_length,
      n_bases,
      clip_qual_left,
      clip_qual_right,
      clip_adapter_left,
      clip_adapter_right,
      name,
      flowgram,
      flow_index,
      bases,
      qualities
    })
    // Advance to next read (with padding)
    let data_len = n_flows * 2 + n_bases + n_bases + n_bases
    let padded_data = sff_pad_to_8(data_len)
    pos = data_start + padded_data
    r = r + 1
  }
  Some(SffFile::{ header, reads })
}

///|
/// Read a fixed-length ASCII string from a byte array.
fn sff_read_string(bytes : Array[Int], pos : Int, length : Int) -> String {
  let mut s = ""
  let mut i = 0
  while i < length && pos + i < bytes.length() {
    let byte_val = bytes[pos + i]
    let c = byte_val.to_char().unwrap_or('?')
    s = s + c.to_string()
    i = i + 1
  }
  s
}

// ===== Example data generators =====

///|
/// Create a small example SFF file with 2 reads for testing.
pub fn sff_create_example_file() -> SffFile {
  let header = SffHeader::new(2, 4, "TACG", "TCAG")
  let read1 = SffRead::new(
    "READ0001",
    "ACGTACGTAC",
    [10, 15, 20, 25, 30, 28, 26, 24, 22, 20],
    [1.05, 0.02, 0.98, 1.01],
    [1, 2, 3, 4, 1, 2, 3, 4, 1, 2],
    1,
    10,
    0,
    0
  )
  let read2 = SffRead::new(
    "READ0002",
    "TGCATGCATG",
    [12, 18, 22, 28, 32, 30, 28, 26, 24, 22],
    [0.99, 1.01, 1.02, 0.98],
    [1, 2, 3, 4, 1, 2, 3, 4, 1, 2],
    1,
    10,
    0,
    0
  )
  SffFile::new(header, [read1, read2])
}

///|
/// Create a single example read for testing.
pub fn sff_create_example_read() -> SffRead {
  SffRead::new(
    "TEST001",
    "ACGTACGT",
    [20, 25, 30, 28, 26, 24, 22, 20],
    [1.0, 0.0, 1.0, 0.0],
    [1, 3, 1, 3, 1, 3, 1, 3],
    1,
    8,
    0,
    0
  )
}

///|
/// Create an example SFF header for testing.
pub fn sff_create_example_header() -> SffHeader {
  SffHeader::new(1, 4, "TACG", "TCAG")
}
