// Bio.SeqIO.TwoBitIO - Parser and writer for the UCSC 2bit binary genome format.
//
// The 2bit format compactly stores genome sequences using 2 bits per nucleotide:
//   A = 00, C = 01, T = 10, G = 11
//
// All multi-byte integers in the file are LITTLE-ENDIAN.
//
// File structure:
//   1. Header (16 bytes): magic, version, sequence count, reserved
//   2. Index: per-sequence (name length, name, offset)
//   3. Per-sequence data: dna size, packed dna, N blocks, mask blocks, reserved

// ============================================================================
// Constants
// ============================================================================

///|
/// Magic number for the 2bit format: 0x1A412743 (stored as Int).
let twobit_magic : Int = 0x1A412743

///|
/// Supported 2bit format version (0).
let twobit_version : Int = 0

// ============================================================================
// TwoBitNBlock - a contiguous region of N (unknown) bases
// ============================================================================

///|
/// A block of N (unknown) bases within a sequence, defined by its 0-based
/// start position and length.
pub struct TwoBitNBlock {
  start : Int
  size : Int
}

///|
/// Construct an N-block at `start` spanning `size` bases.
pub fn TwoBitNBlock::new(start : Int, size : Int) -> TwoBitNBlock {
  TwoBitNBlock::{ start, size }
}

///|
/// 0-based start position of the N-block.
pub fn TwoBitNBlock::start(self : TwoBitNBlock) -> Int {
  self.start
}

///|
/// Number of N bases in the block.
pub fn TwoBitNBlock::size(self : TwoBitNBlock) -> Int {
  self.size
}

///|
/// Exclusive end position (start + size) of the N-block.
pub fn TwoBitNBlock::end(self : TwoBitNBlock) -> Int {
  self.start + self.size
}

///|
/// Human-readable representation.
pub fn TwoBitNBlock::to_string(self : TwoBitNBlock) -> String {
  "NBlock(start=" + self.start.to_string() + ", size=" + self.size.to_string() + ")"
}

// ============================================================================
// TwoBitMaskBlock - a contiguous region of soft-masked (lowercase) bases
// ============================================================================

///|
/// A block of soft-masked (lowercase) bases within a sequence, defined by
/// its 0-based start position and length.
pub struct TwoBitMaskBlock {
  start : Int
  size : Int
}

///|
/// Construct a mask-block at `start` spanning `size` bases.
pub fn TwoBitMaskBlock::new(start : Int, size : Int) -> TwoBitMaskBlock {
  TwoBitMaskBlock::{ start, size }
}

///|
/// 0-based start position of the mask-block.
pub fn TwoBitMaskBlock::start(self : TwoBitMaskBlock) -> Int {
  self.start
}

///|
/// Number of masked bases in the block.
pub fn TwoBitMaskBlock::size(self : TwoBitMaskBlock) -> Int {
  self.size
}

///|
/// Exclusive end position (start + size) of the mask-block.
pub fn TwoBitMaskBlock::end(self : TwoBitMaskBlock) -> Int {
  self.start + self.size
}

///|
/// Human-readable representation.
pub fn TwoBitMaskBlock::to_string(self : TwoBitMaskBlock) -> String {
  "MaskBlock(start=" + self.start.to_string() + ", size=" + self.size.to_string() + ")"
}

// ============================================================================
// TwoBitSequence - a single sequence in a 2bit file
// ============================================================================

///|
/// A single sequence stored in 2bit format, containing the packed DNA bytes,
/// N-blocks (unknown regions), and mask-blocks (lowercase regions).
pub struct TwoBitSequence {
  name : String
  dna_size : Int
  mut packed_dna : Array[Int]
  mut n_blocks : Array[TwoBitNBlock]
  mut mask_blocks : Array[TwoBitMaskBlock]
}

///|
/// Construct a TwoBitSequence with the given name and DNA size. The packed
/// DNA, N-blocks, and mask-blocks start empty and can be populated via
/// `add_n_block`, `add_mask_block`, or by packing a sequence string.
pub fn TwoBitSequence::new(name : String, dna_size : Int) -> TwoBitSequence {
  TwoBitSequence::{
    name,
    dna_size,
    packed_dna: Array::new(),
    n_blocks: Array::new(),
    mask_blocks: Array::new(),
  }
}

///|
/// Sequence name.
pub fn TwoBitSequence::name(self : TwoBitSequence) -> String {
  self.name
}

///|
/// Number of bases in the sequence.
pub fn TwoBitSequence::dna_size(self : TwoBitSequence) -> Int {
  self.dna_size
}

///|
/// Packed DNA bytes (each byte holds 4 bases, 2 bits each).
pub fn TwoBitSequence::packed_dna(self : TwoBitSequence) -> Array[Int] {
  self.packed_dna.copy()
}

///|
/// N-blocks (regions of unknown bases).
pub fn TwoBitSequence::n_blocks(self : TwoBitSequence) -> Array[TwoBitNBlock] {
  self.n_blocks.copy()
}

///|
/// Mask-blocks (regions of soft-masked / lowercase bases).
pub fn TwoBitSequence::mask_blocks(self : TwoBitSequence) -> Array[TwoBitMaskBlock] {
  self.mask_blocks.copy()
}

///|
/// Add an N-block covering `size` bases starting at `start`.
pub fn TwoBitSequence::add_n_block(
  self : TwoBitSequence,
  start : Int,
  size : Int,
) -> Unit {
  self.n_blocks.push(TwoBitNBlock::new(start, size))
}

///|
/// Add a mask-block covering `size` bases starting at `start`.
pub fn TwoBitSequence::add_mask_block(
  self : TwoBitSequence,
  start : Int,
  size : Int,
) -> Unit {
  self.mask_blocks.push(TwoBitMaskBlock::new(start, size))
}

///|
/// Number of N-blocks.
pub fn TwoBitSequence::n_block_count(self : TwoBitSequence) -> Int {
  self.n_blocks.length()
}

///|
/// Number of mask-blocks.
pub fn TwoBitSequence::mask_block_count(self : TwoBitSequence) -> Int {
  self.mask_blocks.length()
}

// ============================================================================
// TwoBitFile - a complete 2bit file
// ============================================================================

///|
/// A complete 2bit file containing one or more sequences.
pub struct TwoBitFile {
  sequences : Array[TwoBitSequence]
  version : Int
}

///|
/// Construct an empty TwoBitFile (version 0).
pub fn TwoBitFile::new() -> TwoBitFile {
  TwoBitFile::{ sequences: Array::new(), version: twobit_version }
}

///|
/// Number of sequences in the file.
pub fn TwoBitFile::n_sequences(self : TwoBitFile) -> Int {
  self.sequences.length()
}

///|
/// All sequences in the file.
pub fn TwoBitFile::sequences(self : TwoBitFile) -> Array[TwoBitSequence] {
  self.sequences.copy()
}

///|
/// Format version (always 0).
pub fn TwoBitFile::version(self : TwoBitFile) -> Int {
  self.version
}

///|
/// Append a sequence to the file.
pub fn TwoBitFile::add_sequence(self : TwoBitFile, seq : TwoBitSequence) -> Unit {
  self.sequences.push(seq)
}

///|
/// Look up a sequence by name, returning `Some(seq)` if found.
pub fn TwoBitFile::get_sequence(
  self : TwoBitFile,
  name : String,
) -> TwoBitSequence? {
  let mut found : TwoBitSequence? = None
  for i in 0..<self.sequences.length() {
    if self.sequences[i].name == name {
      found = Some(self.sequences[i])
      break
    }
  }
  found
}

///|
/// Names of all sequences in the file.
pub fn TwoBitFile::sequence_names(self : TwoBitFile) -> Array[String] {
  let names : Array[String] = Array::new()
  for i in 0..<self.sequences.length() {
    names.push(self.sequences[i].name)
  }
  names
}

// ============================================================================
// Core encoding / decoding
// ============================================================================

///|
/// Convert a nucleotide character to its 2-bit code:
/// A -> 0, C -> 1, T -> 2, G -> 3. Returns -1 for N or any other character.
pub fn twobit_base_to_code(base : Char) -> Int {
  let c = base.to_int()
  if c == 'A'.to_int() {
    0
  } else if c == 'C'.to_int() {
    1
  } else if c == 'T'.to_int() {
    2
  } else if c == 'G'.to_int() {
    3
  } else {
    -1
  }
}

///|
/// Convert a 2-bit code back to a nucleotide character:
/// 0 -> A, 1 -> C, 2 -> T, 3 -> G.
pub fn twobit_code_to_base(code : Int) -> Char {
  if code == 0 {
    'A'
  } else if code == 1 {
    'C'
  } else if code == 2 {
    'T'
  } else {
    'G'
  }
}

///|
/// Pack a DNA string into 2-bit bytes, detecting N regions and lowercase
/// (soft-masked) regions.
///
/// Returns a tuple of `(packed_bytes, n_blocks, mask_blocks)`:
///   - Lowercase letters are upper-cased for packing but recorded as mask-blocks.
///   - N / n characters are recorded as N-blocks and packed as A (code 0).
pub fn twobit_pack_sequence(
  sequence : String,
) -> (Array[Int], Array[TwoBitNBlock], Array[TwoBitMaskBlock]) {
  let packed : Array[Int] = Array::new()
  let n_blocks : Array[TwoBitNBlock] = Array::new()
  let mask_blocks : Array[TwoBitMaskBlock] = Array::new()
  let len = sequence.length()
  let n_bytes = (len + 3) / 4
  // Initialize all packed bytes to 0.
  for _ in 0..<n_bytes {
    packed.push(0)
  }
  // ASCII code constants.
  let code_A = 'A'.to_int()
  let code_C = 'C'.to_int()
  let code_T = 'T'.to_int()
  let code_G = 'G'.to_int()
  let code_N = 'N'.to_int()
  let code_n = 'n'.to_int()
  let code_a = 'a'.to_int()
  let code_z = 'z'.to_int()
  let mut n_start = 0
  let mut n_size = 0
  let mut mask_start = 0
  let mut mask_size = 0
  for i in 0..<len {
    let c = sequence[i].to_int()
    let is_n = c == code_N || c == code_n
    let is_lower = c >= code_a && c <= code_z
    // Determine the 2-bit code for packing.
    let code : Int = if is_n {
      0 // N packed as A
    } else {
      let upper_c = if is_lower { c - 32 } else { c }
      if upper_c == code_A {
        0
      } else if upper_c == code_C {
        1
      } else if upper_c == code_T {
        2
      } else if upper_c == code_G {
        3
      } else {
        0 // unknown base, pack as A
      }
    }
    // Pack code into the correct byte and bit position.
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    packed[byte_idx] = packed[byte_idx] | (code << shift)
    // Update N-block tracking.
    if is_n {
      // Close any open mask-block.
      if mask_size > 0 {
        mask_blocks.push(TwoBitMaskBlock::new(mask_start, mask_size))
        mask_size = 0
      }
      if n_size == 0 {
        n_start = i
      }
      n_size = n_size + 1
    } else {
      // Close any open N-block.
      if n_size > 0 {
        n_blocks.push(TwoBitNBlock::new(n_start, n_size))
        n_size = 0
      }
      if is_lower {
        if mask_size == 0 {
          mask_start = i
        }
        mask_size = mask_size + 1
      } else {
        // Close any open mask-block.
        if mask_size > 0 {
          mask_blocks.push(TwoBitMaskBlock::new(mask_start, mask_size))
          mask_size = 0
        }
      }
    }
  }
  // Close any remaining open blocks.
  if n_size > 0 {
    n_blocks.push(TwoBitNBlock::new(n_start, n_size))
  }
  if mask_size > 0 {
    mask_blocks.push(TwoBitMaskBlock::new(mask_start, mask_size))
  }
  (packed, n_blocks, mask_blocks)
}

///|
/// Unpack 2-bit bytes back into a DNA string.
///
/// Applies N-blocks (replacing regions with 'N') and mask-blocks (converting
/// regions to lowercase) after decoding the packed bases.
pub fn twobit_unpack_sequence(
  packed : Array[Int],
  dna_size : Int,
  n_blocks : Array[TwoBitNBlock],
  mask_blocks : Array[TwoBitMaskBlock],
) -> String {
  // Decode each base from the packed bytes into an array of ASCII codes.
  let codes : Array[Int] = Array::new()
  for i in 0..<dna_size {
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    let code = (packed[byte_idx] >> shift) & 3
    codes.push(twobit_code_to_base(code).to_int())
  }
  // Apply N-blocks: replace bases with 'N'.
  for j in 0..<n_blocks.length() {
    let block = n_blocks[j]
    for k in 0..<block.size {
      let idx = block.start + k
      if idx >= 0 && idx < dna_size {
        codes[idx] = 'N'.to_int()
      }
    }
  }
  // Apply mask-blocks: convert bases to lowercase (add 32 to ASCII code).
  for j in 0..<mask_blocks.length() {
    let block = mask_blocks[j]
    for k in 0..<block.size {
      let idx = block.start + k
      if idx >= 0 && idx < dna_size {
        codes[idx] = codes[idx] + 32
      }
    }
  }
  // Build the result string.
  let sb = StringBuilder::new()
  for i in 0..<codes.length() {
    sb.write_char(Int::unsafe_to_char(codes[i]))
  }
  sb.to_string()
}

// ============================================================================
// Little-endian helpers (private)
// ============================================================================

///|
/// Encode a 32-bit unsigned integer as 4 little-endian bytes.
fn twobit_encode_u32_le(v : Int) -> Array[Int] {
  [v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF, (v >> 24) & 0xFF]
}

///|
/// Read a 32-bit little-endian unsigned integer from a byte array at `pos`.
fn twobit_read_u32_le(bytes : Array[Int], pos : Int) -> Int {
  bytes[pos] | (bytes[pos + 1] << 8) | (bytes[pos + 2] << 16) | (bytes[pos + 3] << 24)
}

// ============================================================================
// Serialization (to / from byte arrays)
// ============================================================================

///|
/// Serialize a TwoBitFile to a byte array in 2bit binary format (little-endian).
pub fn twobit_to_bytes(tbf : TwoBitFile) -> Array[Int] {
  let bytes : Array[Int] = Array::new()
  let seqs = tbf.sequences
  let n_seq = seqs.length()
  // Compute the total index size to determine sequence offsets.
  let mut index_size = 0
  for i in 0..<n_seq {
    index_size = index_size + 1 + seqs[i].name.length() + 4
  }
  // Compute the absolute offset of each sequence's data block.
  let offsets : Array[Int] = Array::new()
  let mut current_offset = 16 + index_size // 16 = header size
  for i in 0..<n_seq {
    offsets.push(current_offset)
    let seq = seqs[i]
    let packed_size = (seq.dna_size + 3) / 4
    let seq_data_size = 4 + packed_size + 4 + 8 * seq.n_blocks.length() + 4 + 8 * seq.mask_blocks.length() + 4
    current_offset = current_offset + seq_data_size
  }
  // ---- Header (16 bytes) ----
  let magic_bytes = twobit_encode_u32_le(twobit_magic)
  for i in 0..<4 {
    bytes.push(magic_bytes[i])
  }
  let ver_bytes = twobit_encode_u32_le(tbf.version)
  for i in 0..<4 {
    bytes.push(ver_bytes[i])
  }
  let cnt_bytes = twobit_encode_u32_le(n_seq)
  for i in 0..<4 {
    bytes.push(cnt_bytes[i])
  }
  // Reserved
  for i in 0..<4 {
    bytes.push(0)
  }
  // ---- Index ----
  for i in 0..<n_seq {
    let seq = seqs[i]
    let name_len = seq.name.length()
    bytes.push(name_len)
    for j in 0..<name_len {
      bytes.push(seq.name[j].to_int())
    }
    let off_bytes = twobit_encode_u32_le(offsets[i])
    for j in 0..<4 {
      bytes.push(off_bytes[j])
    }
  }
  // ---- Per-sequence data ----
  for i in 0..<n_seq {
    let seq = seqs[i]
    // DNA size
    let ds_bytes = twobit_encode_u32_le(seq.dna_size)
    for j in 0..<4 {
      bytes.push(ds_bytes[j])
    }
    // Packed DNA
    for j in 0..<seq.packed_dna.length() {
      bytes.push(seq.packed_dna[j])
    }
    // N-block count
    let nc_bytes = twobit_encode_u32_le(seq.n_blocks.length())
    for j in 0..<4 {
      bytes.push(nc_bytes[j])
    }
    // N-blocks
    for j in 0..<seq.n_blocks.length() {
      let block = seq.n_blocks[j]
      let s_bytes = twobit_encode_u32_le(block.start)
      for k in 0..<4 {
        bytes.push(s_bytes[k])
      }
      let sz_bytes = twobit_encode_u32_le(block.size)
      for k in 0..<4 {
        bytes.push(sz_bytes[k])
      }
    }
    // Mask-block count
    let mc_bytes = twobit_encode_u32_le(seq.mask_blocks.length())
    for j in 0..<4 {
      bytes.push(mc_bytes[j])
    }
    // Mask-blocks
    for j in 0..<seq.mask_blocks.length() {
      let block = seq.mask_blocks[j]
      let s_bytes = twobit_encode_u32_le(block.start)
      for k in 0..<4 {
        bytes.push(s_bytes[k])
      }
      let sz_bytes = twobit_encode_u32_le(block.size)
      for k in 0..<4 {
        bytes.push(sz_bytes[k])
      }
    }
    // Reserved
    for j in 0..<4 {
      bytes.push(0)
    }
  }
  bytes
}

///|
/// Deserialize a TwoBitFile from a byte array in 2bit binary format.
///
/// Reads and validates the magic number, then parses the header, index, and
/// per-sequence data.
pub fn twobit_from_bytes(data : Array[Int]) -> TwoBitFile {
  // Read and validate magic number.
  let magic = twobit_read_u32_le(data, 0)
  if magic != twobit_magic {
    return TwoBitFile::new()
  }
  let version = twobit_read_u32_le(data, 4)
  let seq_count = twobit_read_u32_le(data, 8)
  // Read the index.
  let names : Array[String] = Array::new()
  let offsets : Array[Int] = Array::new()
  let mut pos = 16
  for i in 0..<seq_count {
    let name_len = data[pos]
    pos = pos + 1
    let name_sb = StringBuilder::new()
    for j in 0..<name_len {
      name_sb.write_char(Int::unsafe_to_char(data[pos + j]))
    }
    pos = pos + name_len
    names.push(name_sb.to_string())
    offsets.push(twobit_read_u32_le(data, pos))
    pos = pos + 4
  }
  // Read per-sequence data.
  let sequences : Array[TwoBitSequence] = Array::new()
  for i in 0..<seq_count {
    let p = offsets[i]
    let dna_size = twobit_read_u32_le(data, p)
    let packed_size = (dna_size + 3) / 4
    let packed : Array[Int] = Array::new()
    let mut q = p + 4
    for j in 0..<packed_size {
      packed.push(data[q])
      q = q + 1
    }
    let n_count = twobit_read_u32_le(data, q)
    q = q + 4
    let n_blocks : Array[TwoBitNBlock] = Array::new()
    for j in 0..<n_count {
      let start = twobit_read_u32_le(data, q)
      let size = twobit_read_u32_le(data, q + 4)
      q = q + 8
      n_blocks.push(TwoBitNBlock::new(start, size))
    }
    let mask_count = twobit_read_u32_le(data, q)
    q = q + 4
    let mask_blocks : Array[TwoBitMaskBlock] = Array::new()
    for j in 0..<mask_count {
      let start = twobit_read_u32_le(data, q)
      let size = twobit_read_u32_le(data, q + 4)
      q = q + 8
      mask_blocks.push(TwoBitMaskBlock::new(start, size))
    }
    // Build the sequence and populate its mutable fields.
    let seq = TwoBitSequence::new(names[i], dna_size)
    seq.packed_dna = packed
    seq.n_blocks = n_blocks
    seq.mask_blocks = mask_blocks
    sequences.push(seq)
  }
  TwoBitFile::{ sequences, version }
}

// ============================================================================
// High-level I/O (hex string encoding)
// ============================================================================

///|
/// Convert an integer 0-15 to a lowercase hex digit character.
fn twobit_hex_digit(n : Int) -> Char {
  if n < 10 {
    Int::unsafe_to_char('0'.to_int() + n)
  } else {
    Int::unsafe_to_char('a'.to_int() + (n - 10))
  }
}

///|
/// Convert a hex digit character's ASCII code to its integer value 0-15.
fn twobit_hex_value(c : Int) -> Int {
  let code_0 = '0'.to_int()
  let code_9 = '9'.to_int()
  let code_a = 'a'.to_int()
  let code_f = 'f'.to_int()
  let code_A = 'A'.to_int()
  let code_F = 'F'.to_int()
  if c >= code_0 && c <= code_9 {
    c - code_0
  } else if c >= code_a && c <= code_f {
    c - code_a + 10
  } else if c >= code_A && c <= code_F {
    c - code_A + 10
  } else {
    0
  }
}

///|
/// Serialize a TwoBitFile to a hex string representation for text storage.
/// Each byte is encoded as two lowercase hex characters.
pub fn twobit_write(tbf : TwoBitFile) -> String {
  let bytes = twobit_to_bytes(tbf)
  let sb = StringBuilder::new()
  for i in 0..<bytes.length() {
    let b = bytes[i]
    sb.write_char(twobit_hex_digit((b >> 4) & 0xF))
    sb.write_char(twobit_hex_digit(b & 0xF))
  }
  sb.to_string()
}

///|
/// Deserialize a TwoBitFile from a hex string produced by `twobit_write`.
pub fn twobit_read(data : String) -> TwoBitFile {
  let bytes : Array[Int] = Array::new()
  let len = data.length()
  let mut i = 0
  while i + 1 < len {
    let high = twobit_hex_value(data[i].to_int())
    let low = twobit_hex_value(data[i + 1].to_int())
    bytes.push(high * 16 + low)
    i = i + 2
  }
  twobit_from_bytes(bytes)
}

///|
/// Extract a single sequence as a DNA string by name. Returns `None` if the
/// sequence is not found.
pub fn twobit_get_sequence_string(tbf : TwoBitFile, name : String) -> String? {
  match tbf.get_sequence(name) {
    Some(seq) =>
      Some(
        twobit_unpack_sequence(
          seq.packed_dna,
          seq.dna_size,
          seq.n_blocks,
          seq.mask_blocks,
        ),
      )
    None => None
  }
}

// ============================================================================
// SeqRecord integration
// ============================================================================

///|
/// Convert all sequences in a TwoBitFile to an array of SeqRecord objects.
/// Each sequence's name becomes the SeqRecord id and name.
pub fn twobit_to_seq_records(tbf : TwoBitFile) -> Array[SeqRecord] {
  let records : Array[SeqRecord] = Array::new()
  for i in 0..<tbf.sequences.length() {
    let seq = tbf.sequences[i]
    let dna_str = twobit_unpack_sequence(
      seq.packed_dna,
      seq.dna_size,
      seq.n_blocks,
      seq.mask_blocks,
    )
    let rec = SeqRecord::new(
      Seq::new(dna_str),
      id=seq.name,
      name=seq.name,
      description="",
    )
    records.push(rec)
  }
  records
}

///|
/// Build a TwoBitFile from an array of SeqRecord objects. Each record's
/// sequence is packed, and N / mask regions are detected automatically.
pub fn twobit_from_seq_records(records : Array[SeqRecord]) -> TwoBitFile {
  let tbf = TwoBitFile::new()
  for rec in records {
    let seq_str = rec.seq.to_string()
    let (packed, n_blocks, mask_blocks) = twobit_pack_sequence(seq_str)
    let seq = TwoBitSequence::new(rec.id, seq_str.length())
    seq.packed_dna = packed
    seq.n_blocks = n_blocks
    seq.mask_blocks = mask_blocks
    tbf.add_sequence(seq)
  }
  tbf
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a sample TwoBitFile with two sequences:
///   - "seq1": DNA "ATGCATGCATGC" (12 bases, no N or mask blocks)
///   - "seq2": DNA "ATGCNNNNATGC" (12 bases, one N-block at position 4, size 4)
pub fn twobit_sample_file() -> TwoBitFile {
  let tbf = TwoBitFile::new()
  // seq1: ATGCATGCATGC, no N or mask
  let (packed1, n1, m1) = twobit_pack_sequence("ATGCATGCATGC")
  let seq1 = TwoBitSequence::new("seq1", 12)
  seq1.packed_dna = packed1
  seq1.n_blocks = n1
  seq1.mask_blocks = m1
  tbf.add_sequence(seq1)
  // seq2: ATGCNNNNATGC, N-block at position 4 size 4
  let (packed2, n2, m2) = twobit_pack_sequence("ATGCNNNNATGC")
  let seq2 = TwoBitSequence::new("seq2", 12)
  seq2.packed_dna = packed2
  seq2.n_blocks = n2
  seq2.mask_blocks = m2
  tbf.add_sequence(seq2)
  tbf
}
