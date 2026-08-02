// Bio.SeqIO.NibIO - Reader and writer for the UCSC nib 2-bit binary sequence
// format, mirroring BioPython's Bio.SeqIO.NibIO.
//
// The nib format stores DNA sequences compactly using 2 bits per nucleotide:
//   T = 00 (0), C = 01 (1), A = 10 (2), G = 11 (3)
//
// Each byte packs 4 bases, MSB first: the first base occupies bits 7-6, the
// second bits 5-4, the third bits 3-2, and the fourth bits 1-0.
//
// Note: this implementation uses `Array[Int]` (each element masked to 0-255)
// for the packed data rather than `Array[UInt]`, because the rest of the
// codebase (xdna_io, twobit_io) consistently uses Int with masking for packed
// byte arrays, and UInt bitwise ergonomics are not exercised anywhere in the
// project. Each element still represents a single packed byte.

// ============================================================================
// Constants
// ============================================================================

///|
/// nib format version this implementation reads and writes.
let nib_version : Int = 1

// ============================================================================
// Base encoding / decoding
// ============================================================================

///|
/// Convert a nucleotide character to its nib 2-bit code:
/// T -> 0, C -> 1, A -> 2, G -> 3. Returns -1 for any other character
/// (including N and lowercase letters; use `nib_base_to_code_ci` for
/// case-insensitive lookup).
pub fn nib_base_to_code(base : Char) -> Int {
  let c = base.to_int()
  if c == 'T'.to_int() {
    0
  } else if c == 'C'.to_int() {
    1
  } else if c == 'A'.to_int() {
    2
  } else if c == 'G'.to_int() {
    3
  } else {
    -1
  }
}

///|
/// Case-insensitive variant of `nib_base_to_code`: lowercase a/c/t/g are
/// accepted and mapped to the same codes as their uppercase counterparts.
/// Returns -1 for non-ACGT characters (including N).
pub fn nib_base_to_code_ci(base : Char) -> Int {
  let c = base.to_int()
  let code_a = 'A'.to_int()
  let code_z = 'Z'.to_int()
  let code_a_lower = 'a'.to_int()
  let code_z_lower = 'z'.to_int()
  let upper_c = if c >= code_a_lower && c <= code_z_lower {
    c - 32
  } else if c >= code_a && c <= code_z {
    c
  } else {
    c
  }
  if upper_c == 'T'.to_int() {
    0
  } else if upper_c == 'C'.to_int() {
    1
  } else if upper_c == 'A'.to_int() {
    2
  } else if upper_c == 'G'.to_int() {
    3
  } else {
    -1
  }
}

///|
/// Convert a nib 2-bit code back to a nucleotide character:
/// 0 -> T, 1 -> C, 2 -> A, 3 -> G. Any other code maps to 'N'.
pub fn nib_code_to_base(code : Int) -> Char {
  if code == 0 {
    'T'
  } else if code == 1 {
    'C'
  } else if code == 2 {
    'A'
  } else if code == 3 {
    'G'
  } else {
    'N'
  }
}

// ============================================================================
// Core encoding / decoding
// ============================================================================

///|
/// Pack a DNA string into 2-bit bytes (4 bases per byte, MSB first). The input
/// is treated case-insensitively: lowercase a/c/t/g are upper-cased before
/// packing. Unknown bases (N, X, etc.) are packed as T (code 0); this is a
/// lossy transformation for non-ACGT characters. Returns ceil(length/4) bytes.
pub fn nib_encode(sequence : String) -> Array[Int] {
  let n = sequence.length()
  let nbytes = (n + 3) / 4
  let packed = Array::make(nbytes, 0)
  for i in 0..<n {
    let code_ci = nib_base_to_code_ci(sequence[i].unsafe_to_char())
    let code = if code_ci < 0 { 0 } else { code_ci }
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    packed[byte_idx] = packed[byte_idx] | (code << shift)
  }
  packed
}

///|
/// Unpack 2-bit bytes back into a DNA string of length `length`. Always
/// returns uppercase characters. If `packed` is shorter than
/// ceil(length/4), missing bytes are treated as 0 (T).
pub fn nib_decode(packed : Array[Int], length : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<length {
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    let code = if byte_idx < packed.length() {
      (packed[byte_idx] >> shift) & 3
    } else {
      0
    }
    sb.write_char(nib_code_to_base(code))
  }
  sb.to_string()
}

// ============================================================================
// NibSequence - a single sequence in nib format
// ============================================================================

///|
/// A single sequence stored in nib format. The packed data uses 2 bits per
/// base (4 bases per byte, MSB first). The `is_uppercase` flag records whether
/// the original sequence was entirely uppercase; if false, `nib_to_string`
/// lowercases the decoded bases.
pub struct NibSequence {
  seq_id : String
  length : Int
  packed_data : Array[Int]
  is_uppercase : Bool
}

///|
/// Construct a NibSequence from a DNA string. The sequence is packed using
/// 2-bit encoding; lowercase letters are accepted and packed as their
/// uppercase equivalents. The `is_uppercase` flag is set to false if any
/// character in `sequence` is a lowercase letter, true otherwise. Unknown
/// characters (N, X, etc.) are packed as T (code 0).
pub fn NibSequence::new(seq_id : String, sequence : String) -> NibSequence {
  let packed = nib_encode(sequence)
  let mut upper = true
  let code_a_lower = 'a'.to_int()
  let code_z_lower = 'z'.to_int()
  for i in 0..<sequence.length() {
    let c = sequence[i].to_int()
    if c >= code_a_lower && c <= code_z_lower {
      upper = false
      break
    }
  }
  NibSequence::{
    seq_id,
    length: sequence.length(),
    packed_data: packed,
    is_uppercase: upper,
  }
}

///|
/// Sequence identifier.
pub fn NibSequence::seq_id(self : NibSequence) -> String {
  self.seq_id
}

///|
/// Number of bases in the sequence.
pub fn NibSequence::length(self : NibSequence) -> Int {
  self.length
}

///|
/// A copy of the packed byte array (each element is 0-255).
pub fn NibSequence::packed_data(self : NibSequence) -> Array[Int] {
  self.packed_data.copy()
}

///|
/// Whether the original sequence was entirely uppercase.
pub fn NibSequence::is_uppercase(self : NibSequence) -> Bool {
  self.is_uppercase
}

///|
/// Human-readable summary of the NibSequence.
pub fn NibSequence::to_string(self : NibSequence) -> String {
  let buf = StringBuilder::new()
  buf.write_string("NibSequence(seq_id=")
  buf.write_string(self.seq_id)
  buf.write_string(", length=")
  buf.write_string(self.length.to_string())
  buf.write_string(", packed_bytes=")
  buf.write_string(self.packed_data.length().to_string())
  buf.write_string(", is_uppercase=")
  buf.write_string(self.is_uppercase.to_string())
  buf.write_string(")")
  buf.to_string()
}

// ============================================================================
// High-level operations
// ============================================================================

///|
/// Decode the packed data back into a DNA string, applying the `is_uppercase`
/// flag (lowercasing all bases when false).
pub fn nib_to_string(nib : NibSequence) -> String {
  let decoded = nib_decode(nib.packed_data, nib.length)
  if nib.is_uppercase {
    decoded
  } else {
    let sb = StringBuilder::new()
    for c in decoded {
      let code = c.to_int()
      let code_A = 'A'.to_int()
      let code_Z = 'Z'.to_int()
      if code >= code_A && code <= code_Z {
        sb.write_char(Int::unsafe_to_char(code + 32))
      } else {
        sb.write_char(c)
      }
    }
    sb.to_string()
  }
}

///|
/// Number of bases in the sequence.
pub fn nib_size(nib : NibSequence) -> Int {
  nib.length
}

///|
/// Number of packed bytes (ceil(length / 4)).
pub fn nib_compressed_size(nib : NibSequence) -> Int {
  (nib.length + 3) / 4
}

///|
/// Compression ratio: original size (1 byte per base) divided by packed size
/// (ceil(length/4) bytes). Returns 0.0 for an empty sequence.
pub fn nib_compression_ratio(nib : NibSequence) -> Double {
  let packed = nib_compressed_size(nib)
  if packed == 0 {
    0.0
  } else {
    nib.length.to_double() / packed.to_double()
  }
}

///|
/// Get the single base at `index` as a 1-character string. Returns an empty
/// string if `index` is out of bounds. Respects the `is_uppercase` flag.
pub fn nib_get_base(nib : NibSequence, index : Int) -> String {
  if index < 0 || index >= nib.length {
    return ""
  }
  let byte_idx = index / 4
  let shift = 6 - (index % 4) * 2
  let code = if byte_idx < nib.packed_data.length() {
    (nib.packed_data[byte_idx] >> shift) & 3
  } else {
    0
  }
  let base = nib_code_to_base(code)
  if nib.is_uppercase {
    base.to_string()
  } else {
    let lc = Int::unsafe_to_char(base.to_int() + 32)
    lc.to_string()
  }
}

///|
/// Extract a subsequence of `length` bases starting at `start`. The result is
/// clamped to the available range (so requesting a span past the end returns
/// fewer bases). Respects the `is_uppercase` flag.
pub fn nib_subsequence(nib : NibSequence, start : Int, length : Int) -> String {
  if start < 0 || start >= nib.length || length <= 0 {
    return ""
  }
  let avail = nib.length - start
  let take = if length < avail { length } else { avail }
  let sb = StringBuilder::new()
  for i in 0..<take {
    let pos = start + i
    let byte_idx = pos / 4
    let shift = 6 - (pos % 4) * 2
    let code = if byte_idx < nib.packed_data.length() {
      (nib.packed_data[byte_idx] >> shift) & 3
    } else {
      0
    }
    let base = nib_code_to_base(code)
    if nib.is_uppercase {
      sb.write_char(base)
    } else {
      sb.write_char(Int::unsafe_to_char(base.to_int() + 32))
    }
  }
  sb.to_string()
}

///|
/// Compute the reverse complement of the sequence. A <-> T, C <-> G. The
/// `is_uppercase` flag is preserved. Unknown bases (decoded as 'N') complement
/// to 'N'.
pub fn nib_reverse_complement(nib : NibSequence) -> NibSequence {
  let decoded = nib_decode(nib.packed_data, nib.length)
  let sb = StringBuilder::new()
  for i = nib.length - 1; i >= 0; i = i - 1 {
    let base = decoded[i].unsafe_to_char()
    let comp : Char = match base {
      'A' => 'T'
      'T' => 'A'
      'C' => 'G'
      'G' => 'C'
      _ => 'N'
    }
    if nib.is_uppercase {
      sb.write_char(comp)
    } else {
      sb.write_char(Int::unsafe_to_char(comp.to_int() + 32))
    }
  }
  let rc_str = sb.to_string()
  NibSequence::new(nib.seq_id + "_rc", rc_str)
}

///|
/// Fraction of bases that are G or C (0.0 to 1.0). Returns 0.0 for an empty
/// sequence.
pub fn nib_gc_content(nib : NibSequence) -> Double {
  if nib.length == 0 {
    return 0.0
  }
  let mut gc = 0
  for i in 0..<nib.length {
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    let code = if byte_idx < nib.packed_data.length() {
      (nib.packed_data[byte_idx] >> shift) & 3
    } else {
      0
    }
    // G -> 3, C -> 1
    if code == 3 || code == 1 {
      gc = gc + 1
    }
  }
  gc.to_double() / nib.length.to_double()
}

// ============================================================================
// Hex string I/O
// ============================================================================

///|
/// Convert the packed data to a lowercase hex string (two hex digits per
/// packed byte), suitable for text storage.
pub fn nib_to_hex(nib : NibSequence) -> String {
  let sb = StringBuilder::new()
  for i in 0..<nib.packed_data.length() {
    let b = nib.packed_data[i]
    sb.write_char(nib_hex_digit((b >> 4) & 0xF))
    sb.write_char(nib_hex_digit(b & 0xF))
  }
  sb.to_string()
}

///|
/// Build a NibSequence from a hex string produced by `nib_to_hex`. Each pair
/// of hex characters decodes to one packed byte. The `is_uppercase` flag
/// defaults to true (hex storage does not preserve case information). The
/// `length` argument is the number of bases (which may be less than
/// `4 * hex.length() / 2` when the last byte is only partially used).
pub fn nib_from_hex(seq_id : String, hex : String, length : Int) -> NibSequence {
  let packed : Array[Int] = Array::new()
  let len = hex.length()
  let mut i = 0
  while i + 1 < len {
    let high = nib_hex_value(hex[i].to_int())
    let low = nib_hex_value(hex[i + 1].to_int())
    packed.push(high * 16 + low)
    i = i + 2
  }
  NibSequence::{
    seq_id,
    length,
    packed_data: packed,
    is_uppercase: true,
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a sample NibSequence: id "sample1", sequence "ATGCATGCATGC"
/// (12 bases, 3 packed bytes), all uppercase.
pub fn nib_sample_data() -> NibSequence {
  NibSequence::new("sample1", "ATGCATGCATGC")
}

// ============================================================================
// Private helpers
// ============================================================================

///|
/// Convert an integer 0-15 to a lowercase hex digit character.
fn nib_hex_digit(n : Int) -> Char {
  if n < 10 {
    Int::unsafe_to_char('0'.to_int() + n)
  } else {
    Int::unsafe_to_char('a'.to_int() + (n - 10))
  }
}

///|
/// Convert a hex digit character's ASCII code to its integer value 0-15.
fn nib_hex_value(c : Int) -> Int {
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
