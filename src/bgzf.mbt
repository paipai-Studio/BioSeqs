///|
/// BGZF (Blocked GZIP) decompression module.
/// BGZF is used by BAM format for compressed storage.
/// Each BGZF file consists of multiple gzip-compressed blocks.

///|
/// BGZF/BGZIP parsing errors.
pub suberror BgzfError {
  BgzfError(String)
}

///|
/// BGZF magic number (1f 8b 08 04).
let bgzf_magic : Array[Int] = [0x1f, 0x8b, 0x08, 0x04]

///|
/// GZIP compression method: DEFLATE (8).
let gzip_deflate_method : Int = 8

///|
/// Fixed Huffman code lengths for literals/lengths (0-287).
/// Values 0-143: 8 bits, 144-255: 9 bits, 256-279: 7 bits, 280-287: 8 bits.
let fixed_literal_lengths : FixedArray[Int] = {
  let t = FixedArray::make(288, 0)
  for i in 0..<144 {
    t[i] = 8
  }
  for i in 144..<256 {
    t[i] = 9
  }
  for i in 256..<280 {
    t[i] = 7
  }
  for i in 280..<288 {
    t[i] = 8
  }
  t
}

///|
/// Fixed Huffman code lengths for distances (0-31). All 5 bits.
let fixed_distance_lengths : FixedArray[Int] = {
  let t = FixedArray::make(32, 0)
  for i in 0..<32 {
    t[i] = 5
  }
  t
}

///|
/// Bit reader for DEFLATE stream.
struct BitReader {
  data : String
  pos : Int
  bit_pos : Int
}

///|
/// Create a new bit reader.
fn bit_reader_new(data : String) -> BitReader {
  BitReader::{ data, pos: 0, bit_pos: 0 }
}

///|
/// Read a single bit. DEFLATE packs bits LSB-first within each byte
/// (RFC 1951 section 3.1.1), so bit 0 of a byte is read first.
fn bit_reader_read_bit(br : BitReader) -> (Int, BitReader) {
  let byte = br.data.unsafe_get(br.pos).to_int() & 0xFF
  let bit = (byte >> br.bit_pos) & 1
  let new_bit_pos = br.bit_pos + 1
  let (new_pos, new_bit_pos) = if new_bit_pos >= 8 {
    (br.pos + 1, 0)
  } else {
    (br.pos, new_bit_pos)
  }
  (bit, BitReader::{ data: br.data, pos: new_pos, bit_pos: new_bit_pos })
}

///|
/// Read n bits as a Huffman code (n <= 16). Huffman codes are serialised
/// MSB-first: the first bit received is the most significant bit of the
/// code value, so bits accumulate with the newest bit in the low position.
fn bit_reader_read_bits(br : BitReader, n : Int) -> (Int, BitReader) {
  let mut result = 0
  let mut remaining = n
  let mut reader = br
  while remaining > 0 {
    let (bit, r) = bit_reader_read_bit(reader)
    reader = r
    result = (result << 1) | bit
    remaining = remaining - 1
  }
  (result, reader)
}

///|
/// Read n bits as a plain integer field (BTYPE, HLIT/HDIST/HCLEN, extra
/// length/distance bits, ...). Such fields are packed LSB-first: the first
/// bit received is the least significant bit of the value.
fn bit_reader_read_value(br : BitReader, n : Int) -> (Int, BitReader) {
  let mut result = 0
  let mut bit_index = 0
  let mut reader = br
  while bit_index < n {
    let (bit, r) = bit_reader_read_bit(reader)
    reader = r
    result = result | (bit << bit_index)
    bit_index = bit_index + 1
  }
  (result, reader)
}

///|
/// Skip bits to byte boundary.
fn bit_reader_align(br : BitReader) -> BitReader {
  if br.bit_pos > 0 {
    BitReader::{ data: br.data, pos: br.pos + 1, bit_pos: 0 }
  } else {
    br
  }
}

///|
/// Get remaining bytes in the reader.
fn bit_reader_remaining(br : BitReader) -> Int {
  let bytes_remaining = br.data.length() - br.pos
  if br.bit_pos > 0 {
    bytes_remaining - 1
  } else {
    bytes_remaining
  }
}

///|
/// Huffman decoder using a lookup table.
struct HuffmanDecoder {
  codes : FixedArray[Int]
}

///|
/// Build a Huffman decoder from code lengths.
fn build_huffman_decoder(
  lengths : FixedArray[Int],
  max_symbols : Int,
) -> HuffmanDecoder {
  let bl_count : FixedArray[Int] = FixedArray::make(16, 0)
  let mut i = 0
  while i < max_symbols {
    let len = lengths[i]
    if len > 0 && len < 16 {
      bl_count[len] = bl_count[len] + 1
    }
    i = i + 1
  }

  let next_code : FixedArray[Int] = FixedArray::make(16, 0)
  let mut code = 0
  for bits in 1..<16 {
    code = (code + bl_count[bits - 1]) << 1
    next_code[bits] = code
  }

  let codes : FixedArray[Int] = FixedArray::make(65536, -1)
  let mut i = 0
  while i < max_symbols {
    let len = lengths[i]
    if len > 0 && len <= 16 {
      let c = next_code[len]
      next_code[len] = c + 1
      if len <= 12 {
        let step = 1 << (12 - len)
        let mut j = c
        while j < 4096 {
          codes[j] = i
          j = j + step
        }
      } else {
        codes[c & 0xFFFF] = i
      }
    }
    i = i + 1
  }

  HuffmanDecoder::{ codes, }
}

///|
/// Decode a symbol from bit reader using the Huffman decoder.
fn huffman_decode(decoder : HuffmanDecoder, br : BitReader) -> (Int, BitReader) {
  let (bits12, r1) = bit_reader_read_bits(br, 12)
  if decoder.codes[bits12] != -1 {
    (decoder.codes[bits12], r1)
  } else {
    let (bits14, r2) = bit_reader_read_bits(br, 14)
    if decoder.codes[bits14] != -1 {
      (decoder.codes[bits14], r2)
    } else {
      let (bits16, r3) = bit_reader_read_bits(br, 16)
      (decoder.codes[bits16], r3)
    }
  }
}

///|
/// Length base values for LZ77 matches.
let length_base : FixedArray[Int] = {
  let t = FixedArray::make(30, 0)
  t[0] = 3
  t[1] = 4
  t[2] = 5
  t[3] = 6
  t[4] = 7
  t[5] = 8
  t[6] = 9
  t[7] = 10
  t[8] = 11
  t[9] = 13
  t[10] = 15
  t[11] = 17
  t[12] = 19
  t[13] = 23
  t[14] = 27
  t[15] = 31
  t[16] = 35
  t[17] = 43
  t[18] = 51
  t[19] = 59
  t[20] = 67
  t[21] = 83
  t[22] = 99
  t[23] = 115
  t[24] = 131
  t[25] = 163
  t[26] = 195
  t[27] = 227
  t[28] = 258
  t[29] = 0
  t
}

///|
/// Length extra bits for LZ77 matches.
let length_extra : FixedArray[Int] = {
  let t = FixedArray::make(30, 0)
  t[0] = 0
  t[1] = 0
  t[2] = 0
  t[3] = 0
  t[4] = 0
  t[5] = 0
  t[6] = 0
  t[7] = 0
  t[8] = 1
  t[9] = 1
  t[10] = 1
  t[11] = 1
  t[12] = 2
  t[13] = 2
  t[14] = 2
  t[15] = 2
  t[16] = 3
  t[17] = 3
  t[18] = 3
  t[19] = 3
  t[20] = 4
  t[21] = 4
  t[22] = 4
  t[23] = 4
  t[24] = 5
  t[25] = 5
  t[26] = 5
  t[27] = 5
  t[28] = 0
  t[29] = 0
  t
}

///|
/// Distance base values for LZ77 matches.
let dist_base : FixedArray[Int] = {
  let t = FixedArray::make(30, 0)
  t[0] = 1
  t[1] = 2
  t[2] = 3
  t[3] = 4
  t[4] = 5
  t[5] = 7
  t[6] = 9
  t[7] = 13
  t[8] = 17
  t[9] = 25
  t[10] = 33
  t[11] = 49
  t[12] = 65
  t[13] = 97
  t[14] = 129
  t[15] = 193
  t[16] = 257
  t[17] = 385
  t[18] = 513
  t[19] = 769
  t[20] = 1025
  t[21] = 1537
  t[22] = 2049
  t[23] = 3073
  t[24] = 4097
  t[25] = 6145
  t[26] = 8193
  t[27] = 12289
  t[28] = 16385
  t[29] = 24577
  t
}

///|
/// Distance extra bits for LZ77 matches.
let dist_extra : FixedArray[Int] = {
  let t = FixedArray::make(30, 0)
  t[0] = 0
  t[1] = 0
  t[2] = 0
  t[3] = 0
  t[4] = 1
  t[5] = 1
  t[6] = 2
  t[7] = 2
  t[8] = 3
  t[9] = 3
  t[10] = 4
  t[11] = 4
  t[12] = 5
  t[13] = 5
  t[14] = 6
  t[15] = 6
  t[16] = 7
  t[17] = 7
  t[18] = 8
  t[19] = 8
  t[20] = 9
  t[21] = 9
  t[22] = 10
  t[23] = 10
  t[24] = 11
  t[25] = 11
  t[26] = 12
  t[27] = 12
  t[28] = 13
  t[29] = 13
  t
}

///|
/// Decode raw DEFLATE data.
fn deflate_decode(data : String, max_output : Int) -> String raise BgzfError {
  let output : FixedArray[UInt16] = FixedArray::make(max_output, 0)
  let mut reader = bit_reader_new(data)
  let mut out_pos = 0

  let mut done = false
  while !done {
    if bit_reader_remaining(reader) == 0 {
      done = true
    } else {
      let (b_final, r1) = bit_reader_read_bit(reader)
      let (b_type, r2) = bit_reader_read_value(r1, 2)
      reader = r2

      if b_type == 0 {
        // Stored block: after byte alignment, LEN/NLEN and the raw bytes
        // are plain little-endian/raw bytes (RFC 1951 section 3.2.4), not
        // bit-packed values, so they must be read directly from the stream.
        reader = bit_reader_align(reader)
        let p = reader.pos
        if p + 4 > data.length() {
          raise BgzfError("Truncated stored DEFLATE block header")
        }
        let len = (data.unsafe_get(p).to_int() & 0xFF) |
          ((data.unsafe_get(p + 1).to_int() & 0xFF) << 8)
        let nlen = (data.unsafe_get(p + 2).to_int() & 0xFF) |
          ((data.unsafe_get(p + 3).to_int() & 0xFF) << 8)
        if (len ^ 0xFFFF) != nlen {
          raise BgzfError("Invalid stored block length")
        }
        if p + 4 + len > data.length() || out_pos + len > max_output {
          raise BgzfError("Truncated stored DEFLATE block data")
        }
        let mut i = 0
        while i < len {
          output[out_pos + i] = (data.unsafe_get(p + 4 + i).to_int() & 0xFF)
            .to_uint16()
          i = i + 1
        }
        reader = BitReader::{ data: reader.data, pos: p + 4 + len, bit_pos: 0 }
        out_pos = out_pos + len
      } else if b_type == 1 {
        let lit_decoder = build_huffman_decoder(fixed_literal_lengths, 288)
        let dist_decoder = build_huffman_decoder(fixed_distance_lengths, 32)
        let mut block_done = false
        while !block_done {
          let (sym, r3) = huffman_decode(lit_decoder, reader)
          reader = r3
          if sym == 256 {
            block_done = true
          } else if sym < 256 {
            output[out_pos] = sym.to_uint16()
            out_pos = out_pos + 1
          } else {
            let length_code = sym - 257
            let extra_len_bits = length_extra[length_code]
            let (extra_len, r4) = bit_reader_read_value(reader, extra_len_bits)
            reader = r4
            let length = length_base[length_code] + extra_len

            let (dist_code, r5) = huffman_decode(dist_decoder, reader)
            reader = r5
            let extra_dist_bits = dist_extra[dist_code]
            let (extra_dist, r6) = bit_reader_read_value(reader, extra_dist_bits)
            reader = r6
            let distance = dist_base[dist_code] + extra_dist

            let mut i = 0
            while i < length {
              output[out_pos] = output[out_pos - distance]
              out_pos = out_pos + 1
              i = i + 1
            }
          }
        }
      } else if b_type == 2 {
        let (num_lit_len, r3) = bit_reader_read_value(reader, 5)
        reader = r3
        let num_lit_len = num_lit_len + 257
        let (num_dist, r4) = bit_reader_read_value(reader, 5)
        reader = r4
        let num_dist = num_dist + 1

        let (num_code_len, r5) = bit_reader_read_value(reader, 4)
        reader = r5
        let num_code_len = num_code_len + 4

        let code_length_lengths : FixedArray[Int] = FixedArray::make(19, 0)
        let code_order = [
          16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
        ]
        let mut i = 0
        while i < num_code_len {
          let (len, r6) = bit_reader_read_value(reader, 3)
          reader = r6
          code_length_lengths[code_order[i]] = len
          i = i + 1
        }

        let code_len_decoder = build_huffman_decoder(code_length_lengths, 19)

        let lit_len_lengths : FixedArray[Int] = FixedArray::make(288, 0)
        let mut j = 0
        while j < num_lit_len {
          let (sym, r7) = huffman_decode(code_len_decoder, reader)
          reader = r7
          if sym < 16 {
            lit_len_lengths[j] = sym
            j = j + 1
          } else if sym == 16 {
            let (rep, r8) = bit_reader_read_value(reader, 2)
            reader = r8
            let rep = rep + 3
            let val = if j > 0 { lit_len_lengths[j - 1] } else { 0 }
            let mut k = 0
            while k < rep {
              lit_len_lengths[j] = val
              j = j + 1
              k = k + 1
            }
          } else if sym == 17 {
            let (zeros, r8) = bit_reader_read_value(reader, 3)
            reader = r8
            let zeros = zeros + 3
            let mut k = 0
            while k < zeros {
              lit_len_lengths[j] = 0
              j = j + 1
              k = k + 1
            }
          } else if sym == 18 {
            let (zeros, r8) = bit_reader_read_value(reader, 7)
            reader = r8
            let zeros = zeros + 11
            let mut k = 0
            while k < zeros {
              lit_len_lengths[j] = 0
              j = j + 1
              k = k + 1
            }
          }
        }

        let dist_lengths : FixedArray[Int] = FixedArray::make(32, 0)
        j = 0
        while j < num_dist {
          let (sym, r7) = huffman_decode(code_len_decoder, reader)
          reader = r7
          if sym < 16 {
            dist_lengths[j] = sym
            j = j + 1
          } else if sym == 16 {
            let (rep, r8) = bit_reader_read_value(reader, 2)
            reader = r8
            let rep = rep + 3
            let val = if j > 0 { dist_lengths[j - 1] } else { 0 }
            let mut k = 0
            while k < rep {
              dist_lengths[j] = val
              j = j + 1
              k = k + 1
            }
          } else if sym == 17 {
            let (zeros, r8) = bit_reader_read_value(reader, 3)
            reader = r8
            let zeros = zeros + 3
            let mut k = 0
            while k < zeros {
              dist_lengths[j] = 0
              j = j + 1
              k = k + 1
            }
          } else if sym == 18 {
            let (zeros, r8) = bit_reader_read_value(reader, 7)
            reader = r8
            let zeros = zeros + 11
            let mut k = 0
            while k < zeros {
              dist_lengths[j] = 0
              j = j + 1
              k = k + 1
            }
          }
        }

        let lit_decoder = build_huffman_decoder(lit_len_lengths, 288)
        let dist_decoder = build_huffman_decoder(dist_lengths, 32)

        let mut block_done = false
        while !block_done {
          let (sym, r9) = huffman_decode(lit_decoder, reader)
          reader = r9
          if sym == 256 {
            block_done = true
          } else if sym < 256 {
            output[out_pos] = sym.to_uint16()
            out_pos = out_pos + 1
          } else {
            let length_code = sym - 257
            let extra_len_bits = length_extra[length_code]
            let (extra_len, r10) = bit_reader_read_value(reader, extra_len_bits)
            reader = r10
            let length = length_base[length_code] + extra_len

            let (dist_code, r11) = huffman_decode(dist_decoder, reader)
            reader = r11
            let extra_dist_bits = dist_extra[dist_code]
            let (extra_dist, r12) = bit_reader_read_value(
              reader, extra_dist_bits,
            )
            reader = r12
            let distance = dist_base[dist_code] + extra_dist

            let mut k = 0
            while k < length {
              output[out_pos] = output[out_pos - distance]
              out_pos = out_pos + 1
              k = k + 1
            }
          }
        }
      } else {
        raise BgzfError("Invalid DEFLATE block type")
      }

      if b_final == 1 {
        done = true
      }
    }
  }

  unsafe_from_uint16_fixedarray(output)
}

///|
/// Read a UInt32 (little-endian) from data at position.
fn get_uint32_le(data : String, pos : Int) -> Int {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)) & 0xFFFFFFFF
}

///|
/// Parse and decompress a single BGZF block.
/// Returns (decompressed_data, next_block_pos).
fn decompress_bgzf_block(
  data : String,
  pos : Int,
) -> (String, Int) raise BgzfError {
  let n = data.length()
  if pos + 18 > n {
    raise BgzfError("BGZF block too short")
  }

  for i in 0..<4 {
    if data.unsafe_get(pos + i).to_int() != bgzf_magic[i] {
      raise BgzfError("Invalid BGZF magic number")
    }
  }

  // gzip member layout: [0..4) magic (1f 8b), [2] CM=8, [3] FLG,
  // [4..8) MTIME, [8] XFL, [9] OS, [10..12) XLEN (if FEXTRA), ...
  // The 4-byte BGZF magic already covers CM and FLG; check CM explicitly.
  let method = data.unsafe_get(pos + 2).to_int() & 0xFF
  if method != gzip_deflate_method {
    raise BgzfError("Unsupported compression method")
  }

  let flags = data.unsafe_get(pos + 3).to_int() & 0xFF
  let has_extra = (flags & 0x04) != 0

  let mut current_pos = pos + 10

  if has_extra {
    let extra_len = (data.unsafe_get(current_pos).to_int() & 0xFF) |
      ((data.unsafe_get(current_pos + 1).to_int() & 0xFF) << 8)
    current_pos = current_pos + 2
    current_pos = current_pos + extra_len
  }

  if (flags & 0x08) != 0 {
    while current_pos < n && data.unsafe_get(current_pos) != 0 {
      current_pos = current_pos + 1
    }
    current_pos = current_pos + 1
  }

  if (flags & 0x10) != 0 {
    while current_pos < n && data.unsafe_get(current_pos) != 0 {
      current_pos = current_pos + 1
    }
    current_pos = current_pos + 1
  }

  if (flags & 0x20) != 0 {
    while current_pos < n && data.unsafe_get(current_pos) != 0 {
      current_pos = current_pos + 1
    }
    current_pos = current_pos + 1
  }

  // BSIZE (at offset 16, little-endian u16) holds the total block size minus
  // one. The 8-byte gzip trailer (CRC32 + ISIZE) occupies the last 8 bytes.
  let bsize = (data.unsafe_get(pos + 16).to_int() & 0xFF) |
    ((data.unsafe_get(pos + 17).to_int() & 0xFF) << 8)
  let block_end = pos + bsize + 1
  let trailer_pos = block_end - 8

  if block_end > n || trailer_pos < current_pos {
    raise BgzfError("Invalid BGZF block size")
  }

  let compressed_data = data[current_pos:trailer_pos].to_owned()

  let stored_crc = get_uint32_le(data, trailer_pos)
  let stored_size = get_uint32_le(data, trailer_pos + 4)

  let decompressed = deflate_decode(compressed_data, stored_size)

  let computed_crc = crc32(decompressed)
  if computed_crc != stored_crc {
    raise BgzfError("CRC32 check failed")
  }

  if decompressed.length() != stored_size {
    raise BgzfError("Decompressed size mismatch")
  }

  (decompressed, trailer_pos + 8)
}

///|
/// Decompress entire BGZF data.
pub fn bgzf_decompress(data : String) -> String raise BgzfError {
  let result = StringBuilder::new()
  let mut pos = 0
  let n = data.length()

  while pos < n - 1 {
    if pos + 4 > n {
      break
    }

    let mut is_gzip = true
    for i in 0..<4 {
      if data.unsafe_get(pos + i).to_int() != bgzf_magic[i] {
        is_gzip = false
        break
      }
    }

    if !is_gzip {
      result.write_string(data[pos:n].to_owned())
      break
    }

    let (block_data, new_pos) = decompress_bgzf_block(data, pos)
    result.write_string(block_data)
    pos = new_pos
  }

  result.to_string()
}

///|
/// Parse BAM file from BGZF-compressed data.
pub fn parse_bam_from_bgzf(data : String) -> BamFile raise BgzfError {
  let decompressed = bgzf_decompress(data)
  let result = parse_bam(decompressed) catch {
    BamError(msg) => raise BgzfError(msg)
  }
  result
}

///|
/// Parse BAM file from base64-encoded BGZF-compressed data.
pub fn parse_bam_from_bgzf_base64(
  base64_str : String,
) -> BamFile raise BgzfError {
  let raw_data = base64_decode(base64_str)
  parse_bam_from_bgzf(raw_data)
}

// ============================================================================
// BGZF writing and random-access reading.
// Port of Biopython's Bio.bgzf BgzfWriter / BgzfReader / BgzfBlocks and the
// virtual offset (coffset/uoffset) machinery. Blocks are written as stored
// (uncompressed) DEFLATE blocks, which is fully spec-compliant BGZF and can
// be read by any conformant reader (samtools bgzip, zlib, etc.).
// ============================================================================

///|
/// Maximum uncompressed payload per BGZF block. Each on-disk block is at
/// most 65536 bytes; a stored DEFLATE block costs 31 bytes of overhead
/// (18-byte gzip header with BC extra field + 5 stored-block header +
/// 8-byte trailer).
let bgzf_max_payload : Int = 65536 - 31 // 65505

///|
/// Combine a compressed block start offset and an uncompressed offset
/// within that block into a 64-bit BGZF virtual offset (as used by BAM/
/// tabix indexing). Mirrors Biopython's `make_virtual_offset`.
pub fn bgzf_make_virtual_offset(block_offset : Int, within_block : Int) -> Int64 {
  (block_offset.to_int64() << 16) | (within_block & 0xFFFF).to_int64()
}

///|
/// Split a virtual offset into (compressed block start offset, uncompressed
/// offset within block). Mirrors Biopython's `split_virtual_offset`.
pub fn bgzf_split_virtual_offset(virtual_offset : Int64) -> (Int, Int) {
  let within = (virtual_offset & 0xFFFFL).to_int()
  let block = (virtual_offset >> 16).to_int()
  (block, within)
}

///|
/// Build one BGZF block around an uncompressed payload, using a single
/// final stored DEFLATE block.
fn bgzf_build_block(payload : String) -> String {
  let n = payload.length()
  let total = 31 + n
  let bsize = total - 1
  let nlen = n ^ 0xFFFF
  let crc = crc32(payload) & 0xFFFFFFFF
  let isize = n & 0xFFFFFFFF

  let buf : FixedArray[UInt16] = FixedArray::make(total, 0)

  // gzip member header with FEXTRA + BC subfield
  buf[0] = (0x1f).to_uint16() // ID1
  buf[1] = (0x8b).to_uint16() // ID2
  buf[2] = (0x08).to_uint16() // CM = DEFLATE
  buf[3] = (0x04).to_uint16() // FLG = FEXTRA
  // MTIME 4 bytes = 0, XFL = 0, OS = 0xff
  buf[9] = (0xff).to_uint16()
  buf[10] = (0x06).to_uint16() // XLEN = 6 (lo)
  buf[11] = (0x00).to_uint16() // XLEN (hi)
  buf[12] = (0x42).to_uint16() // SI1 = 'B'
  buf[13] = (0x43).to_uint16() // SI2 = 'C'
  buf[14] = (0x02).to_uint16() // SLEN = 2 (lo)
  buf[15] = (0x00).to_uint16() // SLEN (hi)
  buf[16] = (bsize & 0xff).to_uint16() // BSIZE (lo)
  buf[17] = ((bsize >> 8) & 0xff).to_uint16() // BSIZE (hi)

  // stored DEFLATE block: BFINAL=1, BTYPE=00
  buf[18] = (0x01).to_uint16()
  buf[19] = (n & 0xff).to_uint16() // LEN (lo)
  buf[20] = ((n >> 8) & 0xff).to_uint16() // LEN (hi)
  buf[21] = (nlen & 0xff).to_uint16() // NLEN (lo)
  buf[22] = ((nlen >> 8) & 0xff).to_uint16() // NLEN (hi)

  for i in 0..<n {
    buf[23 + i] = (payload.unsafe_get(i).to_int() & 0xFF).to_uint16()
  }

  // trailer: CRC32 then ISIZE, both little-endian
  let t = 23 + n
  buf[t] = (crc & 0xFF).to_uint16()
  buf[t + 1] = ((crc >> 8) & 0xFF).to_uint16()
  buf[t + 2] = ((crc >> 16) & 0xFF).to_uint16()
  buf[t + 3] = ((crc >> 24) & 0xFF).to_uint16()
  buf[t + 4] = (isize & 0xFF).to_uint16()
  buf[t + 5] = ((isize >> 8) & 0xFF).to_uint16()
  buf[t + 6] = ((isize >> 16) & 0xFF).to_uint16()
  buf[t + 7] = ((isize >> 24) & 0xFF).to_uint16()

  unsafe_from_uint16_fixedarray(buf)
}

///|
/// The canonical 28-byte empty BGZF block that marks end-of-file.
fn bgzf_eof_block() -> String {
  let bytes = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  ]
  let buf : FixedArray[UInt16] = FixedArray::make(28, 0)
  for i in 0..<28 {
    buf[i] = bytes[i].to_uint16()
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// In-memory BGZF writer. Accumulates uncompressed data and emits one
/// BGZF block per at-most-`bgzf_max_payload` chunk, followed by the
/// canonical 28-byte EOF block on close. Mirrors Biopython's BgzfWriter.
pub struct BgzfWriter {
  priv mut buffer : String
  priv output : StringBuilder
  priv mut block_offset : Int
  priv mut closed : Bool
}

///|
pub fn BgzfWriter::new() -> BgzfWriter {
  BgzfWriter::{
    buffer: "",
    output: StringBuilder::new(),
    block_offset: 0,
    closed: false,
  }
}

///|
/// Emit the current buffer as one BGZF block.
fn BgzfWriter::flush_block(self : BgzfWriter) -> Unit {
  if self.buffer.length() == 0 {
    return
  }
  let block = bgzf_build_block(self.buffer)
  self.output.write_string(block)
  self.block_offset = self.block_offset + block.length()
  self.buffer = ""
}

///|
/// Append uncompressed data, flushing full blocks as needed.
pub fn BgzfWriter::write(
  self : BgzfWriter,
  data : String,
) -> Unit raise BgzfError {
  if self.closed {
    raise BgzfError("write to closed BGZF writer")
  }
  self.buffer = self.buffer + data
  while self.buffer.length() > bgzf_max_payload {
    let payload = self.buffer[0:bgzf_max_payload].to_owned()
    self.buffer = self.buffer[bgzf_max_payload:].to_owned()
    let block = bgzf_build_block(payload)
    self.output.write_string(block)
    self.block_offset = self.block_offset + block.length()
  }
}

///|
/// Flush any remaining data and append the EOF marker block.
pub fn BgzfWriter::close(self : BgzfWriter) -> Unit {
  if self.closed {
    return
  }
  self.flush_block()
  let eof = bgzf_eof_block()
  self.output.write_string(eof)
  self.block_offset = self.block_offset + eof.length()
  self.closed = true
}

///|
/// Close the writer and return the complete BGZF byte stream.
pub fn BgzfWriter::getvalue(self : BgzfWriter) -> String {
  self.close()
  self.output.to_string()
}

///|
/// Return the virtual offset at which the next byte would be written:
/// the flushed file position (block start) plus the buffered byte count.
pub fn BgzfWriter::tell(self : BgzfWriter) -> Int64 {
  bgzf_make_virtual_offset(self.block_offset, self.buffer.length())
}

///|
/// One-shot BGZF compression: data in, a complete BGZF stream (with EOF
/// marker) out. Mirrors `bgzf.compress`-style usage.
pub fn bgzf_compress(data : String) -> String raise BgzfError {
  let writer = BgzfWriter::new()
  writer.write(data)
  writer.getvalue()
}

///|
/// Metadata for one on-disk BGZF block, as reported by `bgzf_blocks`.
/// Mirrors Biopython's BgzfBlocks tuples:
/// (start_offset, block_length, data_offset, data_length).
pub struct BgzfBlockInfo {
  start_offset : Int
  block_length : Int
  data_offset : Int
  data_length : Int
} derive(Debug)

///|
/// Scan a BGZF stream and describe every block without decompressing.
pub fn bgzf_blocks(data : String) -> Array[BgzfBlockInfo] raise BgzfError {
  let blocks : Array[BgzfBlockInfo] = Array::new()
  let n = data.length()
  let mut pos = 0
  let mut data_offset = 0
  while pos < n {
    if pos + 18 > n {
      raise BgzfError(
        "Truncated BGZF block at offset " + pos.to_string(),
      )
    }
    for i in 0..<4 {
      if data.unsafe_get(pos + i).to_int() != bgzf_magic[i] {
        raise BgzfError(
          "Invalid BGZF magic at offset " + pos.to_string(),
        )
      }
    }
    let bsize = (data.unsafe_get(pos + 16).to_int() & 0xFF) |
      ((data.unsafe_get(pos + 17).to_int() & 0xFF) << 8)
    let block_len = bsize + 1
    if pos + block_len > n {
      raise BgzfError(
        "BGZF block at offset " +
        pos.to_string() +
        " exceeds file bounds",
      )
    }
    let isize = get_uint32_le(data, pos + block_len - 4)
    blocks.push(
      BgzfBlockInfo::{
        start_offset: pos,
        block_length: block_len,
        data_offset,
        data_length: isize,
      },
    )
    data_offset = data_offset + isize
    pos = pos + block_len
    if isize == 0 {
      // Canonical empty EOF block; nothing can follow meaningfully.
      break
    }
  }
  blocks
}

///|
/// Random-access reader over a BGZF byte stream. Supports virtual-offset
/// seek/tell (as used by BAM/tabix indexing) plus read/readline.
/// Mirrors Biopython's BgzfReader (string-backed, like all of BioSeqs).
pub struct BgzfReader {
  priv data : String
  priv blocks : Array[BgzfBlockInfo]
  priv mut block_index : Int
  priv mut within : Int
  priv mut cached_index : Int
  priv mut cached_text : String
}

///|
pub fn BgzfReader::new(data : String) -> BgzfReader raise BgzfError {
  let blocks = bgzf_blocks(data)
  BgzfReader::{
    data,
    blocks,
    block_index: 0,
    within: 0,
    cached_index: -1,
    cached_text: "",
  }
}

///|
/// Decompress (and cache) the block currently pointed at.
fn BgzfReader::decompressed_current(self : BgzfReader) -> String raise BgzfError {
  if self.cached_index != self.block_index {
    let info = self.blocks[self.block_index]
    let (text, _) = decompress_bgzf_block(self.data, info.start_offset)
    self.cached_index = self.block_index
    self.cached_text = text
  }
  self.cached_text
}

///|
/// Current virtual offset (block start << 16 | offset within block).
pub fn BgzfReader::tell(self : BgzfReader) -> Int64 {
  if self.block_index >= self.blocks.length() {
    // Positioned at end: report the EOF block start.
    let last = self.blocks[self.blocks.length() - 1]
    return bgzf_make_virtual_offset(last.start_offset, 0)
  }
  let info = self.blocks[self.block_index]
  bgzf_make_virtual_offset(info.start_offset, self.within)
}

///|
/// Seek to a virtual offset previously obtained from `tell` or an index.
pub fn BgzfReader::seek(
  self : BgzfReader,
  virtual_offset : Int64,
) -> Unit raise BgzfError {
  let (coffset, uoffset) = bgzf_split_virtual_offset(virtual_offset)
  let mut found = -1
  let mut i = 0
  while i < self.blocks.length() {
    if self.blocks[i].start_offset == coffset {
      found = i
      i = self.blocks.length()
    }
    i = i + 1
  }
  if found < 0 {
    raise BgzfError(
      "No BGZF block starts at compressed offset " + coffset.to_string(),
    )
  }
  if uoffset > self.blocks[found].data_length {
    raise BgzfError(
      "Uncompressed offset " +
      uoffset.to_string() +
      " exceeds block length " +
      self.blocks[found].data_length.to_string(),
    )
  }
  self.block_index = found
  self.within = uoffset
}

///|
/// Read up to `size` bytes (or all remaining data if omitted).
pub fn BgzfReader::read(self : BgzfReader, size? : Int) -> String raise BgzfError {
  let result = StringBuilder::new()
  let mut remaining = size
  while true {
    if self.block_index >= self.blocks.length() {
      break
    }
    let info = self.blocks[self.block_index]
    if info.data_length == 0 {
      // Empty EOF block.
      break
    }
    if remaining is Some(0) {
      break
    }
    if self.within >= info.data_length {
      self.block_index = self.block_index + 1
      self.within = 0
      continue
    }
    let text = self.decompressed_current()
    let avail = info.data_length - self.within
    let take = match remaining {
      Some(r) => if r < avail { r } else { avail }
      None => avail
    }
    result.write_string(text[self.within:self.within + take].to_owned())
    self.within = self.within + take
    if remaining is Some(r) {
      remaining = Some(r - take)
    }
    if self.within >= info.data_length {
      self.block_index = self.block_index + 1
      self.within = 0
    }
  }
  result.to_string()
}

///|
/// Read one line (including the trailing newline when present).
pub fn BgzfReader::readline(self : BgzfReader) -> String raise BgzfError {
  let result = StringBuilder::new()
  while true {
    if self.block_index >= self.blocks.length() {
      break
    }
    let info = self.blocks[self.block_index]
    if info.data_length == 0 {
      break
    }
    let text = self.decompressed_current()
    // Search for newline within the current block remainder.
    let mut nl = -1
    let mut j = self.within
    while j < info.data_length {
      if text.unsafe_get(j).to_int() == '\n'.to_int() {
        nl = j
        j = info.data_length
      }
      j = j + 1
    }
    if nl >= 0 {
      result.write_string(text[self.within:nl + 1].to_owned())
      self.within = nl + 1
      if self.within >= info.data_length {
        self.block_index = self.block_index + 1
        self.within = 0
      }
      break
    } else {
      result.write_string(text[self.within:info.data_length].to_owned())
      self.block_index = self.block_index + 1
      self.within = 0
    }
  }
  result.to_string()
}
