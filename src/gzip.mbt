///|
/// Standard gzip (RFC 1952) decompression, reusing the DEFLATE decoder from
/// `bgzf.mbt`. Handles a single-member stream (the usual `.gz` file layout).

///|
pub fn gzip_decompress(data : String) -> String raise BgzfError {
  let n = data.length()
  if n < 18 {
    raise BgzfError("gzip data too short")
  }
  if data.unsafe_get(0).to_int() != 0x1F || data.unsafe_get(1).to_int() != 0x8B {
    raise BgzfError("Invalid gzip magic number")
  }
  if (data.unsafe_get(2).to_int() & 0xFF) != 8 {
    raise BgzfError("Unsupported gzip compression method")
  }
  let flags = data.unsafe_get(3).to_int() & 0xFF
  let mut pos = 10
  if (flags & 0x04) != 0 {
    let xlen = (data.unsafe_get(pos).to_int() & 0xFF) |
      ((data.unsafe_get(pos + 1).to_int() & 0xFF) << 8)
    pos = pos + 2 + xlen
  }
  if (flags & 0x08) != 0 {
    while pos < n && data.unsafe_get(pos) != 0 {
      pos = pos + 1
    }
    pos = pos + 1
  }
  if (flags & 0x10) != 0 {
    while pos < n && data.unsafe_get(pos) != 0 {
      pos = pos + 1
    }
    pos = pos + 1
  }
  if (flags & 0x02) != 0 {
    pos = pos + 2
  }
  let trailer_pos = n - 8
  if trailer_pos < pos {
    raise BgzfError("Invalid gzip trailer")
  }
  let isize = get_uint32_le(data, trailer_pos + 4)
  let compressed = data[pos:trailer_pos].to_owned()
  deflate_decode(compressed, isize)
}

///|
/// Decompress `data` if it begins with the gzip magic number; return it
/// unchanged otherwise.
pub fn gunzip_if_needed(data : String) -> String raise BgzfError {
  if data.length() >= 2 &&
    data.unsafe_get(0).to_int() == 0x1F &&
    data.unsafe_get(1).to_int() == 0x8B {
    gzip_decompress(data)
  } else {
    data
  }
}