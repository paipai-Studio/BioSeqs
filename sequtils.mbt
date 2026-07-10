///|
/// Bio.SeqUtils - Utility functions for sequence analysis.

///|
/// CheckSum module
let crc32_table : Array[Int] = make_crc32_table()

///|
fn make_crc32_table() -> Array[Int] {
  let table = Array::new()
  let mut i = 0
  while i < 256 {
    let mut c = i & 0xFFFFFFFF
    let mut j = 0
    while j < 8 {
      if (c & 1) != 0 {
        c = ((c >> 1) & 0x7FFFFFFF) ^ 0xEDB88320
      } else {
        c = (c >> 1) & 0x7FFFFFFF
      }
      j = j + 1
    }
    table.push(c & 0xFFFFFFFF)
    i = i + 1
  }
  table
}

///|
pub fn crc32(seq : String) -> Int {
  let mut crc = 0xFFFFFFFF
  let mut i = 0
  while i < seq.length() {
    match seq.get_char(i) {
      Some(c) => {
        let byte = c.to_int() & 0xFF
        let index = (crc ^ byte) & 0xFF
        crc = crc32_table[index] ^ ((crc >> 8) & 0x00FFFFFF)
      }
      None => ()
    }
    i = i + 1
  }
  (crc ^ 0xFFFFFFFF) | 0
}

///|
pub fn gcg(seq : String) -> Int {
  let mut index = 0
  let mut checksum = 0
  let mut i = 0
  while i < seq.length() {
    match seq.get_char(i) {
      Some(c) => {
        let code = c.to_int()
        let upper_code = if code >= 97 && code <= 122 {
          code - 32
        } else {
          code
        }
        index = index + 1
        checksum = checksum + index * upper_code
        if index == 57 {
          index = 0
        }
      }
      None => ()
    }
    i = i + 1
  }
  checksum % 10000
}

///|
pub fn seguid(seq : String) -> String {
  let seq_upper = seq.to_upper()
  let hash_result = sha1_string(seq_upper)
  let chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  let buf = StringBuilder::new()
  let mut i = 0
  while i < hash_result.length() {
    let byte = hash_result[i]
    let idx1 = byte >> 2
    match chars.get_char(idx1) {
      Some(ch) => buf.write_char(ch)
      None => ()
    }
    let rem = byte & 0x03
    i = i + 1
    if i < hash_result.length() {
      let byte2 = hash_result[i]
      let idx2 = (rem << 4) | (byte2 >> 4)
      match chars.get_char(idx2) {
        Some(ch) => buf.write_char(ch)
        None => ()
      }
      let rem2 = byte2 & 0x0F
      i = i + 1
      if i < hash_result.length() {
        let byte3 = hash_result[i]
        let idx3 = (rem2 << 2) | (byte3 >> 6)
        match chars.get_char(idx3) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
        let idx4 = byte3 & 0x3F
        match chars.get_char(idx4) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
        i = i + 1
      } else {
        let idx3 = rem2 << 2
        match chars.get_char(idx3) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
      }
    } else {
      let idx2 = rem << 4
      match chars.get_char(idx2) {
        Some(ch) => buf.write_char(ch)
        None => ()
      }
    }
  }
  buf.to_string()
}

///|
/// SHA1 implementation
pub fn sha1_string(s : String) -> Array[Int] {
  let bytes = string_to_bytes(s)
  sha1(bytes)
}

///|
pub fn string_to_bytes(s : String) -> Array[Int] {
  let result = Array::new()
  let mut i = 0
  while i < s.length() {
    match s.get_char(i) {
      Some(c) => result.push(c.to_int() & 0xFF)
      None => ()
    }
    i = i + 1
  }
  result
}

///|
pub fn bytes_to_base64(bytes : Array[Int]) -> String {
  let chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  let buf = StringBuilder::new()
  let mut i = 0
  while i < bytes.length() {
    let byte = bytes[i]
    let idx1 = byte >> 2
    match chars.get_char(idx1) {
      Some(ch) => buf.write_char(ch)
      None => ()
    }
    let rem = byte & 0x03
    i = i + 1
    if i < bytes.length() {
      let byte2 = bytes[i]
      let idx2 = (rem << 4) | (byte2 >> 4)
      match chars.get_char(idx2) {
        Some(ch) => buf.write_char(ch)
        None => ()
      }
      let rem2 = byte2 & 0x0F
      i = i + 1
      if i < bytes.length() {
        let byte3 = bytes[i]
        let idx3 = (rem2 << 2) | (byte3 >> 6)
        match chars.get_char(idx3) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
        let idx4 = byte3 & 0x3F
        match chars.get_char(idx4) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
        i = i + 1
      } else {
        let idx3 = rem2 << 2
        match chars.get_char(idx3) {
          Some(ch) => buf.write_char(ch)
          None => ()
        }
      }
    } else {
      let idx2 = rem << 4
      match chars.get_char(idx2) {
        Some(ch) => buf.write_char(ch)
        None => ()
      }
    }
  }
  buf.to_string()
}

///|
fn u32(x : Int) -> Int {
  x & 0xFFFFFFFF
}

///|
fn u32_shr(a : Int, shift : Int) -> Int {
  let mask = 0xFFFFFFFF
  let masked = a & mask
  let mut result = 0
  let mut i = 0
  while i < 32 {
    if (masked & (1 << i)) != 0 {
      let pos = i - shift
      if pos >= 0 {
        result = result | (1 << pos)
      }
    }
    i = i + 1
  }
  result & mask
}

///|
fn rotl5(x : Int) -> Int {
  u32((x << 5) | u32_shr(x, 27))
}

///|
fn rotl30(x : Int) -> Int {
  u32((x << 30) | u32_shr(x, 2))
}

///|
fn rotl1(x : Int) -> Int {
  u32((x << 1) | u32_shr(x, 31))
}

///|
pub fn sha1(message : Array[Int]) -> Array[Int] {
  let mut h0 = 0x67452301
  let mut h1 = 0xEFCDAB89
  let mut h2 = 0x98BADCFE
  let mut h3 = 0x10325476
  let mut h4 = 0xC3D2E1F0

  let original_len = message.length()
  let padded_len = ((original_len + 8) / 64 + 1) * 64
  let padded = Array::new()
  let mut i = 0
  while i < padded_len {
    if i < original_len {
      padded.push(message[i])
    } else if i == original_len {
      padded.push(128)
    } else if i < padded_len - 8 {
      padded.push(0)
    } else {
      let len_bits = original_len * 8
      let byte_idx = i - (padded_len - 8)
      padded.push(u32_shr(len_bits, (7 - byte_idx) * 8) & 255)
    }
    i = i + 1
  }

  let mut block_idx = 0
  while block_idx < padded.length() {
    let w = Array::new()
    let mut j = 0
    while j < 80 {
      if j < 16 {
        let idx = block_idx + j * 4
        let b0 = padded[idx] & 255
        let b1 = padded[idx + 1] & 255
        let b2 = padded[idx + 2] & 255
        let b3 = padded[idx + 3] & 255
        w.push(u32((b0 << 24) | (b1 << 16) | (b2 << 8) | b3))
      } else {
        let s = u32(w[j - 3] ^ w[j - 8] ^ (w[j - 14] ^ w[j - 16]))
        w.push(rotl1(s))
      }
      j = j + 1
    }

    let mut a = h0
    let mut b = h1
    let mut c = h2
    let mut d = h3
    let mut e = h4

    let k1 = 0x5A827999
    let k2 = 0x6ED9EBA1
    let k3 = 0x8F1BBCDC
    let k4 = 0xCA62C1D6

    j = 0
    while j < 80 {
      let mut f = 0
      let mut k_val = 0
      if j < 20 {
        f = u32((b & c) | ((b ^ 0xFFFFFFFF) & d))
        k_val = k1
      } else if j < 40 {
        f = u32(b ^ c ^ d)
        k_val = k2
      } else if j < 60 {
        f = u32((b & c) | (b & d) | (c & d))
        k_val = k3
      } else {
        f = u32(b ^ c ^ d)
        k_val = k4
      }

      let wj = w[j]
      let temp = u32(rotl5(a) + f + e + k_val + wj)

      e = d
      d = c
      c = rotl30(b)
      b = a
      a = temp
      j = j + 1
    }

    h0 = u32(h0 + a)
    h1 = u32(h1 + b)
    h2 = u32(h2 + c)
    h3 = u32(h3 + d)
    h4 = u32(h4 + e)

    block_idx = block_idx + 64
  }

  let result = Array::new()
  result.push(u32_shr(h0, 24) & 255)
  result.push(u32_shr(h0, 16) & 255)
  result.push(u32_shr(h0, 8) & 255)
  result.push(h0 & 255)
  result.push(u32_shr(h1, 24) & 255)
  result.push(u32_shr(h1, 16) & 255)
  result.push(u32_shr(h1, 8) & 255)
  result.push(h1 & 255)
  result.push(u32_shr(h2, 24) & 255)
  result.push(u32_shr(h2, 16) & 255)
  result.push(u32_shr(h2, 8) & 255)
  result.push(h2 & 255)
  result.push(u32_shr(h3, 24) & 255)
  result.push(u32_shr(h3, 16) & 255)
  result.push(u32_shr(h3, 8) & 255)
  result.push(h3 & 255)
  result.push(u32_shr(h4, 24) & 255)
  result.push(u32_shr(h4, 16) & 255)
  result.push(u32_shr(h4, 8) & 255)
  result.push(h4 & 255)
  result
}

///|
/// GC content calculation
pub fn gc(seq : String) -> Double {
  let mut gc = 0
  let mut total = 0
  let mut i = 0
  while i < seq.length() {
    match seq.get_char(i) {
      Some(c) => {
        let code = c.to_int()
        let upper_code = if code >= 97 && code <= 122 {
          code - 32
        } else {
          code
        }
        if upper_code == 71 || upper_code == 67 || upper_code == 83 {
          gc = gc + 1
        }
        if upper_code >= 65 && upper_code <= 90 {
          total = total + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  if total == 0 {
    0.0
  } else {
    gc.to_double() / total.to_double() * 100.0
  }
}

///|
/// ProteinAnalysis class
pub struct ProteinAnalysis {
  sequence : String
} derive(Eq)

///|
pub fn ProteinAnalysis::new(sequence : String) -> ProteinAnalysis {
  ProteinAnalysis::{ sequence, }
}

///|
pub fn ProteinAnalysis::count_amino_acids(
  self : ProteinAnalysis,
) -> Map[Char, Int] {
  let counts = Map([], capacity=20)
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) =>
        match counts.get(c) {
          Some(n) => counts[c] = n + 1
          None => counts[c] = 1
        }
      None => ()
    }
    i = i + 1
  }
  counts
}

///|
pub fn ProteinAnalysis::get_amino_acids_percent(
  self : ProteinAnalysis,
) -> Map[Char, Double] {
  let counts = self.count_amino_acids()
  let total = self.sequence.length().to_double()
  let percent = Map([], capacity=counts.length())
  for c, n in counts {
    percent[c] = n.to_double() / total * 100.0
  }
  percent
}

///|
let amino_acid_weights : Map[Char, Double] = Map(
  [
    ('A', 89.0935),
    ('R', 174.2017),
    ('N', 132.1184),
    ('D', 133.1032),
    ('C', 121.1590),
    ('E', 147.1299),
    ('Q', 146.1451),
    ('G', 75.0669),
    ('H', 155.1552),
    ('I', 131.1736),
    ('L', 131.1736),
    ('K', 146.1882),
    ('M', 149.2124),
    ('F', 165.1900),
    ('P', 115.1310),
    ('S', 105.0930),
    ('T', 119.1197),
    ('W', 204.2262),
    ('Y', 181.1894),
    ('V', 117.1469),
  ],
  capacity=20,
)

///|
pub fn ProteinAnalysis::molecular_weight(self : ProteinAnalysis) -> Double {
  let mut weight = 0.0
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) =>
        match amino_acid_weights.get(c) {
          Some(w) => weight = weight + w
          None => ()
        }
      None => ()
    }
    i = i + 1
  }
  weight - 18.0153 * (self.sequence.length() - 1).to_double()
}

///|
pub fn ProteinAnalysis::aromaticity(self : ProteinAnalysis) -> Double {
  let counts = self.count_amino_acids()
  let total = self.sequence.length().to_double()
  let phe = match counts.get('F') {
    Some(n) => n.to_double()
    None => 0.0
  }
  let trp = match counts.get('W') {
    Some(n) => n.to_double()
    None => 0.0
  }
  let tyr = match counts.get('Y') {
    Some(n) => n.to_double()
    None => 0.0
  }
  (phe + trp + tyr) / total
}

///|
let gravy_values : Map[Char, Double] = Map(
  [
    ('A', 1.8),
    ('R', -4.5),
    ('N', -3.5),
    ('D', -3.5),
    ('C', 2.5),
    ('E', -3.5),
    ('Q', -3.5),
    ('G', -0.4),
    ('H', -3.2),
    ('I', 4.5),
    ('L', 3.8),
    ('K', -3.9),
    ('M', 1.9),
    ('F', 2.8),
    ('P', -1.6),
    ('S', -0.8),
    ('T', -0.7),
    ('W', -0.9),
    ('Y', -1.3),
    ('V', 4.2),
  ],
  capacity=20,
)

///|
pub fn ProteinAnalysis::gravy(self : ProteinAnalysis) -> Double {
  let mut score = 0.0
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) =>
        match gravy_values.get(c) {
          Some(v) => score = score + v
          None => ()
        }
      None => ()
    }
    i = i + 1
  }
  score / self.sequence.length().to_double()
}

///|
/// Melting temperature calculations
pub fn tm_wallace(seq : String) -> Double {
  let mut at = 0
  let mut gc = 0
  let mut i = 0
  while i < seq.length() {
    match seq.get_char(i) {
      Some(c) => {
        let code = c.to_int()
        let upper_code = if code >= 97 && code <= 122 {
          code - 32
        } else {
          code
        }
        if upper_code == 65 || upper_code == 84 {
          at = at + 1
        } else if upper_code == 71 || upper_code == 67 {
          gc = gc + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  2.0 * at.to_double() + 4.0 * gc.to_double()
}

///|
pub fn tm_gc(s : String) -> Double {
  let gc_percent = gc(s)
  let n = s.length().to_double()
  81.5 + 0.41 * gc_percent - 600.0 / n + 16.6 * -1.30103
}
