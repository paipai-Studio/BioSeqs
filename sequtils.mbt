///|
/// Bio.SeqUtils - Utility functions for sequence analysis.

///|
/// CheckSum module

pub fn crc32(seq : String) -> Int {
  let mut crc = 0xFFFFFFFF
  let mut i = 0
  while i < seq.length() {
    match seq.get_char(i) {
      Some(c) => {
        let byte = c.to_int() & 0xFF
        let mut c_val = (crc ^ byte) & 0xFFFFFFFF
        let mut j = 0
        while j < 8 {
          let bit = c_val & 1
          c_val = (c_val >> 1) & 0x7FFFFFFF
          if bit == 1 {
            c_val = c_val ^ 0xEDB88320
          }
          j = j + 1
        }
        crc = c_val
      }
      None => ()
    }
    i = i + 1
  }
  (crc ^ 0xFFFFFFFF) | 0
}

///|
pub fn gcg(seq : String) -> Int {
  let seq_upper = seq.to_upper()
  let mut index = 0
  let mut checksum = 0
  let mut i = 0
  while i < seq_upper.length() {
    match seq_upper.get_char(i) {
      Some(c) => {
        index = index + 1
        checksum = checksum + index * c.to_int()
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
  let mut result = ""
  let mut i = 0
  while i < hash_result.length() {
    let byte = hash_result[i]
    let idx1 = byte >> 2
    match chars.get_char(idx1) {
      Some(ch) => result = result + ch.to_string()
      None => ()
    }
    let rem = byte & 0x03
    i = i + 1
    if i < hash_result.length() {
      let byte2 = hash_result[i]
      let idx2 = (rem << 4) | (byte2 >> 4)
      match chars.get_char(idx2) {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
      let rem2 = byte2 & 0x0F
      i = i + 1
      if i < hash_result.length() {
        let byte3 = hash_result[i]
        let idx3 = (rem2 << 2) | (byte3 >> 6)
        match chars.get_char(idx3) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        let idx4 = byte3 & 0x3F
        match chars.get_char(idx4) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        i = i + 1
      } else {
        let idx3 = rem2 << 2
        match chars.get_char(idx3) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        result = result + "="
      }
    } else {
      let idx2 = rem << 4
      match chars.get_char(idx2) {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
      result = result + "=="
    }
  }
  let mut final_result = ""
  let mut j = 0
  while j < result.length() {
    match result.get_char(j) {
      Some(c) => {
        if c != '=' {
          final_result = final_result + c.to_string()
        }
      }
      None => ()
    }
    j = j + 1
  }
  final_result
}

///|
/// SHA1 implementation

pub fn sha1_string(s : String) -> Array[Int] {
  let bytes = string_to_bytes(s)
  sha1(bytes)
}

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

pub fn bytes_to_base64(bytes : Array[Int]) -> String {
  let chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  let mut result = ""
  let mut i = 0
  while i < bytes.length() {
    let byte = bytes[i]
    let idx1 = byte >> 2
    match chars.get_char(idx1) {
      Some(ch) => result = result + ch.to_string()
      None => ()
    }
    let rem = byte & 0x03
    i = i + 1
    if i < bytes.length() {
      let byte2 = bytes[i]
      let idx2 = (rem << 4) | (byte2 >> 4)
      match chars.get_char(idx2) {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
      let rem2 = byte2 & 0x0F
      i = i + 1
      if i < bytes.length() {
        let byte3 = bytes[i]
        let idx3 = (rem2 << 2) | (byte3 >> 6)
        match chars.get_char(idx3) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        let idx4 = byte3 & 0x3F
        match chars.get_char(idx4) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        i = i + 1
      } else {
        let idx3 = rem2 << 2
        match chars.get_char(idx3) {
          Some(ch) => result = result + ch.to_string()
          None => ()
        }
        result = result + "="
      }
    } else {
      let idx2 = rem << 4
      match chars.get_char(idx2) {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
      result = result + "=="
    }
  }
  let mut final_result = ""
  let mut j = 0
  while j < result.length() {
    match result.get_char(j) {
      Some(c) => {
        if c != '=' {
          final_result = final_result + c.to_string()
        }
      }
      None => ()
    }
    j = j + 1
  }
  final_result
}

pub fn sha1(message : Array[Int]) -> Array[Int] {
  let mut h0 = 0x67452301
  let mut h1 = -271733879
  let mut h2 = -1732584194
  let mut h3 = 0x10325476
  let mut h4 = -1009589776

  let original_len = message.length()
  let padded = Array::new()
  let mut i = 0
  while i < original_len {
    padded.push(message[i])
    i = i + 1
  }
  padded.push(128)

  while padded.length() % 64 != 56 {
    padded.push(0)
  }

  let len_bits = original_len * 8
  i = 7
  while i >= 0 {
    padded.push((len_bits >> (i * 8)) & 255)
    i = i - 1
  }

  let mut block_idx = 0
  while block_idx < padded.length() {
    let w = Array::new()
    let mut j = 0
    while j < 16 {
      let idx = block_idx + j * 4
      let b0 = if idx < padded.length() { padded[idx] } else { 0 }
      let b1 = if idx + 1 < padded.length() { padded[idx + 1] } else { 0 }
      let b2 = if idx + 2 < padded.length() { padded[idx + 2] } else { 0 }
      let b3 = if idx + 3 < padded.length() { padded[idx + 3] } else { 0 }
      let word = ((b0 << 24) & 0xFF000000) | ((b1 << 16) & 0x00FF0000) | ((b2 << 8) & 0x0000FF00) | (b3 & 0x000000FF)
      w.push(word)
      j = j + 1
    }

    j = 16
    while j < 80 {
      let s = (w[j - 3] ^ w[j - 8] ^ w[j - 14] ^ w[j - 16]) & 0xFFFFFFFF
      let rot = ((s << 1) & 0xFFFFFFFF) | ((s >> 31) & 1)
      w.push(rot)
      j = j + 1
    }

    let mut a = h0
    let mut b = h1
    let mut c = h2
    let mut d = h3
    let mut e = h4

    j = 0
    while j < 80 {
      let mut f = 0
      let mut k_val = 0
      if j < 20 {
        f = ((b & c) | (((-1 ^ b) & 0xFFFFFFFF) & d)) & 0xFFFFFFFF
        k_val = 0x5A827999
      } else if j < 40 {
        f = (b ^ c ^ d) & 0xFFFFFFFF
        k_val = 0x6ED9EBA1
      } else if j < 60 {
        f = ((b & c) | (b & d) | (c & d)) & 0xFFFFFFFF
        k_val = -1894007588
      } else {
        f = (b ^ c ^ d) & 0xFFFFFFFF
        k_val = -899497514
      }

      let wj = w[j]
      let rot_a = (((a << 5) & 0xFFFFFFFF) | ((a >> 27) & 0x1F)) & 0xFFFFFFFF
      let temp = (rot_a + f + e + k_val + wj) & 0xFFFFFFFF
      e = d
      d = c
      let rot_b = (((b << 30) & 0xFFFFFFFF) | ((b >> 2) & 0x3FFF)) & 0xFFFFFFFF
      c = rot_b
      b = a
      a = temp
      j = j + 1
    }

    h0 = (h0 + a) & 0xFFFFFFFF
    h1 = (h1 + b) & 0xFFFFFFFF
    h2 = (h2 + c) & 0xFFFFFFFF
    h3 = (h3 + d) & 0xFFFFFFFF
    h4 = (h4 + e) & 0xFFFFFFFF

    block_idx = block_idx + 64
  }

  let result = Array::new()
  result.push((h0 >> 24) & 255)
  result.push((h0 >> 16) & 255)
  result.push((h0 >> 8) & 255)
  result.push(h0 & 255)
  result.push((h1 >> 24) & 255)
  result.push((h1 >> 16) & 255)
  result.push((h1 >> 8) & 255)
  result.push(h1 & 255)
  result.push((h2 >> 24) & 255)
  result.push((h2 >> 16) & 255)
  result.push((h2 >> 8) & 255)
  result.push(h2 & 255)
  result.push((h3 >> 24) & 255)
  result.push((h3 >> 16) & 255)
  result.push((h3 >> 8) & 255)
  result.push(h3 & 255)
  result.push((h4 >> 24) & 255)
  result.push((h4 >> 16) & 255)
  result.push((h4 >> 8) & 255)
  result.push(h4 & 255)
  result
}

///|
/// GC content calculation

pub fn gc(seq : String) -> Double {
  let mut gc = 0
  let mut total = 0
  let seq_upper = seq.to_upper()
  let mut i = 0
  while i < seq_upper.length() {
    match seq_upper.get_char(i) {
      Some(c) => {
        let code = c.to_int()
        if code == 71 || code == 67 || code == 83 {
          gc = gc + 1
        }
        if code >= 65 && code <= 90 {
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

pub fn ProteinAnalysis::new(sequence : String) -> ProteinAnalysis {
  ProteinAnalysis::{ sequence }
}

pub fn ProteinAnalysis::count_amino_acids(self : ProteinAnalysis) -> Map[Char, Int] {
  let counts = Map([], capacity=20)
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) => {
        match counts.get(c) {
          Some(n) => counts[c] = n + 1
          None => counts[c] = 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  counts
}

pub fn ProteinAnalysis::get_amino_acids_percent(self : ProteinAnalysis) -> Map[Char, Double] {
  let counts = self.count_amino_acids()
  let total = self.sequence.length().to_double()
  let percent = Map([], capacity=counts.length())
  for c, n in counts {
    percent[c] = n.to_double() / total * 100.0
  }
  percent
}

let amino_acid_weights : Map[Char, Double] = Map([
  ('A', 89.0935), ('R', 174.2017), ('N', 132.1184), ('D', 133.1032),
  ('C', 121.1590), ('E', 147.1299), ('Q', 146.1451), ('G', 75.0669),
  ('H', 155.1552), ('I', 131.1736), ('L', 131.1736), ('K', 146.1882),
  ('M', 149.2124), ('F', 165.1900), ('P', 115.1310), ('S', 105.0930),
  ('T', 119.1197), ('W', 204.2262), ('Y', 181.1894), ('V', 117.1469),
], capacity=20)

pub fn ProteinAnalysis::molecular_weight(self : ProteinAnalysis) -> Double {
  let mut weight = 0.0
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) => {
        match amino_acid_weights.get(c) {
          Some(w) => weight = weight + w
          None => ()
        }
      }
      None => ()
    }
    i = i + 1
  }
  weight - 18.0153 * (self.sequence.length() - 1).to_double()
}

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

let gravy_values : Map[Char, Double] = Map([
  ('A', 1.8), ('R', -4.5), ('N', -3.5), ('D', -3.5), ('C', 2.5),
  ('E', -3.5), ('Q', -3.5), ('G', -0.4), ('H', -3.2), ('I', 4.5),
  ('L', 3.8), ('K', -3.9), ('M', 1.9), ('F', 2.8), ('P', -1.6),
  ('S', -0.8), ('T', -0.7), ('W', -0.9), ('Y', -1.3), ('V', 4.2),
], capacity=20)

pub fn ProteinAnalysis::gravy(self : ProteinAnalysis) -> Double {
  let mut score = 0.0
  let mut i = 0
  while i < self.sequence.length() {
    match self.sequence.get_char(i) {
      Some(c) => {
        match gravy_values.get(c) {
          Some(v) => score = score + v
          None => ()
        }
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
  let seq_upper = seq.to_upper()
  let mut at = 0
  let mut gc = 0
  let mut i = 0
  while i < seq_upper.length() {
    match seq_upper.get_char(i) {
      Some(c) => {
        let code = c.to_int()
        if code == 65 || code == 84 {
          at = at + 1
        } else if code == 71 || code == 67 {
          gc = gc + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  2.0 * at.to_double() + 4.0 * gc.to_double()
}

pub fn tm_gc(s : String) -> Double {
  let gc_percent = gc(s)
  let n = s.length().to_double()
  81.5 + 0.41 * gc_percent - 600.0 / n + 16.6 * (-1.30103)
}