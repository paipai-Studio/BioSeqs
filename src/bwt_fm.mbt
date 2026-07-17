///|
/// Burrows-Wheeler Transform (BWT) and FM-index.
///
/// References:
/// - Burrows M, Wheeler DJ (1994). "A Block-sorting Lossless Data
///   Compression Algorithm." Digital SRC Research Report 124.
/// - Ferragina P, Manzini G (2000). "Opportunistic data structures with
///   applications." Proc. 41st Symposium on Foundations of Computer Science.
/// - Langmead B, Salzberg SL (2012). "Fast gapped-read alignment with
///   Bowtie 2." Nature Methods. 9(4):357-9.

// ---------------------------------------------------------------------------
// Burrows-Wheeler Transform
// ---------------------------------------------------------------------------

///|
/// Sentinel character appended to the input string.
/// Must be lexicographically smaller than any other character.
let bwt_sentinel : Char = '$'

///|
/// Compute the Burrows-Wheeler Transform of a string.
///
/// The input string must not contain the sentinel character '$'.
/// The sentinel is automatically appended.
///
/// Returns the BWT string.
pub fn bwt(s : String) -> String {
  let text = s + "$"
  let n = text.length()
  // Build suffix array (indices of rotation start positions).
  let sa : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    sa.push(i)
  }
  // Sort suffix array by rotation.
  bwt_suffix_sort(sa, text, n)
  // Build BWT: last column of sorted rotations.
  let bwt_buf = StringBuilder::new()
  for i = 0; i < n; i = i + 1 {
    let rot_start = sa[i]
    // Last character of rotation is at (rot_start - 1) mod n.
    let last_idx = if rot_start == 0 { n - 1 } else { rot_start - 1 }
    bwt_buf.write_char(text.unsafe_get(last_idx).unsafe_to_char())
  }
  bwt_buf.to_string()
}

///|
/// Compute the inverse Burrows-Wheeler Transform.
///
/// Given a BWT string, reconstruct the original text (without sentinel).
pub fn bwt_inverse(bwt_str : String) -> String {
  let n = bwt_str.length()
  if n == 0 {
    return ""
  }
  // Build C-table.
  let c_table = bwt_compute_c_table(bwt_str)
  // Find the position of '$' in BWT.
  let mut dollar_pos = 0
  for i = 0; i < n; i = i + 1 {
    if bwt_str.unsafe_get(i).unsafe_to_char() == bwt_sentinel {
      dollar_pos = i
      break
    }
  }
  // Reconstruct text from right to left using LF mapping.
  let buf = StringBuilder::new()
  let mut current = dollar_pos
  for idx = 0; idx < n; idx = idx + 1 {
    let c = bwt_str.unsafe_get(current).unsafe_to_char()
    buf.write_char(c)
    // LF mapping: LF(i) = C[BWT[i]] + rank(BWT[i], i-1)
    // rank(c, i-1) = number of c's in BWT[0..i-1]
    let rank_below = bwt_rank_up_to(bwt_str, c, current - 1)
    let c_count = match c_table.get(c) {
      Some(v) => v
      None => 0
    }
    current = c_count + rank_below
  }
  // The result is reversed and includes '$'.
  let reversed = bwt_reverse_string(buf.to_string())
  // Remove sentinel (should be at the end).
  if reversed.length() > 0 &&
    reversed.unsafe_get(reversed.length() - 1).unsafe_to_char() == bwt_sentinel {
    reversed[0:reversed.length() - 1].to_owned()
  } else {
    reversed
  }
}

///|
/// Sort suffix array indices by their rotation text.
fn bwt_suffix_sort(sa : Array[Int], text : String, n : Int) -> Unit {
  for i = 0; i < n; i = i + 1 {
    let mut min_idx = i
    for j = i + 1; j < n; j = j + 1 {
      if bwt_compare_suffix(sa[j], sa[min_idx], text, n) < 0 {
        min_idx = j
      }
    }
    let tmp = sa[i]
    sa[i] = sa[min_idx]
    sa[min_idx] = tmp
  }
}

///|
/// Compare two suffixes starting at positions a and b.
fn bwt_compare_suffix(a : Int, b : Int, text : String, n : Int) -> Int {
  for i = 0; i < n; i = i + 1 {
    let ca = text.unsafe_get((a + i) % n).unsafe_to_char()
    let cb = text.unsafe_get((b + i) % n).unsafe_to_char()
    if ca < cb {
      return -1
    } else if ca > cb {
      return 1
    }
  }
  0
}

///|
/// Compute rank of character c up to (and including) position pos in s.
/// If pos < 0, returns 0.
fn bwt_rank_up_to(s : String, c : Char, pos : Int) -> Int {
  if pos < 0 {
    return 0
  }
  let mut count = 0
  let limit = if pos >= s.length() { s.length() - 1 } else { pos }
  for i = 0; i <= limit; i = i + 1 {
    if s.unsafe_get(i).unsafe_to_char() == c {
      count = count + 1
    }
  }
  count
}

///|
/// Compute C-table: for each character c, C[c] = number of characters
/// lexicographically smaller than c.
fn bwt_compute_c_table(bwt_str : String) -> Map[Char, Int] {
  // Count occurrences of each character.
  let counts : Map[Char, Int] = Map([], capacity=0)
  for i = 0; i < bwt_str.length(); i = i + 1 {
    let c = bwt_str.unsafe_get(i).unsafe_to_char()
    match counts.get(c) {
      Some(v) => counts[c] = v + 1
      None => counts[c] = 1
    }
  }
  // Collect characters and sort.
  let chars : Array[Char] = Array::new()
  let keys = counts.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    chars.push(keys[i])
  }
  // Bubble sort characters.
  for i = 0; i < chars.length(); i = i + 1 {
    for j = i + 1; j < chars.length(); j = j + 1 {
      if chars[j] < chars[i] {
        let tmp = chars[i]
        chars[i] = chars[j]
        chars[j] = tmp
      }
    }
  }
  // Build C-table.
  let c_table : Map[Char, Int] = Map([], capacity=0)
  let mut cumulative = 0
  for i = 0; i < chars.length(); i = i + 1 {
    c_table[chars[i]] = cumulative
    let cnt = match counts.get(chars[i]) {
      Some(v) => v
      None => 0
    }
    cumulative = cumulative + cnt
  }
  c_table
}

// ---------------------------------------------------------------------------
// FM-index
// ---------------------------------------------------------------------------

///|
/// FM-index for efficient pattern matching on a text.
///
/// Consists of:
/// - BWT string (last column)
/// - C-table (character counts of lexicographically smaller chars)
/// - Occurrence table (full for correctness)
/// - Suffix array (stored explicitly for accurate locate)
pub struct FMIndex {
  bwt_str : String
  c_table : Map[Char, Int]
  occ : Map[Char, Array[Int]]
  sa : Array[Int]
  text_len : Int
}

///|
/// Build an FM-index from a text string.
///
/// The text must not contain '$'. The sentinel is appended automatically.
pub fn FMIndex::new(text : String) -> FMIndex {
  let bwt_str = bwt(text)
  let n = bwt_str.length()
  let c_table = bwt_compute_c_table(bwt_str)
  // Build suffix array.
  let t = text + "$"
  let sa_arr : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    sa_arr.push(i)
  }
  bwt_suffix_sort(sa_arr, t, n)
  // Build full occurrence table:
  // occ[c][i] = number of occurrences of c in bwt[0..i] (inclusive)
  let occ : Map[Char, Array[Int]] = Map([], capacity=0)
  let chars = c_table.keys().collect()
  for ci = 0; ci < chars.length(); ci = ci + 1 {
    let c = chars[ci]
    let arr : Array[Int] = Array::new()
    let mut count = 0
    for i = 0; i < n; i = i + 1 {
      if bwt_str.unsafe_get(i).unsafe_to_char() == c {
        count = count + 1
      }
      arr.push(count)
    }
    occ[c] = arr
  }
  FMIndex::{ bwt_str, c_table, occ, sa: sa_arr, text_len: n }
}

///|
/// Get the BWT string.
pub fn FMIndex::bwt(self : FMIndex) -> String {
  self.bwt_str
}

///|
/// Get the length of the original text (including sentinel).
pub fn FMIndex::len(self : FMIndex) -> Int {
  self.text_len
}

///|
/// Count the number of exact occurrences of `pattern` in the text.
///
/// Uses backward search on the FM-index.
pub fn FMIndex::count(self : FMIndex, pattern : String) -> Int {
  let (l, r) = fm_backward_search_fn(self, pattern)
  if l > r {
    0
  } else {
    r - l + 1
  }
}

///|
/// Find all exact occurrences of `pattern` in the text.
///
/// Returns an array of starting positions (0-based) in the original text.
/// Positions are sorted.
pub fn FMIndex::locate(self : FMIndex, pattern : String) -> Array[Int] {
  let (l, r) = fm_backward_search_fn(self, pattern)
  let result : Array[Int] = Array::new()
  if l > r {
    return result
  }
  for i = l; i <= r; i = i + 1 {
    result.push(self.sa[i])
  }
  // Sort positions.
  for i = 0; i < result.length(); i = i + 1 {
    for j = i + 1; j < result.length(); j = j + 1 {
      if result[j] < result[i] {
        let tmp = result[i]
        result[i] = result[j]
        result[j] = tmp
      }
    }
  }
  result
}

///|
/// LF (Last-to-Front) mapping.
///
/// Given a position `i` in the BWT (last column), returns the
/// corresponding position in the first column.
pub fn FMIndex::lf(self : FMIndex, i : Int) -> Int {
  if i < 0 || i >= self.bwt_str.length() {
    return i
  }
  let c = self.bwt_str.unsafe_get(i).unsafe_to_char()
  let c_val = match self.c_table.get(c) {
    Some(v) => v
    None => 0
  }
  // rank up to and including position i
  let rank = match self.occ.get(c) {
    Some(arr) => if i < arr.length() { arr[i] } else { 0 }
    None => 0
  }
  c_val + rank - 1
}

///|
/// Backward search: find the SA interval [l, r] containing all
/// occurrences of `pattern`.
fn fm_backward_search_fn(fm : FMIndex, pattern : String) -> (Int, Int) {
  let mut l = 0
  let mut r = fm.text_len - 1
  let mut i = pattern.length() - 1
  while i >= 0 && l <= r {
    let c = pattern.unsafe_get(i).unsafe_to_char()
    let c_val = match fm.c_table.get(c) {
      Some(v) => v
      None => 0
    }
    // new_l = C[c] + Occ(c, l-1)
    let occ_l = if l > 0 {
      match fm.occ.get(c) {
        Some(arr) => if l - 1 < arr.length() { arr[l - 1] } else { 0 }
        None => 0
      }
    } else {
      0
    }
    let new_l = c_val + occ_l
    // new_r = C[c] + Occ(c, r) - 1
    let occ_r = match fm.occ.get(c) {
      Some(arr) => if r < arr.length() { arr[r] } else { 0 }
      None => 0
    }
    let new_r = c_val + occ_r - 1
    l = new_l
    r = new_r
    i = i - 1
  }
  (l, r)
}

///|
/// Convert a suffix array index to the original text position.
/// Computed by walking LF mapping until we reach the sentinel.
fn fm_sa_position_fn(fm : FMIndex, sa_idx : Int) -> Int {
  let mut pos = sa_idx
  let mut distance = 0
  while fm.bwt_str.unsafe_get(pos).unsafe_to_char() != bwt_sentinel {
    pos = fm_lf_step_fn(fm, pos)
    distance = distance + 1
  }
  fm.text_len - 1 - distance
}

///|
/// Single LF step (internal helper).
fn fm_lf_step_fn(fm : FMIndex, i : Int) -> Int {
  let c = fm.bwt_str.unsafe_get(i).unsafe_to_char()
  let c_val = match fm.c_table.get(c) {
    Some(v) => v
    None => 0
  }
  let rank = match fm.occ.get(c) {
    Some(arr) => if i < arr.length() { arr[i] } else { 0 }
    None => 0
  }
  c_val + rank - 1
}

// ---------------------------------------------------------------------------
// Suffix array construction (explicit, for reference)
// ---------------------------------------------------------------------------

///|
/// Compute the suffix array of a string.
///
/// The sentinel '$' is appended automatically.
/// Returns an array of indices such that sa[i] is the starting position
/// of the i-th lexicographically smallest suffix.
pub fn suffix_array(text : String) -> Array[Int] {
  let t = text + "$"
  let n = t.length()
  let sa : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    sa.push(i)
  }
  bwt_suffix_sort(sa, t, n)
  sa
}

// ---------------------------------------------------------------------------
// Utility: reverse string (local copy to avoid dependency)
// ---------------------------------------------------------------------------

///|
/// Reverse a string.
fn bwt_reverse_string(s : String) -> String {
  let n = s.length()
  let result = StringBuilder::new()
  let mut i = n - 1
  while i >= 0 {
    result.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  result.to_string()
}
