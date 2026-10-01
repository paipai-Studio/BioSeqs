///|
/// BWA-lite: FM-index based exact and approximate string matching.
///
/// Implements the core of BWA:
///   1. Suffix array construction (naive O(n log n))
///   2. Burrows-Wheeler Transform from suffix array
///   3. FM-index (C array + OCC table)
///   4. Backward search for exact pattern matching
///   5. Inexact matching with backtracking (edit distance ≤ z)
///   6. Locate positions via sampled suffix array
///
/// Algorithm follows Li & Durbin, "Fast and accurate short-read alignment
/// with Burrows-Wheeler transform" (Bioinformatics 2009).

// ============================================================================
// 1. Suffix array construction
// ============================================================================

///|
/// Build suffix array by sorting suffix indices.
/// O(n log n) via comparison-based sort.
/// Returns SA where SA[i] = starting position of the i-th smallest suffix.
pub fn bwa_suffix_array(text : String) -> Array[Int] {
  let n = text.length()
  let sa = Array::makei(n, fn(i) { i })
  sa.sort_by(fn(a, b) {
    // Compare suffixes starting at a and b
    let mut i = a
    let mut j = b
    while i < n && j < n {
      if text[i] < text[j] {
        return -1
      }
      if text[i] > text[j] {
        return 1
      }
      i = i + 1
      j = j + 1
    }
    // Shorter suffix comes first
    if i >= n && j >= n {
      0
    } else if i >= n {
      -1
    } else {
      1
    }
  })
  sa
}

// ============================================================================
// 2. Burrows-Wheeler Transform
// ============================================================================

///|
/// Compute BWT from text and suffix array.
/// BWT[i] = text[SA[i] - 1] if SA[i] > 0, else '$' (sentinel).
pub fn bwa_bwt(text : String, sa : Array[Int]) -> String {
  let n = text.length()
  let buf : StringBuilder = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    if sa[i] == 0 {
      buf.write_char(text[n - 1].unsafe_to_char())
    } else {
      buf.write_char(text[sa[i] - 1].unsafe_to_char())
    }
  }
  buf.to_string()
}

///|
/// Compute BWT directly from text (builds SA internally).
pub fn bwa_bwt_from_text(text : String) -> String {
  let sa = bwa_suffix_array(text)
  bwa_bwt(text, sa)
}

// ============================================================================
// 3. FM-index (C array + OCC table)
// ============================================================================

///|
/// FM-index structure.
pub struct FmIndex {
  bwt : String
  sa : Array[Int]        // suffix array (kept for locate)
  c : Map[Char, Int]     // C[c] = number of chars < c in text
  occ : Map[Char, Array[Int]]  // OCC[c][i] = count of c in BWT[0:i]
  n : Int                // text length
  sample_rate : Int      // SA sampling rate
}

///|
/// Build FM-index from text.
pub fn fm_index_build(text : String) -> FmIndex {
  // Append sentinel to make the FM-index work correctly.
  let ext = text + "$"
  let n = ext.length()
  let sa = bwa_suffix_array(ext)
  let bwt = bwa_bwt(ext, sa)
  // Build C array: C[c] = count of characters lexicographically less than c
  let char_counts : Map[Char, Int] = Map([], capacity=0)
  for i = 0; i < n; i = i + 1 {
    let c = ext[i].unsafe_to_char()
    match char_counts.get(c) {
      Some(v) => char_counts[c] = v + 1
      None => char_counts[c] = 1
    }
  }
  // Get sorted unique characters
  let chars : Array[Char] = Array(capacity=0)
  for k, v in char_counts {
    chars.push(k)
  }
  chars.sort_by(fn(a, b) {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  // C[c] = cumulative count of chars before c
  let c_map : Map[Char, Int] = Map::new()
  let mut cum = 0
  for ch in chars {
    c_map[ch] = cum
    cum = cum + char_counts.get(ch).unwrap_or(0)
  }
  // Build OCC table: OCC[c][i] = number of c in BWT[0:i]
  let occ : Map[Char, Array[Int]] = Map([], capacity=0)
  for ch in chars {
    occ[ch] = Array::make(n + 1, 0)
  }
  for i = 0; i < n; i = i + 1 {
    let ch = bwt[i].unsafe_to_char()
    // Copy previous counts
    for k, v in occ {
      v[i + 1] = v[i]
    }
    match occ.get(ch) {
      Some(arr) => arr[i + 1] = arr[i] + 1
      None => {
        let arr = Array::make(n + 1, 0)
        arr[i + 1] = 1
        occ[ch] = arr
      }
    }
  }
  { bwt, sa, c: c_map, occ, n, sample_rate: 1 }
}

// ============================================================================
// 4. Backward search (exact matching)
// ============================================================================

///|
/// Occ(c, i): number of character c in BWT[0:i).
fn fm_occ(index : FmIndex, c : Char, i : Int) -> Int {
  match index.occ.get(c) {
    Some(arr) => {
      if i < 0 {
        0
      } else if i > index.n {
        arr[index.n]
      } else {
        arr[i]
      }
    }
    None => 0
  }
}

///|
/// LF-mapping: LF(i) = C[BWT[i]] + Occ(BWT[i], i)
fn fm_lf(index : FmIndex, i : Int) -> Int {
  let c = index.bwt[i].unsafe_to_char()
  let c_val = match index.c.get(c) {
    Some(v) => v
    None => 0
  }
  c_val + fm_occ(index, c, i)
}

///|
/// Backward search for exact pattern matching.
/// Returns (left, right) such that SA[left..right] are all occurrences.
/// Returns (-1, -1) if pattern not found.
pub fn fm_backward_search(index : FmIndex, pattern : String) -> (Int, Int) {
  let m = pattern.length()
  if m == 0 {
    return (0, index.n)
  }
  // Initialize: [0, n) covers all suffixes
  let mut lo = 0
  let mut hi = index.n
  // Process pattern from right to left
  for i = m - 1; i >= 0; i = i - 1 {
    let c = pattern[i].unsafe_to_char()
    let c_val = match index.c.get(c) {
      Some(v) => v
      None => return (-1, -1) // character not in text
    }
    lo = c_val + fm_occ(index, c, lo)
    hi = c_val + fm_occ(index, c, hi)
    if lo >= hi {
      return (-1, -1) // no match
    }
  }
  (lo, hi)
}

///|
/// Locate all occurrences of a pattern in the text.
/// Returns array of 0-based positions.
pub fn fm_locate(index : FmIndex, pattern : String) -> Array[Int] {
  let (lo, hi) = fm_backward_search(index, pattern)
  if lo < 0 {
    return []
  }
  let positions : Array[Int] = Array(capacity=0)
  for i = lo; i < hi; i = i + 1 {
    let pos = index.sa[i]
    // Filter out the sentinel position (last position in extended text)
    if pos < index.n - 1 {
      positions.push(pos)
    }
  }
  positions.sort_by(fn(a, b) { a - b })
  positions
}

// ============================================================================
// 5. Inexact matching (backtracking with edit distance)
// ============================================================================

///|
/// A candidate in the backtracking search.
struct BtEntry {
  pos : Int      // position in pattern (from right)
  lo : Int       // current interval [lo, hi)
  hi : Int
  edits : Int    // edits used so far
  cigar : Array[Char] // operations (reversed)
}

///|
/// Inexact backward search with edit distance ≤ max_edits.
/// Returns array of (position, edit_distance, cigar_string).
pub struct BwaHit {
  pos : Int          // 0-based position in text
  edits : Int        // edit distance
  cigar : String     // CIGAR string
}

///|
/// Inexact matching: find all occurrences of pattern with edit distance ≤ z.
pub fn fm_inexact_search(
  index : FmIndex,
  pattern : String,
  max_edits : Int,
) -> Array[BwaHit] {
  let m = pattern.length()
  if m == 0 {
    return []
  }
  let results : Array[BwaHit] = Array(capacity=0)
  let seen : Map[Int, Bool] = Map([], capacity=0)
  // BFS/DFS backtracking
  let stack : Array[BtEntry] = Array(capacity=0)
  stack.push({
    pos: m - 1,
    lo: 0,
    hi: index.n,
    edits: 0,
    cigar: Array(capacity=0),
  })
  while stack.length() > 0 {
    let entry = stack[stack.length() - 1]
    let _ = stack.remove(stack.length() - 1)
    if entry.edits > max_edits {
      continue
    }
    // If we've consumed the entire pattern, record hits
    if entry.pos < 0 {
      for i = entry.lo; i < entry.hi; i = i + 1 {
        let sa_pos = index.sa[i]
        if !seen.get(sa_pos).unwrap_or(false) {
          seen[sa_pos] = true
          let cigar_arr = entry.cigar
          {
          let n = cigar_arr.length()
          for i = 0; i < n / 2; i = i + 1 {
            let tmp = cigar_arr[i]
            cigar_arr[i] = cigar_arr[n - 1 - i]
            cigar_arr[n - 1 - i] = tmp
          }
        }
          results.push({
            pos: sa_pos,
            edits: entry.edits,
            cigar: compress_cigar_bwa(cigar_arr),
          })
        }
      }
      continue
    }
    let c = pattern[entry.pos].unsafe_to_char()
    // 1. Match (or mismatch)
    for ch, occ_arr in index.occ {
      
      let c_val = match index.c.get(ch) {
        Some(v) => v
        None => continue
      }
      let new_lo = c_val + fm_occ(index, ch, entry.lo)
      let new_hi = c_val + fm_occ(index, ch, entry.hi)
      if new_lo >= new_hi {
        continue
      }
      let is_match = ch == c
      let new_edits = if is_match { entry.edits } else { entry.edits + 1 }
      if new_edits > max_edits {
        continue
      }
      let new_cigar : Array[Char] = Array(capacity=0)
      for op in entry.cigar {
        new_cigar.push(op)
      }
      new_cigar.push(if is_match { 'M' } else { 'X' })
      stack.push({
        pos: entry.pos - 1,
        lo: new_lo,
        hi: new_hi,
        edits: new_edits,
        cigar: new_cigar,
      })
    }
    // 2. Deletion (gap in pattern, i.e., character in text not in pattern)
    if entry.edits + 1 <= max_edits {
      let new_cigar : Array[Char] = Array(capacity=0)
      for op in entry.cigar {
        new_cigar.push(op)
      }
      new_cigar.push('D')
      stack.push({
        pos: entry.pos,
        lo: entry.lo,
        hi: entry.hi,
        edits: entry.edits + 1,
        cigar: new_cigar,
      })
    }
    // 3. Insertion (gap in text, i.e., character in pattern not in text)
    if entry.edits + 1 <= max_edits && entry.pos > 0 {
      let new_cigar : Array[Char] = Array(capacity=0)
      for op in entry.cigar {
        new_cigar.push(op)
      }
      new_cigar.push('I')
      stack.push({
        pos: entry.pos - 1,
        lo: entry.lo,
        hi: entry.hi,
        edits: entry.edits + 1,
        cigar: new_cigar,
      })
    }
  }
  // Sort by edit distance ascending
  results.sort_by(fn(a, b) { a.edits - b.edits })
  results
}

///|
/// Compress CIGAR operations (shared with minimap).
fn compress_cigar_bwa(ops : Array[Char]) -> String {
  if ops.length() == 0 {
    return ""
  }
  let buf : StringBuilder = StringBuilder()
  let mut count = 1
  let mut last = ops[0]
  for i = 1; i < ops.length(); i = i + 1 {
    if ops[i] == last {
      count = count + 1
    } else {
      buf.write_string(count.to_string())
      buf.write_char(last)
      count = 1
      last = ops[i]
    }
  }
  buf.write_string(count.to_string())
  buf.write_char(last)
  buf.to_string()
}

// ============================================================================
// 6. High-level BWA interface
// ============================================================================

///|
/// BWA index: wraps FM-index for a reference genome.
pub struct BwaIndex {
  ref_name : String
  ref_seq : String
  fm : FmIndex
}

///|
/// Build BWA index from a reference sequence.
pub fn bwa_index_build(ref_name : String, ref_seq : String) -> BwaIndex {
  let fm = fm_index_build(ref_seq)
  { ref_name, ref_seq, fm }
}

///|
/// BWA-MEM style alignment: find best match for a read.
/// Uses exact matching first, falls back to inexact if needed.
pub struct BwaAlignment {
  ref_name : String
  pos : Int          // 0-based position
  strand : Int       // 0=forward, 1=reverse
  edits : Int
  cigar : String
  score : Int
}

///|
/// Align a single read against the BWA index.
pub fn bwa_align(index : BwaIndex, read : String) -> Array[BwaAlignment] {
  let results : Array[BwaAlignment] = Array(capacity=0)
  // Try exact match first
  let exact = fm_locate(index.fm, read)
  for pos in exact {
    let cigar = read.length().to_string() + "M"
    results.push({
      ref_name: index.ref_name,
      pos,
      strand: 0,
      edits: 0,
      cigar,
      score: read.length() * 2,
    })
  }
  // Try reverse complement
  let rc = reverse_complement(read)
  let exact_rc = fm_locate(index.fm, rc)
  for pos in exact_rc {
    let cigar = read.length().to_string() + "M"
    results.push({
      ref_name: index.ref_name,
      pos,
      strand: 1,
      edits: 0,
      cigar,
      score: read.length() * 2,
    })
  }
  // If no exact match, try inexact with 1-2 edits
  if results.length() == 0 {
    let inexact = fm_inexact_search(index.fm, read, 2)
    for hit in inexact {
      results.push({
        ref_name: index.ref_name,
        pos: hit.pos,
        strand: 0,
        edits: hit.edits,
        cigar: hit.cigar,
        score: read.length() * 2 - hit.edits * 3,
      })
    }
    // Also try reverse complement
    let inexact_rc = fm_inexact_search(index.fm, rc, 2)
    for hit in inexact_rc {
      results.push({
        ref_name: index.ref_name,
        pos: hit.pos,
        strand: 1,
        edits: hit.edits,
        cigar: hit.cigar,
        score: read.length() * 2 - hit.edits * 3,
      })
    }
  }
  // Sort by score descending
  results.sort_by(fn(a, b) {
    if a.score > b.score { -1 } else if a.score < b.score { 1 } else { 0 }
  })
  results
}

///|
/// Convert BWA alignment to SAM line.
pub fn bwa_to_sam(
  aln : BwaAlignment,
  read_name : String,
  read_seq : String,
) -> String {
  let flag = if aln.strand == 1 { 16 } else { 0 }
  let pos = aln.pos + 1 // SAM is 1-based
  let mapq = if aln.edits == 0 { 60 } else { 60 - aln.edits * 10 }
  let mapq = if mapq < 0 { 0 } else { mapq }
  let cigar = if aln.cigar.length() == 0 { "*" } else { aln.cigar }
  let seq = if aln.strand == 1 {
    reverse_complement(read_seq)
  } else {
    read_seq
  }
  let nm_tag = "NM:i:" + aln.edits.to_string()
  read_name + "\t" + flag.to_string() + "\t" + aln.ref_name + "\t" + pos.to_string() + "\t" + mapq.to_string() + "\t" + cigar + "\t*\t0\t0\t" + seq + "\t*\t" + nm_tag
}

///|
/// Align multiple reads and produce SAM output.
pub fn bwa_align_batch(
  index : BwaIndex,
  reads : Array[(String, String)], // (name, seq)
) -> String {
  let buf : StringBuilder = StringBuilder()
  // SAM header
  buf.write_string("@SQ\tSN:" + index.ref_name + "\tLN:" + index.ref_seq.length().to_string() + "\n")
  for (name, seq) in reads {
    let results = bwa_align(index, seq)
    if results.length() == 0 {
      buf.write_string(name + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + seq + "\t*\n")
    } else {
      for r in results {
        buf.write_string(bwa_to_sam(r, name, seq) + "\n")
      }
    }
  }
  buf.to_string()
}
