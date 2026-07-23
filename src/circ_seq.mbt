///|
/// Bio.SeqUtils.CircSeq - Circular DNA sequence operations.
/// Inspired by Biopython's Bio.SeqUtils.CircSeq module.
///
/// Provides tools for working with circular DNA sequences (plasmids,
/// bacterial chromosomes, etc.) where coordinates wrap around the
/// sequence ends. Supports operations like restriction enzyme analysis,
/// PCR primer design, and coordinate transformation on circular sequences.

///|
/// Represents a circular DNA sequence with its properties.
pub struct CircSeq {
  /// The sequence string (linear representation)
  sequence : String
  /// Total length of the circular sequence
  length : Int
  /// GC content (0.0-1.0)
  gc_content : Double
  /// Molecule name (e.g., plasmid name)
  name : String
}

///|
/// Construct a circular sequence.
pub fn CircSeq::new(sequence~ : String) -> CircSeq {
  let length = sequence.length()
  let gc = compute_gc_content(sequence)
  CircSeq::{ sequence, length, gc_content: gc, name: "" }
}

///|
/// Construct with a name.
pub fn CircSeq::with_name(sequence~ : String, name~ : String) -> CircSeq {
  let length = sequence.length()
  let gc = compute_gc_content(sequence)
  CircSeq::{ sequence, length, gc_content: gc, name }
}

fn compute_gc_content(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let mut gc = 0
  let mut i = 0
  while i < n {
    let c = char_to_upper(seq.unsafe_get(i))
    if c == 71 || c == 67 {
      gc = gc + 1
    }
    i = i + 1
  }
  gc.to_double() / n.to_double()
}

///|
/// Get the sequence length.
pub fn CircSeq::circ_length(self : CircSeq) -> Int {
  self.length
}

///|
/// Get the GC content.
pub fn CircSeq::circ_gc_content(self : CircSeq) -> Double {
  self.gc_content
}

///|
/// Get the molecule name.
pub fn CircSeq::circ_name(self : CircSeq) -> String {
  self.name
}

///|
/// Get the raw sequence string.
pub fn CircSeq::circ_sequence(self : CircSeq) -> String {
  self.sequence
}

///|
/// Access a base at a circular position (0-based, wraps around).
pub fn CircSeq::circ_at(self : CircSeq, pos : Int) -> UInt16 {
  let wrapped = modulo_pos(pos, self.length)
  self.sequence.unsafe_get(wrapped)
}

///|
/// Get a sub-sequence from start to end (0-based, inclusive on both ends).
/// Handles circular wrapping.
pub fn CircSeq::circ_subseq(self : CircSeq, start : Int, end : Int) -> String {
  let wrapped_start = modulo_pos(start, self.length)
  let wrapped_end = modulo_pos(end, self.length)

  if wrapped_start <= wrapped_end {
    // Simple linear case
    if wrapped_end + 1 <= self.length {
      return self.sequence.substring(start=wrapped_start, end=wrapped_end + 1)
    }
    // Edge case: wrapped_end == length (shouldn't happen with modulo)
    return self.sequence.substring(start=wrapped_start)
  } else {
    // Wraps around: start > end in linear space
    // e.g., length=10, start=8, end=2 -> get 8,9,0,1,2
    let first_part = self.sequence.substring(start=wrapped_start)
    let second_part = self.sequence.substring(start=0, end=wrapped_end + 1)
    return first_part + second_part
  }
}

///|
/// Circular slice: returns a linear segment starting at pos with given length.
/// Position 0 is the start of the sequence.
pub fn CircSeq::circ_slice(self : CircSeq, pos : Int, len : Int) -> String {
  let mut result = ""
  let mut i = 0
  while i < len {
    result = result + u16_to_str(self.circ_at(pos + i))
    i = i + 1
  }
  result
}

///|
/// Rotate the sequence so that the given position becomes the new origin.
/// Returns a new CircSeq with the rotated sequence.
pub fn CircSeq::circ_rotate(self : CircSeq, new_origin : Int) -> CircSeq {
  let wrapped = modulo_pos(new_origin, self.length)
  if wrapped == 0 {
    return self
  }
  // New sequence: seq[wrapped:] + seq[:wrapped]
  let new_seq = self.sequence.substring(start=wrapped) + self.sequence.substring(start=0, end=wrapped)
  let gc = compute_gc_content(new_seq)
  CircSeq::{ sequence: new_seq, length: self.length, gc_content: gc, name: self.name }
}

///|
/// Reverse complement of the circular sequence.
pub fn CircSeq::circ_reverse_complement(self : CircSeq) -> CircSeq {
  let rc = circ_rev_comp(self.sequence)
  let gc = compute_gc_content(rc)
  CircSeq::{ sequence: rc, length: self.length, gc_content: gc, name: self.name }
}

fn circ_rev_comp(seq : String) -> String {
  let mut result = ""
  let mut i = seq.length() - 1
  while i >= 0 {
    let c = char_to_upper(seq.unsafe_get(i))
    let comp : UInt16 = match c {
      65 => 84   // 'A' => 'T'
      84 => 65   // 'T' => 'A'
      71 => 67   // 'G' => 'C'
      67 => 71   // 'C' => 'G'
      85 => 65   // 'U' => 'A'
      78 => 78   // 'N' => 'N'
      82 => 89   // 'R' => 'Y'
      89 => 82   // 'Y' => 'R'
      83 => 83   // 'S' => 'S'
      87 => 87   // 'W' => 'W'
      75 => 77   // 'K' => 'M'
      77 => 75   // 'M' => 'K'
      66 => 86   // 'B' => 'V'
      86 => 66   // 'V' => 'B'
      68 => 72   // 'D' => 'H'
      72 => 68   // 'H' => 'D'
      _ => 78    // 'N' default
    }
    result = result + u16_to_str(comp)
    i = i - 1
  }
  result
}

///|
/// Find all occurrences of a pattern in the circular sequence.
/// Returns an array of start positions (0-based).
pub fn CircSeq::circ_find_all(self : CircSeq, pattern : String) -> Array[Int] {
  let results : Array[Int] = Array::new()
  let plen = pattern.length()
  if plen == 0 || plen > self.length {
    return results
  }

  // Search the doubled sequence to handle circular wrapping
  let doubled = self.sequence + self.sequence
  let search_len = self.length  // Only search up to the original length to avoid duplicates

  let mut i = 0
  while i <= search_len - plen {
    // Check if pattern matches at position i
    let mut match_found = true
    let mut j = 0
    while j < plen {
      let pos = i + j
      let c = doubled.unsafe_get(pos)
      let p = pattern.unsafe_get(j)
      if c != p {
        match_found = false
        break
      }
      j = j + 1
    }
    if match_found {
      results.push(modulo_pos(i, self.length))
    }
    i = i + 1
  }
  results
}

///|
/// Find the first occurrence of a pattern in the circular sequence.
/// Returns the 0-based position or -1 if not found.
pub fn CircSeq::circ_find(self : CircSeq, pattern : String) -> Int {
  let positions = self.circ_find_all(pattern)
  if positions.length() > 0 {
    return positions[0]
  }
  -1
}

///|
/// Find the reverse complement of a pattern in the circular sequence.
pub fn CircSeq::circ_find_rc(self : CircSeq, pattern : String) -> Array[Int] {
  let rc_pattern = circ_rev_comp(pattern)
  self.circ_find_all(rc_pattern)
}

///|
/// Check if a restriction enzyme site exists in the sequence.
pub fn CircSeq::circ_has_site(self : CircSeq, site : String) -> Bool {
  self.circ_find(site) >= 0
}

///|
/// Find all restriction enzyme cut sites on both strands.
/// Returns an array of (position, strand) tuples.
pub fn CircSeq::circ_find_sites(self : CircSeq, site : String) -> Array[(Int, String)] {
  let results : Array[(Int, String)] = Array::new()
  let fwd = self.circ_find_all(site)
  for pos in fwd {
    results.push((pos, "+"))
  }
  let rc_sites = self.circ_find_rc(site)
  for pos in rc_sites {
    results.push((pos, "-"))
  }
  results
}

///|
/// Simulate a restriction digest on the circular sequence.
/// Returns an array of fragment sizes (in bp).
pub fn CircSeq::circ_digest(self : CircSeq, sites : Array[String]) -> Array[Int] {
  let cut_positions : Array[Int] = Array::new()

  // Collect all cut positions from both strands for all enzymes
  for site in sites {
    let fwd_positions = self.circ_find_all(site)
    for pos in fwd_positions {
      cut_positions.push(pos)
    }
    let rc_positions = self.circ_find_rc(site)
    for pos in rc_positions {
      cut_positions.push(pos)
    }
  }

  // Sort cut positions
  let sorted = sort_positions(cut_positions)

  // Compute fragment sizes
  let fragments : Array[Int] = Array::new()
  if sorted.length() == 0 {
    // No cuts - the whole molecule is one fragment
    fragments.push(self.length)
    return fragments
  }

  // Remove duplicates
  let unique : Array[Int] = Array::new()
  for pos in sorted {
    if unique.length() == 0 || pos != unique[unique.length() - 1] {
      unique.push(pos)
    }
  }

  // Compute fragments
  let n_cuts = unique.length()
  let mut i = 0
  while i < n_cuts {
    let start = unique[i]
    let end = if i == n_cuts - 1 {
      unique[0] + self.length
    } else {
      unique[i + 1]
    }
    let frag_size = end - start
    if frag_size > 0 {
      fragments.push(frag_size)
    }
    i = i + 1
  }

  fragments
}

///|
/// Calculate fragment sizes from a circular digest using specific enzymes.
pub fn CircSeq::circ_digest_enzymes(self : CircSeq, enzymes : Array[String]) -> Array[Int] {
  self.circ_digest(enzymes)
}

///|
/// Find the best primer pair for PCR amplification on the circular sequence.
/// Returns (forward_primer_start, reverse_primer_start, product_length).
/// The forward primer binds to the sense strand, reverse primer to the antisense strand.
pub fn CircSeq::circ_design_primers(
  self : CircSeq,
  target_start : Int,
  target_end : Int,
  primer_len? : Int = 20,
  gc_target? : Double = 0.5,
) -> (Int, Int, Int) {
  let product_len = circular_distance(target_start, target_end, self.length)
  let fwd_start = target_start
  let fwd_seq = self.circ_slice(fwd_start, primer_len)

  // Reverse primer is downstream of target_end, on the opposite strand
  let rev_start = modulo_pos(target_end - primer_len, self.length)
  let rev_seq = self.circ_slice(rev_start, primer_len)

  // Calculate GC content of primers
  let fwd_gc = calc_gc(fwd_seq)
  let rev_gc = calc_gc(rev_seq)

  // Adjust primer length to reach target GC content
  let adj_fwd_len = adjust_primer_length(self, fwd_start, primer_len, gc_target)
  let adj_rev_len = adjust_primer_length(self, rev_start, primer_len, gc_target)

  // Return the positions (primer_start position, not the primer sequence)
  (fwd_start, rev_start, product_len)
}

fn calc_gc(seq : String) -> Double {
  if seq.length() == 0 { return 0.0 }
  let mut gc = 0
  let mut i = 0
  while i < seq.length() {
    let c = char_to_upper(seq.unsafe_get(i))
    if c == 71 || c == 67 { gc = gc + 1 }
    i = i + 1
  }
  gc.to_double() / seq.length().to_double()
}

fn adjust_primer_length(circ : CircSeq, start : Int, target_len : Int, gc_target : Double) -> Int {
  let mut len = target_len
  let max_len = if circ.length < 40 { circ.length } else { 40 }
  while len < max_len {
    let seq = circ.circ_slice(start, len)
    let gc = calc_gc(seq)
    if gc >= gc_target {
      return len
    }
    len = len + 1
  }
  target_len
}

///|
/// Insert a sequence at a specific position in the circular sequence.
pub fn CircSeq::circ_insert(self : CircSeq, pos : Int, insert_seq : String) -> CircSeq {
  let wrapped = modulo_pos(pos, self.length)
  let new_seq = self.sequence.substring(start=0, end=wrapped) + insert_seq + self.sequence.substring(start=wrapped)
  let new_length = new_seq.length()
  let gc = compute_gc_content(new_seq)
  CircSeq::{ sequence: new_seq, length: new_length, gc_content: gc, name: self.name }
}

///|
/// Delete a range from the circular sequence [start, end] inclusive.
pub fn CircSeq::circ_delete(self : CircSeq, start : Int, end : Int) -> CircSeq {
  let new_seq = circular_delete(self.sequence, start, end)
  let new_length = new_seq.length()
  let gc = compute_gc_content(new_seq)
  CircSeq::{ sequence: new_seq, length: new_length, gc_content: gc, name: self.name }
}

fn circular_delete(seq : String, start : Int, end : Int) -> String {
  let wrapped_start = modulo_pos(start, seq.length())
  let wrapped_end = modulo_pos(end, seq.length())

  if wrapped_start <= wrapped_end {
    // Linear delete
    return seq.substring(start=0, end=wrapped_start) + seq.substring(start=wrapped_end + 1)
  } else {
    // Wraps around: delete from wrapped_start to end, and from start to wrapped_end
    return seq.substring(start=wrapped_end + 1, end=wrapped_start)
  }
}

///|
/// Compute the circular distance between two positions.
/// Returns the shortest distance going from start to end in the forward direction.
pub fn circular_distance(start : Int, end : Int, length : Int) -> Int {
  let wrapped_start = modulo_pos(start, length)
  let wrapped_end = modulo_pos(end, length)
  if wrapped_end >= wrapped_start {
    wrapped_end - wrapped_start
  } else {
    length - wrapped_start + wrapped_end
  }
}

///|
/// Map circular coordinates to linear coordinates (no wrapping).
pub fn circ_to_linear(circ_pos : Int, length : Int) -> Int {
  modulo_pos(circ_pos, length)
}

///|
/// Map linear coordinates back to circular coordinates.
pub fn linear_to_circular(linear_pos : Int, length : Int) -> Int {
  modulo_pos(linear_pos, length)
}

///|
/// Check if a position is within a range on the circular sequence.
pub fn is_within_circular_range(pos : Int, range_start : Int, range_end : Int, length : Int) -> Bool {
  let p = modulo_pos(pos, length)
  let rs = modulo_pos(range_start, length)
  let re = modulo_pos(range_end, length)

  if rs <= re {
    // Range doesn't wrap
    return p >= rs && p <= re
  } else {
    // Range wraps around
    return p >= rs || p <= re
  }
}

///|
/// Compute the overlap between two ranges on the circular sequence.
/// Returns the number of overlapping bases.
pub fn circular_overlap(start1 : Int, end1 : Int, start2 : Int, end2 : Int, length : Int) -> Int {
  let s1 = modulo_pos(start1, length)
  let e1 = modulo_pos(end1, length)
  let s2 = modulo_pos(start2, length)
  let e2 = modulo_pos(end2, length)

  // Convert ranges to arrays of positions and find intersection
  let range1 = expand_range(s1, e1, length)
  let range2 = expand_range(s2, e2, length)

  // Count overlap
  let mut count = 0
  for pos in range1 {
    if range2.contains(pos) {
      count = count + 1
    }
  }
  count
}

fn expand_range(start : Int, end : Int, length : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  if start <= end {
    let mut i = start
    while i <= end {
      result.push(modulo_pos(i, length))
      i = i + 1
    }
  } else {
    // Wraps around
    let mut i = start
    while i < length {
      result.push(i)
      i = i + 1
    }
    let mut j = 0
    while j <= end {
      result.push(j)
      j = j + 1
    }
  }
  result
}

// --- Utility Functions ---

fn modulo_pos(pos : Int, length : Int) -> Int {
  if length <= 0 {
    return 0
  }
  let result = pos % length
  if result < 0 {
    result + length
  } else {
    result
  }
}

fn sort_positions(positions : Array[Int]) -> Array[Int] {
  let result = positions.copy()
  let n = result.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      if result[j] < result[i] {
        let temp = result[i]
        result[i] = result[j]
        result[j] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Generate a random circular sequence of given length and GC content.
pub fn circ_random_sequence(length : Int, gc_content : Double, seed? : Int = 42) -> CircSeq {
  let bases = ["A", "T", "G", "C"]
  // Adjust probabilities based on GC content
  let at_prob = (1.0 - gc_content) / 2.0
  let gc_prob = gc_content / 2.0

  let probs = [at_prob, at_prob, gc_prob, gc_prob]

  // Simple LCG pseudo-random number generator (deterministic)
  let mut state = seed
  let mut seq = ""
  let mut i = 0
  while i < length {
    state = (state * 1103515245 + 12345) % 2147483647
    let rand_val = state.to_double() / 2147483647.0
    let base_idx = pick_base(rand_val, probs)
    seq = seq + bases[base_idx]
    i = i + 1
  }
  CircSeq::new(sequence=seq)
}

fn pick_base(rand_val : Double, probs : Array[Double]) -> Int {
  let cumulative = [probs[0], probs[0] + probs[1], probs[0] + probs[1] + probs[2], 1.0]
  if rand_val <= cumulative[0] {
    return 0
  } else if rand_val <= cumulative[1] {
    return 1
  } else if rand_val <= cumulative[2] {
    return 2
  } else {
    return 3
  }
}

///|
/// Simulate a ligation reaction on two circular sequences.
/// Returns the concatenated circular product.
pub fn circ_ligate(linear_fragment : String, vector : CircSeq, insert_pos : Int) -> CircSeq {
  vector.circ_insert(insert_pos, linear_fragment)
}

///|
/// Compute Tm (melting temperature) for a primer on a circular template.
/// Uses the nearest-neighbor method approximation.
pub fn circ_tm(template : CircSeq, start : Int, len : Int, oligo_conc? : Double = 5.0e-8) -> Double {
  let seq = template.circ_slice(start, len)
  let gc = calc_gc(seq)
  let at = 1.0 - gc

  // Nearest-neighbor Tm calculation (simplified)
  // Tm = 64.9 + 41 * (GC - 0.16) for sequences with < 14 bases
  // Tm = 64.9 + 41 * (GC - 0.16) for longer sequences
  // Actually, use the standard formula:
  // Tm = (delta_H / delta_S) - 273.15 + 16.6 * ln([Na+])
  // Simplified: Tm = 4 * GC_length + 2 * AT_length (rough estimate)
  let gc_count = (gc * len.to_double()).round().to_int()
  let at_count = len - gc_count

  // Rough Tm estimate
  let tm = 4.0 * gc_count.to_double() + 2.0 * at_count.to_double()
  tm
}

///|
/// List common restriction enzyme recognition sequences.
pub fn circ_common_restriction_sites() -> Array[(String, String)] {
  [
    ("EcoRI", "GAATTC"),
    ("BamHI", "GGATCC"),
    ("HindIII", "AAGCTT"),
    ("XbaI", "TCTAGA"),
    ("SalI", "GTCGAC"),
    ("PstI", "CTGCAG"),
    ("KpnI", "GGTACC"),
    ("NotI", "GCGGCCGC"),
    ("NcoI", "CCATGG"),
    ("XhoI", "CTCGAG"),
    ("SpeI", "ACTAGT"),
    ("BglII", "AGATCT"),
    ("MluI", "ACGCGT"),
    ("NheI", "GCTAGC"),
    ("PacI", "TTAATTAA"),
    ("PvuII", "CAGCTG"),
    ("ScaI", "AGTACT"),
    ("SmaI", "CCCGGG"),
    ("StuI", "AGGCCT"),
    ("XmaI", "CCCGGG")
  ]
}

///|
/// Find all unique restriction fragments from a digest.
pub fn circ_unique_fragments(digest_fragments : Array[Int]) -> Array[Int] {
  let sorted = sort_positions(digest_fragments)
  let unique : Array[Int] = Array::new()
  for frag in sorted {
    if unique.length() == 0 || frag != unique[unique.length() - 1] {
      unique.push(frag)
    }
  }
  unique
}

///|
/// Find the largest unique fragment in a digest (useful for linearization).
pub fn circ_largest_fragment(digest_fragments : Array[Int]) -> Int {
  let mut largest = 0
  for frag in digest_fragments {
    if frag > largest {
      largest = frag
    }
  }
  largest
}
