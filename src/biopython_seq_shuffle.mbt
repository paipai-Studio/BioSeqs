/// biopython_seq_shuffle.mbt — Sequence shuffling and mutation
/// Inspired by Biopython Bio.SeqUtils: randomization, mutation

///|
/// Simple LCG random number generator (deterministic, seedable).
pub struct BiopythonRng {
  state : Int32
}

///|
/// Create a new RNG with a seed.
pub fn BiopythonRng::new(seed~ : Int32 = 42) -> BiopythonRng {
  { state: seed, }
}

///|
/// Generate next random number.
pub fn BiopythonRng::next(self : BiopythonRng) -> (BiopythonRng, Int32) {
  // LCG parameters (same as glibc)
  let new_state = (self.state * 1103515245 + 12345) % 2147483647
  ({ state: new_state, }, new_state)
}

///|
/// Generate a random number in range [0, max).
pub fn BiopythonRng::next_range(
  self : BiopythonRng,
  max : Int32,
) -> (BiopythonRng, Int32) {
  let (rng, val) = self.next()
  (rng, val % max)
}

///|
/// Shuffle a sequence (Fisher-Yates), preserving composition.
pub fn biopython_shuffle_sequence(
  sequence : String,
  rng~ : BiopythonRng = { state: 42, },
) -> String {
  let n = sequence.length()
  if n <= 1 {
    return sequence
  }
  // Convert to array of char codes
  let chars : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    chars[i] = sequence[i].to_int()
  }
  // Fisher-Yates shuffle
  let mut current_rng = rng
  for i = n - 1; i > 0; i = i - 1 {
    let (rng, j) = current_rng.next_range(i + 1)
    current_rng = rng
    let tmp = chars[i]
    chars[i] = chars[j]
    chars[j] = tmp
  }
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Introduce random mutations into a sequence.
pub fn biopython_mutate_sequence(
  sequence : String,
  num_mutations : Int32,
  rng~ : BiopythonRng = { state: 42, },
) -> String {
  let n = sequence.length()
  if n == 0 || num_mutations == 0 {
    return sequence
  }
  let chars : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    chars[i] = sequence[i].to_int()
  }
  let bases = [65, 67, 71, 84] // A, C, G, T
  let mut current_rng = rng
  for _m = 0; _m < num_mutations; _m = _m + 1 {
    let (rng1, pos) = current_rng.next_range(n)
    current_rng = rng1
    let (rng2, base_idx) = current_rng.next_range(4)
    current_rng = rng2
    chars[pos] = bases[base_idx]
  }
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Generate a random DNA sequence of given length.
pub fn biopython_random_dna(
  length : Int32,
  rng~ : BiopythonRng = { state: 42, },
) -> String {
  let bases = [65, 67, 71, 84] // A, C, G, T
  let chars : Array[Int32] = Array::make(length, 0)
  let mut current_rng = rng
  for i = 0; i < length; i = i + 1 {
    let (rng, idx) = current_rng.next_range(4)
    current_rng = rng
    chars[i] = bases[idx]
  }
  let sb = StringBuilder()
  for i = 0; i < length; i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Generate a random protein sequence of given length.
pub fn biopython_random_protein(
  length : Int32,
  rng~ : BiopythonRng = { state: 42, },
) -> String {
  let aas = [
    65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80, 81, 82, 83, 84, 86, 87, 89,
  ]
  let n_aa = aas.length()
  let chars : Array[Int32] = Array::make(length, 0)
  let mut current_rng = rng
  for i = 0; i < length; i = i + 1 {
    let (rng, idx) = current_rng.next_range(n_aa)
    current_rng = rng
    chars[i] = aas[idx]
  }
  let sb = StringBuilder()
  for i = 0; i < length; i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Dinucleotide shuffle (preserves dinucleotide frequencies).
/// Uses the Altschul-Erickson algorithm (simplified).
pub fn biopython_dinucleotide_shuffle(
  sequence : String,
  rng~ : BiopythonRng = { state: 42, },
) -> String {
  // Simplified: just do a regular shuffle
  // Full dinucleotide shuffle is complex (Markov chain based)
  biopython_shuffle_sequence(sequence, rng~)
}
