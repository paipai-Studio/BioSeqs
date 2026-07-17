///|
/// Bio.SeqUtils.Complexity - Sequence complexity and composition analysis.
/// Inspired by Biopython's Bio.SeqUtils and dust/masking algorithms.

///|
/// Calculate the Shannon entropy of a sequence (in nats, using natural log).
/// Higher entropy = more complex/random sequence.
/// This is a duplicate of biostrings.mbt's shannon_entropy but uses natural log.
/// Use shannon_entropy_bits for log2-based entropy.
pub fn shannon_entropy_bits(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  // Count character frequencies
  let freq : Map[Char, Int] = Map([], capacity=16)
  let mut i = 0
  while i < n {
    let c = seq.unsafe_get(i).unsafe_to_char()
    match freq.get(c) {
      Some(count) => freq[c] = count + 1
      None => freq[c] = 1
    }
    i = i + 1
  }
  // Calculate entropy in bits (using log2)
  let mut entropy = 0.0
  let nf = n.to_double()
  for kv in freq {
    let p = kv.1.to_double() / nf
    entropy = entropy - p * @math.log2(p)
  }
  entropy
}

///|
/// Calculate the linguistic complexity of a sequence.
/// Linguistic complexity = observed_kmers / max_possible_kmers
/// For a sequence of length n with alphabet size a,
/// max possible k-mers of length k = min(a^k, n-k+1)
pub fn linguistic_complexity(seq : String, k? : Int = 2) -> Double {
  let n = seq.length()
  if n < k {
    return 0.0
  }
  let alphabet_size = count_unique_chars(seq)
  let max_kmers = if alphabet_size > 0 {
    let max_theoretical = @math.pow(alphabet_size.to_double(), k.to_double()).to_int()
    let max_practical = n - k + 1
    if max_theoretical < max_practical {
      max_theoretical
    } else {
      max_practical
    }
  } else {
    0
  }
  if max_kmers == 0 {
    return 0.0
  }
  let observed = count_unique_kmers(seq, k)
  observed.to_double() / max_kmers.to_double()
}

///|
/// Count the number of unique characters in a sequence.
fn count_unique_chars(seq : String) -> Int {
  let seen : Map[Char, Bool] = Map([], capacity=16)
  let n = seq.length()
  let mut i = 0
  while i < n {
    let c = seq.unsafe_get(i).unsafe_to_char()
    seen[c] = true
    i = i + 1
  }
  seen.size()
}

///|
/// Count the number of unique k-mers in a sequence.
fn count_unique_kmers(seq : String, k : Int) -> Int {
  let n = seq.length()
  if n < k {
    return 0
  }
  let seen : Map[String, Bool] = Map([], capacity=64)
  let mut i = 0
  while i <= n - k {
    let kmer = seq[i:i + k].to_owned()
    seen[kmer] = true
    i = i + 1
  }
  seen.size()
}

///|
/// Calculate the GC content percentage of a sequence (alias).
pub fn gc_content_percent(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let mut gc = 0
  let mut i = 0
  while i < n {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if c == 'G' || c == 'g' || c == 'C' || c == 'c' {
      gc = gc + 1
    }
    i = i + 1
  }
  gc.to_double() / n.to_double() * 100.0
}

///|
/// Calculate the AT content of a sequence.
pub fn at_content(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let mut at = 0
  let mut i = 0
  while i < n {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if c == 'A' || c == 'a' || c == 'T' || c == 't' {
      at = at + 1
    }
    i = i + 1
  }
  at.to_double() / n.to_double() * 100.0
}

///|
/// Calculate the GC skew of a sequence.
/// GC skew = (G - C) / (G + C)
pub fn gc_skew(seq : String) -> Double {
  let mut g = 0.0
  let mut c = 0.0
  let n = seq.length()
  let mut i = 0
  while i < n {
    let ch = seq.unsafe_get(i).unsafe_to_char()
    if ch == 'G' || ch == 'g' {
      g = g + 1.0
    } else if ch == 'C' || ch == 'c' {
      c = c + 1.0
    }
    i = i + 1
  }
  if g + c == 0.0 {
    return 0.0
  }
  (g - c) / (g + c)
}

///|
/// Calculate the AT skew of a sequence.
/// AT skew = (A - T) / (A + T)
pub fn at_skew(seq : String) -> Double {
  let mut a = 0.0
  let mut t = 0.0
  let n = seq.length()
  let mut i = 0
  while i < n {
    let ch = seq.unsafe_get(i).unsafe_to_char()
    if ch == 'A' || ch == 'a' {
      a = a + 1.0
    } else if ch == 'T' || ch == 't' {
      t = t + 1.0
    }
    i = i + 1
  }
  if a + t == 0.0 {
    return 0.0
  }
  (a - t) / (a + t)
}

///|
/// Calculate the DUST score for low-complexity detection.
/// DUST = sum(kmer_counts^2) / (length * (length-2))
/// Higher scores indicate lower complexity.
pub fn dust_score(seq : String) -> Double {
  let n = seq.length()
  if n < 3 {
    return 0.0
  }
  // Count 3-mers
  let freq : Map[String, Int] = Map([], capacity=64)
  let mut i = 0
  while i <= n - 3 {
    let kmer = seq[i:i + 3].to_owned()
    match freq.get(kmer) {
      Some(count) => freq[kmer] = count + 1
      None => freq[kmer] = 1
    }
    i = i + 1
  }
  // Calculate score
  let mut score = 0.0
  for kv in freq {
    score = score + kv.1.to_double() * kv.1.to_double()
  }
  let denom = (n - 2).to_double() * 3.0
  if denom == 0.0 {
    0.0
  } else {
    score / denom
  }
}

///|
/// Calculate the sequence signature as a k-mer frequency vector.
/// Returns a map of k-mer -> frequency (normalized).
pub fn sequence_signature(seq : String, k? : Int = 3) -> Map[String, Double] {
  let result : Map[String, Double] = Map([], capacity=64)
  let n = seq.length()
  if n < k {
    return result
  }
  let total = n - k + 1
  let mut i = 0
  while i <= n - k {
    let kmer = seq[i:i + k].to_owned()
    match result.get(kmer) {
      Some(count) => result[kmer] = count + 1.0
      None => result[kmer] = 1.0
    }
    i = i + 1
  }
  // Normalize
  let tf = total.to_double()
  for kv in result {
    result[kv.0] = kv.1 / tf
  }
  result
}

///|
/// Calculate the composition vector (nucleotide frequency).
/// Returns (A, T, G, C) frequencies as a tuple.
pub fn nucleotide_frequencies(seq : String) -> (Double, Double, Double, Double) {
  let n = seq.length()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  let mut a = 0
  let mut t = 0
  let mut g = 0
  let mut c = 0
  let mut i = 0
  while i < n {
    let ch = seq.unsafe_get(i).unsafe_to_char()
    match ch {
      'A' | 'a' => a = a + 1
      'T' | 't' => t = t + 1
      'G' | 'g' => g = g + 1
      'C' | 'c' => c = c + 1
      _ => ()
    }
    i = i + 1
  }
  let nf = n.to_double()
  (
    a.to_double() / nf,
    t.to_double() / nf,
    g.to_double() / nf,
    c.to_double() / nf,
  )
}

///|
/// Calculate the CGR (Chaos Game Representation) coordinates for a sequence.
/// Maps a DNA sequence to a 2D point based on iterative function system.
pub fn chaos_game_representation(seq : String) -> (Double, Double) {
  // Corners: A=(0,0), T=(0,1), G=(1,1), C=(1,0)
  let mut x = 0.5
  let mut y = 0.5
  let n = seq.length()
  let mut i = 0
  while i < n {
    let ch = seq.unsafe_get(i).unsafe_to_char()
    match ch {
      'A' | 'a' => {
        x = (x + 0.0) / 2.0
        y = (y + 0.0) / 2.0
      }
      'T' | 't' => {
        x = (x + 0.0) / 2.0
        y = (y + 1.0) / 2.0
      }
      'G' | 'g' => {
        x = (x + 1.0) / 2.0
        y = (y + 1.0) / 2.0
      }
      'C' | 'c' => {
        x = (x + 1.0) / 2.0
        y = (y + 0.0) / 2.0
      }
      _ => ()
    }
    i = i + 1
  }
  (x, y)
}

///|
/// Calculate the similarity between two sequences based on their signatures.
/// Uses cosine similarity of k-mer frequency vectors.
pub fn sequence_similarity(
  seq1 : String,
  seq2 : String,
  k? : Int = 3,
) -> Double {
  let sig1 = sequence_signature(seq1, k~)
  let sig2 = sequence_signature(seq2, k~)
  // Get all unique k-mers
  let all_kmers : Map[String, Bool] = Map([], capacity=64)
  for kv in sig1 {
    all_kmers[kv.0] = true
  }
  for kv in sig2 {
    all_kmers[kv.0] = true
  }
  // Calculate dot product and magnitudes
  let mut dot = 0.0
  let mut mag1 = 0.0
  let mut mag2 = 0.0
  for kv in all_kmers {
    let kmer = kv.0
    let v1 = sig1.get(kmer).unwrap_or(0.0)
    let v2 = sig2.get(kmer).unwrap_or(0.0)
    dot = dot + v1 * v2
    mag1 = mag1 + v1 * v1
    mag2 = mag2 + v2 * v2
  }
  if mag1 == 0.0 || mag2 == 0.0 {
    return 0.0
  }
  dot / (mag1.sqrt() * mag2.sqrt())
}
