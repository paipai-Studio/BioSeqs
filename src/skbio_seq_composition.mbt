/// Sequence composition: dinucleotide, trinucleotide, and codon frequencies.
///
/// scikit-bio equivalent: sequence composition analysis
///
/// These functions compute k-mer frequency distributions for
/// k=2 (dinucleotide), k=3 (trinucleotide/codon), and general k.

///|
/// Compute dinucleotide frequencies.
/// Returns a map from dinucleotide (e.g., "AT") to count.
pub fn skbio_dinucleotide_counts(seq : String) -> Map[String, Int32] {
  let result : Map[String, Int32] = Map([])
  if seq.length() < 2 {
    return result
  }
  for i = 0; i < seq.length() - 1; i = i + 1 {
    let di = seq.substring(start=i, end=i + 2)
    match result.get(di) {
      Some(c) => result[di] = c + 1
      None => result[di] = 1
    }
  }
  result
}

///|
/// Compute dinucleotide frequencies as normalized values (0.0 to 1.0).
pub fn skbio_dinucleotide_frequencies(seq : String) -> Map[String, Double] {
  let counts = skbio_dinucleotide_counts(seq)
  let total = (seq.length() - 1).to_double()
  let result : Map[String, Double] = Map([])
  if total <= 0.0 {
    return result
  }
  let keys = counts.keys().to_array()
  for k in keys {
    result[k] = counts[k].to_double() / total
  }
  result
}

///|
/// Compute trinucleotide frequencies.
/// Returns a map from trinucleotide (e.g., "ATG") to count.
pub fn skbio_trinucleotide_counts(seq : String) -> Map[String, Int32] {
  let result : Map[String, Int32] = Map([])
  if seq.length() < 3 {
    return result
  }
  for i = 0; i < seq.length() - 2; i = i + 1 {
    let tri = seq.substring(start=i, end=i + 3)
    match result.get(tri) {
      Some(c) => result[tri] = c + 1
      None => result[tri] = 1
    }
  }
  result
}

///|
/// Compute trinucleotide frequencies as normalized values.
pub fn skbio_trinucleotide_frequencies(seq : String) -> Map[String, Double] {
  let counts = skbio_trinucleotide_counts(seq)
  let total = (seq.length() - 2).to_double()
  let result : Map[String, Double] = Map([])
  if total <= 0.0 {
    return result
  }
  let keys = counts.keys().to_array()
  for k in keys {
    result[k] = counts[k].to_double() / total
  }
  result
}

///|
/// Compute codon frequencies from a CDS sequence.
/// The sequence length should be a multiple of 3.
/// Returns a map from codon to count.
pub fn skbio_codon_counts(seq : String) -> Map[String, Int32] {
  let result : Map[String, Int32] = Map([])
  let n = seq.length()
  let mut i = 0
  while i + 3 <= n {
    let codon = seq.substring(start=i, end=i + 3)
    match result.get(codon) {
      Some(c) => result[codon] = c + 1
      None => result[codon] = 1
    }
    i = i + 3
  }
  result
}

///|
/// Compute codon frequencies as normalized values (per codon).
pub fn skbio_codon_frequencies(seq : String) -> Map[String, Double] {
  let counts = skbio_codon_counts(seq)
  let total = (seq.length() / 3).to_double()
  let result : Map[String, Double] = Map([])
  if total <= 0.0 {
    return result
  }
  let keys = counts.keys().to_array()
  for k in keys {
    result[k] = counts[k].to_double() / total
  }
  result
}

///|
/// Compute GC content at codon positions 1, 2, 3 separately.
/// Returns (GC1, GC2, GC3) as fractions.
pub fn skbio_gc_by_codon_position(seq : String) -> (Double, Double, Double) {
  let n = seq.length()
  let mut gc1 = 0.0
  let mut gc2 = 0.0
  let mut gc3 = 0.0
  let mut count1 = 0.0
  let mut count2 = 0.0
  let mut count3 = 0.0
  let mut i = 0
  while i + 3 <= n {
    let c1 = seq.unsafe_get(i).unsafe_to_char()
    let c2 = seq.unsafe_get(i + 1).unsafe_to_char()
    let c3 = seq.unsafe_get(i + 2).unsafe_to_char()
    if c1 == 'G' || c1 == 'C' || c1 == 'g' || c1 == 'c' {
      gc1 = gc1 + 1.0
    }
    if c2 == 'G' || c2 == 'C' || c2 == 'g' || c2 == 'c' {
      gc2 = gc2 + 1.0
    }
    if c3 == 'G' || c3 == 'C' || c3 == 'g' || c3 == 'c' {
      gc3 = gc3 + 1.0
    }
    count1 = count1 + 1.0
    count2 = count2 + 1.0
    count3 = count3 + 1.0
    i = i + 3
  }
  let r1 = if count1 > 0.0 { gc1 / count1 } else { 0.0 }
  let r2 = if count2 > 0.0 { gc2 / count2 } else { 0.0 }
  let r3 = if count3 > 0.0 { gc3 / count3 } else { 0.0 }
  (r1, r2, r3)
}

///|
/// Compute general k-mer frequency map.
pub fn skbio_kmer_counts(seq : String, k : Int32) -> Map[String, Int32] {
  let result : Map[String, Int32] = Map([])
  let n = seq.length()
  if k <= 0 || k > n {
    return result
  }
  for i = 0; i <= n - k; i = i + 1 {
    let kmer = seq.substring(start=i, end=i + k)
    match result.get(kmer) {
      Some(c) => result[kmer] = c + 1
      None => result[kmer] = 1
    }
  }
  result
}

///|
/// Compute nucleotide composition (A, C, G, T counts).
/// Returns (A, C, G, T) counts.
pub fn skbio_nucleotide_counts(seq : String) -> (Int32, Int32, Int32, Int32) {
  let mut a : Int32 = 0
  let mut c : Int32 = 0
  let mut g : Int32 = 0
  let mut t : Int32 = 0
  for ch in seq {
    match ch {
      'A' | 'a' => a = a + 1
      'C' | 'c' => c = c + 1
      'G' | 'g' => g = g + 1
      'T' | 't' | 'U' | 'u' => t = t + 1
      _ => ()
    }
  }
  (a, c, g, t)
}

///|
/// Compute the odds ratio for a dinucleotide.
/// Odds(XY) = f(XY) / (f(X) * f(Y))
/// where f is the frequency.
pub fn skbio_dinucleotide_odds_ratio(seq : String, di : String) -> Double {
  if di.length() != 2 {
    return 0.0
  }
  let di_freq = skbio_dinucleotide_frequencies(seq)
  let (a, c, g, t) = skbio_nucleotide_counts(seq)
  let total = (a + c + g + t).to_double()
  if total == 0.0 {
    return 0.0
  }
  let x = di.substring(start=0, end=1)
  let y = di.substring(start=1, end=2)
  let x_count = match x {
    "A" | "a" => a
    "C" | "c" => c
    "G" | "g" => g
    "T" | "t" | "U" | "u" => t
    _ => 0
  }
  let y_count = match y {
    "A" | "a" => a
    "C" | "c" => c
    "G" | "g" => g
    "T" | "t" | "U" | "u" => t
    _ => 0
  }
  let x_freq = x_count.to_double() / total
  let y_freq = y_count.to_double() / total
  let xy_freq = match di_freq.get(di) {
    Some(f) => f
    None => 0.0
  }
  if x_freq > 0.0 && y_freq > 0.0 {
    xy_freq / (x_freq * y_freq)
  } else {
    0.0
  }
}
