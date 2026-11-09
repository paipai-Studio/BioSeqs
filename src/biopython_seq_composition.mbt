/// biopython_seq_composition.mbt — Sequence composition analysis
/// Inspired by Biopython Bio.SeqUtils: GC content, GC skew, base frequency, N-gram frequency

///|
/// Calculate GC content (fraction) of a DNA sequence.
pub fn biopython_gc_content(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let mut gc = 0
  for i = 0; i < n; i = i + 1 {
    let ch = sequence[i].to_int()
    // G=71, C=67, g=103, c=99, S=83 (G or C), s=115
    if ch == 71 || ch == 67 || ch == 103 || ch == 99 || ch == 83 || ch == 115 {
      gc = gc + 1
    }
  }
  gc.to_double() / n.to_double() * 100.0
}

///|
/// Calculate GC skew at each position (cumulative).
/// GC skew = (G - C) / (G + C)
pub fn biopython_gc_skew(sequence : String) -> Array[Double] {
  let n = sequence.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut g_count = 0
  let mut c_count = 0
  for i = 0; i < n; i = i + 1 {
    let ch = sequence[i].to_int()
    if ch == 71 || ch == 103 {
      g_count = g_count + 1
    }
    if ch == 67 || ch == 99 {
      c_count = c_count + 1
    }
    let total = g_count + c_count
    result[i] = if total > 0 {
      (g_count - c_count).to_double() / total.to_double()
    } else {
      0.0
    }
  }
  result
}

///|
/// Calculate AT skew at each position (cumulative).
pub fn biopython_at_skew(sequence : String) -> Array[Double] {
  let n = sequence.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut a_count = 0
  let mut t_count = 0
  for i = 0; i < n; i = i + 1 {
    let ch = sequence[i].to_int()
    if ch == 65 || ch == 97 {
      a_count = a_count + 1
    }
    if ch == 84 || ch == 116 {
      t_count = t_count + 1
    }
    let total = a_count + t_count
    result[i] = if total > 0 {
      (a_count - t_count).to_double() / total.to_double()
    } else {
      0.0
    }
  }
  result
}

///|
/// Base frequency composition: returns map of base -> count.
pub fn biopython_base_frequency(sequence : String) -> Array[(String, Int32)] {
  let bases = ["A", "C", "G", "T", "U", "N"]
  let counts : Array[Int32] = Array::make(bases.length(), 0)
  for i = 0; i < sequence.length(); i = i + 1 {
    let ch = sequence[i:i + 1].to_owned().to_upper()
    for j = 0; j < bases.length(); j = j + 1 {
      if ch == bases[j] {
        counts[j] = counts[j] + 1
        break
      }
    }
  }
  let result : Array[(String, Int32)] = Array(capacity=0)
  for j = 0; j < bases.length(); j = j + 1 {
    if counts[j] > 0 {
      result.push((bases[j], counts[j]))
    }
  }
  result
}

///|
/// N-gram (k-mer) frequency table.
pub fn biopython_ngram_frequency(
  sequence : String,
  n : Int32,
) -> Array[(String, Int32)] {
  let len = sequence.length()
  if len < n {
    return []
  }
  let freqs : Array[(String, Int32)] = Array(capacity=0)
  for i = 0; i <= len - n; i = i + 1 {
    let word = sequence[i:i + n].to_owned()
    let mut found = false
    for j = 0; j < freqs.length(); j = j + 1 {
      if freqs[j].0 == word {
        let (w, c) = freqs[j]
        freqs[j] = (w, c + 1)
        found = true
        break
      }
    }
    if !found {
      freqs.push((word, 1))
    }
  }
  // Sort by count descending
  freqs.sort_by(fn(a, b) {
    if b.1 < a.1 {
      -1
    } else if b.1 > a.1 {
      1
    } else {
      0
    }
  })
  freqs
}

///|
/// Purine-pyrimidine ratio.
pub fn biopython_purine_pyrimidine_ratio(sequence : String) -> Double {
  let mut purines = 0
  let mut pyrimidines = 0
  for i = 0; i < sequence.length(); i = i + 1 {
    let ch = sequence[i].to_int()
    // A=65, G=71, a=97, g=103 (purines)
    if ch == 65 || ch == 71 || ch == 97 || ch == 103 {
      purines = purines + 1
    }
    // C=67, T=84, U=85, c=99, t=116, u=117 (pyrimidines)
    if ch == 67 || ch == 84 || ch == 85 || ch == 99 || ch == 116 || ch == 117 {
      pyrimidines = pyrimidines + 1
    }
  }
  if pyrimidines == 0 {
    0.0
  } else {
    purines.to_double() / pyrimidines.to_double()
  }
}
