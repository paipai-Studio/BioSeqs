// Scikit-bio third-round: sequence simulation and alignment-free distances (items 20-23)
// simulate_seq_on_tree, feature_frequency_profile, composition_vector_distance, return_time_distance

///|
pub fn skbio3_simulate_seq_on_tree(
  root_seq : String,
  branch_length : Double,
  mutation_rate : Double,
) -> String {
  // Simulate sequence evolution along a single branch
  // Jukes-Cantor model: each site has probability p = 3/4(1 - exp(-4*mu*t/3)) of changing
  let p_change = 0.75 * (1.0 - @math.exp(-4.0 * mutation_rate * branch_length / 3.0))
  let n = root_seq.length()
  let sb = StringBuilder()
  for i in 0..<n {
    let c = Int::unsafe_to_char(root_seq[i].to_int())
    // With probability p_change, mutate to a random different nucleotide
    let r = skbio3_pseudo_random()
    if r < p_change {
      // Mutate to one of the 3 other nucleotides
      let bases = ['A', 'C', 'G', 'T']
      let mut chosen = c
      let r2 = skbio3_pseudo_random()
      let idx = (r2 * 3.0).floor().to_int()
      let mut count = 0
      for bi in 0..<4 {
        if bases[bi] != c {
          if count == idx {
            chosen = bases[bi]
            break
          }
          count = count + 1
        }
      }
      let _ = sb.write_char(chosen)
    } else {
      let _ = sb.write_char(c)
    }
  }
  sb.to_string()
}

fn skbio3_pseudo_random() -> Double {
  // Simple LCG pseudo-random number generator for simulation
  // Uses a global mutable state
  skbio3_rng_state[0] = skbio3_rng_state[0] * 1103515245 + 12345
  let r = (skbio3_rng_state[0] % 2147483647).to_double() / 2147483647.0
  if r < 0.0 { 0.0 - r } else { r }
}

let skbio3_rng_state : Array[Int] = [42]

///|
pub fn skbio3_feature_frequency_profile(
  seq : String,
  k : Int,
) -> Map[String, Double] {
  // Feature Frequency Profile (FFP): k-mer frequency distribution
  // Returns normalized k-mer frequencies
  let n = seq.length()
  let profile : Map[String, Double] = Map([])
  if n < k || k <= 0 {
    return profile
  }
  let mut total = 0
  let n_minus_k = n - k + 1
  for i in 0..<n_minus_k {
    let kmer = seq[i:i + k].to_owned()
    let current = profile.get(kmer).unwrap_or(0.0)
    profile[kmer] = current + 1.0
    total = total + 1
  }
  if total > 0 {
    profile.each(fn(key, value) {
      profile[key] = value / total.to_double()
    })
  }
  profile
}

///|
pub fn skbio3_ffp_distance(
  seq1 : String,
  seq2 : String,
  k : Int,
) -> Double {
  // FFP distance: Jensen-Shannon divergence between k-mer profiles
  let p1 = skbio3_feature_frequency_profile(seq1, k)
  let p2 = skbio3_feature_frequency_profile(seq2, k)
  // Collect all k-mers
  let all_kmers : Array[String] = []
  p1.each(fn(key, _) {
    all_kmers.push(key)
  })
  p2.each(fn(key, _) {
    let mut found = false
    for i in 0..<all_kmers.length() {
      if all_kmers[i] == key {
        found = true
        break
      }
    }
    if not(found) {
      all_kmers.push(key)
    }
  })
  // Compute JS divergence
  let mut js_div = 0.0
  for i in 0..<all_kmers.length() {
    let kmer = all_kmers[i]
    let pi = p1.get(kmer).unwrap_or(0.0)
    let qi = p2.get(kmer).unwrap_or(0.0)
    let mi = (pi + qi) / 2.0
    if mi > 0.0 {
      if pi > 0.0 {
        js_div = js_div + 0.5 * pi * @math.ln(pi / mi)
      }
      if qi > 0.0 {
        js_div = js_div + 0.5 * qi * @math.ln(qi / mi)
      }
    }
  }
  js_div.sqrt()
}

///|
pub fn skbio3_composition_vector_distance(
  seq1 : String,
  seq2 : String,
  k : Int,
) -> Double {
  // Composition Vector (CV) distance: remove compositional bias, compare residual vectors
  let p1 = skbio3_feature_frequency_profile(seq1, k)
  let p2 = skbio3_feature_frequency_profile(seq2, k)
  // Collect all k-mers
  let all_kmers : Array[String] = []
  p1.each(fn(key, _) {
    all_kmers.push(key)
  })
  p2.each(fn(key, _) {
    let mut found = false
    for i in 0..<all_kmers.length() {
      if all_kmers[i] == key {
        found = true
        break
      }
    }
    if not(found) {
      all_kmers.push(key)
    }
  })
  // CV: a_i = f_i - product of marginal frequencies (composition-free residual)
  // Simplified: use raw frequencies minus uniform expectation
  let n_kmers = all_kmers.length().to_double()
  let uniform = if n_kmers > 0.0 { 1.0 / n_kmers } else { 0.0 }
  let mut dot_product = 0.0
  let mut norm1 = 0.0
  let mut norm2 = 0.0
  for i in 0..<all_kmers.length() {
    let kmer = all_kmers[i]
    let pi = p1.get(kmer).unwrap_or(0.0) - uniform
    let qi = p2.get(kmer).unwrap_or(0.0) - uniform
    dot_product = dot_product + pi * qi
    norm1 = norm1 + pi * pi
    norm2 = norm2 + qi * qi
  }
  if norm1 <= 0.0 || norm2 <= 0.0 {
    return 1.0
  }
  // Cosine distance
  let cosine = dot_product / (norm1.sqrt() * norm2.sqrt())
  1.0 - cosine
}

///|
pub fn skbio3_return_time_distance(
  seq1 : String,
  seq2 : String,
  k : Int,
) -> Double {
  // Return Time Distance (RTD): based on first return times of k-mers
  // RTD measures the average time for a k-mer to reappear in the sequence
  let rtd1 = skbio3_return_time_vector(seq1, k)
  let rtd2 = skbio3_return_time_vector(seq2, k)
  // Euclidean distance between RTD vectors
  let all_kmers : Array[String] = []
  rtd1.each(fn(key, _) {
    all_kmers.push(key)
  })
  rtd2.each(fn(key, _) {
    let mut found = false
    for i in 0..<all_kmers.length() {
      if all_kmers[i] == key {
        found = true
        break
      }
    }
    if not(found) {
      all_kmers.push(key)
    }
  })
  let mut sum_sq = 0.0
  for i in 0..<all_kmers.length() {
    let kmer = all_kmers[i]
    let r1 = rtd1.get(kmer).unwrap_or(0.0)
    let r2 = rtd2.get(kmer).unwrap_or(0.0)
    let diff = r1 - r2
    sum_sq = sum_sq + diff * diff
  }
  sum_sq.sqrt()
}

fn skbio3_return_time_vector(
  seq : String,
  k : Int,
) -> Map[String, Double] {
  // For each k-mer, compute the mean return time (distance between consecutive occurrences)
  let n = seq.length()
  let rtd : Map[String, Double] = Map([])
  if n < k || k <= 0 {
    return rtd
  }
  // Track positions of each k-mer
  let positions : Map[String, Array[Int]] = Map([])
  let n_minus_k = n - k + 1
  for i in 0..<n_minus_k {
    let kmer = seq[i:i + k].to_owned()
    let pos_list = positions.get(kmer).unwrap_or([])
    pos_list.push(i)
    positions[kmer] = pos_list
  }
  // Compute mean return time for each k-mer
  positions.each(fn(key, pos_list) {
    if pos_list.length() > 1 {
      let mut sum_rts = 0.0
      for j in 1..<pos_list.length() {
        sum_rts = sum_rts + (pos_list[j] - pos_list[j - 1]).to_double()
      }
      rtd[key] = sum_rts / (pos_list.length() - 1).to_double()
    } else {
      rtd[key] = n.to_double()  // single occurrence → large return time
    }
  })
  rtd
}
