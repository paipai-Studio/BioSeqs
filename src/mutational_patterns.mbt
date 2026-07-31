///|
/// MutationalPatterns module - Somatic mutation pattern analysis.
/// Inspired by Bioconductor MutationalPatterns package.
///
/// References:
///   - https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html
///   - Blokzijl F et al. (2018) MutationalPatterns: comprehensive genome-wide analysis
///     of mutational processes. Genome Medicine 10:33.

///|
/// The 6 basic mutation types (pyrimidine-based).
pub enum MutationCategory {
  C_A /// C>A transversion
  C_G /// C>G transversion
  C_T /// C>T transition
  T_A /// T>A transversion
  T_C /// T>C transition
  T_G /// T>G transversion
}

///|
/// Create a C>A mutation category.
pub fn mutation_c_a() -> MutationCategory {
  MutationCategory::C_A
}

///|
/// Create a C>G mutation category.
pub fn mutation_c_g() -> MutationCategory {
  MutationCategory::C_G
}

///|
/// Create a C>T mutation category.
pub fn mutation_c_t() -> MutationCategory {
  MutationCategory::C_T
}

///|
/// Create a T>A mutation category.
pub fn mutation_t_a() -> MutationCategory {
  MutationCategory::T_A
}

///|
/// Create a T>C mutation category.
pub fn mutation_t_c() -> MutationCategory {
  MutationCategory::T_C
}

///|
/// Create a T>G mutation category.
pub fn mutation_t_g() -> MutationCategory {
  MutationCategory::T_G
}

///|
/// MutationCategory to string.
pub fn MutationCategory::to_string(self : MutationCategory) -> String {
  match self {
    MutationCategory::C_A => "C>A"
    MutationCategory::C_G => "C>G"
    MutationCategory::C_T => "C>T"
    MutationCategory::T_A => "T>A"
    MutationCategory::T_C => "T>C"
    MutationCategory::T_G => "T>G"
  }
}

///|
/// All 6 mutation types in standard order.
pub fn all_mutation_types() -> Array[MutationCategory] {
  [
    MutationCategory::C_A, MutationCategory::C_G, MutationCategory::C_T,
    MutationCategory::T_A, MutationCategory::T_C, MutationCategory::T_G,
  ]
}

///|
/// The 16 trinucleotide contexts (3' base × 5' base).
/// Each context is of form "N.N" where N is A, C, G, or T.
pub fn trinucleotide_contexts() -> Array[String] {
  let bases = ["A", "C", "G", "T"]
  let contexts : Array[String] = Array::new()
  let mut i = 0
  while i < 4 {
    let mut j = 0
    while j < 4 {
      contexts.push(bases[i] + "." + bases[j])
      j = j + 1
    }
    i = i + 1
  }
  contexts
}

///|
/// A single somatic mutation record.
pub struct SomaticMutation {
  chromosome : String
  position : Int
  ref_base : String
  alt_base : String
  trinucleotide_context : String /// e.g., "ACT" (ref is middle)
}

///|
/// 96-channel mutation count matrix.
/// Rows: 96 trinucleotide-context × mutation-type combinations
/// Columns: samples
pub struct MutationMatrix {
  /// Sample IDs
  sample_ids : Array[String]
  /// 96 channel labels (e.g., "A[C>A]A", "A[C>A]C", ...)
  channel_labels : Array[String]
  /// Count matrix: 96 rows × n_samples columns
  counts : Array[Array[Int]]
  total_mutations : Array[Int] /// per-sample total
}

///|
/// A mutational signature.
pub struct MutSigProfile {
  id : String
  /// 96-channel probability distribution (sums to 1)
  probabilities : Array[Double]
}

///|
/// Create a SomaticMutation.
pub fn SomaticMutation::new(
  chromosome : String,
  position : Int,
  ref_base : String,
  alt_base : String,
  trinucleotide_context : String,
) -> SomaticMutation {
  SomaticMutation::{
    chromosome,
    position,
    ref_base,
    alt_base,
    trinucleotide_context,
  }
}

///|
/// Get all 96 channel labels in standard order.
/// Format: "X[Y>Z]X" where X is flanking, Y is ref, Z is alt.
pub fn get_channel_labels() -> Array[String] {
  let contexts = trinucleotide_contexts()
  let types = all_mutation_types()
  let labels : Array[String] = Array::new()
  // Order: for each mutation type, iterate over 16 contexts
  let mut ti = 0
  while ti < types.length() {
    let type_str = types[ti].to_string()
    let mut ci = 0
    while ci < contexts.length() {
      let ctx = contexts[ci]
      // ctx is like "A.C" meaning A[C]C
      let parts = ctx.split(".")
      let left = parts.next().unwrap_or("A").to_string()
      let right = parts.next().unwrap_or("A").to_string()
      labels.push(left + "[" + type_str + "]" + right)
      ci = ci + 1
    }
    ti = ti + 1
  }
  labels
}

///|
/// Determine mutation type from ref and alt base.
/// Uses pyrimidine-based convention (C/T as reference).
pub fn get_mutation_type(ref_base : String, alt_base : String) -> MutationCategory {
  let r = ref_base.to_upper()
  let a = alt_base.to_upper()
  // If ref is purine (A/G), reverse complement to get pyrimidine-based
  let (r2, a2) = if r == "A" || r == "G" {
    (complement_base(r), complement_base(a))
  } else {
    (r, a)
  }
  if r2 == "C" && a2 == "A" {
    MutationCategory::C_A
  } else if r2 == "C" && a2 == "G" {
    MutationCategory::C_G
  } else if r2 == "C" && a2 == "T" {
    MutationCategory::C_T
  } else if r2 == "T" && a2 == "A" {
    MutationCategory::T_A
  } else if r2 == "T" && a2 == "C" {
    MutationCategory::T_C
  } else if r2 == "T" && a2 == "G" {
    MutationCategory::T_G
  } else {
    MutationCategory::C_A /// default fallback
  }
}

///|
/// Complement a single base.
fn complement_base(base : String) -> String {
  match base.to_upper() {
    "A" => "T"
    "T" => "A"
    "G" => "C"
    "C" => "G"
    _ => "N"
  }
}

///|
/// Normalize a trinucleotide context to pyrimidine-based convention.
/// If middle base is purine (A/G), reverse complement the whole trimer.
pub fn normalize_context(trinucleotide : String) -> String {
  let upper = trinucleotide.to_upper()
  if upper.length() != 3 {
    return upper
  }
  let middle = upper.substring(start=1, end=2)
  if middle == "A" || middle == "G" {
    // Reverse complement
    let b0 = complement_base(upper.substring(start=0, end=1))
    let b1 = complement_base(upper.substring(start=1, end=2))
    let b2 = complement_base(upper.substring(start=2, end=3))
    b2 + b1 + b0
  } else {
    upper
  }
}

///|
/// Get the 0-based channel index (0-95) for a mutation.
/// Based on the trinucleotide context and mutation type.
pub fn get_channel_index(ref_base : String, alt_base : String, trinucleotide : String) -> Int {
  let norm_ctx = normalize_context(trinucleotide)
  let mtype = get_mutation_type(ref_base, alt_base)
  // Find type index (0-5)
  let types = all_mutation_types()
  let mut type_idx = 0
  let mut ti = 0
  while ti < types.length() {
    if types[ti].to_string() == mtype.to_string() {
      type_idx = ti
      break
    }
    ti = ti + 1
  }
  // Find context index (0-15)
  // norm_ctx is like "ACT", middle is the ref base
  let left = norm_ctx.substring(start=0, end=1)
  let right = norm_ctx.substring(start=2, end=3)
  let bases = ["A", "C", "G", "T"]
  let mut left_idx = 0
  let mut right_idx = 0
  let mut bi = 0
  while bi < 4 {
    if bases[bi] == left {
      left_idx = bi
    }
    if bases[bi] == right {
      right_idx = bi
    }
    bi = bi + 1
  }
  let ctx_idx = left_idx * 4 + right_idx
  type_idx * 16 + ctx_idx
}

///|
/// Build a 96-channel mutation count matrix from a list of mutations.
pub fn build_mutation_matrix(
  mutations : Array[SomaticMutation],
  sample_ids : Array[String],
) -> MutationMatrix {
  let n_channels = 96
  let n_samples = sample_ids.length()
  let channel_labels = get_channel_labels()
  // Initialize count matrix
  let counts : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < n_channels {
    let row : Array[Int] = Array::new()
    let mut j = 0
    while j < n_samples {
      row.push(0)
      j = j + 1
    }
    counts.push(row)
    i = i + 1
  }
  // Build sample index lookup
  let sample_idx : Map[String, Int] = Map::new()
  let mut si = 0
  while si < n_samples {
    sample_idx.set(sample_ids[si], si)
    si = si + 1
  }
  // Count mutations
  let mut mi = 0
  while mi < mutations.length() {
    let mut_record = mutations[mi]
    let sidx = sample_idx.get(mut_record.chromosome).unwrap_or(-1)
    if sidx >= 0 {
      let ch_idx = get_channel_index(
        mut_record.ref_base,
        mut_record.alt_base,
        mut_record.trinucleotide_context,
      )
      if ch_idx >= 0 && ch_idx < n_channels {
        counts[ch_idx][sidx] = counts[ch_idx][sidx] + 1
      }
    }
    mi = mi + 1
  }
  // Calculate per-sample totals
  let totals : Array[Int] = Array::new()
  let mut sj = 0
  while sj < n_samples {
    let mut total = 0
    let mut ci = 0
    while ci < n_channels {
      total = total + counts[ci][sj]
      ci = ci + 1
    }
    totals.push(total)
    sj = sj + 1
  }
  MutationMatrix::{
    sample_ids,
    channel_labels,
    counts,
    total_mutations: totals,
  }
}

///|
/// Get a specific count from the matrix.
pub fn MutationMatrix::get_count(
  self : MutationMatrix,
  channel_idx : Int,
  sample_idx : Int,
) -> Int {
  if channel_idx < 0 || channel_idx >= self.counts.length() {
    return 0
  }
  if sample_idx < 0 || sample_idx >= self.counts[channel_idx].length() {
    return 0
  }
  self.counts[channel_idx][sample_idx]
}

///|
/// Convert counts to relative frequencies (normalize per sample).
pub fn MutationMatrix::to_relative_frequencies(self : MutationMatrix) -> Array[Array[Double]] {
  let n_channels = self.counts.length()
  let n_samples = if n_channels > 0 { self.counts[0].length() } else { 0 }
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_channels {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let total = self.total_mutations[j].to_double()
      let val = if total > 0.0 {
        self.counts[i][j].to_double() / total
      } else {
        0.0
      }
      row.push(val)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Compute cosine similarity between two vectors.
pub fn cosine_similarity(vec1 : Array[Double], vec2 : Array[Double]) -> Double {
  let len = if vec1.length() < vec2.length() { vec1.length() } else { vec2.length() }
  let mut dot = 0.0
  let mut norm1 = 0.0
  let mut norm2 = 0.0
  let mut i = 0
  while i < len {
    dot = dot + vec1[i] * vec2[i]
    norm1 = norm1 + vec1[i] * vec1[i]
    norm2 = norm2 + vec2[i] * vec2[i]
    i = i + 1
  }
  let denom = norm1.sqrt() * norm2.sqrt()
  if denom < 1.0e-10 {
    return 0.0
  }
  dot / denom
}

///|
/// Compute cosine similarity between two samples in a mutation matrix.
pub fn MutationMatrix::sample_cosine_similarity(
  self : MutationMatrix,
  sample_idx1 : Int,
  sample_idx2 : Int,
) -> Double {
  let n_channels = self.counts.length()
  let vec1 : Array[Double] = Array::new()
  let vec2 : Array[Double] = Array::new()
  let mut i = 0
  while i < n_channels {
    vec1.push(self.counts[i][sample_idx1].to_double())
    vec2.push(self.counts[i][sample_idx2].to_double())
    i = i + 1
  }
  cosine_similarity(vec1, vec2)
}

///|
/// Calculate signature exposure (contribution) of a known signature to a sample.
/// Returns the count attributed to the signature.
pub fn calculate_exposure(
  sample_profile : Array[Double],
  signature : MutSigProfile,
) -> Double {
  // Project sample profile onto signature
  let sig = signature.probabilities
  let dot = cosine_similarity(sample_profile, sig)
  // Total mutations * cosine similarity gives approximate exposure
  let mut total = 0.0
  let mut i = 0
  while i < sample_profile.length() {
    total = total + sample_profile[i]
    i = i + 1
  }
  total * dot
}

///|
/// Compute the contribution of multiple signatures to a sample.
/// Returns an array of contribution fractions.
pub fn fit_signatures(
  sample_profile : Array[Double],
  signatures : Array[MutSigProfile],
) -> Array[Double] {
  let n_sigs = signatures.length()
  if n_sigs == 0 {
    return Array::new()
  }
  // Calculate cosine similarity of each signature to the sample
  let sims : Array[Double] = Array::new()
  let mut total_sim = 0.0
  let mut i = 0
  while i < n_sigs {
    let sim = cosine_similarity(sample_profile, signatures[i].probabilities)
    if sim < 0.0 {
      sims.push(0.0)
    } else {
      sims.push(sim)
    }
    total_sim = total_sim + sims[i]
    i = i + 1
  }
  // Normalize to get fractions
  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < n_sigs {
    if total_sim > 0.0 {
      result.push(sims[j] / total_sim)
    } else {
      result.push(1.0 / n_sigs.to_double())
    }
    j = j + 1
  }
  result
}

///|
/// Create a sample mutational signature from a probability vector.
pub fn MutSigProfile::new(
  id : String,
  probabilities : Array[Double],
) -> MutSigProfile {
  MutSigProfile::{ id, probabilities }
}

///|
/// Create known COSMIC signatures (simplified versions of SBS1, SBS2, SBS5).
/// These are approximated for demonstration.
pub fn known_cosmic_signatures() -> Array[MutSigProfile] {
  let sigs : Array[MutSigProfile] = Array::new()
  // SBS1: spontaneous deamination of 5-methylcytosine (C>T at NpCpG)
  let sbs1 : Array[Double] = Array::new()
  let mut i = 0
  while i < 96 {
    // C>T mutations at NCG contexts (channels 32-47, specifically NCG)
    // C>T is type index 2, context index for CG is 1*4+1=5 (C left, G right)
    // So channel = 2*16 + 5 = 37
    if i == 37 {
      sbs1.push(0.8)
    } else if i >= 32 && i < 48 {
      sbs1.push(0.01)
    } else {
      sbs1.push(0.001)
    }
    i = i + 1
  }
  sigs.push(MutSigProfile::new("SBS1", sbs1))
  // SBS5: unknown, flat distribution
  let sbs5 : Array[Double] = Array::new()
  let mut j = 0
  while j < 96 {
    sbs5.push(1.0 / 96.0)
    j = j + 1
  }
  sigs.push(MutSigProfile::new("SBS5", sbs5))
  sigs
}

///|
/// Generate a mutation spectrum summary string.
pub fn MutationMatrix::spectrum_summary(self : MutationMatrix) -> String {
  let types = all_mutation_types()
  let mut s = "=== Mutation Spectrum Summary ===\n".to_string()
  let mut ti = 0
  while ti < types.length() {
    let type_str = types[ti].to_string()
    let mut total = 0
    let mut si = 0
    while si < self.sample_ids.length() {
      let mut ci = 0
      while ci < 16 {
        total = total + self.counts[ti * 16 + ci][si]
        ci = ci + 1
      }
      si = si + 1
    }
    s = s + type_str + ": " + total.to_string() + "\n"
    ti = ti + 1
  }
  s
}

///|
/// Generate sample data for testing.
pub fn mutational_patterns_sample_data() -> Array[SomaticMutation] {
  let mutations : Array[SomaticMutation] = Array::new()
  // Create some test mutations in sample "S1"
  // C>T at ACG context
  mutations.push(SomaticMutation::new("S1", 100, "C", "T", "ACT"))
  mutations.push(SomaticMutation::new("S1", 200, "C", "T", "ACG"))
  mutations.push(SomaticMutation::new("S1", 300, "C", "A", "CCA"))
  mutations.push(SomaticMutation::new("S1", 400, "T", "C", "ATC"))
  mutations.push(SomaticMutation::new("S1", 500, "C", "G", "GCG"))
  // Create some test mutations in sample "S2"
  mutations.push(SomaticMutation::new("S2", 100, "C", "A", "TCA"))
  mutations.push(SomaticMutation::new("S2", 200, "C", "T", "TCG"))
  mutations.push(SomaticMutation::new("S2", 300, "T", "G", "TTG"))
  mutations.push(SomaticMutation::new("S2", 400, "G", "A", "AGT"))
  mutations
}
