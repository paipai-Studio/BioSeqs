/// biotite_seq_profile.mbt — Sequence profile (PSSM)
/// Inspired by biotite.sequence.profile.SequenceProfile

///|
/// Sequence profile (Position-Specific Scoring Matrix).
/// Each position has a score for each possible symbol.
pub struct BiotiteSequenceProfile {
  symbols : Array[String] // alphabet symbols
  scores : Array[Array[Double]] // position -> symbol -> score
  consensus : Array[Int32] // consensus symbol index per position
  n_sequences : Int32 // number of sequences used to build profile
  gap_threshold : Double // fraction of gaps above which position is ignored
}

///|
/// Build a profile from a multiple sequence alignment.
///
/// Parameters:
/// - alignment: array of aligned sequences (same length, with '-' for gaps)
/// - symbols: the alphabet (e.g. ["A","C","G","T"])
/// - pseudocount: added to each count to avoid zero probabilities
pub fn BiotiteSequenceProfile::from_alignment(
  alignment : Array[String],
  symbols : Array[String],
  pseudocount~ : Double = 0.1,
) -> BiotiteSequenceProfile {
  let n_seqs = alignment.length()
  if n_seqs == 0 {
    return {
      symbols,
      scores: [],
      consensus: [],
      n_sequences: 0,
      gap_threshold: 0.5,
    }
  }
  let aln_len = alignment[0].length()
  let n_sym = symbols.length()
  let scores : Array[Array[Double]] = Array::make(
    aln_len,
    Array::make(n_sym, 0.0),
  )
  let consensus : Array[Int32] = Array::make(aln_len, 0)
  for pos = 0; pos < aln_len; pos = pos + 1 {
    let counts : Array[Double] = Array::make(n_sym, pseudocount)
    let mut gap_count = 0
    for seq_idx = 0; seq_idx < n_seqs; seq_idx = seq_idx + 1 {
      let ch = alignment[seq_idx][pos:pos + 1].to_owned()
      if ch == "-" {
        gap_count = gap_count + 1
      } else {
        // Find symbol index
        for sym_idx = 0; sym_idx < n_sym; sym_idx = sym_idx + 1 {
          if ch == symbols[sym_idx] {
            counts[sym_idx] = counts[sym_idx] + 1.0
            break
          }
        }
      }
    }
    // Convert counts to log-odds scores
    let total = biotite_arr_sum(counts)
    let mut best_idx = 0
    let mut best_count = 0.0
    for sym_idx = 0; sym_idx < n_sym; sym_idx = sym_idx + 1 {
      let freq = counts[sym_idx] / total
      // Background frequency: uniform
      let bg = 1.0 / n_sym.to_double()
      if freq > 0.0 {
        scores[pos][sym_idx] = @math.log2(freq / bg)
      } else {
        scores[pos][sym_idx] = -10.0
      }
      if counts[sym_idx] > best_count {
        best_count = counts[sym_idx]
        best_idx = sym_idx
      }
    }
    consensus[pos] = best_idx
  }
  { symbols, scores, consensus, n_sequences: n_seqs, gap_threshold: 0.5, }
}

///|
/// Score a sequence against the profile.
/// Returns the total log-odds score.
pub fn BiotiteSequenceProfile::score(
  self : BiotiteSequenceProfile,
  sequence : String,
) -> Double {
  let n = sequence.length()
  let mut total = 0.0
  for pos = 0; pos < n && pos < self.scores.length(); pos = pos + 1 {
    let ch = sequence[pos:pos + 1].to_owned()
    for sym_idx = 0; sym_idx < self.symbols.length(); sym_idx = sym_idx + 1 {
      if ch == self.symbols[sym_idx] {
        total = total + self.scores[pos][sym_idx]
        break
      }
    }
  }
  total
}

///|
/// Get the consensus sequence.
pub fn BiotiteSequenceProfile::consensus_sequence(
  self : BiotiteSequenceProfile,
) -> String {
  let mut result = ""
  for pos = 0; pos < self.consensus.length(); pos = pos + 1 {
    result = result + self.symbols[self.consensus[pos]]
  }
  result
}

///|
/// Get the score matrix for a specific position.
pub fn BiotiteSequenceProfile::position_scores(
  self : BiotiteSequenceProfile,
  pos : Int32,
) -> Array[Double] {
  if pos >= 0 && pos < self.scores.length() {
    self.scores[pos]
  } else {
    []
  }
}

///|
/// Get the number of positions in the profile.
pub fn BiotiteSequenceProfile::length(self : BiotiteSequenceProfile) -> Int32 {
  self.scores.length()
}

///|
/// Check if a position should be ignored (too many gaps).
pub fn BiotiteSequenceProfile::is_ignored(
  self : BiotiteSequenceProfile,
  pos : Int32,
  gap_fraction : Array[Double],
) -> Bool {
  if pos >= 0 && pos < gap_fraction.length() {
    gap_fraction[pos] > self.gap_threshold
  } else {
    true
  }
}

///|
fn biotite_arr_sum(arr : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < arr.length(); i = i + 1 {
    total = total + arr[i]
  }
  total
}
