///|
/// Needleman-Wunsch global sequence alignment via dynamic programming.
///
/// References:
/// - Needleman SB, Wunsch CD (1970). "A general method applicable to the
///   search for similarities in the amino acid sequence of two proteins."
///   J Mol Biol. 48(3):443-53.
/// - Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis.
///   Cambridge University Press, 1998. Chapter 2.

///|
/// Traceback direction codes for the NW DP matrix.
///
/// Unlike Smith-Waterman, there is no `Stop` variant because global
/// alignment always traces back from (n, m) to (0, 0).
#derive(Eq, Debug)
pub enum NWTrace {
  Diag // match/mismatch: comes from (i-1, j-1)
  Up // gap in seq2 (deletion): comes from (i-1, j)
  Left // gap in seq1 (insertion): comes from (i, j-1)
  Start // only used at (0, 0) to terminate traceback
}

///|
pub extend NWTrace with Eq::{equal, not_equal}

///|
pub extend NWTrace with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Scoring scheme for global alignment.
///
/// Reuses the same structure concept as Smith-Waterman's `Scoring`
/// but defined independently to avoid cross-module coupling.
pub struct NWScoring {
  match_reward : Double
  mismatch_penalty : Double
  gap_penalty : Double
}

///|
/// Create a scoring scheme.
pub fn NWScoring::new(
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> NWScoring {
  NWScoring::{ match_reward, mismatch_penalty, gap_penalty, }
}

///|
/// Default DNA scoring: match=2, mismatch=-1, gap=-1.
pub fn NWScoring::default_dna() -> NWScoring {
  NWScoring::{ match_reward: 2.0, mismatch_penalty: -1.0, gap_penalty: -1.0, }
}

///|
/// Default protein scoring: match=5, mismatch=-4, gap=-3.
pub fn NWScoring::default_protein() -> NWScoring {
  NWScoring::{ match_reward: 5.0, mismatch_penalty: -4.0, gap_penalty: -3.0, }
}

///|
/// Full result of a Needleman-Wunsch alignment.
pub struct NWResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  score_matrix : Array[Array[Double]]
  trace_matrix : Array[Array[NWTrace]]
}

///|
/// Get aligned sequence 1.
pub fn NWResult::aligned1(self : NWResult) -> String {
  self.aligned_seq1
}

///|
/// Get aligned sequence 2.
pub fn NWResult::aligned2(self : NWResult) -> String {
  self.aligned_seq2
}

///|
/// Get match line ('|' for match, '.' for mismatch, ' ' for gap).
pub fn NWResult::match_line(self : NWResult) -> String {
  self.match_line
}

///|
/// Get alignment score.
pub fn NWResult::score(self : NWResult) -> Double {
  self.score
}

///|
/// Get the full score matrix for inspection.
pub fn NWResult::score_matrix(self : NWResult) -> Array[Array[Double]] {
  self.score_matrix
}

///|
/// Get the full traceback matrix for inspection.
pub fn NWResult::trace_matrix(self : NWResult) -> Array[Array[NWTrace]] {
  self.trace_matrix
}

///|
/// Get alignment length (both aligned sequences have the same length).
pub fn NWResult::alignment_length(self : NWResult) -> Int32 {
  self.aligned_seq1.length()
}

///|
/// Count the number of matches in the alignment.
pub fn NWResult::num_matches(self : NWResult) -> Int32 {
  let mut count = 0
  for i = 0; i < self.match_line.length(); i = i + 1 {
    if self.match_line.unsafe_get(i).unsafe_to_char() == '|' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of mismatches (excluding gaps).
pub fn NWResult::num_mismatches(self : NWResult) -> Int32 {
  let mut count = 0
  for i = 0; i < self.match_line.length(); i = i + 1 {
    if self.match_line.unsafe_get(i).unsafe_to_char() == '.' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of gap positions in sequence 1.
pub fn NWResult::num_gaps_seq1(self : NWResult) -> Int32 {
  let mut count = 0
  for i = 0; i < self.aligned_seq1.length(); i = i + 1 {
    if self.aligned_seq1.unsafe_get(i).unsafe_to_char() == '-' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of gap positions in sequence 2.
pub fn NWResult::num_gaps_seq2(self : NWResult) -> Int32 {
  let mut count = 0
  for i = 0; i < self.aligned_seq2.length(); i = i + 1 {
    if self.aligned_seq2.unsafe_get(i).unsafe_to_char() == '-' {
      count = count + 1
    }
  }
  count
}

///|
/// Calculate percent identity (0.0 – 100.0).
pub fn NWResult::percent_identity(self : NWResult) -> Double {
  let aln_len = self.alignment_length()
  if aln_len == 0 {
    return 0.0
  }
  self.num_matches().to_double() / aln_len.to_double() * 100.0
}

///|
/// Format the alignment as a readable block (3 lines).
pub fn NWResult::format(self : NWResult) -> String {
  let buf = StringBuilder()
  buf.write_string(self.aligned_seq1)
  buf.write_char('\n')
  buf.write_string(self.match_line)
  buf.write_char('\n')
  buf.write_string(self.aligned_seq2)
  buf.to_string()
}

///|
/// Format the alignment with line wrapping for long sequences.
/// `width` — number of characters per line.
pub fn NWResult::format_wrapped(self : NWResult, width : Int32) -> String {
  let buf = StringBuilder()
  let total = self.aligned_seq1.length()
  let mut start = 0
  let mut block = 1
  while start < total {
    let end = if start + width > total { total } else { start + width }
    buf.write_string("Block ")
    buf.write_string(block.to_string())
    buf.write_string(" (positions ")
    buf.write_string(start.to_string())
    buf.write_char('-')
    buf.write_string((end - 1).to_string())
    buf.write_string("):\n")
    buf.write_string(self.aligned_seq1[start:end].to_owned())
    buf.write_char('\n')
    buf.write_string(self.match_line[start:end].to_owned())
    buf.write_char('\n')
    buf.write_string(self.aligned_seq2[start:end].to_owned())
    buf.write_char('\n')
    if end < total {
      buf.write_char('\n')
    }
    start = end
    block = block + 1
  }
  buf.to_string()
}

///|
/// Run Needleman-Wunsch global alignment.
///
/// `seq1` — first sequence (string).
/// `seq2` — second sequence (string).
/// `scoring` — scoring scheme (use `NWScoring::default_dna()` etc.).
///
/// Returns an `NWResult` containing aligned sequences, score,
/// and the full DP matrices.
pub fn needleman_wunsch(
  seq1 : String,
  seq2 : String,
  scoring : NWScoring,
) -> NWResult {
  let out = align_core_simple(
    seq1,
    seq2,
    scoring.match_reward,
    scoring.mismatch_penalty,
    scoring.gap_penalty,
    AlignMode::Global,
  )
  let trace : Array[Array[NWTrace]] = out.trace_matrix.map(fn(row) {
    row.map(align_trace_to_nw)
  })
  NWResult::{
    aligned_seq1: out.aligned_seq1,
    aligned_seq2: out.aligned_seq2,
    match_line: out.match_line,
    score: out.score,
    score_matrix: out.score_matrix,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Needleman-Wunsch with default DNA scoring.
pub fn needleman_wunsch_dna(seq1 : String, seq2 : String) -> NWResult {
  needleman_wunsch(seq1, seq2, NWScoring::default_dna())
}

///|
/// Convenience: run Needleman-Wunsch with default protein scoring.
pub fn needleman_wunsch_protein(seq1 : String, seq2 : String) -> NWResult {
  needleman_wunsch(seq1, seq2, NWScoring::default_protein())
}

///|
/// Run Needleman-Wunsch global alignment using a substitution matrix.
///
/// `seq1`, `seq2` — protein sequences (single-letter code strings).
/// `matrix` — substitution matrix (e.g., `@data.blosum62()`).
/// `gap_penalty` — penalty for a gap (typically negative).
///
/// Returns an `NWResult` containing aligned sequences, score,
/// and the full DP matrices.
pub fn needleman_wunsch_with_matrix(
  seq1 : String,
  seq2 : String,
  matrix : @data.SubstitutionMatrix,
  gap_penalty : Double,
) -> NWResult {
  let out = align_core_matrix(
    seq1,
    seq2,
    matrix,
    gap_penalty,
    AlignMode::Global,
  )
  let trace : Array[Array[NWTrace]] = out.trace_matrix.map(fn(row) {
    row.map(align_trace_to_nw)
  })
  NWResult::{
    aligned_seq1: out.aligned_seq1,
    aligned_seq2: out.aligned_seq2,
    match_line: out.match_line,
    score: out.score,
    score_matrix: out.score_matrix,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Needleman-Wunsch with BLOSUM62 matrix.
pub fn needleman_wunsch_blosum62(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> NWResult {
  needleman_wunsch_with_matrix(seq1, seq2, @data.blosum62(), gap_penalty)
}

///|
/// Convenience: run Needleman-Wunsch with PAM250 matrix.
pub fn needleman_wunsch_pam250(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> NWResult {
  needleman_wunsch_with_matrix(seq1, seq2, @data.pam250(), gap_penalty)
}

///|
/// Convert an AlignTrace code to a Needleman-Wunsch NWTrace.
fn align_trace_to_nw(t : AlignTrace) -> NWTrace {
  match t {
    AlignTrace::Stop => NWTrace::Start
    AlignTrace::Diag => NWTrace::Diag
    AlignTrace::Up => NWTrace::Up
    AlignTrace::Left => NWTrace::Left
  }
}
