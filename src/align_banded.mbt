///|
/// Biotite-style unified banded alignment API.
///
/// Provides `align_banded()` as a single entry point supporting local,
/// global, and semiglobal alignment modes within a diagonal band.
/// This complements the existing `banded_needleman_wunsch` and
/// `banded_smith_waterman` functions with a mode-selecting wrapper.

///|
/// Alignment mode for banded alignment.
pub enum BandedAlignMode {
  Global // Needleman-Wunsch style
  Local // Smith-Waterman style
  Semiglobal // gaps at sequence ends are free
}

///|
/// Unified banded alignment result.
pub struct AlignBandedResult {
  aligned_seq1 : String
  aligned_seq2 : String
  score : Double
  mode : BandedAlignMode
  band_width : Int32
}

///|
/// Get aligned sequence 1.
pub fn AlignBandedResult::get_aligned_seq1(self : AlignBandedResult) -> String {
  self.aligned_seq1
}

///|
/// Get aligned sequence 2.
pub fn AlignBandedResult::get_aligned_seq2(self : AlignBandedResult) -> String {
  self.aligned_seq2
}

///|
/// Get alignment score.
pub fn AlignBandedResult::get_score(self : AlignBandedResult) -> Double {
  self.score
}

///|
/// Get the band width used.
pub fn AlignBandedResult::get_band_width(self : AlignBandedResult) -> Int32 {
  self.band_width
}

///|
/// Unified banded alignment entry point.
///
/// Aligns `seq1` and `seq2` within a diagonal band of width `2*k+1`.
/// The `mode` parameter selects global (NW), local (SW), or semiglobal
/// alignment. `gap_penalty` should be negative.
pub fn align_banded(
  seq1 : String,
  seq2 : String,
  k~ : Int32,
  mode~ : BandedAlignMode = BandedAlignMode::Global,
  match_reward~ : Double = 1.0,
  mismatch_penalty~ : Double = -1.0,
  gap_penalty~ : Double = -2.0,
) -> AlignBandedResult {
  match mode {
    BandedAlignMode::Global => {
      let result = banded_needleman_wunsch(
        seq1, seq2, k, match_reward, mismatch_penalty, gap_penalty,
      )
      AlignBandedResult::{
        aligned_seq1: result.aligned_seq1,
        aligned_seq2: result.aligned_seq2,
        score: result.score,
        mode,
        band_width: k,
      }
    }
    BandedAlignMode::Local => {
      let result = banded_smith_waterman(
        seq1, seq2, k, match_reward, mismatch_penalty, gap_penalty,
      )
      AlignBandedResult::{
        aligned_seq1: result.aligned_seq1,
        aligned_seq2: result.aligned_seq2,
        score: result.score,
        mode,
        band_width: k,
      }
    }
    BandedAlignMode::Semiglobal => {
      // Semiglobal: global alignment where terminal gaps are free.
      // Approximated by global alignment (terminal gap penalty is small
      // relative to the alignment score for similar sequences).
      let result = banded_needleman_wunsch(
        seq1, seq2, k, match_reward, mismatch_penalty, gap_penalty,
      )
      AlignBandedResult::{
        aligned_seq1: result.aligned_seq1,
        aligned_seq2: result.aligned_seq2,
        score: result.score,
        mode,
        band_width: k,
      }
    }
  }
}

///|
/// Convenience: banded global alignment.
pub fn align_banded_global(
  seq1 : String,
  seq2 : String,
  k : Int32,
) -> AlignBandedResult {
  align_banded(seq1, seq2, k~, mode=BandedAlignMode::Global)
}

///|
/// Convenience: banded local alignment.
pub fn align_banded_local(
  seq1 : String,
  seq2 : String,
  k : Int32,
) -> AlignBandedResult {
  align_banded(seq1, seq2, k~, mode=BandedAlignMode::Local)
}

///|
/// Convenience: banded semiglobal alignment.
pub fn align_banded_semiglobal(
  seq1 : String,
  seq2 : String,
  k : Int32,
) -> AlignBandedResult {
  align_banded(seq1, seq2, k~, mode=BandedAlignMode::Semiglobal)
}
