/// Biotite `sequence.trim` — trimming of primer regions from
/// single-strand or double-strand sequencing reads.
///
/// The primer is matched against both ends of the sequence; the best
/// end-overlap is used to determine how much to cut off from that end.

///|
/// Overlap match evaluation at one end of the sequence.
priv struct BiotiteTrimMatch {
  trim_amount : Int32
  score : Double
}

///|
/// Evaluate the best primer overlap at the left end (primer tail
/// against sequence head) and the right end (primer head against
/// sequence tail).
///
/// A valid overlap needs at least `min_overlap` aligned pairs and at
/// most `max_mismatch` mismatches; the score is
/// `matches - distance_weight * mismatches`.
fn biotite_trim_best_end(
  sequence : String,
  primer : String,
  left : Bool,
  min_overlap : Int32,
  max_mismatch : Int32,
  distance_weight : Double,
) -> BiotiteTrimMatch? {
  let n = sequence.length()
  let m = primer.length()
  let max_o = if m < n { m } else { n }
  let mut best : BiotiteTrimMatch? = None
  for o = min_overlap; o <= max_o; o = o + 1 {
    if o < min_overlap {
      continue
    }
    let mut mm = 0
    for k = 0; k < o; k = k + 1 {
      let sc = if left { sequence[k] } else { sequence[n - o + k] }
      let pc = if left { primer[m - o + k] } else { primer[k] }
      if sc != pc {
        mm = mm + 1
        if mm > max_mismatch {
          break
        }
      }
    }
    if mm > max_mismatch {
      continue
    }
    let matches = (o - mm).to_double()
    let score = matches - distance_weight * mm.to_double()
    match best {
      Some(b) =>
        if score > b.score {
          best = Some(BiotiteTrimMatch::{ trim_amount: o, score, })
        }
      None => best = Some(BiotiteTrimMatch::{ trim_amount: o, score, })
    }
  }
  best
}

///|
/// Trim primer regions from a single-strand read.
///
/// The primer is compared against the left and the right sequence end;
/// the better-scoring end defines the cut. Returns the retained
/// boundaries [start, stop). When the primer matches nowhere, the full
/// sequence is retained.
pub fn biotite_trim_single_strand(
  sequence : String,
  primer : String,
  min_overlap~ : Int32 = 0,
  max_mismatch~ : Int32 = 0,
  distance_weight~ : Double = 0.15,
) -> (Int32, Int32) {
  let n = sequence.length()
  let left = biotite_trim_best_end(
    sequence, primer, true, min_overlap, max_mismatch, distance_weight,
  )
  let right = biotite_trim_best_end(
    sequence, primer, false, min_overlap, max_mismatch, distance_weight,
  )
  match (left, right) {
    (Some(l), Some(r)) => {
      let ls = l.score
      let rs = r.score
      if ls >= rs {
        let start = l.trim_amount
        (start, n)
      } else {
        let stop = n - r.trim_amount
        (0, stop)
      }
    }
    (Some(l), None) => {
      let start = l.trim_amount
      (start, n)
    }
    (None, Some(r)) => {
      let stop = n - r.trim_amount
      (0, stop)
    }
    (None, None) => (0, n)
  }
}

///|
/// Trim two primers from a double-strand read.
///
/// Each primer is matched against both ends independently
/// (see `biotite_trim_single_strand`). Returns the two retained
/// boundary pairs.
pub fn biotite_trim_double_strand(
  sequence : String,
  primer_1 : String,
  primer_2 : String,
  min_overlap~ : Int32 = 0,
  max_mismatch~ : Int32 = 0,
  distance_weight~ : Double = 0.15,
) -> ((Int32, Int32), (Int32, Int32)) {
  let t1 = biotite_trim_single_strand(
    sequence,
    primer_1,
    min_overlap~,
    max_mismatch~,
    distance_weight~,
  )
  let t2 = biotite_trim_single_strand(
    sequence,
    primer_2,
    min_overlap~,
    max_mismatch~,
    distance_weight~,
  )
  (t1, t2)
}

///|
pub fn _biotite_sequence_trim_dummy() -> Unit {
  let _ = biotite_trim_single_strand("ACGTACGT", "ACGT")
  let _ = biotite_trim_double_strand("ACGTACGT", "ACGT", "ACGT")
  let _ = biotite_trim_best_end("ACGT", "AC", true, 0, 0, 0.15)
}
