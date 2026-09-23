// Detailed coordinate alignment statistics compatible with Biopython
// Bio.Align.Alignment.counts and AlignmentCounts.

///|
/// Error raised while configuring or calculating detailed alignment counts.
pub suberror AlignmentCountsError {
  AlignmentCountsError(String)
}

///|
/// Affine gap scores for each gap side and direction.
///
/// An insertion is a gap in the target row; a deletion is a gap in the query
/// row. A gap of length n contributes one open score and n - 1 extend scores.
pub struct AlignmentGapScores {
  open_left_insertion_score : Double
  extend_left_insertion_score : Double
  open_left_deletion_score : Double
  extend_left_deletion_score : Double
  open_internal_insertion_score : Double
  extend_internal_insertion_score : Double
  open_internal_deletion_score : Double
  extend_internal_deletion_score : Double
  open_right_insertion_score : Double
  extend_right_insertion_score : Double
  open_right_deletion_score : Double
  extend_right_deletion_score : Double
} derive(Eq, Debug)
pub extend AlignmentGapScores with Eq::{equal, not_equal}
pub extend AlignmentGapScores with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn alignment_counts_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn alignment_counts_validate_score(
  value : Double,
  label : String,
) -> Unit raise AlignmentCountsError {
  if !alignment_counts_is_finite(value) {
    raise AlignmentCountsError(label + " must be finite")
  }
}

///|
/// Construct fully position- and direction-specific affine gap scores.
pub fn AlignmentGapScores::create(
  open_left_insertion_score~ : Double,
  extend_left_insertion_score~ : Double,
  open_left_deletion_score~ : Double,
  extend_left_deletion_score~ : Double,
  open_internal_insertion_score~ : Double,
  extend_internal_insertion_score~ : Double,
  open_internal_deletion_score~ : Double,
  extend_internal_deletion_score~ : Double,
  open_right_insertion_score~ : Double,
  extend_right_insertion_score~ : Double,
  open_right_deletion_score~ : Double,
  extend_right_deletion_score~ : Double,
) -> AlignmentGapScores raise AlignmentCountsError {
  let values = [
    open_left_insertion_score, extend_left_insertion_score, open_left_deletion_score,
    extend_left_deletion_score, open_internal_insertion_score, extend_internal_insertion_score,
    open_internal_deletion_score, extend_internal_deletion_score, open_right_insertion_score,
    extend_right_insertion_score, open_right_deletion_score, extend_right_deletion_score,
  ]
  let labels = [
    "open left insertion score", "extend left insertion score", "open left deletion score",
    "extend left deletion score", "open internal insertion score", "extend internal insertion score",
    "open internal deletion score", "extend internal deletion score", "open right insertion score",
    "extend right insertion score", "open right deletion score", "extend right deletion score",
  ]
  let mut index = 0
  while index < values.length() {
    alignment_counts_validate_score(values[index], labels[index])
    index = index + 1
  }
  AlignmentGapScores::{
    open_left_insertion_score,
    extend_left_insertion_score,
    open_left_deletion_score,
    extend_left_deletion_score,
    open_internal_insertion_score,
    extend_internal_insertion_score,
    open_internal_deletion_score,
    extend_internal_deletion_score,
    open_right_insertion_score,
    extend_right_insertion_score,
    open_right_deletion_score,
    extend_right_deletion_score,
  }
}

///|
/// Use one open score and one extend score for every gap class.
pub fn AlignmentGapScores::affine(
  open_score : Double,
  extend_score : Double,
) -> AlignmentGapScores raise AlignmentCountsError {
  alignment_counts_validate_score(open_score, "gap open score")
  alignment_counts_validate_score(extend_score, "gap extend score")
  AlignmentGapScores::{
    open_left_insertion_score: open_score,
    extend_left_insertion_score: extend_score,
    open_left_deletion_score: open_score,
    extend_left_deletion_score: extend_score,
    open_internal_insertion_score: open_score,
    extend_internal_insertion_score: extend_score,
    open_internal_deletion_score: open_score,
    extend_internal_deletion_score: extend_score,
    open_right_insertion_score: open_score,
    extend_right_insertion_score: extend_score,
    open_right_deletion_score: open_score,
    extend_right_deletion_score: extend_score,
  }
}

///|
/// Use direction-specific affine scores independent of gap position.
pub fn AlignmentGapScores::directional(
  insertion_open_score~ : Double,
  insertion_extend_score~ : Double,
  deletion_open_score~ : Double,
  deletion_extend_score~ : Double,
) -> AlignmentGapScores raise AlignmentCountsError {
  alignment_counts_validate_score(insertion_open_score, "insertion open score")
  alignment_counts_validate_score(
    insertion_extend_score, "insertion extend score",
  )
  alignment_counts_validate_score(deletion_open_score, "deletion open score")
  alignment_counts_validate_score(
    deletion_extend_score, "deletion extend score",
  )
  AlignmentGapScores::{
    open_left_insertion_score: insertion_open_score,
    extend_left_insertion_score: insertion_extend_score,
    open_left_deletion_score: deletion_open_score,
    extend_left_deletion_score: deletion_extend_score,
    open_internal_insertion_score: insertion_open_score,
    extend_internal_insertion_score: insertion_extend_score,
    open_internal_deletion_score: deletion_open_score,
    extend_internal_deletion_score: deletion_extend_score,
    open_right_insertion_score: insertion_open_score,
    extend_right_insertion_score: insertion_extend_score,
    open_right_deletion_score: deletion_open_score,
    extend_right_deletion_score: deletion_extend_score,
  }
}

///|
/// Optional composition and gap scoring used by `alignment_counts`.
///
/// A substitution matrix takes precedence over match/mismatch scores.
/// Substitution and gap scores remain unknown unless their corresponding
/// configuration is present.
pub struct AlignmentCountsConfig {
  wildcard : Char?
  substitution_matrix : SubstitutionMatrix?
  match_score : Double?
  mismatch_score : Double?
  gap_scores : AlignmentGapScores?
}

///|
pub fn AlignmentCountsConfig::create(
  wildcard? : Char? = None,
  substitution_matrix? : SubstitutionMatrix? = None,
  match_score? : Double? = None,
  mismatch_score? : Double? = None,
  gap_scores? : AlignmentGapScores? = None,
) -> AlignmentCountsConfig raise AlignmentCountsError {
  match (match_score, mismatch_score) {
    (Some(match_value), Some(mismatch_value)) => {
      alignment_counts_validate_score(match_value, "match score")
      alignment_counts_validate_score(mismatch_value, "mismatch score")
    }
    (None, None) => ()
    _ =>
      raise AlignmentCountsError(
        "match and mismatch scores must be provided together",
      )
  }
  AlignmentCountsConfig::{
    wildcard,
    substitution_matrix,
    match_score,
    mismatch_score,
    gap_scores,
  }
}

///|
pub fn AlignmentCountsConfig::default() -> AlignmentCountsConfig {
  AlignmentCountsConfig::{
    wildcard: None,
    substitution_matrix: None,
    match_score: None,
    mismatch_score: None,
    gap_scores: None,
  }
}

///|
/// Configure wildcard-aware counting without calculating scores.
pub fn AlignmentCountsConfig::with_wildcard(
  wildcard : Char,
) -> AlignmentCountsConfig {
  AlignmentCountsConfig::{
    wildcard: Some(wildcard),
    substitution_matrix: None,
    match_score: None,
    mismatch_score: None,
    gap_scores: None,
  }
}

///|
/// Configure substitution scoring and positive-match counting.
pub fn AlignmentCountsConfig::with_matrix(
  matrix : SubstitutionMatrix,
  wildcard? : Char? = None,
) -> AlignmentCountsConfig {
  AlignmentCountsConfig::{
    wildcard,
    substitution_matrix: Some(matrix),
    match_score: None,
    mismatch_score: None,
    gap_scores: None,
  }
}

///|
/// Configure complete match/mismatch and affine gap scoring.
pub fn AlignmentCountsConfig::scored(
  match_score~ : Double,
  mismatch_score~ : Double,
  gap_scores~ : AlignmentGapScores,
  wildcard? : Char? = None,
) -> AlignmentCountsConfig raise AlignmentCountsError {
  alignment_counts_validate_score(match_score, "match score")
  alignment_counts_validate_score(mismatch_score, "mismatch score")
  AlignmentCountsConfig::{
    wildcard,
    substitution_matrix: None,
    match_score: Some(match_score),
    mismatch_score: Some(mismatch_score),
    gap_scores: Some(gap_scores),
  }
}

///|
/// Configure complete substitution-matrix and affine gap scoring.
pub fn AlignmentCountsConfig::scored_matrix(
  matrix~ : SubstitutionMatrix,
  gap_scores~ : AlignmentGapScores,
  wildcard? : Char? = None,
) -> AlignmentCountsConfig {
  AlignmentCountsConfig::{
    wildcard,
    substitution_matrix: Some(matrix),
    match_score: None,
    mismatch_score: None,
    gap_scores: Some(gap_scores),
  }
}

///|
/// Detailed counts and optional scores for a coordinate alignment.
pub struct AlignmentCounts {
  open_left_insertions : Int
  extend_left_insertions : Int
  open_left_deletions : Int
  extend_left_deletions : Int
  open_internal_insertions : Int
  extend_internal_insertions : Int
  open_internal_deletions : Int
  extend_internal_deletions : Int
  open_right_insertions : Int
  extend_right_insertions : Int
  open_right_deletions : Int
  extend_right_deletions : Int
  aligned : Int
  identities : Int
  mismatches : Int
  positives : Int?
  gap_score : Double?
  substitution_score : Double?
  score : Double?
} derive(Debug)
pub extend AlignmentCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn AlignmentCounts::left_insertions(self : AlignmentCounts) -> Int {
  self.open_left_insertions + self.extend_left_insertions
}

///|
pub fn AlignmentCounts::internal_insertions(self : AlignmentCounts) -> Int {
  self.open_internal_insertions + self.extend_internal_insertions
}

///|
pub fn AlignmentCounts::right_insertions(self : AlignmentCounts) -> Int {
  self.open_right_insertions + self.extend_right_insertions
}

///|
pub fn AlignmentCounts::left_deletions(self : AlignmentCounts) -> Int {
  self.open_left_deletions + self.extend_left_deletions
}

///|
pub fn AlignmentCounts::internal_deletions(self : AlignmentCounts) -> Int {
  self.open_internal_deletions + self.extend_internal_deletions
}

///|
pub fn AlignmentCounts::right_deletions(self : AlignmentCounts) -> Int {
  self.open_right_deletions + self.extend_right_deletions
}

///|
pub fn AlignmentCounts::open_insertions(self : AlignmentCounts) -> Int {
  self.open_left_insertions +
  self.open_internal_insertions +
  self.open_right_insertions
}

///|
pub fn AlignmentCounts::extend_insertions(self : AlignmentCounts) -> Int {
  self.extend_left_insertions +
  self.extend_internal_insertions +
  self.extend_right_insertions
}

///|
pub fn AlignmentCounts::insertions(self : AlignmentCounts) -> Int {
  self.open_insertions() + self.extend_insertions()
}

///|
pub fn AlignmentCounts::open_deletions(self : AlignmentCounts) -> Int {
  self.open_left_deletions +
  self.open_internal_deletions +
  self.open_right_deletions
}

///|
pub fn AlignmentCounts::extend_deletions(self : AlignmentCounts) -> Int {
  self.extend_left_deletions +
  self.extend_internal_deletions +
  self.extend_right_deletions
}

///|
pub fn AlignmentCounts::deletions(self : AlignmentCounts) -> Int {
  self.open_deletions() + self.extend_deletions()
}

///|
pub fn AlignmentCounts::open_left_gaps(self : AlignmentCounts) -> Int {
  self.open_left_insertions + self.open_left_deletions
}

///|
pub fn AlignmentCounts::extend_left_gaps(self : AlignmentCounts) -> Int {
  self.extend_left_insertions + self.extend_left_deletions
}

///|
pub fn AlignmentCounts::left_gaps(self : AlignmentCounts) -> Int {
  self.open_left_gaps() + self.extend_left_gaps()
}

///|
pub fn AlignmentCounts::open_internal_gaps(self : AlignmentCounts) -> Int {
  self.open_internal_insertions + self.open_internal_deletions
}

///|
pub fn AlignmentCounts::extend_internal_gaps(self : AlignmentCounts) -> Int {
  self.extend_internal_insertions + self.extend_internal_deletions
}

///|
pub fn AlignmentCounts::internal_gaps(self : AlignmentCounts) -> Int {
  self.open_internal_gaps() + self.extend_internal_gaps()
}

///|
pub fn AlignmentCounts::open_right_gaps(self : AlignmentCounts) -> Int {
  self.open_right_insertions + self.open_right_deletions
}

///|
pub fn AlignmentCounts::extend_right_gaps(self : AlignmentCounts) -> Int {
  self.extend_right_insertions + self.extend_right_deletions
}

///|
pub fn AlignmentCounts::right_gaps(self : AlignmentCounts) -> Int {
  self.open_right_gaps() + self.extend_right_gaps()
}

///|
pub fn AlignmentCounts::open_gaps(self : AlignmentCounts) -> Int {
  self.open_insertions() + self.open_deletions()
}

///|
pub fn AlignmentCounts::extend_gaps(self : AlignmentCounts) -> Int {
  self.extend_insertions() + self.extend_deletions()
}

///|
pub fn AlignmentCounts::gaps(self : AlignmentCounts) -> Int {
  self.open_gaps() + self.extend_gaps()
}

///|
pub fn AlignmentCounts::summary(self : AlignmentCounts) -> String {
  "AlignmentCounts(aligned=" +
  self.aligned.to_string() +
  ", identities=" +
  self.identities.to_string() +
  ", mismatches=" +
  self.mismatches.to_string() +
  ", gaps=" +
  self.gaps().to_string() +
  ")"
}

///|
struct AlignmentCountsAccumulator {
  mut open_left_insertions : Int
  mut extend_left_insertions : Int
  mut open_left_deletions : Int
  mut extend_left_deletions : Int
  mut open_internal_insertions : Int
  mut extend_internal_insertions : Int
  mut open_internal_deletions : Int
  mut extend_internal_deletions : Int
  mut open_right_insertions : Int
  mut extend_right_insertions : Int
  mut open_right_deletions : Int
  mut extend_right_deletions : Int
  mut aligned : Int
  mut identities : Int
  mut mismatches : Int
  mut positives : Int
  mut matrix_score : Double
}

///|
fn alignment_counts_accumulator() -> AlignmentCountsAccumulator {
  AlignmentCountsAccumulator::{
    open_left_insertions: 0,
    extend_left_insertions: 0,
    open_left_deletions: 0,
    extend_left_deletions: 0,
    open_internal_insertions: 0,
    extend_internal_insertions: 0,
    open_internal_deletions: 0,
    extend_internal_deletions: 0,
    open_right_insertions: 0,
    extend_right_insertions: 0,
    open_right_deletions: 0,
    extend_right_deletions: 0,
    aligned: 0,
    identities: 0,
    mismatches: 0,
    positives: 0,
    matrix_score: 0.0,
  }
}

///|
fn alignment_counts_normalize_coordinates(
  coordinates : Array[Int],
  sequence_length : Int,
  reverse : Bool,
) -> Array[Int] {
  let normalized : Array[Int] = []
  for coordinate in coordinates {
    normalized.push(
      if reverse {
        sequence_length - coordinate
      } else {
        coordinate
      },
    )
  }
  normalized
}

///|
fn alignment_counts_complement(value : Char) -> Char {
  match value {
    'A' => 'T'
    'a' => 't'
    'C' => 'G'
    'c' => 'g'
    'G' => 'C'
    'g' => 'c'
    'T' => 'A'
    't' => 'a'
    'U' => 'A'
    'u' => 'a'
    _ => value
  }
}

///|
fn alignment_counts_residue(
  sequence : String,
  sequence_length : Int,
  position : Int,
  reverse : Bool,
) -> Char {
  if reverse {
    alignment_counts_complement(
      sequence.unsafe_get(sequence_length - position - 1).unsafe_to_char(),
    )
  } else {
    sequence.unsafe_get(position).unsafe_to_char()
  }
}

///|
fn alignment_counts_add_gap(
  counts : AlignmentCountsAccumulator,
  insertion : Bool,
  continuation : Bool,
  start : Int,
  end : Int,
  left : Int,
  right : Int,
  size : Int,
) -> Unit {
  let side = if start == left { 0 } else if end == right { 2 } else { 1 }
  let open = if continuation { 0 } else { 1 }
  let extend = if continuation { size } else { size - 1 }
  if insertion {
    if side == 0 {
      counts.open_left_insertions = counts.open_left_insertions + open
      counts.extend_left_insertions = counts.extend_left_insertions + extend
    } else if side == 1 {
      counts.open_internal_insertions = counts.open_internal_insertions + open
      counts.extend_internal_insertions = counts.extend_internal_insertions +
        extend
    } else {
      counts.open_right_insertions = counts.open_right_insertions + open
      counts.extend_right_insertions = counts.extend_right_insertions + extend
    }
  } else if side == 0 {
    counts.open_left_deletions = counts.open_left_deletions + open
    counts.extend_left_deletions = counts.extend_left_deletions + extend
  } else if side == 1 {
    counts.open_internal_deletions = counts.open_internal_deletions + open
    counts.extend_internal_deletions = counts.extend_internal_deletions + extend
  } else {
    counts.open_right_deletions = counts.open_right_deletions + open
    counts.extend_right_deletions = counts.extend_right_deletions + extend
  }
}

///|
fn alignment_counts_add_residues(
  counts : AlignmentCountsAccumulator,
  first : Char,
  second : Char,
  config : AlignmentCountsConfig,
) -> Unit raise AlignmentCountsError {
  match config.wildcard {
    Some(wildcard) => if first == wildcard || second == wildcard { return }
    None => ()
  }
  if first == second {
    counts.identities = counts.identities + 1
  } else {
    counts.mismatches = counts.mismatches + 1
  }
  match config.substitution_matrix {
    Some(matrix) => {
      let first_text = first.to_string()
      let second_text = second.to_string()
      if !matrix.contains(first_text) {
        raise AlignmentCountsError(
          "target residue '" + first_text + "' is not in substitution matrix",
        )
      }
      if !matrix.contains(second_text) {
        raise AlignmentCountsError(
          "query residue '" + second_text + "' is not in substitution matrix",
        )
      }
      let value = matrix.get_score_case_insensitive(first_text, second_text)
      counts.matrix_score = counts.matrix_score + value.to_double()
      if value > 0 {
        counts.positives = counts.positives + 1
      }
    }
    None => ()
  }
}

///|
fn alignment_counts_add_pair(
  counts : AlignmentCountsAccumulator,
  first_sequence : String,
  first_length : Int,
  first_coordinates : Array[Int],
  first_reverse : Bool,
  second_sequence : String,
  second_length : Int,
  second_coordinates : Array[Int],
  second_reverse : Bool,
  config : AlignmentCountsConfig,
) -> Unit raise AlignmentCountsError {
  if first_coordinates.length() != second_coordinates.length() {
    raise AlignmentCountsError(
      "alignment coordinate rows must have equal length",
    )
  }
  if first_coordinates.length() == 1 {
    raise AlignmentCountsError(
      "an alignment path must be empty or contain at least two points",
    )
  }
  let first = alignment_counts_normalize_coordinates(
    first_coordinates, first_length, first_reverse,
  )
  let second = alignment_counts_normalize_coordinates(
    second_coordinates, second_length, second_reverse,
  )
  if first.length() == 0 {
    return
  }
  let first_left = first[0]
  let first_right = first[first.length() - 1]
  let second_left = second[0]
  let second_right = second[second.length() - 1]
  let mut path = 0
  let mut index = 0
  while index + 1 < first.length() {
    let first_start = first[index]
    let first_end = first[index + 1]
    let second_start = second[index]
    let second_end = second[index + 1]
    let first_step = first_end - first_start
    let second_step = second_end - second_start
    if first_step < 0 || second_step < 0 {
      raise AlignmentCountsError(
        "normalized alignment coordinates must be non-decreasing",
      )
    }
    if first_step == 0 && second_step == 0 {
      // Multiple-alignment row pairs may share a gap inserted for other rows.
    } else if first_step == 0 {
      alignment_counts_add_gap(
        counts,
        true,
        path == 1,
        first_start,
        first_end,
        first_left,
        first_right,
        second_step,
      )
      path = 1
    } else if second_step == 0 {
      alignment_counts_add_gap(
        counts,
        false,
        path == 2,
        second_start,
        second_end,
        second_left,
        second_right,
        first_step,
      )
      path = 2
    } else {
      if first_step != second_step {
        raise AlignmentCountsError(
          "aligned coordinate steps must have equal length",
        )
      }
      path = 0
      counts.aligned = counts.aligned + first_step
      if first_sequence.length() > 0 && second_sequence.length() > 0 {
        let mut offset = 0
        while offset < first_step {
          let first_residue = alignment_counts_residue(
            first_sequence,
            first_length,
            first_start + offset,
            first_reverse,
          )
          let second_residue = alignment_counts_residue(
            second_sequence,
            second_length,
            second_start + offset,
            second_reverse,
          )
          alignment_counts_add_residues(
            counts, first_residue, second_residue, config,
          )
          offset = offset + 1
        }
      }
    }
    index = index + 1
  }
}

///|
fn alignment_counts_gap_score(
  counts : AlignmentCountsAccumulator,
  scores : AlignmentGapScores,
) -> Double {
  counts.open_left_insertions.to_double() * scores.open_left_insertion_score +
  counts.extend_left_insertions.to_double() * scores.extend_left_insertion_score +
  counts.open_left_deletions.to_double() * scores.open_left_deletion_score +
  counts.extend_left_deletions.to_double() * scores.extend_left_deletion_score +
  counts.open_internal_insertions.to_double() *
  scores.open_internal_insertion_score +
  counts.extend_internal_insertions.to_double() *
  scores.extend_internal_insertion_score +
  counts.open_internal_deletions.to_double() *
  scores.open_internal_deletion_score +
  counts.extend_internal_deletions.to_double() *
  scores.extend_internal_deletion_score +
  counts.open_right_insertions.to_double() * scores.open_right_insertion_score +
  counts.extend_right_insertions.to_double() *
  scores.extend_right_insertion_score +
  counts.open_right_deletions.to_double() * scores.open_right_deletion_score +
  counts.extend_right_deletions.to_double() * scores.extend_right_deletion_score
}

///|
fn alignment_counts_finish(
  counts : AlignmentCountsAccumulator,
  config : AlignmentCountsConfig,
) -> AlignmentCounts {
  let positives = match config.substitution_matrix {
    Some(_) => Some(counts.positives)
    None => None
  }
  let substitution_score = match config.substitution_matrix {
    Some(_) => Some(counts.matrix_score)
    None =>
      match (config.match_score, config.mismatch_score) {
        (Some(match_value), Some(mismatch_value)) =>
          Some(
            counts.identities.to_double() * match_value +
            counts.mismatches.to_double() * mismatch_value,
          )
        _ => None
      }
  }
  let gap_score = match config.gap_scores {
    Some(scores) => Some(alignment_counts_gap_score(counts, scores))
    None => None
  }
  let score = match (substitution_score, gap_score) {
    (Some(substitution_value), Some(gap_value)) =>
      Some(substitution_value + gap_value)
    _ => None
  }
  AlignmentCounts::{
    open_left_insertions: counts.open_left_insertions,
    extend_left_insertions: counts.extend_left_insertions,
    open_left_deletions: counts.open_left_deletions,
    extend_left_deletions: counts.extend_left_deletions,
    open_internal_insertions: counts.open_internal_insertions,
    extend_internal_insertions: counts.extend_internal_insertions,
    open_internal_deletions: counts.open_internal_deletions,
    extend_internal_deletions: counts.extend_internal_deletions,
    open_right_insertions: counts.open_right_insertions,
    extend_right_insertions: counts.extend_right_insertions,
    open_right_deletions: counts.open_right_deletions,
    extend_right_deletions: counts.extend_right_deletions,
    aligned: counts.aligned,
    identities: counts.identities,
    mismatches: counts.mismatches,
    positives,
    gap_score,
    substitution_score,
    score,
  }
}

///|
/// Calculate Biopython-style detailed statistics for a pairwise alignment.
pub fn CoordinatePairwiseAlignment::alignment_counts(
  self : CoordinatePairwiseAlignment,
  config? : AlignmentCountsConfig = AlignmentCountsConfig::default(),
) -> AlignmentCounts raise AlignmentCountsError {
  let counts = alignment_counts_accumulator()
  alignment_counts_add_pair(
    counts,
    self.target_sequence,
    self.target_length,
    self.target_coordinates,
    false,
    self.query_sequence,
    self.query_length,
    self.query_coordinates,
    self.is_reverse(),
    config,
  )
  alignment_counts_finish(counts, config)
}

///|
/// Sum detailed statistics over every unordered pair of MSA rows.
pub fn CoordinateMultipleAlignment::alignment_counts(
  self : CoordinateMultipleAlignment,
  config? : AlignmentCountsConfig = AlignmentCountsConfig::default(),
) -> AlignmentCounts raise AlignmentCountsError {
  if self.names.length() != self.sequences.length() ||
    self.names.length() != self.coordinates.length() {
    raise AlignmentCountsError(
      "MSA names, sequences, and coordinates must have equal row counts",
    )
  }
  let counts = alignment_counts_accumulator()
  let mut first = 0
  while first < self.sequences.length() {
    let mut second = first + 1
    while second < self.sequences.length() {
      alignment_counts_add_pair(
        counts,
        self.sequences[first],
        self.sequences[first].length(),
        self.coordinates[first],
        false,
        self.sequences[second],
        self.sequences[second].length(),
        self.coordinates[second],
        false,
        config,
      )
      second = second + 1
    }
    first = first + 1
  }
  alignment_counts_finish(counts, config)
}
