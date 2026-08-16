///|
/// RangedSummarizedExperiment - Genomic-range-aware assay container inspired
/// by Bioconductor SummarizedExperiment.

///|
pub suberror RangedSummarizedExperimentError {
  RangedSummarizedExperimentError(String)
}

///|
/// A SummarizedExperiment whose rows are represented by genomic ranges.
///
/// `row_data` and `row_names` are optional parallel annotations. Empty arrays
/// represent absent annotations; otherwise their length must equal `nrow`.
pub struct RangedSummarizedExperiment {
  experiment : SummarizedExperiment
  row_ranges : GRanges
  row_range_groups : GRangesList?
  row_data : Array[Map[String, String]]
  row_names : Array[String]
}

///|
/// Construct a range-aware experiment and validate all parallel dimensions.
pub fn RangedSummarizedExperiment::new(
  assays~ : Map[String, Array[Array[Double]]],
  row_ranges~ : GRanges,
  col_data~ : Array[Map[String, String]],
  row_data? : Array[Map[String, String]] = [],
  row_names? : Array[String] = [],
  metadata? : Map[String, String] = Map([]),
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  match
    rse_validation_error(assays, row_ranges, col_data, row_data, row_names) {
    Some(message) => raise RangedSummarizedExperimentError(message)
    None => {
      let experiment = summarized_experiment(
        assays,
        rse_ranges_as_tuples(row_ranges),
        col_data,
        metadata,
      )
      RangedSummarizedExperiment::{
        experiment,
        row_ranges,
        row_range_groups: None,
        row_data,
        row_names,
      }
    }
  }
}

///|
/// Construct an experiment whose rows are compound genomic features.
///
/// Each GRangesList element corresponds to one assay row, such as a transcript
/// represented by multiple exons. Feature-level GRangesList metadata is used as
/// row data when explicit row data is absent.
pub fn RangedSummarizedExperiment::new_with_range_groups(
  assays~ : Map[String, Array[Array[Double]]],
  row_range_groups~ : GRangesList,
  col_data~ : Array[Map[String, String]],
  row_data? : Array[Map[String, String]] = [],
  row_names? : Array[String] = [],
  metadata? : Map[String, String] = Map([]),
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  let effective_row_data = if row_data.length() == 0 &&
    row_range_groups.element_metadata().length() > 0 {
    row_range_groups.element_metadata()
  } else {
    row_data
  }
  let bounds = row_range_groups.feature_bounds()
  match
    rse_grouped_validation_error(
      assays, row_range_groups, col_data, effective_row_data, row_names,
    ) {
    Some(message) => raise RangedSummarizedExperimentError(message)
    None => {
      let experiment = summarized_experiment(
        assays,
        rse_ranges_as_tuples(bounds),
        col_data,
        metadata,
      )
      RangedSummarizedExperiment::{
        experiment,
        row_ranges: bounds,
        row_range_groups: Some(row_range_groups),
        row_data: effective_row_data,
        row_names,
      }
    }
  }
}

///|
/// Attach full genomic ranges to an existing SummarizedExperiment.
pub fn RangedSummarizedExperiment::from_experiment(
  experiment~ : SummarizedExperiment,
  row_ranges~ : GRanges,
  row_data? : Array[Map[String, String]] = [],
  row_names? : Array[String] = [],
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  RangedSummarizedExperiment::new(
    assays=rse_copy_assays(experiment),
    row_ranges~,
    col_data=se_col_data(experiment),
    row_data~,
    row_names~,
    metadata=se_metadata(experiment),
  )
}

///|
/// Attach grouped genomic features to an existing SummarizedExperiment.
pub fn RangedSummarizedExperiment::from_experiment_with_range_groups(
  experiment~ : SummarizedExperiment,
  row_range_groups~ : GRangesList,
  row_data? : Array[Map[String, String]] = [],
  row_names? : Array[String] = [],
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  RangedSummarizedExperiment::new_with_range_groups(
    assays=rse_copy_assays(experiment),
    row_range_groups~,
    col_data=se_col_data(experiment),
    row_data~,
    row_names~,
    metadata=se_metadata(experiment),
  )
}

///|
/// Return the underlying SummarizedExperiment.
pub fn RangedSummarizedExperiment::experiment(
  self : RangedSummarizedExperiment,
) -> SummarizedExperiment {
  self.experiment
}

///|
/// Return an assay by name.
pub fn RangedSummarizedExperiment::assay(
  self : RangedSummarizedExperiment,
  name : String,
) -> Array[Array[Double]]? {
  se_assay(self.experiment, name)
}

///|
/// Return all assay names.
pub fn RangedSummarizedExperiment::assay_names(
  self : RangedSummarizedExperiment,
) -> Array[String] {
  se_assay_names(self.experiment)
}

///|
/// Return one genomic range associated with each assay row.
///
/// For grouped rows this returns feature bounds. Use `row_range_groups` to
/// retrieve the exact member ranges.
pub fn RangedSummarizedExperiment::row_ranges(
  self : RangedSummarizedExperiment,
) -> GRanges {
  self.row_ranges
}

///|
/// Return exact grouped row ranges when this is a GRangesList-backed object.
pub fn RangedSummarizedExperiment::row_range_groups(
  self : RangedSummarizedExperiment,
) -> GRangesList? {
  match self.row_range_groups {
    Some(groups) =>
      Some(
        grl_create_unchecked(
          groups.elements(),
          groups.names(),
          groups.element_metadata(),
        ),
      )
    None => None
  }
}

///|
/// Return whether assay rows are represented by compound genomic features.
pub fn RangedSummarizedExperiment::has_grouped_row_ranges(
  self : RangedSummarizedExperiment,
) -> Bool {
  self.row_range_groups is Some(_)
}

///|
/// Return optional per-row annotations.
pub fn RangedSummarizedExperiment::row_data(
  self : RangedSummarizedExperiment,
) -> Array[Map[String, String]] {
  self.row_data
}

///|
/// Return optional row names.
pub fn RangedSummarizedExperiment::row_names(
  self : RangedSummarizedExperiment,
) -> Array[String] {
  self.row_names
}

///|
/// Return per-column sample annotations.
pub fn RangedSummarizedExperiment::col_data(
  self : RangedSummarizedExperiment,
) -> Array[Map[String, String]] {
  se_col_data(self.experiment)
}

///|
/// Return experiment-wide metadata.
pub fn RangedSummarizedExperiment::metadata(
  self : RangedSummarizedExperiment,
) -> Map[String, String] {
  se_metadata(self.experiment)
}

///|
/// Return the number of genomic features.
pub fn RangedSummarizedExperiment::nrow(
  self : RangedSummarizedExperiment,
) -> Int {
  match self.row_range_groups {
    Some(groups) => groups.length()
    None => granges_length(self.row_ranges)
  }
}

///|
/// Return the number of samples.
pub fn RangedSummarizedExperiment::ncol(
  self : RangedSummarizedExperiment,
) -> Int {
  se_ncol(self.experiment)
}

///|
/// Check ranges, annotations, and assay dimensions.
pub fn RangedSummarizedExperiment::is_valid(
  self : RangedSummarizedExperiment,
) -> Bool {
  match self.row_range_groups {
    Some(groups) =>
      rse_grouped_validation_error(
        rse_copy_assays(self.experiment),
        groups,
        se_col_data(self.experiment),
        self.row_data,
        self.row_names,
      )
      is None &&
      rse_ranges_equal(self.row_ranges, groups.feature_bounds())
    None =>
      rse_validation_error(
        rse_copy_assays(self.experiment),
        self.row_ranges,
        se_col_data(self.experiment),
        self.row_data,
        self.row_names,
      )
      is None
  }
}

///|
/// Replace row ranges without changing assays or annotations.
pub fn RangedSummarizedExperiment::with_row_ranges(
  self : RangedSummarizedExperiment,
  row_ranges : GRanges,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  RangedSummarizedExperiment::new(
    assays=rse_copy_assays(self.experiment),
    row_ranges~,
    col_data=se_col_data(self.experiment),
    row_data=self.row_data,
    row_names=self.row_names,
    metadata=se_metadata(self.experiment),
  )
}

///|
/// Replace row ranges with compound features without changing assays.
pub fn RangedSummarizedExperiment::with_row_range_groups(
  self : RangedSummarizedExperiment,
  row_range_groups : GRangesList,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  RangedSummarizedExperiment::new_with_range_groups(
    assays=rse_copy_assays(self.experiment),
    row_range_groups~,
    col_data=se_col_data(self.experiment),
    row_data=self.row_data,
    row_names=self.row_names,
    metadata=se_metadata(self.experiment),
  )
}

///|
/// Select, reorder, or duplicate rows while keeping all row data coordinated.
pub fn RangedSummarizedExperiment::subset_rows(
  self : RangedSummarizedExperiment,
  indices : Array[Int],
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  rse_validate_indices(indices, self.nrow(), "row")
  let row_range_groups = match self.row_range_groups {
    Some(groups) => Some(rse_take_range_groups(groups, indices))
    None => None
  }
  RangedSummarizedExperiment::{
    experiment: se_subset_rows(self.experiment, indices),
    row_ranges: rse_take_ranges(self.row_ranges, indices),
    row_range_groups,
    row_data: rse_take_maps(self.row_data, indices),
    row_names: rse_take_strings(self.row_names, indices),
  }
}

///|
/// Select, reorder, or duplicate columns while keeping sample data coordinated.
pub fn RangedSummarizedExperiment::subset_cols(
  self : RangedSummarizedExperiment,
  indices : Array[Int],
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  rse_validate_indices(indices, self.ncol(), "column")
  RangedSummarizedExperiment::{
    experiment: se_subset_cols(self.experiment, indices),
    row_ranges: self.row_ranges,
    row_range_groups: self.row_range_groups,
    row_data: self.row_data,
    row_names: self.row_names,
  }
}

///|
/// Find all row-range overlaps with a subject GRanges.
///
/// Returned pairs contain zero-based `(row_index, subject_index)` values.
pub fn RangedSummarizedExperiment::find_overlaps(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[(Int, Int)] {
  match self.row_range_groups {
    Some(groups) => return groups.find_overlaps(subject, ignore_strand~)
    None => ()
  }
  let hits : Array[(Int, Int)] = []
  for query_index = 0; query_index < self.nrow(); query_index = query_index + 1 {
    for subject_index = 0
        subject_index < granges_length(subject)
        subject_index = subject_index + 1 {
      if rse_ranges_overlap(
          self.row_ranges,
          query_index,
          subject,
          subject_index,
          ignore_strand,
        ) {
        hits.push((query_index, subject_index))
      }
    }
  }
  hits
}

///|
/// Count subject overlaps for each assay row.
pub fn RangedSummarizedExperiment::count_overlaps(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  let counts = Array::make(self.nrow(), 0)
  for hit in self.find_overlaps(subject, ignore_strand~) {
    counts[hit.0] = counts[hit.0] + 1
  }
  counts
}

///|
/// Return an overlap flag for each assay row.
pub fn RangedSummarizedExperiment::overlaps_any(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Bool] {
  self.count_overlaps(subject, ignore_strand~).map(fn(count) { count > 0 })
}

///|
/// Keep rows that overlap at least one subject range.
pub fn RangedSummarizedExperiment::subset_by_overlaps(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  let indices : Array[Int] = []
  let keep = self.overlaps_any(subject, ignore_strand~)
  for index = 0; index < keep.length(); index = index + 1 {
    if keep[index] {
      indices.push(index)
    }
  }
  self.subset_rows(indices)
}

///|
/// Return the nearest subject index for each row, or `-1` when unavailable.
pub fn RangedSummarizedExperiment::nearest(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  match self.row_range_groups {
    Some(groups) => return groups.nearest(subject, ignore_strand~)
    None => ()
  }
  let nearest_indices = Array::make(self.nrow(), -1)
  for query_index = 0; query_index < self.nrow(); query_index = query_index + 1 {
    let mut best_index = -1
    let mut best_distance = 0
    for subject_index = 0
        subject_index < granges_length(subject)
        subject_index = subject_index + 1 {
      if self.row_ranges.seqnames[query_index] ==
        subject.seqnames[subject_index] &&
        rse_strands_compatible(
          self.row_ranges.strands[query_index],
          subject.strands[subject_index],
          ignore_strand,
        ) {
        let distance = rse_interval_distance(
          self.row_ranges.starts[query_index],
          self.row_ranges.ends[query_index],
          subject.starts[subject_index],
          subject.ends[subject_index],
        )
        if best_index < 0 || distance < best_distance {
          best_index = subject_index
          best_distance = distance
        }
      }
    }
    nearest_indices[query_index] = best_index
  }
  nearest_indices
}

///|
/// Return the nearest subject distance for each row, or `-1` when unavailable.
pub fn RangedSummarizedExperiment::distance_to_nearest(
  self : RangedSummarizedExperiment,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  match self.row_range_groups {
    Some(groups) => return groups.distance_to_nearest(subject, ignore_strand~)
    None => ()
  }
  let distances = Array::make(self.nrow(), -1)
  let nearest_indices = self.nearest(subject, ignore_strand~)
  for query_index = 0; query_index < self.nrow(); query_index = query_index + 1 {
    let subject_index = nearest_indices[query_index]
    if subject_index >= 0 {
      distances[query_index] = rse_interval_distance(
        self.row_ranges.starts[query_index],
        self.row_ranges.ends[query_index],
        subject.starts[subject_index],
        subject.ends[subject_index],
      )
    }
  }
  distances
}

///|
/// Compute per-base coverage for configured sequence lengths.
pub fn RangedSummarizedExperiment::coverage(
  self : RangedSummarizedExperiment,
  sequence_lengths : Map[String, Int],
) -> Map[String, Array[Int]] {
  match self.row_range_groups {
    Some(groups) => groups.coverage(sequence_lengths)
    None => coverage_granges(self.row_ranges, sequence_lengths)
  }
}

///|
/// Shift every genomic row range by a fixed number of bases.
pub fn RangedSummarizedExperiment::shift(
  self : RangedSummarizedExperiment,
  amount : Int,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  match self.row_range_groups {
    Some(groups) => self.with_row_range_groups(groups.shift(amount))
    None => self.with_row_ranges(granges_shift(self.row_ranges, amount))
  }
}

///|
/// Narrow every genomic row range using one-based relative coordinates.
pub fn RangedSummarizedExperiment::narrow(
  self : RangedSummarizedExperiment,
  start : Int,
  end : Int,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  match self.row_range_groups {
    Some(groups) => self.with_row_range_groups(groups.narrow(start, end))
    None => self.with_row_ranges(granges_narrow(self.row_ranges, start, end))
  }
}

///|
/// Resize every genomic row range while fixing its start, end, or center.
pub fn RangedSummarizedExperiment::resize(
  self : RangedSummarizedExperiment,
  width : Int,
  fix : String,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  match self.row_range_groups {
    Some(groups) => self.with_row_range_groups(groups.resize(width, fix))
    None => self.with_row_ranges(granges_resize(self.row_ranges, width, fix))
  }
}

///|
/// Create flanking ranges while preserving assays and annotations.
pub fn RangedSummarizedExperiment::flank(
  self : RangedSummarizedExperiment,
  width : Int,
  start : Bool,
  both : Bool,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  match self.row_range_groups {
    Some(groups) => self.with_row_range_groups(groups.flank(width, start, both))
    None =>
      self.with_row_ranges(granges_flank(self.row_ranges, width, start, both))
  }
}

///|
/// Create strand-aware promoter ranges.
pub fn RangedSummarizedExperiment::promoters(
  self : RangedSummarizedExperiment,
  upstream? : Int = 2000,
  downstream? : Int = 200,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  if upstream < 0 || downstream < 0 {
    raise RangedSummarizedExperimentError(
      "promoter upstream and downstream widths must be non-negative",
    )
  }
  match self.row_range_groups {
    Some(groups) => {
      let promoted = groups.promoters(upstream~, downstream~) catch {
        GRangesListError(message) =>
          raise RangedSummarizedExperimentError(message)
      }
      self.with_row_range_groups(promoted)
    }
    None =>
      self.with_row_ranges(
        granges_promoters(self.row_ranges, upstream, downstream),
      )
  }
}

///|
/// Sort rows by sequence name, start, end, and strand.
pub fn RangedSummarizedExperiment::sort(
  self : RangedSummarizedExperiment,
  decreasing? : Bool = false,
) -> RangedSummarizedExperiment raise RangedSummarizedExperimentError {
  let indices : Array[Int] = []
  for index = 0; index < self.nrow(); index = index + 1 {
    indices.push(index)
  }
  let mut swapped = true
  while swapped {
    swapped = false
    for index = 0; index + 1 < indices.length(); index = index + 1 {
      let comparison = rse_compare_ranges(
        self.row_ranges,
        indices[index],
        indices[index + 1],
      )
      if (!decreasing && comparison > 0) || (decreasing && comparison < 0) {
        let temporary = indices[index]
        indices[index] = indices[index + 1]
        indices[index + 1] = temporary
        swapped = true
      }
    }
  }
  self.subset_rows(indices)
}

///|
/// Return a compact container summary.
pub fn RangedSummarizedExperiment::summary(
  self : RangedSummarizedExperiment,
) -> String {
  let range_label = if self.has_grouped_row_ranges() {
    " range groups x "
  } else {
    " ranges x "
  }
  "RangedSummarizedExperiment(" +
  self.nrow().to_string() +
  range_label +
  self.ncol().to_string() +
  " samples, assays=[" +
  self.assay_names().join(", ") +
  "])"
}

///|
fn rse_grouped_validation_error(
  assays : Map[String, Array[Array[Double]]],
  row_range_groups : GRangesList,
  col_data : Array[Map[String, String]],
  row_data : Array[Map[String, String]],
  row_names : Array[String],
) -> String? {
  match
    grl_validation_error(
      row_range_groups.elements,
      row_range_groups.names,
      row_range_groups.element_metadata,
    ) {
    Some(message) => return Some(message)
    None => ()
  }
  rse_validation_error(
    assays,
    row_range_groups.feature_bounds(),
    col_data,
    row_data,
    row_names,
  )
}

///|
fn rse_validation_error(
  assays : Map[String, Array[Array[Double]]],
  row_ranges : GRanges,
  col_data : Array[Map[String, String]],
  row_data : Array[Map[String, String]],
  row_names : Array[String],
) -> String? {
  let nrow = row_ranges.seqnames.length()
  if row_ranges.starts.length() != nrow ||
    row_ranges.ends.length() != nrow ||
    row_ranges.strands.length() != nrow ||
    row_ranges.widths.length() != nrow {
    return Some("all GRanges fields must have the same length")
  }
  for index = 0; index < nrow; index = index + 1 {
    let expected_width = row_ranges.ends[index] - row_ranges.starts[index] + 1
    if expected_width < 0 || row_ranges.widths[index] != expected_width {
      return Some("invalid genomic range width at row " + index.to_string())
    }
  }
  if row_data.length() != 0 && row_data.length() != nrow {
    return Some("row_data length must equal the number of row ranges")
  }
  if row_names.length() != 0 && row_names.length() != nrow {
    return Some("row_names length must equal the number of row ranges")
  }

  let assay_names = assays.keys().collect()
  for assay_name in assay_names {
    let assay = assays[assay_name]
    if assay.length() != nrow {
      return Some(
        "assay '" + assay_name + "' row count must equal rowRanges length",
      )
    }
    for row_index = 0; row_index < assay.length(); row_index = row_index + 1 {
      if assay[row_index].length() != col_data.length() {
        return Some(
          "assay '" + assay_name + "' column count must equal col_data length",
        )
      }
    }
  }
  None
}

///|
fn rse_ranges_equal(left : GRanges, right : GRanges) -> Bool {
  left.seqnames == right.seqnames &&
  left.starts == right.starts &&
  left.ends == right.ends &&
  left.strands == right.strands &&
  left.widths == right.widths
}

///|
fn rse_ranges_as_tuples(row_ranges : GRanges) -> Array[(String, Int, Int)] {
  let result : Array[(String, Int, Int)] = []
  for index = 0; index < granges_length(row_ranges); index = index + 1 {
    result.push(
      (
        row_ranges.seqnames[index],
        row_ranges.starts[index],
        row_ranges.ends[index],
      ),
    )
  }
  result
}

///|
fn rse_copy_assays(
  experiment : SummarizedExperiment,
) -> Map[String, Array[Array[Double]]] {
  let names = se_assay_names(experiment)
  let assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=names.length(),
  )
  for name in names {
    match se_assay(experiment, name) {
      Some(assay) => assays[name] = assay
      None => ()
    }
  }
  assays
}

///|
fn rse_validate_indices(
  indices : Array[Int],
  upper_bound : Int,
  axis : String,
) -> Unit raise RangedSummarizedExperimentError {
  for index in indices {
    if index < 0 || index >= upper_bound {
      raise RangedSummarizedExperimentError(
        axis + " index out of bounds: " + index.to_string(),
      )
    }
  }
}

///|
fn rse_take_ranges(row_ranges : GRanges, indices : Array[Int]) -> GRanges {
  let seqnames : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let strands : Array[Strand] = []
  let widths : Array[Int] = []
  for index in indices {
    seqnames.push(row_ranges.seqnames[index])
    starts.push(row_ranges.starts[index])
    ends.push(row_ranges.ends[index])
    strands.push(row_ranges.strands[index])
    widths.push(row_ranges.widths[index])
  }
  GRanges::{ seqnames, starts, ends, strands, widths }
}

///|
fn rse_take_range_groups(
  groups : GRangesList,
  indices : Array[Int],
) -> GRangesList {
  let elements = indices.map(fn(index) {
    grl_copy_ranges(groups.elements[index])
  })
  let names = if groups.names.length() == 0 {
    []
  } else {
    indices.map(fn(index) { groups.names[index] })
  }
  let metadata = if groups.element_metadata.length() == 0 {
    []
  } else {
    indices.map(fn(index) { groups.element_metadata[index] })
  }
  grl_create_unchecked(elements, names, metadata)
}

///|
fn rse_take_maps(
  values : Array[Map[String, String]],
  indices : Array[Int],
) -> Array[Map[String, String]] {
  if values.length() == 0 {
    return []
  }
  indices.map(fn(index) { values[index] })
}

///|
fn rse_take_strings(
  values : Array[String],
  indices : Array[Int],
) -> Array[String] {
  if values.length() == 0 {
    return []
  }
  indices.map(fn(index) { values[index] })
}

///|
fn rse_strands_compatible(
  query : Strand,
  subject : Strand,
  ignore_strand : Bool,
) -> Bool {
  ignore_strand ||
  query == Strand::Star ||
  subject == Strand::Star ||
  query == subject
}

///|
fn rse_ranges_overlap(
  query : GRanges,
  query_index : Int,
  subject : GRanges,
  subject_index : Int,
  ignore_strand : Bool,
) -> Bool {
  query.seqnames[query_index] == subject.seqnames[subject_index] &&
  rse_strands_compatible(
    query.strands[query_index],
    subject.strands[subject_index],
    ignore_strand,
  ) &&
  query.starts[query_index] <= subject.ends[subject_index] &&
  query.ends[query_index] >= subject.starts[subject_index]
}

///|
fn rse_interval_distance(
  first_start : Int,
  first_end : Int,
  second_start : Int,
  second_end : Int,
) -> Int {
  if first_end < second_start {
    second_start - first_end - 1
  } else if second_end < first_start {
    first_start - second_end - 1
  } else {
    0
  }
}

///|
fn rse_compare_ranges(row_ranges : GRanges, left : Int, right : Int) -> Int {
  if row_ranges.seqnames[left] < row_ranges.seqnames[right] {
    return -1
  }
  if row_ranges.seqnames[left] > row_ranges.seqnames[right] {
    return 1
  }
  if row_ranges.starts[left] < row_ranges.starts[right] {
    return -1
  }
  if row_ranges.starts[left] > row_ranges.starts[right] {
    return 1
  }
  if row_ranges.ends[left] < row_ranges.ends[right] {
    return -1
  }
  if row_ranges.ends[left] > row_ranges.ends[right] {
    return 1
  }
  rse_strand_rank(row_ranges.strands[left]) -
  rse_strand_rank(row_ranges.strands[right])
}

///|
fn rse_strand_rank(strand : Strand) -> Int {
  match strand {
    Strand::Plus => 0
    Strand::Minus => 1
    Strand::Star => 2
  }
}
