///|
/// GRangesList - grouped genomic ranges inspired by Bioconductor
/// GenomicRanges.

///|
pub suberror GRangesListError {
  GRangesListError(String)
}

///|
/// A collection of GRanges elements representing compound genomic features.
///
/// Typical elements are transcripts whose member ranges are exons. Names and
/// element metadata are optional, but when present they are parallel to the
/// outer list rather than to individual ranges.
pub struct GRangesList {
  elements : Array[GRanges]
  names : Array[String]
  element_metadata : Array[Map[String, String]]
}

///|
/// Construct and validate a grouped genomic-range collection.
pub fn GRangesList::new(
  elements~ : Array[GRanges],
  names? : Array[String] = [],
  element_metadata? : Array[Map[String, String]] = [],
) -> GRangesList raise GRangesListError {
  match grl_validation_error(elements, names, element_metadata) {
    Some(message) => raise GRangesListError(message)
    None => grl_create_unchecked(elements, names, element_metadata)
  }
}

///|
/// Function-style constructor matching the existing `granges` API.
pub fn granges_list(
  elements : Array[GRanges],
  names? : Array[String] = [],
  element_metadata? : Array[Map[String, String]] = [],
) -> GRangesList raise GRangesListError {
  GRangesList::new(elements~, names~, element_metadata~)
}

///|
/// Split a flat GRanges into named elements in first-appearance group order.
pub fn granges_split_as_list(
  ranges : GRanges,
  groups : Array[String],
) -> GRangesList raise GRangesListError {
  match grl_range_validation_error(ranges) {
    Some(message) => raise GRangesListError(message)
    None => ()
  }
  if groups.length() != granges_length(ranges) {
    raise GRangesListError(
      "group vector length must equal the number of genomic ranges",
    )
  }
  let names : Array[String] = []
  for group in groups {
    if !grl_strings_contains(names, group) {
      names.push(group)
    }
  }
  let elements : Array[GRanges] = []
  for name in names {
    let indices : Array[Int] = []
    for index = 0; index < groups.length(); index = index + 1 {
      if groups[index] == name {
        indices.push(index)
      }
    }
    elements.push(grl_take_ranges(ranges, indices))
  }
  grl_create_unchecked(elements, names, [])
}

///|
/// Partition a flat GRanges using explicit element lengths.
pub fn granges_list_from_partition(
  ranges : GRanges,
  element_lengths : Array[Int],
  names? : Array[String] = [],
  element_metadata? : Array[Map[String, String]] = [],
) -> GRangesList raise GRangesListError {
  match grl_range_validation_error(ranges) {
    Some(message) => raise GRangesListError(message)
    None => ()
  }
  if names.length() != 0 && names.length() != element_lengths.length() {
    raise GRangesListError(
      "names length must equal the number of list elements",
    )
  }
  if element_metadata.length() != 0 &&
    element_metadata.length() != element_lengths.length() {
    raise GRangesListError(
      "element metadata length must equal the number of list elements",
    )
  }
  let mut total = 0
  for length in element_lengths {
    if length < 0 {
      raise GRangesListError("element lengths must be non-negative")
    }
    total = total + length
  }
  if total != granges_length(ranges) {
    raise GRangesListError(
      "element lengths must sum to the number of genomic ranges",
    )
  }
  let elements : Array[GRanges] = []
  let mut offset = 0
  for length in element_lengths {
    elements.push(grl_take_contiguous(ranges, offset, length))
    offset = offset + length
  }
  grl_create_unchecked(elements, names, element_metadata)
}

///|
/// Return the number of compound features.
pub fn GRangesList::length(self : GRangesList) -> Int {
  self.elements.length()
}

///|
/// Return optional feature names. An empty array represents an unnamed list.
pub fn GRangesList::names(self : GRangesList) -> Array[String] {
  self.names.copy()
}

///|
/// Return optional feature-level metadata.
pub fn GRangesList::element_metadata(
  self : GRangesList,
) -> Array[Map[String, String]] {
  self.element_metadata.copy()
}

///|
/// Return copies of all GRanges elements.
pub fn GRangesList::elements(self : GRangesList) -> Array[GRanges] {
  self.elements.map(fn(element) { grl_copy_ranges(element) })
}

///|
/// Return the number of ranges in each compound feature.
pub fn GRangesList::element_lengths(self : GRangesList) -> Array[Int] {
  self.elements.map(fn(element) { granges_length(element) })
}

///|
/// Return true when there are no elements or all elements are empty.
pub fn GRangesList::is_empty(self : GRangesList) -> Bool {
  if self.elements.length() == 0 {
    return true
  }
  for element in self.elements {
    if granges_length(element) > 0 {
      return false
    }
  }
  true
}

///|
/// Validate all outer and inner parallel dimensions.
pub fn GRangesList::is_valid(self : GRangesList) -> Bool {
  grl_validation_error(self.elements, self.names, self.element_metadata) is None
}

///|
/// Return a copied GRanges element by zero-based index.
pub fn GRangesList::get(self : GRangesList, index : Int) -> GRanges? {
  if index < 0 || index >= self.elements.length() {
    None
  } else {
    Some(grl_copy_ranges(self.elements[index]))
  }
}

///|
/// Return the first element with a matching name.
pub fn GRangesList::get_by_name(self : GRangesList, name : String) -> GRanges? {
  for index = 0; index < self.names.length(); index = index + 1 {
    if self.names[index] == name {
      return Some(grl_copy_ranges(self.elements[index]))
    }
  }
  None
}

///|
/// Select, reorder, or duplicate outer list elements.
pub fn GRangesList::subset(
  self : GRangesList,
  indices : Array[Int],
) -> GRangesList raise GRangesListError {
  grl_validate_indices(indices, self.elements.length())
  let elements = indices.map(fn(index) { grl_copy_ranges(self.elements[index]) })
  let names = if self.names.length() == 0 {
    []
  } else {
    indices.map(fn(index) { self.names[index] })
  }
  let metadata = if self.element_metadata.length() == 0 {
    []
  } else {
    indices.map(fn(index) { self.element_metadata[index] })
  }
  grl_create_unchecked(elements, names, metadata)
}

///|
/// Concatenate all member ranges into a flat GRanges.
pub fn GRangesList::unlist(self : GRangesList) -> GRanges {
  let result = grl_empty_ranges()
  for element in self.elements {
    grl_append_ranges(result, element)
  }
  result
}

///|
/// Repartition a flat GRanges using this object's element boundaries.
pub fn GRangesList::relist(
  self : GRangesList,
  ranges : GRanges,
) -> GRangesList raise GRangesListError {
  granges_list_from_partition(
    ranges,
    self.element_lengths(),
    names=self.names,
    element_metadata=self.element_metadata,
  )
}

///|
/// Summarize each feature as one bounding range.
///
/// Features spanning multiple sequence names receive `"*"` as the sequence
/// name, and mixed strands receive `Strand::Star`. Empty features are encoded
/// as a zero-width range at 0.
pub fn GRangesList::feature_bounds(self : GRangesList) -> GRanges {
  let result = grl_empty_ranges()
  for element in self.elements {
    if granges_length(element) == 0 {
      result.seqnames.push("")
      result.starts.push(0)
      result.ends.push(-1)
      result.strands.push(Strand::Star)
      result.widths.push(0)
    } else {
      let mut seqname = element.seqnames[0]
      let mut strand = element.strands[0]
      let mut start = element.starts[0]
      let mut end = element.ends[0]
      for index = 1; index < granges_length(element); index = index + 1 {
        if element.seqnames[index] != seqname {
          seqname = "*"
        }
        if element.strands[index] != strand {
          strand = Strand::Star
        }
        if element.starts[index] < start {
          start = element.starts[index]
        }
        if element.ends[index] > end {
          end = element.ends[index]
        }
      }
      result.seqnames.push(seqname)
      result.starts.push(start)
      result.ends.push(end)
      result.strands.push(strand)
      result.widths.push(end - start + 1)
    }
  }
  result
}

///|
/// Append the outer elements of two lists.
pub fn GRangesList::concat(
  self : GRangesList,
  other : GRangesList,
) -> GRangesList {
  let elements = self.elements()
  for element in other.elements {
    elements.push(grl_copy_ranges(element))
  }
  let names : Array[String] = []
  if self.names.length() > 0 || other.names.length() > 0 {
    for index = 0; index < self.length(); index = index + 1 {
      names.push(if self.names.length() > 0 { self.names[index] } else { "" })
    }
    for index = 0; index < other.length(); index = index + 1 {
      names.push(if other.names.length() > 0 { other.names[index] } else { "" })
    }
  }
  let metadata : Array[Map[String, String]] = []
  if self.element_metadata.length() > 0 || other.element_metadata.length() > 0 {
    for index = 0; index < self.length(); index = index + 1 {
      metadata.push(
        if self.element_metadata.length() > 0 {
          self.element_metadata[index]
        } else {
          Map([])
        },
      )
    }
    for index = 0; index < other.length(); index = index + 1 {
      metadata.push(
        if other.element_metadata.length() > 0 {
          other.element_metadata[index]
        } else {
          Map([])
        },
      )
    }
  }
  grl_create_unchecked(elements, names, metadata)
}

///|
/// Concatenate corresponding GRanges elements.
pub fn GRangesList::parallel_concat(
  self : GRangesList,
  other : GRangesList,
) -> GRangesList raise GRangesListError {
  grl_require_parallel(self, other)
  let elements : Array[GRanges] = []
  for index = 0; index < self.length(); index = index + 1 {
    elements.push(
      grl_concat_ranges(self.elements[index], other.elements[index]),
    )
  }
  let names = grl_preferred_names(self, other)
  let metadata = grl_preferred_metadata(self, other)
  grl_create_unchecked(elements, names, metadata)
}

///|
/// Shift every member range while preserving group boundaries.
pub fn GRangesList::shift(self : GRangesList, amount : Int) -> GRangesList {
  grl_map(self, fn(element) { granges_shift(element, amount) })
}

///|
/// Narrow every member range using one-based relative coordinates.
pub fn GRangesList::narrow(
  self : GRangesList,
  start : Int,
  end : Int,
) -> GRangesList {
  grl_map(self, fn(element) { granges_narrow(element, start, end) })
}

///|
/// Resize every member range.
pub fn GRangesList::resize(
  self : GRangesList,
  width : Int,
  fix : String,
) -> GRangesList {
  grl_map(self, fn(element) { granges_resize(element, width, fix) })
}

///|
/// Create flanking ranges for every member range.
pub fn GRangesList::flank(
  self : GRangesList,
  width : Int,
  start : Bool,
  both : Bool,
) -> GRangesList {
  grl_map(self, fn(element) { granges_flank(element, width, start, both) })
}

///|
/// Create strand-aware promoter ranges for every member range.
pub fn GRangesList::promoters(
  self : GRangesList,
  upstream? : Int = 2000,
  downstream? : Int = 200,
) -> GRangesList raise GRangesListError {
  if upstream < 0 || downstream < 0 {
    raise GRangesListError(
      "promoter upstream and downstream widths must be non-negative",
    )
  }
  grl_map(self, fn(element) { granges_promoters(element, upstream, downstream) })
}

///|
/// Reduce overlapping ranges independently within each feature.
pub fn GRangesList::reduce(
  self : GRangesList,
  min_gapwidth? : Int = 0,
) -> GRangesList {
  grl_map(self, fn(element) { granges_reduce(element, min_gapwidth) })
}

///|
/// Disjoin ranges independently within each feature.
pub fn GRangesList::disjoin(self : GRangesList) -> GRangesList {
  grl_map(self, fn(element) { granges_disjoin(element) })
}

///|
/// Sort ranges inside each feature without changing outer feature order.
pub fn GRangesList::sort_ranges(self : GRangesList) -> GRangesList {
  grl_map(self, fn(element) { granges_sort(element) })
}

///|
/// Calculate the sum of member widths for every compound feature.
pub fn GRangesList::total_widths(self : GRangesList) -> Array[Int] {
  self.elements.map(fn(element) {
    let mut total = 0
    for width in element.widths {
      total = total + width
    }
    total
  })
}

///|
/// Compute an element-wise union of two parallel lists.
pub fn GRangesList::parallel_union(
  self : GRangesList,
  other : GRangesList,
) -> GRangesList raise GRangesListError {
  grl_parallel_set_operation(self, other, fn(left, right) {
    granges_union(left, right)
  })
}

///|
/// Compute an element-wise intersection of two parallel lists.
pub fn GRangesList::parallel_intersect(
  self : GRangesList,
  other : GRangesList,
) -> GRangesList raise GRangesListError {
  grl_parallel_set_operation(self, other, fn(left, right) {
    granges_intersect(left, right)
  })
}

///|
/// Compute an element-wise asymmetric difference of two parallel lists.
pub fn GRangesList::parallel_setdiff(
  self : GRangesList,
  other : GRangesList,
) -> GRangesList raise GRangesListError {
  grl_parallel_set_operation(self, other, fn(left, right) {
    granges_setdiff(left, right)
  })
}

///|
/// Find subject ranges overlapping each compound query feature.
///
/// A feature-subject pair is returned once even when multiple member ranges
/// overlap the same subject range.
pub fn GRangesList::find_overlaps(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[(Int, Int)] {
  let hits : Array[(Int, Int)] = []
  for query_index = 0
      query_index < self.elements.length()
      query_index = query_index + 1 {
    for subject_index = 0
        subject_index < granges_length(subject)
        subject_index = subject_index + 1 {
      if grl_element_overlaps_range(
          self.elements[query_index],
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
/// Count distinct subject-range overlaps for each feature.
pub fn GRangesList::count_overlaps(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  let counts = Array::make(self.length(), 0)
  for hit in self.find_overlaps(subject, ignore_strand~) {
    counts[hit.0] = counts[hit.0] + 1
  }
  counts
}

///|
/// Return one overlap flag per compound feature.
pub fn GRangesList::overlaps_any(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Bool] {
  self.count_overlaps(subject, ignore_strand~).map(fn(count) { count > 0 })
}

///|
/// Keep features that overlap at least one subject range.
pub fn GRangesList::subset_by_overlaps(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> GRangesList raise GRangesListError {
  let indices : Array[Int] = []
  let keep = self.overlaps_any(subject, ignore_strand~)
  for index = 0; index < keep.length(); index = index + 1 {
    if keep[index] {
      indices.push(index)
    }
  }
  self.subset(indices)
}

///|
/// Find overlaps between compound features in two grouped collections.
pub fn GRangesList::find_group_overlaps(
  self : GRangesList,
  subject : GRangesList,
  ignore_strand? : Bool = false,
) -> Array[(Int, Int)] {
  let hits : Array[(Int, Int)] = []
  for query_index = 0
      query_index < self.length()
      query_index = query_index + 1 {
    for subject_index = 0
        subject_index < subject.length()
        subject_index = subject_index + 1 {
      if grl_elements_overlap(
          self.elements[query_index],
          subject.elements[subject_index],
          ignore_strand,
        ) {
        hits.push((query_index, subject_index))
      }
    }
  }
  hits
}

///|
/// Count distinct subject features overlapping each query feature.
pub fn GRangesList::count_group_overlaps(
  self : GRangesList,
  subject : GRangesList,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  let counts = Array::make(self.length(), 0)
  for hit in self.find_group_overlaps(subject, ignore_strand~) {
    counts[hit.0] = counts[hit.0] + 1
  }
  counts
}

///|
/// Find overlaps between flat query ranges and grouped subject features.
pub fn granges_find_overlaps_list(
  query : GRanges,
  subject : GRangesList,
  ignore_strand? : Bool = false,
) -> Array[(Int, Int)] {
  let hits : Array[(Int, Int)] = []
  for query_index = 0
      query_index < granges_length(query)
      query_index = query_index + 1 {
    for subject_index = 0
        subject_index < subject.length()
        subject_index = subject_index + 1 {
      if grl_element_overlaps_range(
          subject.elements[subject_index],
          query,
          query_index,
          ignore_strand,
        ) {
        hits.push((query_index, subject_index))
      }
    }
  }
  hits
}

///|
/// Return the nearest subject range for each compound feature.
pub fn GRangesList::nearest(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  let result = Array::make(self.length(), -1)
  for query_index = 0
      query_index < self.length()
      query_index = query_index + 1 {
    let nearest = grl_nearest_for_element(
      self.elements[query_index],
      subject,
      ignore_strand,
    )
    result[query_index] = nearest.0
  }
  result
}

///|
/// Return the minimum genomic distance to a subject range for each feature.
pub fn GRangesList::distance_to_nearest(
  self : GRangesList,
  subject : GRanges,
  ignore_strand? : Bool = false,
) -> Array[Int] {
  let result = Array::make(self.length(), -1)
  for query_index = 0
      query_index < self.length()
      query_index = query_index + 1 {
    let nearest = grl_nearest_for_element(
      self.elements[query_index],
      subject,
      ignore_strand,
    )
    result[query_index] = nearest.1
  }
  result
}

///|
/// Compute per-base coverage across all member ranges.
pub fn GRangesList::coverage(
  self : GRangesList,
  sequence_lengths : Map[String, Int],
) -> Map[String, Array[Int]] {
  coverage_granges(self.unlist(), sequence_lengths)
}

///|
/// Return a compact list summary.
pub fn GRangesList::summary(self : GRangesList) -> String {
  let mut total_ranges = 0
  for length in self.element_lengths() {
    total_ranges = total_ranges + length
  }
  "GRangesList(" +
  self.length().to_string() +
  " elements, " +
  total_ranges.to_string() +
  " ranges)"
}

///|
fn grl_create_unchecked(
  elements : Array[GRanges],
  names : Array[String],
  element_metadata : Array[Map[String, String]],
) -> GRangesList {
  GRangesList::{
    elements: elements.map(fn(element) { grl_copy_ranges(element) }),
    names: names.copy(),
    element_metadata: element_metadata.copy(),
  }
}

///|
fn grl_validation_error(
  elements : Array[GRanges],
  names : Array[String],
  element_metadata : Array[Map[String, String]],
) -> String? {
  if names.length() != 0 && names.length() != elements.length() {
    return Some("names length must equal the number of list elements")
  }
  if element_metadata.length() != 0 &&
    element_metadata.length() != elements.length() {
    return Some(
      "element metadata length must equal the number of list elements",
    )
  }
  for index = 0; index < elements.length(); index = index + 1 {
    match grl_range_validation_error(elements[index]) {
      Some(message) =>
        return Some(
          "invalid GRanges element " + index.to_string() + ": " + message,
        )
      None => ()
    }
  }
  None
}

///|
fn grl_range_validation_error(ranges : GRanges) -> String? {
  let length = ranges.seqnames.length()
  if ranges.starts.length() != length ||
    ranges.ends.length() != length ||
    ranges.strands.length() != length ||
    ranges.widths.length() != length {
    return Some("all GRanges fields must have the same length")
  }
  for index = 0; index < length; index = index + 1 {
    let expected_width = ranges.ends[index] - ranges.starts[index] + 1
    if expected_width < 0 || ranges.widths[index] != expected_width {
      return Some("invalid genomic range width at index " + index.to_string())
    }
  }
  None
}

///|
fn grl_validate_indices(
  indices : Array[Int],
  upper_bound : Int,
) -> Unit raise GRangesListError {
  for index in indices {
    if index < 0 || index >= upper_bound {
      raise GRangesListError(
        "list element index out of bounds: " + index.to_string(),
      )
    }
  }
}

///|
fn grl_empty_ranges() -> GRanges {
  GRanges::{ seqnames: [], starts: [], ends: [], strands: [], widths: [], }
}

///|
fn grl_copy_ranges(ranges : GRanges) -> GRanges {
  GRanges::{
    seqnames: ranges.seqnames.copy(),
    starts: ranges.starts.copy(),
    ends: ranges.ends.copy(),
    strands: ranges.strands.copy(),
    widths: ranges.widths.copy(),
  }
}

///|
fn grl_take_ranges(ranges : GRanges, indices : Array[Int]) -> GRanges {
  let result = grl_empty_ranges()
  for index in indices {
    result.seqnames.push(ranges.seqnames[index])
    result.starts.push(ranges.starts[index])
    result.ends.push(ranges.ends[index])
    result.strands.push(ranges.strands[index])
    result.widths.push(ranges.widths[index])
  }
  result
}

///|
fn grl_take_contiguous(ranges : GRanges, offset : Int, length : Int) -> GRanges {
  let indices : Array[Int] = []
  for index = offset; index < offset + length; index = index + 1 {
    indices.push(index)
  }
  grl_take_ranges(ranges, indices)
}

///|
fn grl_append_ranges(target : GRanges, source : GRanges) -> Unit {
  for index = 0; index < granges_length(source); index = index + 1 {
    target.seqnames.push(source.seqnames[index])
    target.starts.push(source.starts[index])
    target.ends.push(source.ends[index])
    target.strands.push(source.strands[index])
    target.widths.push(source.widths[index])
  }
}

///|
fn grl_concat_ranges(left : GRanges, right : GRanges) -> GRanges {
  let result = grl_copy_ranges(left)
  grl_append_ranges(result, right)
  result
}

///|
fn grl_strings_contains(values : Array[String], target : String) -> Bool {
  for value in values {
    if value == target {
      return true
    }
  }
  false
}

///|
fn grl_map(
  ranges : GRangesList,
  operation : (GRanges) -> GRanges,
) -> GRangesList {
  grl_create_unchecked(
    ranges.elements.map(fn(element) { operation(element) }),
    ranges.names,
    ranges.element_metadata,
  )
}

///|
fn grl_require_parallel(
  left : GRangesList,
  right : GRangesList,
) -> Unit raise GRangesListError {
  if left.length() != right.length() {
    raise GRangesListError(
      "parallel GRangesList operations require equal element counts",
    )
  }
}

///|
fn grl_preferred_names(
  left : GRangesList,
  right : GRangesList,
) -> Array[String] {
  if left.names.length() > 0 {
    left.names
  } else {
    right.names
  }
}

///|
fn grl_preferred_metadata(
  left : GRangesList,
  right : GRangesList,
) -> Array[Map[String, String]] {
  if left.element_metadata.length() > 0 {
    left.element_metadata
  } else {
    right.element_metadata
  }
}

///|
fn grl_parallel_set_operation(
  left : GRangesList,
  right : GRangesList,
  operation : (GRanges, GRanges) -> GRanges,
) -> GRangesList raise GRangesListError {
  grl_require_parallel(left, right)
  let elements : Array[GRanges] = []
  for index = 0; index < left.length(); index = index + 1 {
    elements.push(operation(left.elements[index], right.elements[index]))
  }
  grl_create_unchecked(
    elements,
    grl_preferred_names(left, right),
    grl_preferred_metadata(left, right),
  )
}

///|
fn grl_strands_compatible(
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
fn grl_ranges_overlap(
  query : GRanges,
  query_index : Int,
  subject : GRanges,
  subject_index : Int,
  ignore_strand : Bool,
) -> Bool {
  query.seqnames[query_index] == subject.seqnames[subject_index] &&
  grl_strands_compatible(
    query.strands[query_index],
    subject.strands[subject_index],
    ignore_strand,
  ) &&
  query.starts[query_index] <= subject.ends[subject_index] &&
  query.ends[query_index] >= subject.starts[subject_index]
}

///|
fn grl_element_overlaps_range(
  element : GRanges,
  subject : GRanges,
  subject_index : Int,
  ignore_strand : Bool,
) -> Bool {
  for member_index = 0
      member_index < granges_length(element)
      member_index = member_index + 1 {
    if grl_ranges_overlap(
        element, member_index, subject, subject_index, ignore_strand,
      ) {
      return true
    }
  }
  false
}

///|
fn grl_elements_overlap(
  query : GRanges,
  subject : GRanges,
  ignore_strand : Bool,
) -> Bool {
  for subject_index = 0
      subject_index < granges_length(subject)
      subject_index = subject_index + 1 {
    if grl_element_overlaps_range(query, subject, subject_index, ignore_strand) {
      return true
    }
  }
  false
}

///|
fn grl_interval_distance(
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
fn grl_nearest_for_element(
  element : GRanges,
  subject : GRanges,
  ignore_strand : Bool,
) -> (Int, Int) {
  let mut best_index = -1
  let mut best_distance = -1
  for member_index = 0
      member_index < granges_length(element)
      member_index = member_index + 1 {
    for subject_index = 0
        subject_index < granges_length(subject)
        subject_index = subject_index + 1 {
      if element.seqnames[member_index] == subject.seqnames[subject_index] &&
        grl_strands_compatible(
          element.strands[member_index],
          subject.strands[subject_index],
          ignore_strand,
        ) {
        let distance = grl_interval_distance(
          element.starts[member_index],
          element.ends[member_index],
          subject.starts[subject_index],
          subject.ends[subject_index],
        )
        if best_index < 0 ||
          distance < best_distance ||
          (distance == best_distance && subject_index < best_index) {
          best_index = subject_index
          best_distance = distance
        }
      }
    }
  }
  (best_index, best_distance)
}
