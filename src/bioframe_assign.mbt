/// Assign view labels to intervals and select by labels.
///
/// Bioframe equivalent: `assign_view` + `select_labels`
///
/// This module provides label-based interval annotation: assign each
/// interval a label based on which view region it falls into, then
/// select/drop intervals by label.

///|
/// An annotated interval: original interval + assigned label
pub struct BfAnnotatedInterval {
  interval : BedInterval
  label : String
}

///|
pub fn BfAnnotatedInterval::interval(self : BfAnnotatedInterval) -> BedInterval {
  self.interval
}

///|
pub fn BfAnnotatedInterval::label(self : BfAnnotatedInterval) -> String {
  self.label
}

///|
/// Assign each interval a label based on the view region it falls into.
/// Intervals not contained in any view region get label "NA".
pub fn bf_assign_labels(
  intervals : Array[BedInterval],
  view : BfViewFrame,
) -> Array[BfAnnotatedInterval] {
  let result : Array[BfAnnotatedInterval] = Array(capacity=0)
  for iv in intervals {
    let mut label = "NA"
    for vr in view.regions {
      if iv.chrom == vr.chrom && vr.start <= iv.start && iv.end <= vr.end {
        label = vr.name
        break
      }
    }
    result.push({ interval: iv, label, })
  }
  result
}

///|
/// Select annotated intervals by label.
pub fn bf_select_by_label(
  annotated : Array[BfAnnotatedInterval],
  labels : Array[String],
) -> Array[BfAnnotatedInterval] {
  let label_set : Set[String] = Set::from_array(labels)
  let result : Array[BfAnnotatedInterval] = Array(capacity=0)
  for ai in annotated {
    if label_set.contains(ai.label) {
      result.push(ai)
    }
  }
  result
}

///|
/// Drop annotated intervals by label.
pub fn bf_drop_by_label(
  annotated : Array[BfAnnotatedInterval],
  labels : Array[String],
) -> Array[BfAnnotatedInterval] {
  let label_set : Set[String] = Set::from_array(labels)
  let result : Array[BfAnnotatedInterval] = Array(capacity=0)
  for ai in annotated {
    if !label_set.contains(ai.label) {
      result.push(ai)
    }
  }
  result
}

///|
/// Get all unique labels from annotated intervals.
pub fn bf_get_labels(annotated : Array[BfAnnotatedInterval]) -> Array[String] {
  let label_set : Set[String] = Set([])
  for ai in annotated {
    label_set.add(ai.label)
  }
  label_set.to_array()
}

///|
/// Count intervals per label.
pub fn bf_count_by_label(
  annotated : Array[BfAnnotatedInterval],
) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([])
  for ai in annotated {
    match counts.get(ai.label) {
      Some(c) => counts[ai.label] = c + 1
      None => counts[ai.label] = 1
    }
  }
  counts
}

///|
/// Group annotated intervals by label.
pub fn bf_group_by_label(
  annotated : Array[BfAnnotatedInterval],
) -> Map[String, Array[BedInterval]] {
  let groups : Map[String, Array[BedInterval]] = Map([])
  for ai in annotated {
    match groups.get(ai.label) {
      Some(arr) => arr.push(ai.interval)
      None => groups[ai.label] = [ai.interval]
    }
  }
  groups
}

///|
/// Extract plain intervals from annotated intervals.
pub fn bf_strip_annotations(
  annotated : Array[BfAnnotatedInterval],
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for ai in annotated {
    result.push(ai.interval)
  }
  result
}
