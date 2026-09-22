///|
/// gffutils port — feature interval predicates, generic parent/child traversal
/// (via the `ID`/`Parent` attributes), and merging of overlapping features,
/// operating on BioSeqs' `GFFFeature`/`GFFRecord`.

///|
/// Two features overlap (same `seqid`, 1-based inclusive intervals intersect).
pub fn gffutils_overlaps(f1 : GFFFeature, f2 : GFFFeature) -> Bool {
  f1.seqid == f2.seqid && f1.start <= f2.end && f1.end >= f2.start
}

///|
/// Alias for `gffutils_overlaps`.
pub fn gffutils_intersects(f1 : GFFFeature, f2 : GFFFeature) -> Bool {
  gffutils_overlaps(f1, f2)
}

///|
/// `f1` fully contains `f2` (same `seqid`).
pub fn gffutils_contains(f1 : GFFFeature, f2 : GFFFeature) -> Bool {
  f1.seqid == f2.seqid && f1.start <= f2.start && f1.end >= f2.end
}

///|
/// `f1` lies fully within `f2`.
pub fn gffutils_within(f1 : GFFFeature, f2 : GFFFeature) -> Bool {
  gffutils_contains(f2, f1)
}

///|
/// Distance in bases between two features (0 if they overlap; otherwise the
/// gap length in 1-based inclusive coordinates).
pub fn gffutils_distance(f1 : GFFFeature, f2 : GFFFeature) -> Int {
  if gffutils_overlaps(f1, f2) {
    0
  } else if f1.end < f2.start {
    f2.start - f1.end - 1
  } else {
    f1.start - f2.end - 1
  }
}

///|
fn gffutils_trim(s : String) -> String {
  let mut lo = 0
  let mut hi = s.length() - 1
  while lo <= hi {
    let c = s.unsafe_get(lo).to_int()
    if c == 32 || c == 9 {
      lo = lo + 1
    } else {
      break
    }
  }
  while hi >= lo {
    let c = s.unsafe_get(hi).to_int()
    if c == 32 || c == 9 {
      hi = hi - 1
    } else {
      break
    }
  }
  s[lo:hi + 1].to_owned()
}

///|
/// Parent IDs of a feature, splitting a comma-separated `Parent` attribute.
fn gffutils_parent_ids(f : GFFFeature) -> Array[String] {
  match f.attributes.get("Parent") {
    None => []
    Some(p) => {
      let out : Array[String] = Array::new()
      for part in p.split(",") {
        let t = gffutils_trim(part.to_owned())
        if t.length() > 0 {
          out.push(t)
        }
      }
      out
    }
  }
}

///|
/// Immediate children of a feature: features whose `Parent` attribute lists the
/// given feature ID.
pub fn gffutils_children(record : GFFRecord, feature_id : String) -> Array[GFFFeature] {
  let out : Array[GFFFeature] = Array::new()
  for f in record.features {
    for p in gffutils_parent_ids(f) {
      if p == feature_id {
        out.push(f)
        break
      }
    }
  }
  out
}

///|
/// Immediate parents of a feature: the features listed in its `Parent`
/// attribute.
pub fn gffutils_parents(record : GFFRecord, feature_id : String) -> Array[GFFFeature] {
  let wanted : Map[String, Bool] = Map([], capacity=0)
  for f in record.features {
    match f.get_id() {
      None => ()
      Some(id) =>
        if id == feature_id {
          for p in gffutils_parent_ids(f) {
            wanted.set(p, true)
          }
        }
    }
  }
  let out : Array[GFFFeature] = Array::new()
  for f in record.features {
    match f.get_id() {
      None => ()
      Some(id) =>
        match wanted.get(id) {
          None => ()
          Some(_) => out.push(f)
        }
    }
  }
  out
}

///|
/// All transitive descendants of a feature (breadth-first via the `Parent`
/// attribute).
pub fn gffutils_all_children(record : GFFRecord, feature_id : String) -> Array[GFFFeature] {
  let seen : Map[String, Bool] = Map([], capacity=0)
  let out : Array[GFFFeature] = Array::new()
  let mut frontier = gffutils_children(record, feature_id)
  while frontier.length() > 0 {
    let next : Array[GFFFeature] = Array::new()
    for f in frontier {
      let key = match f.get_id() {
        None => continue
        Some(id) => id
      }
      match seen.get(key) {
        Some(_) => continue
        None => ()
      }
      seen.set(key, true)
      out.push(f)
      for c in gffutils_children(record, key) {
        next.push(c)
      }
    }
    frontier = next
  }
  out
}

///|
/// All transitive ancestors of a feature (via the `Parent` attribute).
pub fn gffutils_all_parents(record : GFFRecord, feature_id : String) -> Array[GFFFeature] {
  let seen : Map[String, Bool] = Map([], capacity=0)
  let out : Array[GFFFeature] = Array::new()
  let mut frontier = gffutils_parents(record, feature_id)
  while frontier.length() > 0 {
    let next : Array[GFFFeature] = Array::new()
    for f in frontier {
      let key = match f.get_id() {
        None => continue
        Some(id) => id
      }
      match seen.get(key) {
        Some(_) => continue
        None => ()
      }
      seen.set(key, true)
      out.push(f)
      for p in gffutils_parents(record, key) {
        next.push(p)
      }
    }
    frontier = next
  }
  out
}


///|
fn gffutils_lt(a : GFFFeature, b : GFFFeature) -> Bool {
  if a.seqid != b.seqid {
    a.seqid < b.seqid
  } else if a.feature_type != b.feature_type {
    a.feature_type < b.feature_type
  } else if a.strand != b.strand {
    a.strand < b.strand
  } else if a.start != b.start {
    a.start < b.start
  } else {
    a.end < b.end
  }
}

///|
fn gffutils_sort_features(record : GFFRecord) -> Array[GFFFeature] {
  let out = record.features.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && gffutils_lt(key, out[j]) {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Merge overlapping (or book-ended) features of the same `seqid`, feature type
/// and strand into larger features, returning the merged features
/// (`gffutils.create_db(...).merge()`). Merged features carry the first
/// feature's source/type strand and an empty attribute set.
pub fn gffutils_merge(record : GFFRecord) -> Array[GFFFeature] {
  let sorted = gffutils_sort_features(record)
  let out : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < sorted.length() {
    let base = sorted[i]
    let mut merged = base
    let mut j = i + 1
    while j < sorted.length() {
      let f = sorted[j]
      if f.seqid == merged.seqid && f.feature_type == merged.feature_type &&
        f.strand == merged.strand && f.start <= merged.end + 1 {
        merged = GFFFeature::{
          seqid: merged.seqid,
          source: merged.source,
          feature_type: merged.feature_type,
          start: merged.start,
          end: if f.end > merged.end { f.end } else { merged.end },
          score: None,
          strand: merged.strand,
          phase: None,
          attributes: Map([], capacity=0),
        }
        j = j + 1
      } else {
        break
      }
    }
    out.push(merged)
    i = j
  }
  out
}