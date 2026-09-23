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
      let out : Array[String] = Array(capacity=0)
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
  let out : Array[GFFFeature] = Array(capacity=0)
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
  let out : Array[GFFFeature] = Array(capacity=0)
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
  let out : Array[GFFFeature] = Array(capacity=0)
  let mut frontier = gffutils_children(record, feature_id)
  while frontier.length() > 0 {
    let next : Array[GFFFeature] = Array(capacity=0)
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
  let out : Array[GFFFeature] = Array(capacity=0)
  let mut frontier = gffutils_parents(record, feature_id)
  while frontier.length() > 0 {
    let next : Array[GFFFeature] = Array(capacity=0)
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
/// Merge features matching the default criteria of `gffutils.FeatureDB.merge()`
/// — same `seqid`, feature type and strand, with `cur.start <= acc.end +
/// min_dist` — into larger features. A merged feature carries the first
/// component's attributes plus a fresh `ID` of the form `<featuretype>_<N>` and
/// a comma-joined `source`; standalone features are returned unchanged.
pub fn gffutils_merge(record : GFFRecord, min_dist? : Int = 1) -> Array[GFFFeature] {
  let sorted = gffutils_sort_features(record)
  let out : Array[GFFFeature] = Array(capacity=0)
  let counters : Map[String, Int] = Map([], capacity=0)
  let mut i = 0
  while i < sorted.length() {
    let first = sorted[i]
    let mut max_end = first.end
    let children : Array[GFFFeature] = Array(capacity=0)
    children.push(first)
    let mut j = i + 1
    while j < sorted.length() {
      let cur = sorted[j]
      if cur.seqid == first.seqid && cur.strand == first.strand &&
        cur.feature_type == first.feature_type && cur.start <= max_end + min_dist {
        if cur.end > max_end {
          max_end = cur.end
        }
        children.push(cur)
        j = j + 1
      } else {
        break
      }
    }
    if children.length() == 1 {
      out.push(first)
    } else {
      let ft = first.feature_type
      let n = match counters.get(ft) {
        None => {
          counters.set(ft, 1)
          1
        }
        Some(c) => {
          counters.set(ft, c + 1)
          c + 1
        }
      }
      let new_id = ft + "_" + n.to_string()
      let seen_src : Map[String, Bool] = Map([], capacity=0)
      let src_sb = StringBuilder()
      let mut first_src = true
      for ch in children {
        match seen_src.get(ch.source) {
          Some(_) => ()
          None => {
            seen_src.set(ch.source, true)
            if !(first_src) {
              src_sb.write_char(',')
            }
            src_sb.write_string(ch.source)
            first_src = false
          }
        }
      }
      let attrs : Map[String, String] = Map([], capacity=0)
      for (k, v) in first.attributes {
        attrs.set(k, v)
      }
      attrs.set("ID", new_id)
      out.push(
        GFFFeature::{
          seqid: first.seqid,
          source: src_sb.to_string(),
          feature_type: ft,
          start: first.start,
          end: max_end,
          score: first.score,
          strand: first.strand,
          phase: first.phase,
          attributes: attrs,
        },
      )
    }
    i = j
  }
  out
}

///|
/// Convert a `GFFFeature` to a BED line (0-based half-open coordinates), with
/// columns `chrom, chromStart, chromEnd, name, score, strand`.
pub fn gffutils_to_bed(feature : GFFFeature) -> String {
  let name = match feature.get_id() {
    None => "."
    Some(id) => id
  }
  let score = match feature.score {
    None => "."
    Some(s) => s.to_string()
  }
  let strand = if feature.strand.length() == 0 { "." } else { feature.strand }
  feature.seqid + "\t" + (feature.start - 1).to_string() + "\t" +
    feature.end.to_string() + "\t" + name + "\t" + score + "\t" + strand
}

///|
/// Convert all features of a `GFFRecord` to BED text (one line per feature).
pub fn gffutils_record_to_bed(record : GFFRecord) -> String {
  let sb = StringBuilder()
  let mut first = true
  for f in record.features {
    if !(first) {
      sb.write_char('\n')
    }
    sb.write_string(gffutils_to_bed(f))
    first = false
  }
  sb.to_string()
}