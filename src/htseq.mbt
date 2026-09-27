///|
/// HTSeq — read counting for RNA-seq (htseq-count).
///
/// Assigns aligned reads from a BAM file to annotation features (e.g. genes or
/// exons) and tallies per-feature counts plus the standard "not counted"
/// categories (`no_feature`, `ambiguous`, `not_aligned`). Modes follow
/// htseq-count: `union` / `intersection-nonempty` (any overlap) and
/// `intersection-strict` (read fully contained in the feature).

///|
/// Result of `htseq_count`.
pub struct HtseqCounts {
  counts : Map[String, Int]
  no_feature : Int
  ambiguous : Int
  not_aligned : Int
}

///|
/// Per-feature read counts.
pub fn HtseqCounts::get_counts(self : HtseqCounts) -> Map[String, Int] {
  self.counts
}

///|
/// Reads overlapping no feature.
pub fn HtseqCounts::get_no_feature(self : HtseqCounts) -> Int {
  self.no_feature
}

///|
/// Reads overlapping multiple features (ambiguous under the current mode).
pub fn HtseqCounts::get_ambiguous(self : HtseqCounts) -> Int {
  self.ambiguous
}

///|
/// Reads that were unmapped.
pub fn HtseqCounts::get_not_aligned(self : HtseqCounts) -> Int {
  self.not_aligned
}

///|
/// Count of one feature, or 0 if absent.
pub fn HtseqCounts::count_of(self : HtseqCounts, feature_id : String) -> Int {
  match self.counts.get(feature_id) {
    None => 0
    Some(c) => c
  }
}

///|
/// A feature's identifier (GFF3 `ID` attribute, else a positional fallback).
fn htseq_feature_id(f : GFFFeature, index : Int) -> String {
  match f.get_id() {
    Some(id) => id
    None => "feature_" + index.to_string()
  }
}

///|
/// Overlap test between a read (1-based inclusive) and a GFF feature.
fn htseq_overlap(
  start_1b : Int,
  end_1b : Int,
  f : GFFFeature,
  mode : String,
) -> Bool {
  if mode == "intersection-strict" {
    start_1b >= f.start && end_1b <= f.end
  } else {
    start_1b <= f.end && end_1b >= f.start
  }
}

///|
/// Count reads per feature. `features` are the annotation features (typically
/// genes or exons, already extracted to one record per feature of interest).
/// `mode` is `"union"` (default), `"intersection-nonempty"` (both = any
/// overlap) or `"intersection-strict"` (read must be fully contained).
pub fn htseq_count(
  bam : BamFile,
  features : Array[GFFFeature],
  mode? : String = "union",
) -> HtseqCounts {
  let header = bam.header
  let counts : Map[String, Int] = Map([], capacity=0)
  let mut no_feature = 0
  let mut ambiguous = 0
  let mut not_aligned = 0

  for rec in bam.records {
    if rec.is_unmapped() {
      not_aligned = not_aligned + 1
      continue
    }
    let ref_name = rec.get_reference_name(header)
    let start_1b = rec.pos + 1
    let end_1b = rec.get_end() + 1

    let hits : Array[Int] = Array(capacity=0)
    let mut k = 0
    while k < features.length() {
      let f = features[k]
      if f.seqid == ref_name && htseq_overlap(start_1b, end_1b, f, mode) {
        hits.push(k)
      }
      k = k + 1
    }

    if hits.length() == 0 {
      no_feature = no_feature + 1
    } else if hits.length() == 1 {
      let fid = htseq_feature_id(features[hits[0]], hits[0])
      let prev = match counts.get(fid) {
        None => 0
        Some(c) => c
      }
      counts.set(fid, prev + 1)
    } else {
      ambiguous = ambiguous + 1
    }
  }

  HtseqCounts::{ counts, no_feature, ambiguous, not_aligned }
}
