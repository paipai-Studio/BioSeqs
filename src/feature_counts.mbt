///|
/// featureCounts module - Read counting over genomic features.
/// Inspired by Bioconductor Rsubread/featureCounts package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/Rsubread.html
///   - Liao Y et al. (2014) featureCounts: an efficient general purpose
///     program for assigning sequence reads to genomic features.
///     Bioinformatics 30(7):923-930.

///|
/// Strand-specificity mode for counting.
pub enum StrandMode {
  Unstranded
  Stranded
  Reversed
} derive(Eq)

///|
/// Create unstranded mode.
pub fn strand_unstranded() -> StrandMode {
  StrandMode::Unstranded
}

///|
/// Create stranded mode.
pub fn strand_stranded() -> StrandMode {
  StrandMode::Stranded
}

///|
/// Create reversely stranded mode.
pub fn strand_reversed() -> StrandMode {
  StrandMode::Reversed
}

///|
/// Feature annotation record for counting.
pub struct FeatureAnnotation {
  /// Chromosome name
  chr : String
  /// Start position (1-based)
  start : Int
  /// End position (1-based, inclusive)
  end_ : Int
  /// Strand ('+', '-', or '*')
  strand : String
  /// Gene ID
  gene_id : String
  /// Transcript ID (optional)
  transcript_id : String
  /// Feature type (exon, gene, transcript, etc.)
  feature_type : String
  /// Unique feature identifier
  feature_id : String
}

///|
/// Create a new FeatureAnnotation.
pub fn FeatureAnnotation::new(
  chr~ : String,
  start~ : Int,
  end_~ : Int,
  strand~ : String,
  gene_id~ : String,
  feature_type? : String = "exon",
  transcript_id? : String = "",
  feature_id? : String = "",
) -> FeatureAnnotation {
  let fid = if feature_id == "" {
    gene_id + ":" + start.to_string() + "-" + end_.to_string()
  } else {
    feature_id
  }
  FeatureAnnotation::{
    chr, start, end_, strand, gene_id,
    transcript_id, feature_type,
    feature_id: fid,
  }
}

///|
/// Get the length of the feature.
pub fn FeatureAnnotation::length(self : FeatureAnnotation) -> Int {
  if self.end_ > self.start { self.end_ - self.start + 1 } else { 0 }
}

///|
/// Check if this feature overlaps with a given range.
pub fn FeatureAnnotation::overlaps(
  self : FeatureAnnotation,
  chr~ : String,
  start~ : Int,
  end_~ : Int,
) -> Bool {
  if self.chr != chr { return false }
  self.start <= end_ && self.end_ >= start
}

///|
/// Get the overlap length between this feature and a range.
pub fn FeatureAnnotation::overlap_length(
  self : FeatureAnnotation,
  chr~ : String,
  start~ : Int,
  end_~ : Int,
) -> Int {
  if self.chr != chr { return 0 }
  let ov_start = if self.start > start { self.start } else { start }
  let ov_end = if self.end_ < end_ { self.end_ } else { end_ }
  if ov_end >= ov_start { ov_end - ov_start + 1 } else { 0 }
}

///|
/// Read alignment record for counting.
pub struct ReadAlignment {
  /// Read ID
  read_id : String
  /// Chromosome
  chr : String
  /// Start position (1-based)
  start : Int
  /// End position (1-based, inclusive)
  end_ : Int
  /// Strand
  strand : String
  /// Mapping quality
  mapq : Int
  /// Number of alignments (for multi-mappers)
  n_alignments : Int
  /// Is paired-end read
  is_paired : Bool
  /// Mate chromosome (for paired-end)
  mate_chr : String
  /// Mate start (for paired-end)
  mate_start : Int
  /// Fragment length (for paired-end)
  fragment_length : Int
}

///|
/// Create a new ReadAlignment.
pub fn ReadAlignment::new(
  read_id~ : String,
  chr~ : String,
  start~ : Int,
  end_~ : Int,
  strand? : String = "+",
  mapq? : Int = 60,
) -> ReadAlignment {
  ReadAlignment::{
    read_id, chr, start, end_, strand, mapq,
    n_alignments: 1,
    is_paired: false,
    mate_chr: "",
    mate_start: 0,
    fragment_length: 0,
  }
}

///|
/// Create a paired-end ReadAlignment.
pub fn ReadAlignment::new_paired(
  read_id~ : String,
  chr~ : String,
  start~ : Int,
  end_~ : Int,
  strand? : String = "+",
  mapq? : Int = 60,
  mate_chr? : String = "",
  mate_start? : Int = 0,
) -> ReadAlignment {
  let frag_len = if mate_start > 0 {
    if mate_start > end_ { mate_start - start + 1 } else { end_ - mate_start + 1 }
  } else {
    0
  }
  ReadAlignment::{
    read_id, chr, start, end_, strand, mapq,
    n_alignments: 1,
    is_paired: true,
    mate_chr: if mate_chr == "" { chr } else { mate_chr },
    mate_start,
    fragment_length: frag_len,
  }
}

///|
/// Get the fragment span for paired-end reads.
pub fn ReadAlignment::fragment_span(self : ReadAlignment) -> (String, Int, Int) {
  if !self.is_paired {
    return (self.chr, self.start, self.end_)
  }
  if self.chr == self.mate_chr {
    let f_start = if self.start < self.mate_start { self.start } else { self.mate_start }
    let f_end = if self.end_ > self.mate_start { self.end_ } else { self.mate_start }
    (self.chr, f_start, f_end)
  } else {
    (self.chr, self.start, self.end_)
  }
}

///|
/// Configuration for featureCounts.
pub struct FeatureCountsConfig {
  /// Minimum overlap required for counting (bp)
  min_overlap : Int
  /// Minimum total overlap (bp) for multi-exon features
  min_overlap_total : Int
  /// Fraction of read that must overlap
  frac_overlap : Double
  /// Minimum mapping quality
  min_mapq : Int
  /// Ignore strand specificity
  ignore_strand : Bool
  /// Strand mode
  strand_mode : StrandMode
  /// Count multi-mapping reads multiple times
  multi_count : Bool
  /// Fractional counting for multi-mappers
  frac_multi_count : Bool
  /// Include chimeric reads
  include_chimeric : Bool
  /// Count fragments (paired-end) or reads
  count_fragments : Bool
}

///|
/// Create default FeatureCountsConfig.
pub fn FeatureCountsConfig::new() -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: 1,
    min_overlap_total: 10,
    frac_overlap: 0.0,
    min_mapq: 0,
    ignore_strand: false,
    strand_mode: StrandMode::Unstranded,
    multi_count: false,
    frac_multi_count: false,
    include_chimeric: false,
    count_fragments: false,
  }
}

///|
/// Set minimum overlap (returns new config).
pub fn FeatureCountsConfig::set_min_overlap(self : FeatureCountsConfig, val~ : Int) -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: val,
    min_overlap_total: self.min_overlap_total,
    frac_overlap: self.frac_overlap,
    min_mapq: self.min_mapq,
    ignore_strand: self.ignore_strand,
    strand_mode: self.strand_mode,
    multi_count: self.multi_count,
    frac_multi_count: self.frac_multi_count,
    include_chimeric: self.include_chimeric,
    count_fragments: self.count_fragments,
  }
}

///|
/// Set strand mode (returns new config).
pub fn FeatureCountsConfig::set_strand_mode(self : FeatureCountsConfig, mode~ : StrandMode) -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: self.min_overlap,
    min_overlap_total: self.min_overlap_total,
    frac_overlap: self.frac_overlap,
    min_mapq: self.min_mapq,
    ignore_strand: self.ignore_strand,
    strand_mode: mode,
    multi_count: self.multi_count,
    frac_multi_count: self.frac_multi_count,
    include_chimeric: self.include_chimeric,
    count_fragments: self.count_fragments,
  }
}

///|
/// Set minimum mapping quality (returns new config).
pub fn FeatureCountsConfig::set_min_mapq(self : FeatureCountsConfig, val~ : Int) -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: self.min_overlap,
    min_overlap_total: self.min_overlap_total,
    frac_overlap: self.frac_overlap,
    min_mapq: val,
    ignore_strand: self.ignore_strand,
    strand_mode: self.strand_mode,
    multi_count: self.multi_count,
    frac_multi_count: self.frac_multi_count,
    include_chimeric: self.include_chimeric,
    count_fragments: self.count_fragments,
  }
}

///|
/// Set count fragments flag (returns new config).
pub fn FeatureCountsConfig::set_count_fragments(self : FeatureCountsConfig, val~ : Bool) -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: self.min_overlap,
    min_overlap_total: self.min_overlap_total,
    frac_overlap: self.frac_overlap,
    min_mapq: self.min_mapq,
    ignore_strand: self.ignore_strand,
    strand_mode: self.strand_mode,
    multi_count: self.multi_count,
    frac_multi_count: self.frac_multi_count,
    include_chimeric: self.include_chimeric,
    count_fragments: val,
  }
}

///|
/// Set multi-count flag (returns new config).
pub fn fcfg_set_multi_count(cfg : FeatureCountsConfig, val~ : Bool) -> FeatureCountsConfig {
  FeatureCountsConfig::{
    min_overlap: cfg.min_overlap,
    min_overlap_total: cfg.min_overlap_total,
    frac_overlap: cfg.frac_overlap,
    min_mapq: cfg.min_mapq,
    ignore_strand: cfg.ignore_strand,
    strand_mode: cfg.strand_mode,
    multi_count: val,
    frac_multi_count: cfg.frac_multi_count,
    include_chimeric: cfg.include_chimeric,
    count_fragments: cfg.count_fragments,
  }
}

///|
/// Feature counting result.
pub struct FeatureCountsResult {
  /// Feature IDs
  feature_ids : Array[String]
  /// Gene IDs
  gene_ids : Array[String]
  /// Per-sample count matrix (sample x feature)
  counts : Array[Array[Double]]
  /// Sample names
  sample_names : Array[String]
  /// Assignment summary
  total_reads : Array[Int]
  /// Assigned reads
  assigned_reads : Array[Int]
  /// Unassigned reads (no feature overlap)
  unassigned_no_features : Array[Int]
  /// Unassigned reads (ambiguous)
  unassigned_ambiguous : Array[Int]
  /// Unassigned reads (multi-mapping)
  unassigned_multi : Array[Int]
  /// Unassigned reads (low quality)
  unassigned_low_quality : Array[Int]
  /// Unassigned reads (overlapping features)
  unassigned_overlap : Array[Int]
}

///|
/// Get number of features.
pub fn FeatureCountsResult::n_features(self : FeatureCountsResult) -> Int {
  self.feature_ids.length()
}

///|
/// Get number of samples.
pub fn FeatureCountsResult::n_samples(self : FeatureCountsResult) -> Int {
  self.sample_names.length()
}

///|
/// Get count for a specific feature and sample.
pub fn FeatureCountsResult::get_count(
  self : FeatureCountsResult,
  feature_id~ : String,
  sample_idx~ : Int,
) -> Double {
  let fi = fc_feature_index(self, feature_id)
  if fi < 0 { return 0.0 }
  if sample_idx < 0 || sample_idx >= self.counts.length() { return 0.0 }
  self.counts[sample_idx][fi]
}

///|
/// Get library size (total assigned counts) for a sample.
pub fn FeatureCountsResult::library_size(
  self : FeatureCountsResult,
  sample_idx~ : Int,
) -> Double {
  if sample_idx < 0 || sample_idx >= self.counts.length() { return 0.0 }
  let mut total = 0.0
  let mut i = 0
  while i < self.counts[sample_idx].length() {
    total = total + self.counts[sample_idx][i]
    i = i + 1
  }
  total
}

///|
/// Get CPM (counts per million) normalized count.
pub fn FeatureCountsResult::cpm(
  self : FeatureCountsResult,
  feature_id~ : String,
  sample_idx~ : Int,
) -> Double {
  let raw = self.get_count(feature_id=feature_id, sample_idx=sample_idx)
  let lib = self.library_size(sample_idx=sample_idx)
  if lib == 0.0 { 0.0 } else { raw / lib * 1_000_000.0 }
}

///|
/// Get total count for a feature across all samples.
pub fn FeatureCountsResult::get_feature_total(
  self : FeatureCountsResult,
  feature_id~ : String,
) -> Double {
  let fi = fc_feature_index(self, feature_id)
  if fi < 0 { return 0.0 }
  let mut total = 0.0
  let mut i = 0
  while i < self.counts.length() {
    total = total + self.counts[i][fi]
    i = i + 1
  }
  total
}

///|
/// Summary of counting results.
pub fn FeatureCountsResult::summary(self : FeatureCountsResult) -> String {
  let mut s = "=== featureCounts Summary ===" + "\n"
  s = s + "Features: " + self.feature_ids.length().to_string() + "\n"
  s = s + "Samples: " + self.sample_names.length().to_string() + "\n"
  let mut i = 0
  while i < self.sample_names.length() {
    let name = self.sample_names[i]
    let total = self.total_reads[i]
    let assigned = self.assigned_reads[i]
    let pct = if total > 0 {
      (assigned.to_double() / total.to_double() * 100.0).to_string() + "%"
    } else { "0.0%" }
    s = s + "  " + name + ": total=" + total.to_string() + ", assigned=" + assigned.to_string() + " (" + pct + ")" + "\n"
    s = s + "    Unassigned: no_features=" + self.unassigned_no_features[i].to_string()
    s = s + ", ambiguous=" + self.unassigned_ambiguous[i].to_string()
    s = s + ", multi=" + self.unassigned_multi[i].to_string()
    s = s + ", low_quality=" + self.unassigned_low_quality[i].to_string() + "\n"
    i = i + 1
  }
  s
}

///|
/// Convert results to string with first 10 features.
pub fn FeatureCountsResult::to_string(self : FeatureCountsResult) -> String {
  let mut s = self.summary()
  s = s + "\n=== Count Matrix (first 10 features) ===" + "\n"
  let n_features = if self.feature_ids.length() < 10 { self.feature_ids.length() } else { 10 }
  let mut i = 0
  while i < n_features {
    let fid = self.feature_ids[i]
    let mut row = "  " + fid
    let mut j = 0
    while j < self.sample_names.length() {
      row = row + " " + self.counts[j][i].to_string()
      j = j + 1
    }
    s = s + row + "\n"
    i = i + 1
  }
  s
}

///|
/// Main feature counting function.
pub fn feature_counts_count(
  features : Array[FeatureAnnotation],
  reads : Array[ReadAlignment],
  sample_names? : Array[String] = [],
  config? : FeatureCountsConfig = FeatureCountsConfig::new(),
) -> FeatureCountsResult {
  let n_features = features.length()
  let n_reads = reads.length()
  let n_samples = if sample_names.length() > 0 { sample_names.length() } else { 1 }

  let feat_ids : Array[String] = Array::new()
  let gene_ids : Array[String] = Array::new()
  let mut i = 0
  while i < n_features {
    feat_ids.push(features[i].feature_id)
    gene_ids.push(features[i].gene_id)
    i = i + 1
  }

  let counts : Array[Array[Double]] = Array::new()
  let mut s = 0
  while s < n_samples {
    let sample_counts : Array[Double] = Array::new()
    let mut f = 0
    while f < n_features {
      sample_counts.push(0.0)
      f = f + 1
    }
    counts.push(sample_counts)
    s = s + 1
  }

  let total_reads : Array[Int] = Array::new()
  let assigned_reads : Array[Int] = Array::new()
  let unassigned_no_features : Array[Int] = Array::new()
  let unassigned_ambiguous : Array[Int] = Array::new()
  let unassigned_multi : Array[Int] = Array::new()
  let unassigned_low_quality : Array[Int] = Array::new()
  let unassigned_overlap : Array[Int] = Array::new()

  let mut si = 0
  while si < n_samples {
    total_reads.push(0)
    assigned_reads.push(0)
    unassigned_no_features.push(0)
    unassigned_ambiguous.push(0)
    unassigned_multi.push(0)
    unassigned_low_quality.push(0)
    unassigned_overlap.push(0)
    si = si + 1
  }

  let sn = if sample_names.length() > 0 {
    sample_names
  } else {
    let def : Array[String] = Array::new()
    def.push("sample_1")
    def
  }

  let mut r = 0
  while r < n_reads {
    let read = reads[r]
    let sample_idx = r % n_samples

    total_reads[sample_idx] = total_reads[sample_idx] + 1

    if read.mapq < config.min_mapq {
      unassigned_low_quality[sample_idx] = unassigned_low_quality[sample_idx] + 1
      r = r + 1
      continue
    }

    let (f_chr, f_start, f_end) = if config.count_fragments && read.is_paired {
      read.fragment_span()
    } else {
      (read.chr, read.start, read.end_)
    }

    let read_len = f_end - f_start + 1

    let overlapping : Array[(Int, Int)] = Array::new()
    let mut fi = 0
    while fi < n_features {
      let feat = features[fi]
      if feat.chr != f_chr { fi = fi + 1; continue }

      let ov_len = feat.overlap_length(chr=f_chr, start=f_start, end_=f_end)
      if ov_len < config.min_overlap { fi = fi + 1; continue }

      if config.frac_overlap > 0.0 {
        let frac = ov_len.to_double() / read_len.to_double()
        if frac < config.frac_overlap { fi = fi + 1; continue }
      }

      if !config.ignore_strand {
        let strand_match = fc_check_strand_match(feat.strand, read.strand, config.strand_mode)
        if !strand_match { fi = fi + 1; continue }
      }

      overlapping.push((fi, ov_len))
      fi = fi + 1
    }

    if overlapping.length() == 0 {
      unassigned_no_features[sample_idx] = unassigned_no_features[sample_idx] + 1
      r = r + 1
      continue
    }

    if overlapping.length() > 1 {
      unassigned_ambiguous[sample_idx] = unassigned_ambiguous[sample_idx] + 1
      r = r + 1
      continue
    }

    let (best_fi, _) = overlapping[0]

    if read.n_alignments > 1 && !config.multi_count {
      unassigned_multi[sample_idx] = unassigned_multi[sample_idx] + 1
      r = r + 1
      continue
    }

    let count_val = if config.frac_multi_count && read.n_alignments > 1 {
      1.0 / read.n_alignments.to_double()
    } else {
      1.0
    }

    counts[sample_idx][best_fi] = counts[sample_idx][best_fi] + count_val
    assigned_reads[sample_idx] = assigned_reads[sample_idx] + 1
    r = r + 1
  }

  FeatureCountsResult::{
    feature_ids: feat_ids,
    gene_ids: gene_ids,
    counts,
    sample_names: sn,
    total_reads,
    assigned_reads,
    unassigned_no_features,
    unassigned_ambiguous,
    unassigned_multi,
    unassigned_low_quality,
    unassigned_overlap,
  }
}

///|
/// Create sample data for testing.
pub fn feature_counts_sample_data() -> (Array[FeatureAnnotation], Array[ReadAlignment], Array[String]) {
  let features : Array[FeatureAnnotation] = Array::new()
  features.push(FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+", gene_id="GeneA",
    feature_id="exon_A1",
  ))
  features.push(FeatureAnnotation::new(
    chr="chr1", start=3000, end_=5000, strand="+", gene_id="GeneA",
    feature_id="exon_A2",
  ))
  features.push(FeatureAnnotation::new(
    chr="chr1", start=6000, end_=8000, strand="-", gene_id="GeneB",
    feature_id="exon_B1",
  ))
  features.push(FeatureAnnotation::new(
    chr="chr2", start=1000, end_=3000, strand="+", gene_id="GeneC",
    feature_id="exon_C1",
  ))

  let reads : Array[ReadAlignment] = Array::new()
  reads.push(ReadAlignment::new(read_id="read1", chr="chr1", start=1500, end_=1600))
  reads.push(ReadAlignment::new(read_id="read2", chr="chr1", start=3500, end_=3600))
  reads.push(ReadAlignment::new(read_id="read3", chr="chr1", start=6500, end_=6600, strand="-"))
  reads.push(ReadAlignment::new(read_id="read4", chr="chr2", start=1500, end_=1600))
  reads.push(ReadAlignment::new(read_id="read5", chr="chr1", start=1500, end_=1600, mapq=30))
  reads.push(ReadAlignment::new(read_id="read6", chr="chr3", start=1000, end_=1100))
  reads.push(ReadAlignment::new(read_id="read7", chr="chr1", start=1700, end_=1800))

  let sample_names : Array[String] = Array::new()
  sample_names.push("sample_1")

  (features, reads, sample_names)
}

// --- Internal helpers ---

fn fc_feature_index(result : FeatureCountsResult, feature_id : String) -> Int {
  let mut i = 0
  while i < result.feature_ids.length() {
    if result.feature_ids[i] == feature_id { return i }
    i = i + 1
  }
  -1
}

fn fc_check_strand_match(feature_strand : String, read_strand : String, mode : StrandMode) -> Bool {
  if feature_strand == "*" { return true }
  match mode {
    StrandMode::Unstranded => true
    StrandMode::Stranded => feature_strand == read_strand
    StrandMode::Reversed => {
      if feature_strand == "+" { read_strand == "-" }
      else { read_strand == "+" }
    }
  }
}
