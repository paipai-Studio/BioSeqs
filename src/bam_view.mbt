/// Filter BAM records by flag, mapping quality, region, and read group.
///
/// pysam/samtools equivalent: `samtools view -f -F -q -r`
/// pysam.AlignmentFile.fetch with flag filters

///|
/// Filter parameters for bam_view
pub struct BamViewFilter {
  /// Require these flags (AND with flag, must equal required_flags)
  require_flags : Int32
  /// Exclude records with any of these flags set
  exclude_flags : Int32
  /// Minimum mapping quality (inclusive)
  min_mapq : Int32
  /// Only include reads mapping to this reference name (None = all)
  contig : String?
  /// Region start (0-based, inclusive), only used when contig is set
  region_start : Int32?
  /// Region end (0-based, exclusive), only used when contig is set
  region_end : Int32?
  /// Only include reads with this read group (None = all)
  read_group : String?
}

///|
pub fn BamViewFilter::new() -> BamViewFilter {
  {
    require_flags: 0,
    exclude_flags: 0,
    min_mapq: 0,
    contig: None,
    region_start: None,
    region_end: None,
    read_group: None,
  }
}

///|
pub fn BamViewFilter::require(
  self : BamViewFilter,
  flags : Int32,
) -> BamViewFilter {
  { ..self, require_flags: self.require_flags.lor(flags), }
}

///|
pub fn BamViewFilter::exclude(
  self : BamViewFilter,
  flags : Int32,
) -> BamViewFilter {
  { ..self, exclude_flags: self.exclude_flags.lor(flags), }
}

///|
pub fn BamViewFilter::min_mapq(
  self : BamViewFilter,
  q : Int32,
) -> BamViewFilter {
  { ..self, min_mapq: q, }
}

///|
pub fn BamViewFilter::contig(
  self : BamViewFilter,
  name : String,
) -> BamViewFilter {
  { ..self, contig: Some(name), }
}

///|
pub fn BamViewFilter::region(
  self : BamViewFilter,
  start : Int32,
  end : Int32,
) -> BamViewFilter {
  { ..self, region_start: Some(start), region_end: Some(end), }
}

///|
pub fn BamViewFilter::read_group(
  self : BamViewFilter,
  rg : String,
) -> BamViewFilter {
  { ..self, read_group: Some(rg), }
}

///|
/// Check if a single record passes the filter.
pub fn bam_view_passes(
  rec : BamRecord,
  filter : BamViewFilter,
  header : BamHeader,
) -> Bool {
  // Flag filter: required flags must be set, excluded flags must not
  if rec.flag.land(filter.require_flags) != filter.require_flags {
    return false
  }
  if rec.flag.land(filter.exclude_flags) != 0 {
    return false
  }
  // Mapping quality
  if rec.mapq < filter.min_mapq {
    return false
  }
  // Contig / region filter
  match filter.contig {
    Some(ctg) => {
      let ref_name = rec.get_reference_name(header)
      if ref_name != ctg {
        return false
      }
      // Region overlap check
      match (filter.region_start, filter.region_end) {
        (Some(rs), Some(re)) => {
          let rec_end = rec.get_end() + 1
          if rec.pos >= re || rec_end <= rs {
            return false
          }
        }
        _ => ()
      }
    }
    None => ()
  }
  // Read group filter
  match filter.read_group {
    Some(rg) =>
      match rec.get_tag("RG") {
        Some(BamTagValue::Str(s)) => if s != rg { return false }
        _ => return false
      }
    None => ()
  }
  true
}

///|
/// Filter BAM records, returning a new BamFile with only passing records.
pub fn bam_view(bam : BamFile, filter : BamViewFilter) -> BamFile {
  let passing : Array[BamRecord] = []
  for rec in bam.records {
    if bam_view_passes(rec, filter, bam.header) {
      passing.push(rec)
    }
  }
  { header: bam.header, records: passing, }
}

///|
/// Count records passing the filter without building a new BamFile.
pub fn bam_view_count(bam : BamFile, filter : BamViewFilter) -> Int32 {
  let mut count : Int32 = 0
  for rec in bam.records {
    if bam_view_passes(rec, filter, bam.header) {
      count = count + 1
    }
  }
  count
}
