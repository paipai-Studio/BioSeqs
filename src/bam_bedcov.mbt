/// Coverage statistics from BED file regions.
///
/// pysam/samtools equivalent: `samtools bedcov`

///|
/// Coverage result for a single BED region
pub struct BedCovResult {
  chrom : String
  start : Int32
  end : Int32
  total_coverage : Int64
  mean_depth : Double
  region_length : Int32
}

///|
pub fn BedCovResult::chrom(self : BedCovResult) -> String {
  self.chrom
}

///|
pub fn BedCovResult::start(self : BedCovResult) -> Int32 {
  self.start
}

///|
pub fn BedCovResult::end(self : BedCovResult) -> Int32 {
  self.end
}

///|
pub fn BedCovResult::total_coverage(self : BedCovResult) -> Int64 {
  self.total_coverage
}

///|
pub fn BedCovResult::mean_depth(self : BedCovResult) -> Double {
  self.mean_depth
}

///|
pub fn BedCovResult::region_length(self : BedCovResult) -> Int32 {
  self.region_length
}

///|
/// A simple BED region
pub struct BedRegion {
  chrom : String
  start : Int32
  end : Int32
}

///|
/// Parse a BED file (3-column: chrom, start, end) into regions.
pub fn bam_parse_bed(content : String) -> Array[BedRegion] {
  let regions : Array[BedRegion] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 || trimmed[0:1].to_owned() == "#" {
      continue
    }
    let fields = trimmed.split("\t").to_array()
    if fields.length() >= 3 {
      let chrom = fields[0].to_owned()
      let start = parse_int(fields[1].to_owned())
      let end = parse_int(fields[2].to_owned())
      regions.push({ chrom, start, end, })
    }
  }
  regions
}

///|
/// Compute coverage for each BED region from a BamFile.
pub fn bam_bedcov(
  bam : BamFile,
  regions : Array[BedRegion],
) -> Array[BedCovResult] {
  let results : Array[BedCovResult] = Array(capacity=0)
  for region in regions {
    let mut total : Int64 = 0
    let reg_len = region.end - region.start
    if reg_len <= 0 {
      results.push({
        chrom: region.chrom,
        start: region.start,
        end: region.end,
        total_coverage: 0,
        mean_depth: 0.0,
        region_length: 0,
      })
      continue
    }
    for rec in bam.records {
      if rec.is_unmapped() {
        continue
      }
      let ref_name = rec.get_reference_name(bam.header)
      if ref_name != region.chrom {
        continue
      }
      let rec_start = rec.pos
      let rec_end = rec.get_end() + 1
      let ov_start = if rec_start > region.start {
        rec_start
      } else {
        region.start
      }
      let ov_end = if rec_end < region.end { rec_end } else { region.end }
      if ov_end > ov_start {
        total = total + (ov_end - ov_start).to_int64()
      }
    }
    let mean = total.to_double() / reg_len.to_double()
    results.push({
      chrom: region.chrom,
      start: region.start,
      end: region.end,
      total_coverage: total,
      mean_depth: mean,
      region_length: reg_len,
    })
  }
  results
}

///|
/// Convenience: parse BED content and compute coverage in one call.
pub fn bam_bedcov_from_bed(
  bam : BamFile,
  bed_content : String,
) -> Array[BedCovResult] {
  let regions = bam_parse_bed(bed_content)
  bam_bedcov(bam, regions)
}
