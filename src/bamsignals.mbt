///|
/// bamsignals module, inspired by Bioconductor bamsignals package.
/// Provides functions for extracting and analyzing signals from BAM files.
/// Supports signal counting and normalization for ChIP-seq and other sequencing data.

/// Count mode for signal extraction
pub enum BamsigCountMode {
  /// Count all reads
  All
  /// Only count reads on specific strand
  Sense
  /// Only count reads on opposite strand
  Antisense
  /// Count paired-end reads as single events
  PairedEnd
} derive(Eq)

/// Normalization method for signal
pub enum BamsigNormMethod {
  /// No normalization
  None
  /// RPM (reads per million)
  RPM
  /// RPKM (reads per kilobase per million)
  RPKM
  /// CPM (counts per million)
  CPM
}

/// Signal extraction parameters
pub struct BamsigParams {
  /// Count mode
  count_mode : BamsigCountMode
  /// Whether to filter duplicates
  filter_dup : Bool
  /// Minimum mapping quality
  min_mapq : Int
  /// Whether to use single-end mode
  single_end : Bool
  /// Extend reads to fragment size
  extend_len : Int
}

/// Create default parameters
pub fn BamsigParams::new() -> BamsigParams {
  BamsigParams::{ count_mode: BamsigCountMode::All, filter_dup: false, min_mapq: 10, single_end: false, extend_len: 150 }
}

/// Set count mode
pub fn BamsigParams::bamsig_set_count_mode(self : BamsigParams, mode : BamsigCountMode) -> BamsigParams {
  BamsigParams::{ count_mode: mode, filter_dup: self.filter_dup, min_mapq: self.min_mapq, single_end: self.single_end, extend_len: self.extend_len }
}

/// Set filter duplicates
pub fn BamsigParams::bamsig_set_filter_dup(self : BamsigParams, filter : Bool) -> BamsigParams {
  BamsigParams::{ count_mode: self.count_mode, filter_dup: filter, min_mapq: self.min_mapq, single_end: self.single_end, extend_len: self.extend_len }
}

/// Set minimum mapping quality
pub fn BamsigParams::bamsig_set_min_mapq(self : BamsigParams, mapq : Int) -> BamsigParams {
  BamsigParams::{ count_mode: self.count_mode, filter_dup: self.filter_dup, min_mapq: mapq, single_end: self.single_end, extend_len: self.extend_len }
}

/// Set extend length
pub fn BamsigParams::bamsig_set_extend(self : BamsigParams, extend_len : Int) -> BamsigParams {
  BamsigParams::{ count_mode: self.count_mode, filter_dup: self.filter_dup, min_mapq: self.min_mapq, single_end: self.single_end, extend_len: extend_len }
}

/// Genomic region for signal extraction
pub struct BamsigRegion {
  /// Chromosome
  chrom : String
  /// Start position (0-based)
  start : Int
  /// End position (0-based, exclusive)
  end : Int
  /// Region ID
  id : String
}

/// Create new region
pub fn BamsigRegion::new(chrom : String, start : Int, end : Int, id : String) -> BamsigRegion {
  BamsigRegion::{ chrom, start, end, id }
}

/// Get chromosome
pub fn BamsigRegion::bamsig_chrom(self : BamsigRegion) -> String {
  self.chrom
}

/// Get start position
pub fn BamsigRegion::bamsig_start(self : BamsigRegion) -> Int {
  self.start
}

/// Get end position
pub fn BamsigRegion::bamsig_end(self : BamsigRegion) -> Int {
  self.end
}

/// Get width
pub fn BamsigRegion::bamsig_width(self : BamsigRegion) -> Int {
  self.end - self.start
}

/// Get ID
pub fn BamsigRegion::bamsig_id(self : BamsigRegion) -> String {
  self.id
}

/// BAM record for signal extraction
pub struct BamsigRecord {
  /// Chromosome
  chrom : String
  /// Start position
  pos : Int
  /// CIGAR string (simplified)
  cigar : String
  /// Strand (true = sense, false = antisense)
  strand : Bool
  /// Mapping quality
  mapq : Int
  /// Is duplicate
  is_dup : Bool
  /// Is paired end
  is_paired : Bool
  /// Mate pair flag
  is_first : Bool
}

/// Create new record
pub fn BamsigRecord::new(chrom : String, pos : Int, cigar : String, strand : Bool, mapq : Int, is_dup : Bool, is_paired : Bool, is_first : Bool) -> BamsigRecord {
  BamsigRecord::{ chrom, pos, cigar, strand, mapq, is_dup, is_paired, is_first }
}

/// Get alignment length from CIGAR
pub fn BamsigRecord::bamsig_align_length(self : BamsigRecord) -> Int {
  parse_cigar_length(self.cigar)
}

/// Get end position
pub fn BamsigRecord::bamsig_end(self : BamsigRecord) -> Int {
  self.pos + self.bamsig_align_length()
}

/// Check if record passes filter
pub fn BamsigRecord::bamsig_passes_filter(self : BamsigRecord, params : BamsigParams) -> Bool {
  if self.is_dup && params.filter_dup {
    return false
  }
  if self.mapq < params.min_mapq {
    return false
  }
  if params.single_end && self.is_paired {
    return false
  }
  true
}

/// Signal counts for regions
pub struct BamsigSignal {
  /// Region IDs
  region_ids : Array[String]
  /// Count matrix (regions x samples)
  counts : Array[Array[Double]]
  /// Normalized counts
  norm_counts : Array[Array[Double]]
  /// Total reads per sample
  total_reads : Array[Double]
}

/// Create new signal
pub fn BamsigSignal::new(region_ids : Array[String], counts : Array[Array[Double]], 
             norm_counts : Array[Array[Double]], total_reads : Array[Double]) -> BamsigSignal {
  BamsigSignal::{ region_ids, counts, norm_counts, total_reads }
}

/// Get number of regions
pub fn BamsigSignal::bamsig_n_regions(self : BamsigSignal) -> Int {
  self.region_ids.length()
}

/// Get number of samples
pub fn BamsigSignal::bamsig_n_samples(self : BamsigSignal) -> Int {
  if self.counts.length() > 0 { self.counts[0].length() } else { 0 }
}

/// Get count for region and sample
pub fn BamsigSignal::bamsig_get_count(self : BamsigSignal, region_idx : Int, sample_idx : Int) -> Double {
  self.counts[region_idx][sample_idx]
}

/// Get normalized count
pub fn BamsigSignal::bamsig_get_norm_count(self : BamsigSignal, region_idx : Int, sample_idx : Int) -> Double {
  self.norm_counts[region_idx][sample_idx]
}

/// Get region ID
pub fn BamsigSignal::bamsig_get_region_id(self : BamsigSignal, idx : Int) -> String {
  self.region_ids[idx]
}

/// Sum counts across samples for a region
pub fn BamsigSignal::bamsig_row_sum(self : BamsigSignal, region_idx : Int) -> Double {
  let row = self.counts[region_idx]
  let mut sum = 0.0
  let mut i = 0
  while i < row.length() {
    sum = sum + row[i]
    i = i + 1
  }
  sum
}

/// Sum counts across regions for a sample
pub fn BamsigSignal::bamsig_col_sum(self : BamsigSignal, sample_idx : Int) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < self.counts.length() {
    sum = sum + self.counts[i][sample_idx]
    i = i + 1
  }
  sum
}

/// Filter regions by count threshold
pub fn BamsigSignal::bamsig_filter(self : BamsigSignal, min_count : Double, min_samples : Int) -> BamsigSignal {
  let new_region_ids : Array[String] = Array::new()
  let new_counts : Array[Array[Double]] = Array::new()
  let new_norm_counts : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < self.region_ids.length() {
    let count_row = self.counts[i]
    let mut passing_samples = 0
    
    let mut j = 0
    while j < count_row.length() {
      if count_row[j] >= min_count {
        passing_samples = passing_samples + 1
      }
      j = j + 1
    }
    
    if passing_samples >= min_samples {
      new_region_ids.push(self.region_ids[i])
      new_counts.push(count_row)
      new_norm_counts.push(self.norm_counts[i])
    }
    
    i = i + 1
  }
  
  BamsigSignal::new(new_region_ids, new_counts, new_norm_counts, self.total_reads)
}

/// Chromatin state analysis result
pub struct BamsigChromState {
  /// State labels
  labels : Array[String]
  /// Colors for each state
  colors : Array[String]
  /// State frequencies
  frequencies : Array[Double]
  /// Active promoter percentage
  active_promoter_pct : Double
  /// Active enhancer percentage
  active_enhancer_pct : Double
  /// Repressed percentage
  repressed_pct : Double
}

/// Create new chromatin state result
pub fn BamsigChromState::new(labels : Array[String], colors : Array[String], frequencies : Array[Double],
             active_promoter_pct : Double, active_enhancer_pct : Double, 
             repressed_pct : Double) -> BamsigChromState {
  BamsigChromState::{ labels, colors, frequencies, active_promoter_pct, active_enhancer_pct, repressed_pct }
}

/// Get number of states
pub fn BamsigChromState::bamsig_n_states(self : BamsigChromState) -> Int {
  self.labels.length()
}

/// Get state label
pub fn BamsigChromState::bamsig_get_label(self : BamsigChromState, idx : Int) -> String {
  self.labels[idx]
}

/// Get state frequency
pub fn BamsigChromState::bamsig_get_frequency(self : BamsigChromState, idx : Int) -> Double {
  self.frequencies[idx]
}

/// Get active promoter percentage
pub fn BamsigChromState::bamsig_active_promoter(self : BamsigChromState) -> Double {
  self.active_promoter_pct
}

/// Get active enhancer percentage
pub fn BamsigChromState::bamsig_active_enhancer(self : BamsigChromState) -> Double {
  self.active_enhancer_pct
}

/// Get repressed percentage
pub fn BamsigChromState::bamsig_repressed(self : BamsigChromState) -> Double {
  self.repressed_pct
}

/// Parse CIGAR string to get alignment length
fn parse_cigar_length(cigar : String) -> Int {
  let mut total = 0
  let mut i = 0
  while i < cigar.length() {
    let mut num_str = ""
    while i < cigar.length() && cigar[i] >= '0' && cigar[i] <= '9' {
      num_str = num_str + cigar[i].to_string()
      i = i + 1
    }
    if num_str.length() > 0 && i < cigar.length() {
      let mut num = 0
      let mut k = 0
      while k < num_str.length() {
        let digit = num_str[k].to_int() - '0'.to_int()
        num = num * 10 + digit
        k = k + 1
      }
      let op = cigar[i]
      if op == 'M' || op == 'X' || op == '=' {
        total = total + num
      }
      i = i + 1
    }
  }
  total
}

/// Create example records for testing
pub fn bamsig_create_example_records() -> Array[BamsigRecord] {
  [
    BamsigRecord::new("chr1", 100, "50M", true, 30, false, false, true),
    BamsigRecord::new("chr1", 200, "45M", false, 25, false, false, true),
    BamsigRecord::new("chr1", 150, "60M", true, 35, false, false, true),
    BamsigRecord::new("chr1", 300, "55M", true, 20, true, false, true),
    BamsigRecord::new("chr2", 100, "40M", false, 28, false, false, true),
    BamsigRecord::new("chr2", 150, "48M", true, 32, false, false, true),
  ]
}

/// Create example regions
pub fn bamsig_create_example_regions() -> Array[BamsigRegion] {
  [
    BamsigRegion::new("chr1", 100, 200, "Region1"),
    BamsigRegion::new("chr1", 150, 250, "Region2"),
    BamsigRegion::new("chr1", 300, 400, "Region3"),
    BamsigRegion::new("chr2", 100, 200, "Region4"),
  ]
}

/// Count signals in regions
pub fn bamsig_count_signals(regions : Array[BamsigRegion], 
                            records_by_chrom : Map[String, Array[BamsigRecord]],
                            params : BamsigParams) -> Array[Array[Double]] {
  let n_regions = regions.length()
  let counts : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_regions {
    let region = regions[i]
    let mut region_count = 0.0
    
    let chrom_records = records_by_chrom.get(region.chrom)
    if chrom_records is Some(records) {
      let recs = records
      let mut j = 0
      while j < recs.length() {
        let rec = recs[j]
        if rec.bamsig_passes_filter(params) {
          let rec_end = rec.pos + rec.bamsig_align_length()
          if rec.pos < region.end && rec_end > region.start {
            let overlap_start = if rec.pos > region.start { rec.pos } else { region.start }
            let overlap_end = if rec_end < region.end { rec_end } else { region.end }
            let _overlap = if overlap_end > overlap_start { overlap_end - overlap_start } else { 0 }
            
            if params.count_mode == BamsigCountMode::All {
              region_count = region_count + 1.0
            } else if params.count_mode == BamsigCountMode::Sense && rec.strand {
              region_count = region_count + 1.0
            } else if params.count_mode == BamsigCountMode::Antisense && !rec.strand {
              region_count = region_count + 1.0
            }
          }
        }
        j = j + 1
      }
    }
    
    counts.push([region_count])
    i = i + 1
  }
  
  counts
}

/// Normalize signals
pub fn bamsig_normalize_signals(counts : Array[Array[Double]], 
                                 total_reads : Array[Double],
                                 norm_method : BamsigNormMethod) -> Array[Array[Double]] {
  let norm_counts : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < counts.length() {
    let row = counts[i]
    let norm_row : Array[Double] = Array::new()
    
    let mut j = 0
    while j < row.length() {
      let c = row[j]
      let total = if j < total_reads.length() { total_reads[j] } else { 1.0 }
      
      let norm_val = match norm_method {
        BamsigNormMethod::None => c
        BamsigNormMethod::RPM => if total > 0.0 { c / total * 1000000.0 } else { 0.0 }
        BamsigNormMethod::CPM => if total > 0.0 { c / total * 1000000.0 } else { 0.0 }
        BamsigNormMethod::RPKM => if total > 0.0 { c / total * 1000000.0 / 1000.0 } else { 0.0 }
      }
      
      norm_row.push(norm_val)
      j = j + 1
    }
    
    norm_counts.push(norm_row)
    i = i + 1
  }
  
  norm_counts
}

/// Analyze chromatin states
pub fn bamsig_analyze_chromatin_states(signal : BamsigSignal, 
                                         promoter_regions : Array[Int],
                                         enhancer_regions : Array[Int],
                                         repressed_regions : Array[Int]) -> BamsigChromState {
  let n_regions = signal.bamsig_n_regions()
  let n_states = 4
  
  let labels = ["Active Promoter", "Active Enhancer", "Repressed", "Inactive"]
  let colors = ["#FF0000", "#00FF00", "#0000FF", "#CCCCCC"]
  let frequencies = Array::make(n_states, 0.0)
  
  let n_promoter = promoter_regions.length().to_double()
  let n_enhancer = enhancer_regions.length().to_double()
  let n_repressed = repressed_regions.length().to_double()
  
  let mut active_promoter = 0.0
  let mut active_enhancer = 0.0
  let mut repressed = 0.0
  
  let mut i = 0
  while i < n_regions {
    let total_count = signal.bamsig_row_sum(i)
    let is_active = total_count > 0.0
    
    let is_promoter = promoter_regions.contains(i)
    let is_enhancer = enhancer_regions.contains(i)
    let is_repressed = repressed_regions.contains(i)
    
    if is_active {
      if is_promoter {
        active_promoter = active_promoter + 1.0
      } else if is_enhancer {
        active_enhancer = active_enhancer + 1.0
      }
    }
    
    if is_repressed {
      repressed = repressed + 1.0
    }
    
    i = i + 1
  }
  
  let _total_annotated = n_promoter + n_enhancer + n_repressed
  let total_promoter_pct = if n_promoter > 0.0 { active_promoter / n_promoter * 100.0 } else { 0.0 }
  let enhancer_pct = if n_enhancer > 0.0 { active_enhancer / n_enhancer * 100.0 } else { 0.0 }
  let repressed_pct = if n_repressed > 0.0 { repressed / n_repressed * 100.0 } else { 0.0 }
  
  frequencies[0] = total_promoter_pct
  frequencies[1] = enhancer_pct
  frequencies[2] = repressed_pct
  frequencies[3] = 100.0 - total_promoter_pct - enhancer_pct - repressed_pct
  
  BamsigChromState::new(labels, colors, frequencies, 
                         total_promoter_pct, enhancer_pct, repressed_pct)
}

/// Create example signal data
pub fn bamsig_create_example_signal() -> BamsigSignal {
  let regions = bamsig_create_example_regions()
  let n_regions = regions.length()
  let _n_samples = 2
  
  let region_ids : Array[String] = Array::new()
  let counts : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_regions {
    region_ids.push(regions[i].id)
    let row : Array[Double] = Array::new()
    row.push((i + 1).to_double() * 50.0)
    row.push((i + 1).to_double() * 45.0)
    counts.push(row)
    i = i + 1
  }
  
  let total_reads = [10000.0, 9500.0]
  let norm_counts = bamsig_normalize_signals(counts, total_reads, BamsigNormMethod::RPM)
  
  BamsigSignal::new(region_ids, counts, norm_counts, total_reads)
}
