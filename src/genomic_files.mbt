///|
/// GenomicFiles module - distributed processing of large genomic files.
/// Inspired by Bioconductor GenomicFiles package.
/// Supports querying BAM, BED, VCF and BigWig files by chromosome/range.

///|
/// File type enumeration for supported genomic file formats.
pub enum FileType {
  BAM
  BED
  VCF
  BigWig
  Unknown
}

///|
/// Convert FileType to string representation.
pub fn FileType::to_string(self : FileType) -> String {
  match self {
    FileType::BAM => "BAM"
    FileType::BED => "BED"
    FileType::VCF => "VCF"
    FileType::BigWig => "BigWig"
    FileType::Unknown => "Unknown"
  }
}

///|
/// GenomicFile - represents a single genomic file with metadata.
/// The file_path field can hold either a file path or inline text content
/// (used for testing without real file I/O).
pub struct GenomicFile {
  mut file_path : String
  file_type : FileType
  index_path : String
  seqnames : Array[String]
  seq_lengths : Map[String, Int]
}

///|
/// Create a new GenomicFile with the given path and type.
/// Initializes empty seqnames and seq_lengths.
pub fn GenomicFile::new(
  file_path : String,
  file_type : FileType,
) -> GenomicFile {
  GenomicFile::{
    file_path,
    file_type,
    index_path: file_path + ".idx",
    seqnames: Array::new(),
    seq_lengths: Map([], capacity=10),
  }
}

///|
/// Add a sequence name and length to the file metadata.
pub fn GenomicFile::add_seqname(
  self : GenomicFile,
  seqname : String,
  length : Int,
) -> Unit {
  self.seqnames.push(seqname)
  self.seq_lengths[seqname] = length
}

///|
/// String representation of GenomicFile.
pub fn GenomicFile::to_string(self : GenomicFile) -> String {
  let mut result = "GenomicFile:\n"
  result = result + "  Path: " + self.file_path + "\n"
  result = result + "  Type: " + self.file_type.to_string() + "\n"
  result = result + "  Index: " + self.index_path + "\n"
  result = result + "  Sequences: " + self.seqnames.length().to_string() + "\n"
  let mut i = 0
  while i < self.seqnames.length() {
    let name = self.seqnames[i]
    let len = self.seq_lengths.get(name)
    let len_str = if len.is_some() { len.unwrap().to_string() } else { "0" }
    result = result + "    " + name + ": " + len_str + "\n"
    i = i + 1
  }
  result
}

///|
/// GenomicFilesCollection - a collection of genomic files for batch processing.
pub struct GenomicFilesCollection {
  files : Array[GenomicFile]
  sample_names : Array[String]
  mut yield_size : Int
}

///|
/// Create a new GenomicFilesCollection.
pub fn GenomicFilesCollection::new(
  files : Array[GenomicFile],
  sample_names : Array[String],
) -> GenomicFilesCollection {
  GenomicFilesCollection::{ files, sample_names, yield_size: 1000 }
}

///|
/// Set the yield size (number of records per batch) for the collection.
pub fn GenomicFilesCollection::set_yield_size(
  self : GenomicFilesCollection,
  yield_size : Int,
) -> Unit {
  self.yield_size = yield_size
}

///|
/// Get the number of files in the collection.
pub fn GenomicFilesCollection::length(self : GenomicFilesCollection) -> Int {
  self.files.length()
}

///|
/// ScanResult - results of scanning a file over a set of ranges.
pub struct ScanResult {
  seqname : String
  ranges : Array[(Int, Int)]
  values : Array[Array[Double]]
}

///|
/// Create a new empty ScanResult for a given seqname.
pub fn ScanResult::new(seqname : String) -> ScanResult {
  ScanResult::{ seqname, ranges: Array::new(), values: Array::new() }
}

///|
/// Add a range and its associated values to the scan result.
pub fn ScanResult::add(
  self : ScanResult,
  range : (Int, Int),
  values : Array[Double],
) -> Unit {
  self.ranges.push(range)
  self.values.push(values)
}

///|
/// String representation of ScanResult.
pub fn ScanResult::to_string(self : ScanResult) -> String {
  let mut result = "ScanResult:\n"
  result = result + "  Seqname: " + self.seqname + "\n"
  result = result + "  Ranges: " + self.ranges.length().to_string() + "\n"
  let mut i = 0
  while i < self.ranges.length() {
    let (start, end) = self.ranges[i]
    result = result +
      "    [" +
      start.to_string() +
      "-" +
      end.to_string() +
      "]: "
    let vals = self.values[i]
    let mut j = 0
    while j < vals.length() {
      if j > 0 {
        result = result + ", "
      }
      result = result + vals[j].to_string()
      j = j + 1
    }
    result = result + "\n"
    i = i + 1
  }
  result
}

///|
/// ReduceResult - aggregated results after applying reduce function.
pub struct ReduceResult {
  summary : Map[String, Double]
  mut total_ranges : Int
  mut total_values : Int
}

///|
/// Create a new empty ReduceResult.
pub fn ReduceResult::new() -> ReduceResult {
  ReduceResult::{
    summary: Map([], capacity=10),
    total_ranges: 0,
    total_values: 0,
  }
}

///|
/// String representation of ReduceResult.
pub fn ReduceResult::to_string(self : ReduceResult) -> String {
  let mut result = "ReduceResult:\n"
  result = result + "  Total ranges: " + self.total_ranges.to_string() + "\n"
  result = result + "  Total values: " + self.total_values.to_string() + "\n"
  result = result + "  Summary:\n"
  let keys = self.summary.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    let val = self.summary.get(key)
    let val_str = if val.is_some() { val.unwrap().to_string() } else { "0" }
    result = result + "    " + key + ": " + val_str + "\n"
    i = i + 1
  }
  result
}

///|
/// Square root using Newton's method.
pub fn genomic_files_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x == 1.0 {
    1.0
  } else {
    let mut z = x
    let mut prev = z
    let mut i = 0
    while i < 30 {
      z = (z + x / z) / 2.0
      if (z - prev).abs() < 0.0000001 {
        break
      }
      prev = z
      i = i + 1
    }
    z
  }
}

///|
/// Check if two intervals [start1, end1] and [start2, end2] overlap.
pub fn genomic_files_is_overlapping(
  start1 : Int,
  end1 : Int,
  start2 : Int,
  end2 : Int,
) -> Bool {
  start1 <= end2 && start2 <= end1
}

///|
/// Parse a BED format line: chrom, start, end, name, score, strand...
/// Returns (chrom, start, end, score). Score defaults to 0.0 if absent.
pub fn genomic_files_parse_bed_line(
  line : String,
) -> (String, Int, Int, Double) {
  let parts = split_by_char(line, 9)
  if parts.length() < 3 {
    return ("", 0, 0, 0.0)
  }
  let chrom = parts[0]
  let start = parse_int(parts[1])
  let end = parse_int(parts[2])
  let score = if parts.length() >= 5 {
    let parsed = parse_double(parts[4])
    if parsed.is_some() {
      parsed.unwrap()
    } else {
      0.0
    }
  } else {
    0.0
  }
  (chrom, start, end, score)
}

///|
/// Parse a VCF format line: chrom, pos, id, ref, alt, qual, filter, info...
/// Returns (chrom, pos, ref, alt). Skips header lines starting with '#'.
pub fn genomic_files_parse_vcf_line(
  line : String,
) -> (String, Int, String, String) {
  if line.length() == 0 {
    return ("", 0, "", "")
  }
  if line.unsafe_get(0).to_int() == '#'.to_int() {
    return ("", 0, "", "")
  }
  let parts = split_by_char(line, 9)
  if parts.length() < 5 {
    return ("", 0, "", "")
  }
  let chrom = parts[0]
  let pos = parse_int(parts[1])
  let ref_allele = parts[3]
  let alt_allele = parts[4]
  (chrom, pos, ref_allele, alt_allele)
}

///|
/// Parse a SAM format line: qname, flag, rname, pos, mapq, cigar, etc.
/// Returns (rname, pos, mapq, cigar, qname) or empty values if header.
pub fn genomic_files_parse_sam_line(
  line : String,
) -> (String, Int, Int, String, String) {
  if line.length() == 0 {
    return ("", 0, 0, "", "")
  }
  if line.unsafe_get(0).to_int() == '@'.to_int() {
    return ("", 0, 0, "", "")
  }
  let parts = split_by_char(line, 9)
  if parts.length() < 6 {
    return ("", 0, 0, "", "")
  }
  let qname = parts[0]
  let rname = parts[2]
  let pos = parse_int(parts[3])
  let mapq = parse_int(parts[4])
  let cigar = parts[5]
  (rname, pos, mapq, cigar, qname)
}

///|
/// Scan a BAM file (or SAM text content) over the given ranges.
/// For each range, computes [read_count, mean_mapq].
/// The file_path is treated as SAM text content if it contains newlines,
/// otherwise a deterministic simulation based on range length is used.
pub fn genomic_files_scan_bam(
  file : GenomicFile,
  ranges : Array[(String, Int, Int)],
) -> ScanResult {
  if ranges.length() == 0 {
    return ScanResult::new("")
  }
  let first_seqname = ranges[0].0
  let result = ScanResult::new(first_seqname)
  let has_content = contains(file.file_path, "\n")
  let mut i = 0
  while i < ranges.length() {
    let (seqname, start, end) = ranges[i]
    let values : Array[Double] = Array::new()
    let mut read_count = 0
    let mut mapq_sum = 0
    if has_content {
      let lines = split_by_char(file.file_path, '\n'.to_int())
      let mut j = 0
      while j < lines.length() {
        let (rname, pos, mapq, _cigar, _qname) = genomic_files_parse_sam_line(
          lines[j],
        )
        if rname.length() > 0 && rname == seqname {
          let read_end = pos + 100
          if genomic_files_is_overlapping(pos, read_end, start, end) {
            read_count = read_count + 1
            mapq_sum = mapq_sum + mapq
          }
        }
        j = j + 1
      }
    } else {
      let range_len = end - start + 1
      read_count = range_len / 50
      if read_count < 1 {
        read_count = 1
      }
      mapq_sum = read_count * 30
    }
    values.push(read_count.to_double())
    let mean_mapq = if read_count > 0 {
      mapq_sum.to_double() / read_count.to_double()
    } else {
      0.0
    }
    values.push(mean_mapq)
    result.add((start, end), values)
    i = i + 1
  }
  result
}

///|
/// Scan a BED file over the given ranges.
/// For each range, computes [feature_count, sum_of_scores, mean_score].
/// The file_path is treated as BED text content.
pub fn genomic_files_scan_bed(
  file : GenomicFile,
  ranges : Array[(String, Int, Int)],
) -> ScanResult {
  if ranges.length() == 0 {
    return ScanResult::new("")
  }
  let first_seqname = ranges[0].0
  let result = ScanResult::new(first_seqname)
  let bed_features : Array[(String, Int, Int, Double)] = Array::new()
  if contains(file.file_path, "\n") {
    let lines = split_by_char(file.file_path, '\n'.to_int())
    let mut j = 0
    while j < lines.length() {
      let line = trim(lines[j])
      if line.length() > 0 && line.unsafe_get(0).to_int() != '#'.to_int() {
        let (chrom, start, end, score) = genomic_files_parse_bed_line(line)
        if chrom.length() > 0 {
          bed_features.push((chrom, start, end, score))
        }
      }
      j = j + 1
    }
  }
  let mut i = 0
  while i < ranges.length() {
    let (seqname, start, end) = ranges[i]
    let values : Array[Double] = Array::new()
    let mut feature_count = 0
    let mut score_sum = 0.0
    let mut k = 0
    while k < bed_features.length() {
      let (chrom, fstart, fend, score) = bed_features[k]
      if chrom == seqname &&
        genomic_files_is_overlapping(start, end, fstart, fend) {
        feature_count = feature_count + 1
        score_sum = score_sum + score
      }
      k = k + 1
    }
    values.push(feature_count.to_double())
    values.push(score_sum)
    let mean_score = if feature_count > 0 {
      score_sum / feature_count.to_double()
    } else {
      0.0
    }
    values.push(mean_score)
    result.add((start, end), values)
    i = i + 1
  }
  result
}

///|
/// Scan a VCF file over the given ranges.
/// For each range, computes [variant_count, mean_qual].
/// The file_path is treated as VCF text content.
pub fn genomic_files_scan_vcf(
  file : GenomicFile,
  ranges : Array[(String, Int, Int)],
) -> ScanResult {
  if ranges.length() == 0 {
    return ScanResult::new("")
  }
  let first_seqname = ranges[0].0
  let result = ScanResult::new(first_seqname)
  let vcf_records : Array[(String, Int, String, String, Double)] = Array::new()
  if contains(file.file_path, "\n") {
    let lines = split_by_char(file.file_path, '\n'.to_int())
    let mut j = 0
    while j < lines.length() {
      let line = lines[j]
      if line.length() > 0 && line.unsafe_get(0).to_int() != '#'.to_int() {
        let (chrom, pos, ref_allele, alt_allele) = genomic_files_parse_vcf_line(
          line,
        )
        if chrom.length() > 0 {
          let parts = split_by_char(line, 9)
          let qual = if parts.length() >= 6 {
            let parsed = parse_double(parts[5])
            if parsed.is_some() {
              parsed.unwrap()
            } else {
              0.0
            }
          } else {
            0.0
          }
          vcf_records.push((chrom, pos, ref_allele, alt_allele, qual))
        }
      }
      j = j + 1
    }
  }
  let mut i = 0
  while i < ranges.length() {
    let (seqname, start, end) = ranges[i]
    let values : Array[Double] = Array::new()
    let mut variant_count = 0
    let mut qual_sum = 0.0
    let mut k = 0
    while k < vcf_records.length() {
      let (chrom, pos, _ref_a, _alt_a, qual) = vcf_records[k]
      if chrom == seqname && pos >= start && pos <= end {
        variant_count = variant_count + 1
        qual_sum = qual_sum + qual
      }
      k = k + 1
    }
    values.push(variant_count.to_double())
    let mean_qual = if variant_count > 0 {
      qual_sum / variant_count.to_double()
    } else {
      0.0
    }
    values.push(mean_qual)
    result.add((start, end), values)
    i = i + 1
  }
  result
}

///|
/// Reduce scan results by range using the specified reduce function.
/// Supported reduce_func values: "sum", "mean", "max", "min".
/// Aggregates values across all scan results for each range index.
pub fn genomic_files_reduce_by_range(
  scan_results : Array[ScanResult],
  reduce_func : String,
) -> ReduceResult {
  let result = ReduceResult::new()
  if scan_results.length() == 0 {
    return result
  }
  let first = scan_results[0]
  let n_ranges = first.ranges.length()
  result.total_ranges = n_ranges
  let mut total_values = 0
  let mut i = 0
  while i < n_ranges {
    let (start, end) = first.ranges[i]
    let key = "range_" +
      i.to_string() +
      "_[" +
      start.to_string() +
      "-" +
      end.to_string() +
      "]"
    let mut aggregated = 0.0
    let mut count = 0
    let mut max_val = 0.0
    let mut min_val = 0.0
    let mut has_value = false
    let mut j = 0
    while j < scan_results.length() {
      if i < scan_results[j].ranges.length() {
        let vals = scan_results[j].values[i]
        let mut k = 0
        while k < vals.length() {
          let v = vals[k]
          if reduce_func == "sum" {
            aggregated = aggregated + v
          } else if reduce_func == "mean" {
            aggregated = aggregated + v
            count = count + 1
          } else if reduce_func == "max" {
            if !has_value || v > max_val {
              max_val = v
            }
          } else if reduce_func == "min" {
            if !has_value || v < min_val {
              min_val = v
            }
          }
          has_value = true
          total_values = total_values + 1
          k = k + 1
        }
      }
      j = j + 1
    }
    let final_val = if reduce_func == "sum" {
      aggregated
    } else if reduce_func == "mean" {
      if count > 0 {
        aggregated / count.to_double()
      } else {
        0.0
      }
    } else if reduce_func == "max" {
      max_val
    } else if reduce_func == "min" {
      min_val
    } else {
      aggregated
    }
    result.summary[key] = final_val
    i = i + 1
  }
  result.total_values = total_values
  result
}

///|
/// Yield reads in batches of yield_size from a BAM/SAM file.
/// Returns an array of batch strings, each containing up to yield_size reads.
pub fn genomic_files_yield_reads(
  file : GenomicFile,
  yield_size : Int,
) -> Array[String] {
  let result : Array[String] = Array::new()
  if yield_size <= 0 {
    return result
  }
  if !contains(file.file_path, "\n") {
    return result
  }
  let lines = split_by_char(file.file_path, '\n'.to_int())
  let current_batch : Array[String] = Array::new()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() > 0 && line.unsafe_get(0).to_int() != '@'.to_int() {
      current_batch.push(line)
      if current_batch.length() >= yield_size {
        let mut batch_str = ""
        let mut k = 0
        while k < current_batch.length() {
          if k > 0 {
            batch_str = batch_str + "\n"
          }
          batch_str = batch_str + current_batch[k]
          k = k + 1
        }
        result.push(batch_str)
        current_batch.clear()
      }
    }
    i = i + 1
  }
  if current_batch.length() > 0 {
    let mut batch_str = ""
    let mut k = 0
    while k < current_batch.length() {
      if k > 0 {
        batch_str = batch_str + "\n"
      }
      batch_str = batch_str + current_batch[k]
      k = k + 1
    }
    result.push(batch_str)
  }
  result
}

///|
/// Count reads in each range for a BAM/SAM file.
/// Returns a Map keyed by "seqname:start-end" with read counts.
pub fn genomic_files_count_reads(
  file : GenomicFile,
  ranges : Array[(String, Int, Int)],
) -> Map[String, Int] {
  let result : Map[String, Int] = Map([], capacity=ranges.length())
  let has_content = contains(file.file_path, "\n")
  let mut i = 0
  while i < ranges.length() {
    let (seqname, start, end) = ranges[i]
    let key = seqname + ":" + start.to_string() + "-" + end.to_string()
    let mut count = 0
    if has_content {
      let lines = split_by_char(file.file_path, '\n'.to_int())
      let mut j = 0
      while j < lines.length() {
        let (rname, pos, _mapq, _cigar, _qname) = genomic_files_parse_sam_line(
          lines[j],
        )
        if rname.length() > 0 && rname == seqname {
          let read_end = pos + 100
          if genomic_files_is_overlapping(pos, read_end, start, end) {
            count = count + 1
          }
        }
        j = j + 1
      }
    } else {
      let range_len = end - start + 1
      count = range_len / 50
      if count < 1 {
        count = 1
      }
    }
    result[key] = count
    i = i + 1
  }
  result
}

///|
/// Compute coverage for each range in a BAM/SAM file.
/// Returns an Array[Double] with mean coverage per range.
pub fn genomic_files_coverage(
  file : GenomicFile,
  ranges : Array[(String, Int, Int)],
) -> Array[Double] {
  let result : Array[Double] = Array::make(ranges.length(), 0.0)
  let has_content = contains(file.file_path, "\n")
  let mut i = 0
  while i < ranges.length() {
    let (seqname, start, end) = ranges[i]
    let range_len = end - start + 1
    let mut total_depth = 0.0
    if has_content {
      let lines = split_by_char(file.file_path, '\n'.to_int())
      let mut j = 0
      while j < lines.length() {
        let (rname, pos, _mapq, _cigar, _qname) = genomic_files_parse_sam_line(
          lines[j],
        )
        if rname.length() > 0 && rname == seqname {
          let read_end = pos + 100
          if genomic_files_is_overlapping(pos, read_end, start, end) {
            let overlap_start = if pos > start { pos } else { start }
            let overlap_end = if read_end < end { read_end } else { end }
            let overlap_len = overlap_end - overlap_start + 1
            if overlap_len > 0 {
              total_depth = total_depth + overlap_len.to_double()
            }
          }
        }
        j = j + 1
      }
    } else {
      let simulated_reads = range_len / 50
      if simulated_reads < 1 {
        total_depth = range_len.to_double()
      } else {
        total_depth = simulated_reads.to_double() * 100.0
      }
    }
    let coverage = if range_len > 0 {
      total_depth / range_len.to_double()
    } else {
      0.0
    }
    result[i] = coverage
    i = i + 1
  }
  result
}

///|
/// Summarize ranges across all files in a collection.
/// Returns an Array[ScanResult], one per file in the collection.
pub fn genomic_files_summarize_ranges(
  collection : GenomicFilesCollection,
  ranges : Array[(String, Int, Int)],
) -> Array[ScanResult] {
  let result : Array[ScanResult] = Array::new()
  let mut i = 0
  while i < collection.files.length() {
    let file = collection.files[i]
    let scan_result = match file.file_type {
      FileType::BAM => genomic_files_scan_bam(file, ranges)
      FileType::BED => genomic_files_scan_bed(file, ranges)
      FileType::VCF => genomic_files_scan_vcf(file, ranges)
      FileType::BigWig => genomic_files_scan_bam(file, ranges)
      FileType::Unknown => genomic_files_scan_bam(file, ranges)
    }
    result.push(scan_result)
    i = i + 1
  }
  result
}

///|
/// Create an example GenomicFile for testing.
/// Returns a BAM file with sample SAM content.
pub fn genomic_files_create_example_file() -> GenomicFile {
  let file = GenomicFile::new("example.bam", FileType::BAM)
  file.add_seqname("chr1", 100000)
  file.add_seqname("chr2", 80000)
  let content = "read1\t0\tchr1\t100\t30\t100M\t*\t0\t0\tACGT\tIIII\n" +
    "read2\t0\tchr1\t150\t40\t100M\t*\t0\t0\tACGT\tIIII\n" +
    "read3\t0\tchr1\t200\t35\t100M\t*\t0\t0\tACGT\tIIII\n" +
    "read4\t0\tchr2\t500\t25\t100M\t*\t0\t0\tACGT\tIIII\n" +
    "read5\t0\tchr2\t550\t45\t100M\t*\t0\t0\tACGT\tIIII"
  file.file_path = content
  file
}

///|
/// Create an example GenomicFilesCollection for testing.
pub fn genomic_files_create_example_collection() -> GenomicFilesCollection {
  let files : Array[GenomicFile] = Array::new()
  let bam_file = genomic_files_create_example_file()
  files.push(bam_file)
  let bed_file = GenomicFile::new(
    "chr1\t100\t200\tfeature1\t10.5\n" +
    "chr1\t150\t250\tfeature2\t20.0\n" +
    "chr2\t500\t600\tfeature3\t15.0",
    FileType::BED,
  )
  files.push(bed_file)
  let vcf_file = GenomicFile::new(
    "chr1\t100\tvar1\tA\tT\t30.0\tPASS\t.\n" +
    "chr1\t180\tvar2\tG\tC\t40.5\tPASS\t.\n" +
    "chr2\t520\tvar3\tC\tG\t25.0\tPASS\t.",
    FileType::VCF,
  )
  files.push(vcf_file)
  let sample_names : Array[String] = ["sample_bam", "sample_bed", "sample_vcf"]
  GenomicFilesCollection::new(files, sample_names)
}

///|
/// Create example ranges for testing.
pub fn genomic_files_create_example_ranges() -> Array[(String, Int, Int)] {
  [("chr1", 100, 200), ("chr1", 150, 250), ("chr2", 500, 600)]
}
