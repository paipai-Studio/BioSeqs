///|
/// Tests for the genomic_files module.
/// Mirrors the Bioconductor GenomicFiles package API.

///|
test "genomic_files_create_file_bam" {
  // FileType variants are read-only outside the package, so obtain a BAM
  // file_type from the example factory function.
  let example = @src.genomic_files_create_example_file()
  let file = @src.GenomicFile::new("test.bam", example.file_type)
  assert_eq(file.file_path, "test.bam")
  assert_eq(file.index_path, "test.bam.idx")
  assert_eq(file.seqnames.length(), 0)
}

///|
test "genomic_files_create_file_bed" {
  // Obtain a BED file_type from the example collection (files[1] is BED).
  let collection = @src.genomic_files_create_example_collection()
  let bed_file_type = collection.files[1].file_type
  let file = @src.GenomicFile::new("regions.bed", bed_file_type)
  assert_eq(file.file_path, "regions.bed")
  assert_eq(file.index_path, "regions.bed.idx")
}

///|
test "genomic_files_add_seqname" {
  let example = @src.genomic_files_create_example_file()
  let file = @src.GenomicFile::new("test.bam", example.file_type)
  file.add_seqname("chr1", 100000)
  file.add_seqname("chr2", 80000)
  assert_eq(file.seqnames.length(), 2)
  assert_eq(file.seqnames[0], "chr1")
  assert_eq(file.seqnames[1], "chr2")
  assert_eq(file.seq_lengths.get("chr1"), Some(100000))
  assert_eq(file.seq_lengths.get("chr2"), Some(80000))
}

///|
test "genomic_files_file_type_to_string" {
  // Only BAM/BED/VCF are obtainable from the example factory functions
  // (the enum is read-only outside the package).
  let example = @src.genomic_files_create_example_file()
  assert_eq(@src.FileType::to_string(example.file_type), "BAM")
  let collection = @src.genomic_files_create_example_collection()
  assert_eq(@src.FileType::to_string(collection.files[1].file_type), "BED")
  assert_eq(@src.FileType::to_string(collection.files[2].file_type), "VCF")
}

///|
test "genomic_files_is_overlapping_positive" {
  assert_true(@src.genomic_files_is_overlapping(100, 200, 150, 250))
  assert_true(@src.genomic_files_is_overlapping(100, 500, 200, 300))
}

///|
test "genomic_files_is_overlapping_negative" {
  assert_true(!@src.genomic_files_is_overlapping(100, 200, 300, 400))
  assert_true(!@src.genomic_files_is_overlapping(500, 600, 100, 200))
}

///|
test "genomic_files_is_overlapping_boundary" {
  // Touching at the boundary counts as overlap.
  assert_true(@src.genomic_files_is_overlapping(100, 200, 200, 300))
  assert_true(@src.genomic_files_is_overlapping(100, 200, 50, 100))
}

///|
test "genomic_files_parse_bed_line_basic" {
  let line = "chr1\t100\t200\tfeature1\t10.5"
  let (chrom, start, end, score) = @src.genomic_files_parse_bed_line(line)
  assert_eq(chrom, "chr1")
  assert_eq(start, 100)
  assert_eq(end, 200)
  assert_true(score > 10.4 && score < 10.6)
}

///|
test "genomic_files_parse_bed_line_no_score" {
  let line = "chr2\t300\t400"
  let (chrom, start, end, score) = @src.genomic_files_parse_bed_line(line)
  assert_eq(chrom, "chr2")
  assert_eq(start, 300)
  assert_eq(end, 400)
  assert_eq(score, 0.0)
}

///|
test "genomic_files_parse_bed_line_too_short" {
  let (chrom, _start, _end, _score) = @src.genomic_files_parse_bed_line("chr1")
  assert_eq(chrom, "")
}

///|
test "genomic_files_parse_vcf_line_basic" {
  let line = "chr1\t100\tvar1\tA\tT\t30.0\tPASS\t."
  let (chrom, pos, ref_a, alt_a) = @src.genomic_files_parse_vcf_line(line)
  assert_eq(chrom, "chr1")
  assert_eq(pos, 100)
  assert_eq(ref_a, "A")
  assert_eq(alt_a, "T")
}

///|
test "genomic_files_parse_vcf_line_header" {
  let (chrom, _pos, _ref_a, _alt_a) = @src.genomic_files_parse_vcf_line(
    "#chrom\tPOS\tID",
  )
  assert_eq(chrom, "")
}

///|
test "genomic_files_sqrt" {
  let result = @src.genomic_files_sqrt(4.0)
  assert_true(result > 1.99 && result < 2.01)
  let result2 = @src.genomic_files_sqrt(0.0)
  assert_eq(result2, 0.0)
  let result3 = @src.genomic_files_sqrt(1.0)
  assert_eq(result3, 1.0)
}

///|
test "genomic_files_create_example_file" {
  let file = @src.genomic_files_create_example_file()
  assert_eq(file.seqnames.length(), 2)
  assert_eq(file.seqnames[0], "chr1")
  assert_eq(file.seqnames[1], "chr2")
  assert_eq(file.seq_lengths.get("chr1"), Some(100000))
  assert_eq(file.seq_lengths.get("chr2"), Some(80000))
  // The file content includes 5 SAM records.
  assert_true(file.file_path.contains("read1"))
  assert_true(file.file_path.contains("read5"))
}

///|
test "genomic_files_create_example_collection" {
  let collection = @src.genomic_files_create_example_collection()
  assert_eq(collection.length(), 3)
  assert_eq(collection.sample_names.length(), 3)
  assert_eq(collection.sample_names[0], "sample_bam")
  assert_eq(collection.sample_names[1], "sample_bed")
  assert_eq(collection.sample_names[2], "sample_vcf")
}

///|
test "genomic_files_scan_bam_with_content" {
  let file = @src.genomic_files_create_example_file()
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 100, 200),
    ("chr2", 500, 600),
  ]
  let result = @src.genomic_files_scan_bam(file, ranges)
  assert_eq(result.ranges.length(), 2)
  assert_eq(result.values[0].length(), 2)
  // Each scan returns [read_count, mean_mapq].
  // chr1:100-200 contains read1 (pos 100), read2 (pos 150), read3 (pos 200).
  assert_eq(result.values[0][0], 3.0)
}

///|
test "genomic_files_scan_bam_simulated" {
  // A path with no newlines triggers the deterministic simulation.
  let example = @src.genomic_files_create_example_file()
  let file = @src.GenomicFile::new("no_real_file.bam", example.file_type)
  let ranges : Array[(String, Int, Int)] = [("chr1", 100, 299)]
  let result = @src.genomic_files_scan_bam(file, ranges)
  assert_eq(result.ranges.length(), 1)
  // range_len = 200, read_count = 200/50 = 4.
  assert_eq(result.values[0][0], 4.0)
}

///|
test "genomic_files_scan_bed_with_content" {
  let content = "chr1\t100\t200\tf1\t10.0\nchr1\t150\t250\tf2\t20.0\nchr2\t500\t600\tf3\t15.0"
  let collection = @src.genomic_files_create_example_collection()
  let file = @src.GenomicFile::new(content, collection.files[1].file_type)
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 100, 250),
    ("chr2", 500, 600),
  ]
  let result = @src.genomic_files_scan_bed(file, ranges)
  assert_eq(result.ranges.length(), 2)
  // First range overlaps both chr1 features -> count = 2.
  assert_eq(result.values[0][0], 2.0)
  // Sum of scores 10 + 20 = 30.
  assert_true(result.values[0][1] > 29.9 && result.values[0][1] < 30.1)
  // chr2:500-600 overlaps one feature -> count = 1.
  assert_eq(result.values[1][0], 1.0)
}

///|
test "genomic_files_scan_vcf_with_content" {
  let content = "chr1\t100\tvar1\tA\tT\t30.0\tPASS\t.\nchr1\t180\tvar2\tG\tC\t40.5\tPASS\t.\nchr2\t520\tvar3\tC\tG\t25.0\tPASS\t."
  let collection = @src.genomic_files_create_example_collection()
  let file = @src.GenomicFile::new(content, collection.files[2].file_type)
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 100, 200),
    ("chr2", 500, 600),
  ]
  let result = @src.genomic_files_scan_vcf(file, ranges)
  assert_eq(result.ranges.length(), 2)
  // chr1:100-200 contains 2 variants.
  assert_eq(result.values[0][0], 2.0)
  // chr2:500-600 contains 1 variant.
  assert_eq(result.values[1][0], 1.0)
}

///|
test "genomic_files_reduce_by_range_sum" {
  // Two scan results, each with one range and two values.
  let r1 = @src.ScanResult::new("chr1")
  r1.add((100, 200), [3.0, 30.0])
  let r2 = @src.ScanResult::new("chr1")
  r2.add((100, 200), [1.0, 40.0])
  let result = @src.genomic_files_reduce_by_range([r1, r2], "sum")
  assert_eq(result.total_ranges, 1)
  let keys = result.summary.keys().collect()
  assert_eq(keys.length(), 1)
  // Sum = (3 + 30) + (1 + 40) = 74.
  let v = result.summary.get(keys[0])
  assert_true(v is Some(_))
  assert_true(v.unwrap() > 73.9 && v.unwrap() < 74.1)
}

///|
test "genomic_files_reduce_by_range_mean" {
  let r1 = @src.ScanResult::new("chr1")
  r1.add((100, 200), [10.0, 20.0])
  let result = @src.genomic_files_reduce_by_range([r1], "mean")
  assert_eq(result.total_ranges, 1)
  let keys = result.summary.keys().collect()
  let v = result.summary.get(keys[0])
  assert_true(v is Some(_))
  // Mean of [10, 20] = 15.
  assert_true(v.unwrap() > 14.9 && v.unwrap() < 15.1)
}

///|
test "genomic_files_yield_reads" {
  let file = @src.genomic_files_create_example_file()
  // 5 reads total; yield in batches of 2 -> 3 batches (2, 2, 1).
  let batches = @src.genomic_files_yield_reads(file, 2)
  assert_eq(batches.length(), 3)
}

///|
test "genomic_files_yield_reads_empty" {
  // Without newline characters the file is treated as a path -> no batches.
  let example = @src.genomic_files_create_example_file()
  let file = @src.GenomicFile::new("path_only.bam", example.file_type)
  let batches = @src.genomic_files_yield_reads(file, 10)
  assert_eq(batches.length(), 0)
}

///|
test "genomic_files_count_reads" {
  let file = @src.genomic_files_create_example_file()
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 100, 200),
    ("chr2", 500, 600),
  ]
  let counts = @src.genomic_files_count_reads(file, ranges)
  // chr1:100-200 overlaps reads at 100, 150, 200.
  let c1 = counts.get("chr1:100-200")
  assert_true(c1 is Some(_))
  assert_eq(c1.unwrap(), 3)
  // chr2:500-600 overlaps reads at 500 and 550.
  let c2 = counts.get("chr2:500-600")
  assert_true(c2 is Some(_))
  assert_eq(c2.unwrap(), 2)
}

///|
test "genomic_files_coverage" {
  let file = @src.genomic_files_create_example_file()
  let ranges : Array[(String, Int, Int)] = [("chr1", 100, 200)]
  let coverage = @src.genomic_files_coverage(file, ranges)
  assert_eq(coverage.length(), 1)
  // Coverage is total_depth / range_len. With three 100-bp reads covering
  // at least part of [100, 200] the coverage is > 0.
  assert_true(coverage[0] > 0.0)
}

///|
test "genomic_files_summarize_ranges" {
  let collection = @src.genomic_files_create_example_collection()
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 100, 200),
    ("chr2", 500, 600),
  ]
  let results = @src.genomic_files_summarize_ranges(collection, ranges)
  // One ScanResult per file in the collection (BAM, BED, VCF).
  assert_eq(results.length(), 3)
  let mut i = 0
  while i < results.length() {
    assert_eq(results[i].ranges.length(), 2)
    i = i + 1
  }
}

///|
test "genomic_files_collection_length" {
  let collection = @src.genomic_files_create_example_collection()
  let bam_file_type = collection.files[0].file_type
  let bed_file_type = collection.files[1].file_type
  let file1 = @src.GenomicFile::new("a.bam", bam_file_type)
  let file2 = @src.GenomicFile::new("b.bed", bed_file_type)
  let new_collection = @src.GenomicFilesCollection::new([file1, file2], [
    "s1", "s2",
  ])
  assert_eq(new_collection.length(), 2)
  assert_eq(new_collection.sample_names.length(), 2)
}

///|
test "genomic_files_collection_set_yield_size" {
  let collection = @src.genomic_files_create_example_collection()
  let bam_file_type = collection.files[0].file_type
  let file = @src.GenomicFile::new("a.bam", bam_file_type)
  let new_collection = @src.GenomicFilesCollection::new([file], ["s1"])
  new_collection.set_yield_size(500)
  assert_eq(new_collection.yield_size, 500)
}
