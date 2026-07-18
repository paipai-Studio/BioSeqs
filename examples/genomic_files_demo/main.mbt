///|
fn main {
  println("=== GenomicFiles Demo ===")

  // 1. Create GenomicFile objects (BAM, BED, VCF)
  println("\n1. Create example GenomicFile (BAM):")
  let bam_file = @src.genomic_files_create_example_file()
  println(bam_file.to_string())

  // 2. Create GenomicFilesCollection
  println("\n2. Create GenomicFilesCollection:")
  let collection = @src.genomic_files_create_example_collection()
  println("  Files in collection: " + collection.length().to_string())
  let mut i = 0
  while i < collection.files.length() {
    let f = collection.files[i]
    println(
      "  " +
      i.to_string() +
      ": " +
      f.file_type.to_string() +
      " (" +
      collection.sample_names[i] +
      ")",
    )
    i = i + 1
  }

  // 3. Create example ranges and scan BAM/BED/VCF
  let ranges = @src.genomic_files_create_example_ranges()
  println("\n3. Example ranges:")
  i = 0
  while i < ranges.length() {
    let (chr, start, end) = ranges[i]
    println("  " + chr + ":" + start.to_string() + "-" + end.to_string())
    i = i + 1
  }

  println("\n4. Scan BAM file over ranges:")
  let bam_scan = @src.genomic_files_scan_bam(bam_file, ranges)
  println(bam_scan.to_string())

  println("\n5. Scan BED file over ranges:")
  let bed_file = collection.files[1]
  let bed_scan = @src.genomic_files_scan_bed(bed_file, ranges)
  println(bed_scan.to_string())

  println("\n6. Scan VCF file over ranges:")
  let vcf_file = collection.files[2]
  let vcf_scan = @src.genomic_files_scan_vcf(vcf_file, ranges)
  println(vcf_scan.to_string())

  // 4. reduce_by_range
  println("\n7. Reduce by range (sum):")
  let scan_results : Array[@src.ScanResult] = [bam_scan, bed_scan, vcf_scan]
  let reduce_result = @src.genomic_files_reduce_by_range(scan_results, "sum")
  println(reduce_result.to_string())

  // 5. yield_reads
  println("\n8. Yield reads (batch size 1000):")
  let batches = @src.genomic_files_yield_reads(bam_file, 1000)
  println("  Number of batches: " + batches.length().to_string())
  let mut j = 0
  while j < batches.length() {
    println(
      "  Batch " +
      j.to_string() +
      ": " +
      batches[j].length().to_string() +
      " chars",
    )
    j = j + 1
  }

  // 6. count_reads
  println("\n9. Count reads per range:")
  let counts = @src.genomic_files_count_reads(bam_file, ranges)
  let count_keys = counts.keys().collect()
  j = 0
  while j < count_keys.length() {
    let key = count_keys[j]
    let val = counts.get(key)
    let val_str = if val.is_some() { val.unwrap().to_string() } else { "0" }
    println("  " + key + ": " + val_str)
    j = j + 1
  }

  // 7. coverage
  println("\n10. Coverage per range:")
  let coverage = @src.genomic_files_coverage(bam_file, ranges)
  let mut k = 0
  while k < coverage.length() {
    let (chr, start, end) = ranges[k]
    println(
      "  " +
      chr +
      ":" +
      start.to_string() +
      "-" +
      end.to_string() +
      ": " +
      coverage[k].to_string(),
    )
    k = k + 1
  }

  // 8. summarize_ranges
  println("\n11. Summarize ranges across collection:")
  let summaries = @src.genomic_files_summarize_ranges(collection, ranges)
  println("  Summaries returned: " + summaries.length().to_string())
  let mut m = 0
  while m < summaries.length() {
    println(
      "  File " +
      m.to_string() +
      " (" +
      collection.files[m].file_type.to_string() +
      "): " +
      summaries[m].ranges.length().to_string() +
      " ranges scanned",
    )
    m = m + 1
  }

  // 9. parse BED line
  println("\n12. Parse BED line:")
  let bed_line = "chr1\t100\t200\tfeature1\t10.5"
  let (bchrom, bstart, bend, bscore) = @src.genomic_files_parse_bed_line(
    bed_line,
  )
  println(
    "  chrom=" +
    bchrom +
    ", start=" +
    bstart.to_string() +
    ", end=" +
    bend.to_string() +
    ", score=" +
    bscore.to_string(),
  )

  // 10. parse VCF line
  println("\n13. Parse VCF line:")
  let vcf_line = "chr1\t100\tvar1\tA\tT\t30.0\tPASS\t."
  let (vchrom, vpos, vref, valt) = @src.genomic_files_parse_vcf_line(vcf_line)
  println(
    "  chrom=" +
    vchrom +
    ", pos=" +
    vpos.to_string() +
    ", ref=" +
    vref +
    ", alt=" +
    valt,
  )

  println("\nGenomicFiles demo completed!")
}
