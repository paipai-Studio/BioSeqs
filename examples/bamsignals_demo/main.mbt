///|
/// bamsignals Demo - Signal Extraction and Normalization from BAM files.
/// Inspired by Bioconductor bamsignals package.

fn main {
  println("=== bamsignals Demo ===")

  // 1. Create example regions
  println("")
  println("1. Create example regions:")
  let regions = @bio.bamsig_create_example_regions()
  let mut i = 0
  while i < regions.length() {
    let r = regions[i]
    println("   " + r.bamsig_id() +
      ": " + r.bamsig_chrom() + ":" + r.bamsig_start().to_string() +
      "-" + r.bamsig_end().to_string() +
      " (width=" + r.bamsig_width().to_string() + ")")
    i = i + 1
  }

  // 2. Create example records
  println("")
  println("2. Create example records:")
  let records = @bio.bamsig_create_example_records()
  let mut j = 0
  while j < records.length() {
    let rec = records[j]
    println("   " + rec.chrom + ":" + rec.pos.to_string() +
      " " + rec.cigar +
      " mapq=" + rec.mapq.to_string() +
      " strand=" + rec.strand.to_string())
    j = j + 1
  }

  // 3. Build Map of records by chromosome
  println("")
  println("3. Build records-by-chromosome map:")
  let chr1_records : Array[@bio.BamsigRecord] = Array::new()
  let chr2_records : Array[@bio.BamsigRecord] = Array::new()
  let mut k = 0
  while k < records.length() {
    let rec = records[k]
    if rec.chrom == "chr1" {
      chr1_records.push(rec)
    } else {
      chr2_records.push(rec)
    }
    k = k + 1
  }
  let records_map : Map[String, Array[@bio.BamsigRecord]] = Map(
    [("chr1", chr1_records), ("chr2", chr2_records)],
    capacity=2,
  )
  println("   Chromosomes in map: " + records_map.length().to_string())

  // 4. Create params and count signals
  println("")
  println("4. Count signals in regions:")
  let params = @bio.BamsigParams::new()
  println("   Params: count_mode=All, filter_dup=false, min_mapq=10")
  let counts = @bio.bamsig_count_signals(regions, records_map, params)
  println("   Signal counts:")
  let mut r = 0
  while r < counts.length() {
    let row = counts[r]
    let region_id = regions[r].bamsig_id()
    println("   " + region_id + ": " + row[0].to_string())
    r = r + 1
  }

  // 5. Create example signal (includes RPM normalization internally)
  println("")
  println("5. Create example signal (with RPM normalization):")
  let signal = @bio.bamsig_create_example_signal()
  println("   Regions: " + signal.bamsig_n_regions().to_string())
  println("   Samples: " + signal.bamsig_n_samples().to_string())
  println("   Signal details (counts and RPM-normalized values):")
  let mut s = 0
  while s < signal.bamsig_n_regions() {
    let region_id = signal.bamsig_get_region_id(s)
    let count = signal.bamsig_get_count(s, 0)
    let norm = signal.bamsig_get_norm_count(s, 0)
    println("   " + region_id +
      ": count=" + count.to_string() +
      ", norm(RPM)=" + norm.to_string())
    s = s + 1
  }

  // 6. Filter signal by count threshold
  println("")
  println("6. Filter signal (min_count=10, min_samples=1):")
  let filtered = signal.bamsig_filter(10.0, 1)
  println("   Regions before: " + signal.bamsig_n_regions().to_string())
  println("   Regions after: " + filtered.bamsig_n_regions().to_string())

  // 7. Row and column sums
  println("")
  println("7. Signal statistics:")
  let mut t = 0
  while t < signal.bamsig_n_regions() {
    let region_id = signal.bamsig_get_region_id(t)
    let row_sum = signal.bamsig_row_sum(t)
    println("   " + region_id + " row_sum=" + row_sum.to_string())
    t = t + 1
  }
  let col_sum = signal.bamsig_col_sum(0)
  println("   Sample 0 col_sum=" + col_sum.to_string())

  println("")
  println("=== bamsignals Demo Completed ===")
}