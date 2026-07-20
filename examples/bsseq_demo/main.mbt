///|
/// bsseq Demo.
fn main {
  println("=== bsseq Demo ===")

  // 1. Create example BS data
  println("\n1. Create example BS data:")
  let bsdata = @src.bsseq_create_example_bsdata()
  println(bsdata.to_string())

  // 2. Demonstrate methylation rate calculation
  println("\n2. Methylation Rate Calculation:")
  let rate = @src.bsseq_get_methylation_rate(8, 2)
  println("  rate(meth=8, unmeth=2) = " + rate.to_string())

  // 3. Demonstrate coverage retrieval
  println("\n3. Coverage retrieval:")
  let coverage = @src.bsseq_get_coverage(bsdata)
  println("  Coverage array length: " + coverage.length().to_string())
  let mut i = 0
  while i < coverage.length() && i < 5 {
    println("    [" + i.to_string() + "] " + coverage[i].to_string())
    i = i + 1
  }

  // 4. Demonstrate methylation proportion retrieval
  println("\n4. Methylation proportion retrieval:")
  let methylation = @src.bsseq_get_methylation(bsdata)
  println("  Methylation array length: " + methylation.length().to_string())
  let mut j = 0
  while j < methylation.length() && j < 5 {
    println("    [" + j.to_string() + "] " + methylation[j].to_string())
    j = j + 1
  }

  // 5. Demonstrate coverage filtering
  println("\n5. Filter by coverage (5..15):")
  let filtered_cov = @src.bsseq_filter_by_coverage(bsdata, 5, 15)
  println(
    "  Records before: " +
    bsdata.n_records().to_string() +
    ", after: " +
    filtered_cov.n_records().to_string(),
  )

  // 6. Demonstrate context filtering (CpG / CHG / CHH)
  println("\n6. Filter by context:")
  let cpg_data = @src.bsseq_filter_by_context(bsdata, "CpG")
  let chg_data = @src.bsseq_filter_by_context(bsdata, "CHG")
  let chh_data = @src.bsseq_filter_by_context(bsdata, "CHH")
  println("  CpG records: " + cpg_data.n_records().to_string())
  println("  CHG records: " + chg_data.n_records().to_string())
  println("  CHH records: " + chh_data.n_records().to_string())

  // 7. Demonstrate methylation smoothing
  println("\n7. Smooth methylation (window=2):")
  let meth_arr : Array[Int] = [8, 9, 5, 3, 7, 6, 4, 2, 8, 9]
  let cov_arr : Array[Int] = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
  let smoothed = @src.bsseq_smooth_methylation(meth_arr, cov_arr, 2)
  println("  Smoothed values: " + smoothed.length().to_string())
  let mut k = 0
  while k < smoothed.length() && k < 5 {
    println("    [" + k.to_string() + "] " + smoothed[k].to_string())
    k = k + 1
  }

  // 8. Demonstrate BSmooth algorithm
  println("\n8. BSmooth algorithm (window=3):")
  let bsmooth_result = @src.bsseq_bsmooth(bsdata, 3)
  println(bsmooth_result.to_string())

  // 9. Demonstrate t-statistic computation
  println("\n9. Compute t-statistics:")
  let group1 : Array[Double] = [
    0.8, 0.9, 0.5, 0.3, 0.7, 0.6, 0.4, 0.2, 0.8, 0.9,
  ]
  let group2 : Array[Double] = [
    0.4, 0.3, 0.5, 0.8, 0.7, 0.2, 0.1, 0.9, 0.5, 0.6,
  ]
  let tstat = @src.bsseq_compute_tstat(group1, group2)
  println(tstat.to_string())

  // 10. Demonstrate DMR detection
  println("\n10. Find DMRs:")
  let positions : Array[Int] = [
    1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
  ]
  let dmr_result = @src.bsseq_find_dmr(tstat, positions, "chr1", 1.5, 500)
  println(dmr_result.to_string())

  // 11. Demonstrate hyper / hypo methylated regions
  println("\n11. Hyper / Hypo methylated regions:")
  let hyper = @src.bsseq_find_hyper_dmr(dmr_result)
  let hypo = @src.bsseq_find_hypo_dmr(dmr_result)
  println("  Hyper DMRs: " + hyper.length().to_string())
  println("  Hypo DMRs: " + hypo.length().to_string())
  let mut h = 0
  while h < hyper.length() {
    println("    hyper: " + hyper[h].to_string())
    h = h + 1
  }
  let mut l = 0
  while l < hypo.length() {
    println("    hypo: " + hypo[l].to_string())
    l = l + 1
  }

  // 12. Demonstrate merging CpGs
  println("\n12. Merge CpGs (max_gap=150):")
  let records = @src.bsseq_create_example_records()
  let merged = @src.bsseq_merge_cpgs(records, 150)
  println(
    "  Records before: " +
    records.length().to_string() +
    ", after merge: " +
    merged.length().to_string(),
  )

  // 13. Demonstrate region methylation computation
  println("\n13. Region methylation:")
  let region_meth = @src.bsseq_region_methylation(bsdata, "chr1", 1000, 1500)
  println("  Mean methylation in chr1:1000-1500 = " + region_meth.to_string())

  // 14. Demonstrate normalization
  println("\n14. Normalize (global):")
  let normalized = @src.bsseq_normalize(bsdata, "global")
  println("  Normalized records: " + normalized.n_records().to_string())
  let norm_none = @src.bsseq_normalize(bsdata, "none")
  println("  Normalize 'none' records: " + norm_none.n_records().to_string())

  // 15. Demonstrate combining samples
  println("\n15. Combine samples:")
  let cov_matrix = @src.bsseq_create_example_cov_matrix()
  println(cov_matrix.to_string())

  // 16. Demonstrate BSRecord creation
  println("\n16. BSRecord creation:")
  let rec = @src.BSRecord::new("chr2", 5000, "+", 12, 3)
  println("  " + rec.to_string())
  println("  coverage=" + rec.coverage().to_string())
  println("  rate=" + rec.rate().to_string())

  // 17. Demonstrate methylation diff and data summary
  println("\n17. Methylation diff and summary:")
  let diff = @src.bsseq_compute_methylation_diff(0.85, 0.42)
  println("  diff(0.85, 0.42) = " + diff.to_string())
  println(@src.bsseq_summary(bsdata))

  println("\nbsseq demo completed!")
}
