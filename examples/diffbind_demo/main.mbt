///|
fn main {
  println("=== DiffBind Demo ===")

  // 1. Create example DBAConfig (2 conditions, 2 replicates each)
  println("\n1. Create example DBAConfig:")
  let config = @src.diffbind_create_example_config()
  println("  Samples: " + config.samples.length().to_string())
  let mut i = 0
  while i < config.samples.length() {
    let s = config.samples[i]
    println(
      "  " +
      s.sample_id +
      ": condition=" +
      s.condition +
      ", replicate=" +
      s.replicate.to_string() +
      ", peaks=" +
      s.peaks.length().to_string(),
    )
    i = i + 1
  }

  // 2. Show example peaks and create new DBAPeak
  println("\n2. Create new DBAPeak via constructor:")
  let new_peak = @src.DBAPeak::new("chr5", 7000, 8000)
  println("  " + new_peak.to_string())

  // 3. peakset - add a new sample
  println("\n3. peakset: add a new sample:")
  let new_sample = @src.DBASample::new("C1", "ConditionC", 1)
  let config2 = @src.diffbind_peakset(config, new_sample)
  println(
    "  After adding: " +
    config2.samples.length().to_string() +
    " samples (was " +
    config.samples.length().to_string() +
    ")",
  )

  // 4. count - count peaks overlapping ranges
  println("\n4. count: count peaks in query ranges:")
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 1500, 3500),
    ("chr2", 1000, 5500),
    ("chr3", 1000, 3000),
  ]
  let counts = @src.diffbind_count(config, ranges)
  let count_keys = counts.keys().collect()
  let mut j = 0
  while j < count_keys.length() {
    let key = count_keys[j]
    let val = counts.get(key)
    let val_str = if val.is_some() { val.unwrap().to_string() } else { "0" }
    println("  " + key + ": " + val_str + " peaks")
    j = j + 1
  }

  // 5. overlap - peak overlap analysis
  println("\n5. overlap: peak consistency analysis (min_overlap=0.5):")
  let overlap_result = @src.diffbind_overlap(config, 0.5)
  println(overlap_result.to_string())

  // 6. consensus_peaks
  println("\n6. find_consensus_peaks:")
  let consensus = @src.diffbind_find_consensus_peaks(overlap_result)
  println("  Consensus peaks: " + consensus.length().to_string())
  let mut k = 0
  while k < consensus.length() {
    println("  " + k.to_string() + ": " + consensus[k].to_string())
    k = k + 1
  }

  // 7. merge_peaks
  println("\n7. merge_peaks:")
  let peaks_a = @src.diffbind_create_example_peaks()
  let peaks_b = @src.diffbind_create_example_peaks()
  let merged = @src.diffbind_merge_peaks(peaks_a, peaks_b)
  println(
    "  Merged " +
    peaks_a.length().to_string() +
    " + " +
    peaks_b.length().to_string() +
    " peaks into " +
    merged.length().to_string() +
    " peaks",
  )
  let mut m = 0
  while m < merged.length() {
    println("  " + m.to_string() + ": " + merged[m].to_string())
    m = m + 1
  }

  // 8. normalize (TMM)
  println("\n8. normalize (TMM):")
  let norm_matrix = @src.diffbind_normalize(config, "TMM")
  println(
    "  Normalized matrix: " +
    norm_matrix.length().to_string() +
    " peaks x " +
    (if norm_matrix.length() > 0 {
      norm_matrix[0].length().to_string()
    } else {
      "0".to_string()
    }) +
    " samples",
  )
  let mut p = 0
  while p < norm_matrix.length() && p < 3 {
    let mut row_str = "  peak " + p.to_string() + ": ["
    let mut q = 0
    while q < norm_matrix[p].length() {
      if q > 0 {
        row_str = row_str + ", "
      }
      row_str = row_str + norm_matrix[p][q].to_string()
      q = q + 1
    }
    row_str = row_str + "]"
    println(row_str)
    p = p + 1
  }

  // 9. analyze - differential binding analysis
  println("\n9. analyze: differential binding (ConditionA vs ConditionB):")
  let result = @src.diffbind_analyze(config, "ConditionA", "ConditionB")
  println(result.to_string())

  // 10. report - significant peaks
  println("\n10. report: significant peaks (FDR < 0.05):")
  let significant = @src.diffbind_report(result, 0.05)
  println("  Significant peaks: " + significant.length().to_string())
  let mut n = 0
  while n < significant.length() {
    println("  " + n.to_string() + ": " + significant[n].to_string())
    n = n + 1
  }

  // 11. volcano_data and ma_data from example result
  println("\n11. Visualization data (from example result):")
  let example_result = @src.diffbind_create_example_result()
  let volcano = @src.diffbind_volcano_data(example_result)
  println("  Volcano data points: " + volcano.length().to_string())
  let mut v = 0
  while v < volcano.length() {
    let (fc, nlp) = volcano[v]
    println(
      "  " +
      v.to_string() +
      ": log2FC=" +
      fc.to_string() +
      ", -log10(p)=" +
      nlp.to_string(),
    )
    v = v + 1
  }
  let ma = @src.diffbind_ma_data(example_result)
  println("  MA data points: " + ma.length().to_string())
  let mut w = 0
  while w < ma.length() {
    let (mean_log2, fc) = ma[w]
    println(
      "  " +
      w.to_string() +
      ": mean_log2=" +
      mean_log2.to_string() +
      ", log2FC=" +
      fc.to_string(),
    )
    w = w + 1
  }

  println("\nDiffBind demo completed!")
}
