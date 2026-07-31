fn main {
  println("=== Bioconductor-style impute Demo ===\n")

  println("--- 1. Data with missing values ---")
  let data = [
    [1.0, 2.0, Double::nan(), 4.0, 5.0],
    [Double::nan(), 7.0, 8.0, Double::nan(), 10.0],
    [11.0, Double::nan(), 13.0, 14.0, Double::nan()],
  ]
  let stats = @bio.impute_na_stats(data)
  println(@bio.impute_na_summary(stats))

  println("\n--- 2. Impute by Row Mean ---")
  let mean_imp = @bio.impute_by_row_mean(data)
  let mut r = 0
  while r < mean_imp.length() {
    let row = mean_imp[r]
    let mut s = "  row[\{r}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      s = s + row[c].to_string()
      c = c + 1
    }
    s = s + "]"
    println(s)
    r = r + 1
  }

  println("\n--- 3. Impute by Column Median ---")
  let med_imp = @bio.impute_by_col_median(data)
  let mut i = 0
  while i < med_imp.length() {
    let row = med_imp[i]
    let mut s = "  row[\{i}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      s = s + row[c].to_string()
      c = c + 1
    }
    s = s + "]"
    println(s)
    i = i + 1
  }

  println("\n--- 4. LOCF (Last Observation Carried Forward) by row ---")
  let locf_imp = @bio.impute_locf(data, by_row=true)
  let mut j = 0
  while j < locf_imp.length() {
    let row = locf_imp[j]
    let mut s = "  row[\{j}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      let v = row[c]
      s = s + (if v.is_nan() { "NA".to_string() } else { v.to_string() })
      c = c + 1
    }
    s = s + "]"
    println(s)
    j = j + 1
  }

  println("\n--- 5. NOCB (Next Observation Carried Back) by row ---")
  let nocb_imp = @bio.impute_nocb(data, by_row=true)
  let mut k = 0
  while k < nocb_imp.length() {
    let row = nocb_imp[k]
    let mut s = "  row[\{k}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      let v = row[c]
      s = s + (if v.is_nan() { "NA".to_string() } else { v.to_string() })
      c = c + 1
    }
    s = s + "]"
    println(s)
    k = k + 1
  }

  println("\n--- 6. KNN Imputation (k=2) ---")
  let p = @bio.KNNImputeParam::with_k(2)
  let knn_imp = @bio.impute_by_knn(data, p)
  let mut m = 0
  while m < knn_imp.length() {
    let row = knn_imp[m]
    let mut s = "  row[\{m}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      s = s + row[c].to_string()
      c = c + 1
    }
    s = s + "]"
    println(s)
    m = m + 1
  }

  println("\nAll imputation methods completed successfully!")
}
