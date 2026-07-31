///|
/// hmisc module demo - Harrell Miscellaneous Statistical Utilities.
///
/// Run: moon run examples/hmisc_demo
fn main {
  println("=== hmisc Module Demo ===")
  println("Harrell Miscellaneous Statistical Utilities")
  println("")

  // Sample data
  println("=== Sample Data ===")
  let data = @bio.hmisc_sample_data()
  let names = @bio.hmisc_sample_names()
  println("Variables: \{names.length()}")
  let mut i = 0
  while i < names.length() {
    println("  \{names[i]}: \{array_to_string(data[i])}")
    i = i + 1
  }
  println("")

  // Descriptive statistics
  println("=== Descriptive Statistics ===")
  i = 0
  while i < names.length() {
    let stats = @bio.hmisc_describe(data[i], name=names[i])
    println(@bio.hmisc_describe_summary(stats))
    println("")
    i = i + 1
  }

  // Pearson correlation
  println("=== Pearson Correlation ===")
  let pearson_corr = @bio.hmisc_rcorr(data, names=names, cor_type=@bio.hmisc_cor_type_pearson())
  println(@bio.hmisc_rcorr_summary(pearson_corr))
  println("")

  // Spearman correlation
  println("=== Spearman Correlation ===")
  let spearman_corr = @bio.hmisc_rcorr(data, names=names, cor_type=@bio.hmisc_cor_type_spearman())
  println(@bio.hmisc_rcorr_summary(spearman_corr))
  println("")

  // Variable clustering
  println("=== Variable Clustering ===")
  let vc_result = @bio.hmisc_varclus(data, names=names, min_cluster_size=2)
  println(@bio.hmisc_varclus_summary(vc_result))
  println("Merge steps: \{vc_result.merges.length()}")
  let mut j = 0
  while j < vc_result.heights.length() {
    let step = vc_result.merges[j]
    println("  Step \{j}: merge \{step[0]} and \{step[1]} at height \{hmisc_fmt(vc_result.heights[j], 4)}")
    j = j + 1
  }
  println("")

  // Somers' d
  println("=== Somers' d Statistic ===")
  let x_data = data[0]
  let g1 : Array[Double] = []
  let g2 : Array[Double] = []
  let mut k = 0
  while k < x_data.length() {
    if x_data[k] <= 5.0 {
      g1.push(x_data[k])
    } else {
      g2.push(x_data[k])
    }
    k = k + 1
  }
  println("Group 1 (X <= 5): \{array_to_string(g1)}")
  println("Group 2 (X > 5): \{array_to_string(g2)}")
  let somers = @bio.hmisc_somers_d(g1, g2)
  println("Somers' d: \{hmisc_fmt(somers.d, 4)}")
  println("  95% CI: [\{hmisc_fmt(somers.lower, 4)}, \{hmisc_fmt(somers.upper, 4)}]")
  println("  n: \{somers.n}")
  println("")

  // Imputation with NaN values
  println("=== Imputation ===")
  let nan : Double = @double.not_a_number
  let messy_data = [
    [1.0, nan, 3.0, 4.0, nan, 6.0, 7.0, nan, 9.0, 10.0],
    [2.1, 3.9, nan, 8.0, 9.9, nan, 14.1, 16.0, nan, 20.0],
    [5.0, 4.0, 6.0, nan, 7.0, 5.0, nan, 4.0, 9.0, 6.0],
  ]
  println("Before imputation (NaN values present):")
  let mut r = 0
  while r < messy_data.length() {
    println("  Var \{names[r]}: \{array_to_string(messy_data[r])}")
    r = r + 1
  }
  let imputed = @bio.hmisc_impute(messy_data)
  println("After imputation (column means):")
  r = 0
  while r < imputed.length() {
    println("  Var \{names[r]}: \{array_to_string(imputed[r])}")
    r = r + 1
  }
  println("")
  println("=== Demo Complete ===")
}

fn array_to_string(arr : Array[Double]) -> String {
  let mut result = "["
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      result = result + ", "
    }
    if arr[i].is_nan() {
      result = result + "NaN"
    } else {
      result = result + arr[i].to_string()
    }
    i = i + 1
  }
  result = result + "]"
  result
}

fn hmisc_fmt(value : Double, decimals : Int) -> String {
  if value.is_nan() {
    return "NaN"
  }
  let rounded = (value * @math.pow(10.0, decimals.to_double())).round() / @math.pow(10.0, decimals.to_double())
  let int_part = rounded.to_int()
  let dec_part = ((rounded - int_part.to_double()) * @math.pow(10.0, decimals.to_double())).round().to_int()
  if decimals == 0 {
    return "\{int_part}"
  }
  let sign = if value < 0.0 { "-" } else { "" }
  let dec_str = pad_zeros(dec_part, decimals)
  "\{sign}\{int_part}.\{dec_str}"
}

fn pad_zeros(value : Int, length : Int) -> String {
  let abs_val = value.abs()
  let s = "\{abs_val}"
  if s.length() >= length {
    return s
  }
  let mut result = ""
  let mut i = 0
  while i < length - s.length() {
    result = result + "0"
    i = i + 1
  }
  result + s
}