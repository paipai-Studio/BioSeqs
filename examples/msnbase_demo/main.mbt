///|
/// MSnbase Demo.
///
/// This example demonstrates core mass-spectrometry data structures
/// inspired by Bioconductor's MSnbase package, including:
///   1. Spectrum objects (MS1, MS2) - creation, filtering, centroiding
///   2. Spectrum processing - TIC normalization, smoothing, baseline corr.
///   3. Chromatogram (XIC) - AUC, FWHM, apex, retention time analysis
///   4. MSnSet quantitative proteomics data structures
///   5. MSnSet processing - log2, imputation, normalization
///   6. Peptide-to-protein summarization (sum/mean/median)
///   7. Sample QC metrics
///
/// Inspired by Bioconductor MSnbase (Gatto & Lilley 2012).

fn pad_right(s : String, n : Int) -> String {
  if s.length() >= n { s } else { s + " ".repeat(n - s.length()) }
}
fn pad_left(s : String, n : Int) -> String {
  if s.length() >= n { s } else { " ".repeat(n - s.length()) + s }
}
fn string_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n { return -1 }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched { return i }
    i = i + 1
  }
  -1
}
fn fmt_fixed(d : Double, decimals : Int) -> String {
  let s = d.to_string()
  let dot_idx = string_index_of(s, ".")
  if dot_idx >= 0 {
    let end = if dot_idx + 1 + decimals < s.length() { dot_idx + 1 + decimals } else { s.length() }
    s.substring(start=0, end=end)
  } else {
    let mut r = s + "."
    while r.length() < s.length() + 1 + decimals { r = r + "0" }
    r
  }
}

fn print_spectrum_info(label : String, sp : @bio.Spectrum) -> Unit {
  println("\{label}:")
  println("  peaks=\{sp.peaks_count()}, ms_level=...")
  println("  rt=" + fmt_fixed(sp.rt(), 2) + "s, TIC=" + fmt_fixed(sp.tic(), 2))
  let (bpi_mz, bpi_int) = sp.base_peak()
  println("  base peak: m/z=" + fmt_fixed(bpi_mz, 4) + ", intensity=" + fmt_fixed(bpi_int, 2))
  if sp.peaks_count() <= 12 {
    let mut line = "  mz = ["
    let mz_arr = sp.mz()
    for i = 0; i < mz_arr.length(); i = i + 1 {
      if i > 0 { line = line + ", " }
      line = line + fmt_fixed(mz_arr[i], 4)
    }
    println(line + "]")
    line = "  int= ["
    let int_arr = sp.intensity()
    for i = 0; i < int_arr.length(); i = i + 1 {
      if i > 0 { line = line + ", " }
      line = line + fmt_fixed(int_arr[i], 1)
    }
    println(line + "]")
  } else {
    println("  (Showing only first 6 / last 6 peaks)")
    let mut line = "  first 6 mz = ["
    for i in 0..<6 {
      if i > 0 { line = line + ", " }
      line = line + fmt_fixed(sp.mz()[i], 4)
    }
    println(line + "]")
    let total = sp.peaks_count()
    line = "  last 6  mz = ["
    for i in (total - 6)..<total {
      if i > total - 6 { line = line + ", " }
      line = line + fmt_fixed(sp.mz()[i], 4)
    }
    println(line + "]")
  }
}

///|
fn main {
  println("=== MSnbase Demo ===")
  println("Mass Spectrometry Data Structures & Processing")
  println("Inspired by Bioconductor MSnbase package")
  println("")

  // =====================================================================
  // Part 1: Spectrum (MS1 & MS2)
  // =====================================================================
  println("=== 1. Creating and inspecting Spectrum objects ===")

  // Build a simple synthetic MS1 spectrum (profile mode)
  let mz_vals : Array[Double] = [
    300.0510, 300.5422, 301.1033, 301.5512, 302.0601, 302.6020,
    303.0988, 303.6200, 304.0980, 304.5612, 305.0410, 305.5882,
    306.0100, 306.6008, 307.0455, 307.5830, 308.0291, 308.6100,
    309.0820, 309.5299, 310.0711, 310.5550, 311.0200, 311.6603,
  ]
  // Corresponding intensities with peaks at 302.06, 305.04, 308.6, 311.02
  let int_vals : Array[Double] = [
    120.0, 180.0, 250.0, 400.0, 980.0, 420.0,
    260.0, 180.0, 220.0, 380.0, 780.0, 390.0,
    200.0, 160.0, 150.0, 170.0, 210.0, 880.0,
    450.0, 260.0, 190.0, 220.0, 920.0, 460.0,
  ]
  let ms1_sp = @bio.Spectrum::new(
    mz_vals, int_vals,
    @bio.ms_level_ms1(),
    @bio.polarity_positive(),
    125.48,
  )
  println("")
  print_spectrum_info("  MS1 Spectrum (profile)", ms1_sp)

  // Find a peak by m/z with ppm tolerance
  let target_peak = 302.0580
  let idx_found = ms1_sp.find_peak(target_peak, ppm=50.0)
  println(
    "\n  Searching for m/z=\{target_peak} (±50 ppm): " +
    (if idx_found >= 0 {
      "found at index=\{idx_found}, " +
      "m/z=" + fmt_fixed(ms1_sp.mz()[idx_found], 4) +
      " (Δ=" + fmt_fixed((ms1_sp.mz()[idx_found] - target_peak) * 1000000.0 / target_peak, 1) + " ppm)"
    } else {
      "NOT FOUND"
    }),
  )

  // Build a simple MS2 spectrum with precursor
  let ms2_mz : Array[Double] = [
    110.0712, 148.1200, 218.1502, 286.1890,
    364.2150, 431.2588, 524.2660, 650.3410,
  ]
  let ms2_int : Array[Double] = [
    850.0, 420.0, 9800.0, 2100.0,
    1500.0, 680.0, 320.0, 180.0,
  ]
  let ms2_sp = @bio.Spectrum::new(
    ms2_mz, ms2_int,
    @bio.ms_level_ms2(),
    @bio.polarity_positive(),
    125.48,
  ).with_precursor(524.2655, 2)
  println("")
  print_spectrum_info("  MS2 Spectrum", ms2_sp)
  println("  precursor m/z=" + fmt_fixed(ms2_sp.precursor_mz, 4))
  println("  precursor charge=\{ms2_sp.precursor_charge}")

  println("\n=== 2. Spectrum Processing ===")

  // (a) Filter by m/z range
  let filtered = ms1_sp.filter_mz_range(301.0, 307.0)
  println("\n  a) Filter m/z [301, 307]:")
  println("     peaks before: \{ms1_sp.peaks_count()}  after: \{filtered.peaks_count()}")
  print_spectrum_info("     Filtered MS1", filtered)

  // (b) TIC normalization (sum of intensities = 1)
  let tic_norm = ms1_sp.normalize_tic()
  println("\n  b) TIC normalized spectrum (sum I = 1):")
  println("     original TIC = " + fmt_fixed(ms1_sp.tic(), 2) + ", new TIC = " + fmt_fixed(tic_norm.tic(), 4))
  // Verify sum of normalized intensities = 1
  let mut sum_check = 0.0
  for k = 0; k < tic_norm.intensity().length(); k = k + 1 {
    let v = tic_norm.intensity()[k]
    sum_check = sum_check + v
  }
  println("     intensity sum after normalize_tic: " + fmt_fixed(sum_check, 6))

  // (c) Moving-average smoothing
  let smoothed = ms1_sp.smooth_moving_average(half_window=2)
  // Print some comparison
  println("\n  c) Smoothed spectrum (MA, half_window=2):")
  println("     sample (index 2-7):")
  for i in 2..<8 {
    println("       mz=" + fmt_fixed(ms1_sp.mz()[i], 4) + " | orig=" + fmt_fixed(ms1_sp.intensity()[i], 1) +
      "  smooth=" + fmt_fixed(smoothed.intensity()[i], 1))
  }

  // (d) Baseline correction
  let baseline_corr = ms1_sp.baseline_correct_minwin(half_window=3)
  println("\n  d) Baseline-corrected (min window, hw=3):")
  for i in 2..<8 {
    println("       mz=" + fmt_fixed(ms1_sp.mz()[i], 4) + " | orig=" + fmt_fixed(ms1_sp.intensity()[i], 1) +
      "  bcorr=" + fmt_fixed(baseline_corr.intensity()[i], 1))
  }

  // (e) Centroiding (pick peaks above SNR)
  let centroided = ms1_sp.centroid_simple(snr_threshold=1.5)
  println("\n  e) Centroided (simple SNR threshold):")
  println("     profile peaks: \{ms1_sp.peaks_count()}  -> centroid peaks: \{centroided.peaks_count()}")
  let mut line = "     centroid m/z  = ["
  let cmz = centroided.mz()
  for i = 0; i < cmz.length(); i = i + 1 {
    if i > 0 { line = line + ", " }
    line = line + fmt_fixed(cmz[i], 4)
  }
  println(line + "]")
  line = "     centroid int. = ["
  let cint = centroided.intensity()
  for i = 0; i < cint.length(); i = i + 1 {
    if i > 0 { line = line + ", " }
    line = line + fmt_fixed(cint[i], 1)
  }
  println(line + "]")

  // =====================================================================
  // Part 2: Chromatogram (XIC)
  // =====================================================================
  println("\n=== 3. Creating and analyzing Chromatogram (XIC) ===")
  // Synthetic elution profile (Gaussian-like peak at ~125s)
  let rt_vals : Array[Double] = [
    118.0, 119.0, 120.0, 121.0, 122.0, 123.0,
    124.0, 124.5, 125.0, 125.5, 126.0, 127.0,
    128.0, 129.0, 130.0, 131.0, 132.0,
  ]
  let xic_int : Array[Double] = [
    200.0,  350.0,  600.0, 1200.0, 2400.0, 5000.0,
    9500.0, 12000.0, 13500.0, 12200.0, 9000.0, 4500.0,
    2000.0,  900.0,  400.0,  250.0,  200.0,
  ]
  let xic = @bio.Chromatogram::new(
    rt_vals, xic_int,
    mz_target=524.2655,
    ppm_tolerance=10.0,
    acquisition_mode="XIC",
  )
  println("  Chromatogram (XIC):")
  println("    points=\{xic.n_points()}, target m/z=" + fmt_fixed(xic.mz_target, 4))
  println("    acquisition mode: \{xic.acquisition_mode}")
  println("    total AUC (trapezoid) = " + fmt_fixed(xic.total_auc(), 1))
  let (apex_rt, apex_int) = xic.apex()
  println("    apex: rt=" + fmt_fixed(apex_rt, 2) + "s, intensity=" + fmt_fixed(apex_int, 1))
  println("    FWHM (approx.) = " + fmt_fixed(xic.fwhm(), 2) + "s")

  // =====================================================================
  // Part 3: MSnSet Quantitative Proteomics
  // =====================================================================
  println("\n=== 4. Building MSnSet Quantitative Proteomics Object ===")
  // 8 peptides (features) x 4 samples
  let peptide_names = [
    "Pep_A1_GeneX", "Pep_A2_GeneX", "Pep_A3_GeneX",
    "Pep_B1_GeneY", "Pep_B2_GeneY",
    "Pep_C1_GeneZ", "Pep_C2_GeneZ",
    "Pep_D1_GeneW",
  ]
  let sample_names_prot = ["Sample_1", "Sample_2", "Sample_3", "Sample_4"]
  // Protein mapping: each peptide -> its gene name
  // Rows: peptides, Cols: samples
  let prot_exprs : Array[Array[Double]] = [
    // Sample_1  Sample_2  Sample_3  Sample_4
    [1200000.0, 1350000.0, 1100000.0, 1500000.0], // A1 GeneX
    [ 900000.0,  980000.0,  850000.0, 1100000.0], // A2 GeneX
    [ 650000.0,       0.0,  600000.0,  720000.0], // A3 GeneX (missing Sample_2)
    [2100000.0, 2200000.0, 2000000.0, 2400000.0], // B1 GeneY
    [1400000.0, 1420000.0, 1350000.0, 1580000.0], // B2 GeneY
    [ 800000.0,  750000.0,  820000.0,  900000.0], // C1 GeneZ
    [ 550000.0,  520000.0,  580000.0,  620000.0], // C2 GeneZ
    [3200000.0, 3000000.0, 3100000.0, 3500000.0], // D1 GeneW
  ]
  // Build feature annotations with protein/gene mapping
  let fd_arr : Array[@bio.MSnFeatureData] = []
  for i = 0; i < peptide_names.length(); i = i + 1 {
    let pname = peptide_names[i]
    let prot = if i < 3 { "GeneX" } else { if i < 5 { "GeneY" } else { if i < 7 { "GeneZ" } else { "GeneW" } } }
    fd_arr.push(@bio.MSnFeatureData::new(pname, [], prot, 0, 0.0, 0.0))
  }
  // Sample annotations
  let sd_arr : Array[@bio.MSnSampleData] = []
  for j = 0; j < sample_names_prot.length(); j = j + 1 {
    let group = if j % 2 == 0 { "Control" } else { "Treatment" }
    sd_arr.push(@bio.MSnSampleData::new(sample_names_prot[j], group, "", 0, ""))
  }
  let msnset = @bio.MSnSet::new(prot_exprs, fd_arr, sd_arr)
  println("  MSnSet built:")
  println("    features (peptides): \{msnset.n_features()}")
  println("    samples: \{msnset.n_samples()}")
  println("    feature names: \{msnset.feature_names().join(", ")}")
  println("    sample names:  \{msnset.sample_names().join(", ")}")

  // Query rows and columns
  let feature_query = "Pep_A1_GeneX"
  let sample_query = "Sample_3"
  println("\n  Query feature '\{feature_query}':  \{msnset.get_feature(feature_query)}")
  let sample_vals = msnset.get_sample(sample_query)
  println("  Query sample '\{sample_query}':  \{sample_vals}")

  println("\n=== 5. MSnSet Processing Pipeline ===")

  // (a) Log2 transform
  let logged = msnset.log2_transform(offset=1.0)
  println("  a) Log2 (offset=1) applied. is_log2=\{logged.is_log2}")
  println("     '\{feature_query}' log2:")
  let mut line = "       ["
  let lf = logged.get_feature(feature_query)
  for i = 0; i < lf.length(); i = i + 1 {
    if i > 0 { line = line + ", " }
    line = line + fmt_fixed(lf[i], 4)
  }
  println(line + "]")

  // (b) Impute missing values (row mean)
  let imputed = msnset.impute_missing(method="mean")
  println("\n  b) Mean imputation of missing (zero/NaN):")
  println("     Feature '\{peptide_names[2]}' (A3 had missing in sample1):")
  println("       orig:    \{prot_exprs[2]}")
  line = "       imputed: ["
  let if2 = imputed.get_feature(peptide_names[2])
  for i = 0; i < if2.length(); i = i + 1 {
    if i > 0 { line = line + ", " }
    line = line + fmt_fixed(if2[i], 1)
  }
  println(line + "]")

  // (c) Sum normalization (TMT-style: each sample sums to 1e6)
  let sum_norm = msnset.normalize_sum(target=1000000.0)
  println("\n  c) Sum normalization (each sample col sums to 1e6):")
  // Compute column sums
  for j in 0..<sample_names_prot.length() {
    let mut cs = 0.0
    for i in 0..<msnset.n_features() {
      cs = cs + sum_norm.exprs()[i][j]
    }
    println("     sample \{sample_names_prot[j]} sum = " + cs.to_int().to_string())
  }

  // (d) Median centering normalization
  let med_norm = msnset.normalize_median_center()
  println("\n  d) Median-centering normalization applied.")

  println("\n=== 6. Peptide-to-Protein Summarization ===")
  println("  Peptides mapped to 4 proteins.")
  println("  3 methods: 'sum' (default), 'mean', 'median'")
  let methods_arr : Array[String] = ["sum", "mean", "median"]
  for midx = 0; midx < methods_arr.length(); midx = midx + 1 {
    let method = methods_arr[midx]
    let prots = msnset.summarize_to_proteins(method=method)
    println("")
    println("  --- method='\{method}' ---")
    println("    protein features: \{prots.n_features()} (\{prots.feature_names().join(", ")})")
    // Print each protein row
    let pnames = prots.feature_names()
    for pidx = 0; pidx < pnames.length(); pidx = pidx + 1 {
      let pname = pnames[pidx]
      let mut line2 = "    " + pad_right(pname, 8) + ": ["
      let pf = prots.get_feature(pname)
      for i = 0; i < pf.length(); i = i + 1 {
        if i > 0 { line2 = line2 + ", " }
        line2 = line2 + fmt_fixed(pf[i], 1)
      }
      println(line2 + "]")
    }
  }

  println("\n=== 7. Sample QC Metrics ===")
  let qc = msnset.sample_qc()
  println("  " + pad_right("sample", 18) + " | miss% |   mean   |  median  |   max")
  println("  " + "-".repeat(18 + 56))
  for idx = 0; idx < qc.length(); idx = idx + 1 {
    let (name, n_miss, mean_v, median_v, max_v) = qc[idx]
    let miss_pct = n_miss.to_double() / msnset.n_features().to_double() * 100.0
    println(
      "  " + pad_right(name, 18) + " | " + fmt_fixed(miss_pct, 1) + "% " +
      "| " + pad_left(fmt_fixed(mean_v, 2), 8) + " | " + pad_left(fmt_fixed(median_v, 2), 8) + " | " + pad_left(fmt_fixed(max_v, 1), 8),
    )
  }

  println("\n=== MSnbase Demo Complete ===")
}
