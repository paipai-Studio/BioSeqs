///|
/// speckle propeller demo: differential proportions analysis for
/// single-cell cell types across biological replicates, including
/// getTransformedProps (logit / arcsin transforms), limma moderated
/// t-tests / ANOVA with robust empirical Bayes variance shrinkage,
/// median-library-size normalisation (normCounts) and beta-binomial
/// parameter estimation (estimateBetaParam(FromCounts)).
///
/// Run: moon run examples/propeller_demo
fn main raise {
  println("=== speckle propeller Demo ===")
  println("")

  // ---- Cell-level data: 3 cell types, 4 samples, 2 groups --------------
  // Counts table: rows are cell types (lexicographic order), columns
  // samples s1..s4; s1,s2 are controls and s3,s4 are treated.
  let types = ["B cell", "Monocyte", "T cell"]
  let samples2 = ["s1", "s2", "s3", "s4"]
  let tab2 = [
    [180.0, 200.0, 150.0, 160.0],
    [100.0, 120.0, 240.0, 210.0],
    [320.0, 300.0, 210.0, 230.0],
  ]
  let (clusters, sample_cells) = expand(tab2, types, samples2)
  let groups = []
  for cell in 0..<clusters.length() {
    groups.push(
      if sample_cells[cell] == "s1" || sample_cells[cell] == "s2" {
        "ctrl"
      } else {
        "treat"
      },
    )
  }

  println("--- 1. getTransformedProps: counts and proportions ---")
  let pl = @bio.get_transformed_props(clusters, sample_cells, transform="logit")
  let header = [lpad("type", 10)]
  for s in pl.sample_names {
    header.push(cpad(s, 9))
  }
  println("  " + join(header))
  for i in 0..<pl.cluster_names.length() {
    let line = [lpad(pl.cluster_names[i], 10)]
    for j in 0..<pl.sample_names.length() {
      line.push(cpad(f0(pl.counts.get(i, j)), 9))
    }
    println("  " + join(line))
  }
  println("  proportions:")
  for i in 0..<pl.cluster_names.length() {
    let line = [lpad(pl.cluster_names[i], 10)]
    for j in 0..<pl.sample_names.length() {
      line.push(cpad(f4(pl.proportions.get(i, j)), 9))
    }
    println("  " + join(line))
  }
  let pasin = @bio.get_transformed_props(
    clusters,
    sample_cells,
    transform="asin",
  )
  println(
    "  arcsin-sqrt transformed (T cell, s1): " + f4(pasin.transformed.get(2, 0)),
  )
  println("")

  println("--- 2. Moderated t-test: propeller_ttest ---")
  // Samples x groups group-membership design; columns are the
  // lexicographically sorted group names ["ctrl", "treat"], so the
  // [1, -1] contrast compares ctrl against treat.
  let design2 = @bio.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [0.0, 1.0],
  ])
  let trows = @bio.propeller_ttest(pl, design2, [1.0, -1.0])
  println("  cluster     ctrl     treat    ratio     t          p          FDR")
  for r in trows {
    println(
      "  " +
      lpad(r.cluster, 10) +
      " " +
      rpad(f3(r.means[0]), 8) +
      rpad(f3(r.means[1]), 9) +
      rpad(f4(r.prop_ratio), 8) +
      rpad(sf3(r.statistic), 11) +
      rpad(e3(r.p_value), 11) +
      e3(r.fdr),
    )
  }
  println("")

  println("--- 3. One-call wrapper (baseline proportions, sorted by p) ---")
  let rows = @bio.propeller(clusters, sample_cells, groups, transform="logit")
  for r in rows {
    let rr = match r.prop_ratio {
      Some(v) => f4(v)
      None => "-"
    }
    println(
      "  " +
      lpad(r.cluster, 10) +
      " [" +
      r.kind +
      "] baseline=" +
      f3(r.baseline_prop) +
      "  ratio=" +
      rr +
      "  stat=" +
      sf3(r.statistic) +
      "  p=" +
      e3(r.p_value) +
      "  FDR=" +
      e3(r.fdr),
    )
  }
  println("")

  // ---- Three groups: moderated ANOVA across 6 samples ------------------
  println("--- 4. Moderated ANOVA: ctrl vs drugA vs drugB ---")
  let samples3 = ["r1", "r2", "r3", "r4", "r5", "r6"]
  let tab3 = [
    [150.0, 140.0, 180.0, 170.0, 120.0, 130.0],
    [90.0, 100.0, 200.0, 190.0, 320.0, 300.0],
    [300.0, 310.0, 220.0, 240.0, 160.0, 170.0],
  ]
  let (clusters3, cells3) = expand(tab3, types, samples3)
  let pl3 = @bio.get_transformed_props(clusters3, cells3)
  // Samples r1,r2 -> ctrl; r3,r4 -> drugA; r5,r6 -> drugB.
  let design3 = @bio.Matrix::from_2d(6, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
  ])
  let arows = @bio.propeller_anova(pl3, design3)
  println("  cluster     ctrl     drugA    drugB    F          p          FDR")
  for r in arows {
    println(
      "  " +
      lpad(r.cluster, 10) +
      " " +
      rpad(f3(r.means[0]), 8) +
      rpad(f3(r.means[1]), 9) +
      rpad(f3(r.means[2]), 9) +
      rpad(f3(r.statistic), 11) +
      rpad(e3(r.p_value), 11) +
      e3(r.fdr),
    )
  }
  println("")

  // ---- normCounts and beta-binomial parameter estimation ---------------
  println("--- 5. normCounts: median library size normalisation ---")
  let counts = @bio.Matrix::from_2d(3, 4, tab2)
  let nc = @bio.propeller_norm_counts(counts)
  let raw = ["  Monocyte raw : "]
  for v in tab2[1] {
    raw.push(rpad(f2(v), 10))
  }
  println(join(raw))
  let normline = ["  Monocyte norm: "]
  for j in 0..<4 {
    normline.push(rpad(f2(nc.get(1, j)), 10))
  }
  println(join(normline))
  let lnc = @bio.propeller_norm_counts(counts, log=true)
  let logline = ["  Monocyte log2: "]
  for j in 0..<4 {
    logline.push(rpad(f4(lnc.get(1, j)), 10))
  }
  println(join(logline))
  println("")

  println("--- 6. Beta parameters of proportions ---")
  let bp = @bio.estimate_beta_params_from_counts(counts)
  for i in 0..<types.length() {
    println(
      "  " +
      lpad(types[i], 10) +
      " alpha=" +
      rpad(f3(bp.alpha[i]), 9) +
      " beta=" +
      rpad(f3(bp.beta[i]), 9) +
      " pi=" +
      rpad(f3(bp.pi[i]), 7) +
      " dispersion=" +
      f4(bp.dispersion[i]),
    )
  }
  // Method of moments for one numeric proportion vector (sample s1).
  let props_s1 = [tab2[0][0] / 600.0, tab2[1][0] / 600.0, tab2[2][0] / 600.0]
  let m = @bio.estimate_beta_param(props_s1)
  println("  s1 proportions: alpha=" + f3(m.alpha) + " beta=" + f3(m.beta))
}

///|
/// Expand a counts table (types x samples) into cell-level cluster and
/// sample vectors, as required by get_transformed_props / propeller.
fn expand(
  tab : Array[Array[Double]],
  type_names : Array[String],
  sample_names : Array[String],
) -> (Array[String], Array[String]) {
  let clusters = []
  let samples = []
  for i in 0..<type_names.length() {
    for j in 0..<sample_names.length() {
      let count = tab[i][j].to_int()
      for _ in 0..<count {
        clusters.push(type_names[i])
        samples.push(sample_names[j])
      }
    }
  }
  (clusters, samples)
}

///|
fn join(parts : Array[String]) -> String {
  let mut out = ""
  for p in parts {
    out += p
  }
  out
}

///|
fn rpad(s : String, width : Int) -> String {
  let mut out = s
  while out.length() < width {
    out += " "
  }
  out
}

///|
fn lpad(s : String, width : Int) -> String {
  let mut out = s
  while out.length() < width {
    out = " " + out
  }
  out
}

///|
fn cpad(s : String, width : Int) -> String {
  let mut out = s
  let mut left = true
  while out.length() < width {
    if left {
      out = " " + out
    } else {
      out += " "
    }
    left = !left
  }
  out
}

///|
/// Fixed-point rendering with ``d`` decimals (d up to 4).
fn fd(x : Double, d : Int) -> String {
  let neg = x < 0.0
  let scale = if d == 0 {
    1.0
  } else if d == 1 {
    10.0
  } else if d == 2 {
    100.0
  } else if d == 3 {
    1000.0
  } else {
    10000.0
  }
  let v = x.abs() * scale + 0.5
  let scaled = v.floor().to_int()
  let whole = scaled / scale.to_int()
  let frac = scaled % scale.to_int()
  let sign = if neg && (whole > 0 || frac > 0) { "-" } else { "" }
  if d == 0 {
    return sign + whole.to_string()
  }
  let mut fs = frac.to_string()
  while fs.length() < d {
    fs = "0" + fs
  }
  sign + whole.to_string() + "." + fs
}

///|
fn f0(x : Double) -> String {
  fd(x, 0)
}

///|
fn f2(x : Double) -> String {
  fd(x, 2)
}

///|
fn f3(x : Double) -> String {
  fd(x, 3)
}

///|
fn f4(x : Double) -> String {
  fd(x, 4)
}

///|
/// Signed fixed-point with 3 decimals.
fn sf3(x : Double) -> String {
  if x > 0.0 {
    "+" + f3(x)
  } else {
    f3(x)
  }
}

///|
/// Scientific notation with 3-digit mantissa, e.g. 8.112e-10.
fn e3(x : Double) -> String {
  if x <= 0.0 {
    return "0"
  }
  let mut exp = 0
  let mut v = x
  while v >= 10.0 {
    v /= 10.0
    exp += 1
  }
  while v < 1.0 {
    v *= 10.0
    exp -= 1
  }
  f3(v) + "e" + exp.to_string()
}
