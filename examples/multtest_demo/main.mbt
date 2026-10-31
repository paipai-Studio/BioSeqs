///|
/// multtest demo: marginal p-value adjustment (nine procedures from
/// multtest::mt.rawp2adjp) plus Westfall-Young permutation maxT / minP
/// adjusted p-values (multtest::mt.maxT / mt.minP engines).
///
/// Run with: moon run examples/multtest_demo
fn round3(x : Double) -> Double {
  (x * 1000.0).round() / 1000.0
}

///|
fn main raise {
  println("== marginal raw-p adjustment (GWAS-style 12 SNPs) ==")
  let rawp = [
    0.0002, 0.0011, 0.004, 0.009, 0.015, 0.028, 0.06, 0.12, 0.23, 0.41, 0.68, 0.95,
  ]
  let res = @bio.mt_rawp2adjp(rawp~, alpha=[0.05, 0.10])
  println("sorted p-values:")
  let mut row = ""
  for p in res.rawp {
    row += " " + round3(p).to_string()
  }
  println(row)
  let names = [
    "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH_0.05",
    "TSBH_0.1",
  ]
  println("")
  println("adjusted p-values and numbers rejected at alpha = 0.05:")
  for name in names {
    let col = res.adjp.get(name).unwrap()
    let (r, _) = @bio.mt_reject(adjp=col, alpha=0.05)
    let mut line = ""
    line += name
    while line.length() < 12 {
      line += " "
    }
    line += ":"
    for x in col {
      line += " " + round3(x).to_string()
    }
    line += "   reject=" + r.to_string()
    println(line)
  }
  println("")
  println(
    "ABH estimated number of true nulls: " + round3(res.h0_abh).to_string(),
  )

  // --- Westfall-Young permutation adjustment ------------------------------
  // m = 4 hypotheses, B = 8 permutations of the (absolute) test statistic.
  println("")
  println("== Westfall-Young permutation adjustment (B = 8) ==")
  let observed = [0.75, 0.85, 0.90, 0.88]
  let null_stats = [
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
    [0.05, 0.15, 0.25, 0.95, 0.05, 0.15, 0.25, 0.35],
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.9, 0.2, 0.2],
  ]
  let adj_maxt = @bio.mt_max_t(observed~, null_stats~)
  let adj_minp = @bio.mt_min_p(observed~, null_stats~)
  let mut obs_line = "observed t :"
  for t in observed {
    obs_line += " " + t.to_string()
  }
  println(obs_line)
  let mut maxt_line = "maxT adjp  :"
  for x in adj_maxt {
    maxt_line += " " + round3(x).to_string()
  }
  println(maxt_line)
  let mut minp_line = "minP adjp  :"
  for x in adj_minp {
    minp_line += " " + round3(x).to_string()
  }
  println(minp_line)
  let (r_maxt, _) = @bio.mt_reject(adjp=adj_maxt, alpha=0.5)
  let (r_minp, _) = @bio.mt_reject(adjp=adj_minp, alpha=0.5)
  println(
    "rejected at 0.5: maxT=" +
    r_maxt.to_string() +
    " minP=" +
    r_minp.to_string(),
  )
}
