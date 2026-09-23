///|
/// scikit-allel `allel.stats.hw` — Hardy–Weinberg statistics at array
/// (per-variant) granularity: observed/expected heterozygosity and the
/// inbreeding coefficient.

///|
/// Observed heterozygosity per variant: `n_het / n_called` over non-missing
/// calls (`allel.heterozygosity_observed`).
pub fn skallel_heterozygosity_observed(
  g : Array[Array[Array[Int]]],
  fill : Double,
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for gi in g {
    let mut het = 0
    let mut called = 0
    for s in gi {
      let a = s[0]
      let b = s[1]
      if a >= 0 && b >= 0 {
        called = called + 1
        if a != b {
          het = het + 1
        }
      }
    }
    if called > 0 {
      out.push(het.to_double() / called.to_double())
    } else {
      out.push(fill)
    }
  }
  out
}

///|
/// Expected heterozygosity per variant under Hardy–Weinberg equilibrium:
/// `1 − Σ af^ploidy` (`allel.heterozygosity_expected`). `fill` is used where
/// allele frequencies do not sum to one.
pub fn skallel_heterozygosity_expected(
  af : Array[Array[Double]],
  ploidy : Int,
  fill : Double,
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for afi in af {
    let mut s = 0.0
    let mut total = 0.0
    for f in afi {
      s = s + @math.pow(f, ploidy.to_double())
      total = total + f
    }
    if total < 1.0 {
      out.push(fill)
    } else {
      out.push(1.0 - s)
    }
  }
  out
}

///|
/// Inbreeding coefficient per variant (`F = 1 − Ho/He`, `allel.inbreeding_coefficient`).
/// `fill` is used where the expected heterozygosity is zero (or no calls).
pub fn skallel_inbreeding_coefficient(
  g : Array[Array[Array[Int]]],
  fill : Double,
) -> Array[Double] {
  let mut max_allele = 0
  for gi in g {
    for s in gi {
      for a in s {
        if a > max_allele {
          max_allele = a
        }
      }
    }
  }
  let counts = skallel_gt_count_alleles(g, max_allele)
  let out : Array[Double] = Array::new()
  for i = 0; i < g.length(); i = i + 1 {
    let gi = g[i]
    let mut het = 0
    let mut called = 0
    for s in gi {
      let a = s[0]
      let b = s[1]
      if a >= 0 && b >= 0 {
        called = called + 1
        if a != b {
          het = het + 1
        }
      }
    }
    let ci = counts[i]
    let mut total = 0
    for c in ci {
      total = total + c
    }
    let mut he = 0.0
    if total > 0 {
      let mut s = 0.0
      for c in ci {
        if c > 0 {
          let f = c.to_double() / total.to_double()
          s = s + f * f
        }
      }
      he = 1.0 - s
    }
    if he > 0.0 {
      let ho = het.to_double() / called.to_double()
      out.push(1.0 - ho / he)
    } else {
      out.push(fill)
    }
  }
  out
}