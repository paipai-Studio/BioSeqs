///|
/// scikit-allel port — population-genetics statistics derived from VCF allele
/// counts (`ac`, per variant × allele) and biallelic genotype-dosage matrices
/// (`gn`, per variant × sample, 0/1/2 with −1 = missing): mean pairwise
/// difference (π / Dxy), Hudson's Fst, Rogers–Huff linkage disequilibrium (r),
/// and Patterson-scaled principal components analysis.

///|
pub struct SkallelPca {
  coords : Array[Array[Double]]
  explained_variance_ratio : Array[Double]
}

///|
fn ac_sum(ac : Array[Int]) -> Int {
  let mut s = 0
  for a in ac {
    s = s + a
  }
  s
}

///|
/// Mean pairwise difference (per variant) within a population, a.k.a. π
/// (`allel.mean_pairwise_difference`). Returns 0.0 for variants with fewer
/// than two haplotypes (scikit-allel's NaN fill).
pub fn skallel_mean_pairwise_difference(ac : Array[Array[Int]]) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for aci in ac {
    let an = ac_sum(aci)
    if an < 2 {
      out.push(0.0)
      continue
    }
    let n_pairs = (an * (an - 1)).to_double() / 2.0
    let mut n_same = 0.0
    for ai in aci {
      n_same = n_same + (ai * (ai - 1)).to_double() / 2.0
    }
    out.push((n_pairs - n_same) / n_pairs)
  }
  out
}

///|
/// Mean pairwise difference between two populations, a.k.a. Dxy
/// (`allel.mean_pairwise_difference_between`). Returns 0.0 where either
/// population has no haplotypes.
pub fn skallel_mean_pairwise_difference_between(
  ac1 : Array[Array[Int]],
  ac2 : Array[Array[Int]],
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  let n = if ac1.length() < ac2.length() { ac1.length() } else { ac2.length() }
  for i = 0; i < n; i = i + 1 {
    let an1 = ac_sum(ac1[i])
    let an2 = ac_sum(ac2[i])
    if an1 <= 0 || an2 <= 0 {
      out.push(0.0)
      continue
    }
    let n_pairs = (an1 * an2).to_double()
    let mut n_same = 0
    let m = if ac1[i].length() < ac2[i].length() {
      ac1[i].length()
    } else {
      ac2[i].length()
    }
    for j = 0; j < m; j = j + 1 {
      n_same = n_same + ac1[i][j] * ac2[i][j]
    }
    out.push((n_pairs - n_same.to_double()) / n_pairs)
  }
  out
}

///|
/// Hudson's estimator of Fst: returns `(num, den)` per variant
/// (`allel.hudson_fst`), where `num = Dxy − (π1 + π2)/2` and `den = Dxy`. The
/// per-variant Fst is `num/den`; the block estimator is `Σnum/Σden`.
pub fn skallel_hudson_fst(
  ac1 : Array[Array[Int]],
  ac2 : Array[Array[Int]],
) -> (Array[Double], Array[Double]) {
  let pi1 = skallel_mean_pairwise_difference(ac1)
  let pi2 = skallel_mean_pairwise_difference(ac2)
  let dxy = skallel_mean_pairwise_difference_between(ac1, ac2)
  let num : Array[Double] = Array::new()
  for i = 0; i < dxy.length(); i = i + 1 {
    num.push(dxy[i] - (pi1[i] + pi2[i]) / 2.0)
  }
  (num, dxy)
}

///|
/// Hudson's Fst averaged over variants (the block estimator): `Σnum / Σden`.
pub fn skallel_hudson_fst_block(
  ac1 : Array[Array[Int]],
  ac2 : Array[Array[Int]],
) -> Double {
  let (num, den) = skallel_hudson_fst(ac1, ac2)
  let mut sn = 0.0
  let mut sd = 0.0
  for i = 0; i < num.length(); i = i + 1 {
    sn = sn + num[i]
    sd = sd + den[i]
  }
  sn / sd
}

///|
fn pearson_r_pairwise(x : Array[Int], y : Array[Int]) -> Double {
  let mut sx = 0.0
  let mut sy = 0.0
  let mut n = 0
  for k = 0; k < x.length(); k = k + 1 {
    if x[k] != -1 && y[k] != -1 {
      sx = sx + x[k].to_double()
      sy = sy + y[k].to_double()
      n = n + 1
    }
  }
  if n < 2 {
    return 0.0
  }
  let mx = sx / n.to_double()
  let my = sy / n.to_double()
  let mut sxx = 0.0
  let mut syy = 0.0
  let mut sxy = 0.0
  for k = 0; k < x.length(); k = k + 1 {
    if x[k] != -1 && y[k] != -1 {
      let dx = x[k].to_double() - mx
      let dy = y[k].to_double() - my
      sxx = sxx + dx * dx
      syy = syy + dy * dy
      sxy = sxy + dx * dy
    }
  }
  let den = (sxx * syy).sqrt()
  if den == 0.0 {
    return 0.0
  }
  sxy / den
}

///|
/// Rogers–Huff linkage disequilibrium: signed correlation coefficient `r` for
/// every pair of variants, in condensed form (upper triangle, row-major)
/// (`allel.rogers_huff_r`). Missing genotypes (−1) are excluded pairwise.
pub fn skallel_rogers_huff_r(gn : Array[Array[Int]]) -> Array[Double] {
  let out : Array[Double] = Array::new()
  let n = gn.length()
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      out.push(pearson_r_pairwise(gn[i], gn[j]))
    }
  }
  out
}

///|
/// Principal components analysis of biallelic genotype dosages with the
/// Patterson et al. (2006) scaling (`allel.pca`). Input `gn` is
/// (n_variants × n_samples) with no missing data. Returns sample coordinates
/// (n_samples × n_components) and the explained-variance ratio per component.
pub fn skallel_pca(gn : Array[Array[Int]], n_components : Int) -> SkallelPca {
  let n_variants = gn.length()
  let n_samples = if n_variants > 0 { gn[0].length() } else { 0 }
  if n_samples == 0 {
    return SkallelPca::{ coords: [], explained_variance_ratio: [] }
  }

  // Patterson scaling: x' = (gn − 2p) / sqrt(2·p·(1−p)), p = mean(gn)/2
  let gn_norm : Array[Array[Double]] = Array::new()
  for v = 0; v < n_variants; v = v + 1 {
    let mut sum = 0.0
    for s = 0; s < n_samples; s = s + 1 {
      sum = sum + gn[v][s].to_double()
    }
    let mean = sum / n_samples.to_double()
    let p = mean / 2.0
    let scale = (2.0 * p * (1.0 - p)).sqrt()
    let row : Array[Double] = Array::new()
    for s = 0; s < n_samples; s = s + 1 {
      if scale > 0.0 {
        row.push((gn[v][s].to_double() - mean) / scale)
      } else {
        row.push(0.0)
      }
    }
    gn_norm.push(row)
  }

  // Gram matrix (n_samples × n_samples): G[i][j] = Σ_v gn_norm[v][i]·gn_norm[v][j]
  let gram : Array[Array[Double]] = Array::new()
  for i = 0; i < n_samples; i = i + 1 {
    let row = Array::make(n_samples, 0.0)
    for j = 0; j < n_samples; j = j + 1 {
      let mut acc = 0.0
      for v = 0; v < n_variants; v = v + 1 {
        acc = acc + gn_norm[v][i] * gn_norm[v][j]
      }
      row[j] = acc
    }
    gram.push(row)
  }

  let (eigenvalues, eigenvectors) = linalg_jacobi(gram, n_samples, 100 * n_samples * n_samples)
  let (evals, evecs) = linalg_sort_eigen(eigenvalues, eigenvectors, n_samples)

  let nc = if n_components < n_samples { n_components } else { n_samples }
  let coords : Array[Array[Double]] = Array::new()
  let mut var_sum = 0.0
  for c = 0; c < nc; c = c + 1 {
    if evals[c] > 0.0 {
      var_sum = var_sum + evals[c]
    }
  }
  let ratio : Array[Double] = Array::new()
  for c = 0; c < nc; c = c + 1 {
    ratio.push(if var_sum > 0.0 { evals[c] / var_sum } else { 0.0 })
  }
  for s = 0; s < n_samples; s = s + 1 {
    let row : Array[Double] = Array::new()
    for c = 0; c < nc; c = c + 1 {
      let scale = if evals[c] > 0.0 { evals[c].sqrt() } else { 0.0 }
      row.push(evecs[s][c] * scale)
    }
    coords.push(row)
  }

  SkallelPca::{ coords, explained_variance_ratio: ratio }
}///|
fn wc_het_count(gv : Array[Array[Int]], subpops : Array[Array[Int]], d : Int) -> Int {
  let mut cnt = 0
  for p = 0; p < subpops.length(); p = p + 1 {
    for si in subpops[p] {
      let a1 = gv[si][0]
      let a2 = gv[si][1]
      if a1 == d && a2 != d && a2 >= 0 {
        cnt = cnt + 1
      } else if a2 == d && a1 != d && a1 >= 0 {
        cnt = cnt + 1
      }
    }
  }
  cnt
}

///|
/// Per-variant Weir & Cockerham (1984) variance components `(a, b, c)` for one
/// variant, following scikit-allel's `_weir_cockerham_fst`.
fn wc_fst_variant(
  gv : Array[Array[Int]],
  subpops : Array[Array[Int]],
  max_allele : Int,
) -> (Array[Double], Array[Double], Array[Double]) {
  let n_alleles = max_allele + 1
  let r = subpops.length()

  // ac[d][p] = count of allele d in subpopulation p
  let ac : Array[Array[Int]] = Array::new()
  for d = 0; d < n_alleles; d = d + 1 {
    ac.push(Array::make(r, 0))
  }
  let an = Array::make(r, 0)
  for p = 0; p < r; p = p + 1 {
    for si in subpops[p] {
      for k = 0; k < 2; k = k + 1 {
        let allele = gv[si][k]
        if allele >= 0 && allele < n_alleles {
          ac[allele][p] = ac[allele][p] + 1
          an[p] = an[p] + 1
        }
      }
    }
  }

  let n = Array::make(r, 0)
  let mut n_total = 0
  for p = 0; p < r; p = p + 1 {
    n[p] = an[p] / 2
    n_total = n_total + n[p]
  }

  let a = Array::make(n_alleles, 0.0)
  let b = Array::make(n_alleles, 0.0)
  let c = Array::make(n_alleles, 0.0)
  if n_total <= 0 || r < 2 {
    return (a, b, c)
  }

  let n_bar = n_total.to_double() / r.to_double()
  let mut sum_n2 = 0
  for p = 0; p < r; p = p + 1 {
    sum_n2 = sum_n2 + n[p] * n[p]
  }
  let n_C =
    (n_total.to_double() - sum_n2.to_double() / n_total.to_double()) /
    (r - 1).to_double()

  let mut an_total = 0
  for p = 0; p < r; p = p + 1 {
    an_total = an_total + an[p]
  }

  for d = 0; d < n_alleles; d = d + 1 {
    let mut ac_total_d = 0
    for p = 0; p < r; p = p + 1 {
      ac_total_d = ac_total_d + ac[d][p]
    }
    let p_bar = ac_total_d.to_double() / an_total.to_double()

    let mut ss_num = 0.0
    for p = 0; p < r; p = p + 1 {
      let p_dp = if an[p] > 0 {
        ac[d][p].to_double() / an[p].to_double()
      } else {
        0.0
      }
      let diff = p_dp - p_bar
      ss_num = ss_num + n[p].to_double() * diff * diff
    }
    let s_squared = ss_num / (n_bar * (r - 1).to_double())

    let h_bar = wc_het_count(gv, subpops, d).to_double() / n_total.to_double()

    a[d] = (n_bar / n_C) *
      (
        s_squared -
        (1.0 / (n_bar - 1.0)) *
        (p_bar * (1.0 - p_bar) - (r - 1).to_double() * s_squared / r.to_double() - h_bar / 4.0)
      )
    b[d] = (n_bar / (n_bar - 1.0)) *
      (
        p_bar * (1.0 - p_bar) -
        (r - 1).to_double() * s_squared / r.to_double() -
        ((2.0 * n_bar - 1.0) * h_bar) / (4.0 * n_bar)
      )
    c[d] = h_bar / 2.0
  }

  (a, b, c)
}

///|
/// Weir & Cockerham (1984) variance components for every variant, from a
/// diploid genotype array (n_variants × n_samples × 2, −1 = missing) and a
/// list of subpopulation sample indices (`allel.weir_cockerham_fst`). Returns
/// `(a, b, c)` each of shape (n_variants × n_alleles).
pub fn skallel_weir_cockerham_fst(
  g : Array[Array[Array[Int]]],
  subpops : Array[Array[Int]],
  max_allele : Int,
) -> (Array[Array[Double]], Array[Array[Double]], Array[Array[Double]]) {
  let a : Array[Array[Double]] = Array::new()
  let b : Array[Array[Double]] = Array::new()
  let c : Array[Array[Double]] = Array::new()
  for v = 0; v < g.length(); v = v + 1 {
    let (av, bv, cv) = wc_fst_variant(g[v], subpops, max_allele)
    a.push(av)
    b.push(bv)
    c.push(cv)
  }
  (a, b, c)
}

///|
/// Per-variant Weir & Cockerham Fst, averaging variance components over
/// alleles: `Σa / (Σa + Σb + Σc)`. Returns 0.0 where the denominator is zero
/// (scikit-allel's NaN).
pub fn skallel_weir_cockerham_fst_variant(
  g : Array[Array[Array[Int]]],
  subpops : Array[Array[Int]],
  max_allele : Int,
) -> Array[Double] {
  let (a, b, c) = skallel_weir_cockerham_fst(g, subpops, max_allele)
  let out : Array[Double] = Array::new()
  for v = 0; v < a.length(); v = v + 1 {
    let mut sa = 0.0
    let mut sb = 0.0
    let mut sc = 0.0
    for d = 0; d < a[v].length(); d = d + 1 {
      sa = sa + a[v][d]
      sb = sb + b[v][d]
      sc = sc + c[v][d]
    }
    let den = sa + sb + sc
    out.push(if den == 0.0 { 0.0 } else { sa / den })
  }
  out
}

///|
/// Weir & Cockerham Fst averaged over all variants and alleles (block
/// estimator): `Σa / (Σa + Σb + Σc)`.
pub fn skallel_weir_cockerham_fst_block(
  g : Array[Array[Array[Int]]],
  subpops : Array[Array[Int]],
  max_allele : Int,
) -> Double {
  let (a, b, c) = skallel_weir_cockerham_fst(g, subpops, max_allele)
  let mut sa = 0.0
  let mut sb = 0.0
  let mut sc = 0.0
  for v = 0; v < a.length(); v = v + 1 {
    for d = 0; d < a[v].length(); d = d + 1 {
      sa = sa + a[v][d]
      sb = sb + b[v][d]
      sc = sc + c[v][d]
    }
  }
  let den = sa + sb + sc
  if den == 0.0 {
    0.0
  } else {
    sa / den
  }
}