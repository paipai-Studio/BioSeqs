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
  let out : Array[Double] = Array(capacity=0)
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
  let out : Array[Double] = Array(capacity=0)
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
  let num : Array[Double] = Array(capacity=0)
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
  let out : Array[Double] = Array(capacity=0)
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
  let gn_norm : Array[Array[Double]] = Array(capacity=0)
  for v = 0; v < n_variants; v = v + 1 {
    let mut sum = 0.0
    for s = 0; s < n_samples; s = s + 1 {
      sum = sum + gn[v][s].to_double()
    }
    let mean = sum / n_samples.to_double()
    let p = mean / 2.0
    let scale = (2.0 * p * (1.0 - p)).sqrt()
    let row : Array[Double] = Array(capacity=0)
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
  let gram : Array[Array[Double]] = Array(capacity=0)
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
  let coords : Array[Array[Double]] = Array(capacity=0)
  let mut var_sum = 0.0
  for c = 0; c < nc; c = c + 1 {
    if evals[c] > 0.0 {
      var_sum = var_sum + evals[c]
    }
  }
  let ratio : Array[Double] = Array(capacity=0)
  for c = 0; c < nc; c = c + 1 {
    ratio.push(if var_sum > 0.0 { evals[c] / var_sum } else { 0.0 })
  }
  for s = 0; s < n_samples; s = s + 1 {
    let row : Array[Double] = Array(capacity=0)
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
  let ac : Array[Array[Int]] = Array(capacity=0)
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
  let a : Array[Array[Double]] = Array(capacity=0)
  let b : Array[Array[Double]] = Array(capacity=0)
  let c : Array[Array[Double]] = Array(capacity=0)
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
  let out : Array[Double] = Array(capacity=0)
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
}///|
/// scikit-allel `admixture.py` — unbiased allele-frequency estimators of
/// heterozygosity and the Patterson F2/F3/F4(D) statistics (Patterson et al.
/// 2012, Appendix A). All take biallelic allele-count arrays (n_variants × 2).

///|
fn alt_freq(ac_i : Array[Int]) -> Double {
  let an = ac_i[0] + ac_i[1]
  if an == 0 {
    0.0
  } else {
    ac_i[1].to_double() / an.to_double()
  }
}

///|
/// Unbiased estimator of `h`, where `2h` is the heterozygosity
/// (`allel.admixture.h_hat`).
pub fn skallel_h_hat(ac : Array[Array[Int]]) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for aci in ac {
    let an = aci[0] + aci[1]
    if an < 2 {
      out.push(0.0)
    } else {
      out.push((aci[0] * aci[1]).to_double() / (an * (an - 1)).to_double())
    }
  }
  out
}

///|
/// Unbiased estimator of F2(A, B), the branch length between two populations
/// (`allel.patterson_f2`).
pub fn skallel_patterson_f2(
  aca : Array[Array[Int]],
  acb : Array[Array[Int]],
) -> Array[Double] {
  let ha = skallel_h_hat(aca)
  let hb = skallel_h_hat(acb)
  let out : Array[Double] = Array(capacity=0)
  let n = if aca.length() < acb.length() { aca.length() } else { acb.length() }
  for i = 0; i < n; i = i + 1 {
    let sa = (aca[i][0] + aca[i][1]).to_double()
    let sb = (acb[i][0] + acb[i][1]).to_double()
    let a = alt_freq(aca[i])
    let b = alt_freq(acb[i])
    let d = a - b
    out.push(d * d - ha[i] / sa - hb[i] / sb)
  }
  out
}

///|
/// Unbiased estimator of F3(C; A, B). Returns `(T, B)`: per-variant
/// un-normalised f3 estimates and the heterozygosity in C (`allel.patterson_f3`).
pub fn skallel_patterson_f3(
  acc : Array[Array[Int]],
  aca : Array[Array[Int]],
  acb : Array[Array[Int]],
) -> (Array[Double], Array[Double]) {
  let hc = skallel_h_hat(acc)
  let t : Array[Double] = Array(capacity=0)
  let b : Array[Double] = Array(capacity=0)
  let mut n = if acc.length() < aca.length() { acc.length() } else { aca.length() }
  if n > acb.length() {
    n = acb.length()
  }
  for i = 0; i < n; i = i + 1 {
    let sc = (acc[i][0] + acc[i][1]).to_double()
    let c = alt_freq(acc[i])
    let a = alt_freq(aca[i])
    let b_ = alt_freq(acb[i])
    t.push((c - a) * (c - b_) - hc[i] / sc)
    b.push(2.0 * hc[i])
  }
  (t, b)
}

///|
/// Unbiased estimator of the four-population statistic D(A, B; C, D), a.k.a.
/// the ABBA-BABA test. Returns `(num, den)` per variant; the averaged D is
/// `Σnum / Σden` (`allel.patterson_d`).
pub fn skallel_patterson_d(
  aca : Array[Array[Int]],
  acb : Array[Array[Int]],
  acc : Array[Array[Int]],
  acd : Array[Array[Int]],
) -> (Array[Double], Array[Double]) {
  let num : Array[Double] = Array(capacity=0)
  let den : Array[Double] = Array(capacity=0)
  let mut n = aca.length()
  n = if n < acb.length() { n } else { acb.length() }
  n = if n < acc.length() { n } else { acc.length() }
  n = if n < acd.length() { n } else { acd.length() }
  for i = 0; i < n; i = i + 1 {
    let a = alt_freq(aca[i])
    let b = alt_freq(acb[i])
    let c = alt_freq(acc[i])
    let d = alt_freq(acd[i])
    num.push((a - b) * (c - d))
    den.push((a + b - 2.0 * a * b) * (c + d - 2.0 * c * d))
  }
  (num, den)
}

///|
/// Averaged D statistic over variants (the block estimator): `Σnum / Σden`.
pub fn skallel_patterson_d_average(
  aca : Array[Array[Int]],
  acb : Array[Array[Int]],
  acc : Array[Array[Int]],
  acd : Array[Array[Int]],
) -> Double {
  let (num, den) = skallel_patterson_d(aca, acb, acc, acd)
  let mut sn = 0.0
  let mut sd = 0.0
  for i = 0; i < num.length(); i = i + 1 {
    sn = sn + num[i]
    sd = sd + den[i]
  }
  if sd == 0.0 {
    0.0
  } else {
    sn / sd
  }
}

///|
/// scikit-allel `selection.py` (Garud H) + windowing — haplotype-diversity and
/// soft-sweep statistics, and a windowed Tajima's D wrapper.

///|
fn skallel_sort_strings_asc(arr : Array[String]) -> Array[String] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] > key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn skallel_sort_doubles_desc(arr : Array[Double]) -> Array[Double] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] < key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Frequencies of distinct haplotype columns in a haplotype array
/// (n_variants × n_haplotypes), sorted descending.
fn distinct_haplotype_frequencies(h : Array[Array[Int]]) -> Array[Double] {
  let n_variants = h.length()
  let n_hap = if n_variants > 0 { h[0].length() } else { 0 }
  if n_hap == 0 {
    return []
  }
  let haps : Array[String] = Array(capacity=0)
  for c = 0; c < n_hap; c = c + 1 {
    let sb = StringBuilder()
    for v = 0; v < n_variants; v = v + 1 {
      sb.write_char((h[v][c] + 48).to_uint16().unsafe_to_char())
    }
    haps.push(sb.to_string())
  }
  let sorted = skallel_sort_strings_asc(haps)
  let freqs : Array[Double] = Array(capacity=0)
  let mut i = 0
  while i < n_hap {
    let mut j = i
    while j < n_hap && sorted[j] == sorted[i] {
      j = j + 1
    }
    freqs.push((j - i).to_double() / n_hap.to_double())
    i = j
  }
  skallel_sort_doubles_desc(freqs)
}

///|
/// Haplotype diversity of a haplotype array (`allel.haplotype_diversity`): the
/// unbiased estimator `(1 − Σf²) · n/(n−1)` over distinct haplotype
/// frequencies.
pub fn skallel_haplotype_diversity(h : Array[Array[Int]]) -> Double {
  let n_hap = if h.length() > 0 { h[0].length() } else { 0 }
  if n_hap <= 1 {
    return 0.0
  }
  let f = distinct_haplotype_frequencies(h)
  let mut sum_sq = 0.0
  for fi in f {
    sum_sq = sum_sq + fi * fi
  }
  (1.0 - sum_sq) * n_hap.to_double() / (n_hap - 1).to_double()
}

///|
/// Garud et al. (2015) H1, H12, H123 and H2/H1 statistics for detecting soft
/// sweeps (`allel.garud_h`).
pub fn skallel_garud_h(h : Array[Array[Int]]) -> (Double, Double, Double, Double) {
  let f = distinct_haplotype_frequencies(h)
  let mut h1 = 0.0
  for fi in f {
    h1 = h1 + fi * fi
  }
  let f0 = if f.length() > 0 { f[0] } else { 0.0 }
  let f1 = if f.length() > 1 { f[1] } else { 0.0 }
  let f2 = if f.length() > 2 { f[2] } else { 0.0 }
  let mut tail2 = 0.0
  for i = 2; i < f.length(); i = i + 1 {
    tail2 = tail2 + f[i] * f[i]
  }
  let mut tail3 = 0.0
  for i = 3; i < f.length(); i = i + 1 {
    tail3 = tail3 + f[i] * f[i]
  }
  let h12 = (f0 + f1) * (f0 + f1) + tail2
  let h123 = (f0 + f1 + f2) * (f0 + f1 + f2) + tail3
  let h2 = h1 - f0 * f0
  (h1, h12, h123, if h1 > 0.0 { h2 / h1 } else { 0.0 })
}

///|
/// Equal-sized index windows `[start, stop)` over `n` items, stepping by
/// `step`; a trailing partial window is dropped.
pub fn skallel_index_windows(n : Int, size : Int, step : Int) -> Array[(Int, Int)] {
  let out : Array[(Int, Int)] = Array(capacity=0)
  let mut ws = 0
  while ws + size <= n {
    out.push((ws, ws + size))
    ws = ws + step
  }
  out
}

///|
/// Tajima's D computed over non-overlapping (or stepped) index windows of
/// biallelic allele counts, reusing `popgen_tajima_d` per window. Returns one
/// D value per window.
pub fn skallel_windowed_tajima_d(
  ac : Array[Array[Int]],
  n_haplotypes : Int,
  size : Int,
  step : Int,
) -> Array[Double] {
  let pi = skallel_mean_pairwise_difference(ac)
  let out : Array[Double] = Array(capacity=0)
  let n = ac.length()
  let mut ws = 0
  while ws + size <= n {
    let mut s = 0
    let mut pi_sum = 0.0
    for i = ws; i < ws + size; i = i + 1 {
      pi_sum = pi_sum + pi[i]
      if ac[i].length() >= 2 && ac[i][1] > 0 && ac[i][0] > 0 {
        s = s + 1
      }
    }
    let res = popgen_tajima_d(pi_sum, s, n_haplotypes, size)
    out.push(res.statistic)
    ws = ws + step
  }
  out
}///|
/// Hudson's Fst (block estimator) over an index window of variants, from
/// precomputed per-variant π and Dxy.
fn hudson_fst_block_from(
  dxy : Array[Double],
  pia : Array[Double],
  pib : Array[Double],
  ws : Int,
  size : Int,
) -> Double {
  let mut sn = 0.0
  let mut sd = 0.0
  for v = ws; v < ws + size; v = v + 1 {
    let num = dxy[v] - (pia[v] + pib[v]) / 2.0
    sn = sn + num
    sd = sd + dxy[v]
  }
  if sd == 0.0 {
    0.0
  } else {
    sn / sd
  }
}

///|
/// Population branching statistic (PBS) comparing allele frequencies among
/// three populations, in index windows (`allel.pbs`). Uses Hudson's Fst.
pub fn skallel_pbs(
  ac1 : Array[Array[Int]],
  ac2 : Array[Array[Int]],
  ac3 : Array[Array[Int]],
  window_size : Int,
  step : Int,
  normed : Bool,
) -> Array[Double] {
  let pi1 = skallel_mean_pairwise_difference(ac1)
  let pi2 = skallel_mean_pairwise_difference(ac2)
  let pi3 = skallel_mean_pairwise_difference(ac3)
  let dxy12 = skallel_mean_pairwise_difference_between(ac1, ac2)
  let dxy13 = skallel_mean_pairwise_difference_between(ac1, ac3)
  let dxy23 = skallel_mean_pairwise_difference_between(ac2, ac3)
  let mut n = ac1.length()
  if ac2.length() < n {
    n = ac2.length()
  }
  if ac3.length() < n {
    n = ac3.length()
  }
  let out : Array[Double] = Array(capacity=0)
  let mut ws = 0
  while ws + window_size <= n {
    let f12 = hudson_fst_block_from(dxy12, pi1, pi2, ws, window_size)
    let f13 = hudson_fst_block_from(dxy13, pi1, pi3, ws, window_size)
    let f23 = hudson_fst_block_from(dxy23, pi2, pi3, ws, window_size)
    let c12 = if f12 < 0.0 { 0.0 } else if f12 > 0.99999 { 0.99999 } else { f12 }
    let c13 = if f13 < 0.0 { 0.0 } else if f13 > 0.99999 { 0.99999 } else { f13 }
    let c23 = if f23 < 0.0 { 0.0 } else if f23 > 0.99999 { 0.99999 } else { f23 }
    let t12 = -@math.ln(1.0 - c12)
    let t13 = -@math.ln(1.0 - c13)
    let t23 = -@math.ln(1.0 - c23)
    let p = (t12 + t13 - t23) / 2.0
    if normed {
      let norm = 1.0 + (t12 + t13 + t23) / 2.0
      out.push(p / norm)
    } else {
      out.push(p)
    }
    ws = ws + step
  }
  out
}
///|
/// Patterson F-statistic between two populations (`allel.patterson_fst`).
/// Returns the per-variant `(num, den)`; the Fst estimate is `Σnum / Σden`.
pub fn skallel_patterson_fst(
  aca : Array[Array[Int]],
  acb : Array[Array[Int]],
) -> (Array[Double], Array[Double]) {
  let num = skallel_patterson_f2(aca, acb)
  let ha = skallel_h_hat(aca)
  let hb = skallel_h_hat(acb)
  let den : Array[Double] = Array(capacity=0)
  let nt = if aca.length() < acb.length() { aca.length() } else { acb.length() }
  for i = 0; i < nt; i = i + 1 {
    den.push(num[i] + ha[i] + hb[i])
  }
  (num, den)
}
