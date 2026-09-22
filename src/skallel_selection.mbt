///|
/// scikit-allel `selection.py` — the EHH / iHS / nSL / standardization pipeline
/// for detecting selection, ported from the `allel.opt.stats` scanning
/// algorithms (shared-suffix-length streaming scans).
///
/// Note: where scikit-allel returns NaN (e.g., EHH never decays below
/// `min_ehh` before the contig edge, or a missing value), this port returns
/// 0.0.

///|
/// Reverse the variant (row) order of a haplotype array.
fn reverse_variants(h : Array[Array[Int]]) -> Array[Array[Int]] {
  let n = h.length()
  let out : Array[Array[Int]] = Array::new()
  let mut i = n - 1
  while i >= 0 {
    out.push(h[i])
    i = i - 1
  }
  out
}

///|
fn reverse_doubles(a : Array[Double]) -> Array[Double] {
  let n = a.length()
  let out : Array[Double] = Array::new()
  let mut i = n - 1
  while i >= 0 {
    out.push(a[i])
    i = i - 1
  }
  out
}

///|
/// Extended haplotype homozygosity (EHH) decay away from the first variant
/// (`allel.ehh_decay`). `ehh[k]` is the fraction of haplotype pairs sharing a
/// prefix of at least `k + 1` variants.
pub fn skallel_ehh_decay(h : Array[Array[Int]]) -> Array[Double] {
  let n_variants = h.length()
  let n_hap = if n_variants > 0 { h[0].length() } else { 0 }
  let n_pairs = n_hap * (n_hap - 1) / 2
  let ehh = Array::make(n_variants, 0.0)
  if n_pairs == 0 {
    return ehh
  }
  let b = Array::make(n_variants + 1, 0)
  for j = 0; j < n_hap; j = j + 1 {
    for k = j + 1; k < n_hap; k = k + 1 {
      let mut spl = 0
      while spl < n_variants && h[spl][j] == h[spl][k] {
        spl = spl + 1
      }
      b[spl] = b[spl] + 1
    }
  }
  let mut suffix = 0
  for t = n_variants; t >= 1; t = t - 1 {
    suffix = suffix + b[t]
    ehh[t - 1] = suffix.to_double() / n_pairs.to_double()
  }
  ehh
}

///|
fn ssl2ihh(
  ssl : Array[Int],
  n_pairs : Int,
  l_max : Int,
  variant_idx : Int,
  gaps : Array[Double],
  min_ehh : Double,
  include_edges : Bool,
) -> Double {
  if n_pairs == 0 {
    return 0.0
  }
  let hh_breaks = Array::make(variant_idx + 1, 0)
  for j = 0; j < n_pairs; j = j + 1 {
    let l = ssl[j]
    if l <= variant_idx {
      hh_breaks[l] = hh_breaks[l] + 1
    }
  }
  let mut n_pairs_ident = n_pairs - hh_breaks[0]
  let mut ehh_prv = n_pairs_ident.to_double() / n_pairs.to_double()
  if ehh_prv <= min_ehh {
    return 0.0
  }
  let mut ihh = 0.0
  for i = 1; i <= variant_idx; i = i + 1 {
    n_pairs_ident = n_pairs_ident - hh_breaks[i]
    let ehh_cur = n_pairs_ident.to_double() / n_pairs.to_double()
    let gap = gaps[variant_idx - i]
    if gap < 0.0 {
      return 0.0
    }
    ihh = ihh + gap * (ehh_cur + ehh_prv) / 2.0
    if ehh_cur <= min_ehh {
      return ihh
    }
    ehh_prv = ehh_cur
  }
  if include_edges {
    ihh
  } else {
    0.0
  }
}

///|
/// Integrated haplotype homozygosity (iHH) scanning backwards from each
/// variant, splitting reference (0) and alternate (1) haplotype pairs
/// (scikit-allel's `ihh01_scan` forward pass).
fn ihh01_scan(
  h : Array[Array[Int]],
  gaps : Array[Double],
  min_ehh : Double,
  min_maf : Double,
  include_edges : Bool,
) -> (Array[Double], Array[Double]) {
  let n_variants = h.length()
  let n_hap = if n_variants > 0 { h[0].length() } else { 0 }
  let n_pairs = n_hap * (n_hap - 1) / 2
  let ssl = Array::make(n_pairs, 0)
  let ssl00 = Array::make(n_pairs, 0)
  let ssl11 = Array::make(n_pairs, 0)
  let vstat0 = Array::make(n_variants, 0.0)
  let vstat1 = Array::make(n_variants, 0.0)

  for i = 0; i < n_variants; i = i + 1 {
    let mut u = 0
    let mut u00 = 0
    let mut u11 = 0
    let mut c0 = 0
    let mut c1 = 0
    let mut l_max_00 = 0
    let mut l_max_11 = 0
    for j = 0; j < n_hap; j = j + 1 {
      let a1 = h[i][j]
      if a1 == 0 {
        c0 = c0 + 1
      } else if a1 == 1 {
        c1 = c1 + 1
      }
      for k = j + 1; k < n_hap; k = k + 1 {
        let a2 = h[i][k]
        let l = if a1 < 0 || a2 < 0 {
          ssl[u] + 1
        } else if a1 == 0 && a2 == 0 {
          ssl[u] + 1
        } else if a1 == 1 && a2 == 1 {
          ssl[u] + 1
        } else {
          0
        }
        ssl[u] = l
        if a1 == 0 && a2 == 0 {
          ssl00[u00] = l
          u00 = u00 + 1
          if l > l_max_00 {
            l_max_00 = l
          }
        } else if a1 == 1 && a2 == 1 {
          ssl11[u11] = l
          u11 = u11 + 1
          if l > l_max_11 {
            l_max_11 = l
          }
        }
        u = u + 1
      }
    }
    let maf = if c0 + c1 == 0 {
      0.0
    } else if c0 < c1 {
      c0.to_double() / (c0 + c1).to_double()
    } else {
      c1.to_double() / (c0 + c1).to_double()
    }
    if maf < min_maf {
      vstat0[i] = 0.0
      vstat1[i] = 0.0
    } else {
      vstat0[i] = ssl2ihh(ssl00, u00, l_max_00, i, gaps, min_ehh, include_edges)
      vstat1[i] = ssl2ihh(ssl11, u11, l_max_11, i, gaps, min_ehh, include_edges)
    }
  }
  (vstat0, vstat1)
}

///|
/// Unstandardized integrated haplotype score (iHS) for each variant
/// (`allel.ihs`), comparing iHH between the reference and alternate alleles.
pub fn skallel_ihs(
  h : Array[Array[Int]],
  pos : Array[Int],
  min_ehh : Double,
  include_edges : Bool,
) -> Array[Double] {
  let gaps : Array[Double] = Array::new()
  for k = 0; k < pos.length() - 1; k = k + 1 {
    gaps.push((pos[k + 1] - pos[k]).to_double())
  }
  let (fwd0, fwd1) = ihh01_scan(h, gaps, min_ehh, 0.0, include_edges)
  let (rev0, rev1) = ihh01_scan(reverse_variants(h), reverse_doubles(gaps), min_ehh, 0.0, include_edges)
  let rev0r = reverse_doubles(rev0)
  let rev1r = reverse_doubles(rev1)
  let out = Array::make(h.length(), 0.0)
  for i = 0; i < h.length(); i = i + 1 {
    let ihh0 = fwd0[i] + rev0r[i]
    let ihh1 = fwd1[i] + rev1r[i]
    out[i] = if ihh0 > 0.0 && ihh1 > 0.0 {
      @math.ln(ihh1 / ihh0)
    } else {
      0.0
    }
  }
  out
}

///|
/// Number of segregating sites by length scan (scikit-allel's `nsl01_scan`
/// forward pass), returning average shared-suffix lengths for reference (0)
/// and alternate (1) haplotype pairs.
fn nsl01_scan(h : Array[Array[Int]]) -> (Array[Double], Array[Double]) {
  let n_variants = h.length()
  let n_hap = if n_variants > 0 { h[0].length() } else { 0 }
  let n_pairs = n_hap * (n_hap - 1) / 2
  let ssl = Array::make(n_pairs, 0)
  let vstat0 = Array::make(n_variants, 0.0)
  let vstat1 = Array::make(n_variants, 0.0)

  for i = 0; i < n_variants; i = i + 1 {
    let mut u = 0
    let mut ssl00_sum = 0
    let mut ssl11_sum = 0
    let mut u00 = 0
    let mut u11 = 0
    for j = 0; j < n_hap; j = j + 1 {
      let a1 = h[i][j]
      for k = j + 1; k < n_hap; k = k + 1 {
        let a2 = h[i][k]
        let l = if a1 < 0 || a2 < 0 {
          ssl[u] + 1
        } else if a1 == 0 && a2 == 0 {
          ssl[u] + 1
        } else if a1 == 1 && a2 == 1 {
          ssl[u] + 1
        } else {
          0
        }
        ssl[u] = l
        if a1 >= 0 && a2 >= 0 {
          if a1 == 0 && a2 == 0 {
            ssl00_sum = ssl00_sum + l
            u00 = u00 + 1
          } else if a1 == 1 && a2 == 1 {
            ssl11_sum = ssl11_sum + l
            u11 = u11 + 1
          }
        }
        u = u + 1
      }
    }
    vstat0[i] = if u00 > 0 {
      ssl00_sum.to_double() / u00.to_double()
    } else {
      0.0
    }
    vstat1[i] = if u11 > 0 {
      ssl11_sum.to_double() / u11.to_double()
    } else {
      0.0
    }
  }
  (vstat0, vstat1)
}

///|
/// Unstandardized number of segregating sites by length (nSL) for each
/// variant, after Ferrer-Admetlla et al. (2014) (`allel.nsl`).
pub fn skallel_nsl(h : Array[Array[Int]]) -> Array[Double] {
  let (fwd0, fwd1) = nsl01_scan(h)
  let (rev0, rev1) = nsl01_scan(reverse_variants(h))
  let rev0r = reverse_doubles(rev0)
  let rev1r = reverse_doubles(rev1)
  let out = Array::make(h.length(), 0.0)
  for i = 0; i < h.length(); i = i + 1 {
    let nsl0 = fwd0[i] + rev0r[i]
    let nsl1 = fwd1[i] + rev1r[i]
    out[i] = if nsl0 > 0.0 && nsl1 > 0.0 {
      @math.ln(nsl1 / nsl0)
    } else {
      0.0
    }
  }
  out
}

///|
/// Sort an integer array in place ascending (copying).
fn sort_ints_asc_copy(arr : Array[Int]) -> Array[Int] {
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
/// Make similar-sized allele-count bins (scikit-allel's
/// `make_similar_sized_bins`).
fn make_similar_sized_bins(x : Array[Int], n : Int) -> Array[Int] {
  let y = sort_ints_asc_copy(x)
  let bins : Array[Int] = Array::new()
  if y.length() == 0 {
    return bins
  }
  bins.push(y[0])
  let step = y.length() / n
  if step > 0 {
    let mut i = step
    while i < y.length() {
      let v = y[i]
      if v > bins[bins.length() - 1] {
        bins.push(v)
      }
      i = i + step
    }
  }
  if bins[bins.length() - 1] != y[y.length() - 1] {
    bins.push(y[y.length() - 1])
  }
  bins
}

///|
/// Centre and scale to unit variance (`allel.standardize`).
pub fn skallel_standardize(score : Array[Double]) -> Array[Double] {
  let n = score.length()
  if n == 0 {
    return []
  }
  let mut mean = 0.0
  for x in score {
    mean = mean + x
  }
  mean = mean / n.to_double()
  let mut variance = 0.0
  for x in score {
    let d = x - mean
    variance = variance + d * d
  }
  let std = (variance / n.to_double()).sqrt()
  let out : Array[Double] = Array::new()
  for x in score {
    if std > 0.0 {
      out.push((x - mean) / std)
    } else {
      out.push(0.0)
    }
  }
  out
}

///|
/// Standardize scores within alternate-allele-count bins
/// (`allel.standardize_by_allele_count`).
pub fn skallel_standardize_by_allele_count(
  score : Array[Double],
  aac : Array[Int],
  n_bins : Int,
) -> Array[Double] {
  let bins = make_similar_sized_bins(aac, n_bins)
  let out : Array[Double] = Array::make(score.length(), 0.0)
  if bins.length() < 2 {
    return out
  }
  for bi = 0; bi < bins.length() - 1; bi = bi + 1 {
    let x1 = bins[bi]
    let x2 = bins[bi + 1]
    // collect scores in this bin
    let mut mean = 0.0
    let mut cnt = 0
    for i = 0; i < score.length(); i = i + 1 {
      let in_bin = if bi == bins.length() - 2 {
        aac[i] >= x1
      } else if bi == 0 {
        aac[i] < x2
      } else {
        aac[i] >= x1 && aac[i] < x2
      }
      if in_bin {
        mean = mean + score[i]
        cnt = cnt + 1
      }
    }
    if cnt == 0 {
      continue
    }
    mean = mean / cnt.to_double()
    let mut variance = 0.0
    for i = 0; i < score.length(); i = i + 1 {
      let in_bin = if bi == bins.length() - 2 {
        aac[i] >= x1
      } else if bi == 0 {
        aac[i] < x2
      } else {
        aac[i] >= x1 && aac[i] < x2
      }
      if in_bin {
        let d = score[i] - mean
        variance = variance + d * d
      }
    }
    let std = (variance / cnt.to_double()).sqrt()
    for i = 0; i < score.length(); i = i + 1 {
      let in_bin = if bi == bins.length() - 2 {
        aac[i] >= x1
      } else if bi == 0 {
        aac[i] < x2
      } else {
        aac[i] >= x1 && aac[i] < x2
      }
      if in_bin {
        out[i] = if std > 0.0 { (score[i] - mean) / std } else { 0.0 }
      }
    }
  }
  out
}///|
/// Integrated haplotype homozygosity (iHH) scanning backwards from each
/// variant over the whole sample (no reference/alternate split), matching
/// scikit-allel's `ihh_scan`.
pub fn skallel_ihh_scan(
  h : Array[Array[Int]],
  gaps : Array[Double],
  min_ehh : Double,
  include_edges : Bool,
) -> Array[Double] {
  let n_variants = h.length()
  let n_hap = if n_variants > 0 { h[0].length() } else { 0 }
  let n_pairs = n_hap * (n_hap - 1) / 2
  let ssl = Array::make(n_pairs, 0)
  let vihh = Array::make(n_variants, 0.0)
  for i = 0; i < n_variants; i = i + 1 {
    let mut u = 0
    let mut l_max = 0
    for j = 0; j < n_hap; j = j + 1 {
      let a1 = h[i][j]
      for k = j + 1; k < n_hap; k = k + 1 {
        let a2 = h[i][k]
        let l = if a1 != a2 && a1 >= 0 && a2 >= 0 {
          0
        } else {
          ssl[u] + 1
        }
        ssl[u] = l
        u = u + 1
        if l > l_max {
          l_max = l
        }
      }
    }
    vihh[i] = ssl2ihh(ssl, n_pairs, l_max, i, gaps, min_ehh, include_edges)
  }
  vihh
}