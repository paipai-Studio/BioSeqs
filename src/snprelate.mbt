///|
/// SNPRelate - SNP genetic relationship analysis on genotype dosage data.
/// Inspired by the Bioconductor SNPRelate package (Zheng & Weir, 2016).
/// Implements SNP statistics, IBS/IBD-KING relatedness, PCA (EIGMIX-style
/// covariance eigen-analysis) and LD pruning. Input matrices hold genotype
/// dosages (rows = SNPs, columns = samples) with values 0, 1, 2 and -1 for
/// missing calls.

///|
/// Per-SNP summary statistics.
pub struct SnpVariantStats {
  allele_frequency : Array[Double]
  maf : Array[Double]
  missing_rate : Array[Double]
}

///|
/// Compute per-SNP allele frequency (alternate), minor allele frequency
/// and missing call rate.
pub fn snpgds_snp_stats(dosages : Array[Array[Int]]) -> SnpVariantStats {
  let n_snp = dosages.length()
  let af : Array[Double] = Array(capacity=n_snp)
  let maf : Array[Double] = Array(capacity=n_snp)
  let miss : Array[Double] = Array(capacity=n_snp)
  for row in dosages {
    let mut called = 0
    let mut alt = 0
    let mut missing = 0
    for dosage in row {
      if dosage < 0 {
        missing = missing + 1
      } else {
        called = called + 2
        alt = alt + dosage
      }
    }
    let p = if called > 0 { alt.to_double() / called.to_double() } else { 0.0 }
    af.push(p)
    maf.push(if p <= 0.5 { p } else { 1.0 - p })
    miss.push(
      if row.length() > 0 {
        missing.to_double() / row.length().to_double()
      } else {
        0.0
      },
    )
  }
  SnpVariantStats::{ allele_frequency: af, maf, missing_rate: miss, }
}

///|
/// Extract one sample's dosage column (across all SNPs).
fn snp_sample_column(dosages : Array[Array[Int]], sample : Int) -> Array[Int] {
  let col : Array[Int] = Array(capacity=dosages.length())
  for row in dosages {
    col.push(row[sample])
  }
  col
}

///|
/// Count per-sample heterozygous loci restricted to loci valid for the
/// given pair. Returns (sum of squared genotype differences, IBS0 count,
/// n valid loci, heterozygote count of sample i, heterozygote count of
/// sample j) following SNPRelate's KING robust accumulator.
fn snp_king_pair_accumulators(
  row_i : Array[Int],
  row_j : Array[Int],
) -> (Int, Int, Int, Int, Int) {
  let mut sum_sq = 0
  let mut ibs0 = 0
  let mut loci = 0
  let mut n_het_i = 0
  let mut n_het_j = 0
  let m = if row_i.length() < row_j.length() {
    row_i.length()
  } else {
    row_j.length()
  }
  let mut k = 0
  while k < m {
    let gi = row_i[k]
    let gj = row_j[k]
    if gi >= 0 && gj >= 0 {
      loci = loci + 1
      let diff = gi - gj
      sum_sq = sum_sq + diff * diff
      if diff == 2 || diff == -2 {
        ibs0 = ibs0 + 1
      }
      if gi == 1 {
        n_het_i = n_het_i + 1
      }
      if gj == 1 {
        n_het_j = n_het_j + 1
      }
    }
    k = k + 1
  }
  (sum_sq, ibs0, loci, n_het_i, n_het_j)
}

///|
/// KING-robust IBD results (Manichaikul et al. 2010), matching
/// snpgdsIBDKing with family=NULL: kinship uses the between-family form
/// phi = 0.5 - SumSq / (4 * min(N1_Aa, N2_Aa)), where the accumulators
/// run over loci valid for both samples.
pub struct SnpKingResult {
  /// Robust kinship coefficients (n x n; diagonal 0.5).
  kinship : Array[Array[Double]]
  /// Proportion of IBS0 loci (n x n; diagonal 0).
  ibs0 : Array[Array[Double]]
}

///|
/// Compute KING-robust kinship and IBS0 proportion matrices.
pub fn snpgds_ibd_king(dosages : Array[Array[Int]]) -> SnpKingResult {
  let n = if dosages.length() > 0 { dosages[0].length() } else { 0 }
  let columns : Array[Array[Int]] = Array(capacity=n)
  let mut s = 0
  while s < n {
    columns.push(snp_sample_column(dosages, s))
    s = s + 1
  }
  let kinship = snp_make_matrix(n, 0.0)
  let ibs0 = snp_make_matrix(n, 0.0)
  let mut i = 0
  while i < n {
    kinship[i][i] = 0.5
    let mut j = i + 1
    while j < n {
      let (sum_sq, ibs0_count, loci, n_het_i, n_het_j) = snp_king_pair_accumulators(
        columns[i],
        columns[j],
      )
      let _ = loci
      let min_het = if n_het_i < n_het_j { n_het_i } else { n_het_j }
      let phi = if min_het > 0 {
        0.5 - sum_sq.to_double() / (4.0 * min_het.to_double())
      } else {
        0.5
      }
      kinship[i][j] = phi
      kinship[j][i] = phi
      let rate = if loci > 0 {
        ibs0_count.to_double() / loci.to_double()
      } else {
        0.0
      }
      ibs0[i][j] = rate
      ibs0[j][i] = rate
      j = j + 1
    }
    i = i + 1
  }
  SnpKingResult::{ kinship, ibs0, }
}

///|
/// KING homozygosity IBD results (KING Homo method): IBD sharing k0/k1
/// estimates based on allele frequencies, matching snpgdsIBDKing with
/// method="homo".
pub struct SnpKingHomoResult {
  /// Probability of zero IBD sharing (k0).
  k0 : Array[Array[Double]]
  /// Probability of one allele IBD sharing (k1).
  k1 : Array[Array[Double]]
}

///|
/// Compute KING homo k0/k1 matrices. Requires per-SNP alternate allele
/// frequencies; pairs with zero frequency information yield 0.
pub fn snpgds_ibd_king_homo(
  dosages : Array[Array[Int]],
  allele_freqs : Array[Double],
) -> SnpKingHomoResult {
  let n = if dosages.length() > 0 { dosages[0].length() } else { 0 }
  let columns : Array[Array[Int]] = Array(capacity=n)
  let mut s = 0
  while s < n {
    columns.push(snp_sample_column(dosages, s))
    s = s + 1
  }
  let k0 = snp_make_matrix(n, 0.0)
  let k1 = snp_make_matrix(n, 0.0)
  let mut sum_afreq = 0.0
  let mut sum_afreq2 = 0.0
  for p in allele_freqs {
    sum_afreq = sum_afreq + p * (1.0 - p)
    sum_afreq2 = sum_afreq2 + p * p * (1.0 - p) * (1.0 - p)
  }
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let (sum_sq, ibs0, _loci, _ni, _nj) = snp_king_pair_accumulators(
        columns[i],
        columns[j],
      )
      let _ = _loci
      let _ = _ni
      let _ = _nj
      let denom1 = 8.0 * sum_afreq
      let theta = if denom1 > 0.0 {
        0.5 - sum_sq.to_double() / denom1
      } else {
        0.0
      }
      let denom2 = 2.0 * sum_afreq2
      let kk0 = if denom2 > 0.0 { ibs0.to_double() / denom2 } else { 0.0 }
      let kk1 = 2.0 - 2.0 * kk0 - 4.0 * theta
      k0[i][j] = kk0
      k0[j][i] = kk0
      k1[i][j] = kk1
      k1[j][i] = kk1
      j = j + 1
    }
    i = i + 1
  }
  SnpKingHomoResult::{ k0, k1, }
}

///|
/// IBS analysis results.
pub struct SnpIbsResult {
  /// IBS proportions (0.5*IBS1 + IBS2) / loci (n x n).
  ibs : Array[Array[Double]]
  /// Genetic distance = 1 - IBS proportion (n x n).
  distance : Array[Array[Double]]
}

///|
/// Compute IBS proportions and genetic distances between samples.
pub fn snpgds_ibs(dosages : Array[Array[Int]]) -> SnpIbsResult {
  let n = if dosages.length() > 0 { dosages[0].length() } else { 0 }
  let ibs = snp_make_matrix(n, 0.0)
  let dist = snp_make_matrix(n, 0.0)
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let mut loci = 0
      let mut ibs2 = 0
      let mut ibs1 = 0
      let mut k = 0
      while k < dosages.length() {
        let gi = dosages[k][i]
        let gj = dosages[k][j]
        if gi >= 0 && gj >= 0 {
          loci = loci + 1
          let diff = (gi - gj).abs()
          if diff == 0 {
            ibs2 = ibs2 + 1
          } else if diff == 1 {
            ibs1 = ibs1 + 1
          }
        }
        k = k + 1
      }
      let proportion = if loci > 0 {
        (0.5 * ibs1.to_double() + ibs2.to_double()) / loci.to_double()
      } else {
        0.0
      }
      ibs[i][j] = proportion
      ibs[j][i] = proportion
      dist[i][j] = 1.0 - proportion
      dist[j][i] = 1.0 - proportion
      j = j + 1
    }
    ibs[i][i] = 1.0
    i = i + 1
  }
  SnpIbsResult::{ ibs, distance: dist, }
}

///|
/// PCA results for SNP data.
pub struct SnpPcaResult {
  /// Leading eigenvalues (descending).
  eigenvalues : Array[Double]
  /// Sample loadings: eigenvector[k][sample].
  eigenvectors : Array[Array[Double]]
  /// Fraction of variance explained per component.
  var_explained : Array[Double]
}

///|
/// SNP PCA following SNPRelate's EIGMIX covariance method: genotypes are
/// centered at 2p, scaled by 1/sqrt(sum 4 p (1-p)) with missing calls
/// treated as the population mean, and the covariance matrix is
/// decomposed with a Jacobi eigen solver.
pub fn snpgds_pca(
  dosages : Array[Array[Int]],
  n_components : Int,
) -> SnpPcaResult {
  let stats = snpgds_snp_stats(dosages)
  let n_snp = dosages.length()
  let n = if n_snp > 0 { dosages[0].length() } else { 0 }
  if n == 0 || n_snp == 0 {
    return SnpPcaResult::{
      eigenvalues: [],
      eigenvectors: [],
      var_explained: [],
    }
  }
  let mut scale_sum = 0.0
  let mut k = 0
  while k < n_snp {
    let p = stats.allele_frequency[k]
    scale_sum = scale_sum + 4.0 * p * (1.0 - p)
    k = k + 1
  }
  let scale = if scale_sum > 0.0 { 1.0 / scale_sum.sqrt() } else { 0.0 }
  // centered scaled genotype matrix g[snp][sample]
  let scaled : Array[Array[Double]] = Array(capacity=n_snp)
  k = 0
  while k < n_snp {
    let two_p = 2.0 * stats.allele_frequency[k]
    let row : Array[Double] = Array(capacity=n)
    let mut s = 0
    while s < n {
      let g = dosages[k][s]
      row.push(if g >= 0 { (g.to_double() - two_p) * scale } else { 0.0 })
      s = s + 1
    }
    scaled.push(row)
    k = k + 1
  }
  // covariance matrix between samples
  let cov = snp_make_matrix(n, 0.0)
  let mut i = 0
  while i < n {
    let mut j = i
    while j < n {
      let mut acc = 0.0
      k = 0
      while k < n_snp {
        acc = acc + scaled[k][i] * scaled[k][j]
        k = k + 1
      }
      cov[i][j] = acc
      cov[j][i] = acc
      j = j + 1
    }
    i = i + 1
  }
  let count = if n_components < n { n_components } else { n }
  let (values, vectors) = snp_eigen_decomposition(cov, count)
  let mut total = 0.0
  for v in values {
    total = total + (if v > 0.0 { v } else { 0.0 })
  }
  let explained : Array[Double] = []
  for v in values {
    explained.push(if total > 0.0 { v / total } else { 0.0 })
  }
  SnpPcaResult::{
    eigenvalues: values,
    eigenvectors: vectors,
    var_explained: explained,
  }
}

///|
/// Pairwise linkage disequilibrium (method "corr": Pearson correlation
/// of genotype dosages over pairwise-complete calls).
pub fn snpgds_ld(
  dosages : Array[Array[Int]],
  snp_i : Int,
  snp_j : Int,
) -> Double raise SeqIOError {
  if snp_i < 0 ||
    snp_j < 0 ||
    snp_i >= dosages.length() ||
    snp_j >= dosages.length() {
    raise SeqIOError("LD SNP index out of range")
  }
  let row_i = dosages[snp_i]
  let row_j = dosages[snp_j]
  let mut n = 0
  let mut sum_i = 0.0
  let mut sum_j = 0.0
  let mut k = 0
  while k < row_i.length() && k < row_j.length() {
    let gi = row_i[k]
    let gj = row_j[k]
    if gi >= 0 && gj >= 0 {
      n = n + 1
      sum_i = sum_i + gi.to_double()
      sum_j = sum_j + gj.to_double()
    }
    k = k + 1
  }
  if n < 2 {
    return 0.0
  }
  let mean_i = sum_i / n.to_double()
  let mean_j = sum_j / n.to_double()
  let mut cov = 0.0
  let mut var_i = 0.0
  let mut var_j = 0.0
  k = 0
  while k < row_i.length() && k < row_j.length() {
    let gi = row_i[k]
    let gj = row_j[k]
    if gi >= 0 && gj >= 0 {
      let di = gi.to_double() - mean_i
      let dj = gj.to_double() - mean_j
      cov = cov + di * dj
      var_i = var_i + di * di
      var_j = var_j + dj * dj
    }
    k = k + 1
  }
  let denom = (var_i * var_j).sqrt()
  if denom > 0.0 {
    cov / denom
  } else {
    0.0
  }
}

///|
/// Result of LD pruning: indices of SNPs retained in file order.
pub struct SnpLdPruneResult {
  kept_indices : Array[Int]
  removed_indices : Array[Int]
}

///|
/// Greedy LD pruning following snpgdsLDpruning(method="corr"): scan SNPs
/// in file order; drop a candidate when its absolute genotype
/// correlation with any previously kept SNP within `slide_max` positions
/// exceeds `threshold`. SNPs with `maf <= maf_threshold` are skipped up
/// front.
pub fn snpgds_ld_pruning(
  dosages : Array[Array[Int]],
  threshold : Double,
  maf_threshold : Double,
  slide_max : Int,
) -> SnpLdPruneResult {
  let stats = snpgds_snp_stats(dosages)
  let kept : Array[Int] = []
  let removed : Array[Int] = []
  let n_snp = dosages.length()
  let mut i = 0
  while i < n_snp {
    let mut to_keep = true
    if stats.maf[i] <= maf_threshold {
      to_keep = false
    }
    if to_keep {
      let mut start = i - slide_max
      if start < 0 {
        start = 0
      }
      let mut j = start
      while j < i {
        let mut kept_j = false
        for kj in kept {
          if kj == j {
            kept_j = true
            break
          }
        }
        if kept_j {
          let r = snpgds_ld(dosages, i, j) catch { _ => 0.0 }
          if r.abs() > threshold {
            to_keep = false
            break
          }
        }
        j = j + 1
      }
    }
    if to_keep {
      kept.push(i)
    } else {
      removed.push(i)
    }
    i = i + 1
  }
  SnpLdPruneResult::{ kept_indices: kept, removed_indices: removed, }
}

///|
fn snp_make_matrix(n : Int, fill : Double) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    out.push(Array::make(n, fill))
    i = i + 1
  }
  out
}

///|
/// Jacobi eigen decomposition of a symmetric matrix. Returns the top
/// `n_components` eigenpairs sorted by descending eigenvalue, with
/// eigenvectors aligned to samples.
fn snp_eigen_decomposition(
  matrix : Array[Array[Double]],
  n_components : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let n = matrix.length()
  if n == 0 {
    return ([], [])
  }
  let a : Array[Array[Double]] = Array(capacity=n)
  let v : Array[Array[Double]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array(capacity=n)
    let ident : Array[Double] = Array(capacity=n)
    let mut j = 0
    while j < n {
      row.push(matrix[i][j])
      ident.push(if i == j { 1.0 } else { 0.0 })
      j = j + 1
    }
    a.push(row)
    v.push(ident)
    i = i + 1
  }
  let mut sweep = 0
  while sweep < 60 {
    let mut max_val = 0.0
    let mut max_i = 0
    let mut max_j = 1
    i = 0
    while i < n {
      let mut j = i + 1
      while j < n {
        if a[i][j].abs() > max_val {
          max_val = a[i][j].abs()
          max_i = i
          max_j = j
        }
        j = j + 1
      }
      i = i + 1
    }
    if max_val < 1.0e-10 {
      break
    }
    let theta = 0.5 *
      @math.atan2(2.0 * a[max_i][max_j], a[max_j][max_j] - a[max_i][max_i])
    let c = @math.cos(theta)
    let s = @math.sin(theta)
    // rotate: a = J^T A J, v = V J
    i = 0
    while i < n {
      let aik_i = a[i][max_i]
      let aik_j = a[i][max_j]
      a[i][max_i] = c * aik_i - s * aik_j
      a[i][max_j] = s * aik_i + c * aik_j
      let vik_i = v[i][max_i]
      let vik_j = v[i][max_j]
      v[i][max_i] = c * vik_i - s * vik_j
      v[i][max_j] = s * vik_i + c * vik_j
      i = i + 1
    }
    i = 0
    while i < n {
      let aki_i = a[max_i][i]
      let akj_i = a[max_j][i]
      a[max_i][i] = c * aki_i - s * akj_i
      a[max_j][i] = s * aki_i + c * akj_i
      i = i + 1
    }
    sweep = sweep + 1
  }
  // collect diagonal and sort descending
  let diag : Array[(Double, Int)] = Array(capacity=0)
  i = 0
  while i < n {
    diag.push((a[i][i], i))
    i = i + 1
  }
  let mut j = 0
  while j < diag.length() - 1 {
    let mut k = j + 1
    while k < diag.length() {
      if diag[j].0 < diag[k].0 {
        let tmp = diag[j]
        diag[j] = diag[k]
        diag[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }
  let count = if n_components < n { n_components } else { n }
  let values : Array[Double] = Array(capacity=count)
  let vectors : Array[Array[Double]] = Array(capacity=count)
  i = 0
  while i < count {
    values.push(diag[i].0)
    let idx = diag[i].1
    let vec : Array[Double] = Array(capacity=n)
    let mut s = 0
    while s < n {
      vec.push(v[s][idx])
      s = s + 1
    }
    vectors.push(vec)
    i = i + 1
  }
  (values, vectors)
}
