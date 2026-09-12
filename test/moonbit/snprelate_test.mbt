///|
/// Tests for SNPRelate-style analyses. The dosage matrix is a 4-SNP x
/// 4-sample toy dataset; expected values were computed with a Python
/// mirror implementing the SNPRelate formulas.

///|
fn snp_test_dosages() -> Array[Array[Int]] {
  [
    [0, 0, 2, 2], // SNP0, p=0.5
    [0, 1, 1, 2], // SNP1, p=0.5
    [1, 1, 0, 1], // SNP2, p=0.375
    [-1, 0, 1, 2], // SNP3, p=0.5, one missing call
  ]
}

///|
fn snp_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
test "snpgds_snp_stats" {
  let stats = @src.snpgds_snp_stats(snp_test_dosages())
  assert_eq(stats.allele_frequency.length(), 4)
  assert_true(snp_close(stats.allele_frequency[0], 0.5, 1.0e-12))
  assert_true(snp_close(stats.allele_frequency[2], 0.375, 1.0e-12))
  assert_true(snp_close(stats.maf[2], 0.375, 1.0e-12))
  // SNP0 is p=0.5 exactly; SNP3 has one missing call of four
  assert_true(snp_close(stats.maf[0], 0.5, 1.0e-12))
  assert_true(snp_close(stats.missing_rate[3], 0.25, 1.0e-12))
  assert_eq(stats.missing_rate[0], 0.0)
}

///|
test "snpgds_ibd_king_pairwise" {
  let king = @src.snpgds_ibd_king(snp_test_dosages())
  assert_eq(king.kinship.length(), 4)
  // diagonal is 0.5
  assert_true(snp_close(king.kinship[0][0], 0.5, 1.0e-12))
  // pair (0,1): SumSq=1, min het=1 -> phi = 0.5 - 1/4 = 0.25
  assert_true(snp_close(king.kinship[0][1], 0.25, 1.0e-12))
  // pair (0,2): SumSq=6, min het=1 -> phi = 0.5 - 1.5 = -1.0
  assert_true(snp_close(king.kinship[0][2], -1.0, 1.0e-12))
  // pair (1,2): SumSq=6, min het=2 -> phi = 0.5 - 0.75 = -0.25
  assert_true(snp_close(king.kinship[1][2], -0.25, 1.0e-12))
  // symmetric
  assert_true(snp_close(king.kinship[2][3], king.kinship[3][2], 0.0))
  // IBS0 rates: pair (0,2) has 1 IBS0 locus of 3 valid
  assert_true(snp_close(king.ibs0[0][2], 1.0 / 3.0, 1.0e-12))
  assert_true(snp_close(king.ibs0[1][2], 0.25, 1.0e-12))
}

///|
test "snpgds_ibd_king_parent_offspring" {
  // A parent-offspring pair never shows IBS0 and KING kinship is ~0.25.
  // Rows are SNPs, columns are samples. Parent is het at every locus;
  // child inherits one allele per locus; samples 2/3 are homozygous
  // reference / alternate at every locus.
  let dosages = [
    [1, 1, 0, 2], // parent, child, hom-ref, hom-alt
    [1, 0, 0, 2],
    [1, 2, 0, 2],
    [1, 1, 0, 2],
    [1, 1, 0, 2],
    [1, 0, 0, 2],
  ]
  let king = @src.snpgds_ibd_king(dosages)
  // parent-child: SumSq = 0+1+1+0+0+1 = 3, parent het=6, child het=3 ->
  // min het = 3 -> phi = 0.5 - 3/12 = 0.25 (the KING parent-offspring
  // expectation)
  assert_true(snp_close(king.kinship[0][1], 0.25, 1.0e-12))
  assert_true(snp_close(king.ibs0[0][1], 0.0, 1.0e-12))
  // unrelated hom-ref vs hom-alt: all IBS0, SumSq = 6*4 = 24, min het = 0
  assert_true(snp_close(king.kinship[2][3], 0.5, 1.0e-12))
  assert_true(snp_close(king.ibs0[2][3], 1.0, 1.0e-12))
}

///|
test "snpgds_ibd_king_homo" {
  let stats = @src.snpgds_snp_stats(snp_test_dosages())
  let king = @src.snpgds_ibd_king_homo(
    snp_test_dosages(),
    stats.allele_frequency,
  )
  assert_eq(king.k0.length(), 4)
  // mirror: sum p(1-p) = 0.984375, sum p^2(1-p)^2 = 0.242431640625
  // pair (0,1): SumSq=1, IBS0=0 -> theta = 0.5 - 1/7.875, k0 = 0
  let theta01 = 0.5 - 1.0 / (8.0 * 0.984375)
  let k1_01 = 2.0 - 4.0 * theta01
  assert_true(snp_close(king.k0[0][1], 0.0, 1.0e-12))
  assert_true(snp_close(king.k1[0][1], k1_01, 1.0e-12))
  // pair (0,2): SumSq=6, IBS0=1 -> k0 = 1/(2*0.242431640625)
  let k0_02 = 1.0 / (2.0 * 0.242431640625)
  let theta02 = 0.5 - 6.0 / (8.0 * 0.984375)
  assert_true(snp_close(king.k0[0][2], k0_02, 1.0e-12))
  assert_true(
    snp_close(king.k1[0][2], 2.0 - 2.0 * k0_02 - 4.0 * theta02, 1.0e-12),
  )
}

///|
test "snpgds_ibs_proportions" {
  let ibs = @src.snpgds_ibs(snp_test_dosages())
  assert_eq(ibs.ibs.length(), 4)
  // pair (0,1): 2 IBS2 + 1 IBS1 over 3 loci -> (0.5 + 2)/3
  assert_true(snp_close(ibs.ibs[0][1], (0.5 + 2.0) / 3.0, 1.0e-12))
  // pair (1,2): 1 IBS2 + 2 IBS1 over 4 loci -> (0.5*2 + 1)/4
  assert_true(snp_close(ibs.ibs[1][2], 2.0 / 4.0, 1.0e-12))
  // distance = 1 - proportion
  assert_true(snp_close(ibs.distance[0][1], 1.0 - (0.5 + 2.0) / 3.0, 1.0e-12))
  assert_true(snp_close(ibs.distance[2][3], 0.375, 1.0e-12))
  // diagonal IBS proportion is 1
  assert_true(snp_close(ibs.ibs[1][1], 1.0, 1.0e-12))
}

///|
test "snpgds_pca_planted_structure" {
  // two sample clusters with opposite genotype structure: PC1 must
  // separate them and eigenvalues must be descending
  let dosages = [
    [0, 0, 0, 2, 2, 2],
    [0, 0, 1, 2, 2, 1],
    [1, 0, 0, 1, 2, 2],
    [2, 1, 0, 0, 0, 1],
  ]
  let pca = @src.snpgds_pca(dosages, 3)
  assert_eq(pca.eigenvalues.length(), 3)
  assert_eq(pca.eigenvectors.length(), 3)
  for k in 0..<2 {
    assert_true(pca.eigenvalues[k] >= pca.eigenvalues[k + 1])
  }
  let mut total = 0.0
  for v in pca.var_explained {
    total = total + v
  }
  assert_true(snp_close(total, 1.0, 1.0e-9))
  // PC1 loadings: same sign within cluster, opposite between clusters
  let pc1 = pca.eigenvectors[0]
  let first_cluster = pc1[0] + pc1[1] + pc1[2]
  let second_cluster = pc1[3] + pc1[4] + pc1[5]
  assert_true(first_cluster * second_cluster < 0.0)
  assert_true(pca.var_explained[0] > 0.4)
}

///|
test "snpgds_pca_handles_missing" {
  let dosages = [[0, 0, 2, 2], [-1, 0, 1, 2], [1, 1, 0, 0]]
  let pca = @src.snpgds_pca(dosages, 2)
  assert_eq(pca.eigenvalues.length(), 2)
  for v in pca.eigenvalues {
    assert_true(v == v && v.abs() <= 1.0e300)
  }
}

///|
test "snpgds_ld_correlation" {
  let dosages = snp_test_dosages()
  // mirror: corr(SNP0, SNP1) = sqrt(0.5) = 0.70710678...
  assert_true(
    snp_close(@src.snpgds_ld(dosages, 0, 1), 0.7071067811865476, 1.0e-12),
  )
  // corr(SNP0, SNP2) = -1/sqrt(3)
  assert_true(
    snp_close(@src.snpgds_ld(dosages, 0, 2), -0.5773502691896258, 1.0e-12),
  )
  // corr(SNP1, SNP2) = 0
  assert_true(snp_close(@src.snpgds_ld(dosages, 1, 2), 0.0, 1.0e-12))
  // out-of-range raises
  let failed = try {
    let _ = @src.snpgds_ld(dosages, 0, 9)
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "snpgds_ld_pruning" {
  // SNP0..SNP3 with corr pairs: (0,1)=0.707, (0,2)=-0.577, (0,3)=0.866,
  // (1,2)=0, (1,3)=0.866, (2,3)=0
  let dosages = snp_test_dosages()
  // threshold 0.9 keeps everything (all |r| below)
  let result_loose = @src.snpgds_ld_pruning(dosages, 0.9, 0.05, 4)
  assert_eq(result_loose.kept_indices, [0, 1, 2, 3])
  assert_eq(result_loose.removed_indices, [])
  // threshold 0.5: SNP0 kept; SNP1 removed (r=0.707); SNP2 removed
  // (|r|=0.577 with SNP0 exceeds the threshold); SNP3 removed (r=0.866
  // with SNP0)
  let result = @src.snpgds_ld_pruning(dosages, 0.5, 0.05, 4)
  assert_eq(result.kept_indices, [0])
  assert_eq(result.removed_indices, [1, 2, 3])
  // high MAF threshold removes all SNPs
  let result_strict = @src.snpgds_ld_pruning(dosages, 0.5, 0.6, 4)
  assert_eq(result_strict.kept_indices, [])
  assert_eq(result_strict.removed_indices, [0, 1, 2, 3])
}

///|
test "snpgds_ld_pruning_sliding_window" {
  // with a slide window of 2 the candidate at index 3 only compares
  // against indices 1-2 (SNP0 is out of the window), so SNP3 survives
  let dosages = snp_test_dosages()
  let result = @src.snpgds_ld_pruning(dosages, 0.5, 0.05, 2)
  // window 2: SNP1 (r with SNP0 = 0.707) removed; SNP2 (|r| with SNP0 =
  // 0.577, inside the window) removed; SNP3 compares against SNP1
  // (removed) and SNP2 (kept, r=0) -> kept
  assert_eq(result.kept_indices, [0, 3])
  assert_eq(result.removed_indices, [1, 2])
}
