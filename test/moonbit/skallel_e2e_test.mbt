///|
/// End-to-end chain test: a small VCF → genotype/haplotype arrays → LD, PCA and
/// iHS. Verifies that the `allel.io` bridge (`skallel_vcf_genotypes` /
/// `skallel_gt_to_*`) feeds the `allel.stats` functions correctly.

///|
test "skallel_e2e_vcf_ld_pca_ihs" {
  let vcf_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/0\nchr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/0\nchr1\t300\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t1/1\nchr1\t400\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\n"
  let vcf = @src.parse_vcf(vcf_text)
  let samples = vcf.get_sample_names()
  let g = @src.skallel_vcf_genotypes(vcf, samples)

  // dosage matrix (missing -> -1)
  let gn = @src.skallel_gt_to_n_alt(g, fill=-1)
  assert_eq((gn[0][0], gn[0][1], gn[0][2]), (0, 2, 0))
  assert_eq((gn[3][0], gn[3][1], gn[3][2]), (0, 1, -1))

  // LD: pairwise Rogers–Huff r (matches the scikit-allel docstring example)
  let r = @src.skallel_rogers_huff_r(gn)
  let exp = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
  assert_eq(r.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    assert_true(assert_close(r[i], exp[i], 1.0e-5))
  }

  // haplotypes -> iHS
  let h = @src.skallel_gt_to_haplotypes(g)
  assert_eq((h[0][0], h[0][1], h[0][2], h[0][3], h[0][4], h[0][5]), (0, 0, 1, 1, 0, 0))
  let ihs = @src.skallel_ihs(h, [100, 200, 300, 400], 0.0, false)
  assert_eq(ihs.length(), 4)

  // PCA on missing-free dosage (fill = 0)
  let gn0 = @src.skallel_gt_to_n_alt(g)
  let pca = @src.skallel_pca(gn0, 3)
  assert_eq(pca.coords.length(), 3)
  assert_eq(pca.coords[0].length(), 3)
  let mut s = 0.0
  for x in pca.explained_variance_ratio {
    s = s + x
  }
  assert_true(assert_close(s, 1.0, 1.0e-9))
}