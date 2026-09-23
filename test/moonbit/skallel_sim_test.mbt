///|
/// Simulation sanity test: a hard-sweep VCF should produce (1) a positive,
/// localized iHS signal near the sweep, and (2) a non-degenerate PCA.

///|
test "skallel_sim_sweep" {
  let vcf_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\nchr1\t1\t.\tA\tT\t.\t.\t.\tGT\t1/0\t1/1\t1/1\t1/0\t0/1\t1/0\nchr1\t2\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/1\t1/1\t1/0\t1/0\nchr1\t3\t.\tA\tT\t.\t.\t.\tGT\t1/0\t1/0\t0/0\t1/0\t0/0\t0/0\nchr1\t4\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t0/1\t1/0\t1/0\t0/0\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/0\t1/1\nchr1\t6\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t1/0\t0/1\t0/1\nchr1\t7\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t1/1\t1/1\t0/0\nchr1\t8\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/0\t0/0\t0/0\t1/1\t0/0\nchr1\t9\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/0\t1/1\t0/0\t0/1\t1/1\nchr1\t10\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/1\t0/0\t1/1\t1/0\nchr1\t11\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/1\t1/1\t1/1\nchr1\t12\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t0/1\t1/1\t1/0\t0/0\nchr1\t13\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/0\t0/0\t1/1\nchr1\t14\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1\t0/0\t1/0\nchr1\t15\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t1/0\t1/0\t0/0\nchr1\t16\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1\t1/1\t0/1\nchr1\t17\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/0\t0/1\t0/1\t1/1\t0/0\nchr1\t18\t.\tA\tT\t.\t.\t.\tGT\t1/0\t0/1\t0/1\t0/1\t1/0\t0/0\nchr1\t19\t.\tA\tT\t.\t.\t.\tGT\t1/0\t0/1\t0/0\t1/1\t1/0\t1/0\nchr1\t20\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/0\t0/1\t1/0\t1/0\t0/0\n"
  let vcf = @src.parse_vcf(vcf_text)
  let samples = vcf.get_sample_names()
  let g = @src.skallel_vcf_genotypes(vcf, samples)
  let h = @src.skallel_gt_to_haplotypes(g)
  let gn = @src.skallel_gt_to_n_alt(g)

  let pos : Array[Int] = Array(capacity=0)
  for v = 0; v < 20; v = v + 1 {
    pos.push(v + 1)
  }

  // iHS: a clear positive peak within the sweep region (variants 3..12)
  let ihs = @src.skallel_ihs(h, pos, 0.0, false)
  assert_eq(ihs.length(), 20)
  let mut max_v = 0.0
  let mut max_i = 0
  for i = 0; i < 20; i = i + 1 {
    if ihs[i] > max_v {
      max_v = ihs[i]
      max_i = i
    }
  }
  assert_true(max_v > 0.5)
  assert_true(max_i >= 2 && max_i <= 12)

  // nSL: at least one positive peak (sweep region is non-degenerate)
  let nsl = @src.skallel_nsl(h)
  let mut nsl_max = 0.0
  for x in nsl {
    if x > nsl_max {
      nsl_max = x
    }
  }
  assert_true(nsl_max > 0.3)

  // PCA: non-degenerate (first component explains < 90% of variance)
  let pca = @src.skallel_pca(gn, 5)
  assert_true(pca.explained_variance_ratio[0] < 0.9)
}