///|
/// Tests for the skallel VCF→array bridge. Values are hand-computed.

///|
test "skallel_vcf_bridge" {
  let content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\nchr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\nchr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t./.\t1/1\n"
  let vcf = @src.parse_vcf(content)
  let samples = vcf.get_sample_names()
  assert_eq(samples.length(), 2)
  assert_eq(samples[0], "S1")

  // GT parse
  assert_eq(@src.skallel_parse_gt("0/1").length(), 2)
  assert_eq(@src.skallel_parse_gt("1|0")[1], 0)
  assert_eq(@src.skallel_parse_gt("./.")[0], -1)

  let g = @src.skallel_vcf_genotypes(vcf, samples)
  assert_eq(g.length(), 3)
  assert_eq((g[0][0][0], g[0][0][1]), (0, 0))
  assert_eq((g[0][1][0], g[0][1][1]), (0, 1))
  assert_eq((g[2][0][0], g[2][0][1]), (-1, -1))
  assert_eq((g[2][1][0], g[2][1][1]), (1, 1))

  // to_n_alt (missing -> 0 by default)
  let gn = @src.skallel_gt_to_n_alt(g)
  assert_eq((gn[0][0], gn[0][1]), (0, 1))
  assert_eq((gn[1][0], gn[1][1]), (1, 2))
  assert_eq((gn[2][0], gn[2][1]), (0, 2))
  // missing -> -1
  assert_eq(@src.skallel_gt_to_n_alt(g, fill=-1)[2][0], -1)

  // haplotypes
  let h = @src.skallel_gt_to_haplotypes(g)
  assert_eq((h[0][0], h[0][1], h[0][2], h[0][3]), (0, 0, 0, 1))
  assert_eq((h[1][0], h[1][1], h[1][2], h[1][3]), (0, 1, 1, 1))
  assert_eq((h[2][0], h[2][1], h[2][2], h[2][3]), (-1, -1, 1, 1))

  // count_alleles
  let ac = @src.skallel_gt_count_alleles(g, 1)
  assert_eq((ac[0][0], ac[0][1]), (3, 1))
  assert_eq((ac[1][0], ac[1][1]), (1, 3))
  assert_eq((ac[2][0], ac[2][1]), (0, 2))

  // count_alleles_subpop (sample 0 only)
  let ac0 = @src.skallel_gt_count_alleles_subpop(g, [0], 1)
  assert_eq((ac0[0][0], ac0[0][1]), (2, 0))
  assert_eq((ac0[1][0], ac0[1][1]), (1, 1))
  assert_eq((ac0[2][0], ac0[2][1]), (0, 0))
}