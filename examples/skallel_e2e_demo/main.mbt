///|
/// scikit-allel end-to-end demo: a simulated hard-sweep VCF (20 variants,
/// 6 samples, 12 haplotypes) → genotype/haplotype arrays → iHS / nSL / PCA,
/// through the `allel.io` bridge. The sweep carries a conserved ALT haplotype
/// (sample S1) over variants 5–15, giving a positive iHS/nSL signal.

///|
fn fmt_i(row : Array[Int]) -> String {
  let sb = StringBuilder()
  let mut first = true
  for x in row {
    if not(first) {
      sb.write_char(' ')
    }
    sb.write_string(x.to_string())
    first = false
  }
  sb.to_string()
}

///|
fn fmt_d(row : Array[Double]) -> String {
  let sb = StringBuilder()
  let mut first = true
  for x in row {
    if not(first) {
      sb.write_char(' ')
    }
    sb.write_string(x.to_string())
    first = false
  }
  sb.to_string()
}

///|
fn main {
  let vcf_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\nchr1\t1\t.\tA\tT\t.\t.\t.\tGT\t1/0\t1/1\t1/1\t1/0\t0/1\t1/0\nchr1\t2\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/1\t1/1\t1/0\t1/0\nchr1\t3\t.\tA\tT\t.\t.\t.\tGT\t1/0\t1/0\t0/0\t1/0\t0/0\t0/0\nchr1\t4\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t0/1\t1/0\t1/0\t0/0\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/0\t1/1\nchr1\t6\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t1/0\t0/1\t0/1\nchr1\t7\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t1/1\t1/1\t0/0\nchr1\t8\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/0\t0/0\t0/0\t1/1\t0/0\nchr1\t9\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/0\t1/1\t0/0\t0/1\t1/1\nchr1\t10\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/1\t0/0\t1/1\t1/0\nchr1\t11\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/1\t1/1\t1/1\nchr1\t12\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t0/1\t1/1\t1/0\t0/0\nchr1\t13\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/0\t0/0\t1/1\nchr1\t14\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1\t0/0\t1/0\nchr1\t15\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t1/0\t1/0\t0/0\nchr1\t16\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1\t1/1\t0/1\nchr1\t17\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/0\t0/1\t0/1\t1/1\t0/0\nchr1\t18\t.\tA\tT\t.\t.\t.\tGT\t1/0\t0/1\t0/1\t0/1\t1/0\t0/0\nchr1\t19\t.\tA\tT\t.\t.\t.\tGT\t1/0\t0/1\t0/0\t1/1\t1/0\t1/0\nchr1\t20\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/0\t0/1\t1/0\t1/0\t0/0\n"
  let vcf = @src.parse_vcf(vcf_text)
  let samples = vcf.get_sample_names()
  let g = @src.skallel_vcf_genotypes(vcf, samples)
  let h = @src.skallel_gt_to_haplotypes(g)
  let gn = @src.skallel_gt_to_n_alt(g, fill=-1)

  let pos2 : Array[Int] = Array(capacity=0)
  for v = 0; v < 20; v = v + 1 {
    pos2.push(v + 1)
  }

  println("samples = \{samples}")
  println("\n-- iHS (unstandardized; sweep peak expected near variants 3–11) --")
  println(fmt_d(@src.skallel_ihs(h, pos2, 0.0, false)))
  println("\n-- nSL (unstandardized) --")
  println(fmt_d(@src.skallel_nsl(h)))
  println("\n-- LD (Rogers–Huff r, condensed) --")
  println(fmt_d(@src.skallel_rogers_huff_r(gn)))

  let pca = @src.skallel_pca(gn, 5)
  println("\n-- PCA explained-variance ratio --")
  println(fmt_d(pca.explained_variance_ratio))
  println("-- PCA coordinates (samples x components; S1 = sweep sample) --")
  for row in pca.coords {
    println(fmt_d(row))
  }
  // --- SFS / Hardy–Weinberg / Patterson Fst ---
  let ac = @src.skallel_gt_count_alleles(g, 1)
  let dac : Array[Int] = Array(capacity=0)
  for aci in ac {
    dac.push(aci[1])
  }
  println("\n-- SFS unfolded (ALT allele count, n=12 chromosomes) --")
  println(fmt_i(@src.skallel_sfs(dac, Some(12))))
  println("-- SFS folded (minor allele count) --")
  println(fmt_i(@src.skallel_sfs_folded(ac, Some(12))))

  println("\n-- observed heterozygosity Ho (per variant) --")
  println(fmt_d(@src.skallel_heterozygosity_observed(g, -1.0)))
  println("-- inbreeding coefficient F (per variant; -999 = He=0) --")
  println(fmt_d(@src.skallel_inbreeding_coefficient(g, -999.0)))

  println("\n-- Patterson Fst (S1-S3 vs S4-S6) --")
  let ac1 = @src.skallel_gt_count_alleles_subpop(g, [0, 1, 2], 1)
  let ac2 = @src.skallel_gt_count_alleles_subpop(g, [3, 4, 5], 1)
  let (num, den) = @src.skallel_patterson_fst(ac1, ac2)
  let mut sn = 0.0
  let mut sd = 0.0
  for i = 0; i < num.length(); i = i + 1 {
    sn = sn + num[i]
    sd = sd + den[i]
  }
  println("averaged Fst = " + (sn / sd).to_string())

}