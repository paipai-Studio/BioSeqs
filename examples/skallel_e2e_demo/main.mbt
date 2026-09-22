///|
/// scikit-allel end-to-end demo: a small VCF → genotype/haplotype arrays →
/// linkage disequilibrium (LD), PCA and iHS, through the `allel.io` bridge.

///|
fn fmt_i(row : Array[Int]) -> String {
  let sb = StringBuilder::new()
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
  let sb = StringBuilder::new()
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
  let vcf_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/0\nchr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/0\nchr1\t300\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t1/1\nchr1\t400\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\n"
  let vcf = @src.parse_vcf(vcf_text)
  let samples = vcf.get_sample_names()
  let g = @src.skallel_vcf_genotypes(vcf, samples)

  println("samples = \{samples}")

  let gn = @src.skallel_gt_to_n_alt(g, fill=-1)
  println("\n-- dosage matrix gn (missing = -1) --")
  for row in gn {
    println(fmt_i(row))
  }

  let r = @src.skallel_rogers_huff_r(gn)
  println("\n-- linkage disequilibrium (Rogers–Huff r, condensed) --")
  println(fmt_d(r))

  let h = @src.skallel_gt_to_haplotypes(g)
  println("\n-- haplotype array h --")
  for row in h {
    println(fmt_i(row))
  }

  let ihs = @src.skallel_ihs(h, [100, 200, 300, 400], 0.0, false)
  println("\n-- iHS (unstandardized) --")
  println(fmt_d(ihs))

  let gn0 = @src.skallel_gt_to_n_alt(g)
  let pca = @src.skallel_pca(gn0, 3)
  println("\n-- PCA explained-variance ratio --")
  println(fmt_d(pca.explained_variance_ratio))
  println("-- PCA coordinates (samples x components) --")
  for row in pca.coords {
    println(fmt_d(row))
  }
}