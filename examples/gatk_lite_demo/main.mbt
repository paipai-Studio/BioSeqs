///|
/// GATK-lite 变异检出演示: pileup -> genotype likelihood -> active region -> VCF.
fn main {
  println(
    "=== GATK-lite 变异检出 (genotype likelihood + active region) ===",
  )
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGT"
  println("参考 \{reference.length()}bp (非周期)")
  // 10 条 read 覆盖 [10, 30)
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let sliced = reference[10:30].to_owned()
    // 位置 15(读内偏移 5): G -> C 杂合
    let mutated = if i < 5 {
      sliced[0:5].to_owned() + "C" + sliced[6:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "r" + i.to_string(), seq: mutated, pos: 10, })
  }
  // 1. pileup
  println("\n1. 位点 pileup (位置 16, 1-based)")
  let site = @src.gatk_pileup_site(reference, reads, 15)
  println("  depth=\{site.depth} ref=\{site.ref_base}")
  for pair in site.base_counts {
    println("  \{pair.0}: \{pair.1}")
  }
  // 2. genotype likelihood
  println("\n2. 基因型后验")
  let gts = @src.gatk_genotype_likelihoods(site, 'C')
  for g in gts {
    println("  \{g.gt}: posterior=\{g.posterior} GQ=\{g.log10_error}")
  }
  // 3. 变异检出
  println("\n3. 全流程检出")
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  println("  " + @src.gatk_variants_summary(variants))
  for v in variants {
    println(
      "  \{v.chrom}:\{v.pos} \{v.ref}>\{v.alt} GT=\{v.gt} DP=\{v.dp} AD=\{v.ad_ref},\{v.ad_alt} GQ=\{v.gq} AF=\{v.af}",
    )
  }
  // 4. indel 场景
  println("\n4. indel 检出 (3bp 缺失)")
  let del_reads : Array[@src.GatkRead] = []
  for i = 0; i < 8; i = i + 1 {
    let sliced = reference[0:20].to_owned()
    let seq = if i < 6 {
      sliced[0:5].to_owned() + sliced[8:].to_owned()
    } else {
      sliced
    }
    del_reads.push({ name: "d" + i.to_string(), seq, pos: 0, })
  }
  let del_variants = @src.gatk_call_variants("chr1", reference, del_reads)
  println("  " + @src.gatk_variants_summary(del_variants))
  for v in del_variants {
    println(
      "  \{v.chrom}:\{v.pos} \{v.ref}>\{v.alt} GT=\{v.gt} AD=\{v.ad_ref},\{v.ad_alt}",
    )
  }
  // 5. VCF 输出与 parse_vcf 闭环
  println("\n5. VCF 4.2 输出 (与 parse_vcf 闭环)")
  let vcf_text = @src.gatk_write_vcf("chr1", reference, variants)
  let lines = vcf_text.split("\n").to_array()
  for l in lines {
    let s = l.to_owned()
    if s.length() > 0 && (s.has_prefix("#") || s.has_prefix("chr1")) {
      println("  " + s)
    }
  }
  let vcf = @src.parse_vcf(vcf_text)
  println(
    "  parse_vcf 解析: \{vcf.len()} 条记录, SNP=\{vcf.records[0].is_snp()}",
  )
  println("\n完成: pileup -> 似然 -> 检出 -> VCF 全链路本地完成.")
}
