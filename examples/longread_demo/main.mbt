///|
/// 长读长深化演示: minimap 预设比对 -> medaka-lite 变异检出.
fn main {
  println("=== 长读长管线: minimap 预设 + medaka-lite 变异检出 ===")
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGTCCATTGCAGGT"
  println("参考 \{reference.length()}bp")
  let refs : Array[@src.MinimapRefSeq] = [
    @src.MinimapRefSeq::new("ctg1", reference),
  ]
  // 1. 三种预设
  println("\n1. 预设参数 (对齐 minimap2 -x)")
  for preset in [@src.MapOnt, @src.MapPb, @src.MapSr] {
    let (k, w, band, skip) = preset.params()
    println(
      "  \{preset.to_string()}: k=\{k} w=\{w} band=\{band} chain_skip=\{skip}",
    )
  }
  // 2. map-ont 比对
  println("\n2. map-ont 长读比对")
  let index = @src.minimap_index_build_preset(refs, @src.MapOnt)
  let region = reference[10:50].to_owned()
  // 模拟 ONT 噪声: 2 处错配
  let noisy = region[0:5].to_owned() +
    "T" +
    region[6:12].to_owned() +
    "G" +
    region[13:].to_owned()
  let results = @src.minimap_align_preset(
    index,
    noisy,
    "ONT_READ_001",
    @src.MapOnt,
  )
  for r in results {
    println(
      "  \{r.query_name} -> \{r.ref_name}:\{r.ref_start}-\{r.ref_end} strand=\{r.strand} mapq=\{r.mapq} cigar=\{r.cigar}",
    )
  }
  // 3. medaka-lite 变异检出
  println("\n3. medaka-lite 变异检出 (ONT 错误模型)")
  // 单一杂合 SNP 场景(简化共识模型不处理缺失下游的连锁错位)
  let reads : Array[@src.MedakaRead] = []
  for i = 0; i < 5; i = i + 1 {
    let sliced = reference[0:40].to_owned()
    let seq = if i < 3 {
      // 杂合 SNP(位置 15: G -> A)
      sliced[0:15].to_owned() + "A" + sliced[16:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "lr" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.medaka_call_variants("ctg1", reference, reads)
  println("  " + @src.medaka_variants_summary(variants))
  for v in variants {
    println(
      "  \{v.chrom}:\{v.pos} \{v.ref}>\{v.alt} GT=\{v.gt} DP=\{v.dp} AD=\{v.ad_ref},\{v.ad_alt} GQ=\{v.gq}",
    )
  }
  // 4. HiFi 对照
  println("\n4. HiFi 模型对照 (相同证据更高 GQ)")
  let v_hifi = @src.medaka_call_variants(
    "ctg1",
    reference,
    reads,
    model=@src.Hifi,
  )
  for v in v_hifi {
    println(
      "  \{v.chrom}:\{v.pos} \{v.ref}>\{v.alt} GQ=\{v.gq} (ONT 时更低)",
    )
  }
  // 5. VCF 输出
  println("\n5. VCF 输出 (复用 gatk_write_vcf)")
  let vcf_text = @src.gatk_write_vcf(
    "ctg1",
    reference,
    variants,
    sample="LONGREAD",
  )
  let lines = vcf_text.split("\n").to_array()
  for l in lines {
    let s = l.to_owned()
    if s.has_prefix("ctg1") {
      println("  " + s)
    }
  }
  println(
    "\n完成: 长读预设比对 -> 高错误率共识变异检出 -> VCF 全链路本地完成.",
  )
}
