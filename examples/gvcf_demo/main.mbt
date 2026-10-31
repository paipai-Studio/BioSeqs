///|
/// gVCF 参考块 / 多样本合并 / 联合基因分型演示.
///
/// 流程: 两个样本各自的变异 -> 单样本 gVCF(变异记录 + 参考块) -> 合并 ->
/// 联合基因分型 -> 多样本 VCF + JSON.

///|
fn gvcf_demo_ref() -> String {
  "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGT"
}

///|
fn gvcf_demo_read(
  reference : String,
  pos : Int32,
  offset : Int32?,
  alt : Char,
  name : String,
) -> @src.GatkRead {
  let end = if pos + 20 < reference.length() {
    pos + 20
  } else {
    reference.length()
  }
  let mut seq = reference[pos:end].to_owned()
  match offset {
    Some(off) =>
      if off >= 0 && off < seq.length() {
        seq = seq[0:off].to_owned() + alt.to_string() + seq[off + 1:].to_owned()
      }
    None => ()
  }
  { name, seq, pos, }
}

///|
fn gvcf_demo_sample_variants(
  reference : String,
  start : Int32,
  offset : Int32,
  alt : Char,
  sample : String,
) -> Array[@src.GatkVariant] {
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(offset) } else { None }
    reads.push(
      gvcf_demo_read(reference, start, off, alt, sample + "_r" + i.to_string()),
    )
  }
  @src.gatk_call_variants("chr1", reference, reads)
}

///|
fn gvcf_demo_show(callset : @src.GvcfCallset) -> Unit {
  for record in callset.records {
    if record.is_block {
      println(
        "  [block]   \{record.chrom}:\{record.pos}-\{record.block_end} ref=\{record.ref_allele} <NON_REF>",
      )
    } else {
      let alts = StringBuilder()
      for i, alt in record.alts {
        if i > 0 {
          alts.write_char(',')
        }
        alts.write_string(alt)
      }
      let mut sample_note = ""
      for s in callset.samples {
        match record.samples.get(s) {
          Some(g) => sample_note = sample_note + " " + s + "=" + g.gt
          None => ()
        }
      }
      println(
        "  [variant] \{record.chrom}:\{record.pos} ref=\{record.ref_allele} alt=\{alts.to_string()}\{sample_note}",
      )
    }
  }
}

///|
fn main {
  println("=== gVCF 参考块 / 多样本合并 / 联合基因分型 ===")
  let reference = gvcf_demo_ref()
  println("参考长度: \{reference.length()}")

  let variants_s1 = gvcf_demo_sample_variants(reference, 10, 5, 'C', "S1")
  let variants_s2 = gvcf_demo_sample_variants(reference, 20, 4, 'C', "S2")
  println("\n样本 S1 检出变异: \{variants_s1.length()}")
  println("样本 S2 检出变异: \{variants_s2.length()}")

  let cs1 = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    variants_s1,
    sample="S1",
  ) catch {
    @src.GvcfError(message) => {
      println("S1 gVCF 构建失败: " + message)
      return
    }
  }
  let cs2 = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    variants_s2,
    sample="S2",
  ) catch {
    @src.GvcfError(message) => {
      println("S2 gVCF 构建失败: " + message)
      return
    }
  }
  println("\n1. 单样本 gVCF (S1)")
  gvcf_demo_show(cs1)
  println("  摘要: \{@src.gvcf_summary(cs1)}")
  println("\n   单样本 gVCF (S2)")
  gvcf_demo_show(cs2)
  println("  摘要: \{@src.gvcf_summary(cs2)}")

  let merged = @src.gvcf_merge([cs1, cs2]) catch {
    @src.GvcfError(message) => {
      println("合并失败: " + message)
      return
    }
  }
  println(
    "\n2. 多样本合并(未覆盖到变异的样本按参考块判为 0/0)",
  )
  for record in merged.records {
    let mut note = ""
    for s in merged.samples {
      match record.samples.get(s) {
        Some(g) => note = note + " " + s + "=" + g.gt
        None => ()
      }
    }
    println("  \{record.chrom}:\{record.pos} ref=\{record.ref_allele} \{note}")
  }

  let joint = @src.gvcf_joint_genotype(merged)
  println("\n3. 联合基因分型(按各样本 PL 重算 GT)")
  for record in joint.records {
    let mut note = ""
    for s in joint.samples {
      match record.samples.get(s) {
        Some(g) =>
          note = note + " " + s + "=" + g.gt + "(GQ=" + g.gq.to_string() + ")"
        None => ()
      }
    }
    println("  \{record.chrom}:\{record.pos} ref=\{record.ref_allele} \{note}")
  }

  println("\n4. 多样本 VCF 输出")
  println(@src.gvcf_write_vcf(joint))

  println("5. JSON 序列化 / roundtrip")
  let json = @src.gvcf_callset_to_json(joint)
  println("  JSON 长度 = \{json.length()} 字节")
  match @src.gvcf_callset_from_json(json) {
    Some(restored) =>
      println(
        "  还原: samples=\{restored.samples.length()}, records=\{restored.records.length()}",
      )
    None => println("  roundtrip 失败")
  }

  println(
    "\n演示完成: gVCF 参考块 + 多样本合并 + 联合基因分型闭环",
  )
}
