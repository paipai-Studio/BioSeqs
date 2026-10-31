///|
/// mutect_lite_demo — Mutect2 风格肿瘤-正常体细胞变异检出演示:
/// Part A: somatic SNV / germline / pop_af 的 LOD 检出 -> TMB ->
///         somatic VCF -> 突变签名 96 通道衔接;
/// Part B: indel 检出与简化推断的移位伪影说明.

///|
fn ml_demo_pick_base(rng : @src.Rng) -> String {
  let i = rng.randint(4)
  match i {
    0 => "A"
    1 => "C"
    2 => "G"
    _ => "T"
  }
}

///|
fn ml_demo_other(base : String) -> String {
  match base {
    "A" => "T"
    "C" => "G"
    "G" => "A"
    _ => "C"
  }
}

///|
fn ml_demo_char_at(s : String, pos : Int32) -> String {
  s.unsafe_get(pos).unsafe_to_char().to_string()
}

///|
fn ml_demo_replace(s : String, pos : Int32, b : String) -> String {
  let head = if pos > 0 { s[:pos].to_owned() } else { "" }
  let tail = if pos + 1 < s.length() { s[pos + 1:].to_owned() } else { "" }
  head + b + tail
}

///|
fn ml_demo_insert(s : String, anchor : Int32, ins : String) -> String {
  s[:anchor + 1].to_owned() + ins + s[anchor + 1:].to_owned()
}

///|
fn main {
  let rng = @src.Rng::new(20261031UL)
  let sb = StringBuilder()
  for _i = 0; _i < 120; _i = _i + 1 {
    sb.write_string(ml_demo_pick_base(rng))
  }
  let reference = sb.to_string()
  println("1. 合成参考: 120bp (确定性 RNG seed=20261031)")
  // ---- Part A: SNV / germline / pop_af ----
  let snv1_pos : Int32 = 40
  let snv1_alt = ml_demo_other(ml_demo_char_at(reference, snv1_pos))
  let snv2_pos : Int32 = 90
  let snv2_alt = ml_demo_other(ml_demo_char_at(reference, snv2_pos))
  let germ_pos : Int32 = 70
  let germ_alt = ml_demo_other(ml_demo_char_at(reference, germ_pos))
  let pop_pos : Int32 = 30
  let pop_alt = ml_demo_other(ml_demo_char_at(reference, pop_pos))
  let tumor_reads : Array[@src.GatkRead] = Array(capacity=20)
  for i = 0; i < 20; i = i + 1 {
    let mut seq = reference
    if i < 8 {
      seq = ml_demo_replace(seq, snv1_pos, snv1_alt)
    }
    if i < 6 {
      seq = ml_demo_replace(seq, snv2_pos, snv2_alt)
    }
    if i < 10 {
      seq = ml_demo_replace(seq, germ_pos, germ_alt)
    }
    if i < 6 {
      seq = ml_demo_replace(seq, pop_pos, pop_alt)
    }
    tumor_reads.push(@src.GatkRead::{ name: "t\{i}", seq, pos: 0, })
  }
  let normal_reads : Array[@src.GatkRead] = Array(capacity=20)
  for i = 0; i < 20; i = i + 1 {
    let seq = if i < 10 {
      ml_demo_replace(reference, germ_pos, germ_alt)
    } else {
      reference
    }
    normal_reads.push(@src.GatkRead::{ name: "n\{i}", seq, pos: 0, })
  }
  let tumor = @src.MutectSample::{ name: "TUMOR", reads: tumor_reads, }
  let normal = @src.MutectSample::{ name: "NORMAL", reads: normal_reads, }
  println("2. Part A 注入: somatic SNV x2 + germline x1 + pop_af x1")
  let pop_af : Map[String, Double] = Map([], capacity=0)
  let pop_key = "chr1:" +
    (pop_pos + 1).to_string() +
    ":" +
    ml_demo_char_at(reference, pop_pos) +
    ">" +
    pop_alt
  pop_af.set(pop_key, 0.05)
  let calls = @src.mutect_call("chr1", reference, tumor, normal, pop_af~)
  println("3. Part A 检出 (pos ref>alt | tumor_af | lod_t | lod_n | filter):")
  for v in calls {
    println(
      "  \{v.pos} \{v.ref_allele}>\{v.alt} | \{v.tumor_af} | \{v.lod_t} | \{v.lod_n} | \{v.filter}",
    )
  }
  let summary = @src.mutect_summary(calls, callable_mb=0.12)
  println(
    "4. \{summary.n_pass} (SNV=\{summary.n_snv}, indel=\{summary.n_indel}), TMB=\{summary.tmb_per_mb} /Mb",
  )
  println(@src.mutect_summary_line(summary))
  let vcf = @src.mutect_write_vcf(calls, "TUMOR", "NORMAL")
  println("5. somatic VCF (数据行, 共 \{vcf.length()} 字节):")
  let mut shown = 0
  for line in vcf.split("\n") {
    if line != "" && line[0:1].to_owned() != "#" && shown < 4 {
      println("  " + line.to_owned())
      shown = shown + 1
    }
  }
  let muts = @src.mutect_to_somatic_mutations(
    reference,
    calls,
    sample_name="TUMOR",
  )
  let matrix = @src.build_mutation_matrix(muts, ["TUMOR"])
  println(
    "6. 突变签名衔接: PASS SNV \{muts.length()} 条进入 96 通道矩阵",
  )
  println(matrix.spectrum_summary())
  // ---- Part B: indel 与移位伪影 ----
  println(
    "7. Part B: indel 检出(独立样本, 锚点 105 后插入 GGA, 6/20 条)",
  )
  let ins_anchor : Int32 = 105
  let ins_tumor_reads : Array[@src.GatkRead] = Array(capacity=20)
  for i = 0; i < 20; i = i + 1 {
    let seq = if i < 6 {
      ml_demo_insert(reference, ins_anchor, "GGA")
    } else {
      reference
    }
    ins_tumor_reads.push(@src.GatkRead::{ name: "t\{i}", seq, pos: 0, })
  }
  let ins_tumor = @src.MutectSample::{ name: "TUMOR", reads: ins_tumor_reads, }
  let ins_calls = @src.mutect_call("chr1", reference, ins_tumor, normal)
  println("  检出 \{ins_calls.length()} 条候选, 其中插入类:")
  for v in ins_calls {
    if v.ref_allele.length() == 1 && v.alt.length() > 1 {
      println(
        "  \{v.pos} \{v.ref_allele}>\{v.alt} | \{v.tumor_af} | \{v.filter}",
      )
    }
  }
  println(
    "  注: 其余候选为插入读段下游移位产生的伪影 SNV(简化逐位点推断的",
  )
  println(
    "  已知局限, Mutect2 以局部组装解决; 本模块以 SNV 优先清洗移除部分伪插入)",
  )
  println(
    "演示完成: LOD 检出 -> 过滤 -> TMB -> VCF -> 突变签名闭环",
  )
}
