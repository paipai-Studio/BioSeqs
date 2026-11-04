// ===========================================================================
// mutect_lite_test.mbt — Mutect2 风格体细胞变异检出测试.
// ===========================================================================

///|
/// 80bp 确定性参考("ACGT" 循环 20 次): pos 40 为 'A', pos 55 为 'T'.
fn ml_reference() -> String {
  let sb = StringBuilder()
  for i = 0; i < 80; i = i + 1 {
    let b = match i % 4 {
      0 => "A"
      1 => "C"
      2 => "G"
      _ => "T"
    }
    sb.write_string(b)
  }
  sb.to_string()
}

///|
/// 把 s 的 pos(0-based) 处碱基替换为 b.
fn ml_replace(s : String, pos : Int32, b : Char) -> String {
  let head = if pos > 0 { s[:pos].to_owned() } else { "" }
  let tail = if pos + 1 < s.length() { s[pos + 1:].to_owned() } else { "" }
  head + b.to_string() + tail
}

///|
/// n 条全参考读段, 其中前 n_alt 条在 alt_pos(0-based) 替换为 alt_base.
fn ml_make_reads(
  reference : String,
  n : Int32,
  alt_pos : Int32,
  alt_base : Char,
  n_alt : Int32,
) -> Array[@src.GatkRead] {
  let reads : Array[@src.GatkRead] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let seq = if i < n_alt {
      ml_replace(reference, alt_pos, alt_base)
    } else {
      reference
    }
    reads.push(@src.GatkRead::{ name: "r\{i}", seq, pos: 0, })
  }
  reads
}

///|
/// n 条全参考读段, 其中前 n_ins 条在锚点 anchor(0-based) 之后插入 ins_seq.
fn ml_make_insert_reads(
  reference : String,
  n : Int32,
  anchor : Int32,
  ins_seq : String,
  n_ins : Int32,
) -> Array[@src.GatkRead] {
  let reads : Array[@src.GatkRead] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let seq = if i < n_ins {
      reference[:anchor + 1].to_owned() +
      ins_seq +
      reference[anchor + 1:].to_owned()
    } else {
      reference
    }
    reads.push(@src.GatkRead::{ name: "r\{i}", seq, pos: 0, })
  }
  reads
}

///|
/// n 条全参考读段, 其中前 n_del 条删除锚点 anchor 之后 dl 个碱基.
fn ml_make_delete_reads(
  reference : String,
  n : Int32,
  anchor : Int32,
  dl : Int32,
  n_del : Int32,
) -> Array[@src.GatkRead] {
  let reads : Array[@src.GatkRead] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let seq = if i < n_del {
      reference[:anchor + 1].to_owned() + reference[anchor + 1 + dl:].to_owned()
    } else {
      reference
    }
    reads.push(@src.GatkRead::{ name: "r\{i}", seq, pos: 0, })
  }
  reads
}

///|
test "mutect: somatic SNV 检出 PASS" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  inspect(calls.length(), content="1")
  let v = calls[0]
  inspect(v.pos, content="41")
  inspect(v.ref_allele, content="A")
  inspect(v.alt, content="T")
  inspect(v.filter, content="PASS")
  inspect(v.tumor_alt, content="8")
  inspect(v.tumor_ref, content="12")
  inspect(v.normal_alt, content="0")
  assert_eq(v.tumor_af, 0.4)
  assert_true(v.lod_t > 5.3)
  assert_true(v.lod_n < 0.000001)
}

///|
test "mutect: germline 过滤(normal 含 alt 证据)" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 10),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, 40, 'T', 10),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  inspect(calls.length(), content="1")
  inspect(calls[0].filter, content="germline")
  assert_true(calls[0].lod_n > 2.2)
}

///|
test "mutect: weak evidence 过滤(低 LOD)" {
  let ref_ = ml_reference()
  // 3/20 alt: af=0.15 但 LOD 不足以支持 somatic
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 3),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  inspect(calls.length(), content="1")
  inspect(calls[0].filter, content="weak_evidence")
  assert_true(calls[0].lod_t < 5.3)
}

///|
test "mutect: 人群 AF 过滤(pop_af 表)" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let pop_af : Map[String, Double] = Map([], capacity=0)
  pop_af.set("chr1:41:A>T", 0.1)
  let calls = @src.mutect_call("chr1", ref_, tumor, normal, pop_af~)
  inspect(calls.length(), content="1")
  inspect(calls[0].filter, content="pop_af")
}

///|
test "mutect: 插入检出" {
  let ref_ = ml_reference()
  // 锚点 50('G') 后插入 "GGT", 5/20 条读段
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_insert_reads(ref_, 20, 50, "GGT", 5),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  // 简化 indel 推断把插入报在插入物内部锚点(此处为 T>TGTT),
  // 断言存在 PASS 的插入候选(ref 单碱基, alt = 锚 + 3bp 插入物)
  let mut found = false
  for v in calls {
    if v.ref_allele.length() == 1 && v.alt.length() == 4 {
      assert_eq(v.filter, "PASS")
      assert_eq(v.tumor_alt, 5)
      assert_true(v.lod_t > 5.3)
      found = true
    }
  }
  assert_eq(found, true)
}

///|
test "mutect: 缺失检出" {
  let ref_ = ml_reference()
  // 锚点 55('T') 后删 3 碱基, 5/20 条读段
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_delete_reads(ref_, 20, 55, 3, 5),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  let mut found = false
  for v in calls {
    if v.ref_allele.length() == 4 && v.alt.length() == 1 && v.alt == "T" {
      assert_eq(v.filter, "PASS")
      found = true
    }
  }
  assert_eq(found, true)
}

///|
test "mutect: 汇总与 TMB" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  let summary = @src.mutect_summary(calls, callable_mb=2.0)
  inspect(summary.n_candidates, content="1")
  inspect(summary.n_pass, content="1")
  inspect(summary.n_snv, content="1")
  inspect(summary.n_indel, content="0")
  assert_eq(summary.tmb_per_mb, 0.5)
  let pass = @src.mutect_pass(calls)
  inspect(pass.length(), content="1")
}

///|
test "mutect: somatic VCF 输出格式" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  let vcf = @src.mutect_write_vcf(calls, "TUMOR", "NORMAL")
  assert_eq(vcf.contains("##fileformat=VCFv4.2"), true)
  assert_eq(vcf.contains("##INFO=<ID=TLOD"), true)
  assert_eq(vcf.contains("##FILTER=<ID=germline"), true)
  assert_eq(vcf.contains("FORMAT\tTUMOR\tNORMAL"), true)
  assert_eq(vcf.contains("chr1\t41\t.\tA\tT"), true)
  assert_eq(vcf.contains("TLOD="), true)
  assert_eq(vcf.contains("0/1:12,8:0.4"), true)
  assert_eq(vcf.contains("0/0:20,0:0"), true)
}

///|
test "mutect: 衔接 mutational_patterns 突变签名" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{
    name: "NORMAL",
    reads: ml_make_reads(ref_, 20, -1, ' ', 0),
  }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  let muts = @src.mutect_to_somatic_mutations(ref_, calls, sample_name="TUMOR")
  inspect(muts.length(), content="1")
  inspect(muts[0].chromosome, content="TUMOR")
  inspect(muts[0].position, content="41")
  inspect(muts[0].ref_base, content="A")
  inspect(muts[0].alt_base, content="T")
  // pos 40 上下文: ref[39]='T', ref[40]='A', ref[41]='C'
  inspect(muts[0].trinucleotide_context, content="TAC")
  let matrix = @src.build_mutation_matrix(muts, ["TUMOR"])
  inspect(matrix.total_mutations[0], content="1")
}

///|
test "mutect: JSON 导出与 tumor-only 模式" {
  let ref_ = ml_reference()
  let tumor = @src.MutectSample::{
    name: "TUMOR",
    reads: ml_make_reads(ref_, 20, 40, 'T', 8),
  }
  let normal = @src.MutectSample::{ name: "NORMAL", reads: [], }
  let calls = @src.mutect_call("chr1", ref_, tumor, normal)
  inspect(calls.length(), content="1")
  inspect(calls[0].filter, content="PASS")
  inspect(calls[0].normal_alt, content="0")
  let json = @src.mutect_json(calls)
  assert_eq(json.contains("\"filter\""), true)
  assert_eq(json.contains("\"lod_t\""), true)
  let line = @src.mutect_summary_line(@src.mutect_summary(calls))
  assert_eq(line.contains("pass=1"), true)
}
