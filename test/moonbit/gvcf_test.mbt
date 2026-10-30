// ===========================================================================
// gvcf_test.mbt — gVCF 参考块 / 多样本合并 / 联合基因分型 测试.
// ===========================================================================

///|
fn gvcf_het_reads(
  reference : String,
  start : Int32,
  offset : Int32,
  alt : Char,
  n : Int32,
  alt_share : Int32,
) -> Array[@src.GatkRead] {
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < n; i = i + 1 {
    let off : Int32? = if i < alt_share { Some(offset) } else { None }
    reads.push(
      gatk_make_read(reference, start, 20, off, alt, "r" + i.to_string()),
    )
  }
  reads
}

///|
fn gvcf_sample_a(reference : String) -> Array[@src.GatkVariant] {
  // 覆盖 [10, 30): 位置 16 处杂合.
  let reads = gvcf_het_reads(reference, 10, 5, 'C', 10, 5)
  @src.gatk_call_variants("chr1", reference, reads)
}

///|
fn gvcf_sample_b(reference : String) -> Array[@src.GatkVariant] {
  // 覆盖 [20, 40): 位置 25 处杂合 (ref T -> C).
  let reads = gvcf_het_reads(reference, 20, 4, 'C', 10, 5)
  @src.gatk_call_variants("chr1", reference, reads)
}

// ===========================================================================
// 建库: 变异记录 + 参考块
// ===========================================================================

///|
test "gvcf_build_emits_variant_and_reference_blocks" {
  let reference = gatk_ref()
  let variants = gvcf_sample_a(reference)
  assert_eq(variants.length(), 1)
  let callset = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    variants,
    sample="S1",
  ) catch {
    @src.GvcfError(message) => abort("build failed: " + message)
  }
  assert_eq(callset.chrom, "chr1")
  assert_eq(callset.samples.length(), 1)
  // 前块 [1,15] + 变异(16) + 后块 [17,60].
  assert_eq(callset.records.length(), 3)
  let first = callset.records[0]
  assert_true(first.is_block)
  assert_eq(first.pos, 1)
  assert_eq(first.block_end, 15)
  assert_eq(first.alts, ["<NON_REF>"])
  let variant = callset.records[1]
  assert_true(!variant.is_block)
  assert_eq(variant.pos, 16)
  assert_eq(variant.ref_allele, "G")
  assert_eq(variant.alts, ["C"])
  assert_eq(variant.samples.get("S1").unwrap().gt, "0/1")
  let last = callset.records[2]
  assert_true(last.is_block)
  assert_eq(last.pos, 17)
  assert_eq(last.block_end, reference.length())
}

///|
test "gvcf_build_rejects_empty_reference" {
  let raised = try {
    ignore(@src.gvcf_build_from_variants("chr1", "", []))
    false
  } catch {
    @src.GvcfError(_) => true
  }
  assert_true(raised)
}

// ===========================================================================
// 基因型似然
// ===========================================================================

///|
test "gvcf_pl_from_ad_genotypes_by_depth" {
  assert_eq(@src.gvcf_gt_from_pl(@src.gvcf_pl_from_ad(20, 0)), "0/0")
  assert_eq(@src.gvcf_gt_from_pl(@src.gvcf_pl_from_ad(10, 10)), "0/1")
  assert_eq(@src.gvcf_gt_from_pl(@src.gvcf_pl_from_ad(1, 20)), "1/1")
  // 最小 PL 归一为 0.
  let pl = @src.gvcf_pl_from_ad(10, 10)
  assert_eq(pl.length(), 3)
  assert_true(@src.gvcf_gq_from_pl(pl) > 0)
}

// ===========================================================================
// 多样本合并 + 联合基因分型
// ===========================================================================

///|
test "gvcf_merge_assigns_hom_ref_to_unaffected_sample" {
  let reference = gatk_ref()
  let cs_a = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_a(reference),
    sample="S1",
  ) catch {
    @src.GvcfError(_) => abort("build A failed")
  }
  let cs_b = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_b(reference),
    sample="S2",
  ) catch {
    @src.GvcfError(_) => abort("build B failed")
  }
  let merged = @src.gvcf_merge([cs_a, cs_b]) catch {
    @src.GvcfError(message) => abort("merge failed: " + message)
  }
  assert_eq(merged.samples, ["S1", "S2"])
  assert_eq(merged.records.length(), 2)
  let at16 = merged.records[0]
  assert_eq(at16.pos, 16)
  assert_eq(at16.samples.get("S1").unwrap().gt, "0/1")
  assert_eq(at16.samples.get("S2").unwrap().gt, "0/0")
  let at25 = merged.records[1]
  assert_eq(at25.pos, 25)
  assert_eq(at25.samples.get("S1").unwrap().gt, "0/0")
  assert_eq(at25.samples.get("S2").unwrap().gt, "0/1")
}

///|
test "gvcf_merge_rejects_empty_input" {
  let empty : Array[@src.GvcfCallset] = []
  let raised = try {
    ignore(@src.gvcf_merge(empty))
    false
  } catch {
    @src.GvcfError(_) => true
  }
  assert_true(raised)
}

///|
test "gvcf_joint_genotype_recomputes_from_pl" {
  let reference = gatk_ref()
  let cs_a = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_a(reference),
    sample="S1",
  ) catch {
    @src.GvcfError(_) => abort("build A failed")
  }
  let cs_b = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_b(reference),
    sample="S2",
  ) catch {
    @src.GvcfError(_) => abort("build B failed")
  }
  let merged = @src.gvcf_merge([cs_a, cs_b]) catch {
    @src.GvcfError(_) => abort("merge failed")
  }
  let joint = @src.gvcf_joint_genotype(merged)
  assert_eq(joint.samples, ["S1", "S2"])
  let at16 = joint.records[0]
  assert_eq(at16.samples.get("S1").unwrap().gt, "0/1")
  assert_eq(at16.samples.get("S2").unwrap().gt, "0/0")
}

// ===========================================================================
// 解析 / VCF 输出 / JSON
// ===========================================================================

///|
fn gvcf_demo_text() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
  "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=4\tGT:DP:GQ:PL\t0/0:0:60:0,60,90\n" +
  "chr1\t5\t.\tA\tC\t.\t.\t.\tGT:DP:GQ:PL\t0/1:10:30:30,0,30\n"
}

///|
test "gvcf_parse_reads_blocks_and_samples" {
  match @src.gvcf_parse(gvcf_demo_text()) {
    Some(callset) => {
      assert_eq(callset.chrom, "chr1")
      assert_eq(callset.samples, ["S1"])
      assert_eq(callset.records.length(), 2)
      let block = callset.records[0]
      assert_true(block.is_block)
      assert_eq(block.block_end, 4)
      let variant = callset.records[1]
      assert_true(!variant.is_block)
      assert_eq(variant.pos, 5)
      assert_eq(variant.alts, ["C"])
      let geno = variant.samples.get("S1").unwrap()
      assert_eq(geno.gt, "0/1")
      assert_eq(geno.pl, [30, 0, 30])
      assert_eq(geno.dp, 10)
    }
    None => assert_true(false)
  }
}

///|
test "gvcf_parse_rejects_non_vcf" {
  assert_true(@src.gvcf_parse("") is None)
  assert_true(@src.gvcf_parse("hello world") is None)
}

///|
test "gvcf_json_roundtrip_preserves_records" {
  let reference = gatk_ref()
  let callset = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_a(reference),
    sample="S1",
  ) catch {
    @src.GvcfError(_) => abort("build failed")
  }
  let json = @src.gvcf_callset_to_json(callset)
  assert_true(json.contains("\"is_block\":true"))
  match @src.gvcf_callset_from_json(json) {
    Some(restored) => {
      assert_eq(restored.chrom, callset.chrom)
      assert_eq(restored.samples, callset.samples)
      assert_eq(restored.records.length(), callset.records.length())
      assert_eq(restored.records[1].ref_allele, callset.records[1].ref_allele)
      assert_eq(restored.records[1].samples.get("S1").unwrap().gt, "0/1")
    }
    None => assert_true(false)
  }
}

///|
test "gvcf_write_vcf_has_multi_sample_columns" {
  let reference = gatk_ref()
  let cs_a = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_a(reference),
    sample="S1",
  ) catch {
    @src.GvcfError(_) => abort("build A failed")
  }
  let cs_b = @src.gvcf_build_from_variants(
    "chr1",
    reference,
    gvcf_sample_b(reference),
    sample="S2",
  ) catch {
    @src.GvcfError(_) => abort("build B failed")
  }
  let joint = @src.gvcf_joint_genotype(
    @src.gvcf_merge([cs_a, cs_b]) catch {
      @src.GvcfError(_) => abort("merge failed")
    },
  )
  let vcf = @src.gvcf_write_vcf(joint)
  assert_true(
    vcf.contains(
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    ),
  )
  assert_true(vcf.contains("0/1"))
  assert_true(@src.gvcf_summary(joint).contains("samples=2"))
}

// ===========================================================================
// ALT 索引重映射(同位点不同 ALT 的多样本合并)
// ===========================================================================

///|
fn gvcf_variant_record(
  pos : Int32,
  ref_allele : String,
  alts : Array[String],
  sample : String,
  gt : String,
  pl : Array[Int32],
) -> @src.GvcfRecord {
  @src.GvcfRecord::{
    chrom: "chr1",
    pos,
    ref_allele,
    alts,
    is_block: false,
    block_end: pos,
    samples: Map([(sample, @src.GvcfGenotype::{ gt, pl, dp: 10, gq: 30, })]),
  }
}

///|
test "gvcf_merge_remaps_alt_indices_per_sample" {
  // S1: pos 10 A->C 杂合; S2: pos 10 A->G 杂合(同位点不同 ALT).
  let cs_a = @src.GvcfCallset::{
    chrom: "chr1",
    samples: ["S1"],
    records: [gvcf_variant_record(10, "A", ["C"], "S1", "0/1", [30, 0, 30])],
  }
  let cs_b = @src.GvcfCallset::{
    chrom: "chr1",
    samples: ["S2"],
    records: [gvcf_variant_record(10, "A", ["G"], "S2", "0/1", [30, 0, 30])],
  }
  let merged = @src.gvcf_merge([cs_a, cs_b]) catch {
    @src.GvcfError(message) => abort("merge failed: " + message)
  }
  assert_eq(merged.records.length(), 1)
  let record = merged.records[0]
  assert_eq(record.alts, ["C", "G"])
  // S1 的 C 映射到合并索引 1; S2 的 G 映射到合并索引 2(而非被错读为 C).
  assert_eq(record.samples.get("S1").unwrap().gt, "0/1")
  assert_eq(record.samples.get("S2").unwrap().gt, "0/2")
  // 等位数变化(1 -> 2)时 PL 基因型空间无法对应, 整字段置缺失.
  assert_eq(record.samples.get("S1").unwrap().pl.length(), 0)
  assert_eq(record.samples.get("S2").unwrap().pl.length(), 0)
  // 联合分型不再按错位的 PL 把 S2 判成 REF/C (0/1).
  let joint = @src.gvcf_joint_genotype(merged)
  assert_eq(joint.records[0].samples.get("S1").unwrap().gt, "0/1")
  assert_eq(joint.records[0].samples.get("S2").unwrap().gt, "0/2")
}

///|
test "gvcf_merge_remaps_pl_when_alt_sets_match" {
  // 三个样本同位点: S1 ALT=[C] 1/1; S2 ALT=[A] 0/1; S3 ALT=[C,A] 1/2.
  let cs_a = @src.GvcfCallset::{
    chrom: "chr1",
    samples: ["S1"],
    records: [gvcf_variant_record(20, "T", ["C"], "S1", "1/1", [90, 30, 0])],
  }
  let cs_b = @src.GvcfCallset::{
    chrom: "chr1",
    samples: ["S2"],
    records: [gvcf_variant_record(20, "T", ["A"], "S2", "0/1", [30, 0, 30])],
  }
  let cs_c = @src.GvcfCallset::{
    chrom: "chr1",
    samples: ["S3"],
    records: [
      gvcf_variant_record(20, "T", ["C", "A"], "S3", "1/2", [
        90, 30, 40, 20, 0, 60,
      ]),
    ],
  }
  let merged = @src.gvcf_merge([cs_a, cs_b, cs_c]) catch {
    @src.GvcfError(message) => abort("merge failed: " + message)
  }
  let record = merged.records[0]
  assert_eq(record.alts, ["A", "C"])
  // S1: C -> 合并索引 2, GT 1/1 -> 2/2; 等位数 1 != 2, PL 置缺失.
  assert_eq(record.samples.get("S1").unwrap().gt, "2/2")
  assert_eq(record.samples.get("S1").unwrap().pl.length(), 0)
  // S2: A -> 合并索引 1, GT 0/1 保持.
  assert_eq(record.samples.get("S2").unwrap().gt, "0/1")
  // S3: 原 ALT [C,A] 与合并 [A,C] 等位数一致, PL 双射重排:
  // 原空间(0,C=1,A=2) -> 新空间(0,A=1,C=2).
  let geno = record.samples.get("S3").unwrap()
  assert_eq(geno.gt, "1/2")
  assert_eq(geno.pl, [90, 20, 60, 30, 0, 40])
  // 联合分型按重排后的 6 项 PL 重选基因型, 最小项对应 C/A 杂合(1/2).
  let joint = @src.gvcf_joint_genotype(merged)
  assert_eq(joint.records[0].samples.get("S3").unwrap().gt, "1/2")
  assert_eq(joint.records[0].samples.get("S1").unwrap().gt, "2/2")
}

///|
test "gvcf_gt_from_pl_supports_multiallelic" {
  assert_eq(@src.gvcf_gt_from_pl([0, 60, 90]), "0/0")
  assert_eq(@src.gvcf_gt_from_pl([30, 0, 30]), "0/1")
  assert_eq(@src.gvcf_gt_from_pl([90, 30, 0]), "1/1")
  assert_eq(@src.gvcf_gt_from_pl([90, 20, 60, 30, 0, 40]), "1/2")
  assert_eq(@src.gvcf_gt_from_pl([90, 30, 40, 0, 20, 60]), "0/2")
  assert_eq(@src.gvcf_gt_from_pl([90, 30, 40, 20, 60, 0]), "2/2")
  assert_eq(@src.gvcf_gt_from_pl([10, 20]), "./.")
  assert_eq(@src.gvcf_gt_from_pl([1, 2, 3, 4]), "./.")
}
