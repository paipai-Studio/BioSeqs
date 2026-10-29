// ===========================================================================
// gatk_lite_test.mbt — GATK 风格变异检出测试
// ===========================================================================

///|
/// 参考 60bp, 非周期(避免 indel 错位巧合).
fn gatk_ref() -> String {
  "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGT"
}

///|
/// 从参考切片构造读段, 并在指定读内偏移注入替换碱基.
fn gatk_make_read(
  reference : String,
  pos : Int32,
  len : Int32,
  offset : Int32?,
  alt : Char,
  name : String,
) -> @src.GatkRead {
  let end = if pos + len < reference.length() {
    pos + len
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
test "gatk pileup site counts" {
  let reference = gatk_ref()
  // 10 条 read 覆盖 [10, 30), 位置 15 处 5 条 alt(读内偏移 5)
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let site = @src.gatk_pileup_site(reference, reads, 15)
  assert_eq(site.depth, 10)
  assert_eq(site.ref_base, 'G')
  let mut g_count = 0
  let mut c_count = 0
  for pair in site.base_counts {
    if pair.0 == 'G' {
      g_count = pair.1
    }
    if pair.0 == 'C' {
      c_count = pair.1
    }
  }
  assert_eq(g_count, 5)
  assert_eq(c_count, 5)
  assert_eq(site.insertions.length(), 0)
  assert_eq(site.deletions.length(), 0)
}

///|
test "gatk pileup no coverage" {
  let reference = gatk_ref()
  let site = @src.gatk_pileup_site(reference, [], 5)
  assert_eq(site.depth, 0)
  assert_eq(site.base_counts.length(), 0)
}

///|
test "gatk genotype likelihoods heterozygous" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let site = @src.gatk_pileup_site(reference, reads, 15)
  let gts = @src.gatk_genotype_likelihoods(site, 'C')
  assert_eq(gts[0].gt, "0/1")
  assert_true(gts[0].posterior > 0.5)
  let mut total = 0.0
  for g in gts {
    total = total + g.posterior
  }
  assert_true(total > 0.999 && total < 1.001)
  assert_true(gts[0].log10_error > 0.0)
}

///|
test "gatk genotype likelihoods homozygous alt" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 8; i = i + 1 {
    reads.push(
      gatk_make_read(reference, 0, 20, Some(5), 'G', "r" + i.to_string()),
    )
  }
  // 位置 5: 参考碱基 A("ATCGGA" 位置 5)
  let site = @src.gatk_pileup_site(reference, reads, 5)
  assert_eq(site.ref_base, 'A')
  let gts = @src.gatk_genotype_likelihoods(site, 'G')
  assert_eq(gts[0].gt, "1/1")
  assert_true(gts[0].log10_error > 20.0)
}

///|
test "gatk call heterozygous snp" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 1)
  let v = variants[0]
  assert_eq(v.chrom, "chr1")
  assert_eq(v.pos, 16)
  assert_eq(v.ref, "G")
  assert_eq(v.alt, "C")
  assert_eq(v.gt, "0/1")
  assert_eq(v.dp, 10)
  assert_eq(v.ad_ref, 5)
  assert_eq(v.ad_alt, 5)
  assert_true(v.gq >= 20)
  assert_true(v.af > 0.49 && v.af < 0.51)
}

///|
test "gatk call homozygous snp" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    // 位置 7(读内偏移 7): 参考碱基 T("ACGTACGT" 位置 7), alt A
    reads.push(
      gatk_make_read(reference, 0, 20, Some(7), 'A', "r" + i.to_string()),
    )
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 1)
  assert_eq(variants[0].gt, "1/1")
  assert_eq(variants[0].ref, "T")
  assert_eq(variants[0].alt, "A")
  assert_eq(variants[0].ad_alt, 10)
}

///|
test "gatk no variant on pure reference" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    reads.push(
      gatk_make_read(reference, 10, 20, None, 'N', "r" + i.to_string()),
    )
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 0)
}

///|
test "gatk low allele frequency filtered" {
  let reference = gatk_ref()
  // 10 条中仅 1 条 alt(10% < 阈值 20%)
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i == 0 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 0)
}

///|
test "gatk low depth filtered" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = [
    gatk_make_read(reference, 0, 20, Some(5), 'G', "r0"),
    gatk_make_read(reference, 0, 20, Some(5), 'G', "r1"),
  ]
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  assert_eq(variants.length(), 0)
}

///|
test "gatk active regions detection" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  // 位置 10-14 集中错配: 读内偏移 10-14 全替换
  for i = 0; i < 3; i = i + 1 {
    let sliced = reference[0:20].to_owned()
    // 偏移 10-14 替换为 T(参考该处 ACGTA)
    let mutated = sliced[0:10].to_owned() + "TTTTT" + sliced[15:].to_owned()
    reads.push({ name: "r" + i.to_string(), seq: mutated, pos: 0, })
  }
  let regions = @src.gatk_active_regions(reference, reads, window=10)
  assert_true(regions.length() > 0)
  let (s, e) = regions[0]
  assert_true(s <= 10)
  assert_true(e > 10)
}

///|
test "gatk assemble consensus majority vote" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 5; i = i + 1 {
    let off : Int32? = if i < 4 { Some(2) } else { None }
    // 位置 2: 参考 G(2%4=2), alt T
    reads.push(gatk_make_read(reference, 0, 15, off, 'T', "r" + i.to_string()))
  }
  let consensus = @src.gatk_assemble_consensus(reference, reads, 0, 15)
  assert_eq(consensus[2:3].to_owned(), "T")
}

///|
test "gatk call insertion variant" {
  // 非周期参考(避免周期序列的 indel 错位巧合)
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCAT"
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 8; i = i + 1 {
    let seq = if i < 6 {
      // 读内偏移 5 后插入 GG(参考位置 5 后)
      let sliced = reference[0:20].to_owned()
      sliced[0:6].to_owned() + "GG" + sliced[6:].to_owned()
    } else {
      reference[0:20].to_owned()
    }
    reads.push({ name: "r" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let mut found_ins = false
  for v in variants {
    if v.alt.length() > v.ref.length() {
      found_ins = true
      assert_true(v.alt.length() - v.ref.length() == 2)
    }
  }
  assert_true(found_ins)
}

///|
test "gatk call deletion variant" {
  // 非周期参考(避免周期序列的 indel 错位巧合)
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCAT"
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 8; i = i + 1 {
    let seq = if i < 6 {
      // 读内偏移 5-7 删除(参考位置 5,6,7 的 ATT)
      let sliced = reference[0:20].to_owned()
      sliced[0:5].to_owned() + sliced[8:].to_owned()
    } else {
      reference[0:20].to_owned()
    }
    reads.push({ name: "r" + i.to_string(), seq, pos: 0, })
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let mut found_del = false
  for v in variants {
    if v.ref.length() > v.alt.length() {
      found_del = true
      assert_true(v.ref.length() - v.alt.length() >= 1)
    }
  }
  assert_true(found_del)
}

///|
test "gatk write vcf and parse_vcf roundtrip" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let vcf_text = @src.gatk_write_vcf("chr1", reference, variants)
  let vcf = @src.parse_vcf(vcf_text)
  assert_eq(vcf.len(), 1)
  let rec = vcf.records[0]
  assert_eq(rec.chrom, "chr1")
  assert_eq(rec.pos, 16)
  assert_eq(rec.ref_seq, "G")
  assert_eq(rec.alt[0], "C")
  assert_true(rec.is_snp())
  assert_eq(rec.info.get("DP").unwrap(), "10")
  let gts = rec.samples.get("SAMPLE").unwrap()
  assert_true(gts.length() > 0)
  assert_true(gts[0].contains("0/1"))
}

///|
test "gatk variants summary" {
  let reference = gatk_ref()
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let off : Int32? = if i < 5 { Some(5) } else { None }
    reads.push(gatk_make_read(reference, 10, 20, off, 'C', "r" + i.to_string()))
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let s = @src.gatk_variants_summary(variants)
  assert_true(s.contains("gatk-lite"))
  assert_true(s.contains("1 SNP"))
}
