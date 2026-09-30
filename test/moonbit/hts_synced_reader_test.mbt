///|
fn mk_variant(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alts : Array[String]
) -> @src.SyncVariant {
  @src.SyncVariant::make(chrom, pos, ref_allele, alts, ".")
}

test "synced_reader regions init string basic" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000,chr2:5000-6000")
  assert_eq(r.nseqs(), 2)
  assert_eq(r.nregs_for_seq("chr1"), 1)
  assert_eq(r.nregs_for_seq("chr2"), 1)
}

test "synced_reader regions init whole chrom" {
  let r = @src.bcf_sr_regions_init_string("chr1,chr2,chr3")
  assert_eq(r.nseqs(), 3)
  assert_eq(r.nregs_for_seq("chr1"), 1)
  assert_eq(r.nregs_for_seq("chr2"), 1)
  assert_eq(r.nregs_for_seq("chr3"), 1)
}

test "synced_reader regions init single pos" {
  let r = @src.bcf_sr_regions_init_string("chr1:500")
  assert_eq(r.nseqs(), 1)
  assert_eq(r.nregs_for_seq("chr1"), 1)
  let regs = r.get_regs("chr1")
  assert_eq(regs.length(), 1)
  // 1-based 500 → 0-based 499
  assert_eq(regs[0].0, 499)
  assert_eq(regs[0].1, 499)
}

test "synced_reader regions init range" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000")
  let regs = r.get_regs("chr1")
  assert_eq(regs.length(), 1)
  // 1-based → 0-based
  assert_eq(regs[0].0, 999)
  assert_eq(regs[0].1, 1999)
}

test "synced_reader regions merge overlapping" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000,chr1:1500-3000")
  let regs = r.get_regs("chr1")
  // Should merge into one interval
  assert_eq(regs.length(), 1)
  assert_eq(regs[0].0, 999)
  assert_eq(regs[0].1, 2999)
}

test "synced_reader regions merge adjacent" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000,chr1:2000-3000")
  let regs = r.get_regs("chr1")
  // 2000 and 2001 are adjacent (0-based: 1999 and 2000), should merge
  assert_eq(regs.length(), 1)
}

test "synced_reader regions no merge gap" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000,chr1:3000-4000")
  let regs = r.get_regs("chr1")
  // Gap between 2000 and 3000, should NOT merge
  assert_eq(regs.length(), 2)
}

test "synced_reader regions overlap yes" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000")
  // 0-based: query [999, 1999]
  assert_eq(r.overlap("chr1", 1500, 1600), 0)
  assert_eq(r.overlap("chr1", 999, 999), 0)
  assert_eq(r.overlap("chr1", 1999, 1999), 0)
}

test "synced_reader regions overlap no" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000")
  assert_eq(r.overlap("chr1", 500, 998), -1)
  assert_eq(r.overlap("chr1", 2000, 3000), -1)
}

test "synced_reader regions overlap wrong seq" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-2000")
  assert_eq(r.overlap("chr2", 1500, 1600), -1)
}

test "synced_reader regions weird name" {
  let r = @src.bcf_sr_regions_init_string("{chr-1}:1000-2000")
  assert_eq(r.nseqs(), 1)
  assert_eq(r.nregs_for_seq("chr-1"), 1)
}

test "synced_reader regions open ended" {
  let r = @src.bcf_sr_regions_init_string("chr1:1000-")
  let regs = r.get_regs("chr1")
  assert_eq(regs.length(), 1)
  assert_eq(regs[0].0, 999)
  // End should be max_csi_coor
  assert_true(regs[0].1 > 1000000000)
}

test "synced_reader classify allele snp" {
  assert_eq(@src.classify_allele("A", "T"), @src.sr_snp)
  assert_eq(@src.classify_allele("G", "C"), @src.sr_snp)
}

test "synced_reader classify allele indel" {
  assert_eq(@src.classify_allele("A", "AT"), @src.sr_indel)
  assert_eq(@src.classify_allele("AT", "A"), @src.sr_indel)
  assert_eq(@src.classify_allele("A", "ATT"), @src.sr_indel)
}

test "synced_reader classify allele ref" {
  assert_eq(@src.classify_allele("A", "."), @src.sr_ref)
  assert_eq(@src.classify_allele("A", ""), @src.sr_ref)
}

test "synced_reader classify variant" {
  let vt = @src.classify_variant("A", ["T"])
  assert_true((vt & @src.sr_snp) != 0)
  let vt2 = @src.classify_variant("A", ["AT"])
  assert_true((vt2 & @src.sr_indel) != 0)
  let vt3 = @src.classify_variant("A", ["T", "AT"])
  assert_true((vt3 & @src.sr_snp) != 0)
  assert_true((vt3 & @src.sr_indel) != 0)
}

test "synced_reader pair compatible snp snp" {
  let t = @src.sr_snp
  assert_true(@src.pair_compatible(t, t, @src.bcf_sr_pair_snps))
}

test "synced_reader pair compatible indel indel" {
  let t = @src.sr_indel
  assert_true(@src.pair_compatible(t, t, @src.bcf_sr_pair_indels))
}

test "synced_reader pair not compatible snp indel" {
  let ts = @src.sr_snp
  let ti = @src.sr_indel
  assert_false(@src.pair_compatible(ts, ti, @src.bcf_sr_pair_snps))
}

test "synced_reader pair compatible any" {
  assert_true(
    @src.pair_compatible(@src.sr_snp, @src.sr_indel, @src.bcf_sr_pair_any)
  )
}

test "synced_reader pair compatible exact" {
  let t = @src.sr_snp
  assert_true(
    @src.pair_compatible(t, t, @src.bcf_sr_pair_exact)
  )
  assert_false(
    @src.pair_compatible(@src.sr_snp, @src.sr_indel, @src.bcf_sr_pair_exact)
  )
}

test "synced_reader init and add" {
  let sr = @src.bcf_sr_init()
  assert_eq(sr.nreaders(), 0)
  let sr2 = sr.add_reader("file1.vcf", [])
  assert_eq(sr2.nreaders(), 1)
}

test "synced_reader next_line single file" {
  let sr = @src.bcf_sr_init()
  let variants = [
    mk_variant("chr1", 100, "A", ["T"]),
    mk_variant("chr1", 200, "G", ["C"]),
  ]
  let sr2 = sr.add_reader("file1.vcf", variants)
  let (sr3, count1) = sr2.next_line()
  assert_eq(count1, 1)
  assert_true(sr3.has_line_at(0))
  let line1 = sr3.get_line(0)
  assert_eq(line1.unwrap_or(mk_variant("", 0, "", [])).pos, 100)
  let (sr4, count2) = sr3.next_line()
  assert_eq(count2, 1)
  let line2 = sr4.get_line(0)
  assert_eq(line2.unwrap_or(mk_variant("", 0, "", [])).pos, 200)
  let (_, count3) = sr4.next_line()
  assert_eq(count3, 0)
}

test "synced_reader next_line two files sync" {
  let sr = @src.bcf_sr_init()
  let v1 = [
    mk_variant("chr1", 100, "A", ["T"]),
    mk_variant("chr1", 200, "G", ["C"]),
    mk_variant("chr1", 300, "C", ["A"]),
  ]
  let v2 = [
    mk_variant("chr1", 100, "A", ["T"]),
    mk_variant("chr1", 300, "C", ["A"]),
  ]
  let sr2 = sr.add_reader("f1.vcf", v1).add_reader("f2.vcf", v2)
  // pos 100: both files have it
  let (sr3, c1) = sr2.next_line()
  assert_eq(c1, 2)
  assert_true(sr3.has_line_at(0))
  assert_true(sr3.has_line_at(1))
  // pos 200: only file1
  let (sr4, c2) = sr3.next_line()
  assert_eq(c2, 1)
  assert_true(sr4.has_line_at(0))
  assert_false(sr4.has_line_at(1))
  // pos 300: both files
  let (sr5, c3) = sr4.next_line()
  assert_eq(c3, 2)
  assert_true(sr5.has_line_at(0))
  assert_true(sr5.has_line_at(1))
  // done
  let (_, c4) = sr5.next_line()
  assert_eq(c4, 0)
}

test "synced_reader next_line different chroms" {
  let sr = @src.bcf_sr_init()
  let v1 = [
    mk_variant("chr1", 100, "A", ["T"]),
    mk_variant("chr2", 50, "G", ["C"]),
  ]
  let v2 = [
    mk_variant("chr2", 50, "G", ["C"]),
    mk_variant("chr3", 10, "C", ["A"]),
  ]
  let sr2 = sr.add_reader("f1.vcf", v1).add_reader("f2.vcf", v2)
  // chr1:100 only in f1
  let (sr3, c1) = sr2.next_line()
  assert_eq(c1, 1)
  // chr2:50 in both
  let (sr4, c2) = sr3.next_line()
  assert_eq(c2, 2)
  // chr3:10 only in f2
  let (sr5, c3) = sr4.next_line()
  assert_eq(c3, 1)
  // done
  let (_, c4) = sr5.next_line()
  assert_eq(c4, 0)
}

test "synced_reader empty readers" {
  let sr = @src.bcf_sr_init()
  let (_, count) = sr.next_line()
  assert_eq(count, 0)
}

test "synced_reader empty variant list" {
  let sr = @src.bcf_sr_init()
  let sr2 = sr.add_reader("empty.vcf", [])
  let (_, count) = sr2.next_line()
  assert_eq(count, 0)
}

test "synced_reader get_line none" {
  let sr = @src.bcf_sr_init()
  let sr2 = sr.add_reader("f.vcf", [mk_variant("chr1", 100, "A", ["T"])])
  // Before next_line, get_line should return None
  assert_true(sr2.get_line(0).is_none())
}
