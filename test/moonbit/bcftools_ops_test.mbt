///|
/// Tests for bcftools-style variant operations.

fn make_test_vcf() -> @src.VcfFile {
  let text =
    #|##fileformat=VCFv4.2
    #|##FILTER=<ID=LowQual,Description="Low quality">
    #|##FILTER=<ID=PASS,Description="All filters passed">
    #|##contig=<ID=chr1,length=1000>
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
    #|chr1	100	.	A	G	50	PASS	DP=10	GT	0/0	0/1
    #|chr1	200	.	AT	A	30	LowQual	DP=5	GT	0/1	1/1
    #|chr1	300	.	C	A,T	99	PASS	DP=20;AF=0.5	GT	0/1	1/2
    #|chr1	400	.	G	GGT	99	PASS	DP=15	GT	0/1	./.
    #|chr1	500	.	TT	T	.	PASS	DP=8	GT	0/0	0/0
  @src.parse_vcf(text)
}

test "vcf_filter_pass" {
  let vcf = make_test_vcf()
  let filtered = @src.vcf_filter_pass(vcf)
  // 4 records have PASS, 1 has LowQual
  assert_eq(filtered.records.length(), 4)
}

test "vcf_filter_min_qual" {
  let vcf = make_test_vcf()
  let filtered = @src.vcf_filter_min_qual(vcf, 40.0)
  // qual >= 40: 100(50), 300(99), 400(99) → 3 records (500 has no qual)
  assert_eq(filtered.records.length(), 3)
}

test "vcf_filter_tag" {
  let vcf = make_test_vcf()
  let filtered = @src.vcf_filter_tag(vcf, "LowQual")
  assert_eq(filtered.records.length(), 1)
  assert_eq(filtered.records[0].pos, 200)
}

test "vcf_exclude_tag" {
  let vcf = make_test_vcf()
  let filtered = @src.vcf_exclude_tag(vcf, "LowQual")
  assert_eq(filtered.records.length(), 4)
}

test "vcf_view_snps" {
  let vcf = make_test_vcf()
  let snps = @src.vcf_view_snps(vcf)
  // pos 100 (A>G) and 300 (C>A,T) — 300 is multi-allelic but each alt is 1bp
  // is_snp checks ref.length==1 && alt.length==1 && alt[0].length==1
  // 300 has alt.length==2 so not SNP
  // 500: ref=TT (len 2), not SNP
  // 400: ref=G, alt=GGT → insertion, not SNP
  assert_eq(snps.records.length(), 1)
  assert_eq(snps.records[0].pos, 100)
}

test "vcf_view_indels" {
  let vcf = make_test_vcf()
  let indels = @src.vcf_view_indels(vcf)
  // pos 200 (AT>A): deletion
  // pos 400 (G>GGT): insertion
  // pos 500 (TT>T): deletion
  assert_eq(indels.records.length(), 3)
}

test "vcf_norm_split" {
  let vcf = make_test_vcf()
  let split = @src.vcf_norm_split(vcf)
  // pos 300 has 2 alt alleles → split into 2
  // total: 4 + 2 = 6
  assert_eq(split.records.length(), 6)
  // Check the split records
  let r300 = split.records[2]
  assert_eq(r300.pos, 300)
  assert_eq(r300.alt.length(), 1)
  assert_eq(r300.alt[0], "A")
  let r300b = split.records[3]
  assert_eq(r300b.pos, 300)
  assert_eq(r300b.alt[0], "T")
}

test "vcf_norm_merge" {
  // Create a VCF with two records at same pos/ref that can be merged
  let text =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	100	.	A	G	50	PASS	DP=10
    #|chr1	100	.	A	T	50	PASS	DP=10
    #|chr1	200	.	C	A	30	PASS	DP=5
  let vcf = @src.parse_vcf(text)
  let merged = @src.vcf_norm_merge(vcf)
  // Two records at chr1:100:A should merge into one multi-allelic
  assert_eq(merged.records.length(), 2)
  assert_eq(merged.records[0].alt.length(), 2)
  assert_eq(merged.records[0].alt[0], "G")
  assert_eq(merged.records[0].alt[1], "T")
}

test "vcf_isec" {
  let vcf_a = make_test_vcf()
  // Make B with only pos 100 and 300
  let text_b =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	100	.	A	G	50	PASS	DP=10
    #|chr1	300	.	C	A,T	99	PASS	DP=20
  let vcf_b = @src.parse_vcf(text_b)
  let isec = @src.vcf_isec(vcf_a, vcf_b)
  assert_eq(isec.records.length(), 2)
  assert_eq(isec.records[0].pos, 100)
  assert_eq(isec.records[1].pos, 300)
}

test "vcf_isec_complement" {
  let vcf_a = make_test_vcf()
  let text_b =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	100	.	A	G	50	PASS	DP=10
    #|chr1	300	.	C	A,T	99	PASS	DP=20
  let vcf_b = @src.parse_vcf(text_b)
  let comp = @src.vcf_isec_complement(vcf_a, vcf_b)
  // A has 5 records, 2 in B → 3 not in B
  assert_eq(comp.records.length(), 3)
}

test "vcf_union" {
  let vcf_a = make_test_vcf()
  let text_b =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	100	.	A	G	50	PASS	DP=10
    #|chr1	600	.	G	C	40	PASS	DP=12
  let vcf_b = @src.parse_vcf(text_b)
  let union = @src.vcf_union(vcf_a, vcf_b)
  // A has 5, B has 2 (1 shared) → 6 unique
  assert_eq(union.records.length(), 6)
}

test "vcf_sort" {
  // Create unsorted VCF
  let text =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	300	.	C	A	30	PASS	DP=5
    #|chr1	100	.	A	G	50	PASS	DP=10
    #|chr2	50	.	G	T	20	PASS	DP=3
    #|chr1	200	.	T	C	40	PASS	DP=8
  let vcf = @src.parse_vcf(text)
  let sorted = @src.vcf_sort(vcf)
  assert_eq(sorted.records[0].chrom, "chr1")
  assert_eq(sorted.records[0].pos, 100)
  assert_eq(sorted.records[1].chrom, "chr1")
  assert_eq(sorted.records[1].pos, 200)
  assert_eq(sorted.records[2].chrom, "chr1")
  assert_eq(sorted.records[2].pos, 300)
  assert_eq(sorted.records[3].chrom, "chr2")
  assert_eq(sorted.records[3].pos, 50)
}

test "vcf_query_basic" {
  let vcf = make_test_vcf()
  let query = @src.vcf_query_basic(vcf)
  let lines = query.to_string().split("\n").to_array()
  // Last line is empty (trailing \n)
  assert_eq(lines.length(), 6) // 5 records + 1 empty
  assert_true(lines[0].contains("chr1	100	A	G"))
}

test "vcf_query_info" {
  let vcf = make_test_vcf()
  let dp = @src.vcf_query_info(vcf, "DP")
  let lines = dp.to_string().split("\n").to_array()
  // All 5 records have DP
  assert_eq(lines.length(), 6) // 5 + 1 empty
  assert_eq(lines[0], "10")
  assert_eq(lines[1], "5")
}

test "vcf_stats" {
  let vcf = make_test_vcf()
  let stats = @src.vcf_stats(vcf)
  assert_eq(stats.total_records, 5)
  assert_eq(stats.snps, 1)
  assert_eq(stats.insertions, 1) // pos 400 G>GGT
  assert_eq(stats.deletions, 2) // pos 200 AT>A, pos 500 TT>T
  assert_eq(stats.multi_allelic, 1) // pos 300 C>A,T
  assert_eq(stats.pass_records, 4)
  assert_eq(stats.failed_records, 1)
}

test "vcf_concat" {
  let text_a =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	300	.	C	A	30	PASS	DP=5
    #|chr1	100	.	A	G	50	PASS	DP=10
  let text_b =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr2	50	.	G	T	20	PASS	DP=3
    #|chr1	200	.	T	C	40	PASS	DP=8
  let vcf_a = @src.parse_vcf(text_a)
  let vcf_b = @src.parse_vcf(text_b)
  let result = @src.vcf_concat([vcf_a, vcf_b])
  // 4 records, sorted
  assert_eq(result.records.length(), 4)
  assert_eq(result.records[0].pos, 100)
  assert_eq(result.records[3].chrom, "chr2")
}

// ============================================================
// bcftools annotate tests
// ============================================================

test "vcf_annotate_set_info" {
  let vcf = make_test_vcf()
  let annotated = @src.vcf_annotate_set_info(vcf, "AF", "0.01")
  for rec in annotated.records {
    match rec.info.get("AF") {
      Some(_) => ()
      None => fail("expected Some")
    }
    assert_eq(rec.info["AF"], "0.01")
  }
}

test "vcf_annotate_remove_info" {
  let vcf = make_test_vcf()
  // Record 3 has DP=20;AF=0.5
  match vcf.records[2].info.get("AF") {
      Some(_) => ()
      None => fail("expected Some")
    }
  let cleaned = @src.vcf_annotate_remove_info(vcf, "AF")
  for rec in cleaned.records {
    match rec.info.get("AF") {
      None => ()
      Some(_) => fail("expected None")
    }
  }
}

test "vcf_annotate_rename_info" {
  let vcf = make_test_vcf()
  let renamed = @src.vcf_annotate_rename_info(vcf, "DP", "DEPTH")
  for rec in renamed.records {
    match rec.info.get("DP") {
      None => ()
      Some(_) => fail("expected None")
    }
    match rec.info.get("DEPTH") {
      Some(_) => ()
      None => fail("expected Some")
    }
  }
}

test "vcf_annotate_set_id" {
  let vcf = make_test_vcf()
  let annotated = @src.vcf_annotate_set_id(vcf, "rs12345")
  for rec in annotated.records {
    assert_eq(rec.id.length(), 1)
    assert_eq(rec.id[0], "rs12345")
  }
  // Clear with "."
  let cleared = @src.vcf_annotate_set_id(vcf, ".")
  for rec in cleared.records {
    assert_eq(rec.id.length(), 0)
  }
}

test "vcf_annotate_set_filter" {
  let vcf = make_test_vcf()
  let annotated = @src.vcf_annotate_set_filter(vcf, "LowQual")
  for rec in annotated.records {
    assert_eq(rec.filter.length(), 1)
    assert_eq(rec.filter[0], "LowQual")
  }
  // Clear with "."
  let cleared = @src.vcf_annotate_set_filter(vcf, ".")
  for rec in cleared.records {
    assert_eq(rec.filter.length(), 0)
  }
}

test "vcf_annotate_set_qual" {
  let vcf = make_test_vcf()
  let annotated = @src.vcf_annotate_set_qual(vcf, 99.0)
  for rec in annotated.records {
    match rec.qual {
      Some(q) => assert_eq(q, 99.0)
      None => fail("qual should be set")
    }
  }
}

test "vcf_annotate_info_if_filter" {
  let vcf = make_test_vcf()
  // Record 2 has LowQual filter
  let annotated = @src.vcf_annotate_info_if_filter(vcf, "LowQual", "CONSEQUENCE", "low")
  // Record 2 should have the new info
  assert_eq(annotated.records[1].info["CONSEQUENCE"], "low")
  // Other records should NOT have it
  match annotated.records[0].info.get("CONSEQUENCE") {
      None => ()
      Some(_) => fail("expected None")
    }
}

test "vcf_annotate_from" {
  let target = make_test_vcf()
  // Source has matching record at chr1:100 with extra INFO
  let source_text =
    #|##fileformat=VCFv4.2
    #|#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    #|chr1	100	.	A	G	50	PASS	DP=10;GENE=BRCA1;IMPACT=HIGH
    #|chr1	200	.	AT	A	30	PASS	DP=5;GENE=TP53
  let source = @src.parse_vcf(source_text)
  let annotated = @src.vcf_annotate_from(target, source)
  // Record at chr1:100 should now have GENE and IMPACT
  assert_eq(annotated.records[0].info["GENE"], "BRCA1")
  assert_eq(annotated.records[0].info["IMPACT"], "HIGH")
  // Record at chr1:200 should have GENE
  assert_eq(annotated.records[1].info["GENE"], "TP53")
  // Record at chr1:300 has no match in source → no new annotations
  match annotated.records[2].info.get("GENE") {
      None => ()
      Some(_) => fail("expected None")
    }
}
