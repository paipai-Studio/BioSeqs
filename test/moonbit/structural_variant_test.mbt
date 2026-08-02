///|
/// Tests for the StructuralVariantAnnotation module.
///
/// Covers SV types, breakend parsing, VCF line parsing, partner pairing,
// gene overlap annotation, filtering, summary, and sample data validation.

// ===========================================================================
// 1. SvType enum
// ===========================================================================

///|
test "sv_type_labels" {
  assert_eq(@src.sv_type_label(@src.sv_type_del()), "DEL")
  assert_eq(@src.sv_type_label(@src.sv_type_dup()), "DUP")
  assert_eq(@src.sv_type_label(@src.sv_type_inv()), "INV")
  assert_eq(@src.sv_type_label(@src.sv_type_ins()), "INS")
  assert_eq(@src.sv_type_label(@src.sv_type_tra()), "TRA")
  assert_eq(@src.sv_type_label(@src.sv_type_bnd()), "BND")
}

///|
test "sv_parse_type_known" {
  assert_true(@src.sv_type_label(@src.sv_parse_type("DEL")) == "DEL")
  assert_true(@src.sv_type_label(@src.sv_parse_type("DUP")) == "DUP")
  assert_true(@src.sv_type_label(@src.sv_parse_type("INV")) == "INV")
  assert_true(@src.sv_type_label(@src.sv_parse_type("INS")) == "INS")
  assert_true(@src.sv_type_label(@src.sv_parse_type("TRA")) == "TRA")
  assert_true(@src.sv_type_label(@src.sv_parse_type("BND")) == "BND")
}

///|
test "sv_parse_type_case_insensitive" {
  assert_true(@src.sv_type_label(@src.sv_parse_type("del")) == "DEL")
  assert_true(@src.sv_type_label(@src.sv_parse_type("Dup")) == "DUP")
}

///|
test "sv_parse_type_unknown" {
  // Unknown types default to BND.
  assert_true(@src.sv_type_label(@src.sv_parse_type("UNKNOWN")) == "BND")
}

// ===========================================================================
// 2. SvBreakend
// ===========================================================================

///|
test "sv_breakend_creation" {
  let b = @src.SvBreakend::new("chr1", 1000, "+", "chr2", 2000, "-")
  assert_eq(b.chrom(), "chr1")
  assert_eq(b.pos(), 1000)
  assert_eq(b.strand(), "+")
  assert_eq(b.mate_chrom(), "chr2")
  assert_eq(b.mate_pos(), 2000)
  assert_eq(b.mate_strand(), "-")
}

///|
test "sv_breakend_has_mate" {
  let b = @src.SvBreakend::new("chr1", 1000, "+", "chr2", 2000, "-")
  assert_true(b.has_mate())
  let b2 = @src.SvBreakend::new("chr1", 1000, "+", "", 0, "")
  assert_false(b2.has_mate())
}

///|
test "sv_breakend_is_inter_chromosomal" {
  let b1 = @src.SvBreakend::new("chr1", 1000, "+", "chr2", 2000, "-")
  assert_true(b1.is_inter_chromosomal())
  let b2 = @src.SvBreakend::new("chr1", 1000, "+", "chr1", 2000, "-")
  assert_false(b2.is_inter_chromosomal())
}

// ===========================================================================
// 3. SvRecord
// ===========================================================================

///|
test "sv_record_creation" {
  let rec = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 500.0, "PASS",
  )
  assert_eq(rec.id(), "sv1")
  assert_eq(rec.chrom(), "chr1")
  assert_eq(rec.pos(), 1000)
  assert_eq(rec.end(), 2000)
  assert_eq(rec.svlen(), -1000)
  assert_true(rec.quality() == 500.0)
  assert_eq(rec.filter(), "PASS")
}

///|
test "sv_record_is_bnd" {
  let bnd = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_bnd(), 1000, 0, None, 500.0, "PASS",
  )
  assert_true(bnd.is_bnd())
  let del = @src.SvRecord::new(
    "sv2", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 500.0, "PASS",
  )
  assert_false(del.is_bnd())
}

///|
test "sv_record_size" {
  let del = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 0.0, "PASS",
  )
  assert_eq(del.size(), 1000)
}

///|
test "sv_record_is_inter_chromosomal" {
  let be = @src.SvBreakend::new("chr1", 1000, "+", "chr2", 2000, "-")
  let bnd = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_bnd(), 1000, 0, Some(be), 0.0, "PASS",
  )
  assert_true(bnd.is_inter_chromosomal())
  let del = @src.SvRecord::new(
    "sv2", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 0.0, "PASS",
  )
  assert_false(del.is_inter_chromosomal())
}

// ===========================================================================
// 4. BND ALT field parsing
// ===========================================================================

///|
test "sv_parse_bnd_alt_right_forward" {
  // t[p[ : extends right, forward strand
  let alt = "A[chr2:2000["
  let be = @src.sv_parse_bnd_alt(alt, "chr1", 1000)
  assert_true(be is Some(_))
  match be {
    Some(b) => {
      assert_eq(b.chrom(), "chr1")
      assert_eq(b.pos(), 1000)
      assert_eq(b.mate_chrom(), "chr2")
      assert_eq(b.mate_pos(), 2000)
    }
    None => assert_true(false, msg="expected breakend")
  }
}

///|
test "sv_parse_bnd_alt_right_reverse" {
  // t]p] : extends right, reverse strand
  let alt = "A]chr2:3000]"
  let be = @src.sv_parse_bnd_alt(alt, "chr1", 1000)
  assert_true(be is Some(_))
  match be {
    Some(b) => {
      assert_eq(b.mate_chrom(), "chr2")
      assert_eq(b.mate_pos(), 3000)
    }
    None => assert_true(false, msg="expected breakend")
  }
}

///|
test "sv_parse_bnd_alt_left_reverse" {
  // ]p]t : extends left, reverse strand
  let alt = "]chr2:4000]A"
  let be = @src.sv_parse_bnd_alt(alt, "chr1", 1000)
  assert_true(be is Some(_))
  match be {
    Some(b) => {
      assert_eq(b.mate_chrom(), "chr2")
      assert_eq(b.mate_pos(), 4000)
    }
    None => assert_true(false, msg="expected breakend")
  }
}

///|
test "sv_parse_bnd_alt_left_forward" {
  // [p[t : extends left, forward strand
  let alt = "[chr2:5000[A"
  let be = @src.sv_parse_bnd_alt(alt, "chr1", 1000)
  assert_true(be is Some(_))
  match be {
    Some(b) => {
      assert_eq(b.mate_chrom(), "chr2")
      assert_eq(b.mate_pos(), 5000)
    }
    None => assert_true(false, msg="expected breakend")
  }
}

///|
test "sv_parse_bnd_alt_invalid" {
  let be = @src.sv_parse_bnd_alt("INVALID", "chr1", 1000)
  assert_true(be is None)
}

// ===========================================================================
// 5. VCF parsing
// ===========================================================================

///|
test "sv_parse_vcf_line_del" {
  let line = "chr1\t1000\tsv_del_1\tA\t<DEL>\t500.0\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000"
  let rec = @src.sv_parse_vcf_line(line)
  assert_true(rec is Some(_))
  match rec {
    Some(r) => {
      assert_eq(r.id(), "sv_del_1")
      assert_eq(r.chrom(), "chr1")
      assert_eq(r.pos(), 1000)
      assert_eq(r.end(), 2000)
      assert_eq(r.svlen(), -1000)
      assert_true(r.quality() == 500.0)
      assert_eq(r.filter(), "PASS")
      assert_true(@src.sv_type_label(r.sv_type()) == "DEL")
    }
    None => assert_true(false, msg="expected record")
  }
}

///|
test "sv_parse_vcf_line_dup" {
  let line = "chr1\t3000\tsv_dup_1\tA\t<DUP>\t450.0\tPASS\tSVTYPE=DUP;END=5000;SVLEN=2000"
  let rec = @src.sv_parse_vcf_line(line)
  match rec {
    Some(r) => {
      assert_eq(r.id(), "sv_dup_1")
      assert_eq(r.end(), 5000)
      assert_eq(r.svlen(), 2000)
      assert_true(@src.sv_type_label(r.sv_type()) == "DUP")
    }
    None => assert_true(false, msg="expected record")
  }
}

///|
test "sv_parse_vcf_line_bnd" {
  let line = "chr1\t6000\tsv_bnd_1\tA\tA]chr2:7000]\t600.0\tPASS\tSVTYPE=BND"
  let rec = @src.sv_parse_vcf_line(line)
  match rec {
    Some(r) => {
      assert_eq(r.id(), "sv_bnd_1")
      assert_true(r.is_bnd())
      match r.breakend() {
        Some(b) => {
          assert_eq(b.mate_chrom(), "chr2")
          assert_eq(b.mate_pos(), 7000)
        }
        None => assert_true(false, msg="expected breakend")
      }
    }
    None => assert_true(false, msg="expected record")
  }
}

///|
test "sv_parse_vcf_skips_header" {
  let line = "##fileformat=VCFv4.2"
  assert_true(@src.sv_parse_vcf_line(line) is None)
  let line2 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  assert_true(@src.sv_parse_vcf_line(line2) is None)
}

///|
test "sv_parse_vcf_multiple_records" {
  let content = @src.sv_sample_vcf()
  let records = @src.sv_parse_vcf(content)
  // Should have 6 SV records (DEL, DUP, INV, INS, 2x BND).
  assert_eq(records.length(), 6)
}

///|
test "sv_parse_vcf_empty" {
  let records = @src.sv_parse_vcf("")
  assert_eq(records.length(), 0)
}

// ===========================================================================
// 6. Partner pairing
// ===========================================================================

///|
test "sv_find_partners_sample" {
  let records = @src.sv_sample_records()
  let pairs = @src.sv_find_partners(records)
  // The two BND records (sv_bnd_1 and sv_bnd_2) should be partners.
  assert_eq(pairs.length(), 1)
  let (id_a, id_b) = pairs[0]
  assert_true(
    (id_a == "sv_bnd_1" && id_b == "sv_bnd_2") ||
      (id_a == "sv_bnd_2" && id_b == "sv_bnd_1"),
  )
}

///|
test "sv_are_partners_mutual" {
  let records = @src.sv_sample_records()
  // Find the two BND records.
  let bnd1 = records[4]
  let bnd2 = records[5]
  assert_true(@src.sv_are_partners(bnd1, bnd2))
  assert_true(@src.sv_are_partners(bnd2, bnd1))
}

///|
test "sv_are_partners_non_bnd" {
  let del = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 0.0, "PASS",
  )
  let dup = @src.SvRecord::new(
    "sv2", "chr1", 3000, @src.sv_type_dup(), 5000, 2000, None, 0.0, "PASS",
  )
  assert_false(@src.sv_are_partners(del, dup))
}

// ===========================================================================
// 7. Gene region annotation
// ===========================================================================

///|
test "sv_gene_region_creation" {
  let g = @src.SvGeneRegion::new("BRCA1", "chr1", 800, 2500, "+")
  assert_eq(g.gene(), "BRCA1")
  assert_eq(g.chrom(), "chr1")
  assert_eq(g.start(), 800)
  assert_eq(g.end(), 2500)
  assert_eq(g.strand(), "+")
}

///|
test "sv_annotate_overlapping_genes" {
  let records = @src.sv_sample_records()
  let genes = @src.sv_sample_genes()
  let annotations = @src.sv_annotate(records, genes)
  assert_eq(annotations.length(), 6)
  // The DEL at chr1:1000-2000 should overlap BRCA1 (chr1:800-2500).
  let del_ann = annotations[0]
  assert_true(del_ann.genes().contains("BRCA1"))
}

///|
test "sv_annotate_has_gene_overlap" {
  let records = @src.sv_sample_records()
  let genes = @src.sv_sample_genes()
  let annotations = @src.sv_annotate(records, genes)
  // At least one SV should overlap a gene.
  let mut any_overlap = false
  for ann in annotations {
    if ann.has_gene_overlap() {
      any_overlap = true
    }
  }
  assert_true(any_overlap)
}

///|
test "sv_annotation_summary" {
  let records = @src.sv_sample_records()
  let genes = @src.sv_sample_genes()
  let annotations = @src.sv_annotate(records, genes)
  let summary = @src.sv_annotation_summary(annotations[0])
  assert_true(summary.contains("sv_del_1"))
  assert_true(summary.contains("DEL"))
}

// ===========================================================================
// 8. Filtering
// ===========================================================================

///|
test "sv_filter_by_quality" {
  let records = @src.sv_sample_records()
  let filtered = @src.sv_filter_by_quality(records, 400.0)
  // Records with QUAL >= 400: sv_del_1(500), sv_dup_1(450), sv_bnd_1(600), sv_bnd_2(600) = 4
  assert_eq(filtered.length(), 4)
}

///|
test "sv_filter_by_size" {
  let records = @src.sv_sample_records()
  let filtered = @src.sv_filter_by_size(records, 500)
  // Records with size >= 500: DEL(1000), DUP(2000), INV(1000) = 3
  // INS has SVLEN=500 (>= 500), so 4.
  assert_eq(filtered.length(), 4)
}

///|
test "sv_filter_by_type" {
  let records = @src.sv_sample_records()
  let dels = @src.sv_filter_by_type(records, @src.sv_type_del())
  assert_eq(dels.length(), 1)
  let bnds = @src.sv_filter_by_type(records, @src.sv_type_bnd())
  assert_eq(bnds.length(), 2)
}

// ===========================================================================
// 9. Summary and counting
// ===========================================================================

///|
test "sv_count_by_type" {
  let records = @src.sv_sample_records()
  let counts = @src.sv_count_by_type(records)
  assert_eq(counts["DEL"], 1)
  assert_eq(counts["DUP"], 1)
  assert_eq(counts["INV"], 1)
  assert_eq(counts["INS"], 1)
  assert_eq(counts["BND"], 2)
}

///|
test "sv_record_summary" {
  let rec = @src.sv_sample_records()[0]
  let summary = @src.sv_record_summary(rec)
  assert_true(summary.contains("sv_del_1"))
  assert_true(summary.contains("DEL"))
  assert_true(summary.contains("chr1"))
  assert_true(summary.contains("1000"))
}

// ===========================================================================
// 10. Sample data validation
// ===========================================================================

///|
test "sv_sample_vcf_has_header" {
  let vcf = @src.sv_sample_vcf()
  assert_true(vcf.contains("##fileformat=VCFv4.2"))
  assert_true(vcf.contains("#CHROM"))
}

///|
test "sv_sample_records_count" {
  let records = @src.sv_sample_records()
  assert_eq(records.length(), 6)
}

///|
test "sv_sample_genes_count" {
  let genes = @src.sv_sample_genes()
  assert_eq(genes.length(), 5)
  assert_eq(genes[0].gene(), "BRCA1")
  assert_eq(genes[4].gene(), "KRAS")
}

///|
test "sv_sample_records_types" {
  let records = @src.sv_sample_records()
  let types = @src.sv_count_by_type(records)
  assert_eq(types["DEL"], 1)
  assert_eq(types["DUP"], 1)
  assert_eq(types["INV"], 1)
  assert_eq(types["INS"], 1)
  assert_eq(types["BND"], 2)
}

// ===========================================================================
// 11. Edge cases
// ===========================================================================

///|
test "sv_parse_vcf_line_missing_svtype" {
  let line = "chr1\t1000\tsv1\tA\t<DEL>\t500.0\tPASS\tEND=2000"
  let rec = @src.sv_parse_vcf_line(line)
  // Without SVTYPE, defaults to BND.
  match rec {
    Some(r) => assert_true(@src.sv_type_label(r.sv_type()) == "BND")
    None => assert_true(false, msg="expected record")
  }
}

///|
test "sv_filter_empty_records" {
  let records : Array[@src.SvRecord] = []
  let filtered = @src.sv_filter_by_quality(records, 100.0)
  assert_eq(filtered.length(), 0)
}

///|
test "sv_find_partners_no_bnd" {
  let del = @src.SvRecord::new(
    "sv1", "chr1", 1000, @src.sv_type_del(), 2000, -1000, None, 0.0, "PASS",
  )
  let dup = @src.SvRecord::new(
    "sv2", "chr1", 3000, @src.sv_type_dup(), 5000, 2000, None, 0.0, "PASS",
  )
  let pairs = @src.sv_find_partners([del, dup])
  assert_eq(pairs.length(), 0)
}

///|
test "sv_annotate_no_genes" {
  let records = @src.sv_sample_records()
  let genes : Array[@src.SvGeneRegion] = []
  let annotations = @src.sv_annotate(records, genes)
  for ann in annotations {
    assert_false(ann.has_gene_overlap())
  }
}
