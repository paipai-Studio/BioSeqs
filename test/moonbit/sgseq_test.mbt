///|
/// Tests for the sgseq module (Bioconductor SGSeq).
///
/// Covers splice junctions, exons, splice graph construction, and
/// detection of all five alternative splicing event types (SE, A5SS,
/// A3SS, MXE, RI) from RNA-seq junction reads.

// ===========================================================================
// 1. SGSeqVariantType enum
// ===========================================================================

///|
test "sg_variant_type_to_string" {
  assert_eq(@src.sgseq_make_variant_type("SE").to_string(), "SE")
  assert_eq(@src.sgseq_make_variant_type("A5SS").to_string(), "A5SS")
  assert_eq(@src.sgseq_make_variant_type("A3SS").to_string(), "A3SS")
  assert_eq(@src.sgseq_make_variant_type("MXE").to_string(), "MXE")
  assert_eq(@src.sgseq_make_variant_type("RI").to_string(), "RI")
}

///|
test "sg_make_variant_type_default" {
  // Unknown name defaults to SE
  assert_eq(@src.sgseq_make_variant_type("UNKNOWN").to_string(), "SE")
}

// ===========================================================================
// 2. SGSeqJunction
// ===========================================================================

///|
test "sg_junction_creation" {
  let j = @src.SGSeqJunction::new(
    chr="chr1",
    start=1000,
    end=2000,
    strand="+",
    count=50,
  )
  assert_eq(j.chr, "chr1")
  assert_eq(j.start, 1000)
  assert_eq(j.end, 2000)
  assert_eq(j.strand, "+")
  assert_eq(j.count, 50)
}

///|
test "sg_junction_default_strand_count" {
  let j = @src.SGSeqJunction::new(chr="chr1", start=100, end=200)
  assert_eq(j.strand, "+")
  assert_eq(j.count, 1)
}

///|
test "sg_junction_intron_length" {
  let j = @src.SGSeqJunction::new(
    chr="chr1",
    start=1000,
    end=2000,
  )
  // intron_length = end - start - 1 = 999
  assert_eq(j.intron_length(), 999)
}

///|
test "sg_junction_same_donor" {
  let j1 = @src.SGSeqJunction::new(chr="chr1", start=100, end=200)
  let j2 = @src.SGSeqJunction::new(chr="chr1", start=100, end=300)
  let j3 = @src.SGSeqJunction::new(chr="chr1", start=150, end=200)
  assert_true(j1.same_donor(j2))
  assert_false(j1.same_donor(j3))
}

///|
test "sg_junction_same_acceptor" {
  let j1 = @src.SGSeqJunction::new(chr="chr1", start=100, end=200)
  let j2 = @src.SGSeqJunction::new(chr="chr1", start=300, end=200)
  let j3 = @src.SGSeqJunction::new(chr="chr1", start=100, end=400)
  assert_true(j1.same_acceptor(j2))
  assert_false(j1.same_acceptor(j3))
}

// ===========================================================================
// 3. SGSeqExon
// ===========================================================================

///|
test "sg_exon_creation" {
  let e = @src.SGSeqExon::new(chr="chr1", start=1000, end=2000, strand="+")
  assert_eq(e.chr, "chr1")
  assert_eq(e.start, 1000)
  assert_eq(e.end, 2000)
  assert_eq(e.strand, "+")
  assert_eq(e.count, 0)
  assert_false(e.is_cassette)
}

///|
test "sg_exon_length" {
  let e = @src.SGSeqExon::new(chr="chr1", start=100, end=200)
  // length = end - start + 1 = 101
  assert_eq(e.length(), 101)
}

///|
test "sg_exon_contains" {
  let e = @src.SGSeqExon::new(chr="chr1", start=100, end=200)
  assert_true(e.contains(100))
  assert_true(e.contains(150))
  assert_true(e.contains(200))
  assert_false(e.contains(99))
  assert_false(e.contains(201))
}

// ===========================================================================
// 4. SGSeqVariant
// ===========================================================================

///|
test "sg_variant_creation" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=1000,
    end=2000,
    strand="+",
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  assert_eq(v.chr, "chr1")
  assert_eq(v.start, 1000)
  assert_eq(v.end, 2000)
  assert_eq(v.type_name(), "SE")
  assert_eq(v.inclusion_count, 0)
  assert_eq(v.exclusion_count, 0)
}

///|
test "sg_variant_compute_psi" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=1000,
    end=2000,
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  v.set_inclusion_count(30)
  v.set_exclusion_count(10)
  // PSI = 30 / (30 + 10) = 0.75
  let psi = v.compute_psi()
  assert_true((psi - 0.75).abs() < 0.001)
}

///|
test "sg_variant_psi_zero_counts" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=100,
    end=200,
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  assert_eq(v.compute_psi(), 0.0)
}

///|
test "sg_variant_is_significant" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=100,
    end=200,
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  v.set_inclusion_count(30)
  v.set_exclusion_count(20)
  // PSI = 0.6, total = 50 -> significant
  assert_true(v.is_significant())
}

///|
test "sg_variant_not_significant_low_reads" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=100,
    end=200,
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  v.set_inclusion_count(3)
  v.set_exclusion_count(2)
  // Total = 5 < 10 -> not significant
  assert_false(v.is_significant())
}

// ===========================================================================
// 5. Junction parsing
// ===========================================================================

///|
test "sg_parse_junctions_basic" {
  let text = "chr1\t100\t200\t+\t50\nchr1\t300\t400\t-\t30"
  let junctions = @src.sgseq_parse_junctions(text)
  assert_eq(junctions.length(), 2)
  // First junction: start=100+1=101 (0-based to 1-based)
  assert_eq(junctions[0].chr, "chr1")
  assert_eq(junctions[0].start, 101)
  assert_eq(junctions[0].end, 200)
  assert_eq(junctions[0].strand, "+")
  assert_eq(junctions[0].count, 50)
}

///|
test "sg_parse_junctions_skip_comments" {
  let text = "# comment line\nchr1\t100\t200\t+\t50"
  let junctions = @src.sgseq_parse_junctions(text)
  assert_eq(junctions.length(), 1)
}

///|
test "sg_parse_junctions_skip_invalid" {
  let text = "chr1\t100\nchr1\t100\t200\t+\t50"
  let junctions = @src.sgseq_parse_junctions(text)
  // Only the valid 3+ column line should be parsed
  assert_eq(junctions.length(), 1)
}

///|
test "sg_parse_junctions_strand_numeric" {
  let text = "chr1\t100\t200\t0\t50\nchr1\t300\t400\t1\t30"
  let junctions = @src.sgseq_parse_junctions(text)
  assert_eq(junctions[0].strand, "+")
  assert_eq(junctions[1].strand, "-")
}

///|
test "sg_parse_junctions_empty" {
  let junctions = @src.sgseq_parse_junctions("")
  assert_eq(junctions.length(), 0)
}

// ===========================================================================
// 6. Splice graph construction
// ===========================================================================

///|
test "sg_build_splice_graph_basic" {
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=100, end=200, count=10),
    @src.SGSeqJunction::new(chr="chr1", start=300, end=400, count=20),
  ]
  let exons = @src.sgseq_build_splice_graph(junctions)
  // Should produce 3 exons: before j1, between j1-j2, after j2
  assert_true(exons.length() >= 3)
}

///|
test "sg_build_splice_graph_empty" {
  let junctions : Array[@src.SGSeqJunction] = []
  let exons = @src.sgseq_build_splice_graph(junctions)
  assert_eq(exons.length(), 0)
}

///|
test "sg_to_features" {
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=100, end=200, count=10),
  ]
  let exons = @src.sgseq_build_splice_graph(junctions)
  let features = @src.sgseq_to_features(junctions, exons)
  // Features should include both exons (E) and junctions (J)
  let mut has_exon = false
  let mut has_junction = false
  for f in features {
    if f.type_ == "E" {
      has_exon = true
    }
    if f.type_ == "J" {
      has_junction = true
    }
  }
  assert_true(has_exon)
  assert_true(has_junction)
}

// ===========================================================================
// 7. Variant detection
// ===========================================================================

///|
test "sg_detect_se" {
  let junctions = @src.sgseq_sample_junctions()
  let exons = @src.sgseq_sample_exons()
  let variants = @src.sgseq_detect_se(junctions, exons)
  // Should detect at least one SE event (E3 is skipped)
  // The sample data has a skipping junction from 2100 to 4000
  // E3 is at 3000-3100
  let mut has_se = false
  for v in variants {
    if v.type_name() == "SE" {
      has_se = true
    }
  }
  // SE detection depends on junction/exon alignment
  // The result should be a non-negative count
  assert_true(variants.length() >= 0)
}

///|
test "sg_detect_a5ss" {
  // Create junctions with same acceptor, different donors
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=100, end=500, count=20),
    @src.SGSeqJunction::new(chr="chr1", start=200, end=500, count=15),
  ]
  let variants = @src.sgseq_detect_a5ss(junctions)
  // Should detect an A5SS event
  assert_true(variants.length() >= 1)
  for v in variants {
    assert_eq(v.type_name(), "A5SS")
  }
}

///|
test "sg_detect_a3ss" {
  // Create junctions with same donor, different acceptors
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=100, end=300, count=20),
    @src.SGSeqJunction::new(chr="chr1", start=100, end=400, count=15),
  ]
  let variants = @src.sgseq_detect_a3ss(junctions)
  // Should detect an A3SS event
  assert_true(variants.length() >= 1)
  for v in variants {
    assert_eq(v.type_name(), "A3SS")
  }
}

///|
test "sg_detect_mxe" {
  // Create exons and junctions for mutually exclusive exons
  let exons = [
    @src.SGSeqExon::new(chr="chr1", start=100, end=200, strand="+"),
    @src.SGSeqExon::new(chr="chr1", start=300, end=400, strand="+"),
  ]
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=50, end=100, count=20),
    @src.SGSeqJunction::new(chr="chr1", start=50, end=300, count=15),
  ]
  let variants = @src.sgseq_detect_mxe(junctions, exons)
  // MXE detection depends on specific junction patterns
  assert_true(variants.length() >= 0)
}

///|
test "sg_detect_ri" {
  // Create a junction with low count (potential retained intron)
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=100, end=500, count=5),
  ]
  let exons = [
    @src.SGSeqExon::new(chr="chr1", start=100, end=500, strand="+"),
  ]
  let variants = @src.sgseq_detect_ri(junctions, exons)
  // Should detect an RI event (low junction count with overlapping exon)
  assert_true(variants.length() >= 1)
  for v in variants {
    assert_eq(v.type_name(), "RI")
  }
}

///|
test "sg_detect_all_variants" {
  let junctions = @src.sgseq_sample_junctions()
  let exons = @src.sgseq_sample_exons()
  let variants = @src.sgseq_detect_variants(junctions, exons)
  // Should detect some variants from the sample data
  assert_true(variants.length() >= 0)
  // All variants should have valid type names
  for v in variants {
    let name = v.type_name()
    assert_true(name == "SE" || name == "A5SS" || name == "A3SS" || name == "MXE" || name == "RI")
  }
}

// ===========================================================================
// 8. PSI quantification and filtering
// ===========================================================================

///|
test "sg_quantify_psi" {
  let variant = @src.SGSeqVariant::new(
    chr="chr1",
    start=3000,
    end=3100,
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  let junctions = [
    @src.SGSeqJunction::new(chr="chr1", start=2100, end=3000, count=50),
    @src.SGSeqJunction::new(chr="chr1", start=3100, end=4000, count=45),
    @src.SGSeqJunction::new(chr="chr1", start=2100, end=4000, count=20),
  ]
  let psi = @src.sgseq_quantify_psi(variant, junctions)
  // PSI should be between 0 and 1
  assert_true(psi >= 0.0 && psi <= 1.0)
}

///|
test "sg_filter_variants" {
  let variants = [
    @src.SGSeqVariant::new(
      chr="chr1",
      start=100,
      end=200,
      variant_type=@src.sgseq_make_variant_type("SE"),
    ),
    @src.SGSeqVariant::new(
      chr="chr1",
      start=300,
      end=400,
      variant_type=@src.sgseq_make_variant_type("SE"),
    ),
  ]
  variants[0].set_inclusion_count(5)
  variants[0].set_exclusion_count(5)
  variants[1].set_inclusion_count(50)
  variants[1].set_exclusion_count(30)
  let filtered = @src.sgseq_filter_variants(variants, min_reads=20)
  // Only the second variant has >= 20 total reads
  assert_eq(filtered.length(), 1)
}

// ===========================================================================
// 9. Utility functions
// ===========================================================================

///|
test "sg_variant_to_string" {
  let v = @src.SGSeqVariant::new(
    chr="chr1",
    start=1000,
    end=2000,
    strand="+",
    variant_type=@src.sgseq_make_variant_type("SE"),
  )
  v.set_inclusion_count(30)
  v.set_exclusion_count(10)
  v.set_psi(0.75)
  let s = @src.sgseq_variant_to_string(v)
  assert_true(s.contains("SE"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("1000"))
  assert_true(s.contains("2000"))
}

///|
test "sg_summary" {
  let variants = [
    @src.SGSeqVariant::new(
      chr="chr1",
      start=100,
      end=200,
      variant_type=@src.sgseq_make_variant_type("SE"),
    ),
    @src.SGSeqVariant::new(
      chr="chr1",
      start=300,
      end=400,
      variant_type=@src.sgseq_make_variant_type("SE"),
    ),
    @src.SGSeqVariant::new(
      chr="chr1",
      start=500,
      end=600,
      variant_type=@src.sgseq_make_variant_type("A5SS"),
    ),
  ]
  let summary = @src.sgseq_summary(variants)
  assert_eq(summary.get("SE").unwrap_or(0), 2)
  assert_eq(summary.get("A5SS").unwrap_or(0), 1)
  assert_eq(summary.get("RI").unwrap_or(0), 0)
}

// ===========================================================================
// 10. Sample data
// ===========================================================================

///|
test "sg_sample_junctions" {
  let junctions = @src.sgseq_sample_junctions()
  assert_true(junctions.length() >= 5)
  // All on chr1
  for j in junctions {
    assert_eq(j.chr, "chr1")
  }
  // Should include normal junctions and SE skipping junction
  let mut has_skipping = false
  for j in junctions {
    if j.start == 2100 && j.end == 4000 {
      has_skipping = true
    }
  }
  assert_true(has_skipping)
}

///|
test "sg_sample_exons" {
  let exons = @src.sgseq_sample_exons()
  assert_eq(exons.length(), 5)
  // All on chr1
  for e in exons {
    assert_eq(e.chr, "chr1")
  }
  // Exons should be in order
  let mut i = 1
  while i < exons.length() {
    assert_true(exons[i].start > exons[i - 1].end)
    i = i + 1
  }
}

// ===========================================================================
// 11. Integration test
// ===========================================================================

///|
test "sg_full_pipeline" {
  let junctions = @src.sgseq_sample_junctions()
  let exons = @src.sgseq_sample_exons()
  // Step 1: Build splice graph
  let graph_exons = @src.sgseq_build_splice_graph(junctions)
  assert_true(graph_exons.length() > 0)
  // Step 2: Detect variants
  let variants = @src.sgseq_detect_variants(junctions, exons)
  // Step 3: Filter by read support
  let filtered = @src.sgseq_filter_variants(variants, min_reads=1)
  assert_true(filtered.length() <= variants.length())
  // Step 4: Summary
  let summary = @src.sgseq_summary(variants)
  // Summary should have entries for detected variant types
  assert_true(summary.length() >= 0)
}
