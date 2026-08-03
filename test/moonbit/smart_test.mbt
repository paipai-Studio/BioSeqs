///|
/// Tests for the SMART domain database format parser (smart.mbt).

// ============================================================================
// SmartResult Construction Tests
// ============================================================================

test "smart_result_construction" {
  let result = @src.SmartResult::new("test_protein")
  assert_eq(result.sequence_id, "test_protein")
  assert_eq(result.sequence_length, 0)
  assert_eq(result.domains.length(), 0)
  assert_eq(result.pfam_domains.length(), 0)
  assert_eq(result.signal_peptides.length(), 0)
  assert_eq(result.transmembrane.length(), 0)
}

// ============================================================================
// SmartDomain Tests
// ============================================================================

test "smart_domain_construction" {
  let domain = @src.SmartDomain::new("SM00001", "ABC_membrane", 1, 100, evalue=1.5e-20)
  assert_eq(domain.domain_id, "SM00001")
  assert_eq(domain.domain_name, "ABC_membrane")
  assert_eq(domain.start, 1)
  assert_eq(domain.end, 100)
  assert_true(domain.evalue > 0.0)
}

// ============================================================================
// SMART Parsing Tests
// ============================================================================

test "smart_parse_empty" {
  let result = @src.smart_parse("")
  assert_eq(result.sequence_id, "")
  assert_eq(result.sequence_length, 0)
  assert_eq(result.domains.length(), 0)
}

test "smart_parse_sample" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  assert_eq(result.sequence_id, "Protein123")
  assert_true(result.domains.length() >= 3)
}

test "smart_parse_domain_type_classification" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  assert_true(result.pfam_domains.length() >= 1)
  assert_true(result.signal_peptides.length() >= 1)
  assert_true(result.transmembrane.length() >= 1)
}

// ============================================================================
// SMART Query Tests
// ============================================================================

test "smart_find_domains_by_name" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  let matches = @src.smart_find_domains(result, "ABC")
  assert_true(matches.length() >= 1)
}

test "smart_filter_by_evalue" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  let significant = @src.smart_filter_evalue(result, 1.0e-10)
  assert_true(significant.length() >= 1)
}

test "smart_best_domain" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  match @src.smart_best_domain(result) {
    Some(best) => assert_true(best.evalue >= 0.0)
    None => assert_true(false)
  }
}

// ============================================================================
// SMART Summary Tests
// ============================================================================

test "smart_summary" {
  let result = @src.SmartResult::new("test_protein")
  let summary = @src.smart_summary(result)
  assert_true(summary.contains("test_protein"))
}

test "smart_total_domains" {
  let content = @src.smart_sample()
  let result = @src.smart_parse(content)
  let total = @src.smart_total_domains(result)
  assert_true(total >= 3)
}
