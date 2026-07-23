///|
/// Tests for ExPASy module.

test "ExPASyRecord creation" {
  let record = @src.ExPASyRecord::new("P04637", "Swiss-Prot")
  assert_eq(record.id, "P04637")
  assert_eq(record.database, "Swiss-Prot")
}

test "ExPASyRecord add field" {
  let record = @src.ExPASyRecord::new("P04637", "Swiss-Prot")
  let record2 = record.add_field("description", "Tumor protein p53")
  
  let desc = record2.get_field("description")
  assert_true(desc is Some(_))
  assert_eq(desc.unwrap(), "Tumor protein p53")
}

test "ExPASyEntry creation" {
  let entry = @src.ExPASyEntry::new("P04637", "TP53")
  assert_eq(entry.accession, "P04637")
  assert_eq(entry.name, "TP53")
}

test "ExPASyEntry add keyword" {
  let entry = @src.ExPASyEntry::new("P04637", "TP53")
  let entry2 = entry.add_keyword("Tumor suppressor")
  
  assert_eq(entry2.keywords.length(), 1)
  assert_eq(entry2.keywords[0], "Tumor suppressor")
}

test "EnzymeEntry creation" {
  let enzyme = @src.EnzymeEntry::new("1.1.1.1", "Alcohol dehydrogenase")
  assert_eq(enzyme.ec_number, "1.1.1.1")
  assert_eq(enzyme.name, "Alcohol dehydrogenase")
}

test "enzyme_parse_ec" {
  let (class, subclass, subsubclass, serial) = @src.enzyme_parse_ec("1.1.1.1")
  assert_eq(class, "1")
  assert_eq(subclass, "1")
  assert_eq(subsubclass, "1")
  assert_eq(serial, "1")
}

test "expasy_get_prosite_ids" {
  let ids = @src.expasy_get_prosite_ids("PCNG")
  assert_true(ids.length() > 0)
}

test "expasy_get_swissprot_entry" {
  let entry = @src.expasy_get_swissprot_entry("P04637")
  assert_true(entry is Some(_))
}

test "expasy_get_enzyme" {
  let enzyme = @src.expasy_get_enzyme("1.1.1.1")
  assert_true(enzyme is Some(_))
}

test "expasy_analyze_protein" {
  let results = @src.expasy_analyze_protein("AVG")
  assert_true(results.contains("molecular_weight"))
  assert_true(results.contains("gravy"))
}

test "create_example_expasy_entry" {
  let entry = @src.create_example_expasy_entry()
  assert_eq(entry.accession, "P04637")
  assert_eq(entry.name, "TP53")
}

test "create_example_enzyme_entry" {
  let enzyme = @src.create_example_enzyme_entry()
  assert_eq(enzyme.ec_number, "1.1.1.1")
  assert_eq(enzyme.name, "Alcohol dehydrogenase")
}