///|
/// Tests for Bio.Geo SOFT format parser.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "geo_entity_type_to_string" {
  assert_eq(@src.GeoEntityType::to_string(@src.geo_platform()), "Platform")
  assert_eq(@src.GeoEntityType::to_string(@src.geo_sample()), "Sample")
  assert_eq(@src.GeoEntityType::to_string(@src.geo_series()), "Series")
  assert_eq(@src.GeoEntityType::to_string(@src.geo_dataset()), "Dataset")
  assert_eq(@src.GeoEntityType::to_string(@src.geo_unknown()), "Unknown")
}

test "geo_attribute_creation" {
  let a = @src.GeoAttribute::new("title", "Test Sample")
  assert_eq(a.name(), "title")
  assert_eq(a.value(), "Test Sample")
}

test "geo_column_def_creation" {
  let c = @src.GeoColumnDef::new("ID", "Probe Set ID")
  assert_eq(c.name(), "ID")
  assert_eq(c.description(), "Probe Set ID")
}

test "geo_data_table_creation" {
  let cols = [
    @src.GeoColumnDef::new("ID", "Probe ID"),
    @src.GeoColumnDef::new("VALUE", "Expression"),
  ]
  let rows = [
    ["1007_s_at", "8.5"],
    ["117_at", "7.2"],
  ]
  let t = @src.GeoDataTable::new(cols, rows)
  assert_eq(t.n_columns(), 2)
  assert_eq(t.n_rows(), 2)
}

test "geo_data_table_column_index" {
  let cols = [
    @src.GeoColumnDef::new("ID", "Probe ID"),
    @src.GeoColumnDef::new("VALUE", "Expression"),
  ]
  let t = @src.GeoDataTable::new(cols, [])
  assert_eq(t.column_index("ID"), 0)
  assert_eq(t.column_index("VALUE"), 1)
  assert_eq(t.column_index("MISSING"), -1)
}

test "geo_data_table_get_column" {
  let cols = [
    @src.GeoColumnDef::new("ID", "Probe ID"),
    @src.GeoColumnDef::new("VALUE", "Expression"),
  ]
  let rows = [
    ["1007_s_at", "8.5"],
    ["117_at", "7.2"],
    ["1294_at", "9.1"],
  ]
  let t = @src.GeoDataTable::new(cols, rows)
  let ids = t.get_column("ID")
  assert_eq(ids.length(), 3)
  assert_eq(ids[0], "1007_s_at")
  assert_eq(ids[2], "1294_at")
  let vals = t.get_column("VALUE")
  assert_eq(vals[0], "8.5")
  assert_eq(vals[1], "7.2")
}

test "geo_data_table_get_column_missing" {
  let cols = [@src.GeoColumnDef::new("ID", "Probe ID")]
  let t = @src.GeoDataTable::new(cols, [])
  let result = t.get_column("MISSING")
  assert_eq(result.length(), 0)
}

test "geo_record_creation" {
  let attrs = [@src.GeoAttribute::new("title", "Test")]
  let cols = [@src.GeoColumnDef::new("ID", "Probe ID")]
  let table = @src.GeoDataTable::new(cols, [["probe1"]])
  let r = @src.GeoRecord::new(@src.geo_sample(), "GSM1", attrs, table)
  assert_eq(r.accession(), "GSM1")
  assert_eq(r.attributes().length(), 1)
  assert_eq(r.title(), "Test")
}

test "geo_record_get_attribute" {
  let attrs = [
    @src.GeoAttribute::new("title", "Test Sample"),
    @src.GeoAttribute::new("source", "brain"),
  ]
  let table = @src.GeoDataTable::new([], [])
  let r = @src.GeoRecord::new(@src.geo_sample(), "GSM1", attrs, table)
  assert_true(r.get_attribute("title").unwrap() == "Test Sample")
  assert_true(r.get_attribute("source").unwrap() == "brain")
  assert_true(r.get_attribute("missing") is None)
}

test "geo_file_creation" {
  let table = @src.GeoDataTable::new([], [])
  let r1 = @src.GeoRecord::new(@src.geo_platform(), "GPL1", [], table)
  let r2 = @src.GeoRecord::new(@src.geo_sample(), "GSM1", [], table)
  let f = @src.GeoFile::new([r1, r2])
  assert_eq(f.n_records(), 2)
}

test "geo_file_get_by_type" {
  let table = @src.GeoDataTable::new([], [])
  let r1 = @src.GeoRecord::new(@src.geo_platform(), "GPL1", [], table)
  let r2 = @src.GeoRecord::new(@src.geo_sample(), "GSM1", [], table)
  let r3 = @src.GeoRecord::new(@src.geo_sample(), "GSM2", [], table)
  let f = @src.GeoFile::new([r1, r2, r3])
  assert_eq(f.get_by_type(@src.geo_platform()).length(), 1)
  assert_eq(f.get_by_type(@src.geo_sample()).length(), 2)
  assert_eq(f.get_by_type(@src.geo_series()).length(), 0)
}

test "geo_file_get_by_accession" {
  let table = @src.GeoDataTable::new([], [])
  let r1 = @src.GeoRecord::new(@src.geo_platform(), "GPL1", [], table)
  let r2 = @src.GeoRecord::new(@src.geo_sample(), "GSM1", [], table)
  let f = @src.GeoFile::new([r1, r2])
  assert_true(f.get_by_accession("GPL1").unwrap().accession() == "GPL1")
  assert_true(f.get_by_accession("GSM1").unwrap().accession() == "GSM1")
  assert_true(f.get_by_accession("MISSING") is None)
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "geo_sample_text_has_platform" {
  let text = @src.geo_sample_text()
  assert_true(text.contains("^PLATFORM"))
  assert_true(text.contains("GPL1"))
}

test "geo_sample_text_has_samples" {
  let text = @src.geo_sample_text()
  assert_true(text.contains("^SAMPLE"))
  assert_true(text.contains("GSM1"))
  assert_true(text.contains("GSM2"))
}

test "geo_sample_text_has_attributes" {
  let text = @src.geo_sample_text()
  assert_true(text.contains("!Platform_title"))
  assert_true(text.contains("!Sample_title"))
  assert_true(text.contains("!Sample_source"))
}

test "geo_sample_text_has_columns" {
  let text = @src.geo_sample_text()
  assert_true(text.contains("#ID"))
  assert_true(text.contains("#VALUE"))
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

test "geo_parse_record_count" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  assert_eq(f.n_records(), 3)
}

test "geo_parse_platform" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let platforms = f.get_by_type(@src.geo_platform())
  assert_eq(platforms.length(), 1)
  let p = platforms[0]
  assert_eq(p.accession(), "GPL1")
  assert_true(p.title().contains("Test Affymetrix"))
}

test "geo_parse_samples" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let samples = f.get_by_type(@src.geo_sample())
  assert_eq(samples.length(), 2)
  assert_eq(samples[0].accession(), "GSM1")
  assert_eq(samples[1].accession(), "GSM2")
}

test "geo_parse_sample_attributes" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let gsm1 = f.get_by_accession("GSM1").unwrap()
  assert_true(gsm1.get_attribute("Sample_title").unwrap() == "Control Sample 1")
  assert_true(gsm1.get_attribute("Sample_source").unwrap() == "brain tissue")
  assert_true(gsm1.get_attribute("Sample_characteristics").unwrap() == "normal")
  // title() convenience method should also work
  assert_true(gsm1.title() == "Control Sample 1")
}

test "geo_parse_data_table" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let p = f.get_by_accession("GPL1").unwrap()
  let table = p.data_table()
  assert_eq(table.n_columns(), 2)
  assert_eq(table.n_rows(), 3)
}

test "geo_parse_data_table_column_values" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let p = f.get_by_accession("GPL1").unwrap()
  let ids = p.data_table().get_column("ID")
  assert_eq(ids.length(), 3)
  assert_eq(ids[0], "1007_s_at")
  assert_eq(ids[2], "1294_at")
  let vals = p.data_table().get_column("VALUE")
  assert_eq(vals[0], "8.5")
  assert_eq(vals[2], "9.1")
}

test "geo_parse_get_by_accession" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let gsm2 = f.get_by_accession("GSM2").unwrap()
  assert_eq(gsm2.accession(), "GSM2")
  assert_true(gsm2.title().contains("Treatment"))
}

test "geo_parse_empty_input" {
  let f = @src.geo_parse("")
  assert_eq(f.n_records(), 0)
}

test "geo_parse_minimal_record" {
  let text = "^SAMPLE = GSM1\n!Sample_title = Test\n"
  let f = @src.geo_parse(text)
  assert_eq(f.n_records(), 1)
  let r = f.get_by_accession("GSM1").unwrap()
  assert_eq(r.title(), "Test")
  assert_eq(r.data_table().n_rows(), 0)
}

test "geo_parse_record_without_data_table" {
  let text = "^SERIES = GSE1\n!Series_title = Test Series\n!Series_summary = A test\n"
  let f = @src.geo_parse(text)
  assert_eq(f.n_records(), 1)
  let r = f.get_by_accession("GSE1").unwrap()
  assert_eq(r.attributes().length(), 2)
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

test "geo_attribute_to_string" {
  let a = @src.GeoAttribute::new("title", "Test")
  let s = a.to_string()
  assert_true(s.contains("title"))
  assert_true(s.contains("Test"))
}

test "geo_data_table_to_string" {
  let cols = [
    @src.GeoColumnDef::new("ID", "Probe"),
    @src.GeoColumnDef::new("VALUE", "Expression"),
  ]
  let rows = [["p1", "1.0"], ["p2", "2.0"]]
  let t = @src.GeoDataTable::new(cols, rows)
  let s = t.to_string()
  assert_true(s.contains("DataTable"))
  assert_true(s.contains("ID"))
  assert_true(s.contains("VALUE"))
}

test "geo_record_to_string" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let r = f.get_by_accession("GSM1").unwrap()
  let s = r.to_string()
  assert_true(s.contains("GeoRecord"))
  assert_true(s.contains("GSM1"))
  assert_true(s.contains("Sample"))
}

test "geo_file_to_string" {
  let text = @src.geo_sample_text()
  let f = @src.geo_parse(text)
  let s = f.to_string()
  assert_true(s.contains("GEO SOFT"))
  assert_true(s.contains("3 records"))
}
