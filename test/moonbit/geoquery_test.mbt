///|
test "geoquery_parse_geo_series_matrix" {
  let content = "!Series_geo_accession = GSE12345\n" +
    "!Platform_geo_accession = GPL96\n" +
    "!Sample_geo_accession = GSM1\n" +
    "!Sample_geo_accession = GSM2\n" +
    "ID_REF\tGSM1\tGSM2\n" +
    "gene1\t10.5\t20.3\n" +
    "gene2\t15.2\t25.8\n"

  let geo = @src.parse_geo_series_matrix(content)

  assert_eq(geo.series_id, "GSE12345")
  assert_eq(geo.platform_id, "GPL96")
  assert_eq(geo.sample_ids.length(), 2)
  assert_eq(geo.gene_names.length(), 2)
}

///|
test "geoquery_parse_soft_format" {
  let content = "!Series_geo_accession = GSE12345\n" +
    "!Series_title = Test Series\n" +
    "!Series_summary = A test series\n" +
    "^SAMPLE=GSM1\n" +
    "!Sample_platform_id = GPL96\n" +
    "^SAMPLE=GSM2\n"

  let series = @src.parse_soft_format(content)

  assert_eq(series.gse_id, "GSE12345")
  assert_eq(series.title, "Test Series")
  assert_eq(series.sample_ids.length(), 2)
}

///|
test "geoquery_geo_data_to_expression_set" {
  let data = [[10.0, 20.0], [15.0, 25.0]]
  let geo = @src.geo_data_new(
    "GSE12345",
    ["GSM1", "GSM2"],
    "GPL96",
    data,
    ["gene1", "gene2"],
    Map([]),
  )

  let eset = @src.geo_data_to_expression_set(geo)

  assert_eq(@src.expression_set_nrow(eset), 2)
  assert_eq(@src.expression_set_ncol(eset), 2)
}

///|
test "geoquery_geo_data_log2_transform" {
  let data = [[8.0, 16.0], [32.0, 64.0]]
  let geo = @src.geo_data_new(
    "GSE12345",
    ["GSM1", "GSM2"],
    "GPL96",
    data,
    ["gene1", "gene2"],
    Map([]),
  )

  let transformed = @src.geo_data_log2_transform(geo)

  assert_eq(transformed.data[0][0], 3.0)
  assert_eq(transformed.data[0][1], 4.0)
}

///|
test "geoquery_geo_data_summary" {
  let data = [[10.0, 20.0]]
  let geo = @src.geo_data_new(
    "GSE12345",
    ["GSM1", "GSM2"],
    "GPL96",
    data,
    ["gene1"],
    Map([]),
  )

  let summary = @src.geo_data_summary(geo)

  assert_true(summary.contains("GEOData Summary"))
  assert_true(summary.contains("Series ID: GSE12345"))
}
