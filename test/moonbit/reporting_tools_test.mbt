///|
/// Test file for ReportingTools module.

test "report_new" {
  let doc = @src.ReportDocument::new("Test Report")
  assert_eq(doc.get_title(), "Test Report")
  assert_eq(doc.get_author(), "BioSeqs")
  assert_eq(doc.get_n_sections(), 0)
}

test "report_set_author" {
  let doc = @src.ReportDocument::new("Test")
  let doc2 = doc.set_author("User")
  assert_eq(doc2.get_author(), "User")
}

test "report_add_text" {
  let doc = @src.ReportDocument::new("Test")
  let doc2 = doc.add_text("sec1", "Introduction", "This is a test report.")
  assert_eq(doc2.get_n_sections(), 1)
  let sec = doc2.get_section("sec1")
  assert_eq(sec.title, "Introduction")
}

test "report_add_table" {
  let doc = @src.ReportDocument::new("Test")
  let table = @src.ReportTable::new("table1", "Results", ["Gene", "log2FC", "p-value"], [["TP53", "BRCA1"], ["2.5", "-3.0"], ["0.001", "0.005"]], "Differential expression results")
  let doc2 = doc.add_table("tab1", table)
  assert_eq(doc2.get_n_tables(), 1)
  assert_eq(doc2.get_n_sections(), 1)
}

test "report_add_plot" {
  let doc = @src.ReportDocument::new("Test")
  let doc2 = doc.add_plot("plot1", "Volcano Plot", "  +\n o\n----\n")
  assert_eq(doc2.get_n_sections(), 1)
}

test "report_get_section" {
  let doc = @src.ReportDocument::new("Test")
  let doc2 = doc.add_text("sec1", "Section 1", "Content 1")
  let sec = doc2.get_section("sec1")
  assert_eq(sec.section_id, "sec1")
  assert_eq(sec.content, "Content 1")
}

test "report_get_section_not_found" {
  let doc = @src.ReportDocument::new("Test")
  let sec = doc.get_section("nonexistent")
  assert_eq(sec.section_id, "")
}

test "report_column_new" {
  let col = @src.ReportColumn::new("Gene", ["TP53", "BRCA1"], false)
  assert_eq(col.name, "Gene")
  assert_eq(col.values.length(), 2)
  assert_eq(col.numeric, false)
}

test "report_table_new" {
  let table = @src.ReportTable::new("t1", "Test Table", ["A", "B"], [["x", "y"], ["1", "2"]], "A test caption")
  assert_eq(table.get_n_rows(), 2)
  assert_eq(table.get_n_columns(), 2)
  assert_eq(table.get_column_names(), ["A", "B"])
}

test "report_table_to_ascii" {
  let table = @src.ReportTable::new("t1", "Test", ["Name", "Value"], [["Item1", "Item2"], ["100", "200"]], "Caption")
  let ascii = table.to_ascii()
  assert_true(ascii.contains("Test"))
  assert_true(ascii.contains("Name"))
  assert_true(ascii.contains("Item1"))
}

test "report_render" {
  let doc = @src.ReportDocument::new("My Report")
  let doc2 = doc.add_text("intro", "Introduction", "This is the introduction.")
  let doc3 = doc2.add_text("methods", "Methods", "We used method X.")
  let rendered = doc3.render()
  assert_true(rendered.contains("My Report"))
  assert_true(rendered.contains("Introduction"))
  assert_true(rendered.contains("Methods"))
}

test "report_summary" {
  let doc = @src.ReportDocument::new("Test")
  let doc2 = doc.add_text("sec1", "Section", "Content")
  let summary = doc2.summary()
  assert_true(summary.contains("Report Summary"))
  assert_true(summary.contains("Sections: 1"))
}

test "report_multiple_sections" {
  let doc = @src.ReportDocument::new("Multi")
  let doc2 = doc.add_text("s1", "S1", "C1")
  let doc3 = doc2.add_text("s2", "S2", "C2")
  let doc4 = doc3.add_text("s3", "S3", "C3")
  assert_eq(doc4.get_n_sections(), 3)
}

test "report_table_from_columns" {
  let col1 = @src.ReportColumn::new("Gene", ["TP53", "BRCA1"], false)
  let col2 = @src.ReportColumn::new("log2FC", ["2.5", "-3.0"], true)
  let table = @src.ReportTable::from_columns("t1", "Results", [col1, col2], "DE results")
  assert_eq(table.get_n_columns(), 2)
  assert_eq(table.get_n_rows(), 2)
}

test "report_section_type_to_string" {
  let text = @src.report_section_text()
  assert_eq(text.to_string(), "text")
  let table = @src.report_section_table()
  assert_eq(table.to_string(), "table")
  let plot = @src.report_section_plot()
  assert_eq(plot.to_string(), "plot")
}

test "report_render_with_table" {
  let doc = @src.ReportDocument::new("Test with Table")
  let table = @src.ReportTable::new("results", "Results", ["Gene", "Value"], [["A", "B"], ["1", "2"]], "")
  let doc2 = doc.add_table("tab", table)
  let rendered = doc2.render()
  assert_true(rendered.contains("Results"))
}

test "report_empty_render" {
  let doc = @src.ReportDocument::new("Empty Report")
  let rendered = doc.render()
  assert_true(rendered.contains("Empty Report"))
  assert_true(rendered.contains("Sections: 0"))
}

test "report_padding_func" {
  // Test that the pad logic works through to_ascii
  let table = @src.ReportTable::new("t1", "Pad Test", ["Col1", "Col2"], [["Short", "LongerValue"], ["1", "2"]], "")
  let ascii = table.to_ascii()
  assert_true(ascii.contains("Short"))
  assert_true(ascii.contains("LongerValue"))
}