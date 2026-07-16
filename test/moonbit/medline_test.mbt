///|
test "medline_record_new" {
  let authors = ["Smith J", "Johnson A"]
  let mesh_terms = ["Genetics", "Molecular Biology"]

  let record = @src.MedlineRecord::new(
    "12345678",
    "Title of article",
    "Abstract text...",
    authors,
    "Journal of Biology",
    "10",
    "2",
    "100-105",
    "2024/03/15",
    mesh_terms,
    [],
    [],
    "eng",
    ["Journal Article"],
    "10.1234/jbio.2024.10.2.100",
    "PMC1234567",
  )

  assert_eq(record.pmid, "12345678")
  assert_eq(record.title, "Title of article")
  assert_eq(record.authors.length(), 2)
  assert_eq(record.mesh_terms.length(), 2)
}

///|
test "medline_parse" {
  let text = "PMID- 12345678\nTI  Title of the article\nAB  This is the abstract text.\nAU  Smith J\nAU  Johnson A\nJT  Journal of Biology\nVI  10\nIP  2\nPG  100-105\nDP  2024/03/15\nMH  Genetics\nMH  Molecular Biology\nLA  eng\nPT  Journal Article\nAID doi:10.1234/jbio.2024.10.2.100"

  let record = @src.parse_medline(text)

  assert_true(record.is_some())
  assert_eq(record.unwrap().pmid, "12345678")
  assert_eq(record.unwrap().title, "Title of the article")
  assert_eq(record.unwrap().authors.length(), 2)
  assert_eq(record.unwrap().mesh_terms.length(), 2)
  assert_eq(record.unwrap().doi, "10.1234/jbio.2024.10.2.100")
}

///|
test "medline_parse_multiple" {
  let text = "PMID- 12345678\nTI  Article 1\nAU  Smith J\n\nPMID- 87654321\nTI  Article 2\nAU  Johnson A"

  let records = @src.parse_medline_records(text)

  assert_eq(records.length(), 2)
  assert_eq(records[0].pmid, "12345678")
  assert_eq(records[1].pmid, "87654321")
}

///|
test "medline_format_citation_apa" {
  let authors = ["Smith J", "Johnson A"]
  let record = @src.MedlineRecord::new(
    "12345678",
    "Title of article",
    "",
    authors,
    "Journal of Biology",
    "10",
    "2",
    "100-105",
    "2024/03/15",
    [],
    [],
    [],
    "eng",
    [],
    "",
    "",
  )

  let citation = @src.format_citation_apa(record)

  assert_true(citation.contains("Smith J"))
  assert_true(citation.contains("(2024)"))
  assert_true(citation.contains("Title of article"))
}

///|
test "medline_get_first_author" {
  let authors = ["Smith J", "Johnson A"]
  let record = @src.MedlineRecord::new(
    "12345678",
    "Title",
    "",
    authors,
    "",
    "",
    "",
    "",
    "",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  assert_eq(@src.get_first_author(record), "Smith J")
}

///|
test "medline_has_abstract" {
  let record_with = @src.MedlineRecord::new(
    "12345678",
    "Title",
    "This is an abstract.",
    [],
    "",
    "",
    "",
    "",
    "",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  let record_without = @src.MedlineRecord::new(
    "87654321",
    "Title",
    "",
    [],
    "",
    "",
    "",
    "",
    "",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  assert_true(@src.has_abstract(record_with))
  assert_true(!@src.has_abstract(record_without))
}

///|
test "medline_get_publication_year" {
  let record = @src.MedlineRecord::new(
    "12345678",
    "Title",
    "",
    [],
    "",
    "",
    "",
    "",
    "2024/03/15",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  assert_eq(@src.get_publication_year(record), "2024")
}

///|
test "medline_filter_by_mesh" {
  let record1 = @src.MedlineRecord::new(
    "12345678",
    "Article 1",
    "",
    [],
    "",
    "",
    "",
    "",
    "2024",
    ["Genetics", "Cancer"],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  let record2 = @src.MedlineRecord::new(
    "87654321",
    "Article 2",
    "",
    [],
    "",
    "",
    "",
    "",
    "2024",
    ["Biochemistry"],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  let records = [record1, record2]
  let filtered = @src.filter_by_mesh(records, "Genetics")

  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].pmid, "12345678")
}

///|
test "medline_count_by_year" {
  let record1 = @src.MedlineRecord::new(
    "1",
    "A",
    "",
    [],
    "",
    "",
    "",
    "",
    "2022",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )
  let record2 = @src.MedlineRecord::new(
    "2",
    "B",
    "",
    [],
    "",
    "",
    "",
    "",
    "2023",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )
  let record3 = @src.MedlineRecord::new(
    "3",
    "C",
    "",
    [],
    "",
    "",
    "",
    "",
    "2023",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )
  let record4 = @src.MedlineRecord::new(
    "4",
    "D",
    "",
    [],
    "",
    "",
    "",
    "",
    "2024",
    [],
    [],
    [],
    "",
    [],
    "",
    "",
  )

  let records = [record1, record2, record3, record4]
  let counts = @src.count_by_year(records)

  assert_eq(counts.length(), 3)
  assert_eq(counts[0].0, "2022")
  assert_eq(counts[0].1, 1)
  assert_eq(counts[1].0, "2023")
  assert_eq(counts[1].1, 2)
}
