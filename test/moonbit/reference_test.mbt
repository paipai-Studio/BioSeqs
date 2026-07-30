///|
/// Test file for reference module.

test "reference_create_basic" {
  let r = @src.BioReference::new(title="Test Paper", authors="Smith J", journal="Nature", year="2024")
  assert_eq!(r.title(), "Test Paper")
  assert_eq!(r.authors(), "Smith J")
  assert_eq!(r.journal(), "Nature")
  assert_eq!(r.year(), "2024")
  assert_eq!(r.pubmed_id(), "")
  assert_eq!(r.doi(), "")
  assert_eq!(r.reference_type(), "journal article")
}

test "reference_with_pubmed" {
  let r = @src.BioReference::with_pubmed("CRISPR Advances", "Doudna J", "Science", "2020", "12345678")
  assert_eq!(r.pubmed_id(), "12345678")
  assert_eq!(r.title(), "CRISPR Advances")
  let citation = r.citation()
  assert_eq!(citation.contains("12345678"), true)
  assert_eq!(citation.contains("Doudna J"), true)
}

test "reference_with_doi" {
  let r = @src.BioReference::with_doi("Protein Folding", "Jones A", "Cell", "2023", "10.1000/test")
  assert_eq!(r.doi(), "10.1000/test")
  assert_eq!(r.citation().contains("10.1000/test"), true)
}

test "reference_setters" {
  let r = @src.BioReference::new()
  r.set_title("New Title")
  r.set_authors("New Author")
  r.set_journal("New Journal")
  r.set_year("2025")
  r.set_pubmed_id("99999")
  r.set_doi("10.999/test")
  r.set_type("book")
  assert_eq!(r.title(), "New Title")
  assert_eq!(r.authors(), "New Author")
  assert_eq!(r.journal(), "New Journal")
  assert_eq!(r.year(), "2025")
  assert_eq!(r.pubmed_id(), "99999")
  assert_eq!(r.doi(), "10.999/test")
  assert_eq!(r.reference_type(), "book")
}

test "reference_locations" {
  let r = @src.BioReference::new(title="Test", authors="A")
  r.add_location(1, 100)
  r.add_location(200, 350)
  assert_eq!(r.get_n_locations(), 2)
  let (s1, e1) = r.get_location(0)
  assert_eq!(s1, 1)
  assert_eq!(e1, 100)
  let (s2, e2) = r.get_location(1)
  assert_eq!(s2, 200)
  assert_eq!(e2, 350)
  // Out of bounds returns (0,0)
  let (s3, e3) = r.get_location(5)
  assert_eq!(s3, 0)
  assert_eq!(e3, 0)
}

test "reference_comment" {
  let r = @src.BioReference::new(title="Test", authors="A")
  r.set_comment("Important discovery")
  assert_eq!(r.get_comment(), "Important discovery")
}

test "reference_citation" {
  let r = @src.BioReference::with_pubmed("Genome Study", "Smith J, Jones A", "Science", "2022", "12345")
  let citation = r.citation()
  assert_eq!(citation.contains("Smith J"), true)
  assert_eq!(citation.contains("Genome Study"), true)
  assert_eq!(citation.contains("Science"), true)
  assert_eq!(citation.contains("2022"), true)
  assert_eq!(citation.contains("PMID: 12345"), true)
}

test "reference_list_basic" {
  let list = @src.ReferenceList::new()
  assert_eq!(list.count(), 0)
  let r = @src.BioReference::new(title="Test", authors="A")
  list.add(r)
  assert_eq!(list.count(), 1)
  let retrieved = list.get(0)
  assert_eq!(retrieved.title(), "Test")
}

test "reference_list_by_author" {
  let list = @src.ReferenceList::new()
  list.add(@src.BioReference::new(title="Paper 1", authors="Smith J"))
  list.add(@src.BioReference::new(title="Paper 2", authors="Jones A"))
  list.add(@src.BioReference::new(title="Paper 3", authors="Smith J, Lee K"))
  let smith_papers = list.by_author("Smith J")
  assert_eq!(smith_papers.length(), 2)
  let jones_papers = list.by_author("Jones A")
  assert_eq!(jones_papers.length(), 1)
  let doe_papers = list.by_author("Doe X")
  assert_eq!(doe_papers.length(), 0)
}

test "reference_list_by_year" {
  let list = @src.ReferenceList::new()
  list.add(@src.BioReference::new(title="Paper 1", year="2020"))
  list.add(@src.BioReference::new(title="Paper 2", year="2022"))
  list.add(@src.BioReference::new(title="Paper 3", year="2020"))
  let y2020 = list.by_year("2020")
  assert_eq!(y2020.length(), 2)
  let y2022 = list.by_year("2022")
  assert_eq!(y2022.length(), 1)
}

test "reference_list_with_pubmed" {
  let list = @src.ReferenceList::new()
  list.add(@src.BioReference::with_pubmed("P1", "A", "J", "2020", "11111"))
  list.add(@src.BioReference::new(title="P2", authors="B"))
  list.add(@src.BioReference::with_pubmed("P3", "C", "J", "2022", "22222"))
  let with_pm = list.with_pubmed()
  assert_eq!(with_pm.length(), 2)
}

test "reference_list_bibliography" {
  let list = @src.ReferenceList::new()
  list.add(@src.BioReference::new(title="Paper 1", authors="Smith J", journal="Nature", year="2020"))
  list.add(@src.BioReference::new(title="Paper 2", authors="Jones A", journal="Science", year="2022"))
  let bib = list.bibliography()
  assert_eq!(bib.contains("1. Smith J"), true)
  assert_eq!(bib.contains("2. Jones A"), true)
  assert_eq!(bib.contains("Nature"), true)
  assert_eq!(bib.contains("Science"), true)
}

test "reference_sample_data" {
  let list = @src.reference_sample_data()
  assert_eq!(list.count(), 3)
  let ref1 = list.get(0)
  assert_eq!(ref1.title(), "The Human Genome: A Complete Sequence")
  assert_eq!(ref1.pubmed_id(), "36189102")
  assert_eq!(ref1.get_n_locations(), 1)
}
