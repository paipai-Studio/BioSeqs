///|
/// reference module example - Literature references for SeqRecord.
///
/// Run: moon run examples/reference_demo
fn main {
  println("=== Bio.Reference Module Demo ===")
  println("Literature Reference Objects for SeqRecord")
  println("")

  // Create basic reference
  println("=== Basic Reference ===")
  let r1 = @bio.BioReference::new(
    title="The Human Genome Sequence",
    authors="Nurk S, Koren S, Rhie A",
    journal="Science",
    year="2022",
  )
  println("Title: " + r1.title())
  println("Authors: " + r1.authors())
  println("Journal: " + r1.journal())
  println("Year: " + r1.year())
  println("Citation: " + r1.citation())
  println("")

  // Reference with PubMed
  println("=== Reference with PubMed ID ===")
  let r2 = @bio.BioReference::with_pubmed(
    "CRISPR-Cas9 Advances",
    "Doudna JA, Charpentier E",
    "Science",
    "2014",
    "24700447",
  )
  println("PubMed ID: " + r2.pubmed_id())
  println("Citation: " + r2.citation())
  println("")

  // Reference with DOI
  println("=== Reference with DOI ===")
  let r3 = @bio.BioReference::with_doi(
    "AlphaFold Protein Structure Prediction",
    "Jumper J, Evans R, Pritzel A",
    "Nature",
    "2021",
    "10.1038/s41586-021-03819-2",
  )
  println("DOI: " + r3.doi())
  println("Citation: " + r3.citation())
  println("")

  // Setters and locations
  println("=== Setters and Locations ===")
  let r4 = @bio.BioReference::new()
  r4.set_title("COVID-19 Vaccine Development")
  r4.set_authors("Kariko K, Weissman C")
  r4.set_journal("Cell")
  r4.set_year("2020")
  r4.set_pubmed_id("33064427")
  r4.set_type("book chapter")
  r4.add_location(1, 5000)
  r4.add_location(10000, 15000)
  r4.set_comment("mRNA vaccine technology breakthrough")
  println("Title: " + r4.title())
  println("Type: " + r4.reference_type())
  println("Locations: " + r4.get_n_locations().to_string())
  println("Comment: " + r4.get_comment())
  println("Citation: " + r4.citation())
  println("")

  // Reference List
  println("=== Reference List ===")
  let list = @bio.ReferenceList::new()
  list.add(r1)
  list.add(r2)
  list.add(r3)
  list.add(r4)
  println("Total references: " + list.count().to_string())
  println("")

  // Query by author
  println("=== Query by Author ===")
  let doudna_papers = list.by_author("Doudna")
  println("Doudna papers: " + doudna_papers.length().to_string())
  println("")

  // Query by year
  println("=== Query by Year ===")
  let papers_2022 = list.by_year("2022")
  println("2022 papers: " + papers_2022.length().to_string())
  println("")

  // Query with PubMed
  println("=== Query with PubMed ID ===")
  let with_pm = list.with_pubmed()
  println("Papers with PubMed IDs: " + with_pm.length().to_string())
  println("")

  // Bibliography
  println("=== Full Bibliography ===")
  println(list.bibliography())

  // Sample data
  println("=== Sample Data ===")
  let sample = @bio.reference_sample_data()
  println("Sample references: " + sample.count().to_string())
  let mut i = 0
  while i < sample.count() {
    println("  " + (i + 1).to_string() + ". " + sample.get(i).title() + " (" + sample.get(i).year() + ")")
    i = i + 1
  }
}
