///|
fn main {
  println("=== Medline/PubMed Demo ===")

  let medline_text = "PMID- 12345678\nTI  - Bioinformatics approaches for gene expression analysis\nAB  - This study demonstrates the application of bioinformatics methods for analyzing gene expression data using RNA-seq technology.\nAU  - Smith J\nAU  - Johnson A\nAU  - Williams B\nJT  - Bioinformatics Journal\nVI  - 20\nIP  - 5\nPG  - 100-105\nDP  - 2024/03/15\nMH  - Gene Expression\nMH  - Bioinformatics\nMH  - Computational Biology\nKW  - RNA-seq\nKW  - Differential Expression\nAD  - Department of Biology, University of XYZ\nLA  - eng\nPT  - Journal Article\nPT  - Research Support Non-US Gov\nAID - doi:10.1234/bioinfo.2024.20.5.100\nAID - PMC1234567"

  println("\n1. Parse Single Medline Record:")
  println("===============================")
  let record = @bio.parse_medline(medline_text)
  if record.is_some() {
    let r = record.unwrap()
    println("PMID: " + r.pmid)
    println("Title: " + r.title)
    println("Abstract: " + r.abstract_text[0:50].to_owned() + "...")
    println("Journal: " + r.journal)
    println("Volume/Issue: " + r.volume + "(" + r.issue + ")")
    println("Pages: " + r.pages)
    println("Publication Date: " + r.publication_date)
    println("Language: " + r.language)
    println("DOI: " + r.doi)
    println("PMC: " + r.pmc)
    println("Affiliations: " + r.affiliations.length().to_string())

    println("\nAuthors (" + r.authors.length().to_string() + "):")
    let mut i = 0
    while i < r.authors.length() {
      println("  - " + r.authors[i])
      i = i + 1
    }

    println("\nMeSH Terms (" + r.mesh_terms.length().to_string() + "):")
    i = 0
    while i < r.mesh_terms.length() {
      println("  - " + r.mesh_terms[i])
      i = i + 1
    }

    println("\nKeywords (" + r.keywords.length().to_string() + "):")
    i = 0
    while i < r.keywords.length() {
      println("  - " + r.keywords[i])
      i = i + 1
    }

    println("\nPublication Types (" + r.pub_types.length().to_string() + "):")
    i = 0
    while i < r.pub_types.length() {
      println("  - " + r.pub_types[i])
      i = i + 1
    }
  }

  println("\n2. APA Citation Format:")
  println("=======================")
  if record.is_some() {
    let citation = @bio.format_citation_apa(record.unwrap())
    println(citation)
  }

  println("\n3. Record Metadata:")
  println("===================")
  if record.is_some() {
    let r = record.unwrap()
    println("First Author: " + @bio.get_first_author(r))
    println("Has Abstract: " + @bio.has_abstract(r).to_string())
    println("Publication Year: " + @bio.get_publication_year(r))
    println("Author Count: " + @bio.get_author_count(r).to_string())
  }

  let multiple_text = "PMID- 12345678\nTI  - Gene Expression Analysis Using RNA-seq\nAU  - Smith J\nAU  - Johnson A\nDP  - 2022/06/15\nMH  - Gene Expression\nMH  - RNA-seq\n\nPMID- 87654321\nTI  - Cancer Genetics and Genomics\nAU  - Johnson A\nAU  - Davis C\nDP  - 2023/01/20\nMH  - Cancer\nMH  - Genetics\nMH  - Genomics\n\nPMID- 11112222\nTI  - Computational Biology Approaches\nAU  - Williams B\nDP  - 2023/08/10\nMH  - Computational Biology\nMH  - Bioinformatics\n\nPMID- 33334444\nTI  - Advanced RNA-seq Technology\nAU  - Davis C\nAU  - Wilson E\nDP  - 2024/02/28\nMH  - RNA-seq\nMH  - Next-Generation Sequencing"

  println("\n4. Parse Multiple Records:")
  println("==========================")
  let records = @bio.parse_medline_records(multiple_text)
  println("Total records parsed: " + records.length().to_string())

  let mut j = 0
  while j < records.length() {
    println("  PMID: " + records[j].pmid + " - " + records[j].title)
    j = j + 1
  }

  println("\n5. Filter by MeSH Term:")
  println("=======================")

  let genetics_records = @bio.filter_by_mesh(records, "Genetics")
  println(
    "Records with 'Genetics' MeSH: " + genetics_records.length().to_string(),
  )
  j = 0
  while j < genetics_records.length() {
    println(
      "  - " + genetics_records[j].pmid + ": " + genetics_records[j].title,
    )
    j = j + 1
  }

  let rnaseq_records = @bio.filter_by_mesh(records, "RNA-seq")
  println(
    "\nRecords with 'RNA-seq' MeSH: " + rnaseq_records.length().to_string(),
  )
  j = 0
  while j < rnaseq_records.length() {
    println("  - " + rnaseq_records[j].pmid + ": " + rnaseq_records[j].title)
    j = j + 1
  }

  println("\n6. Filter by Author:")
  println("====================")
  let smith_records = @bio.filter_by_author(records, "Smith")
  println("Records by 'Smith': " + smith_records.length().to_string())
  j = 0
  while j < smith_records.length() {
    println("  - " + smith_records[j].pmid + ": " + smith_records[j].title)
    j = j + 1
  }

  println("\n7. Count Publications by Year:")
  println("===============================")
  let year_counts = @bio.count_by_year(records)
  j = 0
  while j < year_counts.length() {
    println(
      "  " +
      year_counts[j].0 +
      ": " +
      year_counts[j].1.to_string() +
      " publications",
    )
    j = j + 1
  }

  println("\n8. Extract All PMIDs:")
  println("=====================")
  let pmids = @bio.extract_pmids(records)
  j = 0
  while j < pmids.length() {
    println("  " + pmids[j])
    j = j + 1
  }

  println("\n=== Demo Complete ===")
}
