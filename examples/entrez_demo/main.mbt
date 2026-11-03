///|
/// NCBI Entrez database access demo inspired by Biopython Bio.Entrez.
///
/// Run: moon run examples/entrez_demo/main.mbt
fn main {
  println("=== Bio.Entrez Demo ===")

  println("\n1. ESearch - Search PubMed:")
  let results = @bio.esearch(@bio.entrez_db_pubmed(), "CRISPR")
  for result in results {
    println(
      "   ID: \{result.id}, Title: \{result.title}, Score: \{result.score}",
    )
  }

  println("\n2. EFetch - Fetch PubMed article:")
  let xml = @bio.efetch(@bio.entrez_db_pubmed(), "37566892")
  println("   XML: \{xml[0:80]}...")

  println("\n3. Parse PubMed article:")
  match @bio.parse_pubmed(xml) {
    Some(article) => {
      println("   PMID: \{article.pmid}")
      println("   Title: \{article.title}")
      println("   Journal: \{article.journal}")
    }
    None => println("   Failed to parse")
  }

  println("\n4. ESearch - Search Gene:")
  let gene_results = @bio.esearch(@bio.entrez_db_gene(), "TP53")
  for result in gene_results {
    println("   ID: \{result.id}, Title: \{result.title}")
  }

  println("\n5. EFetch - Fetch Gene record:")
  let gene_xml = @bio.efetch(@bio.entrez_db_gene(), "7157")
  println("   XML: \{gene_xml[0:80]}...")

  println("\n6. Parse Gene record:")
  match @bio.parse_gene(gene_xml) {
    Some(gene) => {
      println("   GeneID: \{gene.gene_id}")
      println("   Symbol: \{gene.symbol}")
      println("   Description: \{gene.description}")
      println("   Chromosome: \{gene.chromosome}")
    }
    None => println("   Failed to parse")
  }

  println("\n7. Get gene info:")
  match @bio.get_gene_info("TP53") {
    Some(gene) => {
      println("   GeneID: \{gene.gene_id}")
      println("   Symbol: \{gene.symbol}")
      println("   Organism: \{gene.organism}")
    }
    None => println("   Gene not found")
  }

  println("\n8. Get taxonomy:")
  match @bio.get_taxonomy("Homo sapiens") {
    Some(tax) => {
      println("   TaxID: \{tax.tax_id}")
      println("   Scientific Name: \{tax.scientific_name}")
      println("   Common Name: \{tax.common_name}")
      println("   Rank: \{tax.rank}")
    }
    None => println("   Taxonomy not found")
  }

  println("\n9. EGQuery - Global query:")
  let global = @bio.egquery("cancer")
  for kv in global {
    println("   \{kv.0}: \{kv.1}")
  }

  println("\n10. EInfo - Database info:")
  let info = @bio.einfo()
  for kv in info {
    println("   \{kv.0}: \{kv.1}")
  }
}
