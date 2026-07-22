///|
/// Taxonomy example: NCBI taxonomy database operations.
///
/// Run: moon run examples/taxonomy_demo/main.mbt
fn main {
  println("=== Bio.Taxonomy Demo ===\n")

  let db = @bio.create_example_taxonomy()

  println("1. Database Statistics:")
  println("   Total taxa: \{db.count_taxa()}")
  println("   Species count: \{db.get_all_species().length()}")
  println("   Genera count: \{db.get_all_genera().length()}")

  println("\n2. Query by TaxID:")
  match db.get_scientific_name("9606") {
    Some(name) => println("   TaxID 9606: \{name}")
    None => println("   TaxID 9606: Not found")
  }

  println("\n3. Query by Name:")
  match db.get_taxid_by_name("Homo sapiens") {
    Some(id) => println("   Homo sapiens: tax_id = \{id}")
    None => println("   Homo sapiens: Not found")
  }

  println("\n4. Lineage of Homo sapiens (9606):")
  let lineage = db.get_lineage("9606")
  let mut i = 0
  while i < lineage.length() {
    println("   \{i + 1}. \{lineage[i]}")
    i = i + 1
  }

  println("\n5. Common Ancestor:")
  match db.get_common_ancestor("9606", "9599") {
    Some(ancestor_id) => {
      match db.get_scientific_name(ancestor_id) {
        Some(name) => println("   Common ancestor of human and chimp: \{name} (tax_id: \{ancestor_id})")
        None => println("   Common ancestor: tax_id \{ancestor_id}")
      }
    }
    None => println("   No common ancestor found")
  }

  println("\n6. Taxonomic Distance:")
  let dist = db.get_distance("9606", "9599")
  println("   Distance between human and chimp: \{dist}")

  let dist2 = db.get_distance("9606", "562")
  println("   Distance between human and E. coli: \{dist2}")

  println("\n7. Ancestor Check:")
  let is_ancestor = db.is_ancestor("9604", "9606")
  println("   Is Hominidae (9604) ancestor of Homo sapiens (9606)? \{is_ancestor}")

  println("\n8. Common Names:")
  let common_names = db.get_common_names("9606")
  let mut names_str = "   Human common names: "
  let mut j = 0
  while j < common_names.length() {
    if j > 0 {
      names_str = names_str + ", "
    }
    names_str = names_str + common_names[j]
    j = j + 1
  }
  println(names_str)

  println("\n9. Synonyms:")
  let synonyms = db.get_synonyms("9606")
  let mut syn_str = "   Human synonyms: "
  let mut k = 0
  while k < synonyms.length() {
    if k > 0 {
      syn_str = syn_str + ", "
    }
    syn_str = syn_str + synonyms[k]
    k = k + 1
  }
  println(syn_str)

  println("\n10. Get all taxa by rank 'species':")
  let species = db.get_all_species()
  let mut l = 0
  while l < species.length() {
    println("   - \{species[l].scientific_name}")
    l = l + 1
  }

  println("\n11. Parse nodes.dmp format:")
  let nodes_content = "1\t|\t0\t|\tno rank\t|\n131567\t|\t1\t|\tno rank\t|\n2759\t|\t131567\t|\tsuperkingdom\t|\n"
  let parsed_db = @bio.parse_nodes_dmp(@bio.TaxonomyDatabase::new(), nodes_content)
  println("   Parsed \{parsed_db.count_taxa()} taxa from nodes.dmp")

  println("\n=== Demo Complete ===")
}