///|
/// Taxonomy module tests

test "TaxonomyDatabase::new" {
  let db = @src.TaxonomyDatabase::new()
  assert_eq(db.count_taxa(), 0)
}

test "TaxonomyDatabase::add_taxon" {
  let db = @src.TaxonomyDatabase::new()
  let taxon = @src.Taxon::new("1", "0", "no rank").set_scientific_name("root")
  ignore(db.add_taxon(taxon))
  assert_eq(db.count_taxa(), 1)
}

test "TaxonomyDatabase::get_taxon" {
  let db = @src.TaxonomyDatabase::new()
  let taxon = @src.Taxon::new("1", "0", "no rank").set_scientific_name("root")
  ignore(db.add_taxon(taxon))
  match db.get_taxon("1") {
    Some(t) => assert_eq(t.scientific_name, "root")
    None => assert_true(false)
  }
}

test "TaxonomyDatabase::get_taxid_by_name" {
  let db = @src.create_example_taxonomy()
  match db.get_taxid_by_name("Homo sapiens") {
    Some(id) => assert_eq(id, "9606")
    None => assert_true(false)
  }
}

test "TaxonomyDatabase::get_lineage" {
  let db = @src.create_example_taxonomy()
  let lineage = db.get_lineage("9606")
  assert_true(lineage.length() > 0)
  assert_eq(lineage[lineage.length() - 1], "Homo sapiens")
}

test "TaxonomyDatabase::get_ancestors" {
  let db = @src.create_example_taxonomy()
  let ancestors = db.get_ancestors("9606")
  assert_true(ancestors.length() > 0)
}

test "TaxonomyDatabase::get_common_ancestor" {
  let db = @src.create_example_taxonomy()
  match db.get_common_ancestor("9606", "9599") {
    Some(id) => assert_eq(id, "9604")
    None => assert_true(false)
  }
}

test "TaxonomyDatabase::is_ancestor" {
  let db = @src.create_example_taxonomy()
  assert_true(db.is_ancestor("9604", "9606"))
  assert_true(!db.is_ancestor("9606", "9604"))
}

test "TaxonomyDatabase::get_distance" {
  let db = @src.create_example_taxonomy()
  let dist = db.get_distance("9606", "9599")
  assert_true(dist > 0)
}

test "TaxonomyDatabase::get_taxa_by_rank" {
  let db = @src.create_example_taxonomy()
  let species = db.get_taxa_by_rank("species")
  assert_true(species.length() >= 2)
}

test "parse_nodes_dmp" {
  let content = "1\t|\t0\t|\tno rank\t|\n131567\t|\t1\t|\tno rank\t|\n2759\t|\t131567\t|\tsuperkingdom\t|\n"
  let db = @src.parse_nodes_dmp(@src.TaxonomyDatabase::new(), content)
  assert_eq(db.count_taxa(), 3)
}

test "parse_names_dmp" {
  let db = @src.TaxonomyDatabase::new()
  ignore(db.add_taxon(@src.Taxon::new("1", "0", "no rank")))
  let content = "1\t|\troot\t|\t\t|\tscientific name\t|\n"
  let db2 = @src.parse_names_dmp(db, content)
  match db2.get_scientific_name("1") {
    Some(name) => assert_eq(name, "root")
    None => assert_true(false)
  }
}

test "create_example_taxonomy" {
  let db = @src.create_example_taxonomy()
  assert_true(db.count_taxa() > 10)
}

test "get_all_species" {
  let db = @src.create_example_taxonomy()
  let species = db.get_all_species()
  assert_true(species.length() >= 2)
}

test "get_common_names" {
  let db = @src.create_example_taxonomy()
  let names = db.get_common_names("9606")
  assert_true(contains_string(names, "human"))
}

fn contains_string(arr : Array[String], s : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == s {
      return true
    }
    i = i + 1
  }
  false
}