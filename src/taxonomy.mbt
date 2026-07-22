///|
/// Bio.Taxonomy - NCBI taxonomy database support.
/// Inspired by Biopython Bio.Taxonomy.
/// Provides NCBI taxonomy database parsing and taxonomy tree operations.

///|
pub struct Taxon {
  tax_id : String
  parent_tax_id : String
  rank : String
  scientific_name : String
  common_names : Array[String]
  synonyms : Array[String]
  lineage : Array[String]
}

///|
pub fn Taxon::new(tax_id : String, parent_tax_id : String, rank : String) -> Taxon {
  Taxon::{
    tax_id,
    parent_tax_id,
    rank,
    scientific_name: "",
    common_names: Array::new(),
    synonyms: Array::new(),
    lineage: Array::new()
  }
}

///|
pub fn Taxon::set_scientific_name(self : Taxon, name : String) -> Taxon {
  Taxon::{
    tax_id: self.tax_id,
    parent_tax_id: self.parent_tax_id,
    rank: self.rank,
    scientific_name: name,
    common_names: self.common_names,
    synonyms: self.synonyms,
    lineage: self.lineage
  }
}

///|
pub fn Taxon::add_common_name(self : Taxon, name : String) -> Taxon {
  let new_names = self.common_names
  new_names.push(name)
  Taxon::{
    tax_id: self.tax_id,
    parent_tax_id: self.parent_tax_id,
    rank: self.rank,
    scientific_name: self.scientific_name,
    common_names: new_names,
    synonyms: self.synonyms,
    lineage: self.lineage
  }
}

///|
pub fn Taxon::add_synonym(self : Taxon, name : String) -> Taxon {
  let new_synonyms = self.synonyms
  new_synonyms.push(name)
  Taxon::{
    tax_id: self.tax_id,
    parent_tax_id: self.parent_tax_id,
    rank: self.rank,
    scientific_name: self.scientific_name,
    common_names: self.common_names,
    synonyms: new_synonyms,
    lineage: self.lineage
  }
}

///|
pub fn Taxon::set_lineage(self : Taxon, lineage : Array[String]) -> Taxon {
  Taxon::{
    tax_id: self.tax_id,
    parent_tax_id: self.parent_tax_id,
    rank: self.rank,
    scientific_name: self.scientific_name,
    common_names: self.common_names,
    synonyms: self.synonyms,
    lineage: lineage
  }
}

///|
pub struct TaxonomyDatabase {
  taxa : @hashmap.HashMap[String, Taxon]
  name_to_taxid : @hashmap.HashMap[String, String]
}

///|
pub fn TaxonomyDatabase::new() -> TaxonomyDatabase {
  TaxonomyDatabase::{
    taxa: @hashmap.HashMap([]),
    name_to_taxid: @hashmap.HashMap([])
  }
}

///|
pub fn TaxonomyDatabase::add_taxon(self : TaxonomyDatabase, taxon : Taxon) -> TaxonomyDatabase {
  self.taxa.set(taxon.tax_id, taxon)
  if taxon.scientific_name != "" {
    self.name_to_taxid.set(to_lowercase(taxon.scientific_name), taxon.tax_id)
  }
  self
}

///|
pub fn TaxonomyDatabase::get_taxon(self : TaxonomyDatabase, tax_id : String) -> Taxon? {
  if self.taxa.contains(tax_id) {
    Some(self.taxa.get(tax_id).unwrap())
  } else {
    None
  }
}

///|
pub fn TaxonomyDatabase::get_taxid_by_name(self : TaxonomyDatabase, name : String) -> String? {
  if self.name_to_taxid.contains(to_lowercase(name)) {
    Some(self.name_to_taxid.get(to_lowercase(name)).unwrap())
  } else {
    None
  }
}

///|
pub fn TaxonomyDatabase::get_ancestors(self : TaxonomyDatabase, tax_id : String) -> Array[Taxon] {
  let result : Array[Taxon] = Array::new()
  let mut current_id = tax_id
  
  while current_id != "" && current_id != "1" && current_id != "0" {
    match self.get_taxon(current_id) {
      Some(taxon) => {
        if taxon.parent_tax_id != "" && taxon.parent_tax_id != "1" && taxon.parent_tax_id != "0" {
          match self.get_taxon(taxon.parent_tax_id) {
            Some(parent) => {
              result.push(parent)
              current_id = parent.tax_id
            }
            None => break
          }
        } else {
          break
        }
      }
      None => break
    }
  }
  
  result
}

///|
pub fn TaxonomyDatabase::get_descendants(self : TaxonomyDatabase, tax_id : String) -> Array[Taxon] {
  let result : Array[Taxon] = Array::new()
  ignore(find_descendants(self, tax_id, result))
  result
}

///|
fn find_descendants(db : TaxonomyDatabase, parent_id : String, result : Array[Taxon]) -> Unit {
  let entries = db.taxa.iter().collect()
  let mut i = 0
  while i < entries.length() {
    let (id, taxon) = entries[i]
    if taxon.parent_tax_id == parent_id {
      result.push(taxon)
      ignore(find_descendants(db, id, result))
    }
    i = i + 1
  }
}

///|
pub fn TaxonomyDatabase::get_lineage(self : TaxonomyDatabase, tax_id : String) -> Array[String] {
  let ancestors = self.get_ancestors(tax_id)
  let result : Array[String] = Array::new()
  
  let mut i = ancestors.length() - 1
  while i >= 0 {
    result.push(ancestors[i].scientific_name)
    i = i - 1
  }
  
  match self.get_taxon(tax_id) {
    Some(taxon) => result.push(taxon.scientific_name)
    None => ()
  }
  
  result
}

///|
pub fn TaxonomyDatabase::get_rank(self : TaxonomyDatabase, tax_id : String) -> String? {
  match self.get_taxon(tax_id) {
    Some(taxon) => Some(taxon.rank)
    None => None
  }
}

///|
pub fn TaxonomyDatabase::get_scientific_name(self : TaxonomyDatabase, tax_id : String) -> String? {
  match self.get_taxon(tax_id) {
    Some(taxon) => Some(taxon.scientific_name)
    None => None
  }
}

///|
pub fn TaxonomyDatabase::get_common_names(self : TaxonomyDatabase, tax_id : String) -> Array[String] {
  match self.get_taxon(tax_id) {
    Some(taxon) => taxon.common_names
    None => Array::new()
  }
}

///|
pub fn TaxonomyDatabase::get_synonyms(self : TaxonomyDatabase, tax_id : String) -> Array[String] {
  match self.get_taxon(tax_id) {
    Some(taxon) => taxon.synonyms
    None => Array::new()
  }
}

///|
pub fn parse_nodes_dmp(db : TaxonomyDatabase, content : String) -> TaxonomyDatabase {
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line = trim(lines[i])
    if line != "" {
      let fields = split_by_char(line, '\t'.to_int())
      let pipe_fields = Array::new()
      let mut j = 0
      while j < fields.length() {
        let f = trim(fields[j])
        if f != "|" {
          pipe_fields.push(f)
        }
        j = j + 1
      }
      if pipe_fields.length() >= 3 {
        let tax_id = pipe_fields[0]
        let parent_tax_id = pipe_fields[1]
        let rank = if pipe_fields.length() > 2 { pipe_fields[2] } else { "" }
        
        let taxon = Taxon::new(tax_id, parent_tax_id, rank)
        ignore(db.add_taxon(taxon))
      }
    }
    i = i + 1
  }
  db
}

///|
pub fn parse_names_dmp(db : TaxonomyDatabase, content : String) -> TaxonomyDatabase {
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line = trim(lines[i])
    if line != "" {
      let fields = split_by_char(line, '\t'.to_int())
      let pipe_fields = Array::new()
      let mut j = 0
      while j < fields.length() {
        let f = trim(fields[j])
        if f != "|" {
          pipe_fields.push(f)
        }
        j = j + 1
      }
      if pipe_fields.length() >= 4 {
        let tax_id = pipe_fields[0]
        let name = pipe_fields[1]
        let name_class = if pipe_fields.length() > 3 { pipe_fields[3] } else { "" }
        
        if db.taxa.contains(tax_id) {
          let taxon = db.taxa.get(tax_id).unwrap()
          
          if name_class == "scientific name" {
            let updated = taxon.set_scientific_name(name)
            db.taxa.set(tax_id, updated)
          } else if name_class == "common name" {
            let updated = taxon.add_common_name(name)
            db.taxa.set(tax_id, updated)
          } else if name_class == "synonym" {
            let updated = taxon.add_synonym(name)
            db.taxa.set(tax_id, updated)
          }
        }
      }
    }
    i = i + 1
  }
  db
}

///|


///|
pub fn parse_taxonomy_database(nodes_content : String, names_content : String) -> TaxonomyDatabase {
  let db = TaxonomyDatabase::new()
  let db_with_nodes = parse_nodes_dmp(db, nodes_content)
  parse_names_dmp(db_with_nodes, names_content)
}

///|
pub fn create_example_taxonomy() -> TaxonomyDatabase {
  let db = TaxonomyDatabase::new()
  
  let root = Taxon::new("1", "0", "no rank").set_scientific_name("root")
  ignore(db.add_taxon(root))
  
  let cellular_organisms = Taxon::new("131567", "1", "no rank").set_scientific_name("cellular organisms")
  ignore(db.add_taxon(cellular_organisms))
  
  let eukaryota = Taxon::new("2759", "131567", "superkingdom").set_scientific_name("Eukaryota")
  ignore(db.add_taxon(eukaryota))
  
  let metazoa = Taxon::new("33208", "2759", "kingdom").set_scientific_name("Metazoa")
  ignore(db.add_taxon(metazoa))
  
  let chordata = Taxon::new("7711", "33208", "phylum").set_scientific_name("Chordata")
  ignore(db.add_taxon(chordata))
  
  let mammalia = Taxon::new("40674", "7711", "class").set_scientific_name("Mammalia")
  ignore(db.add_taxon(mammalia))
  
  let primates = Taxon::new("9443", "40674", "order").set_scientific_name("Primates")
  ignore(db.add_taxon(primates))
  
  let hominidae = Taxon::new("9604", "9443", "family").set_scientific_name("Hominidae")
  ignore(db.add_taxon(hominidae))
  
  let homo = Taxon::new("9605", "9604", "genus").set_scientific_name("Homo")
  ignore(db.add_taxon(homo))
  
  let sapiens = Taxon::new("9606", "9605", "species").set_scientific_name("Homo sapiens")
  let sapiens_with_common = sapiens.add_common_name("human")
  let sapiens_with_synonym = sapiens_with_common.add_synonym("Homo sapiens sapiens")
  ignore(db.add_taxon(sapiens_with_synonym))
  
  let pan = Taxon::new("9598", "9604", "genus").set_scientific_name("Pan")
  ignore(db.add_taxon(pan))
  
  let troglodytes = Taxon::new("9599", "9598", "species").set_scientific_name("Pan troglodytes")
  let troglodytes_with_common = troglodytes.add_common_name("chimpanzee")
  ignore(db.add_taxon(troglodytes_with_common))
  
  let bacteria = Taxon::new("2", "131567", "superkingdom").set_scientific_name("Bacteria")
  ignore(db.add_taxon(bacteria))
  
  let proteobacteria = Taxon::new("1224", "2", "phylum").set_scientific_name("Proteobacteria")
  ignore(db.add_taxon(proteobacteria))
  
  let escherichia = Taxon::new("561", "1224", "genus").set_scientific_name("Escherichia")
  ignore(db.add_taxon(escherichia))
  
  let coli = Taxon::new("562", "561", "species").set_scientific_name("Escherichia coli")
  let coli_with_common = coli.add_common_name("E. coli")
  ignore(db.add_taxon(coli_with_common))
  
  db
}

///|
pub fn TaxonomyDatabase::get_taxa_by_rank(self : TaxonomyDatabase, rank : String) -> Array[Taxon] {
  let result : Array[Taxon] = Array::new()
  let entries = self.taxa.iter().collect()
  let mut i = 0
  while i < entries.length() {
    let (_, taxon) = entries[i]
    if taxon.rank == rank {
      result.push(taxon)
    }
    i = i + 1
  }
  result
}

///|
pub fn TaxonomyDatabase::get_all_species(self : TaxonomyDatabase) -> Array[Taxon] {
  self.get_taxa_by_rank("species")
}

///|
pub fn TaxonomyDatabase::get_all_genera(self : TaxonomyDatabase) -> Array[Taxon] {
  self.get_taxa_by_rank("genus")
}

///|
pub fn TaxonomyDatabase::count_taxa(self : TaxonomyDatabase) -> Int {
  let entries = self.taxa.iter().collect()
  entries.length()
}

///|
pub fn TaxonomyDatabase::is_ancestor(self : TaxonomyDatabase, ancestor_id : String, descendant_id : String) -> Bool {
  let ancestors = self.get_ancestors(descendant_id)
  let mut i = 0
  while i < ancestors.length() {
    if ancestors[i].tax_id == ancestor_id {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn TaxonomyDatabase::get_common_ancestor(self : TaxonomyDatabase, tax_id1 : String, tax_id2 : String) -> String? {
  let ancestors1 = self.get_ancestors(tax_id1)
  let ancestors2 = self.get_ancestors(tax_id2)
  
  let mut i = 0
  while i < ancestors1.length() {
    let mut j = 0
    while j < ancestors2.length() {
      if ancestors1[i].tax_id == ancestors2[j].tax_id {
        return Some(ancestors1[i].tax_id)
      }
      j = j + 1
    }
    i = i + 1
  }
  
  None
}

///|
pub fn TaxonomyDatabase::get_distance(self : TaxonomyDatabase, tax_id1 : String, tax_id2 : String) -> Int {
  match self.get_common_ancestor(tax_id1, tax_id2) {
    Some(ancestor_id) => {
      let dist1 = get_distance_to_ancestor(self, tax_id1, ancestor_id)
      let dist2 = get_distance_to_ancestor(self, tax_id2, ancestor_id)
      dist1 + dist2
    }
    None => -1
  }
}

///|
fn get_distance_to_ancestor(db : TaxonomyDatabase, tax_id : String, ancestor_id : String) -> Int {
  let ancestors = db.get_ancestors(tax_id)
  let mut distance = 0
  let mut i = 0
  while i < ancestors.length() {
    distance = distance + 1
    if ancestors[i].tax_id == ancestor_id {
      return distance
    }
    i = i + 1
  }
  -1
}