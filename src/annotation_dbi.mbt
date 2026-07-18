///|
pub struct AnnotationDb {
  gene_info : Map[String, GeneInfo]
  id_maps : Map[String, Map[String, String]]
  go_terms : Map[String, Array[AnnotationGOTerm]]
  kegg_terms : Map[String, Array[AnnotationKEGGTerm]]
  org_db : String
}

///|
pub struct GeneInfo {
  gene_id : String
  symbol : String
  description : String
  chr : String
  start : Int
  end : Int
  strand : String
  biotype : String
}

///|
pub fn GeneInfo::new(
  gene_id : String,
  symbol : String,
  description : String,
  chr : String,
  start : Int,
  end : Int,
  strand : String,
  biotype : String,
) -> GeneInfo {
  GeneInfo::{ gene_id, symbol, description, chr, start, end, strand, biotype }
}

///|
pub struct AnnotationGOTerm {
  go_id : String
  go_name : String
  evidence_code : String
  ontology : String
}

///|
pub fn AnnotationGOTerm::new(
  go_id : String,
  go_name : String,
  evidence_code : String,
  ontology : String,
) -> AnnotationGOTerm {
  AnnotationGOTerm::{ go_id, go_name, evidence_code, ontology }
}

///|
pub struct AnnotationKEGGTerm {
  kegg_id : String
  kegg_name : String
  pathway_type : String
}

///|
pub fn AnnotationKEGGTerm::new(
  kegg_id : String,
  kegg_name : String,
  pathway_type : String,
) -> AnnotationKEGGTerm {
  AnnotationKEGGTerm::{ kegg_id, kegg_name, pathway_type }
}

///|
pub fn create_annotation_db(org_db : String) -> AnnotationDb {
  AnnotationDb::{
    gene_info: Map([], capacity=1024),
    id_maps: Map([], capacity=16),
    go_terms: Map([], capacity=1024),
    kegg_terms: Map([], capacity=1024),
    org_db,
  }
}

///|
pub fn add_gene_info(db : AnnotationDb, gene_info : GeneInfo) -> Unit {
  db.gene_info[gene_info.gene_id] = gene_info
}

///|
pub fn add_id_map(
  db : AnnotationDb,
  from_type : String,
  mappings : Map[String, String],
) -> Unit {
  db.id_maps[from_type] = mappings
}

///|
pub fn add_go_terms(
  db : AnnotationDb,
  gene_id : String,
  terms : Array[AnnotationGOTerm],
) -> Unit {
  db.go_terms[gene_id] = terms
}

///|
pub fn add_kegg_terms(
  db : AnnotationDb,
  gene_id : String,
  terms : Array[AnnotationKEGGTerm],
) -> Unit {
  db.kegg_terms[gene_id] = terms
}

///|
pub fn annotation_dbi_get_gene_info(
  db : AnnotationDb,
  gene_id : String,
) -> GeneInfo? {
  if db.gene_info.contains(gene_id) {
    Some(db.gene_info[gene_id])
  } else {
    None
  }
}

///|
pub fn map_ids(
  db : AnnotationDb,
  ids : Array[String],
  from_type : String,
  to_type : String,
) -> Map[String, String] {
  let result = Map([], capacity=ids.length())

  if db.id_maps.contains(from_type) {
    let mapping = db.id_maps[from_type]
    let mut i = 0
    while i < ids.length() {
      let id = ids[i]
      if mapping.contains(id) {
        result[id] = mapping[id]
      }
      i = i + 1
    }
  }

  result
}

///|
pub fn select(
  db : AnnotationDb,
  keys : Array[String],
  columns : Array[String],
) -> Array[Map[String, String]] {
  let result : Array[Map[String, String]] = []

  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    if db.gene_info.contains(key) {
      let gene = db.gene_info[key]
      let row = Map([], capacity=columns.length())

      let mut j = 0
      while j < columns.length() {
        let col = columns[j]
        if col == "gene_id" {
          row[col] = gene.gene_id
        } else if col == "symbol" {
          row[col] = gene.symbol
        } else if col == "description" {
          row[col] = gene.description
        } else if col == "chr" {
          row[col] = gene.chr
        } else if col == "start" {
          row[col] = gene.start.to_string()
        } else if col == "end" {
          row[col] = gene.end.to_string()
        } else if col == "strand" {
          row[col] = gene.strand
        } else if col == "biotype" {
          row[col] = gene.biotype
        }
        j = j + 1
      }

      result.push(row)
    }
    i = i + 1
  }

  result
}

///|
pub fn get_go_annotations(
  db : AnnotationDb,
  gene_ids : Array[String],
) -> Map[String, Array[AnnotationGOTerm]] {
  let result = Map([], capacity=gene_ids.length())

  let mut i = 0
  while i < gene_ids.length() {
    let gene_id = gene_ids[i]
    if db.go_terms.contains(gene_id) {
      result[gene_id] = db.go_terms[gene_id]
    }
    i = i + 1
  }

  result
}

///|
pub fn get_kegg_annotations(
  db : AnnotationDb,
  gene_ids : Array[String],
) -> Map[String, Array[AnnotationKEGGTerm]] {
  let result = Map([], capacity=gene_ids.length())

  let mut i = 0
  while i < gene_ids.length() {
    let gene_id = gene_ids[i]
    if db.kegg_terms.contains(gene_id) {
      result[gene_id] = db.kegg_terms[gene_id]
    }
    i = i + 1
  }

  result
}

///|
pub fn annotation_dbi_keys(db : AnnotationDb) -> Array[String] {
  db.gene_info.keys().collect()
}

///|
pub fn annotation_dbi_columns(db : AnnotationDb) -> Array[String] {
  [
    "gene_id", "symbol", "description", "chr", "start", "end", "strand", "biotype",
  ]
}

///|
pub fn nrow(db : AnnotationDb) -> Int {
  db.gene_info.length()
}

///|
pub fn ncol(db : AnnotationDb) -> Int {
  8
}

///|
pub fn head(db : AnnotationDb, n : Int) -> Array[GeneInfo] {
  let result : Array[GeneInfo] = []
  let keys = db.gene_info.keys().collect()

  let mut i = 0
  while i < keys.length() && i < n {
    let key = keys[i]
    result.push(db.gene_info[key])
    i = i + 1
  }

  result
}

///|
pub fn summarize_db(db : AnnotationDb) -> String {
  let mut summary = "Annotation Database Summary:\n"
  summary = summary + "  Organism: " + db.org_db + "\n"
  summary = summary +
    "  Number of genes: " +
    db.gene_info.length().to_string() +
    "\n"
  summary = summary +
    "  Number of ID maps: " +
    db.id_maps.length().to_string() +
    "\n"
  summary = summary +
    "  Number of GO annotations: " +
    db.go_terms.length().to_string() +
    "\n"
  summary = summary +
    "  Number of KEGG annotations: " +
    db.kegg_terms.length().to_string() +
    "\n"
  summary
}
