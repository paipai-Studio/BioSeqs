///|
pub struct BioMartDataset {
  name : String
  description : String
  host : String
  dataset : String
  mappings : Array[(String, String)]
  annotations : Array[(String, Array[(String, String)])]
}

///|
pub fn BioMartDataset::new(
  name : String,
  description : String,
  host : String,
  dataset : String,
) -> BioMartDataset {
  BioMartDataset::{
    name,
    description,
    host,
    dataset,
    mappings: Array::new(),
    annotations: Array::new(),
  }
}

///|
pub fn BioMartDataset::add_mapping(
  self : BioMartDataset,
  from_id : String,
  to_id : String,
) -> BioMartDataset {
  let new_mappings = self.mappings.copy()
  new_mappings.push((from_id, to_id))
  BioMartDataset::{
    name: self.name,
    description: self.description,
    host: self.host,
    dataset: self.dataset,
    mappings: new_mappings,
    annotations: self.annotations.copy(),
  }
}

///|
pub fn BioMartDataset::add_annotation(
  self : BioMartDataset,
  gene_id : String,
  annotations : Array[(String, String)],
) -> BioMartDataset {
  let new_annotations = self.annotations.copy()
  new_annotations.push((gene_id, annotations))
  BioMartDataset::{
    name: self.name,
    description: self.description,
    host: self.host,
    dataset: self.dataset,
    mappings: self.mappings.copy(),
    annotations: new_annotations,
  }
}

///|
pub fn BioMartDataset::get_mapping(
  self : BioMartDataset,
  from_id : String,
) -> String? {
  let mut i = 0
  while i < self.mappings.length() {
    if self.mappings[i].0 == from_id {
      return Some(self.mappings[i].1)
    }
    i = i + 1
  }
  None
}

///|
pub fn BioMartDataset::map_ids(
  self : BioMartDataset,
  ids : Array[String],
) -> Array[(String, String?)] {
  let results : Array[(String, String?)] = Array::new()
  let mut i = 0
  while i < ids.length() {
    results.push((ids[i], self.get_mapping(ids[i])))
    i = i + 1
  }
  results
}

///|
pub fn BioMartDataset::get_annotation(
  self : BioMartDataset,
  gene_id : String,
  attribute : String,
) -> String? {
  let mut i = 0
  while i < self.annotations.length() {
    if self.annotations[i].0 == gene_id {
      let annots = self.annotations[i].1
      let mut j = 0
      while j < annots.length() {
        if annots[j].0 == attribute {
          return Some(annots[j].1)
        }
        j = j + 1
      }
      return None
    }
    i = i + 1
  }
  None
}

///|
pub fn BioMartDataset::get_gene_name(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "gene_name")
}

///|
pub fn BioMartDataset::get_gene_description(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "description")
}

///|
pub fn BioMartDataset::get_chromosome(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "chromosome")
}

///|
pub fn BioMartDataset::get_start(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "start")
}

///|
pub fn BioMartDataset::get_end(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "end")
}

///|
pub fn BioMartDataset::get_strand(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "strand")
}

///|
pub fn BioMartDataset::get_biotype(
  self : BioMartDataset,
  gene_id : String,
) -> String? {
  self.get_annotation(gene_id, "biotype")
}

///|
pub fn BioMartDataset::get_external_db(
  self : BioMartDataset,
  gene_id : String,
  db_name : String,
) -> String? {
  self.get_annotation(gene_id, db_name)
}

///|
pub fn BioMartDataset::query(
  self : BioMartDataset,
  gene_ids : Array[String],
  attributes : Array[String],
) -> Array[Array[String]] {
  let results : Array[Array[String]] = Array::new()

  let header : Array[String] = Array::new()
  header.push("gene_id")
  let mut i = 0
  while i < attributes.length() {
    header.push(attributes[i])
    i = i + 1
  }
  results.push(header)

  let mut j = 0
  while j < gene_ids.length() {
    let row : Array[String] = Array::new()
    row.push(gene_ids[j])

    let mut k = 0
    while k < attributes.length() {
      let val = self.get_annotation(gene_ids[j], attributes[k])
      row.push(if val is Some(_) { val.unwrap() } else { "" })
      k = k + 1
    }
    results.push(row)
    j = j + 1
  }
  results
}

///|
pub fn create_human_mart() -> BioMartDataset {
  let mart = BioMartDataset::new(
    "hsapiens_gene_ensembl", "Human genes (GRCh38.p14)", "https://www.ensembl.org/biomart/martservice",
    "hsapiens_gene_ensembl",
  )

  let mart = mart.add_mapping("ENSG00000130203", "TP53")
  let mart = mart.add_mapping("ENSG00000141510", "BRCA1")
  let mart = mart.add_mapping("ENSG00000139618", "BRCA2")
  let mart = mart.add_mapping("ENSG00000160072", "EGFR")
  let mart = mart.add_mapping("ENSG00000157764", "KRAS")

  let mart = mart.add_annotation("ENSG00000130203", [
    ("gene_name", "TP53"),
    ("description", "Tumor protein p53"),
    ("chromosome", "17"),
    ("start", "7661779"),
    ("end", "7687550"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
    ("hgnc_symbol", "TP53"),
    ("uniprot", "P04637"),
  ])
  let mart = mart.add_annotation("ENSG00000141510", [
    ("gene_name", "BRCA1"),
    ("description", "BRCA1 DNA repair associated"),
    ("chromosome", "17"),
    ("start", "43044295"),
    ("end", "43125370"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
    ("hgnc_symbol", "BRCA1"),
    ("uniprot", "P38398"),
  ])
  let mart = mart.add_annotation("ENSG00000139618", [
    ("gene_name", "BRCA2"),
    ("description", "BRCA2 DNA repair associated"),
    ("chromosome", "13"),
    ("start", "32315086"),
    ("end", "32400268"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
    ("hgnc_symbol", "BRCA2"),
    ("uniprot", "P51587"),
  ])
  let mart = mart.add_annotation("ENSG00000160072", [
    ("gene_name", "EGFR"),
    ("description", "Epidermal growth factor receptor"),
    ("chromosome", "7"),
    ("start", "55019037"),
    ("end", "55201656"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
    ("hgnc_symbol", "EGFR"),
    ("uniprot", "P00533"),
  ])
  mart.add_annotation("ENSG00000157764", [
    ("gene_name", "KRAS"),
    ("description", "KRAS proto-oncogene, GTPase"),
    ("chromosome", "12"),
    ("start", "25245307"),
    ("end", "25254677"),
    ("strand", "1"),
    ("biotype", "protein_coding"),
    ("hgnc_symbol", "KRAS"),
    ("uniprot", "P01116"),
  ])
}

///|
pub fn create_mouse_mart() -> BioMartDataset {
  let mart = BioMartDataset::new(
    "mmusculus_gene_ensembl", "Mouse genes (GRCm39)", "https://www.ensembl.org/biomart/martservice",
    "mmusculus_gene_ensembl",
  )

  let mart = mart.add_mapping("ENSMUSG00000059552", "Trp53")
  let mart = mart.add_mapping("ENSMUSG00000015529", "Brca1")
  let mart = mart.add_mapping("ENSMUSG00000035740", "Egfr")

  let mart = mart.add_annotation("ENSMUSG00000059552", [
    ("gene_name", "Trp53"),
    ("description", "Tumor protein p53"),
    ("chromosome", "11"),
    ("start", "30968930"),
    ("end", "30998764"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
  ])
  let mart = mart.add_annotation("ENSMUSG00000015529", [
    ("gene_name", "Brca1"),
    ("description", "BRCA1 DNA repair associated"),
    ("chromosome", "11"),
    ("start", "98253072"),
    ("end", "98339596"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
  ])
  mart.add_annotation("ENSMUSG00000035740", [
    ("gene_name", "Egfr"),
    ("description", "Epidermal growth factor receptor"),
    ("chromosome", "10"),
    ("start", "126095626"),
    ("end", "126266787"),
    ("strand", "-1"),
    ("biotype", "protein_coding"),
  ])
}
