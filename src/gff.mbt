///|
/// Bio.GFF - GFF format parsing and gene annotation features.
/// Inspired by Biopython Bio.GFF.
/// Provides GFF3 format parsing and gene annotation feature extraction.

///|
pub struct GFFFeature {
  seqid : String
  source : String
  feature_type : String
  start : Int
  end : Int
  score : Double?
  strand : String
  phase : Int?
  attributes : @hashmap.HashMap[String, String]
}

///|
pub fn GFFFeature::new(
  seqid : String,
  source : String,
  feature_type : String,
  start : Int,
  end : Int,
  strand : String
) -> GFFFeature {
  GFFFeature::{
    seqid,
    source,
    feature_type,
    start,
    end,
    score: None,
    strand,
    phase: None,
    attributes: @hashmap.HashMap([])
  }
}

///|
pub fn GFFFeature::set_score(self : GFFFeature, score : Double) -> GFFFeature {
  GFFFeature::{
    seqid: self.seqid,
    source: self.source,
    feature_type: self.feature_type,
    start: self.start,
    end: self.end,
    score: Some(score),
    strand: self.strand,
    phase: self.phase,
    attributes: self.attributes
  }
}

///|
pub fn GFFFeature::set_phase(self : GFFFeature, phase : Int) -> GFFFeature {
  GFFFeature::{
    seqid: self.seqid,
    source: self.source,
    feature_type: self.feature_type,
    start: self.start,
    end: self.end,
    score: self.score,
    strand: self.strand,
    phase: Some(phase),
    attributes: self.attributes
  }
}

///|
pub fn GFFFeature::add_attribute(self : GFFFeature, key : String, value : String) -> GFFFeature {
  let new_attrs = self.attributes
  new_attrs.set(key, value)
  GFFFeature::{
    seqid: self.seqid,
    source: self.source,
    feature_type: self.feature_type,
    start: self.start,
    end: self.end,
    score: self.score,
    strand: self.strand,
    phase: self.phase,
    attributes: new_attrs
  }
}

///|
pub fn GFFFeature::get_attribute(self : GFFFeature, key : String) -> String? {
  if self.attributes.contains(key) {
    Some(self.attributes.get(key).unwrap())
  } else {
    None
  }
}

///|
pub fn GFFFeature::get_id(self : GFFFeature) -> String? {
  self.get_attribute("ID")
}

///|
pub fn GFFFeature::get_parent(self : GFFFeature) -> String? {
  self.get_attribute("Parent")
}

///|
pub fn GFFFeature::get_name(self : GFFFeature) -> String? {
  self.get_attribute("Name")
}

///|
pub fn GFFFeature::length(self : GFFFeature) -> Int {
  self.end - self.start + 1
}

///|
pub fn GFFFeature::is_coding(self : GFFFeature) -> Bool {
  self.feature_type == "CDS" || self.feature_type == "mRNA" || self.feature_type == "gene"
}

///|
pub fn GFFFeature::is_on_forward_strand(self : GFFFeature) -> Bool {
  self.strand == "+"
}

///|
pub fn GFFFeature::is_on_reverse_strand(self : GFFFeature) -> Bool {
  self.strand == "-"
}

///|
pub struct GFFRecord {
  features : Array[GFFFeature]
  sequences : @hashmap.HashMap[String, String]
  version : String
}

///|
pub fn GFFRecord::new() -> GFFRecord {
  GFFRecord::{
    features: Array::new(),
    sequences: @hashmap.HashMap([]),
    version: "3"
  }
}

///|
pub fn GFFRecord::add_feature(self : GFFRecord, feature : GFFFeature) -> GFFRecord {
  let new_features = self.features
  new_features.push(feature)
  GFFRecord::{
    features: new_features,
    sequences: self.sequences,
    version: self.version
  }
}

///|
pub fn GFFRecord::add_sequence(self : GFFRecord, seqid : String, sequence : String) -> GFFRecord {
  self.sequences.set(seqid, sequence)
  self
}

///|
pub fn GFFRecord::set_version(self : GFFRecord, version : String) -> GFFRecord {
  GFFRecord::{
    features: self.features,
    sequences: self.sequences,
    version: version
  }
}

///|
pub fn GFFRecord::get_sequence(self : GFFRecord, seqid : String) -> String? {
  if self.sequences.contains(seqid) {
    Some(self.sequences.get(seqid).unwrap())
  } else {
    None
  }
}

///|
pub fn GFFRecord::get_features_by_type(self : GFFRecord, feature_type : String) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    if self.features[i].feature_type == feature_type {
      result.push(self.features[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn GFFRecord::get_features_by_seqid(self : GFFRecord, seqid : String) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    if self.features[i].seqid == seqid {
      result.push(self.features[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn GFFRecord::get_features_by_id(self : GFFRecord, id : String) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    match self.features[i].get_id() {
      Some(feature_id) => {
        if feature_id == id {
          result.push(self.features[i])
        }
      }
      None => ()
    }
    i = i + 1
  }
  result
}

///|
pub fn GFFRecord::get_child_features(self : GFFRecord, parent_id : String) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    match self.features[i].get_parent() {
      Some(parent) => {
        let parents = split_by_char(parent, ','.to_int())
        let mut j = 0
        while j < parents.length() {
          let p = trim(parents[j])
          if p == parent_id {
            result.push(self.features[i])
            break
          }
          j = j + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  result
}

///|
pub fn GFFRecord::get_genes(self : GFFRecord) -> Array[GFFFeature] {
  self.get_features_by_type("gene")
}

///|
pub fn GFFRecord::get_mrnas(self : GFFRecord) -> Array[GFFFeature] {
  self.get_features_by_type("mRNA")
}

///|
pub fn GFFRecord::get_cds_features(self : GFFRecord) -> Array[GFFFeature] {
  self.get_features_by_type("CDS")
}

///|
pub fn GFFRecord::get_exons(self : GFFRecord) -> Array[GFFFeature] {
  self.get_features_by_type("exon")
}

///|
pub fn GFFRecord::get_transcripts_for_gene(self : GFFRecord, gene_id : String) -> Array[GFFFeature] {
  self.get_child_features(gene_id)
}

///|
pub fn GFFRecord::get_cds_for_transcript(self : GFFRecord, transcript_id : String) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let children = self.get_child_features(transcript_id)
  
  let mut i = 0
  while i < children.length() {
    if children[i].feature_type == "CDS" {
      result.push(children[i])
    } else if children[i].feature_type == "exon" {
      let exon_id = match children[i].get_id() {
        Some(id) => id
        None => ""
      }
      let exon_children = self.get_child_features(exon_id)
      let mut j = 0
      while j < exon_children.length() {
        if exon_children[j].feature_type == "CDS" {
          result.push(exon_children[j])
        }
        j = j + 1
      }
    }
    i = i + 1
  }
  
  result
}

///|
pub fn parse_attributes(attr_str : String) -> @hashmap.HashMap[String, String] {
  let attrs = @hashmap.HashMap([])
  let pairs = split_by_char(attr_str, ';'.to_int())
  
  let mut i = 0
  while i < pairs.length() {
    let pair = trim(pairs[i])
    if pair != "" {
      let parts = split_by_char(pair, '='.to_int())
      if parts.length() >= 2 {
        let key = trim(parts[0])
        let value = unescape_attribute_value(trim(parts[1]))
        attrs.set(key, value)
      }
    }
    i = i + 1
  }
  
  attrs
}

///|
fn unescape_attribute_value(value : String) -> String {
  let result = string_replace(value, "%3B", ";")
  let result2 = string_replace(result, "%3D", "=")
  let result3 = string_replace(result2, "%25", "%")
  let result4 = string_replace(result3, "&amp;", "&")
  let result5 = string_replace(result4, "&lt;", "<")
  string_replace(result5, "&gt;", ">")
}

///|
fn string_replace(s : String, from : String, to : String) -> String {
  let n = s.length()
  let m = from.length()
  if m == 0 || m > n {
    return s
  }
  
  let mut result = ""
  let mut i = 0
  while i < n {
    let match_found = starts_with(s, i, from)
    if match_found {
      result = result + to
      i = i + m
    } else {
      result = result + char_at(s, i)
      i = i + 1
    }
  }
  result
}

///|
pub fn bio_parse_gff(content : String) -> GFFRecord {
  let mut record = GFFRecord::new()
  
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line = trim(lines[i])
    
    if line == "" || starts_with(line, 0, "#") {
      if starts_with(line, 0, "##gff-version") {
        let parts = split_by_char(line, ' '.to_int())
        if parts.length() >= 2 {
          record = record.set_version(trim(parts[1]))
        }
      }
      i = i + 1
      continue
    }
    
    let fields = split_by_char(line, '\t'.to_int())
    if fields.length() >= 8 {
      let seqid = trim(fields[0])
      let source = trim(fields[1])
      let feature_type = trim(fields[2])
      let start = parse_int(trim(fields[3]))
      let end = parse_int(trim(fields[4]))
      let strand = trim(fields[6])
      
      let mut feature = GFFFeature::new(seqid, source, feature_type, start, end, strand)
      
      if fields.length() > 5 && trim(fields[5]) != "." {
        match parse_double(trim(fields[5])) {
          Some(score) => feature = feature.set_score(score)
          None => ()
        }
      }
      
      if fields.length() > 7 && trim(fields[7]) != "." {
        let phase = parse_int(trim(fields[7]))
        feature = feature.set_phase(phase)
      }
      
      if fields.length() > 8 {
        let attrs = parse_attributes(trim(fields[8]))
        let entries = attrs.iter().collect()
        let mut j = 0
        while j < entries.length() {
          let (key, value) = entries[j]
          feature = feature.add_attribute(key, value)
          j = j + 1
        }
      }
      
      ignore(record.add_feature(feature))
    }
    
    i = i + 1
  }
  
  record
}

///|
pub fn GFFRecord::to_gff(self : GFFRecord) -> String {
  let mut result = "##gff-version ".to_string() + self.version + "\n"
  
  let mut i = 0
  while i < self.features.length() {
    let f = self.features[i]
    result = result + f.seqid + "\t"
    result = result + f.source + "\t"
    result = result + f.feature_type + "\t"
    result = result + f.start.to_string() + "\t"
    result = result + f.end.to_string() + "\t"
    
    match f.score {
      Some(s) => result = result + s.to_string() + "\t"
      None => result = result + ".\t"
    }
    
    result = result + f.strand + "\t"
    
    match f.phase {
      Some(p) => result = result + p.to_string() + "\t"
      None => result = result + ".\t"
    }
    
    let attrs = f.attributes.iter().collect()
    let mut j = 0
    while j < attrs.length() {
      let (key, value) = attrs[j]
      if j > 0 {
        result = result + ";"
      }
      result = result + key + "=" + escape_attribute_value(value)
      j = j + 1
    }
    result = result + "\n"
    
    i = i + 1
  }
  
  result
}

///|
fn escape_attribute_value(value : String) -> String {
  let result = string_replace(value, "%", "%25")
  let result2 = string_replace(result, ";", "%3B")
  let result3 = string_replace(result2, "=", "%3D")
  let result4 = string_replace(result3, "&", "&amp;")
  let result5 = string_replace(result4, "<", "&lt;")
  string_replace(result5, ">", "&gt;")
}

///|
pub fn create_example_gff() -> GFFRecord {
  let record = GFFRecord::new()
  
  let gene = GFFFeature::new("chr1", "Ensembl", "gene", 10000, 15000, "+")
    .add_attribute("ID", "gene:ENSG00000130203")
    .add_attribute("Name", "TP53")
    .add_attribute("biotype", "protein_coding")
  ignore(record.add_feature(gene))
  
  let mrna = GFFFeature::new("chr1", "Ensembl", "mRNA", 10100, 14800, "+")
    .add_attribute("ID", "transcript:ENST00000335137")
    .add_attribute("Parent", "gene:ENSG00000130203")
    .add_attribute("biotype", "protein_coding")
    .add_attribute("transcript_id", "ENST00000335137")
  ignore(record.add_feature(mrna))
  
  let exon1 = GFFFeature::new("chr1", "Ensembl", "exon", 10100, 10200, "+")
    .add_attribute("ID", "exon:ENST00000335137:1")
    .add_attribute("Parent", "transcript:ENST00000335137")
  ignore(record.add_feature(exon1))
  
  let exon2 = GFFFeature::new("chr1", "Ensembl", "exon", 10500, 10600, "+")
    .add_attribute("ID", "exon:ENST00000335137:2")
    .add_attribute("Parent", "transcript:ENST00000335137")
  ignore(record.add_feature(exon2))
  
  let exon3 = GFFFeature::new("chr1", "Ensembl", "exon", 11000, 11500, "+")
    .add_attribute("ID", "exon:ENST00000335137:3")
    .add_attribute("Parent", "transcript:ENST00000335137")
  ignore(record.add_feature(exon3))
  
  let cds1 = GFFFeature::new("chr1", "Ensembl", "CDS", 10100, 10200, "+")
    .add_attribute("ID", "cds:ENST00000335137:1")
    .add_attribute("Parent", "transcript:ENST00000335137")
    .set_phase(0)
  ignore(record.add_feature(cds1))
  
  let cds2 = GFFFeature::new("chr1", "Ensembl", "CDS", 10500, 10600, "+")
    .add_attribute("ID", "cds:ENST00000335137:2")
    .add_attribute("Parent", "transcript:ENST00000335137")
    .set_phase(0)
  ignore(record.add_feature(cds2))
  
  let cds3 = GFFFeature::new("chr1", "Ensembl", "CDS", 11000, 11450, "+")
    .add_attribute("ID", "cds:ENST00000335137:3")
    .add_attribute("Parent", "transcript:ENST00000335137")
    .set_phase(0)
  ignore(record.add_feature(cds3))
  
  let gene2 = GFFFeature::new("chr17", "Ensembl", "gene", 43044295, 43125370, "-")
    .add_attribute("ID", "gene:ENSG00000141510")
    .add_attribute("Name", "BRCA1")
    .add_attribute("biotype", "protein_coding")
  ignore(record.add_feature(gene2))
  
  record
}

///|
pub fn GFFRecord::count_features(self : GFFRecord) -> Int {
  self.features.length()
}

///|
pub fn GFFRecord::get_unique_seqids(self : GFFRecord) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : @hashmap.HashMap[String, Bool] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < self.features.length() {
    let seqid = self.features[i].seqid
    if !seen.contains(seqid) {
      seen.set(seqid, true)
      result.push(seqid)
    }
    i = i + 1
  }
  
  result
}

///|
pub fn GFFRecord::get_unique_feature_types(self : GFFRecord) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : @hashmap.HashMap[String, Bool] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < self.features.length() {
    let ftype = self.features[i].feature_type
    if !seen.contains(ftype) {
      seen.set(ftype, true)
      result.push(ftype)
    }
    i = i + 1
  }
  
  result
}

///|
pub fn GFFRecord::get_features_in_range(self : GFFRecord, seqid : String, start : Int, end : Int) -> Array[GFFFeature] {
  let result : Array[GFFFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    let f = self.features[i]
    if f.seqid == seqid && f.start <= end && f.end >= start {
      result.push(f)
    }
    i = i + 1
  }
  result
}

///|
pub fn GFFRecord::calculate_gene_length(self : GFFRecord, gene_id : String) -> Int {
  let cds_features = self.get_cds_for_transcript(gene_id)
  let mut total_length = 0
  let mut i = 0
  while i < cds_features.length() {
    total_length = total_length + cds_features[i].length()
    i = i + 1
  }
  total_length
}