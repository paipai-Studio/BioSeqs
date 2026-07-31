///|
/// Bio.SeqIO.UniprotIO - UniProt XML format parsing.
/// Inspired by Biopython's Bio.SeqIO.UniprotIO module.
/// Parses UniProt Knowledgebase XML format entries into SeqRecord structures.
///
/// UniProt XML format structure:
///   <uniprot>
///     <entry>
///       <accession>P12345</accession>
///       <protein><recommendedName><fullName>...</fullName></recommendedName></protein>
///       <gene><name>GENE1</name></gene>
///       <organism><name>Homo sapiens</name></organism>
///       <sequence length="123" mass="12345" modified="...">SEQUENCE</sequence>
///       <feature type="...">...</feature>
///       <reference>...</reference>
///     </entry>
///   </uniprot>

///|
/// A parsed UniProt XML entry with annotations.
pub struct UniprotEntry {
  accession : String
  protein_name : String
  gene_names : Array[String]
  organism : String
  sequence : String
  mut seq_length : Int
  mut seq_mass : Double
  ec_numbers : Array[String]
  mut subcellular_location : Array[String]
  mut keywords : Array[String]
  mut db_references : Array[DbReference]
  mut features : Array[UniprotFeature]
  mut references : Array[UniprotReference]
  comments : Array[String]
}

///|
/// A database cross-reference in a UniProt entry.
pub struct DbReference {
  database : String
  id : String
  properties : Map[String, String]
}

///|
/// A feature annotation in a UniProt entry.
pub struct UniprotFeature {
  feature_type : String
  description : String
  location : String
  qualifiers : Map[String, String]
}

///|
/// A reference in a UniProt entry.
pub struct UniprotReference {
  key : String
  title : String
  authors : Array[String]
  journal : String
  pubmed_id : String
  doi : String
}

///|
/// Create a new DbReference.
pub fn DbReference::new(
  database : String,
  id : String,
  properties? : Map[String, String] = Map([], capacity=8),
) -> DbReference {
  DbReference::{ database, id, properties }
}

///|
/// Create a new UniprotFeature.
pub fn UniprotFeature::new(
  feature_type : String,
  description : String,
  location : String,
  qualifiers? : Map[String, String] = Map([], capacity=8),
) -> UniprotFeature {
  UniprotFeature::{ feature_type, description, location, qualifiers }
}

///|
/// Create a new UniprotReference.
pub fn UniprotReference::new(
  key : String,
  title : String,
  authors : Array[String],
  journal : String,
  pubmed_id : String,
  doi : String,
) -> UniprotReference {
  UniprotReference::{ key, title, authors, journal, pubmed_id, doi }
}

///|
/// Create a new UniprotEntry.
pub fn UniprotEntry::new(
  accession : String,
  protein_name : String,
  gene_names : Array[String],
  organism : String,
  sequence : String,
) -> UniprotEntry {
  UniprotEntry::{
    accession,
    protein_name,
    gene_names,
    organism,
    sequence,
    seq_length: sequence.length(),
    seq_mass: 0.0,
    ec_numbers: Array::new(),
    subcellular_location: Array::new(),
    keywords: Array::new(),
    db_references: Array::new(),
    features: Array::new(),
    references: Array::new(),
    comments: Array::new(),
  }
}

///|
/// Parse a string to a float (Double).
fn u_parse_float(s : String) -> Double {
  if s.length() == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  let mut has_dot = false
  let mut decimal_place = 1.0
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      has_dot = true
      i = i + 1
      continue
    }
    if c >= zero && c <= zero + 9 {
      let digit = (c - zero).to_double()
      if !has_dot {
        result = result * 10.0 + digit
      } else {
        decimal_place = decimal_place / 10.0
        result = result + digit * decimal_place
      }
    }
    i = i + 1
  }
  result * sign
}

///|
/// Find the index of a substring in a string, returns -1 if not found.
fn u_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n {
    return -1
  }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Find the index of a substring from a starting position.
fn u_index_of_from(s : String, sub : String, start : Int) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || start > n - kn {
    return -1
  }
  let mut i = start
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Extract text content between XML tags.
fn extract_tag_content(xml : String, tag : String) -> String {
  let open_tag = "<" + tag + ">"
  let close_tag = "</" + tag + ">"
  let start = u_index_of(xml, open_tag)
  if start < 0 {
    return ""
  }
  let start_idx = start + open_tag.length()
  let end = u_index_of_from(xml, close_tag, start_idx)
  if end < 0 {
    return ""
  }
  let end_idx = end
  xml[start_idx:end_idx].trim().to_owned()
}

///|
/// Extract an attribute value from an opening tag string (e.g., "<dbReference type=\"PDB\" id=\"1ABC\">").
fn extract_attribute_from_str(tag_str : String, attr : String) -> String {
  let attr_pattern = attr + "=\""
  let attr_start = u_index_of(tag_str, attr_pattern)
  if attr_start < 0 {
    return ""
  }
  let as_idx = attr_start + attr_pattern.length()
  let quote_end = u_index_of_from(tag_str, "\"", as_idx)
  if quote_end < 0 {
    return ""
  }
  tag_str[as_idx:quote_end].to_owned()
}

///|
/// Extract text content between the first matching open/close tag pair.
fn extract_first_tag(xml : String, tag : String) -> String {
  let open = "<" + tag
  let start = u_index_of(xml, open)
  if start < 0 {
    return ""
  }
  let s = start
  let tag_end = u_index_of_from(xml, ">", s)
  if tag_end < 0 {
    return ""
  }
  let close_tag = "</" + tag + ">"
  let cstart = u_index_of_from(xml, close_tag, tag_end)
  if cstart < 0 {
    return ""
  }
  xml[tag_end + 1:cstart].trim().to_owned()
}

///|
/// Extract all occurrences of a tag's content.
fn extract_all_tags(xml : String, tag : String) -> Array[String] {
  let results : Array[String] = Array::new()
  let open_tag = "<" + tag
  let close_tag = "</" + tag + ">"
  let mut search_start = 0
  while true {
    let s = u_index_of_from(xml, open_tag, search_start)
    if s < 0 {
      break
    }
    let si = s
    let tag_end = u_index_of_from(xml, ">", si)
    if tag_end < 0 {
      break
    }
    let close = u_index_of_from(xml, close_tag, tag_end)
    if close < 0 {
      break
    }
    results.push(xml[tag_end + 1:close].trim().to_owned())
    search_start = close + close_tag.length()
  }
  results
}

///|
/// Extract all occurrences of a tag with their opening tag string.
/// Returns pairs of (opening_tag_string, inner_content).
fn extract_all_tags_with_attrs(xml : String, tag : String) -> Array[(String, String)] {
  let results : Array[(String, String)] = Array::new()
  let open_tag = "<" + tag
  let close_tag = "</" + tag + ">"
  let mut search_start = 0
  while true {
    let s = u_index_of_from(xml, open_tag, search_start)
    if s < 0 {
      break
    }
    let si = s
    let tag_end = u_index_of_from(xml, ">", si)
    if tag_end < 0 {
      break
    }
    let close = u_index_of_from(xml, close_tag, tag_end)
    if close < 0 {
      break
    }
    let open_str = xml[si:tag_end + 1].to_owned()
    let inner = xml[tag_end + 1:close].trim().to_owned()
    results.push((open_str, inner))
    search_start = close + close_tag.length()
  }
  results
}

///|
/// Extract all gene names from a UniProt entry XML.
fn extract_gene_names(entry_xml : String) -> Array[String] {
  let names : Array[String] = Array::new()
  let gene_blocks = extract_all_tags(entry_xml, "gene")
  for gb in gene_blocks {
    let names_in_gene = extract_all_tags(gb, "name")
    for n in names_in_gene {
      if n.length() > 0 {
        names.push(n)
      }
    }
    let orf_names = extract_all_tags(gb, "orfName")
    for n in orf_names {
      if n.length() > 0 {
        names.push(n)
      }
    }
  }
  names
}

///|
/// Extract organism name from entry XML.
fn extract_organism(entry_xml : String) -> String {
  let org_block = extract_first_tag(entry_xml, "organism")
  if org_block.length() == 0 {
    return ""
  }
  let names = extract_all_tags(org_block, "name")
  if names.length() > 0 {
    names[0]
  } else {
    ""
  }
}

///|
/// Extract protein name from entry XML.
fn extract_protein_name(entry_xml : String) -> String {
  let protein_block = extract_first_tag(entry_xml, "protein")
  if protein_block.length() == 0 {
    return ""
  }
  // Try recommendedName first
  let rec_block = extract_first_tag(protein_block, "recommendedName")
  if rec_block.length() > 0 {
    let full_name = extract_tag_content(rec_block, "fullName")
    if full_name.length() > 0 {
      return full_name
    }
  }
  // Try submittedName
  let sub_block = extract_first_tag(protein_block, "submittedName")
  if sub_block.length() > 0 {
    let full_name = extract_tag_content(sub_block, "fullName")
    if full_name.length() > 0 {
      return full_name
    }
  }
  ""
}

///|
/// Extract sequence and its attributes from entry XML.
fn extract_sequence(entry_xml : String) -> (String, Int, Double) {
  let seq_block = extract_first_tag(entry_xml, "sequence")
  if seq_block.length() == 0 {
    return ("", 0, 0.0)
  }
  let seq = seq_block.replace(old="\n", new="").replace(old=" ", new="").trim().to_owned()
  // Extract mass attribute
  let mass_attr = extract_attribute(entry_xml, "sequence", "mass")
  let mass = if mass_attr.length() > 0 { u_parse_float(mass_attr) } else { 0.0 }
  (seq, seq.length(), mass)
}

///|
/// Extract an XML attribute value from a tag.
fn extract_attribute(xml : String, tag : String, attr : String) -> String {
  let open = "<" + tag
  let start = u_index_of(xml, open)
  if start < 0 {
    return ""
  }
  let s = start
  let tag_end = u_index_of_from(xml, ">", s)
  if tag_end < 0 {
    return ""
  }
  let tag_str = xml[s:tag_end].to_owned()
  let attr_pattern = attr + "=\""
  let attr_start = u_index_of(tag_str, attr_pattern)
  if attr_start < 0 {
    return ""
  }
  let as_idx = attr_start + attr_pattern.length()
  let quote_end = u_index_of_from(tag_str, "\"", as_idx)
  if quote_end < 0 {
    return ""
  }
  tag_str[as_idx:quote_end].to_owned()
}

///|
/// Extract database cross-references from entry XML.
fn extract_db_references(entry_xml : String) -> Array[DbReference] {
  let refs : Array[DbReference] = Array::new()
  let blocks = extract_all_tags_with_attrs(entry_xml, "dbReference")
  let mut bi = 0
  while bi < blocks.length() {
    let block = blocks[bi]
    let open_str = block.0
    let inner = block.1
    let type_attr = extract_attribute_from_str(open_str, "type")
    let id_attr = extract_attribute_from_str(open_str, "id")
    if type_attr.length() > 0 {
      let props : Map[String, String] = Map([], capacity=8)
      // Try to extract properties
      let prop_blocks = extract_all_tags(inner, "property")
      for pb in prop_blocks {
        let ptype = extract_attribute(pb, "property", "type")
        let pval = extract_attribute(pb, "property", "value")
        if ptype.length() > 0 && pval.length() > 0 {
          props.set(ptype, pval)
        }
      }
      refs.push(DbReference::new(type_attr, id_attr, properties=props))
    }
    bi = bi + 1
  }
  refs
}

///|
/// Extract features from entry XML.
fn extract_features(entry_xml : String) -> Array[UniprotFeature] {
  let features : Array[UniprotFeature] = Array::new()
  let blocks = extract_all_tags_with_attrs(entry_xml, "feature")
  let mut bi = 0
  while bi < blocks.length() {
    let block = blocks[bi]
    let open_str = block.0
    let inner = block.1
    let ftype = extract_attribute_from_str(open_str, "type")
    let desc = extract_attribute_from_str(open_str, "description")
    let loc = extract_tag_content(inner, "location")
    let mut loc_str = ""
    if loc.length() > 0 {
      // Extract position info from location block
      let pos = extract_tag_content(loc, "position")
      if pos.length() > 0 {
        let pos_attr = extract_attribute(pos, "position", "value")
        loc_str = pos_attr
      } else {
        let begin = extract_tag_content(loc, "begin")
        let begin_val = extract_attribute(begin, "begin", "position")
        let end = extract_tag_content(loc, "end")
        let end_val = extract_attribute(end, "end", "position")
        if begin_val.length() > 0 && end_val.length() > 0 {
          loc_str = begin_val + "-" + end_val
        } else if begin_val.length() > 0 {
          loc_str = begin_val
        }
      }
    }
    if ftype.length() > 0 {
      features.push(UniprotFeature::new(ftype, desc, loc_str))
    }
    bi = bi + 1
  }
  features
}

///|
/// Extract references from entry XML.
fn extract_references(entry_xml : String) -> Array[UniprotReference] {
  let refs : Array[UniprotReference] = Array::new()
  let blocks = extract_all_tags(entry_xml, "reference")
  for b in blocks {
    let key = extract_attribute(b, "reference", "key")
    let citation = extract_first_tag(b, "citation")
    let title = extract_tag_content(citation, "title")
    let journal = extract_tag_content(citation, "journal")
    let authors : Array[String] = Array::new()
    let person_blocks = extract_all_tags(citation, "person")
    for pb in person_blocks {
      let person_name = extract_attribute(pb, "person", "name")
      if person_name.length() > 0 {
        authors.push(person_name)
      }
    }
    let pubmed_block = extract_first_tag(b, "dbReference")
    let pubmed_id = if pubmed_block.length() > 0 {
      let ptype = extract_attribute(pubmed_block, "dbReference", "type")
      if ptype == "PubMed" {
        extract_attribute(pubmed_block, "dbReference", "id")
      } else {
        ""
      }
    } else {
      ""
    }
    if key.length() > 0 {
      refs.push(UniprotReference::new(key, title, authors, journal, pubmed_id, ""))
    }
  }
  refs
}

///|
/// Extract keywords from entry XML.
fn extract_keywords(entry_xml : String) -> Array[String] {
  let words : Array[String] = Array::new()
  let blocks = extract_all_tags(entry_xml, "keyword")
  for b in blocks {
    let kw_id = extract_attribute(b, "keyword", "id")
    if kw_id.length() > 0 {
      words.push(kw_id)
    }
  }
  words
}

///|
/// Extract subcellular location from entry XML comments.
fn extract_subcellular(entry_xml : String) -> Array[String] {
  let locations : Array[String] = Array::new()
  let comment_blocks = extract_all_tags(entry_xml, "comment")
  for cb in comment_blocks {
    let comment_type = extract_attribute(cb, "comment", "type")
    if comment_type == "subcellular location" {
      let text_blocks = extract_all_tags(cb, "text")
      for tb in text_blocks {
        if tb.length() > 0 {
          locations.push(tb.trim().to_owned())
        }
      }
    }
  }
  locations
}

///|
/// Parse a single UniProt entry from XML string.
pub fn parse_uniprot_entry(entry_xml : String) -> UniprotEntry {
  // Extract accession
  let accession = extract_tag_content(entry_xml, "accession")
  // Extract protein name
  let protein_name = extract_protein_name(entry_xml)
  // Extract gene names
  let gene_names = extract_gene_names(entry_xml)
  // Extract organism
  let organism = extract_organism(entry_xml)
  // Extract sequence
  let (sequence, _, _) = extract_sequence(entry_xml)
  // Create base entry
  let entry = UniprotEntry::new(accession, protein_name, gene_names, organism, sequence)
  // Extract sequence with mass
  let (_, seq_len, seq_mass) = extract_sequence(entry_xml)
  entry.seq_length = seq_len
  entry.seq_mass = seq_mass
  // Extract dbReferences
  entry.db_references = extract_db_references(entry_xml)
  // Extract features
  entry.features = extract_features(entry_xml)
  // Extract references
  entry.references = extract_references(entry_xml)
  // Extract keywords
  entry.keywords = extract_keywords(entry_xml)
  // Extract subcellular location
  entry.subcellular_location = extract_subcellular(entry_xml)
  entry
}

///|
/// Parse UniProt XML content containing one or more entries.
pub fn parse_uniprot_xml(content : String) -> Array[UniprotEntry] {
  let entries : Array[UniprotEntry] = Array::new()
  let mut search_start = 0
  while true {
    let entry_start = u_index_of_from(content, "<entry>", search_start)
    if entry_start < 0 {
      break
    }
    let es = entry_start
    let entry_end = u_index_of_from(content, "</entry>", es)
    if entry_end < 0 {
      break
    }
    let entry_xml = content[es:entry_end + "</entry>".length()].to_owned()
    entries.push(parse_uniprot_entry(entry_xml))
    search_start = entry_end + "</entry>".length()
  }
  entries
}

///|
/// Convert a UniprotEntry to a SeqRecord-compatible map.
/// Returns a tuple of (id, description, sequence, annotations).
pub fn uniprot_to_seqrecord(entry : UniprotEntry) -> (String, String, String, Map[String, String]) {
  let id = entry.accession
  let description = entry.protein_name
  let annotations : Map[String, String] = Map([], capacity=16)
  annotations.set("gene_names", entry.gene_names.join(", "))
  annotations.set("organism", entry.organism)
  annotations.set("keywords", entry.keywords.join(", "))
  annotations.set("seq_length", entry.seq_length.to_string())
  annotations.set("seq_mass", entry.seq_mass.to_string())
  (id, description, entry.sequence, annotations)
}

///|
/// Create a sample UniProt XML entry for testing.
pub fn uniprot_sample_xml() -> String {
  "<uniprot>\n" +
  "  <entry>\n" +
  "    <accession>P12345</accession>\n" +
  "    <protein>\n" +
  "      <recommendedName>\n" +
  "        <fullName>Sample protein kinase</fullName>\n" +
  "      </recommendedName>\n" +
  "    </protein>\n" +
  "    <gene>\n" +
  "      <name>SampleGene</name>\n" +
  "      <orfName>SG1</orfName>\n" +
  "    </gene>\n" +
  "    <organism>\n" +
  "      <name>Homo sapiens</name>\n" +
  "      <name>human</name>\n" +
  "    </organism>\n" +
  "    <sequence length=\"31\" mass=\"3567.89\" modified=\"2023-01-01\">MSSSVSYTGN GQQVRVYLH EDPEPTDET AEVIQKLYH KKD</sequence>\n" +
  "    <keyword id=\"KW-0001\">Transferase</keyword>\n" +
  "    <keyword id=\"KW-0002\">Phosphoprotein</keyword>\n" +
  "    <comment type=\"subcellular location\">\n" +
  "      <text>Cytoplasm</text>\n" +
  "    </comment>\n" +
  "    <feature type=\"Chain\" description=\"Sample chain\">\n" +
  "      <location>\n" +
  "        <position value=\"1\"/>\n" +
  "      </location>\n" +
  "    </feature>\n" +
  "    <feature type=\"Region\" description=\"Catalytic region\">\n" +
  "      <location>\n" +
  "        <begin position=\"5\"/>\n" +
  "        <end position=\"25\"/>\n" +
  "      </location>\n" +
  "    </feature>\n" +
  "    <dbReference type=\"PDB\" id=\"1ABC\"/>\n" +
  "    <dbReference type=\"Ensembl\" id=\"ENSG00000000001\"/>\n" +
  "    <reference key=\"1\">\n" +
  "      <citation>\n" +
  "        <title>Sample paper on protein function</title>\n" +
  "        <journal>J. Biol. Chem.</journal>\n" +
  "        <person name=\"Smith J.\"/>\n" +
  "        <person name=\"Doe A.\"/>\n" +
  "      </citation>\n" +
  "      <dbReference type=\"PubMed\" id=\"12345678\"/>\n" +
  "    </reference>\n" +
  "  </entry>\n" +
  "</uniprot>"
}
