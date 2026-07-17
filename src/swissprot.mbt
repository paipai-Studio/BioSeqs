///|
/// SwissProt - Swiss-Prot/UniProt record parsing.
/// Inspired by Biopython Bio.SwissProt.
/// Parses UniProtKB/Swiss-Prot flat file format.

///|
/// A Swiss-Prot/UniProt record with annotations.
pub struct SwissProtRecord {
  entry_name : String
  accession_numbers : Array[String]
  protein_name : String
  gene_names : Array[String]
  organism : String
  organism_id : String
  taxonomy : Array[String]
  sequence : String
  seq_length : Int
  seq_mass : Double
  keywords : Array[String]
  function : String
  subcellular_location : Array[String]
  tissue_specificity : String
  ec_numbers : Array[String]
  domains : Array[String]
  references : Array[SwissProtReference]
  comments : Array[(String, String)]
  features : Array[SwissProtFeature]
  date_created : String
  date_modified : String
  version : Int
}

///|
/// A reference in a Swiss-Prot record.
pub struct SwissProtReference {
  number : Int
  title : String
  authors : Array[String]
  journal : String
  pubmed_id : String
  doi : String
}

///|
/// A feature annotation in a Swiss-Prot record.
pub struct SwissProtFeature {
  feature_type : String
  description : String
  start : Int
  end : Int
  qualifiers : Map[String, String]
}

///|
/// Create a new SwissProtReference.
pub fn SwissProtReference::new(
  number : Int,
  title : String,
  authors : Array[String],
  journal : String,
  pubmed_id : String,
  doi : String,
) -> SwissProtReference {
  SwissProtReference::{ number, title, authors, journal, pubmed_id, doi }
}

///|
/// Create a new SwissProtFeature.
pub fn SwissProtFeature::new(
  feature_type : String,
  description : String,
  start : Int,
  end : Int,
  qualifiers : Map[String, String],
) -> SwissProtFeature {
  SwissProtFeature::{ feature_type, description, start, end, qualifiers }
}

///|
/// Parse a Swiss-Prot/UniProt record from flat file text.
pub fn parse_swissprot(text : String) -> SwissProtRecord? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut entry_name = ""
  let accession_numbers : Array[String] = Array::new()
  let mut protein_name = ""
  let gene_names : Array[String] = Array::new()
  let mut organism = ""
  let mut organism_id = ""
  let taxonomy : Array[String] = Array::new()
  let mut sequence = ""
  let mut seq_length = 0
  let mut seq_mass = 0.0
  let keywords : Array[String] = Array::new()
  let mut function = ""
  let subcellular_location : Array[String] = Array::new()
  let mut tissue_specificity = ""
  let ec_numbers : Array[String] = Array::new()
  let domains : Array[String] = Array::new()
  let references : Array[SwissProtReference] = Array::new()
  let comments : Array[(String, String)] = Array::new()
  let features : Array[SwissProtFeature] = Array::new()
  let mut date_created = ""
  let mut date_modified = ""
  let mut version = 0

  let mut i = 0
  let mut in_sequence = false
  let mut current_ref_num = 0
  let mut current_ref_title = ""
  let current_ref_authors : Array[String] = Array::new()
  let mut current_ref_journal = ""
  let mut current_ref_pubmed = ""
  let mut current_ref_doi = ""
  let mut in_comment = false
  let mut comment_topic = ""
  let mut comment_text = ""

  while i < lines.length() {
    let line = lines[i]
    if line.length() < 2 {
      i = i + 1
      continue
    }

    let code = line[0:2].to_string()

    if code == "ID" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let parts = split_by_char(content, ' '.to_int())
      if parts.length() > 0 {
        entry_name = parts[0]
      }
    } else if code == "AC" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let acs = split_by_char(content, ';'.to_int())
      for ac in acs {
        let trimmed = ac.trim().to_string()
        if trimmed.length() > 0 {
          accession_numbers.push(trimmed)
        }
      }
    } else if code == "DE" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      let recname_pos = find_substring(content, "RecName: Full=")
      if recname_pos >= 0 {
        let start_pos = recname_pos + 13
        let end_pos = find_substring_from(content, ";", start_pos)
        if end_pos > start_pos {
          protein_name = content[start_pos:end_pos].to_string()
        }
      }
      let ec_pos = find_substring(content, "EC=")
      if ec_pos >= 0 {
        let ec_start = ec_pos + 3
        let ec_end = find_substring_from(content, ";", ec_start)
        if ec_end > ec_start {
          let ec = content[ec_start:ec_end].to_string()
          ec_numbers.push(ec)
        }
      }
    } else if code == "GN" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      let name_pos = find_substring(content, "Name=")
      if name_pos >= 0 {
        let name_start = name_pos + 5
        let name_end = find_substring_from(content, ";", name_start)
        if name_end > name_start {
          let gene = content[name_start:name_end].to_string()
          gene_names.push(gene)
        }
      }
    } else if code == "OS" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let paren_pos = find_substring(content, " (")
      let dot_pos = find_substring(content, ".")
      let end_pos = if paren_pos >= 0 {
        paren_pos
      } else if dot_pos >= 0 {
        dot_pos
      } else {
        content.length()
      }
      organism = content[:end_pos].trim().to_string()
    } else if code == "OX" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      let taxid_pos = find_substring(content, "NCBI_TaxID=")
      if taxid_pos >= 0 {
        let id_start = taxid_pos + 11
        let id_end = find_substring_from(content, ";", id_start)
        if id_end > id_start {
          organism_id = content[id_start:id_end].to_string()
        }
      }
    } else if code == "OC" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let parts = split_by_char(content, ';'.to_int())
      for part in parts {
        let trimmed = part.trim().to_string()
        if trimmed.length() > 0 {
          let dot_pos = find_substring(trimmed, ".")
          if dot_pos >= 0 {
            taxonomy.push(trimmed[:dot_pos].to_string())
          } else {
            taxonomy.push(trimmed)
          }
        }
      }
    } else if code == "KW" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let parts = split_by_char(content, ';'.to_int())
      for part in parts {
        let trimmed = part.trim().to_string()
        if trimmed.length() > 0 {
          let dot_pos = find_substring(trimmed, ".")
          if dot_pos >= 0 {
            keywords.push(trimmed[:dot_pos].to_string())
          } else {
            keywords.push(trimmed)
          }
        }
      }
    } else if code == "FT" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      if content.length() >= 30 {
        let ft_type = content[:15].trim().to_string()
        let ft_desc = content[15:30].trim().to_string()
        let ft_pos = content[30:].trim().to_string()
        if ft_type.length() > 0 && ft_pos.length() > 0 {
          let dotdot = find_substring(ft_pos, "..")
          if dotdot >= 0 {
            let start_str = ft_pos[:dotdot].to_string()
            let end_str = ft_pos[dotdot + 2:].to_string()
            let mut s = 0
            let mut e = 0
            s = parse_int(start_str)
            e = parse_int(end_str)
            let feat = SwissProtFeature::new(ft_type, ft_desc, s, e, Map([]))
            features.push(feat)
          }
        }
      }
    } else if code == "CC" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      if find_substring(content, "-!- FUNCTION:") >= 0 {
        in_comment = true
        comment_topic = "FUNCTION"
        comment_text = ""
      } else if find_substring(content, "-!- SUBCELLULAR LOCATION:") >= 0 {
        in_comment = true
        comment_topic = "SUBCELLULAR_LOCATION"
        comment_text = ""
      } else if find_substring(content, "-!- TISSUE SPECIFICITY:") >= 0 {
        in_comment = true
        comment_topic = "TISSUE_SPECIFICITY"
        comment_text = ""
      } else if find_substring(content, "-!-") >= 0 {
        if in_comment &&
          comment_topic == "FUNCTION" &&
          comment_text.length() > 0 {
          function = comment_text.trim().to_string()
        } else if in_comment &&
          comment_topic == "TISSUE_SPECIFICITY" &&
          comment_text.length() > 0 {
          tissue_specificity = comment_text.trim().to_string()
        }
        in_comment = false
      } else if in_comment {
        let trimmed = content.trim().to_string()
        if trimmed.length() > 0 {
          comment_text = comment_text + " " + trimmed
        }
        if comment_topic == "SUBCELLULAR_LOCATION" {
          let locs = split_by_char(trimmed, '.'.to_int())
          for loc in locs {
            let loc_trim = loc.trim().to_string()
            if loc_trim.length() > 0 {
              subcellular_location.push(loc_trim)
            }
          }
        }
      }
    } else if code == "RN" {
      if current_ref_num > 0 {
        let ref = SwissProtReference::new(
          current_ref_num,
          current_ref_title,
          copy_string_array(current_ref_authors),
          current_ref_journal,
          current_ref_pubmed,
          current_ref_doi,
        )
        references.push(ref)
      }
      current_ref_num = current_ref_num + 1
      current_ref_title = ""
      clear_string_array(current_ref_authors)
      current_ref_journal = ""
      current_ref_pubmed = ""
      current_ref_doi = ""
    } else if code == "RT" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      if content.length() > 0 {
        current_ref_title = current_ref_title + " " + content
      }
    } else if code == "RA" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      let authors_list = split_by_char(content, ','.to_int())
      for author in authors_list {
        let trimmed = author.trim().to_string()
        if trimmed.length() > 0 {
          current_ref_authors.push(trimmed)
        }
      }
    } else if code == "RL" {
      let content = if line.length() > 5 {
        line[5:].trim().to_string()
      } else {
        ""
      }
      current_ref_journal = current_ref_journal + " " + content
    } else if code == "RX" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      let pubmed_pos = find_substring(content, "PUBMED=")
      if pubmed_pos >= 0 {
        let pm_start = pubmed_pos + 7
        let pm_end = find_substring_from(content, ";", pm_start)
        if pm_end > pm_start {
          current_ref_pubmed = content[pm_start:pm_end].to_string()
        }
      }
      let doi_pos = find_substring(content, "DOI=")
      if doi_pos >= 0 {
        let doi_start = doi_pos + 4
        let doi_end = find_substring_from(content, ";", doi_start)
        if doi_end > doi_start {
          current_ref_doi = content[doi_start:doi_end].to_string()
        }
      }
    } else if code == "DT" {
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      if find_substring(content, "CREATED") >= 0 {
        date_created = content[:11].trim().to_string()
      } else if find_substring(content, "SEQUENCE UPDATED") >= 0 {
        date_modified = content[:11].trim().to_string()
      } else if find_substring(content, "VERSION") >= 0 {
        let v_pos = find_substring(content, "v")
        if v_pos >= 0 {
          let v_str = content[v_pos + 1:].trim().to_string()
          version = parse_int(v_str)
        }
      }
    } else if code == "SQ" {
      in_sequence = true
      let content = if line.length() > 5 { line[5:].to_string() } else { "" }
      let aa_pos = find_substring(content, "AA;")
      if aa_pos >= 0 {
        let before_aa = content[:aa_pos].trim().to_string()
        let space_pos = find_substring_rev(before_aa, " ")
        let num_str = if space_pos >= 0 {
          before_aa[space_pos + 1:].trim().to_string()
        } else {
          before_aa
        }
        seq_length = parse_int(num_str)
      }
      let mass_pos = find_substring(content, "Mass:")
      if mass_pos >= 0 {
        let mass_str = content[mass_pos + 5:].trim().to_string()
        let mass_num = split_by_char(mass_str, ' '.to_int())
        if mass_num.length() > 0 {
          match parse_double(mass_num[0]) {
            Some(v) => seq_mass = v
            None => seq_mass = 0.0
          }
        }
      }
    } else if code == "  " && in_sequence {
      let content = line.trim().to_string()
      let cleaned = content.replace(old=" ", new="")
      sequence = sequence + cleaned
    } else if code == "//" {
      if current_ref_num > 0 {
        let ref = SwissProtReference::new(
          current_ref_num,
          current_ref_title,
          copy_string_array(current_ref_authors),
          current_ref_journal,
          current_ref_pubmed,
          current_ref_doi,
        )
        references.push(ref)
      }
      if in_comment && comment_topic == "FUNCTION" && comment_text.length() > 0 {
        function = comment_text.trim().to_string()
      } else if in_comment &&
        comment_topic == "TISSUE_SPECIFICITY" &&
        comment_text.length() > 0 {
        tissue_specificity = comment_text.trim().to_string()
      }
      in_sequence = false
      in_comment = false
    }

    i = i + 1
  }

  if entry_name.length() == 0 && accession_numbers.length() == 0 {
    return None
  }

  Some(SwissProtRecord::{
    entry_name,
    accession_numbers,
    protein_name,
    gene_names,
    organism,
    organism_id,
    taxonomy,
    sequence,
    seq_length,
    seq_mass,
    keywords,
    function,
    subcellular_location,
    tissue_specificity,
    ec_numbers,
    domains,
    references,
    comments,
    features,
    date_created,
    date_modified,
    version,
  })
}

///|
/// Copy a string array.
fn copy_string_array(arr : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  for s in arr {
    result.push(s)
  }
  result
}

///|
/// Clear a string array.
fn clear_string_array(arr : Array[String]) -> Unit {
  while arr.length() > 0 {
    let _ = arr.pop()
  }
}

///|
/// Find substring position in a string.
fn find_substring(s : String, sub : String) -> Int {
  let n = s.length()
  let m = sub.length()
  if m == 0 || m > n {
    return -1
  }
  let mut i = 0
  while i <= n - m {
    let mut j = 0
    while j < m && s.unsafe_get(i + j) == sub.unsafe_get(j) {
      j = j + 1
    }
    if j == m {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Find substring position from a start index.
fn find_substring_from(s : String, sub : String, start : Int) -> Int {
  let n = s.length()
  let m = sub.length()
  if m == 0 || m > n || start >= n {
    return -1
  }
  let mut i = start
  while i <= n - m {
    let mut j = 0
    while j < m && s.unsafe_get(i + j) == sub.unsafe_get(j) {
      j = j + 1
    }
    if j == m {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Find last occurrence of substring.
fn find_substring_rev(s : String, sub : String) -> Int {
  let n = s.length()
  let m = sub.length()
  if m == 0 || m > n {
    return -1
  }
  let mut i = n - m
  while i >= 0 {
    let mut j = 0
    while j < m && s.unsafe_get(i + j) == sub.unsafe_get(j) {
      j = j + 1
    }
    if j == m {
      return i
    }
    i = i - 1
  }
  -1
}

///|
/// Get the primary accession number.
pub fn SwissProtRecord::primary_accession(self : SwissProtRecord) -> String? {
  if self.accession_numbers.length() > 0 {
    Some(self.accession_numbers[0])
  } else {
    None
  }
}

///|
/// Get the primary gene name.
pub fn SwissProtRecord::primary_gene(self : SwissProtRecord) -> String? {
  if self.gene_names.length() > 0 {
    Some(self.gene_names[0])
  } else {
    None
  }
}

///|
/// Get features of a specific type.
pub fn SwissProtRecord::get_features(
  self : SwissProtRecord,
  feature_type : String,
) -> Array[SwissProtFeature] {
  let result : Array[SwissProtFeature] = Array::new()
  for f in self.features {
    if f.feature_type == feature_type {
      result.push(f)
    }
  }
  result
}

///|
/// Create a sample/test SwissProt record for demonstration.
pub fn sample_swissprot_record() -> SwissProtRecord {
  let seq = "MSRSLLLRFLLFLLLLPPLPMLMIPGKPKLKVNLSLKDKISAFIRQLRQEGKPVLDFYLTGDQVEKVKGSFKEVKQELSKNPNAELLVLLENHLKPEDVSYLVGDSLLEILAKQGVNLLVGQPYEGKLLFTKKYIDSVLEKIKQKGVQVVGDLVLGHTLEALLRQANVTIMIGANPNSVGSLFTDFKNVTIGTDFQEIYPSQKTLYNKALNDYIAEFGGNVYSILTKSQELQKQFTENAKILMVGDVSFGHSSAFVKSLGVDIIVDIDVRLVGDSTVTEQFENFKRFKEVNPDDVLVVFAYDDGSLAQFVEKANVPVVQVDNGKLINQKIKQFLDQNNPNIYSVDELRILNKDVYIIRSGGKIDYSTFKELNQKLEKAFEKSIAPYSDPVKEQNIKQYFSYLQSKGVNQISFDDTKVKADFVIQPK"

  let features = [
    SwissProtFeature::new("DOMAIN", "GATase", 1, 180, Map([])),
    SwissProtFeature::new("DOMAIN", "ATPP synthetase", 181, 400, Map([])),
    SwissProtFeature::new("ACT_SITE", "Catalytic", 50, 50, Map([])),
    SwissProtFeature::new("BINDING", "ATP", 120, 125, Map([])),
  ]

  let refs = [
    SwissProtReference::new(
      1,
      "The complete amino acid sequence of human glutamine synthetase.",
      ["Smith A", "Johnson B", "Williams C"],
      "J. Biol. Chem. 265:1234-1240(1990)",
      "2315312",
      "10.1016/S0021-9258(19)38921-3",
    ),
  ]

  SwissProtRecord::{
    entry_name: "GLNA_HUMAN",
    accession_numbers: ["P15104", "Q9Y2K2"],
    protein_name: "Glutamine synthetase",
    gene_names: ["GLUL", "GLNS"],
    organism: "Homo sapiens",
    organism_id: "9606",
    taxonomy: [
      "Eukaryota", "Metazoa", "Chordata", "Craniata", "Vertebrata", "Euteleostomi",
      "Mammalia", "Eutheria", "Euarchontoglires", "Primates", "Haplorrhini", "Catarrhini",
      "Hominoidea", "Hominidae", "Homo",
    ],
    sequence: seq,
    seq_length: 373,
    seq_mass: 42128.5,
    keywords: [
      "ATP-binding", "Glutamine synthesis", "Mitochondrion", "Transferase",
    ],
    function: "Catalyzes the synthesis of glutamine from glutamate and ammonia in an ATP-dependent reaction.",
    subcellular_location: ["Cytoplasm", "Mitochondrion matrix"],
    tissue_specificity: "Expressed in brain, liver, and kidney.",
    ec_numbers: ["6.3.1.2"],
    domains: ["GATase", "ATPP synthetase"],
    references: refs,
    comments: [("FUNCTION", "Catalyzes glutamine synthesis")],
    features,
    date_created: "01-NOV-1986",
    date_modified: "15-JUL-2020",
    version: 188,
  }
}
