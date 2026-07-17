///|
/// Medline/PubMed record parsing, inspired by Biopython Bio.Medline.

///|
pub struct MedlineRecord {
  pmid : String
  title : String
  abstract_text : String
  authors : Array[String]
  journal : String
  volume : String
  issue : String
  pages : String
  publication_date : String
  mesh_terms : Array[String]
  keywords : Array[String]
  affiliations : Array[String]
  language : String
  pub_types : Array[String]
  doi : String
  pmc : String
}

///|
pub fn MedlineRecord::new(
  pmid : String,
  title : String,
  abstract_text : String,
  authors : Array[String],
  journal : String,
  volume : String,
  issue : String,
  pages : String,
  publication_date : String,
  mesh_terms : Array[String],
  keywords : Array[String],
  affiliations : Array[String],
  language : String,
  pub_types : Array[String],
  doi : String,
  pmc : String,
) -> MedlineRecord {
  MedlineRecord::{
    pmid,
    title,
    abstract_text,
    authors,
    journal,
    volume,
    issue,
    pages,
    publication_date,
    mesh_terms,
    keywords,
    affiliations,
    language,
    pub_types,
    doi,
    pmc,
  }
}

///|
pub fn parse_medline(text : String) -> MedlineRecord? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut pmid = ""
  let mut title = ""
  let mut abstract_text = ""
  let authors : Array[String] = Array::new()
  let mut journal = ""
  let mut volume = ""
  let mut issue = ""
  let mut pages = ""
  let mut publication_date = ""
  let mesh_terms : Array[String] = Array::new()
  let keywords : Array[String] = Array::new()
  let affiliations : Array[String] = Array::new()
  let mut language = ""
  let pub_types : Array[String] = Array::new()
  let mut doi = ""
  let mut pmc = ""

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()

    if line.length() >= 4 {
      let tag = line[0:4].trim().to_owned()
      let content_start = if line.length() >= 5 &&
        line.unsafe_get(4).to_int() == '-'.to_int() {
        6
      } else {
        4
      }
      let content = if line.length() > content_start {
        line[content_start:].trim().to_owned()
      } else {
        ""
      }

      if tag == "PMID" {
        pmid = content
      } else if tag == "TI" {
        title = content
      } else if tag == "AB" {
        abstract_text = content
      } else if tag == "AU" {
        authors.push(content)
      } else if tag == "JT" {
        journal = content
      } else if tag == "VI" {
        volume = content
      } else if tag == "IP" {
        issue = content
      } else if tag == "PG" {
        pages = content
      } else if tag == "DP" {
        publication_date = content
      } else if tag == "MH" {
        mesh_terms.push(content)
      } else if tag == "KW" {
        keywords.push(content)
      } else if tag == "AD" {
        affiliations.push(content)
      } else if tag == "LA" {
        language = content
      } else if tag == "PT" {
        pub_types.push(content)
      } else if tag == "AID" {
        if starts_with(content, 0, "doi:") {
          doi = content[4:].trim().to_owned()
        } else if starts_with(content, 0, "PMC") {
          pmc = content
        }
      }
    }

    i = i + 1
  }

  if pmid == "" {
    None
  } else {
    Some(
      MedlineRecord::new(
        pmid, title, abstract_text, authors, journal, volume, issue, pages, publication_date,
        mesh_terms, keywords, affiliations, language, pub_types, doi, pmc,
      ),
    )
  }
}

///|
pub fn parse_medline_records(text : String) -> Array[MedlineRecord] {
  let records : Array[MedlineRecord] = Array::new()
  let lines = split_by_char(text, '\n'.to_int())

  let mut current_block = ""
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()
    if line == "" {
      if current_block.length() > 0 && starts_with(current_block, 0, "PMID") {
        let parsed = parse_medline(current_block)
        if parsed.is_some() {
          records.push(parsed.unwrap())
        }
        current_block = ""
      }
    } else if current_block == "" {
      current_block = line
    } else {
      current_block = current_block + "\n" + line
    }
    i = i + 1
  }

  if current_block.length() > 0 && starts_with(current_block, 0, "PMID") {
    let parsed = parse_medline(current_block)
    if parsed.is_some() {
      records.push(parsed.unwrap())
    }
  }

  records
}

///|
pub fn format_citation_apa(record : MedlineRecord) -> String {
  let author_part = if record.authors.length() > 0 {
    let mut authors_str = ""
    let mut j = 0
    while j < record.authors.length() {
      if j > 0 {
        authors_str = authors_str + ", "
      }
      authors_str = authors_str + record.authors[j]
      j = j + 1
    }
    authors_str
  } else {
    "Anonymous"
  }

  let year = if record.publication_date.length() >= 4 {
    record.publication_date[0:4].to_owned()
  } else {
    ""
  }

  let journal_part = if record.journal != "" {
    " " + record.journal + "."
  } else {
    ""
  }

  let volume_part = if record.volume != "" { " " + record.volume } else { "" }

  let issue_part = if record.issue != "" {
    "(" + record.issue + ")"
  } else {
    ""
  }

  let pages_part = if record.pages != "" { ":" + record.pages } else { "" }

  author_part +
  " (" +
  year +
  "). " +
  record.title +
  "." +
  journal_part +
  volume_part +
  issue_part +
  pages_part
}

///|
pub fn get_first_author(record : MedlineRecord) -> String {
  if record.authors.length() > 0 {
    record.authors[0]
  } else {
    ""
  }
}

///|
pub fn has_abstract(record : MedlineRecord) -> Bool {
  record.abstract_text.length() > 0
}

///|
pub fn get_publication_year(record : MedlineRecord) -> String {
  if record.publication_date.length() >= 4 {
    record.publication_date[0:4].to_owned()
  } else {
    ""
  }
}

///|
pub fn get_author_count(record : MedlineRecord) -> Int {
  record.authors.length()
}

///|
pub fn filter_by_author(
  records : Array[MedlineRecord],
  author_name : String,
) -> Array[MedlineRecord] {
  let filtered : Array[MedlineRecord] = Array::new()
  let mut i = 0
  while i < records.length() {
    let mut j = 0
    while j < records[i].authors.length() {
      if contains(records[i].authors[j], author_name) {
        filtered.push(records[i])
        break
      }
      j = j + 1
    }
    i = i + 1
  }
  filtered
}

///|
pub fn extract_pmids(records : Array[MedlineRecord]) -> Array[String] {
  let pmids : Array[String] = Array::new()
  let mut i = 0
  while i < records.length() {
    pmids.push(records[i].pmid)
    i = i + 1
  }
  pmids
}

///|
pub fn filter_by_mesh(
  records : Array[MedlineRecord],
  mesh_term : String,
) -> Array[MedlineRecord] {
  let filtered : Array[MedlineRecord] = Array::new()
  let mut i = 0
  while i < records.length() {
    let mut j = 0
    while j < records[i].mesh_terms.length() {
      if contains(records[i].mesh_terms[j], mesh_term) {
        filtered.push(records[i])
        break
      }
      j = j + 1
    }
    i = i + 1
  }
  filtered
}

///|
fn contains(s : String, substr : String) -> Bool {
  if substr.length() == 0 {
    return true
  }
  let n = s.length()
  let m = substr.length()
  let mut i = 0
  while i <= n - m {
    let mut j = 0
    while j < m && s.unsafe_get(i + j) == substr.unsafe_get(j) {
      j = j + 1
    }
    if j == m {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn count_by_year(records : Array[MedlineRecord]) -> Array[(String, Int)] {
  let year_counts : Array[(String, Int)] = Array::new()

  let mut i = 0
  while i < records.length() {
    let year = get_publication_year(records[i])
    if year != "" {
      let mut found = false
      let mut j = 0
      while j < year_counts.length() {
        if year_counts[j].0 == year {
          year_counts[j] = (year, year_counts[j].1 + 1)
          found = true
          break
        }
        j = j + 1
      }
      if !found {
        year_counts.push((year, 1))
      }
    }
    i = i + 1
  }

  let mut j = 0
  while j < year_counts.length() - 1 {
    let mut k = j + 1
    while k < year_counts.length() {
      if year_counts[j].0 > year_counts[k].0 {
        let tmp = year_counts[j]
        year_counts[j] = year_counts[k]
        year_counts[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  year_counts
}
