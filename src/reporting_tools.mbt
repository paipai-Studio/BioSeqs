///|
/// ReportingTools module - Report generation for bioinformatics analysis results.
/// Inspired by Bioconductor ReportingTools package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/ReportingTools.html
///   - Hunnewell P, Jiao Y, Lapoint S, Davis S. (2024) ReportingTools: An interface for
///     generating detailed reports on DNA microarray analysis. R package version 2.34.0.

///|
/// Report section types.
pub enum ReportSectionType {
  Text
  Table
  Plot
  Section
  Header
}

///|
/// Create a Text section type.
pub fn report_section_text() -> ReportSectionType {
  ReportSectionType::Text
}

///|
/// Create a Table section type.
pub fn report_section_table() -> ReportSectionType {
  ReportSectionType::Table
}

///|
/// Create a Plot section type.
pub fn report_section_plot() -> ReportSectionType {
  ReportSectionType::Plot
}

///|
/// Create a Section section type.
pub fn report_section_section() -> ReportSectionType {
  ReportSectionType::Section
}

///|
/// Create a Header section type.
pub fn report_section_header() -> ReportSectionType {
  ReportSectionType::Header
}

///|
/// A single report section.
pub struct ReportSection {
  section_id : String
  section_type : ReportSectionType
  title : String
  content : String
  children : Array[ReportSection]
  level : Int
}

///|
/// A report table column.
pub struct ReportColumn {
  name : String
  values : Array[String]
  numeric : Bool
}

///|
/// A report table.
pub struct ReportTable {
  table_id : String
  title : String
  columns : Array[ReportColumn]
  caption : String
}

///|
/// A complete report document.
pub struct ReportDocument {
  title : String
  author : String
  sections : Array[ReportSection]
  tables : Array[ReportTable]
  timestamp : String
  version : String
}

///|
/// Section type to string.
pub fn ReportSectionType::to_string(self : ReportSectionType) -> String {
  match self {
    ReportSectionType::Text => "text"
    ReportSectionType::Table => "table"
    ReportSectionType::Plot => "plot"
    ReportSectionType::Section => "section"
    ReportSectionType::Header => "header"
  }
}

///|
/// Create a new report document.
pub fn ReportDocument::new(title : String) -> ReportDocument {
  ReportDocument::{
    title,
    author: "BioSeqs",
    sections: Array::new(),
    tables: Array::new(),
    timestamp: "",
    version: "1.0.0",
  }
}

///|
/// Set author.
pub fn ReportDocument::set_author(doc : ReportDocument, author : String) -> ReportDocument {
  ReportDocument::{
    title: doc.title,
    author,
    sections: doc.sections,
    tables: doc.tables,
    timestamp: doc.timestamp,
    version: doc.version,
  }
}

///|
/// Add a text section to the report.
pub fn ReportDocument::add_text(doc : ReportDocument, id : String, title : String, content : String) -> ReportDocument {
  let sections = doc.sections.copy()
  sections.push(ReportSection::{
    section_id: id,
    section_type: ReportSectionType::Text,
    title,
    content,
    children: Array::new(),
    level: 1,
  })
  ReportDocument::{
    title: doc.title,
    author: doc.author,
    sections,
    tables: doc.tables,
    timestamp: doc.timestamp,
    version: doc.version,
  }
}

///|
/// Add a table section to the report.
pub fn ReportDocument::add_table(doc : ReportDocument, id : String, table : ReportTable) -> ReportDocument {
  let sections = doc.sections.copy()
  let tables = doc.tables.copy()
  tables.push(table)
  sections.push(ReportSection::{
    section_id: id,
    section_type: ReportSectionType::Table,
    title: table.title,
    content: table.table_id,
    children: Array::new(),
    level: 1,
  })
  ReportDocument::{
    title: doc.title,
    author: doc.author,
    sections,
    tables,
    timestamp: doc.timestamp,
    version: doc.version,
  }
}

///|
/// Add a plot section to the report.
pub fn ReportDocument::add_plot(doc : ReportDocument, id : String, title : String, plot_ascii : String) -> ReportDocument {
  let sections = doc.sections.copy()
  sections.push(ReportSection::{
    section_id: id,
    section_type: ReportSectionType::Plot,
    title,
    content: plot_ascii,
    children: Array::new(),
    level: 1,
  })
  ReportDocument::{
    title: doc.title,
    author: doc.author,
    sections,
    tables: doc.tables,
    timestamp: doc.timestamp,
    version: doc.version,
  }
}

///|
/// Get the title.
pub fn ReportDocument::get_title(doc : ReportDocument) -> String {
  doc.title
}

///|
/// Get the author.
pub fn ReportDocument::get_author(doc : ReportDocument) -> String {
  doc.author
}

///|
/// Get the number of sections.
pub fn ReportDocument::get_n_sections(doc : ReportDocument) -> Int {
  doc.sections.length()
}

///|
/// Get the number of tables.
pub fn ReportDocument::get_n_tables(doc : ReportDocument) -> Int {
  doc.tables.length()
}

///|
/// Create a new report column.
pub fn ReportColumn::new(name : String, values : Array[String], numeric : Bool) -> ReportColumn {
  ReportColumn::{ name, values, numeric }
}

///|
/// Create a new report table from string arrays.
pub fn ReportTable::new(table_id : String, title : String, column_names : Array[String], column_data : Array[Array[String]], caption : String) -> ReportTable {
  let columns : Array[ReportColumn] = Array::new()
  let n = column_names.length()
  let mut i = 0
  while i < n {
    columns.push(ReportColumn::{
      name: column_names[i],
      values: column_data[i],
      numeric: false,
    })
    i = i + 1
  }
  ReportTable::{ table_id, title, columns, caption }
}

///|
/// Create a report table from mixed columns.
pub fn ReportTable::from_columns(table_id : String, title : String, columns : Array[ReportColumn], caption : String) -> ReportTable {
  ReportTable::{ table_id, title, columns, caption }
}

///|
/// Get column names.
pub fn ReportTable::get_column_names(table : ReportTable) -> Array[String] {
  let names : Array[String] = Array::new()
  let mut i = 0
  while i < table.columns.length() {
    names.push(table.columns[i].name)
    i = i + 1
  }
  names
}

///|
/// Get number of rows in a table.
pub fn ReportTable::get_n_rows(table : ReportTable) -> Int {
  if table.columns.length() == 0 { 0 } else { table.columns[0].values.length() }
}

///|
/// Get number of columns.
pub fn ReportTable::get_n_columns(table : ReportTable) -> Int {
  table.columns.length()
}

fn rt_pad_right(s : String, width : Int) -> String {
  if s.length() >= width {
    s.substring(start=0, end=width).to_string()
  } else {
    let mut result = s
    let mut i = s.length()
    while i < width {
      result = result + " "
      i = i + 1
    }
    result
  }
}

fn rt_repeat_char(c : String, n : Int) -> String {
  let mut result = ""
  let mut i = 0
  while i < n {
    result = result + c
    i = i + 1
  }
  result
}

///|
/// Format table as ASCII.
pub fn ReportTable::to_ascii(table : ReportTable) -> String {
  let names = table.get_column_names()
  let n_cols = names.length()
  let n_rows = table.get_n_rows()

  // Calculate column widths based on max of name length and data value lengths
  let widths : Array[Int] = Array::new()
  let mut ci = 0
  while ci < n_cols {
    let name_len = names[ci].length()
    let mut max_val_len = name_len
    if ci < table.columns.length() {
      let values = table.columns[ci].values
      let mut vi = 0
      while vi < values.length() {
        if values[vi].length() > max_val_len {
          max_val_len = values[vi].length()
        }
        vi = vi + 1
      }
    }
    let w = if max_val_len > 10 { max_val_len } else { 10 }
    widths.push(w)
    ci = ci + 1
  }

  let mut result = table.title + "\n"
  if table.caption.length() > 0 {
    result = result + table.caption + "\n"
  }
  result = result + "\n"

  // Header
  let mut hi = 0
  while hi < n_cols {
    let name = if hi < n_cols { names[hi] } else { "" }
    let padded = rt_pad_right(name, widths[hi])
    result = result + padded + " | "
    hi = hi + 1
  }
  result = result + "\n"

  // Separator
  let mut si = 0
  while si < n_cols {
    let dashes = rt_repeat_char("-", widths[si])
    result = result + dashes + "-+-"
    si = si + 1
  }
  result = result + "\n"

  // Rows
  let mut ri = 0
  while ri < n_rows {
    let mut ci2 = 0
    while ci2 < n_cols {
      let val = if ri < table.columns[ci2].values.length() { table.columns[ci2].values[ri] } else { "" }
      let padded = rt_pad_right(val, widths[ci2])
      result = result + padded + " | "
      ci2 = ci2 + 1
    }
    result = result + "\n"
    ri = ri + 1
  }

  result
}

///|
/// Get a report section by id.
pub fn ReportDocument::get_section(doc : ReportDocument, section_id : String) -> ReportSection {
  let mut i = 0
  while i < doc.sections.length() {
    if doc.sections[i].section_id == section_id {
      return doc.sections[i]
    }
    i = i + 1
  }
  // Return empty section if not found
  ReportSection::{
    section_id: "",
    section_type: ReportSectionType::Text,
    title: "",
    content: "Section not found",
    children: Array::new(),
    level: 0,
  }
}

///|
/// Render the full report as text/ASCII.
pub fn ReportDocument::render(doc : ReportDocument) -> String {
  let mut result = ""
  result = result + "========================================\n"
  result = result + doc.title + "\n"
  result = result + "Author: " + doc.author + "\n"
  result = result + "Version: " + doc.version + "\n"
  result = result + "Sections: " + doc.sections.length().to_string() + "\n"
  result = result + "========================================\n\n"

  let mut i = 0
  while i < doc.sections.length() {
    let sec = doc.sections[i]
    result = result + "### " + sec.title + " (" + sec.section_type.to_string() + ") ###\n"
    match sec.section_type {
      ReportSectionType::Text => {
        result = result + sec.content + "\n"
      }
      ReportSectionType::Table => {
        // Find the corresponding table
        let mut j = 0
        while j < doc.tables.length() {
          if doc.tables[j].table_id == sec.content {
            result = result + doc.tables[j].to_ascii() + "\n"
          }
          j = j + 1
        }
      }
      ReportSectionType::Plot => {
        result = result + sec.content + "\n"
      }
      _ => {
        result = result + sec.content + "\n"
      }
    }
    result = result + "\n---\n\n"
    i = i + 1
  }

  result
}

///|
/// Generate a report summary.
pub fn ReportDocument::summary(doc : ReportDocument) -> String {
  let mut s = "=== Report Summary ===\n".to_string()
  s = s + "Title: " + doc.title + "\n"
  s = s + "Author: " + doc.author + "\n"
  s = s + "Sections: " + doc.sections.length().to_string() + "\n"
  s = s + "Tables: " + doc.tables.length().to_string() + "\n"
  let mut total_rows = 0
  let mut i = 0
  while i < doc.tables.length() {
    total_rows = total_rows + doc.tables[i].get_n_rows()
    i = i + 1
  }
  s = s + "Total table rows: " + total_rows.to_string() + "\n"
  s
}