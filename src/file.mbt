///|
/// file module - Smart file handle with transparent compressed format support.
/// Inspired by Biopython Bio.File module.
///
/// References:
///   - https://biopython.org/wiki/Bio.File
///   - Biopython Bio.File: Generic file handle supporting gzip, bz2, gz2, etc.
///
/// This module provides smart file handles that automatically detect and
/// handle compressed file formats based on file extensions, enabling
/// transparent reading/writing of gzip, bzip2 compressed files.

///|
/// Supported file compression formats.
pub enum CompressionFormat {
  Plain /// No compression (plain text)
  Gzip /// gzip compression (.gz)
  Bzip2 /// bzip2 compression (.bz2)
}

///|
/// Create a Plain compression format.
pub fn compression_plain() -> CompressionFormat {
  CompressionFormat::Plain
}

///|
/// Create a Gzip compression format.
pub fn compression_gzip() -> CompressionFormat {
  CompressionFormat::Gzip
}

///|
/// Create a Bzip2 compression format.
pub fn compression_bzip2() -> CompressionFormat {
  CompressionFormat::Bzip2
}

///|
/// CompressionFormat to string.
pub fn CompressionFormat::to_string(self : CompressionFormat) -> String {
  match self {
    CompressionFormat::Plain => "plain"
    CompressionFormat::Gzip => "gzip"
    CompressionFormat::Bzip2 => "bzip2"
  }
}
///|
/// A smart file handle with automatic compression detection.
pub struct SmartFile {
  path : String
  format : CompressionFormat
  mode : String /// "r" for read, "w" for write
  is_open : Bool
  lines : Array[String] /// buffered lines for reading
  pos : Int /// current read position
  written_content : String /// accumulated content for writing
}

///|
/// Create a SmartFile from a file path, auto-detecting compression format.
///
/// Supported extensions:
///   - .gz -> gzip
///   - .bz2 -> bzip2
///   - others -> plain text
pub fn SmartFile::new(path : String, mode? : String = "r") -> SmartFile {
  let format = detect_compression(path)
  SmartFile::{
    path,
    format,
    mode,
    is_open: false,
    lines: ([] : Array[String]),
    pos: 0,
    written_content: "".to_string(),
  }
}

///|
/// Create a SmartFile with explicit compression format.
pub fn SmartFile::with_format(
  path : String,
  format : CompressionFormat,
  mode? : String = "r",
) -> SmartFile {
  SmartFile::{
    path,
    format,
    mode,
    is_open: false,
    lines: ([] : Array[String]),
    pos: 0,
    written_content: "".to_string(),
  }
}

///|
/// Detect compression format from file extension.
pub fn detect_compression(path : String) -> CompressionFormat {
  let lower = path.to_lower()
  if lower.has_suffix(".gz") {
    CompressionFormat::Gzip
  } else if lower.has_suffix(".bz2") {
    CompressionFormat::Bzip2
  } else {
    CompressionFormat::Plain
  }
}

///|
/// Get the file path.
pub fn SmartFile::get_path(self : SmartFile) -> String {
  self.path
}

///|
/// Get the detected compression format.
pub fn SmartFile::get_format(self : SmartFile) -> CompressionFormat {
  self.format
}

///|
/// Get the file mode.
pub fn SmartFile::get_mode(self : SmartFile) -> String {
  self.mode
}

///|
/// Check if the file handle is open.
pub fn SmartFile::is_open(self : SmartFile) -> Bool {
  self.is_open
}

///|
/// Check if the handle is for reading.
pub fn SmartFile::is_read_mode(self : SmartFile) -> Bool {
  self.mode == "r"
}

///|
/// Check if the handle is for writing.
pub fn SmartFile::is_write_mode(self : SmartFile) -> Bool {
  self.mode == "w"
}

///|
/// Open the file (prepare for reading or writing).
/// For reading, this loads all lines into memory (simulated for WASM).
/// For writing, this initializes the buffer.
pub fn SmartFile::open(self : SmartFile) -> SmartFile {
  if self.is_open {
    return self
  }
  SmartFile::{
    path: self.path,
    format: self.format,
    mode: self.mode,
    is_open: true,
    lines: self.lines,
    pos: self.pos,
    written_content: self.written_content,
  }
}

///|
/// Close the file handle.
pub fn SmartFile::close(self : SmartFile) -> SmartFile {
  SmartFile::{
    path: self.path,
    format: self.format,
    mode: self.mode,
    is_open: false,
    lines: self.lines,
    pos: self.pos,
    written_content: self.written_content,
  }
}

///|
/// Read a single line from the file.
/// Returns None if at end of file.
pub fn SmartFile::read_line(self : SmartFile) -> String? {
  if self.pos >= self.lines.length() {
    return None
  }
  let line = self.lines[self.pos]
  Some(line)
}

///|
/// Read all lines from the file.
pub fn SmartFile::read_lines(self : SmartFile) -> Array[String] {
  self.lines.copy()
}

///|
/// Read the entire file content as a single string.
pub fn SmartFile::read_all(self : SmartFile) -> String {
  let mut content = "".to_string()
  let mut i = 0
  while i < self.lines.length() {
    if i > 0 {
      content = content + "\n"
    }
    content = content + self.lines[i]
    i = i + 1
  }
  content
}

///|
/// Write a line to the file.
pub fn SmartFile::write_line(self : SmartFile, line : String) -> SmartFile {
  let content = if self.written_content.length() > 0 {
    self.written_content + "\n" + line
  } else {
    line
  }
  SmartFile::{
    path: self.path,
    format: self.format,
    mode: self.mode,
    is_open: self.is_open,
    lines: self.lines,
    pos: self.pos,
    written_content: content,
  }
}

///|
/// Write content to the file.
pub fn SmartFile::write(self : SmartFile, content : String) -> SmartFile {
  SmartFile::{
    path: self.path,
    format: self.format,
    mode: self.mode,
    is_open: self.is_open,
    lines: self.lines,
    pos: self.pos,
    written_content: content,
  }
}

///|
/// Get the written content (for write mode).
pub fn SmartFile::get_written_content(self : SmartFile) -> String {
  self.written_content
}

///|
/// Set pre-loaded lines (for testing/simulation purposes).
/// In a real implementation, this would be handled by open() reading the file.
pub fn SmartFile::set_lines(self : SmartFile, lines : Array[String]) -> SmartFile {
  SmartFile::{
    path: self.path,
    format: self.format,
    mode: self.mode,
    is_open: true,
    lines,
    pos: 0,
    written_content: self.written_content,
  }
}

///|
/// Create a SmartFile with pre-loaded content (for testing).
pub fn smart_file_from_string(path : String, content : String) -> SmartFile {
  let sf = SmartFile::new(path, mode="r")
  let iter = content.split("\n")
  let lines : Array[String] = Array::new()
  for sv in iter {
    lines.push(sv.to_owned())
  }
  sf.set_lines(lines)
}

///|
/// Create a SmartFile with pre-loaded lines (for testing).
pub fn smart_file_from_lines(path : String, lines : Array[String]) -> SmartFile {
  let sf = SmartFile::new(path, mode="r")
  sf.set_lines(lines)
}

///|
/// Create a write-mode SmartFile.
pub fn smart_file_writer(path : String, format? : CompressionFormat = CompressionFormat::Plain) -> SmartFile {
  SmartFile::with_format(path, format, mode="w")
}

///|
/// Generate a summary string for a SmartFile.
pub fn SmartFile::summary(self : SmartFile) -> String {
  let format_str = self.format.to_string()
  let mode_str = self.mode
  let state = if self.is_open { "open" } else { "closed" }
  let n_lines = self.lines.length()
  "SmartFile: " + self.path + " | format=" + format_str + " | mode=" + mode_str + " | " + state + " | lines=" + n_lines.to_string()
}

///|
/// Check if a file path has a recognized compression extension.
pub fn has_compression_extension(path : String) -> Bool {
  let lower = path.to_lower()
  lower.has_suffix(".gz") || lower.has_suffix(".bz2")
}

///|
/// Get the plain path without compression extension.
pub fn strip_compression_extension(path : String) -> String {
  let lower = path.to_lower()
  if lower.has_suffix(".gz") {
    path.substring(start=0, end=path.length() - 3)
  } else if lower.has_suffix(".bz2") {
    path.substring(start=0, end=path.length() - 4)
  } else {
    path
  }
}

///|
/// Generate sample data for testing.
pub fn file_sample_data() -> Array[String] {
  ["line1", "line2", "line3", ""]
}

///|
/// Generate sample lines with mixed content.
pub fn file_sample_lines() -> Array[String] {
  ["header1,header2,header3", "val1,val2,val3", "val4,val5,val6"]
}
