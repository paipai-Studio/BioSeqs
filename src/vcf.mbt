///|
/// VCF format parsing module, matching pysam VCF specification.

///|
/// VCF header line type.
pub enum VcfHeaderType {
  Meta
  Fixed
} derive(Eq)

///|
/// VCF meta header record (## lines).
pub struct VcfMetaHeader {
  key : String
  value : String
  tags : @hashmap.HashMap[String, String]
}

///|
/// VCF fixed header record (# line).
pub struct VcfFixedHeader {
  columns : Array[String]
}

///|
/// VCF record.
pub struct VcfRecord {
  chrom : String
  pos : Int
  id : Array[String]
  ref : String
  alt : Array[String]
  qual : Double?
  filter : Array[String]
  info : @hashmap.HashMap[String, String]
  format : Array[String]?
  samples : @hashmap.HashMap[String, Array[String]]
}

///|
/// VCF file representation.
pub struct VcfFile {
  meta_headers : Array[VcfMetaHeader]
  fixed_header : VcfFixedHeader?
  records : Array[VcfRecord]
}

///|
/// Parse INFO field into HashMap.
fn parse_info(info_str : String) -> @hashmap.HashMap[String, String] {
  let result = @hashmap.HashMap([])
  if info_str == "." || info_str.length() == 0 {
    return result
  }

  let entries = parse_semicolon_list(info_str)
  let mut i = 0
  while i < entries.length() {
    let entry = entries[i]
    let eq_pos = entry.find("=")
    if eq_pos.is_some() {
      let pos = eq_pos.unwrap()
      let key = entry[0:pos].to_owned()
      let value = if pos + 1 < entry.length() {
        entry[pos + 1:entry.length()].to_owned()
      } else {
        ""
      }
      result.set(key, value)
    } else {
      result.set(entry, "")
    }
    i = i + 1
  }

  result
}

///|
/// Parse FORMAT field into array.
fn parse_format(format_str : String) -> Array[String]? {
  if format_str == "." || format_str.length() == 0 {
    return None
  }
  Some(parse_comma_list(format_str))
}

///|
/// Parse VCF meta header line (##...).
fn parse_vcf_meta_header(line : String) -> VcfMetaHeader? {
  if line.length() < 3 ||
    line.unsafe_get(0).to_int() != '#'.to_int() ||
    line.unsafe_get(1).to_int() != '#'.to_int() {
    return None
  }

  let content = line[2:line.length()].to_owned()
  let eq_pos = content.find("=")
  if eq_pos.is_none() {
    return None
  }

  let pos = eq_pos.unwrap()
  let key = content[0:pos].to_owned()
  let value_part = content[pos + 1:content.length()].to_owned()

  let tags = @hashmap.HashMap([])
  let mut value = value_part

  if value_part.length() > 0 &&
    value_part.unsafe_get(0).to_int() == '<'.to_int() {
    let n = value_part.length()
    if n > 1 && value_part.unsafe_get(n - 1).to_int() == '>'.to_int() {
      let inner = value_part[1:n - 1].to_owned()
      let fields = parse_comma_list(inner)
      let mut i = 0
      while i < fields.length() {
        let field = fields[i]
        let feq_pos = field.find("=")
        if feq_pos.is_some() {
          let fpos = feq_pos.unwrap()
          let fkey = field[0:fpos].to_owned()
          let fvalue = if fpos + 1 < field.length() {
            field[fpos + 1:field.length()].to_owned()
          } else {
            ""
          }
          tags.set(fkey, fvalue)
        }
        i = i + 1
      }
      value = inner
    }
  }

  Some(VcfMetaHeader::{ key, value, tags })
}

///|
/// Parse VCF fixed header line (#...).
fn parse_vcf_fixed_header(line : String) -> VcfFixedHeader? {
  if line.length() < 2 ||
    line.unsafe_get(0).to_int() != '#'.to_int() ||
    line.unsafe_get(1).to_int() == '#'.to_int() {
    return None
  }

  let content = line[1:line.length()].to_owned()
  let columns = split_by_tab(content)
  Some(VcfFixedHeader::{ columns, })
}

///|
/// Parse VCF record line.
fn parse_vcf_record(line : String, sample_names : Array[String]) -> VcfRecord? {
  let parts = split_by_tab(line)
  if parts.length() < 8 {
    return None
  }

  let chrom = parts[0]
  let pos = parse_int(parts[1])
  let id = parse_comma_list(parts[2])
  let ref = parts[3]
  let alt = parse_comma_list(parts[4])
  let qual = parse_double(parts[5])
  let filter = parse_comma_list(parts[6])
  let info = parse_info(parts[7])

  let format = if parts.length() > 8 { parse_format(parts[8]) } else { None }

  let samples = @hashmap.HashMap([])
  if parts.length() > 9 && format.is_some() {
    let format_cols = format.unwrap()
    let mut i = 9
    while i < parts.length() {
      let sample_name = if i - 9 < sample_names.length() {
        sample_names[i - 9]
      } else {
        (i - 9).to_string()
      }
      let sample_data = parse_comma_list(parts[i])
      samples.set(sample_name, sample_data)
      i = i + 1
    }
  }

  Some(VcfRecord::{
    chrom,
    pos,
    id,
    ref,
    alt,
    qual,
    filter,
    info,
    format,
    samples,
  })
}

///|
/// Parse VCF content from string.
pub fn parse_vcf(content : String) -> VcfFile {
  let meta_headers = Array::new()
  let mut fixed_header : VcfFixedHeader? = None
  let records = Array::new()

  let lines = content.split("\n").to_array()
  let mut i = 0
  let n = lines.length()
  while i < n {
    let line = lines[i].to_owned()
    let line_len = line.length()
    if line_len == 0 {
      i = i + 1
      continue
    }

    if line.unsafe_get(0).to_int() == '#'.to_int() {
      if line.length() > 1 && line.unsafe_get(1).to_int() == '#'.to_int() {
        let meta_opt = parse_vcf_meta_header(line)
        if meta_opt.is_some() {
          meta_headers.push(meta_opt.unwrap())
        }
      } else {
        let fixed_opt = parse_vcf_fixed_header(line)
        if fixed_opt.is_some() {
          fixed_header = fixed_opt
        }
      }
    } else {
      let sample_names = match fixed_header {
        Some(h) => {
          let names = Array::new()
          let mut j = 9
          while j < h.columns.length() {
            names.push(h.columns[j])
            j = j + 1
          }
          names
        }
        None => []
      }
      let record_opt = parse_vcf_record(line, sample_names)
      if record_opt.is_some() {
        records.push(record_opt.unwrap())
      }
    }
    i = i + 1
  }

  VcfFile::{ meta_headers, fixed_header, records }
}

///|
/// Get number of records in VCF file.
pub fn VcfFile::len(self : VcfFile) -> Int {
  self.records.length()
}

///|
/// Get number of meta headers.
pub fn VcfFile::num_meta_headers(self : VcfFile) -> Int {
  self.meta_headers.length()
}

///|
/// Get sample names from fixed header.
pub fn VcfFile::get_sample_names(self : VcfFile) -> Array[String] {
  match self.fixed_header {
    Some(h) => {
      let result = Array::new()
      let mut i = 9
      while i < h.columns.length() {
        result.push(h.columns[i])
        i = i + 1
      }
      result
    }
    None => []
  }
}

///|
/// Get VCF version from meta headers.
pub fn VcfFile::get_version(self : VcfFile) -> String? {
  let mut i = 0
  while i < self.meta_headers.length() {
    let header = self.meta_headers[i]
    if header.key == "fileformat" {
      return Some(header.value)
    }
    i = i + 1
  }
  None
}

///|
/// Check if variant is a SNP (single nucleotide polymorphism).
pub fn VcfRecord::is_snp(self : VcfRecord) -> Bool {
  self.ref.length() == 1 && self.alt.length() == 1 && self.alt[0].length() == 1
}

///|
/// Check if variant is an insertion.
pub fn VcfRecord::is_insertion(self : VcfRecord) -> Bool {
  if self.alt.length() == 0 {
    return false
  }
  self.ref.length() < self.alt[0].length()
}

///|
/// Check if variant is a deletion.
pub fn VcfRecord::is_deletion(self : VcfRecord) -> Bool {
  if self.alt.length() == 0 {
    return false
  }
  self.ref.length() > self.alt[0].length()
}

///|
/// Get genotype for a sample.
pub fn VcfRecord::get_genotype(
  self : VcfRecord,
  sample_name : String,
) -> Array[String]? {
  match self.samples.get(sample_name) {
    Some(data) => Some(data.clone())
    None => None
  }
}
