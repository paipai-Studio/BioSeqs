/// Chromosome sizes utilities.
///
/// Bioframe equivalent: `fetch_chromsizes`, `make_chromsizes`
///
/// Build chromosome size tables from various sources.

///|
/// Chromosome size entry
pub struct BfChromSize {
  name : String
  length : Int32
}

///|
pub fn BfChromSize::name(self : BfChromSize) -> String {
  self.name
}

///|
pub fn BfChromSize::length(self : BfChromSize) -> Int32 {
  self.length
}

///|
/// Make chromsizes from an array of (name, length) pairs.
pub fn bf_make_chromsizes(sizes : Array[(String, Int32)]) -> Array[BfChromSize] {
  let result : Array[BfChromSize] = Array(capacity=0)
  for s in sizes {
    result.push({ name: s.0, length: s.1, })
  }
  result
}

///|
/// Fetch chromsizes from BAM header references.
/// Uses BamHeader reference info to build chromsizes.
pub fn bf_fetch_chromsizes_from_bam(
  refs : Array[(String, Int32)],
) -> Array[BfChromSize] {
  bf_make_chromsizes(refs)
}

///|
/// Split a line by tab and space, returning Array[String].
fn bf_split_fields(line : StringView) -> Array[String] {
  let fields : Array[String] = []
  let mut sb = StringBuilder()
  for c in line {
    if c == '\t' || c == ' ' {
      let cur = sb.to_string()
      if cur.length() > 0 {
        fields.push(cur)
      }
      sb = StringBuilder()
    } else {
      sb.write_char(c)
    }
  }
  let cur = sb.to_string()
  if cur.length() > 0 {
    fields.push(cur)
  }
  fields
}

///|
/// Fetch chromsizes from FASTA index (fai) content.
/// FAI format: name\tlength\toffset\tlinebases\tlinewidth
pub fn bf_fetch_chromsizes_from_fai(content : String) -> Array[BfChromSize] {
  let result : Array[BfChromSize] = Array(capacity=0)
  for line in content.split("\n") {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    let fields = bf_split_fields(trimmed)
    if fields.length() >= 2 {
      let name = fields[0]
      let length = parse_int(fields[1])
      result.push({ name, length, })
    }
  }
  result
}

///|
/// Fetch chromsizes from a chrom.sizes file content.
/// Format: name\tlength per line
pub fn bf_fetch_chromsizes_from_file(content : String) -> Array[BfChromSize] {
  let result : Array[BfChromSize] = Array(capacity=0)
  for line in content.split("\n") {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    let fields = bf_split_fields(trimmed)
    if fields.length() >= 2 {
      let name = fields[0]
      let length = parse_int(fields[1])
      result.push({ name, length, })
    }
  }
  result
}

///|
/// Convert chromsizes to (name, length) pairs.
pub fn bf_chromsizes_to_pairs(
  chromsizes : Array[BfChromSize],
) -> Array[(String, Int32)] {
  let result : Array[(String, Int32)] = Array(capacity=0)
  for cs in chromsizes {
    result.push((cs.name, cs.length))
  }
  result
}

///|
/// Get total genome size from chromsizes.
pub fn bf_genome_size(chromsizes : Array[BfChromSize]) -> Int64 {
  let mut total : Int64 = 0
  for cs in chromsizes {
    total = total + cs.length.to_int64()
  }
  total
}

///|
/// Sort chromsizes by chromosome name.
pub fn bf_sort_chromsizes(
  chromsizes : Array[BfChromSize],
) -> Array[BfChromSize] {
  let result = chromsizes.copy()
  result.sort_by(fn(a, b) { compare_chrom_names(a.name, b.name) })
  result
}

///|
fn compare_chrom_names(a : String, b : String) -> Int32 {
  if a < b {
    -1
  } else if a > b {
    1
  } else {
    0
  }
}

///|
/// Write chromsizes to chrom.sizes file format string.
pub fn bf_write_chromsizes(chromsizes : Array[BfChromSize]) -> String {
  let sb = StringBuilder()
  for cs in chromsizes {
    sb.write_string(cs.name)
    sb.write_char('\t')
    sb.write_string(cs.length.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}
