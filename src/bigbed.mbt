// Bio.Align.bigbed - indexed UCSC BigBed records and interval queries.
//
// The implementation follows the BigBed version 4 binary layout used by
// Biopython: a chromosome B+ tree, binary BED data blocks, and an R-tree
// interval index. Writers use standards-compliant zlib stored blocks so the
// output remains deterministic without requiring a native compression library.

///|
pub suberror BigBedError {
  BigBedError(String)
}

///|
pub struct BigBedField {
  as_type : String
  name : String
  comment : String
} derive(Eq, Debug)

///|
pub struct BigBedSchema {
  name : String
  comment : String
  fields : Array[BigBedField]
} derive(Eq, Debug)

///|
pub struct BigBedTarget {
  name : String
  length : Int
} derive(Eq, Debug)

///|
pub struct BigBedRecord {
  chrom : String
  chrom_start : Int
  chrom_end : Int
  name : String
  score : Int
  strand : String
  thick_start : Int
  thick_end : Int
  item_rgb : String
  declared_block_count : Int
  block_sizes : Array[Int]
  block_starts : Array[Int]
  extra_fields : Array[String]
} derive(Eq, Debug)

///|
pub struct BigBedWriteConfig {
  bed_columns : Int
  items_per_slot : Int
  block_size : Int
  compress : Bool
  schema_name : String
  schema_comment : String
  custom_fields : Array[BigBedField]
} derive(Eq, Debug)

///|
pub struct BigBedHeader {
  little_endian : Bool
  version : Int
  zoom_levels : Int
  chromosome_tree_offset : Int
  full_data_offset : Int
  full_index_offset : Int
  field_count : Int
  defined_field_count : Int
  auto_sql_offset : Int
  total_summary_offset : Int
  uncompress_buffer_size : Int
  extra_indices_offset : Int
} derive(Eq, Debug)

///|
pub struct BigBedBlockIndex {
  start_chrom : Int
  start_base : Int
  end_chrom : Int
  end_base : Int
  data_offset : Int
  data_size : Int
} derive(Eq, Debug)

///|
pub struct BigBedSummary {
  record_count : Int
  target_count : Int
  data_block_count : Int
  covered_bases : Int
  minimum_span : Int
  maximum_span : Int
  compressed : Bool
} derive(Eq, Debug)

///|
pub struct BigBedFile {
  header : BigBedHeader
  schema : BigBedSchema
  targets : Array[BigBedTarget]
  source : Array[Int]
  blocks : Array[BigBedBlockIndex]
  record_count : Int
}

///|
enum BigBedEndian {
  BigBedLittleEndian
  BigBedBigEndian
} derive(Eq)

///|
let bigbed_magic : Int = 0x8789F2EB

///|
let bigbed_version : Int = 4

///|
let bigbed_chrom_tree_magic : Int = 0x78CA8C91

///|
let bigbed_rtree_magic : Int = 0x2468ACE0

///|
let bigbed_standard_names : Array[String] = [
  "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd",
  "reserved", "blockCount", "blockSizes", "chromStarts",
]

///|
let bigbed_standard_types : Array[String] = [
  "string", "uint", "uint", "string", "uint", "char[1]", "uint", "uint", "uint",
  "int", "int[blockCount]", "int[blockCount]",
]

///|
let bigbed_standard_comments : Array[String] = [
  "Reference sequence chromosome or scaffold", "Start position in chromosome", "End position in chromosome",
  "Name of item.", "Score (0-1000)", "+ or - for strand", "Start of where display should be thick",
  "End of where display should be thick", "Used as itemRgb", "Number of blocks",
  "Comma separated list of block sizes", "Start positions relative to chromStart",
]

///|
fn bigbed_fail(message : String) -> Unit raise BigBedError {
  raise BigBedError(message)
}

///|
fn bigbed_is_supported_type(as_type : String) -> Bool {
  let base = if contains(as_type, "[") {
    let mut end = 0
    while end < as_type.length() &&
          as_type.unsafe_get(end).to_int() != '['.to_int() {
      end = end + 1
    }
    as_type[0:end].to_owned()
  } else {
    as_type
  }
  base == "int" ||
  base == "uint" ||
  base == "short" ||
  base == "ushort" ||
  base == "byte" ||
  base == "ubyte" ||
  base == "float" ||
  base == "char" ||
  base == "string" ||
  base == "lstring"
}

///|
fn bigbed_validate_text(
  value : String,
  label : String,
) -> Unit raise BigBedError {
  for i in 0..<value.length() {
    let code = value.unsafe_get(i).to_int()
    if code == 0 ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      bigbed_fail(label + " contains a forbidden control character")
    }
  }
}

///|
pub fn BigBedField::create(
  as_type : String,
  name : String,
  comment? : String = "",
) -> BigBedField raise BigBedError {
  if !bigbed_is_supported_type(as_type) {
    bigbed_fail("Unsupported AutoSQL field type '" + as_type + "'")
  }
  if name.length() == 0 {
    bigbed_fail("AutoSQL field name must not be empty")
  }
  bigbed_validate_text(name, "AutoSQL field name")
  BigBedField::{ as_type, name, comment }
}

///|
fn bigbed_standard_fields(count : Int) -> Array[BigBedField] {
  let fields : Array[BigBedField] = []
  for i in 0..<count {
    fields.push(BigBedField::{
      as_type: bigbed_standard_types[i],
      name: bigbed_standard_names[i],
      comment: bigbed_standard_comments[i],
    })
  }
  fields
}

///|
pub fn BigBedSchema::create(
  name : String,
  comment : String,
  fields : Array[BigBedField],
) -> BigBedSchema raise BigBedError {
  if name.length() == 0 {
    bigbed_fail("AutoSQL table name must not be empty")
  }
  if fields.length() < 3 {
    bigbed_fail("AutoSQL schema must contain at least three fields")
  }
  let seen : Map[String, Bool] = Map([])
  for field in fields {
    if seen.contains(field.name) {
      bigbed_fail("Duplicate AutoSQL field '" + field.name + "'")
    }
    seen[field.name] = true
  }
  BigBedSchema::{ name, comment, fields: fields.copy() }
}

///|
pub fn BigBedSchema::default(
  bed_columns? : Int = 12,
) -> BigBedSchema raise BigBedError {
  if bed_columns < 3 || bed_columns > 12 {
    bigbed_fail("BED column count must be between 3 and 12")
  }
  BigBedSchema::{
    name: "bed",
    comment: "Browser Extensible Data",
    fields: bigbed_standard_fields(bed_columns),
  }
}

///|
pub fn BigBedTarget::create(
  name : String,
  length : Int,
) -> BigBedTarget raise BigBedError {
  if name.length() == 0 {
    bigbed_fail("BigBed target name must not be empty")
  }
  bigbed_validate_text(name, "BigBed target name")
  if length <= 0 {
    bigbed_fail("BigBed target length must be positive")
  }
  BigBedTarget::{ name, length }
}

///|
fn bigbed_validate_blocks(
  span : Int,
  sizes : Array[Int],
  starts : Array[Int],
) -> Unit raise BigBedError {
  if sizes.length() != starts.length() {
    bigbed_fail("BED block sizes and starts must have equal lengths")
  }
  if sizes.length() == 0 {
    return
  }
  if starts[0] != 0 {
    bigbed_fail("The first BED block must start at chromStart")
  }
  let mut previous_end = 0
  for i in 0..<sizes.length() {
    if sizes[i] <= 0 {
      bigbed_fail("BED block sizes must be positive")
    }
    if starts[i] < previous_end {
      bigbed_fail("BED blocks must be sorted and non-overlapping")
    }
    let block_end = starts[i] + sizes[i]
    if block_end > span {
      bigbed_fail("BED block extends beyond chromEnd")
    }
    previous_end = block_end
  }
  if previous_end != span {
    bigbed_fail("The last BED block must end at chromEnd")
  }
}

///|
pub fn BigBedRecord::create(
  chrom : String,
  chrom_start : Int,
  chrom_end : Int,
  name? : String = ".",
  score? : Int = 0,
  strand? : String = ".",
  thick_start? : Int = -1,
  thick_end? : Int = -1,
  item_rgb? : String = "0",
  block_count? : Int = -1,
  block_sizes? : Array[Int] = [],
  block_starts? : Array[Int] = [],
  extra_fields? : Array[String] = [],
) -> BigBedRecord raise BigBedError {
  if chrom.length() == 0 {
    bigbed_fail("BED chromosome must not be empty")
  }
  bigbed_validate_text(chrom, "BED chromosome")
  if chrom_start < 0 || chrom_end <= chrom_start {
    bigbed_fail("BED coordinates must define a positive half-open interval")
  }
  if score < 0 || score > 1000 {
    bigbed_fail("BED score must be between 0 and 1000")
  }
  if strand != "+" && strand != "-" && strand != "." {
    bigbed_fail("BED strand must be '+', '-', or '.'")
  }
  bigbed_validate_text(name, "BED name")
  bigbed_validate_text(item_rgb, "BED itemRgb")
  for value in extra_fields {
    bigbed_validate_text(value, "BigBed custom field")
  }
  let normalized_thick_start = if thick_start < 0 {
    chrom_start
  } else {
    thick_start
  }
  let normalized_thick_end = if thick_end < 0 { chrom_end } else { thick_end }
  if normalized_thick_start < chrom_start ||
    normalized_thick_start > normalized_thick_end ||
    normalized_thick_end > chrom_end {
    bigbed_fail("BED thick interval must lie inside the record interval")
  }
  let normalized_block_count = if block_count < 0 {
    block_sizes.length()
  } else {
    block_count
  }
  if normalized_block_count < 0 {
    bigbed_fail("BED block count must not be negative")
  }
  if block_sizes.length() > 0 && block_sizes.length() != normalized_block_count {
    bigbed_fail("BED block count does not match block sizes")
  }
  if block_starts.length() > 0 &&
    block_starts.length() != normalized_block_count {
    bigbed_fail("BED block count does not match block starts")
  }
  if block_starts.length() > 0 && block_sizes.length() == 0 {
    bigbed_fail("BED block starts require corresponding block sizes")
  }
  if block_sizes.length() > 0 && block_starts.length() > 0 {
    bigbed_validate_blocks(chrom_end - chrom_start, block_sizes, block_starts)
  }
  BigBedRecord::{
    chrom,
    chrom_start,
    chrom_end,
    name,
    score,
    strand,
    thick_start: normalized_thick_start,
    thick_end: normalized_thick_end,
    item_rgb,
    declared_block_count: normalized_block_count,
    block_sizes: block_sizes.copy(),
    block_starts: block_starts.copy(),
    extra_fields: extra_fields.copy(),
  }
}

///|
pub fn BigBedWriteConfig::create(
  bed_columns? : Int = 12,
  items_per_slot? : Int = 512,
  block_size? : Int = 256,
  compress? : Bool = true,
  schema_name? : String = "bed",
  schema_comment? : String = "Browser Extensible Data",
  custom_fields? : Array[BigBedField] = [],
) -> BigBedWriteConfig raise BigBedError {
  if bed_columns < 3 || bed_columns > 12 {
    bigbed_fail("BED column count must be between 3 and 12")
  }
  if items_per_slot <= 0 || items_per_slot > 65535 {
    bigbed_fail("BigBed itemsPerSlot must be between 1 and 65535")
  }
  if block_size < 2 || block_size > 65535 {
    bigbed_fail("BigBed blockSize must be between 2 and 65535")
  }
  let fields = bigbed_standard_fields(bed_columns)
  for custom in custom_fields {
    for standard in fields {
      if standard.name == custom.name {
        bigbed_fail(
          "Custom AutoSQL field duplicates standard field '" + custom.name + "'",
        )
      }
    }
  }
  ignore(
    BigBedSchema::create(
      schema_name,
      schema_comment,
      {
        let all = fields.copy()
        for field in custom_fields {
          all.push(field)
        }
        all
      },
    ),
  )
  BigBedWriteConfig::{
    bed_columns,
    items_per_slot,
    block_size,
    compress,
    schema_name,
    schema_comment,
    custom_fields: custom_fields.copy(),
  }
}

///|
pub fn BigBedWriteConfig::default() -> BigBedWriteConfig {
  BigBedWriteConfig::{
    bed_columns: 12,
    items_per_slot: 512,
    block_size: 256,
    compress: true,
    schema_name: "bed",
    schema_comment: "Browser Extensible Data",
    custom_fields: [],
  }
}

///|
pub fn BigBedRecord::span(self : BigBedRecord) -> Int {
  self.chrom_end - self.chrom_start
}

///|
pub fn BigBedRecord::block_count(self : BigBedRecord) -> Int {
  self.declared_block_count
}

///|
pub fn BigBedRecord::query_size(self : BigBedRecord) -> Int {
  if self.block_sizes.length() == 0 ||
    self.block_sizes.length() != self.block_starts.length() {
    self.span()
  } else {
    let mut total = 0
    for size in self.block_sizes {
      total = total + size
    }
    total
  }
}

///|
pub fn BigBedRecord::overlaps(
  self : BigBedRecord,
  chromosome : String,
  start : Int,
  end : Int,
) -> Bool {
  self.chrom == chromosome && self.chrom_start < end && start < self.chrom_end
}

///|
pub fn BigBedRecord::coordinates(self : BigBedRecord) -> Array[Array[Int]] {
  let target : Array[Int] = []
  let query : Array[Int] = []
  let query_size = self.query_size()
  if self.block_sizes.length() == 0 ||
    self.block_sizes.length() != self.block_starts.length() {
    target.push(self.chrom_start)
    target.push(self.chrom_end)
    query.push(0)
    query.push(query_size)
  } else {
    let mut query_position = 0
    for i in 0..<self.block_sizes.length() {
      let block_start = self.chrom_start + self.block_starts[i]
      let block_end = block_start + self.block_sizes[i]
      if i == 0 {
        target.push(block_start)
        query.push(query_position)
      } else {
        target.push(block_start)
        query.push(query_position)
      }
      query_position = query_position + self.block_sizes[i]
      target.push(block_end)
      query.push(query_position)
    }
  }
  if self.strand == "-" {
    for i in 0..<query.length() {
      query[i] = query_size - query[i]
    }
  }
  [target, query]
}

///|
pub fn BigBedRecord::annotation(
  self : BigBedRecord,
  schema : BigBedSchema,
  key : String,
  bed_columns : Int,
) -> String? {
  for i in 0..<self.extra_fields.length() {
    let field_index = bed_columns + i
    if field_index < schema.fields.length() &&
      schema.fields[field_index].name == key {
      return Some(self.extra_fields[i])
    }
  }
  None
}

///|
fn bigbed_escape_auto_sql(value : String) -> String {
  let builder = StringBuilder::new()
  for i in 0..<value.length() {
    let code = value.unsafe_get(i).to_int()
    if code == '"'.to_int() {
      builder.write_string("'")
    } else if code != '\n'.to_int() && code != '\r'.to_int() {
      builder.write_char(value.unsafe_get(i).unsafe_to_char())
    }
  }
  builder.to_string()
}

///|
pub fn BigBedSchema::to_auto_sql(self : BigBedSchema) -> String {
  let output = StringBuilder::new()
  output.write_string("table " + self.name + "\n")
  output.write_string("\"" + bigbed_escape_auto_sql(self.comment) + "\"\n")
  output.write_string("(\n")
  for field in self.fields {
    output.write_string(
      "    " +
      field.as_type +
      " " +
      field.name +
      "; \"" +
      bigbed_escape_auto_sql(field.comment) +
      "\"\n",
    )
  }
  output.write_string(")\n")
  output.to_string()
}

///|
fn bigbed_find_byte(text : String, code : Int, start : Int) -> Int {
  for i in start..<text.length() {
    if text.unsafe_get(i).to_int() == code {
      return i
    }
  }
  -1
}

///|
fn bigbed_quoted_value(text : String) -> String raise BigBedError {
  let first = bigbed_find_byte(text, '"'.to_int(), 0)
  if first < 0 {
    bigbed_fail("Expected quoted AutoSQL text")
  }
  let second = bigbed_find_byte(text, '"'.to_int(), first + 1)
  if second < 0 {
    bigbed_fail("Unterminated quoted AutoSQL text")
  }
  text[first + 1:second].to_owned()
}

///|
pub fn bigbed_parse_auto_sql(text : String) -> BigBedSchema raise BigBedError {
  let lines = text.split("\n").to_array()
  let mut table_name = ""
  let mut table_comment = ""
  let fields : Array[BigBedField] = []
  let mut saw_table = false
  let mut in_fields = false
  for raw in lines {
    let line = trim(raw.to_owned())
    if line.length() == 0 {
      continue
    }
    if !saw_table {
      let words = split_by_char(line, ' '.to_int())
      if words.length() < 2 || words[0] != "table" {
        bigbed_fail("AutoSQL declaration must start with 'table'")
      }
      table_name = words[1]
      saw_table = true
      continue
    }
    if table_comment.length() == 0 && !in_fields {
      if line == "(" {
        in_fields = true
      } else {
        table_comment = bigbed_quoted_value(line)
      }
      continue
    }
    if line == "(" {
      in_fields = true
      continue
    }
    if line == ")" {
      in_fields = false
      continue
    }
    if in_fields {
      let semicolon = bigbed_find_byte(line, ';'.to_int(), 0)
      if semicolon < 0 {
        bigbed_fail("AutoSQL field declaration is missing ';'")
      }
      let definition = trim(line[0:semicolon].to_owned())
      let words = split_by_char(definition, ' '.to_int())
      let compact : Array[String] = []
      for word in words {
        if word.length() > 0 {
          compact.push(word)
        }
      }
      if compact.length() != 2 {
        bigbed_fail("AutoSQL field requires a type and name")
      }
      fields.push(
        BigBedField::create(
          compact[0],
          compact[1],
          comment=bigbed_quoted_value(line),
        ),
      )
    }
  }
  if !saw_table || in_fields || fields.length() == 0 {
    bigbed_fail("Incomplete AutoSQL declaration")
  }
  BigBedSchema::create(table_name, table_comment, fields)
}

///|
fn bigbed_check_range(
  data : Array[Int],
  position : Int,
  size : Int,
  label : String,
) -> Unit raise BigBedError {
  if position < 0 || size < 0 || position > data.length() - size {
    bigbed_fail("Truncated BigBed " + label)
  }
}

///|
fn bigbed_read_u16(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
) -> Int raise BigBedError {
  bigbed_check_range(data, position, 2, "16-bit integer")
  if endian == BigBedLittleEndian {
    data[position] | (data[position + 1] << 8)
  } else {
    (data[position] << 8) | data[position + 1]
  }
}

///|
fn bigbed_read_u32(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
) -> Int raise BigBedError {
  bigbed_check_range(data, position, 4, "32-bit integer")
  if endian == BigBedLittleEndian {
    data[position] |
    (data[position + 1] << 8) |
    (data[position + 2] << 16) |
    (data[position + 3] << 24)
  } else {
    (data[position] << 24) |
    (data[position + 1] << 16) |
    (data[position + 2] << 8) |
    data[position + 3]
  }
}

///|
fn bigbed_read_u64_local(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
) -> Int raise BigBedError {
  let low_position = if endian == BigBedLittleEndian {
    position
  } else {
    position + 4
  }
  let high_position = if endian == BigBedLittleEndian {
    position + 4
  } else {
    position
  }
  let low = bigbed_read_u32(data, low_position, endian)
  let high = bigbed_read_u32(data, high_position, endian)
  if high != 0 || low < 0 {
    bigbed_fail("BigBed offset or count exceeds the local Int range")
  }
  low
}

///|
fn bigbed_push_u16(output : Array[Int], value : Int) -> Unit {
  output.push(value & 0xFF)
  output.push((value >> 8) & 0xFF)
}

///|
fn bigbed_push_u32(output : Array[Int], value : Int) -> Unit {
  output.push(value & 0xFF)
  output.push((value >> 8) & 0xFF)
  output.push((value >> 16) & 0xFF)
  output.push((value >> 24) & 0xFF)
}

///|
fn bigbed_push_u64_local(output : Array[Int], value : Int) -> Unit {
  bigbed_push_u32(output, value)
  bigbed_push_u32(output, 0)
}

///|
fn bigbed_patch_u16(output : Array[Int], position : Int, value : Int) -> Unit {
  output[position] = value & 0xFF
  output[position + 1] = (value >> 8) & 0xFF
}

///|
fn bigbed_patch_u32(output : Array[Int], position : Int, value : Int) -> Unit {
  output[position] = value & 0xFF
  output[position + 1] = (value >> 8) & 0xFF
  output[position + 2] = (value >> 16) & 0xFF
  output[position + 3] = (value >> 24) & 0xFF
}

///|
fn bigbed_patch_u64_local(
  output : Array[Int],
  position : Int,
  value : Int,
) -> Unit {
  bigbed_patch_u32(output, position, value)
  bigbed_patch_u32(output, position + 4, 0)
}

///|
fn bigbed_push_text(output : Array[Int], value : String) -> Unit {
  for i in 0..<value.length() {
    output.push(value.unsafe_get(i).to_int() & 0xFF)
  }
}

///|
fn bigbed_read_text(
  data : Array[Int],
  start : Int,
  end : Int,
) -> String raise BigBedError {
  bigbed_check_range(data, start, end - start, "text")
  let builder = StringBuilder::new()
  for i in start..<end {
    builder.write_char(Int::unsafe_to_char(data[i] & 0xFF))
  }
  builder.to_string()
}

///|
fn bigbed_read_c_string(
  data : Array[Int],
  start : Int,
  limit : Int,
) -> (String, Int) raise BigBedError {
  if start < 0 || limit > data.length() || start >= limit {
    bigbed_fail("Invalid BigBed string bounds")
  }
  let mut end = start
  while end < limit && data[end] != 0 {
    end = end + 1
  }
  if end >= limit {
    bigbed_fail("Unterminated BigBed string")
  }
  (bigbed_read_text(data, start, end), end + 1)
}

///|
fn bigbed_parse_nonnegative_int(
  value : String,
  label : String,
) -> Int raise BigBedError {
  if value.length() == 0 {
    bigbed_fail("Missing integer for " + label)
  }
  let mut result = 0
  for i in 0..<value.length() {
    let code = value.unsafe_get(i).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      bigbed_fail("Invalid integer '" + value + "' for " + label)
    }
    let digit = code - '0'.to_int()
    if result > 214748364 || (result == 214748364 && digit > 7) {
      bigbed_fail("Integer for " + label + " exceeds the local Int range")
    }
    result = result * 10 + digit
  }
  result
}

///|
fn bigbed_parse_csv_ints(
  value : String,
  label : String,
) -> Array[Int] raise BigBedError {
  let values : Array[Int] = []
  for word in split_by_char(value, ','.to_int()) {
    if word.length() > 0 {
      values.push(bigbed_parse_nonnegative_int(word, label))
    }
  }
  values
}

///|
fn bigbed_csv(values : Array[Int]) -> String {
  let output = StringBuilder::new()
  for value in values {
    output.write_string(value.to_string())
    output.write_string(",")
  }
  output.to_string()
}

///|
fn bigbed_schema_for_config(
  config : BigBedWriteConfig,
) -> BigBedSchema raise BigBedError {
  let fields = bigbed_standard_fields(config.bed_columns)
  for field in config.custom_fields {
    fields.push(field)
  }
  BigBedSchema::create(config.schema_name, config.schema_comment, fields)
}

///|
fn bigbed_validate_schema(
  schema : BigBedSchema,
  defined_fields : Int,
) -> Unit raise BigBedError {
  if schema.fields.length() < defined_fields {
    bigbed_fail("AutoSQL schema has fewer fields than the BigBed header")
  }
  for i in 0..<defined_fields {
    if schema.fields[i].name != bigbed_standard_names[i] {
      bigbed_fail(
        "Expected AutoSQL field '" +
        bigbed_standard_names[i] +
        "' at column " +
        (i + 1).to_string(),
      )
    }
  }
}

///|
fn bigbed_target_index(targets : Array[BigBedTarget], name : String) -> Int {
  for i in 0..<targets.length() {
    if targets[i].name == name {
      return i
    }
  }
  -1
}

///|
fn bigbed_validate_records(
  targets : Array[BigBedTarget],
  records : Array[BigBedRecord],
  config : BigBedWriteConfig,
) -> Unit raise BigBedError {
  if targets.length() == 0 {
    bigbed_fail("BigBed requires at least one target sequence")
  }
  let seen : Map[String, Bool] = Map([])
  for target in targets {
    ignore(BigBedTarget::create(target.name, target.length))
    if seen.contains(target.name) {
      bigbed_fail("Duplicate BigBed target '" + target.name + "'")
    }
    seen[target.name] = true
  }
  let mut previous_chrom = -1
  let mut previous_start = -1
  for i in 0..<records.length() {
    let record = records[i]
    let chrom_id = bigbed_target_index(targets, record.chrom)
    if chrom_id < 0 {
      bigbed_fail(
        "Record " +
        i.to_string() +
        " refers to unknown target '" +
        record.chrom +
        "'",
      )
    }
    if record.chrom_end > targets[chrom_id].length {
      bigbed_fail(
        "Record " + i.to_string() + " extends beyond its target sequence",
      )
    }
    if chrom_id < previous_chrom ||
      (chrom_id == previous_chrom && record.chrom_start < previous_start) {
      bigbed_fail("BigBed records must be sorted by target and start")
    }
    if config.bed_columns >= 10 && record.declared_block_count <= 0 {
      bigbed_fail("BED10-BED12 records require at least one block")
    }
    if config.bed_columns >= 11 &&
      record.block_sizes.length() != record.declared_block_count {
      bigbed_fail("BED11-BED12 block count must match blockSizes")
    }
    if config.bed_columns >= 12 &&
      record.block_starts.length() != record.declared_block_count {
      bigbed_fail("BED12 block count must match chromStarts")
    }
    if record.extra_fields.length() != config.custom_fields.length() {
      bigbed_fail(
        "Record " + i.to_string() + " has an incorrect number of custom fields",
      )
    }
    previous_chrom = chrom_id
    previous_start = record.chrom_start
  }
}

///|
fn bigbed_rest_fields(
  record : BigBedRecord,
  bed_columns : Int,
) -> Array[String] {
  let fields : Array[String] = []
  if bed_columns >= 4 {
    fields.push(record.name)
  }
  if bed_columns >= 5 {
    fields.push(record.score.to_string())
  }
  if bed_columns >= 6 {
    fields.push(record.strand)
  }
  if bed_columns >= 7 {
    fields.push(record.thick_start.to_string())
  }
  if bed_columns >= 8 {
    fields.push(record.thick_end.to_string())
  }
  if bed_columns >= 9 {
    fields.push(record.item_rgb)
  }
  if bed_columns >= 10 {
    fields.push(record.declared_block_count.to_string())
  }
  if bed_columns >= 11 {
    fields.push(bigbed_csv(record.block_sizes))
  }
  if bed_columns >= 12 {
    fields.push(bigbed_csv(record.block_starts))
  }
  for value in record.extra_fields {
    fields.push(value)
  }
  fields
}

///|
fn bigbed_join_tabs(fields : Array[String]) -> String {
  let output = StringBuilder::new()
  for i in 0..<fields.length() {
    if i > 0 {
      output.write_string("\t")
    }
    output.write_string(fields[i])
  }
  output.to_string()
}

///|
fn bigbed_record_bytes(
  record : BigBedRecord,
  chrom_id : Int,
  bed_columns : Int,
) -> Array[Int] {
  let output : Array[Int] = []
  bigbed_push_u32(output, chrom_id)
  bigbed_push_u32(output, record.chrom_start)
  bigbed_push_u32(output, record.chrom_end)
  bigbed_push_text(
    output,
    bigbed_join_tabs(bigbed_rest_fields(record, bed_columns)),
  )
  output.push(0)
  output
}

///|
fn bigbed_adler32(data : Array[Int]) -> Int {
  let mut a = 1
  let mut b = 0
  for byte in data {
    a = (a + (byte & 0xFF)) % 65521
    b = (b + a) % 65521
  }
  (b << 16) | a
}

///|
struct BigBedBitReader {
  data : Array[Int]
  bit_position : Int
  bit_limit : Int
}

///|
struct BigBedHuffman {
  counts : Array[Int]
  symbols : Array[Int]
}

///|
let bigbed_length_base : Array[Int] = [
  3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83,
  99, 115, 131, 163, 195, 227, 258,
]

///|
let bigbed_length_extra : Array[Int] = [
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5,
  5, 0,
]

///|
let bigbed_distance_base : Array[Int] = [
  1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769,
  1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
]

///|
let bigbed_distance_extra : Array[Int] = [
  0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
  12, 12, 13, 13,
]

///|
fn bigbed_read_deflate_bits(
  reader : BigBedBitReader,
  count : Int,
) -> (Int, BigBedBitReader) raise BigBedError {
  if count < 0 || count > 24 {
    bigbed_fail("Invalid DEFLATE bit count")
  }
  if reader.bit_position + count > reader.bit_limit {
    bigbed_fail("Truncated DEFLATE stream in BigBed data block")
  }
  let mut value = 0
  for i in 0..<count {
    let position = reader.bit_position + i
    let bit = (reader.data[position / 8] >> (position % 8)) & 1
    value = value | (bit << i)
  }
  (
    value,
    BigBedBitReader::{
      data: reader.data,
      bit_position: reader.bit_position + count,
      bit_limit: reader.bit_limit,
    },
  )
}

///|
fn bigbed_align_deflate_reader(
  reader : BigBedBitReader,
) -> BigBedBitReader raise BigBedError {
  let aligned = (reader.bit_position + 7) / 8 * 8
  if aligned > reader.bit_limit {
    bigbed_fail("Truncated DEFLATE byte alignment")
  }
  BigBedBitReader::{
    data: reader.data,
    bit_position: aligned,
    bit_limit: reader.bit_limit,
  }
}

///|
fn bigbed_build_huffman(
  lengths : Array[Int],
  label : String,
) -> BigBedHuffman raise BigBedError {
  let counts : Array[Int] = Array::make(16, 0)
  for length in lengths {
    if length < 0 || length > 15 {
      bigbed_fail("Invalid " + label + " Huffman code length")
    }
    counts[length] = counts[length] + 1
  }
  if counts[0] == lengths.length() {
    bigbed_fail("Empty " + label + " Huffman tree")
  }
  let mut remaining = 1
  for bits in 1..<16 {
    remaining = (remaining << 1) - counts[bits]
    if remaining < 0 {
      bigbed_fail("Oversubscribed " + label + " Huffman tree")
    }
  }
  let offsets : Array[Int] = Array::make(16, 0)
  for bits in 1..<15 {
    offsets[bits + 1] = offsets[bits] + counts[bits]
  }
  let next = offsets.copy()
  let symbols : Array[Int] = Array::make(lengths.length() - counts[0], 0)
  for symbol in 0..<lengths.length() {
    let length = lengths[symbol]
    if length > 0 {
      symbols[next[length]] = symbol
      next[length] = next[length] + 1
    }
  }
  BigBedHuffman::{ counts, symbols }
}

///|
fn bigbed_decode_huffman_symbol(
  tree : BigBedHuffman,
  input : BigBedBitReader,
) -> (Int, BigBedBitReader) raise BigBedError {
  let mut reader = input
  let mut code = 0
  let mut first = 0
  let mut index = 0
  for length in 1..<16 {
    let (bit, next) = bigbed_read_deflate_bits(reader, 1)
    reader = next
    code = code | bit
    let count = tree.counts[length]
    if code < first + count {
      let symbol_index = index + code - first
      if symbol_index < 0 || symbol_index >= tree.symbols.length() {
        bigbed_fail("Invalid canonical Huffman symbol index")
      }
      return (tree.symbols[symbol_index], reader)
    }
    index = index + count
    first = (first + count) << 1
    code = code << 1
  }
  bigbed_fail("Invalid Huffman code in BigBed DEFLATE stream")
  (0, reader)
}

///|
fn bigbed_fixed_huffman() -> (BigBedHuffman, BigBedHuffman) raise BigBedError {
  let literal_lengths : Array[Int] = Array::make(288, 0)
  for i in 0..<144 {
    literal_lengths[i] = 8
  }
  for i in 144..<256 {
    literal_lengths[i] = 9
  }
  for i in 256..<280 {
    literal_lengths[i] = 7
  }
  for i in 280..<288 {
    literal_lengths[i] = 8
  }
  let distance_lengths : Array[Int] = Array::make(32, 5)
  (
    bigbed_build_huffman(literal_lengths, "fixed literal/length"),
    bigbed_build_huffman(distance_lengths, "fixed distance"),
  )
}

///|
fn bigbed_dynamic_huffman(
  input : BigBedBitReader,
) -> (BigBedHuffman, BigBedHuffman, BigBedBitReader) raise BigBedError {
  let mut reader = input
  let (literal_bits, next) = bigbed_read_deflate_bits(reader, 5)
  reader = next
  let literal_count = literal_bits + 257
  let (distance_bits, next) = bigbed_read_deflate_bits(reader, 5)
  reader = next
  let distance_count = distance_bits + 1
  let (code_bits, next) = bigbed_read_deflate_bits(reader, 4)
  reader = next
  let code_count = code_bits + 4
  if literal_count > 286 || distance_count > 32 {
    bigbed_fail("Invalid dynamic DEFLATE tree dimensions")
  }
  let code_order = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
  ]
  let code_lengths : Array[Int] = Array::make(19, 0)
  for i in 0..<code_count {
    let (length, next) = bigbed_read_deflate_bits(reader, 3)
    reader = next
    code_lengths[code_order[i]] = length
  }
  let code_tree = bigbed_build_huffman(code_lengths, "code-length")
  let total = literal_count + distance_count
  let lengths : Array[Int] = Array::make(total, 0)
  let mut position = 0
  while position < total {
    let (symbol, next) = bigbed_decode_huffman_symbol(code_tree, reader)
    reader = next
    if symbol < 16 {
      lengths[position] = symbol
      position = position + 1
    } else if symbol == 16 {
      if position == 0 {
        bigbed_fail("DEFLATE repeat code has no previous code length")
      }
      let (extra, next) = bigbed_read_deflate_bits(reader, 2)
      reader = next
      let repeat = extra + 3
      if position + repeat > total {
        bigbed_fail("DEFLATE code-length repeat exceeds its tree")
      }
      let previous = lengths[position - 1]
      for _ in 0..<repeat {
        lengths[position] = previous
        position = position + 1
      }
    } else if symbol == 17 {
      let (extra, next) = bigbed_read_deflate_bits(reader, 3)
      reader = next
      let repeat = extra + 3
      if position + repeat > total {
        bigbed_fail("DEFLATE zero repeat exceeds its tree")
      }
      position = position + repeat
    } else if symbol == 18 {
      let (extra, next) = bigbed_read_deflate_bits(reader, 7)
      reader = next
      let repeat = extra + 11
      if position + repeat > total {
        bigbed_fail("DEFLATE long zero repeat exceeds its tree")
      }
      position = position + repeat
    } else {
      bigbed_fail("Invalid DEFLATE code-length symbol")
    }
  }
  let literal_lengths : Array[Int] = []
  for i in 0..<literal_count {
    literal_lengths.push(lengths[i])
  }
  if literal_lengths[256] == 0 {
    bigbed_fail("Dynamic DEFLATE tree has no end-of-block symbol")
  }
  let distance_lengths : Array[Int] = []
  for i in 0..<distance_count {
    distance_lengths.push(lengths[literal_count + i])
  }
  (
    bigbed_build_huffman(literal_lengths, "dynamic literal/length"),
    bigbed_build_huffman(distance_lengths, "dynamic distance"),
    reader,
  )
}

///|
fn bigbed_inflate_huffman_block(
  input : BigBedBitReader,
  literal_tree : BigBedHuffman,
  distance_tree : BigBedHuffman,
  output : Array[Int],
  maximum_output : Int,
) -> BigBedBitReader raise BigBedError {
  let mut reader = input
  let mut done = false
  while !done {
    let (symbol, next) = bigbed_decode_huffman_symbol(literal_tree, reader)
    reader = next
    if symbol < 256 {
      if output.length() >= maximum_output {
        bigbed_fail("BigBed data block exceeds the declared buffer size")
      }
      output.push(symbol)
    } else if symbol == 256 {
      done = true
    } else {
      if symbol < 257 || symbol > 285 {
        bigbed_fail("Invalid DEFLATE length symbol")
      }
      let length_index = symbol - 257
      let (length_extra, next) = bigbed_read_deflate_bits(
        reader,
        bigbed_length_extra[length_index],
      )
      reader = next
      let length = bigbed_length_base[length_index] + length_extra
      let (distance_symbol, next) = bigbed_decode_huffman_symbol(
        distance_tree, reader,
      )
      reader = next
      if distance_symbol < 0 || distance_symbol >= 30 {
        bigbed_fail("Invalid DEFLATE distance symbol")
      }
      let (distance_extra, next) = bigbed_read_deflate_bits(
        reader,
        bigbed_distance_extra[distance_symbol],
      )
      reader = next
      let distance = bigbed_distance_base[distance_symbol] + distance_extra
      if distance <= 0 || distance > output.length() || distance > 32768 {
        bigbed_fail("Invalid LZ77 distance in BigBed DEFLATE stream")
      }
      if output.length() + length > maximum_output {
        bigbed_fail("BigBed data block exceeds the declared buffer size")
      }
      for _ in 0..<length {
        output.push(output[output.length() - distance])
      }
    }
  }
  reader
}

///|
fn bigbed_zlib_store(data : Array[Int]) -> Array[Int] {
  let output : Array[Int] = [0x78, 0x01]
  let mut position = 0
  if data.length() == 0 {
    output.push(1)
    bigbed_push_u16(output, 0)
    bigbed_push_u16(output, 0xFFFF)
  }
  while position < data.length() {
    let remaining = data.length() - position
    let size = if remaining > 65535 { 65535 } else { remaining }
    let final_block = position + size == data.length()
    output.push(if final_block { 1 } else { 0 })
    bigbed_push_u16(output, size)
    bigbed_push_u16(output, size ^ 0xFFFF)
    for i in position..<(position + size) {
      output.push(data[i])
    }
    position = position + size
  }
  let checksum = bigbed_adler32(data)
  output.push((checksum >> 24) & 0xFF)
  output.push((checksum >> 16) & 0xFF)
  output.push((checksum >> 8) & 0xFF)
  output.push(checksum & 0xFF)
  output
}

///|
fn bigbed_zlib_unstore(
  data : Array[Int],
  maximum_output : Int,
) -> Array[Int] raise BigBedError {
  if data.length() < 8 {
    bigbed_fail("Truncated zlib data block")
  }
  let cmf = data[0]
  let flags = data[1]
  if (cmf & 0x0F) != 8 || cmf >> 4 > 7 || ((cmf << 8) | flags) % 31 != 0 {
    bigbed_fail("Invalid zlib header in BigBed data block")
  }
  if (flags & 0x20) != 0 {
    bigbed_fail("Preset zlib dictionaries are not supported")
  }
  let payload_end = data.length() - 4
  let output : Array[Int] = []
  let mut reader = BigBedBitReader::{
    data,
    bit_position: 16,
    bit_limit: payload_end * 8,
  }
  let mut done = false
  while !done {
    let (final_block, next) = bigbed_read_deflate_bits(reader, 1)
    reader = next
    let (block_type, next) = bigbed_read_deflate_bits(reader, 2)
    reader = next
    done = final_block == 1
    if block_type == 0 {
      reader = bigbed_align_deflate_reader(reader)
      let (size, next) = bigbed_read_deflate_bits(reader, 16)
      reader = next
      let (inverse, next) = bigbed_read_deflate_bits(reader, 16)
      reader = next
      if (size ^ 0xFFFF) != inverse {
        bigbed_fail("Invalid stored DEFLATE block length")
      }
      if output.length() + size > maximum_output {
        bigbed_fail("BigBed data block exceeds the declared buffer size")
      }
      for _ in 0..<size {
        let (byte, next) = bigbed_read_deflate_bits(reader, 8)
        reader = next
        output.push(byte)
      }
    } else if block_type == 1 {
      let (literal_tree, distance_tree) = bigbed_fixed_huffman()
      reader = bigbed_inflate_huffman_block(
        reader, literal_tree, distance_tree, output, maximum_output,
      )
    } else if block_type == 2 {
      let (literal_tree, distance_tree, next) = bigbed_dynamic_huffman(reader)
      reader = bigbed_inflate_huffman_block(
        next, literal_tree, distance_tree, output, maximum_output,
      )
    } else {
      bigbed_fail("Reserved DEFLATE block type in BigBed data block")
    }
  }
  reader = bigbed_align_deflate_reader(reader)
  if reader.bit_position != payload_end * 8 {
    bigbed_fail("Unexpected bytes after final DEFLATE block")
  }
  let stored_checksum = (data[payload_end] << 24) |
    (data[payload_end + 1] << 16) |
    (data[payload_end + 2] << 8) |
    data[payload_end + 3]
  if bigbed_adler32(output) != stored_checksum {
    bigbed_fail("BigBed zlib Adler-32 checksum mismatch")
  }
  output
}

///|
fn bigbed_append(output : Array[Int], source : Array[Int]) -> Unit {
  for byte in source {
    output.push(byte)
  }
}

///|
fn bigbed_patch_header(
  output : Array[Int],
  chromosome_tree_offset : Int,
  full_data_offset : Int,
  full_index_offset : Int,
  field_count : Int,
  defined_field_count : Int,
  auto_sql_offset : Int,
  total_summary_offset : Int,
  uncompress_buffer_size : Int,
) -> Unit {
  bigbed_patch_u32(output, 0, bigbed_magic)
  bigbed_patch_u16(output, 4, bigbed_version)
  bigbed_patch_u16(output, 6, 0)
  bigbed_patch_u64_local(output, 8, chromosome_tree_offset)
  bigbed_patch_u64_local(output, 16, full_data_offset)
  bigbed_patch_u64_local(output, 24, full_index_offset)
  bigbed_patch_u16(output, 32, field_count)
  bigbed_patch_u16(output, 34, defined_field_count)
  bigbed_patch_u64_local(output, 36, auto_sql_offset)
  bigbed_patch_u64_local(output, 44, total_summary_offset)
  bigbed_patch_u32(output, 52, uncompress_buffer_size)
  bigbed_patch_u64_local(output, 56, 0)
}

///|
fn bigbed_write_tree_key(
  output : Array[Int],
  key : String,
  key_size : Int,
) -> Unit {
  for i in 0..<key_size {
    output.push(
      if i < key.length() {
        key.unsafe_get(i).to_int() & 0xFF
      } else {
        0
      },
    )
  }
}

///|
fn bigbed_tree_height(item_count : Int, block_size : Int) -> Int {
  let mut height = 1
  let mut capacity = block_size
  while capacity < item_count {
    capacity = if capacity > item_count / block_size {
      item_count
    } else {
      capacity * block_size
    }
    height = height + 1
  }
  height
}

///|
fn bigbed_tree_capacity(
  item_count : Int,
  block_size : Int,
  height : Int,
) -> Int {
  let mut capacity = 1
  for _ in 0..<height {
    capacity = if capacity > item_count / block_size {
      item_count
    } else {
      capacity * block_size
    }
  }
  capacity
}

///|
fn bigbed_write_chrom_node(
  output : Array[Int],
  targets : Array[BigBedTarget],
  order : Array[Int],
  start : Int,
  end : Int,
  key_size : Int,
  block_size : Int,
  height : Int,
) -> Unit {
  let count = end - start
  if height == 1 {
    output.push(1)
    output.push(0)
    bigbed_push_u16(output, count)
    for position in start..<end {
      let chrom_id = order[position]
      let target = targets[chrom_id]
      bigbed_write_tree_key(output, target.name, key_size)
      bigbed_push_u32(output, chrom_id)
      bigbed_push_u32(output, target.length)
    }
    return
  }
  let child_capacity = bigbed_tree_capacity(count, block_size, height - 1)
  let child_count = (count + child_capacity - 1) / child_capacity
  output.push(0)
  output.push(0)
  bigbed_push_u16(output, child_count)
  let patch_positions : Array[Int] = []
  let child_starts : Array[Int] = []
  let child_ends : Array[Int] = []
  let mut cursor = start
  for child in 0..<child_count {
    let remaining_items = end - cursor
    let remaining_children = child_count - child
    let child_size = (remaining_items + remaining_children - 1) /
      remaining_children
    let child_end = cursor + child_size
    bigbed_write_tree_key(output, targets[order[cursor]].name, key_size)
    patch_positions.push(output.length())
    bigbed_push_u64_local(output, 0)
    child_starts.push(cursor)
    child_ends.push(child_end)
    cursor = child_end
  }
  for child in 0..<child_count {
    bigbed_patch_u64_local(output, patch_positions[child], output.length())
    bigbed_write_chrom_node(
      output,
      targets,
      order,
      child_starts[child],
      child_ends[child],
      key_size,
      block_size,
      height - 1,
    )
  }
}

///|
fn bigbed_write_chrom_tree(
  output : Array[Int],
  targets : Array[BigBedTarget],
  configured_block_size : Int,
) -> Unit {
  let mut key_size = 1
  for target in targets {
    if target.name.length() > key_size {
      key_size = target.name.length()
    }
  }
  let order : Array[Int] = []
  for i in 0..<targets.length() {
    order.push(i)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if targets[left].name < targets[right].name {
      -1
    } else if targets[left].name > targets[right].name {
      1
    } else {
      0
    }
  })
  bigbed_push_u32(output, bigbed_chrom_tree_magic)
  bigbed_push_u32(output, configured_block_size)
  bigbed_push_u32(output, key_size)
  bigbed_push_u32(output, 8)
  bigbed_push_u64_local(output, targets.length())
  bigbed_push_u64_local(output, 0)
  bigbed_write_chrom_node(
    output,
    targets,
    order,
    0,
    targets.length(),
    key_size,
    configured_block_size,
    bigbed_tree_height(targets.length(), configured_block_size),
  )
}

///|
fn bigbed_block_overlaps(
  block : BigBedBlockIndex,
  chrom_id : Int,
  start : Int,
  end : Int,
) -> Bool {
  let starts_before_query_end = block.start_chrom < chrom_id ||
    (block.start_chrom == chrom_id && block.start_base < end)
  let query_starts_before_block_end = chrom_id < block.end_chrom ||
    (chrom_id == block.end_chrom && start < block.end_base)
  starts_before_query_end && query_starts_before_block_end
}

///|
fn bigbed_rtree_end(
  blocks : Array[BigBedBlockIndex],
  start : Int,
  end : Int,
) -> (Int, Int) {
  let mut end_chrom = blocks[start].end_chrom
  let mut end_base = blocks[start].end_base
  for i in (start + 1)..<end {
    let block = blocks[i]
    if block.end_chrom > end_chrom ||
      (block.end_chrom == end_chrom && block.end_base > end_base) {
      end_chrom = block.end_chrom
      end_base = block.end_base
    }
  }
  (end_chrom, end_base)
}

///|
fn bigbed_write_rtree_node(
  output : Array[Int],
  blocks : Array[BigBedBlockIndex],
  start : Int,
  end : Int,
  block_size : Int,
  height : Int,
) -> Unit {
  let count = end - start
  if height == 1 {
    output.push(1)
    output.push(0)
    bigbed_push_u16(output, count)
    for i in start..<end {
      let block = blocks[i]
      bigbed_push_u32(output, block.start_chrom)
      bigbed_push_u32(output, block.start_base)
      bigbed_push_u32(output, block.end_chrom)
      bigbed_push_u32(output, block.end_base)
      bigbed_push_u64_local(output, block.data_offset)
      bigbed_push_u64_local(output, block.data_size)
    }
    return
  }
  let child_capacity = bigbed_tree_capacity(count, block_size, height - 1)
  let child_count = (count + child_capacity - 1) / child_capacity
  output.push(0)
  output.push(0)
  bigbed_push_u16(output, child_count)
  let patch_positions : Array[Int] = []
  let child_starts : Array[Int] = []
  let child_ends : Array[Int] = []
  let mut cursor = start
  for child in 0..<child_count {
    let remaining_items = end - cursor
    let remaining_children = child_count - child
    let child_size = (remaining_items + remaining_children - 1) /
      remaining_children
    let child_end = cursor + child_size
    let (end_chrom, end_base) = bigbed_rtree_end(blocks, cursor, child_end)
    bigbed_push_u32(output, blocks[cursor].start_chrom)
    bigbed_push_u32(output, blocks[cursor].start_base)
    bigbed_push_u32(output, end_chrom)
    bigbed_push_u32(output, end_base)
    patch_positions.push(output.length())
    bigbed_push_u64_local(output, 0)
    child_starts.push(cursor)
    child_ends.push(child_end)
    cursor = child_end
  }
  for child in 0..<child_count {
    bigbed_patch_u64_local(output, patch_positions[child], output.length())
    bigbed_write_rtree_node(
      output,
      blocks,
      child_starts[child],
      child_ends[child],
      block_size,
      height - 1,
    )
  }
}

///|
fn bigbed_write_rtree(
  output : Array[Int],
  blocks : Array[BigBedBlockIndex],
  configured_block_size : Int,
  item_count : Int,
  items_per_slot : Int,
  end_file_offset : Int,
) -> Unit {
  bigbed_push_u32(output, bigbed_rtree_magic)
  bigbed_push_u32(output, configured_block_size)
  bigbed_push_u64_local(output, item_count)
  if blocks.length() == 0 {
    for _ in 0..<4 {
      bigbed_push_u32(output, 0)
    }
  } else {
    bigbed_push_u32(output, blocks[0].start_chrom)
    bigbed_push_u32(output, blocks[0].start_base)
    let (end_chrom, end_base) = bigbed_rtree_end(blocks, 0, blocks.length())
    bigbed_push_u32(output, end_chrom)
    bigbed_push_u32(output, end_base)
  }
  bigbed_push_u64_local(output, end_file_offset)
  bigbed_push_u32(output, items_per_slot)
  bigbed_push_u32(output, 0)
  if blocks.length() == 0 {
    output.push(1)
    output.push(0)
    bigbed_push_u16(output, 0)
  } else {
    bigbed_write_rtree_node(
      output,
      blocks,
      0,
      blocks.length(),
      configured_block_size,
      bigbed_tree_height(blocks.length(), configured_block_size),
    )
  }
}

///|
pub fn bigbed_write(
  targets : Array[BigBedTarget],
  records : Array[BigBedRecord],
  config? : BigBedWriteConfig = BigBedWriteConfig::default(),
) -> Array[Int] raise BigBedError {
  bigbed_validate_records(targets, records, config)
  let schema = bigbed_schema_for_config(config)
  let output : Array[Int] = Array::make(64, 0)
  let auto_sql_offset = output.length()
  bigbed_push_text(output, schema.to_auto_sql())
  output.push(0)
  let total_summary_offset = output.length()
  for _ in 0..<40 {
    output.push(0)
  }
  let chromosome_tree_offset = output.length()
  bigbed_write_chrom_tree(output, targets, config.block_size)
  let full_data_offset = output.length()
  bigbed_push_u64_local(output, records.length())
  let blocks : Array[BigBedBlockIndex] = []
  let mut maximum_uncompressed = 0
  let mut record_index = 0
  while record_index < records.length() {
    let block_start_index = record_index
    let block_end_index = if record_index + config.items_per_slot <
      records.length() {
      record_index + config.items_per_slot
    } else {
      records.length()
    }
    let raw : Array[Int] = []
    while record_index < block_end_index {
      let record = records[record_index]
      let chrom_id = bigbed_target_index(targets, record.chrom)
      bigbed_append(
        raw,
        bigbed_record_bytes(record, chrom_id, config.bed_columns),
      )
      record_index = record_index + 1
    }
    if raw.length() > maximum_uncompressed {
      maximum_uncompressed = raw.length()
    }
    let payload = if config.compress { bigbed_zlib_store(raw) } else { raw }
    let data_offset = output.length()
    bigbed_append(output, payload)
    let first = records[block_start_index]
    let mut end_chrom = bigbed_target_index(targets, first.chrom)
    let mut end_base = first.chrom_end
    for i in (block_start_index + 1)..<block_end_index {
      let chrom_id = bigbed_target_index(targets, records[i].chrom)
      if chrom_id > end_chrom ||
        (chrom_id == end_chrom && records[i].chrom_end > end_base) {
        end_chrom = chrom_id
        end_base = records[i].chrom_end
      }
    }
    blocks.push(BigBedBlockIndex::{
      start_chrom: bigbed_target_index(targets, first.chrom),
      start_base: first.chrom_start,
      end_chrom,
      end_base,
      data_offset,
      data_size: payload.length(),
    })
  }
  let full_index_offset = output.length()
  bigbed_write_rtree(
    output,
    blocks,
    config.block_size,
    records.length(),
    config.items_per_slot,
    full_index_offset,
  )
  bigbed_push_u32(output, bigbed_magic)
  bigbed_patch_header(
    output,
    chromosome_tree_offset,
    full_data_offset,
    full_index_offset,
    schema.fields.length(),
    config.bed_columns,
    auto_sql_offset,
    total_summary_offset,
    if config.compress {
      maximum_uncompressed
    } else {
      0
    },
  )
  output
}

///|
fn bigbed_detect_endian(data : Array[Int]) -> BigBedEndian raise BigBedError {
  if data.length() < 4 {
    bigbed_fail("BigBed input is shorter than its magic number")
  }
  let little = data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24)
  if little == bigbed_magic {
    BigBedLittleEndian
  } else {
    let big = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3]
    if big == bigbed_magic {
      BigBedBigEndian
    } else {
      bigbed_fail("Input does not contain the BigBed magic number")
      BigBedLittleEndian
    }
  }
}

///|
fn bigbed_read_header(
  data : Array[Int],
) -> (BigBedHeader, BigBedEndian) raise BigBedError {
  bigbed_check_range(data, 0, 64, "header")
  let endian = bigbed_detect_endian(data)
  let version = bigbed_read_u16(data, 4, endian)
  if version != bigbed_version {
    bigbed_fail(
      "Unsupported BigBed version " +
      version.to_string() +
      "; expected version 4",
    )
  }
  let defined = bigbed_read_u16(data, 34, endian)
  if defined < 3 || defined > 12 {
    bigbed_fail("BigBed defined field count must be between 3 and 12")
  }
  let field_count = bigbed_read_u16(data, 32, endian)
  if field_count < defined {
    bigbed_fail("BigBed field count is smaller than its BED field count")
  }
  let header = BigBedHeader::{
    little_endian: endian == BigBedLittleEndian,
    version,
    zoom_levels: bigbed_read_u16(data, 6, endian),
    chromosome_tree_offset: bigbed_read_u64_local(data, 8, endian),
    full_data_offset: bigbed_read_u64_local(data, 16, endian),
    full_index_offset: bigbed_read_u64_local(data, 24, endian),
    field_count,
    defined_field_count: defined,
    auto_sql_offset: bigbed_read_u64_local(data, 36, endian),
    total_summary_offset: bigbed_read_u64_local(data, 44, endian),
    uncompress_buffer_size: bigbed_read_u32(data, 52, endian),
    extra_indices_offset: bigbed_read_u64_local(data, 56, endian),
  }
  if header.chromosome_tree_offset < 64 ||
    header.full_data_offset <= header.chromosome_tree_offset ||
    header.full_index_offset <= header.full_data_offset ||
    header.full_index_offset >= data.length() {
    bigbed_fail("BigBed header contains invalid section offsets")
  }
  (header, endian)
}

///|
fn bigbed_read_key(
  data : Array[Int],
  position : Int,
  key_size : Int,
) -> String raise BigBedError {
  bigbed_check_range(data, position, key_size, "B+ tree key")
  let mut end = position
  while end < position + key_size && data[end] != 0 {
    end = end + 1
  }
  bigbed_read_text(data, position, end)
}

///|
fn bigbed_read_chrom_node(
  data : Array[Int],
  position : Int,
  key_size : Int,
  endian : BigBedEndian,
  targets : Array[BigBedTarget?],
  depth : Int,
) -> Unit raise BigBedError {
  if depth > 64 {
    bigbed_fail("BigBed chromosome B+ tree exceeds 64 levels")
  }
  bigbed_check_range(data, position, 4, "B+ tree node")
  let is_leaf = data[position] != 0
  let count = bigbed_read_u16(data, position + 2, endian)
  let mut cursor = position + 4
  if is_leaf {
    for _ in 0..<count {
      let name = bigbed_read_key(data, cursor, key_size)
      cursor = cursor + key_size
      let chrom_id = bigbed_read_u32(data, cursor, endian)
      let chrom_size = bigbed_read_u32(data, cursor + 4, endian)
      cursor = cursor + 8
      if chrom_id < 0 || chrom_id >= targets.length() {
        bigbed_fail("Chromosome B+ tree contains an invalid target ID")
      }
      if targets[chrom_id] is Some(_) {
        bigbed_fail("Chromosome B+ tree contains a duplicate target ID")
      }
      targets[chrom_id] = Some(BigBedTarget::create(name, chrom_size))
    }
  } else {
    let children : Array[Int] = []
    for _ in 0..<count {
      ignore(bigbed_read_key(data, cursor, key_size))
      cursor = cursor + key_size
      children.push(bigbed_read_u64_local(data, cursor, endian))
      cursor = cursor + 8
    }
    for child in children {
      bigbed_read_chrom_node(data, child, key_size, endian, targets, depth + 1)
    }
  }
}

///|
fn bigbed_read_chrom_tree(
  data : Array[Int],
  offset : Int,
  endian : BigBedEndian,
) -> Array[BigBedTarget] raise BigBedError {
  bigbed_check_range(data, offset, 32, "chromosome B+ tree header")
  if bigbed_read_u32(data, offset, endian) != bigbed_chrom_tree_magic {
    bigbed_fail("Invalid BigBed chromosome B+ tree magic number")
  }
  let key_size = bigbed_read_u32(data, offset + 8, endian)
  let value_size = bigbed_read_u32(data, offset + 12, endian)
  let item_count = bigbed_read_u64_local(data, offset + 16, endian)
  if key_size <= 0 || value_size != 8 || item_count <= 0 {
    bigbed_fail("Invalid BigBed chromosome B+ tree dimensions")
  }
  let slots : Array[BigBedTarget?] = Array::make(item_count, None)
  bigbed_read_chrom_node(data, offset + 32, key_size, endian, slots, 0)
  let targets : Array[BigBedTarget] = []
  for slot in slots {
    match slot {
      Some(target) => targets.push(target)
      None => bigbed_fail("Chromosome B+ tree is missing a target ID")
    }
  }
  targets
}

///|
fn bigbed_read_rtree_node(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
  blocks : Array[BigBedBlockIndex],
  depth : Int,
) -> Unit raise BigBedError {
  if depth > 64 {
    bigbed_fail("BigBed R-tree exceeds 64 levels")
  }
  bigbed_check_range(data, position, 4, "R-tree node")
  let is_leaf = data[position] != 0
  let count = bigbed_read_u16(data, position + 2, endian)
  let mut cursor = position + 4
  if is_leaf {
    for _ in 0..<count {
      bigbed_check_range(data, cursor, 32, "R-tree leaf")
      blocks.push(BigBedBlockIndex::{
        start_chrom: bigbed_read_u32(data, cursor, endian),
        start_base: bigbed_read_u32(data, cursor + 4, endian),
        end_chrom: bigbed_read_u32(data, cursor + 8, endian),
        end_base: bigbed_read_u32(data, cursor + 12, endian),
        data_offset: bigbed_read_u64_local(data, cursor + 16, endian),
        data_size: bigbed_read_u64_local(data, cursor + 24, endian),
      })
      cursor = cursor + 32
    }
  } else {
    let children : Array[Int] = []
    for _ in 0..<count {
      bigbed_check_range(data, cursor, 24, "R-tree branch")
      children.push(bigbed_read_u64_local(data, cursor + 16, endian))
      cursor = cursor + 24
    }
    for child in children {
      bigbed_read_rtree_node(data, child, endian, blocks, depth + 1)
    }
  }
}

///|
fn bigbed_read_rtree(
  data : Array[Int],
  offset : Int,
  endian : BigBedEndian,
) -> (Array[BigBedBlockIndex], Int) raise BigBedError {
  bigbed_check_range(data, offset, 48, "R-tree header")
  if bigbed_read_u32(data, offset, endian) != bigbed_rtree_magic {
    bigbed_fail("Invalid BigBed R-tree magic number")
  }
  let indexed_count = bigbed_read_u64_local(data, offset + 8, endian)
  let blocks : Array[BigBedBlockIndex] = []
  bigbed_read_rtree_node(data, offset + 48, endian, blocks, 0)
  for block in blocks {
    bigbed_check_range(
      data,
      block.data_offset,
      block.data_size,
      "indexed data block",
    )
  }
  (blocks, indexed_count)
}

///|
fn bigbed_read_schema(
  data : Array[Int],
  header : BigBedHeader,
) -> BigBedSchema raise BigBedError {
  if header.auto_sql_offset == 0 {
    return BigBedSchema::default(bed_columns=header.defined_field_count)
  }
  if header.auto_sql_offset < 64 ||
    header.auto_sql_offset >= header.chromosome_tree_offset {
    bigbed_fail("BigBed AutoSQL offset is outside the metadata section")
  }
  let (text, _) = bigbed_read_c_string(
    data,
    header.auto_sql_offset,
    header.chromosome_tree_offset,
  )
  let schema = bigbed_parse_auto_sql(text)
  if schema.fields.length() != header.field_count {
    bigbed_fail("AutoSQL field count does not match the BigBed header")
  }
  bigbed_validate_schema(schema, header.defined_field_count)
  schema
}

///|
pub fn bigbed_parse(data : Array[Int]) -> BigBedFile raise BigBedError {
  let (header, endian) = bigbed_read_header(data)
  let trailing = bigbed_read_u32(data, data.length() - 4, endian)
  if trailing != bigbed_magic {
    bigbed_fail("BigBed trailing magic number is missing")
  }
  let schema = bigbed_read_schema(data, header)
  let targets = bigbed_read_chrom_tree(
    data,
    header.chromosome_tree_offset,
    endian,
  )
  let record_count = bigbed_read_u64_local(
    data,
    header.full_data_offset,
    endian,
  )
  let (blocks, indexed_count) = bigbed_read_rtree(
    data,
    header.full_index_offset,
    endian,
  )
  if indexed_count != record_count {
    bigbed_fail("BigBed data and R-tree record counts disagree")
  }
  BigBedFile::{
    header,
    schema,
    targets,
    source: data.copy(),
    blocks,
    record_count,
  }
}

///|
fn bigbed_decode_block(
  file : BigBedFile,
  block : BigBedBlockIndex,
) -> Array[Int] raise BigBedError {
  let payload : Array[Int] = []
  for i in block.data_offset..<(block.data_offset + block.data_size) {
    payload.push(file.source[i])
  }
  if file.header.uncompress_buffer_size > 0 {
    bigbed_zlib_unstore(payload, file.header.uncompress_buffer_size)
  } else {
    payload
  }
}

///|
fn bigbed_record_from_words(
  file : BigBedFile,
  chrom_id : Int,
  chrom_start : Int,
  chrom_end : Int,
  words : Array[String],
) -> BigBedRecord raise BigBedError {
  if chrom_id < 0 || chrom_id >= file.targets.length() {
    bigbed_fail("Binary BED record contains an invalid chromosome ID")
  }
  let expected_words = file.header.field_count - 3
  if words.length() != expected_words {
    bigbed_fail(
      "Binary BED record contains " +
      words.length().to_string() +
      " fields after its coordinates; expected " +
      expected_words.to_string(),
    )
  }
  let bed_columns = file.header.defined_field_count
  let name = if bed_columns >= 4 { words[0] } else { "." }
  let score = if bed_columns >= 5 {
    bigbed_parse_nonnegative_int(words[1], "BED score")
  } else {
    0
  }
  let strand = if bed_columns >= 6 { words[2] } else { "." }
  let thick_start = if bed_columns >= 7 {
    bigbed_parse_nonnegative_int(words[3], "BED thickStart")
  } else {
    chrom_start
  }
  let thick_end = if bed_columns >= 8 {
    bigbed_parse_nonnegative_int(words[4], "BED thickEnd")
  } else {
    chrom_end
  }
  let item_rgb = if bed_columns >= 9 { words[5] } else { "0" }
  let declared_blocks = if bed_columns >= 10 {
    bigbed_parse_nonnegative_int(words[6], "BED blockCount")
  } else {
    0
  }
  let block_sizes = if bed_columns >= 11 {
    bigbed_parse_csv_ints(words[7], "BED blockSizes")
  } else {
    []
  }
  let block_starts = if bed_columns >= 12 {
    bigbed_parse_csv_ints(words[8], "BED chromStarts")
  } else {
    []
  }
  if bed_columns >= 11 && block_sizes.length() != declared_blocks {
    bigbed_fail("BED blockCount does not match blockSizes")
  }
  if bed_columns >= 12 && block_starts.length() != declared_blocks {
    bigbed_fail("BED blockCount does not match chromStarts")
  }
  let extra : Array[String] = []
  for i in (bed_columns - 3)..<words.length() {
    extra.push(words[i])
  }
  BigBedRecord::create(
    file.targets[chrom_id].name,
    chrom_start,
    chrom_end,
    name~,
    score~,
    strand~,
    thick_start~,
    thick_end~,
    item_rgb~,
    block_count=declared_blocks,
    block_sizes~,
    block_starts~,
    extra_fields=extra,
  )
}

///|
fn bigbed_decode_records(
  file : BigBedFile,
  block : BigBedBlockIndex,
) -> Array[BigBedRecord] raise BigBedError {
  let data = bigbed_decode_block(file, block)
  let endian = if file.header.little_endian {
    BigBedLittleEndian
  } else {
    BigBedBigEndian
  }
  let records : Array[BigBedRecord] = []
  let mut position = 0
  while position < data.length() {
    bigbed_check_range(data, position, 13, "binary BED record")
    let chrom_id = bigbed_read_u32(data, position, endian)
    let chrom_start = bigbed_read_u32(data, position + 4, endian)
    let chrom_end = bigbed_read_u32(data, position + 8, endian)
    let (rest, next) = bigbed_read_c_string(data, position + 12, data.length())
    let words = if rest.length() == 0 {
      []
    } else {
      split_by_char(rest, '\t'.to_int())
    }
    records.push(
      bigbed_record_from_words(file, chrom_id, chrom_start, chrom_end, words),
    )
    position = next
  }
  records
}

///|
pub fn BigBedFile::records(
  self : BigBedFile,
) -> Array[BigBedRecord] raise BigBedError {
  let records : Array[BigBedRecord] = []
  for block in self.blocks {
    for record in bigbed_decode_records(self, block) {
      records.push(record)
    }
  }
  if records.length() != self.record_count {
    bigbed_fail("Decoded BigBed record count does not match its header")
  }
  records
}

///|
pub fn BigBedFile::search(
  self : BigBedFile,
  chromosome : String,
  start? : Int = 0,
  end? : Int = -1,
) -> Array[BigBedRecord] raise BigBedError {
  let chrom_id = bigbed_target_index(self.targets, chromosome)
  if chrom_id < 0 {
    bigbed_fail("Unknown BigBed target '" + chromosome + "'")
  }
  let resolved_end = if end < 0 { self.targets[chrom_id].length } else { end }
  if start < 0 ||
    resolved_end <= start ||
    resolved_end > self.targets[chrom_id].length {
    bigbed_fail("BigBed search requires a valid positive half-open interval")
  }
  let results : Array[BigBedRecord] = []
  for block in self.blocks {
    if bigbed_block_overlaps(block, chrom_id, start, resolved_end) {
      for record in bigbed_decode_records(self, block) {
        if record.overlaps(chromosome, start, resolved_end) {
          results.push(record)
        }
      }
    }
  }
  results
}

///|
pub fn BigBedFile::find_by_name(
  self : BigBedFile,
  name : String,
) -> Array[BigBedRecord] raise BigBedError {
  let results : Array[BigBedRecord] = []
  for record in self.records() {
    if record.name == name {
      results.push(record)
    }
  }
  results
}

///|
pub fn BigBedFile::target(self : BigBedFile, name : String) -> BigBedTarget? {
  let index = bigbed_target_index(self.targets, name)
  if index < 0 {
    None
  } else {
    Some(self.targets[index])
  }
}

///|
pub fn BigBedFile::is_compressed(self : BigBedFile) -> Bool {
  self.header.uncompress_buffer_size > 0
}

///|
pub fn BigBedFile::summary(
  self : BigBedFile,
) -> BigBedSummary raise BigBedError {
  let mut covered_bases = 0
  let mut minimum_span = 0
  let mut maximum_span = 0
  let records = self.records()
  for i in 0..<records.length() {
    let span = records[i].span()
    covered_bases = covered_bases + span
    if i == 0 || span < minimum_span {
      minimum_span = span
    }
    if span > maximum_span {
      maximum_span = span
    }
  }
  BigBedSummary::{
    record_count: self.record_count,
    target_count: self.targets.length(),
    data_block_count: self.blocks.length(),
    covered_bases,
    minimum_span,
    maximum_span,
    compressed: self.is_compressed(),
  }
}

///|
pub fn BigBedFile::to_bed(self : BigBedFile) -> String raise BigBedError {
  let output = StringBuilder::new()
  let bed_columns = self.header.defined_field_count
  for record in self.records() {
    output.write_string(record.chrom)
    output.write_string("\t")
    output.write_string(record.chrom_start.to_string())
    output.write_string("\t")
    output.write_string(record.chrom_end.to_string())
    let rest = bigbed_rest_fields(record, bed_columns)
    for field in rest {
      output.write_string("\t")
      output.write_string(field)
    }
    output.write_string("\n")
  }
  output.to_string()
}

///|
pub fn bigbed_example_data() -> (Array[BigBedTarget], Array[BigBedRecord]) {
  let targets = [
    BigBedTarget::{ name: "chr1", length: 1000000 },
    BigBedTarget::{ name: "chr2", length: 500000 },
  ]
  let records = [
    BigBedRecord::{
      chrom: "chr1",
      chrom_start: 100,
      chrom_end: 280,
      name: "transcript-A",
      score: 960,
      strand: "+",
      thick_start: 120,
      thick_end: 260,
      item_rgb: "40,120,220",
      declared_block_count: 3,
      block_sizes: [50, 40, 30],
      block_starts: [0, 80, 150],
      extra_fields: [],
    },
    BigBedRecord::{
      chrom: "chr1",
      chrom_start: 250,
      chrom_end: 410,
      name: "transcript-B",
      score: 870,
      strand: "-",
      thick_start: 270,
      thick_end: 390,
      item_rgb: "220,80,80",
      declared_block_count: 2,
      block_sizes: [60, 50],
      block_starts: [0, 110],
      extra_fields: [],
    },
    BigBedRecord::{
      chrom: "chr2",
      chrom_start: 1000,
      chrom_end: 1120,
      name: "transcript-C",
      score: 700,
      strand: "+",
      thick_start: 1000,
      thick_end: 1120,
      item_rgb: "0",
      declared_block_count: 2,
      block_sizes: [45, 50],
      block_starts: [0, 70],
      extra_fields: [],
    },
  ]
  (targets, records)
}
