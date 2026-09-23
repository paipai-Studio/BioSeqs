// Bio.Align.bigmaf - indexed UCSC BigMaf multiple alignments.
//
// BigMaf is a bed3+1 BigBed file whose custom lstring field contains a
// semicolon-delimited MAF block. The binary container is implemented by
// bigbed.mbt; this module owns strict MAF block validation and conversion.

///|
pub suberror BigMafError {
  BigMafError(String)
}

///|
pub struct BigMafInsertion {
  left_status : String
  left_count : Int
  right_status : String
  right_count : Int
} derive(Eq, Debug)
pub extend BigMafInsertion with Eq::{equal, not_equal}
pub extend BigMafInsertion with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafComponent {
  source : String
  start : Int
  size : Int
  strand : String
  source_size : Int
  text : String
  insertion : BigMafInsertion?
  quality : String?
} derive(Eq, Debug)
pub extend BigMafComponent with Eq::{equal, not_equal}
pub extend BigMafComponent with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafEmptyComponent {
  source : String
  start : Int
  size : Int
  strand : String
  source_size : Int
  status : String
} derive(Eq, Debug)
pub extend BigMafEmptyComponent with Eq::{equal, not_equal}
pub extend BigMafEmptyComponent with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafBlock {
  score : Double?
  pass_number : Int?
  components : Array[BigMafComponent]
  empty_components : Array[BigMafEmptyComponent]
  comments : Array[String]
} derive(Eq, Debug)
pub extend BigMafBlock with Eq::{equal, not_equal}
pub extend BigMafBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafWriteConfig {
  compress : Bool
  block_size : Int
  items_per_slot : Int
} derive(Eq, Debug)
pub extend BigMafWriteConfig with Eq::{equal, not_equal}
pub extend BigMafWriteConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafSummary {
  alignment_count : Int
  target_count : Int
  component_count : Int
  empty_component_count : Int
  aligned_columns : Int
  covered_reference_bases : Int
  compressed : Bool
} derive(Eq, Debug)
pub extend BigMafSummary with Eq::{equal, not_equal}
pub extend BigMafSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigMafFile {
  reference : String
  bed : BigBedFile
}

///|
fn bigmaf_fail(message : String) -> Unit raise BigMafError {
  raise BigMafError(message)
}

///|
fn bigmaf_is_space(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\n'.to_int() ||
  code == '\r'.to_int()
}

///|
fn bigmaf_words(text : String) -> Array[String] {
  let words : Array[String] = []
  let mut start = -1
  for index in 0..<text.length() {
    if bigmaf_is_space(text.unsafe_get(index).to_int()) {
      if start >= 0 {
        words.push(text[start:index].to_owned())
        start = -1
      }
    } else if start < 0 {
      start = index
    }
  }
  if start >= 0 {
    words.push(text[start:text.length()].to_owned())
  }
  words
}

///|
fn bigmaf_validate_token(
  value : String,
  label : String,
) -> Unit raise BigMafError {
  if value.length() == 0 {
    bigmaf_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if bigmaf_is_space(code) || code == ';'.to_int() || code == 0 {
      bigmaf_fail(label + " contains a forbidden delimiter")
    }
  }
}

///|
fn bigmaf_validate_comment(value : String) -> Unit raise BigMafError {
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == ';'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() ||
      code == 0 {
      bigmaf_fail("BigMaf comments cannot contain line delimiters")
    }
  }
}

///|
fn bigmaf_parse_int(value : String, label : String) -> Int raise BigMafError {
  if value.length() == 0 {
    bigmaf_fail("Missing integer for " + label)
  }
  let mut first = 0
  if value.unsafe_get(0).to_int() == '-'.to_int() {
    if value.length() == 1 {
      bigmaf_fail("Invalid integer for " + label + ": " + value)
    }
    first = 1
  }
  for index in first..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      bigmaf_fail("Invalid integer for " + label + ": " + value)
    }
  }
  parse_int(value)
}

///|
fn bigmaf_parse_double(
  value : String,
  label : String,
) -> Double raise BigMafError {
  match parse_double(value) {
    Some(number) =>
      if number.is_nan() || number.abs() > 1.0e300 {
        raise BigMafError("Non-finite number for " + label)
      } else {
        number
      }
    None => raise BigMafError("Invalid number for " + label + ": " + value)
  }
}

///|
fn bigmaf_non_gap_count(text : String) -> Int {
  let mut count = 0
  for index in 0..<text.length() {
    if text.unsafe_get(index).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn bigmaf_validate_alignment_text(text : String) -> Unit raise BigMafError {
  bigmaf_validate_token(text, "MAF aligned sequence")
  for index in 0..<text.length() {
    let code = text.unsafe_get(index).to_int()
    if code < 33 || code > 126 {
      bigmaf_fail("MAF aligned sequence contains a non-printable character")
    }
  }
}

///|
fn bigmaf_insertion_status(status : String) -> Bool {
  status == "C" ||
  status == "I" ||
  status == "N" ||
  status == "n" ||
  status == "M" ||
  status == "T"
}

///|
fn bigmaf_empty_status(status : String) -> Bool {
  status == "C" || status == "I" || status == "M" || status == "n"
}

///|
fn bigmaf_source_parts(source : String) -> (String, String) raise BigMafError {
  let mut dot = -1
  for index in 0..<source.length() {
    if source.unsafe_get(index).to_int() == '.'.to_int() {
      dot = index
      break
    }
  }
  if dot <= 0 || dot + 1 >= source.length() {
    bigmaf_fail(
      "BigMaf reference component source must be '<reference>.<chromosome>'",
    )
  }
  (source[0:dot].to_owned(), source[dot + 1:source.length()].to_owned())
}

///|
pub fn BigMafInsertion::create(
  left_status : String,
  left_count : Int,
  right_status : String,
  right_count : Int,
) -> BigMafInsertion raise BigMafError {
  if !bigmaf_insertion_status(left_status) ||
    !bigmaf_insertion_status(right_status) {
    bigmaf_fail("MAF insertion status must be C, I, N, n, M, or T")
  }
  if left_count < 0 || right_count < 0 {
    bigmaf_fail("MAF insertion counts must not be negative")
  }
  BigMafInsertion::{ left_status, left_count, right_status, right_count, }
}

///|
pub fn BigMafComponent::create(
  source : String,
  start : Int,
  size : Int,
  strand : String,
  source_size : Int,
  text : String,
  insertion? : BigMafInsertion? = None,
  quality? : String? = None,
) -> BigMafComponent raise BigMafError {
  bigmaf_validate_token(source, "MAF source")
  if start < 0 || size <= 0 || source_size <= 0 || start + size > source_size {
    bigmaf_fail("MAF component coordinates exceed the source sequence")
  }
  if strand != "+" && strand != "-" {
    bigmaf_fail("MAF component strand must be '+' or '-'")
  }
  bigmaf_validate_alignment_text(text)
  if bigmaf_non_gap_count(text) != size {
    bigmaf_fail("MAF component size does not match its non-gap sequence length")
  }
  match quality {
    Some(value) => {
      if value.length() != text.length() {
        bigmaf_fail("MAF quality string must match the aligned sequence width")
      }
      for index in 0..<value.length() {
        let quality_code = value.unsafe_get(index).to_int()
        let sequence_gap = text.unsafe_get(index).to_int() == '-'.to_int()
        let quality_gap = quality_code == '-'.to_int()
        let valid_score = (
            quality_code >= '0'.to_int() && quality_code <= '9'.to_int()
          ) ||
          quality_code == 'F'.to_int()
        if sequence_gap != quality_gap || (!quality_gap && !valid_score) {
          bigmaf_fail(
            "MAF quality must use 0-9/F and preserve aligned gap columns",
          )
        }
      }
    }
    None => ()
  }
  match insertion {
    Some(value) =>
      ignore(
        BigMafInsertion::create(
          value.left_status,
          value.left_count,
          value.right_status,
          value.right_count,
        ),
      )
    None => ()
  }
  BigMafComponent::{
    source,
    start,
    size,
    strand,
    source_size,
    text,
    insertion,
    quality,
  }
}

///|
pub fn BigMafComponent::with_insertion(
  self : BigMafComponent,
  insertion : BigMafInsertion,
) -> BigMafComponent raise BigMafError {
  BigMafComponent::create(
    self.source,
    self.start,
    self.size,
    self.strand,
    self.source_size,
    self.text,
    insertion=Some(insertion),
    quality=self.quality,
  )
}

///|
pub fn BigMafComponent::with_quality(
  self : BigMafComponent,
  quality : String,
) -> BigMafComponent raise BigMafError {
  BigMafComponent::create(
    self.source,
    self.start,
    self.size,
    self.strand,
    self.source_size,
    self.text,
    insertion=self.insertion,
    quality=Some(quality),
  )
}

///|
pub fn BigMafEmptyComponent::create(
  source : String,
  start : Int,
  size : Int,
  strand : String,
  source_size : Int,
  status : String,
) -> BigMafEmptyComponent raise BigMafError {
  bigmaf_validate_token(source, "MAF empty component source")
  if start < 0 || size < 0 || source_size <= 0 || start + size > source_size {
    bigmaf_fail("MAF empty component coordinates exceed the source sequence")
  }
  if strand != "+" && strand != "-" {
    bigmaf_fail("MAF empty component strand must be '+' or '-'")
  }
  if !bigmaf_empty_status(status) {
    bigmaf_fail("MAF empty component status must be C, I, M, or n")
  }
  BigMafEmptyComponent::{ source, start, size, strand, source_size, status, }
}

///|
pub fn BigMafBlock::create(
  components : Array[BigMafComponent],
  score? : Double? = None,
  pass_number? : Int? = None,
  empty_components? : Array[BigMafEmptyComponent] = [],
  comments? : Array[String] = [],
) -> BigMafBlock raise BigMafError {
  if components.length() == 0 {
    bigmaf_fail("BigMaf blocks require at least one aligned component")
  }
  match score {
    Some(value) =>
      if value.is_nan() || value.abs() > 1.0e300 {
        bigmaf_fail("BigMaf block score must be finite")
      }
    None => ()
  }
  match pass_number {
    Some(value) =>
      if value <= 0 {
        bigmaf_fail("MAF pass annotation must be positive")
      }
    None => ()
  }
  let width = components[0].text.length()
  for component in components {
    ignore(
      BigMafComponent::create(
        component.source,
        component.start,
        component.size,
        component.strand,
        component.source_size,
        component.text,
        insertion=component.insertion,
        quality=component.quality,
      ),
    )
    if component.text.length() != width {
      bigmaf_fail("All MAF components in a block must have equal width")
    }
  }
  for empty in empty_components {
    ignore(
      BigMafEmptyComponent::create(
        empty.source,
        empty.start,
        empty.size,
        empty.strand,
        empty.source_size,
        empty.status,
      ),
    )
  }
  for comment in comments {
    bigmaf_validate_comment(comment)
  }
  BigMafBlock::{
    score,
    pass_number,
    components: components.copy(),
    empty_components: empty_components.copy(),
    comments: comments.copy(),
  }
}

///|
pub fn BigMafWriteConfig::create(
  compress? : Bool = true,
  block_size? : Int = 256,
  items_per_slot? : Int = 512,
) -> BigMafWriteConfig raise BigMafError {
  if block_size < 2 || block_size > 65535 {
    bigmaf_fail("BigMaf blockSize must be between 2 and 65535")
  }
  if items_per_slot <= 0 || items_per_slot > 65535 {
    bigmaf_fail("BigMaf itemsPerSlot must be between 1 and 65535")
  }
  BigMafWriteConfig::{ compress, block_size, items_per_slot, }
}

///|
pub fn BigMafWriteConfig::default() -> BigMafWriteConfig {
  BigMafWriteConfig::{ compress: true, block_size: 256, items_per_slot: 512, }
}

///|
pub fn BigMafComponent::sequence(self : BigMafComponent) -> String {
  let output = StringBuilder()
  for index in 0..<self.text.length() {
    if self.text.unsafe_get(index).to_int() != '-'.to_int() {
      output.write_string(self.text[index:index + 1].to_owned())
    }
  }
  output.to_string()
}

///|
pub fn BigMafComponent::forward_sequence(self : BigMafComponent) -> String {
  let sequence = self.sequence()
  if self.strand == "+" {
    sequence
  } else {
    reverse_complement(sequence)
  }
}

///|
pub fn BigMafComponent::forward_interval(self : BigMafComponent) -> (Int, Int) {
  if self.strand == "+" {
    (self.start, self.start + self.size)
  } else {
    (self.source_size - self.start - self.size, self.source_size - self.start)
  }
}

///|
pub fn BigMafEmptyComponent::forward_interval(
  self : BigMafEmptyComponent,
) -> (Int, Int) {
  if self.strand == "+" {
    (self.start, self.start + self.size)
  } else {
    (self.source_size - self.start - self.size, self.source_size - self.start)
  }
}

///|
pub fn BigMafComponent::column_to_source(
  self : BigMafComponent,
  column : Int,
) -> Int? raise BigMafError {
  if column < 0 || column >= self.text.length() {
    bigmaf_fail("MAF alignment column is out of range")
  }
  if self.text.unsafe_get(column).to_int() == '-'.to_int() {
    return None
  }
  let mut offset = 0
  for index in 0..<column {
    if self.text.unsafe_get(index).to_int() != '-'.to_int() {
      offset = offset + 1
    }
  }
  if self.strand == "+" {
    Some(self.start + offset)
  } else {
    Some(self.source_size - self.start - offset - 1)
  }
}

///|
pub fn BigMafComponent::source_to_column(
  self : BigMafComponent,
  position : Int,
) -> Int? {
  let (forward_start, forward_end) = self.forward_interval()
  if position < forward_start || position >= forward_end {
    return None
  }
  let wanted = if self.strand == "+" {
    position - self.start
  } else {
    self.source_size - self.start - position - 1
  }
  let mut offset = 0
  for column in 0..<self.text.length() {
    if self.text.unsafe_get(column).to_int() != '-'.to_int() {
      if offset == wanted {
        return Some(column)
      }
      offset = offset + 1
    }
  }
  None
}

///|
pub fn BigMafBlock::aligned_columns(self : BigMafBlock) -> Int {
  self.components[0].text.length()
}

///|
pub fn BigMafBlock::reference_component(self : BigMafBlock) -> BigMafComponent {
  self.components[0]
}

///|
pub fn BigMafBlock::component(
  self : BigMafBlock,
  source : String,
) -> BigMafComponent? {
  for component in self.components {
    if component.source == source {
      return Some(component)
    }
  }
  None
}

///|
pub fn BigMafBlock::sources(self : BigMafBlock) -> Array[String] {
  let sources : Array[String] = []
  for component in self.components {
    sources.push(component.source)
  }
  sources
}

///|
pub fn BigMafBlock::map_reference_to(
  self : BigMafBlock,
  reference_position : Int,
  source : String,
) -> Int? {
  match self.components[0].source_to_column(reference_position) {
    Some(column) =>
      match self.component(source) {
        Some(component) =>
          component.column_to_source(column) catch {
            _ => None
          }
        None => None
      }
    None => None
  }
}

///|
pub fn BigMafBlock::pairwise_identity(self : BigMafBlock) -> Double {
  let mut matches = 0
  let mut comparisons = 0
  for column in 0..<self.aligned_columns() {
    for left in 0..<self.components.length() {
      let left_code = self.components[left].text.unsafe_get(column).to_int()
      if left_code == '-'.to_int() {
        continue
      }
      for right in (left + 1)..<self.components.length() {
        let right_code = self.components[right].text.unsafe_get(column).to_int()
        if right_code != '-'.to_int() {
          comparisons = comparisons + 1
          if left_code == right_code {
            matches = matches + 1
          }
        }
      }
    }
  }
  if comparisons == 0 {
    0.0
  } else {
    matches.to_double() / comparisons.to_double()
  }
}

///|
fn bigmaf_write_component(
  output : StringBuilder,
  component : BigMafComponent,
  delimiter : String,
) -> Unit {
  output.write_string(
    "s " +
    component.source +
    " " +
    component.start.to_string() +
    " " +
    component.size.to_string() +
    " " +
    component.strand +
    " " +
    component.source_size.to_string() +
    " " +
    component.text +
    delimiter,
  )
  match component.insertion {
    Some(insertion) =>
      output.write_string(
        "i " +
        component.source +
        " " +
        insertion.left_status +
        " " +
        insertion.left_count.to_string() +
        " " +
        insertion.right_status +
        " " +
        insertion.right_count.to_string() +
        delimiter,
      )
    None => ()
  }
  match component.quality {
    Some(quality) =>
      output.write_string("q " + component.source + " " + quality + delimiter)
    None => ()
  }
}

///|
fn bigmaf_format_block(block : BigMafBlock, delimiter : String) -> String {
  let output = StringBuilder()
  for comment in block.comments {
    output.write_string("# " + comment + delimiter)
  }
  output.write_string("a")
  match block.score {
    Some(score) => output.write_string(" score=" + score.to_string())
    None => ()
  }
  match block.pass_number {
    Some(pass_number) => output.write_string(" pass=" + pass_number.to_string())
    None => ()
  }
  output.write_string(delimiter)
  for component in block.components {
    bigmaf_write_component(output, component, delimiter)
  }
  for empty in block.empty_components {
    output.write_string(
      "e " +
      empty.source +
      " " +
      empty.start.to_string() +
      " " +
      empty.size.to_string() +
      " " +
      empty.strand +
      " " +
      empty.source_size.to_string() +
      " " +
      empty.status +
      delimiter,
    )
  }
  output.to_string()
}

///|
pub fn BigMafBlock::to_embedded_maf(self : BigMafBlock) -> String {
  bigmaf_format_block(self, ";")
}

///|
pub fn BigMafBlock::to_maf(self : BigMafBlock) -> String {
  bigmaf_format_block(self, "\n") + "\n"
}

///|
fn bigmaf_parse_score_line(
  words : Array[String],
) -> (Double?, Int?) raise BigMafError {
  let mut score : Double? = None
  let mut pass_number : Int? = None
  for index in 1..<words.length() {
    let field = words[index]
    let mut equal = -1
    for position in 0..<field.length() {
      if field.unsafe_get(position).to_int() == '='.to_int() {
        equal = position
        break
      }
    }
    if equal <= 0 || equal + 1 >= field.length() {
      bigmaf_fail("Malformed MAF alignment annotation '" + field + "'")
    }
    let key = field[0:equal].to_owned()
    let value = field[equal + 1:field.length()].to_owned()
    if key == "score" {
      if score is Some(_) {
        bigmaf_fail("Duplicate MAF score annotation")
      }
      score = Some(bigmaf_parse_double(value, "MAF score"))
    } else if key == "pass" {
      if pass_number is Some(_) {
        bigmaf_fail("Duplicate MAF pass annotation")
      }
      let parsed = bigmaf_parse_int(value, "MAF pass")
      if parsed <= 0 {
        bigmaf_fail("MAF pass annotation must be positive")
      }
      pass_number = Some(parsed)
    } else {
      bigmaf_fail("Unknown MAF alignment annotation '" + key + "'")
    }
  }
  (score, pass_number)
}

///|
pub fn bigmaf_parse_block(text : String) -> BigMafBlock raise BigMafError {
  if text.length() == 0 {
    bigmaf_fail("BigMaf block text must not be empty")
  }
  let lines = text.split(";").to_array()
  let components : Array[BigMafComponent] = []
  let empty_components : Array[BigMafEmptyComponent] = []
  let comments : Array[String] = []
  let mut score : Double? = None
  let mut pass_number : Int? = None
  let mut saw_alignment = false
  let mut last_component = -1
  for raw in lines {
    let line = trim(raw.to_owned())
    if line.length() == 0 {
      continue
    }
    if line.unsafe_get(0).to_int() == '#'.to_int() {
      comments.push(trim(line[1:line.length()].to_owned()))
      continue
    }
    let words = bigmaf_words(line)
    if words.length() == 0 {
      continue
    }
    if words[0] == "a" {
      if saw_alignment {
        bigmaf_fail("BigMaf payload contains more than one MAF block")
      }
      let (found_score, found_pass) = bigmaf_parse_score_line(words)
      score = found_score
      pass_number = found_pass
      saw_alignment = true
      last_component = -1
    } else if words[0] == "s" {
      if !saw_alignment {
        bigmaf_fail("MAF sequence line appears before the alignment line")
      }
      if words.length() != 7 {
        bigmaf_fail("MAF sequence lines must contain seven fields")
      }
      let component = BigMafComponent::create(
        words[1],
        bigmaf_parse_int(words[2], "MAF component start"),
        bigmaf_parse_int(words[3], "MAF component size"),
        words[4],
        bigmaf_parse_int(words[5], "MAF source size"),
        words[6],
      )
      components.push(component)
      last_component = components.length() - 1
    } else if words[0] == "i" {
      if last_component < 0 || words.length() != 6 {
        bigmaf_fail(
          "MAF insertion lines must follow an s line and have six fields",
        )
      }
      if words[1] != components[last_component].source {
        bigmaf_fail("MAF insertion source does not match the preceding s line")
      }
      if components[last_component].insertion is Some(_) {
        bigmaf_fail("Duplicate MAF insertion line for one component")
      }
      let insertion = BigMafInsertion::create(
        words[2],
        bigmaf_parse_int(words[3], "MAF left insertion count"),
        words[4],
        bigmaf_parse_int(words[5], "MAF right insertion count"),
      )
      components[last_component] = components[last_component].with_insertion(
        insertion,
      )
    } else if words[0] == "q" {
      if last_component < 0 || words.length() != 3 {
        bigmaf_fail(
          "MAF quality lines must follow an s line and have three fields",
        )
      }
      if words[1] != components[last_component].source {
        bigmaf_fail("MAF quality source does not match the preceding s line")
      }
      if components[last_component].quality is Some(_) {
        bigmaf_fail("Duplicate MAF quality line for one component")
      }
      components[last_component] = components[last_component].with_quality(
        words[2],
      )
    } else if words[0] == "e" {
      if !saw_alignment || words.length() != 7 {
        bigmaf_fail("MAF empty component lines must contain seven fields")
      }
      empty_components.push(
        BigMafEmptyComponent::create(
          words[1],
          bigmaf_parse_int(words[2], "MAF empty component start"),
          bigmaf_parse_int(words[3], "MAF empty component size"),
          words[4],
          bigmaf_parse_int(words[5], "MAF empty source size"),
          words[6],
        ),
      )
      last_component = -1
    } else {
      bigmaf_fail("Unexpected MAF line type '" + words[0] + "'")
    }
  }
  if !saw_alignment {
    bigmaf_fail("BigMaf payload is missing its MAF alignment line")
  }
  BigMafBlock::create(
    components,
    score~,
    pass_number~,
    empty_components~,
    comments~,
  )
}

///|
pub fn bigmaf_schema() -> BigBedSchema {
  BigBedSchema::{
    name: "bedMaf",
    comment: "Bed3 with MAF block",
    fields: [
      BigBedField::{
        as_type: "string",
        name: "chrom",
        comment: "Reference sequence chromosome or scaffold",
      },
      BigBedField::{
        as_type: "uint",
        name: "chromStart",
        comment: "Start position in chromosome",
      },
      BigBedField::{
        as_type: "uint",
        name: "chromEnd",
        comment: "End position in chromosome",
      },
      BigBedField::{
        as_type: "lstring",
        name: "mafBlock",
        comment: "MAF block",
      },
    ],
  }
}

///|
fn bigmaf_target_index(
  targets : Array[BigBedTarget],
  chromosome : String,
) -> Int {
  for index in 0..<targets.length() {
    if targets[index].name == chromosome {
      return index
    }
  }
  -1
}

///|
fn bigmaf_wrap_bigbed(
  operation : () -> Array[Int] raise BigBedError,
) -> Array[Int] raise BigMafError {
  operation() catch {
    BigBedError(message) => raise BigMafError("BigBed: " + message)
  }
}

///|
pub fn bigmaf_write(
  reference : String,
  targets : Array[BigBedTarget],
  blocks : Array[BigMafBlock],
  config? : BigMafWriteConfig = BigMafWriteConfig::default(),
) -> Array[Int] raise BigMafError {
  bigmaf_validate_token(reference, "BigMaf reference")
  for index in 0..<reference.length() {
    if reference.unsafe_get(index).to_int() == '.'.to_int() {
      bigmaf_fail("BigMaf reference name cannot contain '.'")
    }
  }
  if targets.length() == 0 {
    bigmaf_fail("BigMaf requires at least one reference target")
  }
  for index in 0..<targets.length() {
    let target = targets[index]
    ignore(BigBedTarget::create(target.name, target.length)) catch {
      BigBedError(message) => bigmaf_fail("BigBed: " + message)
    }
    for previous in 0..<index {
      if targets[previous].name == target.name {
        bigmaf_fail("Duplicate BigMaf target '" + target.name + "'")
      }
    }
  }
  let records : Array[BigBedRecord] = []
  let target_order : Array[Int] = []
  for block in blocks {
    let validated = BigMafBlock::create(
      block.components,
      score=block.score,
      pass_number=block.pass_number,
      empty_components=block.empty_components,
      comments=block.comments,
    )
    let component = validated.components[0]
    let (block_reference, chromosome) = bigmaf_source_parts(component.source)
    if block_reference != reference {
      bigmaf_fail(
        "Reference component '" +
        component.source +
        "' does not match reference '" +
        reference +
        "'",
      )
    }
    if component.strand != "+" {
      bigmaf_fail("BigMaf reference components must use the forward strand")
    }
    let target_index = bigmaf_target_index(targets, chromosome)
    if target_index < 0 {
      bigmaf_fail("BigMaf block refers to unknown target '" + chromosome + "'")
    }
    if component.source_size != targets[target_index].length {
      bigmaf_fail(
        "BigMaf reference component source size does not match its target",
      )
    }
    let record = BigBedRecord::create(
      chromosome,
      component.start,
      component.start + component.size,
      extra_fields=[validated.to_embedded_maf()],
    ) catch {
      BigBedError(message) => {
        bigmaf_fail("BigBed: " + message)
        BigBedRecord::{
          chrom: "",
          chrom_start: 0,
          chrom_end: 1,
          name: ".",
          score: 0,
          strand: ".",
          thick_start: 0,
          thick_end: 1,
          item_rgb: "0",
          declared_block_count: 0,
          block_sizes: [],
          block_starts: [],
          extra_fields: [],
        }
      }
    }
    records.push(record)
    target_order.push(target_index)
  }
  let order : Array[Int] = []
  for index in 0..<records.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if target_order[left] < target_order[right] {
      -1
    } else if target_order[left] > target_order[right] {
      1
    } else if records[left].chrom_start < records[right].chrom_start {
      -1
    } else if records[left].chrom_start > records[right].chrom_start {
      1
    } else if records[left].chrom_end < records[right].chrom_end {
      -1
    } else if records[left].chrom_end > records[right].chrom_end {
      1
    } else {
      0
    }
  })
  let sorted : Array[BigBedRecord] = []
  for index in order {
    sorted.push(records[index])
  }
  let field = BigBedField::{
    as_type: "lstring",
    name: "mafBlock",
    comment: "MAF block",
  }
  let bed_config = BigBedWriteConfig::create(
    bed_columns=3,
    items_per_slot=config.items_per_slot,
    block_size=config.block_size,
    compress=config.compress,
    schema_name="bedMaf",
    schema_comment="Bed3 with MAF block",
    custom_fields=[field],
  ) catch {
    BigBedError(message) => {
      bigmaf_fail("BigBed: " + message)
      BigBedWriteConfig::default()
    }
  }
  bigmaf_wrap_bigbed(fn() -> Array[Int] raise BigBedError {
    bigbed_write(targets, sorted, config=bed_config)
  })
}

///|
fn bigmaf_validate_schema(bed : BigBedFile) -> Unit raise BigMafError {
  if bed.header.defined_field_count != 3 ||
    bed.header.field_count != 4 ||
    bed.schema.name != "bedMaf" ||
    bed.schema.fields.length() != 4 ||
    bed.schema.fields[3].name != "mafBlock" ||
    bed.schema.fields[3].as_type != "lstring" {
    bigmaf_fail(
      "BigMaf requires the standard bedMaf bed3+1 AutoSQL declaration",
    )
  }
}

///|
fn bigmaf_decode_record(
  bed : BigBedFile,
  record : BigBedRecord,
  expected_reference : String,
) -> (String, BigMafBlock) raise BigMafError {
  if record.extra_fields.length() != 1 {
    bigmaf_fail("BigMaf binary records require one mafBlock field")
  }
  let block = bigmaf_parse_block(record.extra_fields[0])
  let component = block.components[0]
  let (reference, chromosome) = bigmaf_source_parts(component.source)
  if expected_reference.length() > 0 && reference != expected_reference {
    bigmaf_fail("BigMaf records use inconsistent reference prefixes")
  }
  if component.strand != "+" {
    bigmaf_fail("BigMaf reference components must use the forward strand")
  }
  if chromosome != record.chrom ||
    component.start != record.chrom_start ||
    component.start + component.size != record.chrom_end {
    bigmaf_fail("BigMaf BED coordinates disagree with the embedded MAF block")
  }
  match bed.target(chromosome) {
    Some(target) =>
      if component.source_size != target.length {
        bigmaf_fail(
          "BigMaf reference component source size disagrees with its target",
        )
      }
    None => bigmaf_fail("BigMaf record refers to an unknown target")
  }
  (reference, block)
}

///|
pub fn bigmaf_parse(data : Array[Int]) -> BigMafFile raise BigMafError {
  let bed = bigbed_parse(data) catch {
    BigBedError(message) => {
      bigmaf_fail("BigBed: " + message)
      BigBedFile::{
        header: BigBedHeader::{
          little_endian: true,
          version: 4,
          zoom_levels: 0,
          chromosome_tree_offset: 0,
          full_data_offset: 0,
          full_index_offset: 0,
          field_count: 0,
          defined_field_count: 0,
          auto_sql_offset: 0,
          total_summary_offset: 0,
          uncompress_buffer_size: 0,
          extra_indices_offset: 0,
        },
        schema: bigmaf_schema(),
        targets: [],
        source: [],
        blocks: [],
        record_count: 0,
      }
    }
  }
  bigmaf_validate_schema(bed)
  let records = bed.records() catch {
    BigBedError(message) => {
      bigmaf_fail("BigBed: " + message)
      []
    }
  }
  let mut reference = ""
  for record in records {
    let (found, _) = bigmaf_decode_record(bed, record, reference)
    if reference.length() == 0 {
      reference = found
    }
  }
  BigMafFile::{ reference, bed, }
}

///|
pub fn BigMafFile::blocks(
  self : BigMafFile,
) -> Array[BigMafBlock] raise BigMafError {
  let records = self.bed.records() catch {
    BigBedError(message) => {
      bigmaf_fail("BigBed: " + message)
      []
    }
  }
  let blocks : Array[BigMafBlock] = []
  for record in records {
    let (_, block) = bigmaf_decode_record(self.bed, record, self.reference)
    blocks.push(block)
  }
  blocks
}

///|
fn bigmaf_query_chromosome(
  file : BigMafFile,
  chromosome : String,
) -> String raise BigMafError {
  if file.reference.length() == 0 {
    return chromosome
  }
  let prefix = file.reference + "."
  if chromosome.has_prefix(prefix) {
    let suffix = chromosome[prefix.length():chromosome.length()].to_owned()
    if suffix.length() == 0 {
      bigmaf_fail("BigMaf query chromosome must not be empty")
    }
    suffix
  } else {
    chromosome
  }
}

///|
pub fn BigMafFile::search(
  self : BigMafFile,
  chromosome : String,
  start? : Int = 0,
  end? : Int = -1,
) -> Array[BigMafBlock] raise BigMafError {
  let resolved = bigmaf_query_chromosome(self, chromosome)
  let records = self.bed.search(resolved, start~, end~) catch {
    BigBedError(message) => {
      bigmaf_fail("BigBed: " + message)
      []
    }
  }
  let blocks : Array[BigMafBlock] = []
  for record in records {
    let (_, block) = bigmaf_decode_record(self.bed, record, self.reference)
    blocks.push(block)
  }
  blocks
}

///|
pub fn BigMafFile::target(
  self : BigMafFile,
  chromosome : String,
) -> BigBedTarget? {
  let resolved = bigmaf_query_chromosome(self, chromosome) catch {
    _ => return None
  }
  self.bed.target(resolved)
}

///|
pub fn BigMafFile::is_compressed(self : BigMafFile) -> Bool {
  self.bed.is_compressed()
}

///|
pub fn BigMafFile::to_maf(self : BigMafFile) -> String raise BigMafError {
  let output = StringBuilder()
  output.write_string("##maf version=1\n\n")
  for block in self.blocks() {
    output.write_string(block.to_maf())
  }
  output.to_string()
}

///|
pub fn BigMafFile::summary(
  self : BigMafFile,
) -> BigMafSummary raise BigMafError {
  let blocks = self.blocks()
  let mut component_count = 0
  let mut empty_component_count = 0
  let mut aligned_columns = 0
  let mut covered_reference_bases = 0
  for block in blocks {
    component_count = component_count + block.components.length()
    empty_component_count = empty_component_count +
      block.empty_components.length()
    aligned_columns = aligned_columns + block.aligned_columns()
    covered_reference_bases = covered_reference_bases + block.components[0].size
  }
  BigMafSummary::{
    alignment_count: blocks.length(),
    target_count: self.bed.targets.length(),
    component_count,
    empty_component_count,
    aligned_columns,
    covered_reference_bases,
    compressed: self.is_compressed(),
  }
}

///|
pub fn bigmaf_example_data() -> (
  String,
  Array[BigBedTarget],
  Array[BigMafBlock],
) raise BigMafError {
  let reference = "hg38"
  let targets = [
    BigBedTarget::create("chr7", 1000) catch {
      BigBedError(message) => {
        bigmaf_fail("BigBed: " + message)
        BigBedTarget::{ name: "chr7", length: 1000, }
      }
    },
    BigBedTarget::create("chr12", 800) catch {
      BigBedError(message) => {
        bigmaf_fail("BigBed: " + message)
        BigBedTarget::{ name: "chr12", length: 800, }
      }
    },
  ]
  let reference_one = BigMafComponent::create(
    "hg38.chr7",
    100,
    8,
    "+",
    1000,
    "ACGT-ACGT",
    quality=Some("9999-FFFF"),
  )
  let chimp_one = BigMafComponent::create(
    "panTro6.chr6", 200, 8, "+", 900, "ACGT-ACGT",
  )
  let mouse_one = BigMafComponent::create(
    "mm39.chr5",
    50,
    8,
    "-",
    700,
    "AC-TTACGT",
    insertion=Some(BigMafInsertion::create("C", 0, "I", 12)),
  )
  let empty = BigMafEmptyComponent::create(
    "canFam6.chr6", 300, 20, "+", 850, "M",
  )
  let first = BigMafBlock::create(
    [reference_one, chimp_one, mouse_one],
    score=Some(23262.0),
    pass_number=Some(1),
    empty_components=[empty],
    comments=["conserved reference block"],
  )
  let second = BigMafBlock::create(
    [
      BigMafComponent::create("hg38.chr7", 300, 6, "+", 1000, "TTAAGA"),
      BigMafComponent::create("panTro6.chr6", 410, 6, "+", 900, "TTAAGA"),
      BigMafComponent::create("mm39.chr5", 120, 5, "+", 700, "TT-AGA"),
    ],
    score=Some(5062.0),
  )
  let third = BigMafBlock::create(
    [
      BigMafComponent::create("hg38.chr12", 40, 7, "+", 800, "GCA-GCTG"),
      BigMafComponent::create("panTro6.chr10", 70, 7, "+", 600, "GCA-GCTG"),
    ],
    score=Some(6636.0),
  )
  (reference, targets, [first, second, third])
}
