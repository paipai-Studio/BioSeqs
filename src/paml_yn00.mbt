// Biopython-compatible Bio.Phylo.PAML.yn00 support.
//
// This module parses YN00 control files and result text. It deliberately does
// not launch the external PAML executable.

///|
/// Error raised for malformed YN00 configuration, output, or queries.
pub suberror Yn00Error {
  Yn00Error(String)
}

///|
/// One validated YN00 option retained in canonical text form.
pub struct Yn00Option {
  name : String
  value : String
} derive(Eq, Debug)

///|
/// A YN00 control file.
pub struct Yn00Control {
  sequence_file : String
  output_file : String
  options : Array[Yn00Option]
} derive(Eq, Debug)

///|
/// Nei-Gojobori (1986) estimates for one sequence pair.
pub struct Yn00Ng86Estimate {
  omega : Double
  dn : Double
  ds : Double
} derive(Eq, Debug)

///|
/// Yang-Nielsen (2000) estimates and standard errors for one sequence pair.
pub struct Yn00YangNielsenEstimate {
  synonymous_sites : Double
  nonsynonymous_sites : Double
  time : Double
  kappa : Double
  omega : Double
  dn : Double
  dn_standard_error : Double
  ds : Double
  ds_standard_error : Double
} derive(Eq, Debug)

///|
/// An estimate from LWL85, modified LWL85, or LPB93.
///
/// Undefined PAML values such as `nan`, `-nan`, and `inf` are represented by
/// `None` instead of leaking non-finite values into downstream calculations.
pub struct Yn00CountingEstimate {
  method_name : String
  ds : Double?
  dn : Double?
  omega : Double?
  synonymous_sites : Double?
  nonsynonymous_sites : Double?
  rho : Double?
} derive(Eq, Debug)

///|
/// All five YN00 estimates for one unordered sequence pair.
pub struct Yn00PairResult {
  first : String
  second : String
  ng86 : Yn00Ng86Estimate
  yang_nielsen : Yn00YangNielsenEstimate
  lwl85 : Yn00CountingEstimate
  lwl85_modified : Yn00CountingEstimate
  lpb93 : Yn00CountingEstimate
} derive(Eq, Debug)

///|
/// A symmetric matrix for one method/statistic combination.
pub struct Yn00DistanceMatrix {
  names : Array[String]
  values : Array[Array[Double?]]
} derive(Eq, Debug)

///|
/// Parsed YN00 pairwise results.
pub struct Yn00Results {
  alignment_file : String
  sequence_count : Int
  codon_count : Int
  sequence_names : Array[String]
  pairs : Array[Yn00PairResult]
} derive(Eq, Debug)

///|
priv struct Yn00PairBuilder {
  first_index : Int
  second_index : Int
  mut comparison_seen : Bool
  mut ng86 : Yn00Ng86Estimate?
  mut yang_nielsen : Yn00YangNielsenEstimate?
  mut lwl85 : Yn00CountingEstimate?
  mut lwl85_modified : Yn00CountingEstimate?
  mut lpb93 : Yn00CountingEstimate?
}

///|
fn yn00_fail(message : String) -> Unit raise Yn00Error {
  raise Yn00Error(message)
}

///|
fn yn00_parse_int(text : String, label : String) -> Int raise Yn00Error {
  codeml_parse_int(text, label) catch {
    CodemlError(_) => {
      yn00_fail(label + " must be an integer")
      0
    }
  }
}

///|
fn yn00_parse_double(text : String, label : String) -> Double raise Yn00Error {
  codeml_parse_double(text, label) catch {
    CodemlError(_) => {
      yn00_fail(label + " must be a finite number")
      0.0
    }
  }
}

///|
fn yn00_copy_options(values : Array[Yn00Option]) -> Array[Yn00Option] {
  let copied : Array[Yn00Option] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn yn00_copy_strings(values : Array[String]) -> Array[String] {
  let copied : Array[String] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn yn00_copy_pairs(values : Array[Yn00PairResult]) -> Array[Yn00PairResult] {
  let copied : Array[Yn00PairResult] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn yn00_copy_optional_matrix(
  values : Array[Array[Double?]],
) -> Array[Array[Double?]] {
  let copied : Array[Array[Double?]] = []
  for row in values {
    let copied_row : Array[Double?] = []
    for value in row {
      copied_row.push(value)
    }
    copied.push(copied_row)
  }
  copied
}

///|
fn yn00_allowed_option(name : String) -> Bool {
  match name {
    "verbose" | "icode" | "weighting" | "commonf3x4" | "ndata" => true
    _ => false
  }
}

///|
fn yn00_validate_option(name : String, value : Int) -> Unit raise Yn00Error {
  match name {
    "verbose" | "weighting" | "commonf3x4" =>
      if value < 0 || value > 1 {
        yn00_fail("YN00 " + name + " must be 0 or 1")
      }
    "icode" =>
      if value < 0 || value > 10 {
        yn00_fail("YN00 icode must be between 0 and 10")
      }
    "ndata" => if value <= 0 { yn00_fail("YN00 ndata must be positive") }
    _ => ()
  }
}

///|
/// Parse and strictly validate YN00 control text.
pub fn yn00_parse_control(text : String) -> Yn00Control raise Yn00Error {
  let mut sequence_file = ""
  let mut output_file = ""
  let options : Array[Yn00Option] = []
  let seen : Array[String] = []
  let lines = codeml_lines(text)
  for line_index in 0..<lines.length() {
    let line = lines[line_index]
    let comment = codeml_index_of_char(line, '*'.to_int(), 0)
    let uncommented = codeml_trim(
      if comment >= 0 {
        line[0:comment].to_owned()
      } else {
        line
      },
    )
    if uncommented.length() == 0 {
      continue
    }
    let equals = codeml_index_of_char(uncommented, '='.to_int(), 0)
    if equals < 0 ||
      codeml_index_of_char(uncommented, '='.to_int(), equals + 1) >= 0 {
      yn00_fail("Malformed YN00 control line " + (line_index + 1).to_string())
    }
    let name = codeml_trim(uncommented[0:equals].to_owned())
    let value = codeml_trim(
      uncommented[equals + 1:uncommented.length()].to_owned(),
    )
    if name.length() == 0 || value.length() == 0 {
      yn00_fail("Malformed YN00 control line " + (line_index + 1).to_string())
    }
    for previous in seen {
      if previous == name {
        yn00_fail("Duplicate YN00 control key: " + name)
      }
    }
    seen.push(name)
    match name {
      "seqfile" => sequence_file = value
      "outfile" => output_file = value
      _ => {
        if !yn00_allowed_option(name) {
          yn00_fail("Unknown YN00 option: " + name)
        }
        let parsed = yn00_parse_int(value, "YN00 " + name)
        yn00_validate_option(name, parsed)
        options.push(Yn00Option::{ name, value: parsed.to_string(), })
      }
    }
  }
  if sequence_file.length() == 0 {
    yn00_fail("YN00 control is missing seqfile")
  }
  if output_file.length() == 0 {
    yn00_fail("YN00 control is missing outfile")
  }
  Yn00Control::{ sequence_file, output_file, options, }
}

///|
/// Write canonical YN00 control text.
pub fn yn00_write_control(control : Yn00Control) -> String {
  let output : Array[String] = [
    "seqfile = " + control.sequence_file,
    "outfile = " + control.output_file,
  ]
  for option in control.options {
    output.push(option.name + " = " + option.value)
  }
  output.join("\n") + "\n"
}

///|
pub fn Yn00Option::name(self : Yn00Option) -> String {
  self.name
}

///|
pub fn Yn00Option::value(self : Yn00Option) -> String {
  self.value
}

///|
pub fn Yn00Control::sequence_file(self : Yn00Control) -> String {
  self.sequence_file
}

///|
pub fn Yn00Control::output_file(self : Yn00Control) -> String {
  self.output_file
}

///|
pub fn Yn00Control::options(self : Yn00Control) -> Array[Yn00Option] {
  yn00_copy_options(self.options)
}

///|
pub fn Yn00Control::option(self : Yn00Control, name : String) -> String? {
  for option in self.options {
    if option.name == name {
      return Some(option.value)
    }
  }
  None
}

///|
fn yn00_pair_position(first_index : Int, second_index : Int) -> Int {
  let larger = if first_index > second_index {
    first_index
  } else {
    second_index
  }
  let smaller = if first_index > second_index {
    second_index
  } else {
    first_index
  }
  larger * (larger - 1) / 2 + smaller
}

///|
fn yn00_builders(sequence_count : Int) -> Array[Yn00PairBuilder] {
  let output : Array[Yn00PairBuilder] = []
  for second_index in 1..<sequence_count {
    for first_index in 0..<second_index {
      output.push(Yn00PairBuilder::{
        first_index,
        second_index,
        comparison_seen: false,
        ng86: None,
        yang_nielsen: None,
        lwl85: None,
        lwl85_modified: None,
        lpb93: None,
      })
    }
  }
  output
}

///|
fn yn00_assign_name(
  names : Array[String],
  index : Int,
  name : String,
) -> Unit raise Yn00Error {
  if index < 0 || index >= names.length() {
    yn00_fail("YN00 sequence index is outside the reported dimensions")
  }
  if name.length() == 0 {
    yn00_fail("YN00 sequence names must not be empty")
  }
  if names[index].length() > 0 && names[index] != name {
    yn00_fail(
      "YN00 sequence index " +
      (index + 1).to_string() +
      " has inconsistent names",
    )
  }
  for previous_index in 0..<names.length() {
    if previous_index != index &&
      names[previous_index].length() > 0 &&
      names[previous_index] == name {
      yn00_fail("Duplicate YN00 sequence name: " + name)
    }
  }
  names[index] = name
}

///|
fn yn00_comparison_header(
  line : String,
) -> (Int, String, Int, String)? raise Yn00Error {
  let value = codeml_trim(line)
  let first_open = codeml_index_of_char(value, '('.to_int(), 0)
  if first_open <= 0 {
    return None
  }
  let marker = codeml_index_of(value, ") vs.", first_open + 1)
  if marker < 0 {
    return None
  }
  let second_open = codeml_index_of_char(value, '('.to_int(), marker + 5)
  let second_close = codeml_last_index_of_char(value, ')'.to_int())
  if second_open < 0 || second_close <= second_open {
    yn00_fail("Malformed YN00 pair comparison header")
  }
  let first_number = yn00_parse_int(
    codeml_trim(value[0:first_open].to_owned()),
    "YN00 first sequence index",
  )
  let second_number = yn00_parse_int(
    codeml_trim(value[marker + 5:second_open].to_owned()),
    "YN00 second sequence index",
  )
  let first_name = value[first_open + 1:marker].to_owned()
  let second_name = value[second_open + 1:second_close].to_owned()
  Some((first_number - 1, first_name, second_number - 1, second_name))
}

///|
fn yn00_measure_token(line : String, label : String) -> String? {
  let mut from = 0
  while from < line.length() {
    let position = codeml_index_of(line, label, from)
    if position < 0 {
      return None
    }
    let before_ok = position == 0 ||
      codeml_is_space(line.unsafe_get(position - 1).to_int()) ||
      line.unsafe_get(position - 1).to_int() == '('.to_int()
    let mut equals = position + label.length()
    while equals < line.length() &&
          codeml_is_space(line.unsafe_get(equals).to_int()) {
      equals = equals + 1
    }
    if before_ok &&
      equals < line.length() &&
      line.unsafe_get(equals).to_int() == '='.to_int() {
      let mut start = equals + 1
      while start < line.length() &&
            codeml_is_space(line.unsafe_get(start).to_int()) {
        start = start + 1
      }
      let mut end = start
      while end < line.length() &&
            !codeml_is_space(line.unsafe_get(end).to_int()) &&
            line.unsafe_get(end).to_int() != ')'.to_int() {
        end = end + 1
      }
      if start == end {
        return None
      }
      return Some(line[start:end].to_owned())
    }
    from = position + label.length()
  }
  None
}

///|
fn yn00_optional_measure(
  token : String,
  label : String,
) -> Double? raise Yn00Error {
  let lower = codeml_trim(token).to_lower()
  if codeml_contains(lower, "nan") ||
    codeml_contains(lower, "inf") ||
    codeml_contains(lower, "#ind") {
    None
  } else {
    Some(yn00_parse_double(lower, label))
  }
}

///|
fn yn00_required_measure(
  line : String,
  label : String,
  method_name : String,
) -> Double? raise Yn00Error {
  match yn00_measure_token(line, label) {
    Some(token) =>
      yn00_optional_measure(token, "YN00 " + method_name + " " + label)
    None => {
      yn00_fail("YN00 " + method_name + " result is missing " + label)
      None
    }
  }
}

///|
fn yn00_counting_estimate(
  line : String,
  method_name : String,
) -> Yn00CountingEstimate raise Yn00Error {
  let ds = yn00_required_measure(line, "dS", method_name)
  let dn = yn00_required_measure(line, "dN", method_name)
  let omega = yn00_required_measure(line, "w", method_name)
  let synonymous_sites : Double? = if method_name == "LPB93" {
    None
  } else {
    yn00_required_measure(line, "S", method_name)
  }
  let nonsynonymous_sites : Double? = if method_name == "LPB93" {
    None
  } else {
    yn00_required_measure(line, "N", method_name)
  }
  let rho : Double? = if method_name == "LWL85m" {
    yn00_required_measure(line, "rho", method_name)
  } else {
    None
  }
  Yn00CountingEstimate::{
    method_name,
    ds,
    dn,
    omega,
    synonymous_sites,
    nonsynonymous_sites,
    rho,
  }
}

///|
fn yn00_parse_counting_section(
  lines : Array[String],
  start : Int,
  names : Array[String],
  builders : Array[Yn00PairBuilder],
) -> Unit raise Yn00Error {
  let mut current : Int? = None
  for line_index in start..<lines.length() {
    let trimmed = codeml_trim(lines[line_index])
    match yn00_comparison_header(trimmed) {
      Some((first_index, first_name, second_index, second_name)) => {
        if first_index <= second_index ||
          first_index >= names.length() ||
          second_index < 0 {
          yn00_fail(
            "YN00 pair header line " +
            (line_index + 1).to_string() +
            " has invalid sequence indices",
          )
        }
        yn00_assign_name(names, first_index, first_name)
        yn00_assign_name(names, second_index, second_name)
        let position = yn00_pair_position(first_index, second_index)
        let pair = builders[position]
        if pair.comparison_seen {
          yn00_fail("Duplicate YN00 pair comparison header")
        }
        pair.comparison_seen = true
        current = Some(position)
        continue
      }
      None => ()
    }
    let method_name = if codeml_starts_with(trimmed, "LWL85m:") {
      "LWL85m"
    } else if codeml_starts_with(trimmed, "LWL85:") {
      "LWL85"
    } else if codeml_starts_with(trimmed, "LPB93:") {
      "LPB93"
    } else {
      ""
    }
    if method_name.length() == 0 {
      continue
    }
    let position = match current {
      Some(value) => value
      None => {
        yn00_fail("YN00 counting estimate appears before a pair header")
        0
      }
    }
    let pair = builders[position]
    let estimate = yn00_counting_estimate(trimmed, method_name)
    match method_name {
      "LWL85" =>
        if pair.lwl85 is Some(_) {
          yn00_fail("Duplicate YN00 LWL85 estimate")
        } else {
          pair.lwl85 = Some(estimate)
        }
      "LWL85m" =>
        if pair.lwl85_modified is Some(_) {
          yn00_fail("Duplicate YN00 LWL85m estimate")
        } else {
          pair.lwl85_modified = Some(estimate)
        }
      "LPB93" =>
        if pair.lpb93 is Some(_) {
          yn00_fail("Duplicate YN00 LPB93 estimate")
        } else {
          pair.lpb93 = Some(estimate)
        }
      _ => ()
    }
  }
}

///|
fn yn00_parse_yang_nielsen_section(
  lines : Array[String],
  start : Int,
  end : Int,
  sequence_count : Int,
  builders : Array[Yn00PairBuilder],
) -> Unit raise Yn00Error {
  for line_index in start..<end {
    let trimmed = codeml_trim(lines[line_index])
    let words = codeml_words(trimmed)
    if words.length() < 2 {
      continue
    }
    let first_number = yn00_parse_int(words[0], "YN00 first sequence index") catch {
      Yn00Error(_) => continue
    }
    let second_number = yn00_parse_int(words[1], "YN00 second sequence index") catch {
      Yn00Error(_) => continue
    }
    let first_index = first_number - 1
    let second_index = second_number - 1
    if first_index <= second_index ||
      first_index >= sequence_count ||
      second_index < 0 {
      yn00_fail(
        "YN00 Yang-Nielsen row " +
        (line_index + 1).to_string() +
        " has invalid sequence indices",
      )
    }
    let values = codeml_scan_numbers(trimmed)
    if values.length() != 11 {
      yn00_fail(
        "YN00 Yang-Nielsen row " +
        (line_index + 1).to_string() +
        " must contain nine estimates",
      )
    }
    let position = yn00_pair_position(first_index, second_index)
    let pair = builders[position]
    if pair.yang_nielsen is Some(_) {
      yn00_fail("Duplicate YN00 Yang-Nielsen pair estimate")
    }
    pair.yang_nielsen = Some(Yn00YangNielsenEstimate::{
      synonymous_sites: values[2],
      nonsynonymous_sites: values[3],
      time: values[4],
      kappa: values[5],
      omega: values[6],
      dn: values[7],
      dn_standard_error: values[8],
      ds: values[9],
      ds_standard_error: values[10],
    })
  }
}

///|
fn yn00_parse_ng_remainder(
  remainder : String,
  expected : Int,
) -> Array[Yn00Ng86Estimate] raise Yn00Error {
  let estimates : Array[Yn00Ng86Estimate] = []
  let mut cursor = 0
  while cursor < remainder.length() {
    while cursor < remainder.length() &&
          codeml_is_space(remainder.unsafe_get(cursor).to_int()) {
      cursor = cursor + 1
    }
    if cursor == remainder.length() {
      break
    }
    let open = codeml_index_of_char(remainder, '('.to_int(), cursor)
    if open < 0 {
      yn00_fail("Malformed YN00 NG86 matrix row")
    }
    let close = codeml_index_of_char(remainder, ')'.to_int(), open + 1)
    if close < 0 {
      yn00_fail("Malformed YN00 NG86 matrix row")
    }
    let omega_text = codeml_trim(remainder[cursor:open].to_owned())
    let inside = remainder[open + 1:close].to_owned()
    let values = codeml_scan_numbers(inside)
    if omega_text.length() == 0 || values.length() != 2 {
      yn00_fail("Malformed YN00 NG86 matrix estimate")
    }
    estimates.push(Yn00Ng86Estimate::{
      omega: yn00_parse_double(omega_text, "YN00 NG86 omega"),
      dn: values[0],
      ds: values[1],
    })
    cursor = close + 1
  }
  if estimates.length() != expected {
    yn00_fail("YN00 NG86 matrix row has the wrong number of estimates")
  }
  estimates
}

///|
fn yn00_parse_ng86_section(
  lines : Array[String],
  start : Int,
  end : Int,
  names : Array[String],
  builders : Array[Yn00PairBuilder],
) -> Unit raise Yn00Error {
  for sequence_index in 0..<names.length() {
    let name = names[sequence_index]
    let mut found = false
    for line_index in start..<end {
      let trimmed = codeml_trim(lines[line_index])
      if !codeml_starts_with(trimmed, name) {
        continue
      }
      let remainder = trimmed[name.length():trimmed.length()].to_owned()
      if sequence_index == 0 {
        if codeml_trim(remainder).length() == 0 {
          found = true
          break
        }
        continue
      }
      let estimates = yn00_parse_ng_remainder(remainder, sequence_index) catch {
        Yn00Error(_) => continue
      }
      for previous_index in 0..<sequence_index {
        let position = yn00_pair_position(sequence_index, previous_index)
        let pair = builders[position]
        if pair.ng86 is Some(_) {
          yn00_fail("Duplicate YN00 NG86 pair estimate")
        }
        pair.ng86 = Some(estimates[previous_index])
      }
      found = true
      break
    }
    if !found {
      yn00_fail("YN00 NG86 matrix is missing sequence row: " + name)
    }
  }
}

///|
fn yn00_find_sections(lines : Array[String]) -> (Int, Int, Int) raise Yn00Error {
  let mut section_a = -1
  let mut section_b = -1
  let mut section_c = -1
  for index in 0..<lines.length() {
    let line = lines[index]
    if codeml_contains(line, "(A) Nei-Gojobori (1986) method") {
      if section_a >= 0 {
        yn00_fail("Duplicate YN00 NG86 section")
      }
      section_a = index
    } else if codeml_contains(line, "(B) Yang & Nielsen (2000) method") {
      if section_b >= 0 {
        yn00_fail("Duplicate YN00 Yang-Nielsen section")
      }
      section_b = index
    } else if codeml_contains(line, "(C) LWL85, LPB93 & LWLm methods") {
      if section_c >= 0 {
        yn00_fail("Duplicate YN00 counting-method section")
      }
      section_c = index
    }
  }
  if section_a < 0 || section_b <= section_a || section_c <= section_b {
    yn00_fail("YN00 results do not contain all three method sections")
  }
  (section_a, section_b, section_c)
}

///|
fn yn00_metadata(
  lines : Array[String],
  section_a : Int,
) -> (String, Int, Int) raise Yn00Error {
  let mut alignment_file = ""
  let mut sequence_count = -1
  let mut codon_count = -1
  for index in 0..<section_a {
    let trimmed = codeml_trim(lines[index])
    if codeml_starts_with(trimmed, "YN00 ") {
      if alignment_file.length() > 0 {
        yn00_fail("Duplicate YN00 alignment header")
      }
      alignment_file = codeml_trim(trimmed[5:trimmed.length()].to_owned())
    }
    if codeml_starts_with(trimmed, "ns =") && codeml_contains(trimmed, "ls =") {
      let values = codeml_scan_numbers(trimmed)
      if values.length() != 2 {
        yn00_fail("Malformed YN00 alignment dimensions")
      }
      sequence_count = values[0].to_int()
      codon_count = values[1].to_int()
    }
  }
  if alignment_file.length() == 0 {
    yn00_fail("YN00 results do not contain an alignment header")
  }
  if sequence_count < 2 || codon_count <= 0 {
    yn00_fail("YN00 results contain invalid alignment dimensions")
  }
  (alignment_file, sequence_count, codon_count)
}

///|
fn yn00_finish_pairs(
  names : Array[String],
  builders : Array[Yn00PairBuilder],
) -> Array[Yn00PairResult] raise Yn00Error {
  for name in names {
    if name.length() == 0 {
      yn00_fail("YN00 results do not identify every sequence")
    }
  }
  let pairs : Array[Yn00PairResult] = []
  for builder in builders {
    if !builder.comparison_seen {
      yn00_fail("YN00 results are missing a pair comparison")
    }
    let ng86 = match builder.ng86 {
      Some(value) => value
      None => {
        yn00_fail("YN00 pair is missing its NG86 estimate")
        Yn00Ng86Estimate::{ omega: 0.0, dn: 0.0, ds: 0.0, }
      }
    }
    let yang_nielsen = match builder.yang_nielsen {
      Some(value) => value
      None => {
        yn00_fail("YN00 pair is missing its Yang-Nielsen estimate")
        Yn00YangNielsenEstimate::{
          synonymous_sites: 0.0,
          nonsynonymous_sites: 0.0,
          time: 0.0,
          kappa: 0.0,
          omega: 0.0,
          dn: 0.0,
          dn_standard_error: 0.0,
          ds: 0.0,
          ds_standard_error: 0.0,
        }
      }
    }
    let lwl85 = match builder.lwl85 {
      Some(value) => value
      None => {
        yn00_fail("YN00 pair is missing its LWL85 estimate")
        Yn00CountingEstimate::{
          method_name: "LWL85",
          ds: None,
          dn: None,
          omega: None,
          synonymous_sites: None,
          nonsynonymous_sites: None,
          rho: None,
        }
      }
    }
    let lwl85_modified = match builder.lwl85_modified {
      Some(value) => value
      None => {
        yn00_fail("YN00 pair is missing its LWL85m estimate")
        Yn00CountingEstimate::{
          method_name: "LWL85m",
          ds: None,
          dn: None,
          omega: None,
          synonymous_sites: None,
          nonsynonymous_sites: None,
          rho: None,
        }
      }
    }
    let lpb93 = match builder.lpb93 {
      Some(value) => value
      None => {
        yn00_fail("YN00 pair is missing its LPB93 estimate")
        Yn00CountingEstimate::{
          method_name: "LPB93",
          ds: None,
          dn: None,
          omega: None,
          synonymous_sites: None,
          nonsynonymous_sites: None,
          rho: None,
        }
      }
    }
    pairs.push(Yn00PairResult::{
      first: names[builder.first_index],
      second: names[builder.second_index],
      ng86,
      yang_nielsen,
      lwl85,
      lwl85_modified,
      lpb93,
    })
  }
  pairs
}

///|
/// Parse a complete YN00 result artifact.
pub fn yn00_parse_results(text : String) -> Yn00Results raise Yn00Error {
  if codeml_trim(text).length() == 0 {
    yn00_fail("YN00 result text must not be empty")
  }
  let lines = codeml_lines(text)
  let (section_a, section_b, section_c) = yn00_find_sections(lines)
  let (alignment_file, sequence_count, codon_count) = yn00_metadata(
    lines, section_a,
  )
  let names = Array::make(sequence_count, "")
  let builders = yn00_builders(sequence_count)
  yn00_parse_counting_section(lines, section_c + 1, names, builders)
  yn00_parse_yang_nielsen_section(
    lines,
    section_b + 1,
    section_c,
    sequence_count,
    builders,
  )
  yn00_parse_ng86_section(lines, section_a + 1, section_b, names, builders)
  let pairs = yn00_finish_pairs(names, builders)
  let expected_pairs = sequence_count * (sequence_count - 1) / 2
  if pairs.length() != expected_pairs {
    yn00_fail("YN00 pair count does not match the reported sequence count")
  }
  Yn00Results::{
    alignment_file,
    sequence_count,
    codon_count,
    sequence_names: names,
    pairs,
  }
}

///|
pub fn Yn00Ng86Estimate::omega(self : Yn00Ng86Estimate) -> Double {
  self.omega
}

///|
pub fn Yn00Ng86Estimate::dn(self : Yn00Ng86Estimate) -> Double {
  self.dn
}

///|
pub fn Yn00Ng86Estimate::ds(self : Yn00Ng86Estimate) -> Double {
  self.ds
}

///|
pub fn Yn00YangNielsenEstimate::synonymous_sites(
  self : Yn00YangNielsenEstimate,
) -> Double {
  self.synonymous_sites
}

///|
pub fn Yn00YangNielsenEstimate::nonsynonymous_sites(
  self : Yn00YangNielsenEstimate,
) -> Double {
  self.nonsynonymous_sites
}

///|
pub fn Yn00YangNielsenEstimate::time(self : Yn00YangNielsenEstimate) -> Double {
  self.time
}

///|
pub fn Yn00YangNielsenEstimate::kappa(self : Yn00YangNielsenEstimate) -> Double {
  self.kappa
}

///|
pub fn Yn00YangNielsenEstimate::omega(self : Yn00YangNielsenEstimate) -> Double {
  self.omega
}

///|
pub fn Yn00YangNielsenEstimate::dn(self : Yn00YangNielsenEstimate) -> Double {
  self.dn
}

///|
pub fn Yn00YangNielsenEstimate::dn_standard_error(
  self : Yn00YangNielsenEstimate,
) -> Double {
  self.dn_standard_error
}

///|
pub fn Yn00YangNielsenEstimate::ds(self : Yn00YangNielsenEstimate) -> Double {
  self.ds
}

///|
pub fn Yn00YangNielsenEstimate::ds_standard_error(
  self : Yn00YangNielsenEstimate,
) -> Double {
  self.ds_standard_error
}

///|
pub fn Yn00CountingEstimate::method_name(self : Yn00CountingEstimate) -> String {
  self.method_name
}

///|
pub fn Yn00CountingEstimate::ds(self : Yn00CountingEstimate) -> Double? {
  self.ds
}

///|
pub fn Yn00CountingEstimate::dn(self : Yn00CountingEstimate) -> Double? {
  self.dn
}

///|
pub fn Yn00CountingEstimate::omega(self : Yn00CountingEstimate) -> Double? {
  self.omega
}

///|
pub fn Yn00CountingEstimate::synonymous_sites(
  self : Yn00CountingEstimate,
) -> Double? {
  self.synonymous_sites
}

///|
pub fn Yn00CountingEstimate::nonsynonymous_sites(
  self : Yn00CountingEstimate,
) -> Double? {
  self.nonsynonymous_sites
}

///|
pub fn Yn00CountingEstimate::rho(self : Yn00CountingEstimate) -> Double? {
  self.rho
}

///|
pub fn Yn00PairResult::first(self : Yn00PairResult) -> String {
  self.first
}

///|
pub fn Yn00PairResult::second(self : Yn00PairResult) -> String {
  self.second
}

///|
pub fn Yn00PairResult::ng86(self : Yn00PairResult) -> Yn00Ng86Estimate {
  self.ng86
}

///|
pub fn Yn00PairResult::yang_nielsen(
  self : Yn00PairResult,
) -> Yn00YangNielsenEstimate {
  self.yang_nielsen
}

///|
pub fn Yn00PairResult::lwl85(self : Yn00PairResult) -> Yn00CountingEstimate {
  self.lwl85
}

///|
pub fn Yn00PairResult::lwl85_modified(
  self : Yn00PairResult,
) -> Yn00CountingEstimate {
  self.lwl85_modified
}

///|
pub fn Yn00PairResult::lpb93(self : Yn00PairResult) -> Yn00CountingEstimate {
  self.lpb93
}

///|
pub fn Yn00Results::alignment_file(self : Yn00Results) -> String {
  self.alignment_file
}

///|
pub fn Yn00Results::sequence_count(self : Yn00Results) -> Int {
  self.sequence_count
}

///|
pub fn Yn00Results::codon_count(self : Yn00Results) -> Int {
  self.codon_count
}

///|
pub fn Yn00Results::sequence_names(self : Yn00Results) -> Array[String] {
  yn00_copy_strings(self.sequence_names)
}

///|
pub fn Yn00Results::pairs(self : Yn00Results) -> Array[Yn00PairResult] {
  yn00_copy_pairs(self.pairs)
}

///|
/// Look up one pair without regard to name order.
pub fn Yn00Results::pair(
  self : Yn00Results,
  first : String,
  second : String,
) -> Yn00PairResult? {
  if first == second {
    return None
  }
  for pair in self.pairs {
    if (pair.first == first && pair.second == second) ||
      (pair.first == second && pair.second == first) {
      return Some(pair)
    }
  }
  None
}

///|
fn yn00_pair_statistic(
  pair : Yn00PairResult,
  method_name : String,
  statistic : String,
) -> Double? raise Yn00Error {
  match method_name {
    "NG86" =>
      match statistic {
        "omega" | "w" => Some(pair.ng86.omega)
        "dN" => Some(pair.ng86.dn)
        "dS" => Some(pair.ng86.ds)
        _ => {
          yn00_fail("Unknown NG86 statistic: " + statistic)
          None
        }
      }
    "YN00" =>
      match statistic {
        "S" => Some(pair.yang_nielsen.synonymous_sites)
        "N" => Some(pair.yang_nielsen.nonsynonymous_sites)
        "t" => Some(pair.yang_nielsen.time)
        "kappa" => Some(pair.yang_nielsen.kappa)
        "omega" | "w" => Some(pair.yang_nielsen.omega)
        "dN" => Some(pair.yang_nielsen.dn)
        "dN SE" => Some(pair.yang_nielsen.dn_standard_error)
        "dS" => Some(pair.yang_nielsen.ds)
        "dS SE" => Some(pair.yang_nielsen.ds_standard_error)
        _ => {
          yn00_fail("Unknown YN00 statistic: " + statistic)
          None
        }
      }
    "LWL85" | "LWL85m" | "LPB93" => {
      let estimate = match method_name {
        "LWL85" => pair.lwl85
        "LWL85m" => pair.lwl85_modified
        _ => pair.lpb93
      }
      match statistic {
        "omega" | "w" => estimate.omega
        "dN" => estimate.dn
        "dS" => estimate.ds
        "S" => estimate.synonymous_sites
        "N" => estimate.nonsynonymous_sites
        "rho" => estimate.rho
        _ => {
          yn00_fail("Unknown " + method_name + " statistic: " + statistic)
          None
        }
      }
    }
    _ => {
      yn00_fail("Unknown YN00 method: " + method_name)
      None
    }
  }
}

///|
/// Return a symmetric matrix for a method/statistic combination.
///
/// Diagonal cells and undefined PAML estimates are represented by `None`.
pub fn Yn00Results::matrix(
  self : Yn00Results,
  method_name : String,
  statistic : String,
) -> Yn00DistanceMatrix raise Yn00Error {
  let values : Array[Array[Double?]] = []
  for _ in 0..<self.sequence_count {
    values.push(Array::make(self.sequence_count, None))
  }
  for pair in self.pairs {
    let value = yn00_pair_statistic(pair, method_name, statistic)
    let mut first_index = -1
    let mut second_index = -1
    for index in 0..<self.sequence_names.length() {
      if self.sequence_names[index] == pair.first {
        first_index = index
      }
      if self.sequence_names[index] == pair.second {
        second_index = index
      }
    }
    values[first_index][second_index] = value
    values[second_index][first_index] = value
  }
  Yn00DistanceMatrix::{ names: yn00_copy_strings(self.sequence_names), values, }
}

///|
/// Mean of all defined off-diagonal pair estimates.
pub fn Yn00Results::mean(
  self : Yn00Results,
  method_name : String,
  statistic : String,
) -> Double raise Yn00Error {
  let mut total = 0.0
  let mut count = 0
  for pair in self.pairs {
    match yn00_pair_statistic(pair, method_name, statistic) {
      Some(value) => {
        total = total + value
        count = count + 1
      }
      None => ()
    }
  }
  if count == 0 {
    yn00_fail(
      "YN00 " + method_name + " " + statistic + " has no defined pair estimates",
    )
  }
  total / count.to_double()
}

///|
pub fn Yn00DistanceMatrix::names(self : Yn00DistanceMatrix) -> Array[String] {
  yn00_copy_strings(self.names)
}

///|
pub fn Yn00DistanceMatrix::values(
  self : Yn00DistanceMatrix,
) -> Array[Array[Double?]] {
  yn00_copy_optional_matrix(self.values)
}

///|
pub fn Yn00DistanceMatrix::value(
  self : Yn00DistanceMatrix,
  first : String,
  second : String,
) -> Double? {
  let mut first_index = -1
  let mut second_index = -1
  for index in 0..<self.names.length() {
    if self.names[index] == first {
      first_index = index
    }
    if self.names[index] == second {
      second_index = index
    }
  }
  if first_index < 0 || second_index < 0 {
    None
  } else {
    self.values[first_index][second_index]
  }
}

///|
/// Compact control fixture adapted from Biopython's official YN00 data.
pub fn yn00_example_control_text() -> String {
  "seqfile = alignment.phylip * coding alignment\n" +
  "outfile = yn00.out\n" +
  "verbose = 1\n" +
  "icode = 0\n" +
  "weighting = 0\n" +
  "commonf3x4 = 0\n" +
  "ndata = 1\n"
}

///|
/// Compact result fixture preserving all five YN00 methods.
pub fn yn00_example_results_text() -> String {
  "YN00 alignment.phylip\n" +
  "ns = 3 ls = 74\n" +
  "\n" +
  "(A) Nei-Gojobori (1986) method\n" +
  "Nei & Gojobori 1986. dN/dS (dN, dS)\n" +
  "Homo_sapie\n" +
  "Pan_troglo -1.0000 (0.0000 0.0207)\n" +
  "Gorilla_go 0.5000 (0.0100 0.0200)-1.0000 (0.0000 0.0421)\n" +
  "\n" +
  "(B) Yang & Nielsen (2000) method\n" +
  "seq. seq. S N t kappa omega dN +- SE dS +- SE\n" +
  "2 1 67.3 154.7 0.0136 3.6564 0.0000 0.0000 +- 0.0000 0.0150 +- 0.0151\n" +
  "3 1 67.4 154.6 0.0200 3.7000 0.5000 0.0100 +- 0.0050 0.0200 +- 0.0100\n" +
  "3 2 67.3 154.7 0.0276 3.6564 0.0000 0.0000 +- 0.0000 0.0303 +- 0.0215\n" +
  "\n" +
  "(C) LWL85, LPB93 & LWLm methods\n" +
  "2 (Pan_troglo) vs. 1 (Homo_sapie)\n" +
  "LWL85: dS = 0.0227 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0\n" +
  "LWL85m: dS = -nan dN = -nan w = -nan S = -nan N = -nan (rho = -nan)\n" +
  "LPB93: dS = 0.0129 dN = 0.0000 w = 0.0000\n" +
  "3 (Gorilla_go) vs. 1 (Homo_sapie)\n" +
  "LWL85: dS = 0.0200 dN = 0.0100 w = 0.5000 S = 45.0 N = 177.0\n" +
  "LWL85m: dS = 0.0180 dN = 0.0110 w = 0.6111 S = 53.9 N = 168.1 (rho = 0.507)\n" +
  "LPB93: dS = 0.0190 dN = 0.0100 w = 0.5263\n" +
  "3 (Gorilla_go) vs. 2 (Pan_troglo)\n" +
  "LWL85: dS = 0.0463 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0\n" +
  "LWL85m: dS = -nan dN = -nan w = -nan S = -nan N = -nan (rho = -nan)\n" +
  "LPB93: dS = 0.0264 dN = 0.0000 w = 0.0000\n"
}
