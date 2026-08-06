// Biopython-compatible Bio.Phylo.PAML.codeml support.
//
// This module parses CODEML control files and CODONML/AAML result text. It
// deliberately does not launch the external PAML executable; callers can
// generate a control file, run CODEML separately, and parse its output here.

///|
/// Error raised for malformed CODEML configuration, output, or comparisons.
pub suberror CodemlError {
  CodemlError(String)
}

///|
/// One validated CODEML option retained in canonical text form.
pub struct CodemlOption {
  name : String
  value : String
} derive(Eq, Debug)

///|
/// A CODEML control file.
pub struct CodemlControl {
  sequence_file : String
  output_file : String
  tree_file : String
  options : Array[CodemlOption]
} derive(Eq, Debug)

///|
/// One named floating-point parameter.
pub struct CodemlNamedValue {
  name : String
  value : Double
} derive(Eq, Debug)

///|
/// Per-gene kappa and omega estimates from a multi-gene model.
pub struct CodemlGeneParameter {
  gene : Int
  kappa : Double
  omega : Double
} derive(Eq, Debug)

///|
/// One CODEML site class.
pub struct CodemlSiteClass {
  index : Int
  proportion : Double
  omega : Double?
  branch_types : Array[CodemlNamedValue]
} derive(Eq, Debug)

///|
/// One row from CODEML's branch dN/dS table.
///
/// Fields are optional because PAML may print `nan` for an unidentifiable
/// branch quantity.
pub struct CodemlBranchResult {
  branch : String
  time : Double?
  nonsynonymous_sites : Double?
  synonymous_sites : Double?
  omega : Double?
  dn : Double?
  ds : Double?
  n_times_dn : Double?
  s_times_ds : Double?
} derive(Eq, Debug)

///|
/// One positively selected site reported by NEB or BEB.
pub struct CodemlPositiveSite {
  analysis_method : String
  position : Int
  amino_acid : String
  probability : Double
  significance : String
  posterior_mean_omega : Double?
  standard_error : Double?
} derive(Eq, Debug)

///|
/// One fitted CODEML NSsites model or one multi-gene result section.
pub struct CodemlModelResult {
  number : Int
  gene : Int
  description : String
  ln_likelihood : Double
  parameter_count : Int
  parameter_list : String
  standard_errors : String
  tree_length : Double?
  tree : String
  dn_tree : String
  ds_tree : String
  omega_tree : String
  kappa : Double?
  omegas : Array[Double]
  dn_tree_length : Double?
  ds_tree_length : Double?
  rates : Array[Double]
  gene_parameters : Array[CodemlGeneParameter]
  named_parameters : Array[CodemlNamedValue]
  site_classes : Array[CodemlSiteClass]
  branches : Array[CodemlBranchResult]
  positive_sites : Array[CodemlPositiveSite]
} derive(Eq, Debug)

///|
/// One maximum-likelihood pairwise codon comparison.
pub struct CodemlPairwiseResult {
  first : String
  second : String
  ln_likelihood : Double
  time : Double
  synonymous_sites : Double
  nonsynonymous_sites : Double
  omega : Double
  dn : Double
  ds : Double
} derive(Eq, Debug)

///|
/// One symmetric amino-acid distance matrix entry.
pub struct CodemlDistance {
  kind : String
  first : String
  second : String
  value : Double
} derive(Eq, Debug)

///|
/// Parsed CODEML or AAML results.
pub struct CodemlResults {
  program : String
  version : String
  model : String
  codon_frequency_model : String
  site_class_model : String
  sequence_count : Int
  site_count : Int
  ln_l_max : Double?
  models : Array[CodemlModelResult]
  pairwise : Array[CodemlPairwiseResult]
  distances : Array[CodemlDistance]
} derive(Eq, Debug)

///|
/// Likelihood-ratio comparison between nested CODEML models.
pub struct CodemlLikelihoodRatio {
  statistic : Double
  degrees_of_freedom : Int
  p_value : Double
  alpha : Double
  significant : Bool
} derive(Eq, Debug)

///|
priv struct CodemlPendingPair {
  first : String
  second : String
  mut ln_likelihood : Double?
}

///|
fn codeml_fail(message : String) -> Unit raise CodemlError {
  raise CodemlError(message)
}

///|
fn codeml_is_space(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn codeml_trim(value : String) -> String {
  let mut start = 0
  let mut end = value.length()
  while start < end && codeml_is_space(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  while end > start && codeml_is_space(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn codeml_starts_with(value : String, prefix : String) -> Bool {
  if prefix.length() > value.length() {
    return false
  }
  for index in 0..<prefix.length() {
    if value.unsafe_get(index).to_int() != prefix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn codeml_ends_with(value : String, suffix : String) -> Bool {
  if suffix.length() > value.length() {
    return false
  }
  let offset = value.length() - suffix.length()
  for index in 0..<suffix.length() {
    if value.unsafe_get(offset + index).to_int() !=
      suffix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn codeml_index_of_char(value : String, wanted : Int, start : Int) -> Int {
  let mut index = start
  while index < value.length() {
    if value.unsafe_get(index).to_int() == wanted {
      return index
    }
    index = index + 1
  }
  -1
}

///|
fn codeml_last_index_of_char(value : String, wanted : Int) -> Int {
  let mut index = value.length() - 1
  while index >= 0 {
    if value.unsafe_get(index).to_int() == wanted {
      return index
    }
    index = index - 1
  }
  -1
}

///|
fn codeml_contains(value : String, needle : String) -> Bool {
  if needle.length() == 0 {
    return true
  }
  if needle.length() > value.length() {
    return false
  }
  for start in 0..<=(value.length() - needle.length()) {
    let mut matched = true
    for offset in 0..<needle.length() {
      if value.unsafe_get(start + offset).to_int() !=
        needle.unsafe_get(offset).to_int() {
        matched = false
        break
      }
    }
    if matched {
      return true
    }
  }
  false
}

///|
fn codeml_index_of(value : String, needle : String, from : Int) -> Int {
  if needle.length() == 0 {
    return from
  }
  if from < 0 || from + needle.length() > value.length() {
    return -1
  }
  for start in from..<=(value.length() - needle.length()) {
    let mut matched = true
    for offset in 0..<needle.length() {
      if value.unsafe_get(start + offset).to_int() !=
        needle.unsafe_get(offset).to_int() {
        matched = false
        break
      }
    }
    if matched {
      return start
    }
  }
  -1
}

///|
fn codeml_lines(text : String) -> Array[String] {
  let output : Array[String] = []
  for view in text.split("\n") {
    let owned = view.to_owned()
    if owned.length() > 0 &&
      owned.unsafe_get(owned.length() - 1).to_int() == '\r'.to_int() {
      output.push(owned[0:owned.length() - 1].to_owned())
    } else {
      output.push(owned)
    }
  }
  output
}

///|
fn codeml_words(value : String) -> Array[String] {
  let words : Array[String] = []
  let mut index = 0
  while index < value.length() {
    while index < value.length() &&
          codeml_is_space(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    let start = index
    while index < value.length() &&
          !codeml_is_space(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if start < index {
      words.push(value[start:index].to_owned())
    }
  }
  words
}

///|
fn codeml_is_digit(code : Int) -> Bool {
  code >= '0'.to_int() && code <= '9'.to_int()
}

///|
fn codeml_parse_int(text : String, label : String) -> Int raise CodemlError {
  let value = codeml_trim(text)
  if value.length() == 0 {
    codeml_fail(label + " must be an integer")
  }
  let mut index = 0
  let mut sign = 1
  if value.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1
    index = 1
  } else if value.unsafe_get(0).to_int() == '+'.to_int() {
    index = 1
  }
  if index == value.length() {
    codeml_fail(label + " must be an integer")
  }
  let mut result = 0
  while index < value.length() {
    let code = value.unsafe_get(index).to_int()
    if !codeml_is_digit(code) {
      codeml_fail(label + " must be an integer")
    }
    let digit = code - '0'.to_int()
    if result > 214748364 || (result == 214748364 && digit > 7) {
      codeml_fail(label + " exceeds the supported integer range")
    }
    result = result * 10 + digit
    index = index + 1
  }
  sign * result
}

///|
fn codeml_is_number(text : String) -> Bool {
  let value = codeml_trim(text)
  if value.length() == 0 {
    return false
  }
  let mut index = 0
  if value.unsafe_get(index).to_int() == '-'.to_int() ||
    value.unsafe_get(index).to_int() == '+'.to_int() {
    index = index + 1
  }
  let mut digits = 0
  while index < value.length() &&
        codeml_is_digit(value.unsafe_get(index).to_int()) {
    digits = digits + 1
    index = index + 1
  }
  if index < value.length() && value.unsafe_get(index).to_int() == '.'.to_int() {
    index = index + 1
    while index < value.length() &&
          codeml_is_digit(value.unsafe_get(index).to_int()) {
      digits = digits + 1
      index = index + 1
    }
  }
  if digits == 0 {
    return false
  }
  if index < value.length() &&
    (
      value.unsafe_get(index).to_int() == 'e'.to_int() ||
      value.unsafe_get(index).to_int() == 'E'.to_int()
    ) {
    index = index + 1
    if index < value.length() &&
      (
        value.unsafe_get(index).to_int() == '-'.to_int() ||
        value.unsafe_get(index).to_int() == '+'.to_int()
      ) {
      index = index + 1
    }
    let exponent_start = index
    while index < value.length() &&
          codeml_is_digit(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if exponent_start == index {
      return false
    }
  }
  index == value.length()
}

///|
fn codeml_parse_double(
  text : String,
  label : String,
) -> Double raise CodemlError {
  let value = codeml_trim(text)
  if !codeml_is_number(value) {
    codeml_fail(label + " must be a finite number")
  }
  let normalized = if value.unsafe_get(0).to_int() == '+'.to_int() {
    value[1:value.length()].to_owned()
  } else {
    value
  }
  let result = match parse_double(normalized) {
    Some(number) => number
    None => {
      codeml_fail(label + " must be a finite number")
      0.0
    }
  }
  if result.is_nan() || result.abs() > 1.0e300 {
    codeml_fail(label + " must be a finite number")
  }
  result
}

///|
fn codeml_parse_optional_double(text : String) -> Double? raise CodemlError {
  let value = codeml_trim(text)
  if value.to_lower() == "nan" || value == "-nan" {
    None
  } else {
    Some(codeml_parse_double(value, "CODEML branch value"))
  }
}

///|
fn codeml_scan_numbers(text : String) -> Array[Double] {
  let values : Array[Double] = []
  let mut index = 0
  while index < text.length() {
    let code = text.unsafe_get(index).to_int()
    let starts = codeml_is_digit(code) ||
      code == '.'.to_int() ||
      (
        (code == '-'.to_int() || code == '+'.to_int()) &&
        index + 1 < text.length() &&
        (
          codeml_is_digit(text.unsafe_get(index + 1).to_int()) ||
          text.unsafe_get(index + 1).to_int() == '.'.to_int()
        )
      )
    if !starts {
      index = index + 1
      continue
    }
    let start = index
    if code == '-'.to_int() || code == '+'.to_int() {
      index = index + 1
    }
    let mut digits = 0
    while index < text.length() &&
          codeml_is_digit(text.unsafe_get(index).to_int()) {
      digits = digits + 1
      index = index + 1
    }
    if index < text.length() && text.unsafe_get(index).to_int() == '.'.to_int() {
      index = index + 1
      while index < text.length() &&
            codeml_is_digit(text.unsafe_get(index).to_int()) {
        digits = digits + 1
        index = index + 1
      }
    }
    if digits == 0 {
      index = start + 1
      continue
    }
    if index < text.length() &&
      (
        text.unsafe_get(index).to_int() == 'e'.to_int() ||
        text.unsafe_get(index).to_int() == 'E'.to_int()
      ) {
      let exponent_mark = index
      index = index + 1
      if index < text.length() &&
        (
          text.unsafe_get(index).to_int() == '-'.to_int() ||
          text.unsafe_get(index).to_int() == '+'.to_int()
        ) {
        index = index + 1
      }
      let exponent_start = index
      while index < text.length() &&
            codeml_is_digit(text.unsafe_get(index).to_int()) {
        index = index + 1
      }
      if exponent_start == index {
        index = exponent_mark
      }
    }
    let token = text[start:index].to_owned()
    if codeml_is_number(token) {
      match
        parse_double(
          if token.unsafe_get(0).to_int() == '+'.to_int() {
            token[1:token.length()].to_owned()
          } else {
            token
          },
        ) {
        Some(value) =>
          if !value.is_nan() && value.abs() <= 1.0e300 {
            values.push(value)
          }
        None => ()
      }
    }
  }
  values
}

///|
fn codeml_scan_fixed_five_numbers(text : String) -> Array[Double] {
  let values : Array[Double] = []
  let mut index = 0
  while index < text.length() {
    if !codeml_is_digit(text.unsafe_get(index).to_int()) {
      index = index + 1
      continue
    }
    let start = index
    let mut whole_digits = 0
    while index < text.length() &&
          codeml_is_digit(text.unsafe_get(index).to_int()) &&
          whole_digits < 3 {
      whole_digits = whole_digits + 1
      index = index + 1
    }
    if index >= text.length() || text.unsafe_get(index).to_int() != '.'.to_int() {
      index = start + 1
      continue
    }
    index = index + 1
    let fraction_start = index
    while index < text.length() &&
          codeml_is_digit(text.unsafe_get(index).to_int()) &&
          index - fraction_start < 5 {
      index = index + 1
    }
    if index - fraction_start != 5 {
      index = start + 1
      continue
    }
    match parse_double(text[start:index].to_owned()) {
      Some(value) => values.push(value)
      None => ()
    }
  }
  values
}

///|
fn codeml_copy_doubles(values : Array[Double]) -> Array[Double] {
  let copied : Array[Double] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_copy_options(values : Array[CodemlOption]) -> Array[CodemlOption] {
  let copied : Array[CodemlOption] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_copy_named_values(
  values : Array[CodemlNamedValue],
) -> Array[CodemlNamedValue] {
  let copied : Array[CodemlNamedValue] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_clone_site_class(value : CodemlSiteClass) -> CodemlSiteClass {
  CodemlSiteClass::{
    index: value.index,
    proportion: value.proportion,
    omega: value.omega,
    branch_types: codeml_copy_named_values(value.branch_types),
  }
}

///|
fn codeml_copy_site_classes(
  values : Array[CodemlSiteClass],
) -> Array[CodemlSiteClass] {
  let copied : Array[CodemlSiteClass] = []
  for value in values {
    copied.push(codeml_clone_site_class(value))
  }
  copied
}

///|
fn codeml_copy_gene_parameters(
  values : Array[CodemlGeneParameter],
) -> Array[CodemlGeneParameter] {
  let copied : Array[CodemlGeneParameter] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_copy_branches(
  values : Array[CodemlBranchResult],
) -> Array[CodemlBranchResult] {
  let copied : Array[CodemlBranchResult] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_copy_positive_sites(
  values : Array[CodemlPositiveSite],
) -> Array[CodemlPositiveSite] {
  let copied : Array[CodemlPositiveSite] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_clone_model(value : CodemlModelResult) -> CodemlModelResult {
  CodemlModelResult::{
    number: value.number,
    gene: value.gene,
    description: value.description,
    ln_likelihood: value.ln_likelihood,
    parameter_count: value.parameter_count,
    parameter_list: value.parameter_list,
    standard_errors: value.standard_errors,
    tree_length: value.tree_length,
    tree: value.tree,
    dn_tree: value.dn_tree,
    ds_tree: value.ds_tree,
    omega_tree: value.omega_tree,
    kappa: value.kappa,
    omegas: codeml_copy_doubles(value.omegas),
    dn_tree_length: value.dn_tree_length,
    ds_tree_length: value.ds_tree_length,
    rates: codeml_copy_doubles(value.rates),
    gene_parameters: codeml_copy_gene_parameters(value.gene_parameters),
    named_parameters: codeml_copy_named_values(value.named_parameters),
    site_classes: codeml_copy_site_classes(value.site_classes),
    branches: codeml_copy_branches(value.branches),
    positive_sites: codeml_copy_positive_sites(value.positive_sites),
  }
}

///|
fn codeml_copy_models(
  values : Array[CodemlModelResult],
) -> Array[CodemlModelResult] {
  let copied : Array[CodemlModelResult] = []
  for value in values {
    copied.push(codeml_clone_model(value))
  }
  copied
}

///|
fn codeml_copy_pairwise(
  values : Array[CodemlPairwiseResult],
) -> Array[CodemlPairwiseResult] {
  let copied : Array[CodemlPairwiseResult] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_copy_distances(
  values : Array[CodemlDistance],
) -> Array[CodemlDistance] {
  let copied : Array[CodemlDistance] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn codeml_is_allowed_option(name : String) -> Bool {
  match name {
    "noisy"
    | "verbose"
    | "runmode"
    | "seqtype"
    | "CodonFreq"
    | "ndata"
    | "clock"
    | "aaDist"
    | "aaRatefile"
    | "model"
    | "NSsites"
    | "icode"
    | "Mgene"
    | "fix_kappa"
    | "kappa"
    | "fix_omega"
    | "omega"
    | "fix_alpha"
    | "alpha"
    | "Malpha"
    | "ncatG"
    | "getSE"
    | "RateAncestor"
    | "Small_Diff"
    | "cleandata"
    | "fix_blength"
    | "method"
    | "rho"
    | "fix_rho" => true
    _ => false
  }
}

///|
fn codeml_is_integer_option(name : String) -> Bool {
  match name {
    "noisy"
    | "verbose"
    | "runmode"
    | "seqtype"
    | "CodonFreq"
    | "ndata"
    | "clock"
    | "aaDist"
    | "model"
    | "icode"
    | "Mgene"
    | "fix_kappa"
    | "fix_omega"
    | "fix_alpha"
    | "Malpha"
    | "ncatG"
    | "getSE"
    | "RateAncestor"
    | "cleandata"
    | "fix_blength"
    | "method"
    | "fix_rho" => true
    _ => false
  }
}

///|
fn codeml_is_double_option(name : String) -> Bool {
  name == "kappa" ||
  name == "omega" ||
  name == "alpha" ||
  name == "Small_Diff" ||
  name == "rho"
}

///|
fn codeml_validate_path(
  value : String,
  label : String,
) -> Unit raise CodemlError {
  if value.length() == 0 {
    codeml_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() || code == 0 {
      codeml_fail(label + " contains an invalid character")
    }
  }
}

///|
fn codeml_validate_option(
  name : String,
  value : String,
) -> Unit raise CodemlError {
  if !codeml_is_allowed_option(name) {
    codeml_fail("Unknown CODEML option '" + name + "'")
  }
  if value.length() == 0 {
    codeml_fail("CODEML option '" + name + "' has an empty value")
  }
  if codeml_is_integer_option(name) {
    ignore(codeml_parse_int(value, "CODEML option '" + name + "'"))
  } else if codeml_is_double_option(name) {
    ignore(codeml_parse_double(value, "CODEML option '" + name + "'"))
  } else if name == "NSsites" {
    let values = codeml_words(value)
    if values.length() == 0 {
      codeml_fail("CODEML NSsites must contain at least one class")
    }
    for item in values {
      let site = codeml_parse_int(item, "CODEML NSsites class")
      if site < 0 {
        codeml_fail("CODEML NSsites classes must be non-negative")
      }
    }
  }
}

///|
fn codeml_has_option(options : Array[CodemlOption], name : String) -> Bool {
  for option in options {
    if option.name == name {
      return true
    }
  }
  false
}

///|
/// Parse and validate CODEML control text.
///
/// PAML comments begin with `*` and extend to the end of the line.
pub fn codeml_parse_control(text : String) -> CodemlControl raise CodemlError {
  let mut sequence_file = ""
  let mut output_file = ""
  let mut tree_file = ""
  let options : Array[CodemlOption] = []
  let lines = codeml_lines(text)
  for line_index in 0..<lines.length() {
    let raw = lines[line_index]
    let comment = codeml_index_of_char(raw, '*'.to_int(), 0)
    let uncommented = if comment >= 0 { raw[0:comment].to_owned() } else { raw }
    let line = codeml_trim(uncommented)
    if line.length() == 0 {
      continue
    }
    let equals = codeml_index_of_char(line, '='.to_int(), 0)
    if equals <= 0 {
      codeml_fail(
        "Malformed CODEML control line " + (line_index + 1).to_string(),
      )
    }
    let name = codeml_trim(line[0:equals].to_owned())
    let value = codeml_trim(line[equals + 1:line.length()].to_owned())
    if name.length() == 0 || value.length() == 0 {
      codeml_fail(
        "Malformed CODEML control line " + (line_index + 1).to_string(),
      )
    }
    if name == "seqfile" {
      if sequence_file.length() > 0 {
        codeml_fail("Duplicate CODEML seqfile entry")
      }
      codeml_validate_path(value, "CODEML seqfile")
      sequence_file = value
    } else if name == "outfile" {
      if output_file.length() > 0 {
        codeml_fail("Duplicate CODEML outfile entry")
      }
      codeml_validate_path(value, "CODEML outfile")
      output_file = value
    } else if name == "treefile" {
      if tree_file.length() > 0 {
        codeml_fail("Duplicate CODEML treefile entry")
      }
      codeml_validate_path(value, "CODEML treefile")
      tree_file = value
    } else {
      if codeml_has_option(options, name) {
        codeml_fail("Duplicate CODEML option '" + name + "'")
      }
      codeml_validate_option(name, value)
      options.push(CodemlOption::{ name, value })
    }
  }
  if sequence_file.length() == 0 {
    codeml_fail("CODEML control is missing seqfile")
  }
  if output_file.length() == 0 {
    codeml_fail("CODEML control is missing outfile")
  }
  if tree_file.length() == 0 {
    codeml_fail("CODEML control is missing treefile")
  }
  CodemlControl::{ sequence_file, output_file, tree_file, options }
}

///|
/// Write a canonical CODEML control file.
pub fn codeml_write_control(control : CodemlControl) -> String {
  let output = StringBuilder::new()
  output.write_string("seqfile = " + control.sequence_file + "\n")
  output.write_string("outfile = " + control.output_file + "\n")
  output.write_string("treefile = " + control.tree_file + "\n")
  for option in control.options {
    output.write_string(option.name + " = " + option.value + "\n")
  }
  output.to_string()
}

///|
pub fn CodemlControl::sequence_file(self : CodemlControl) -> String {
  self.sequence_file
}

///|
pub fn CodemlControl::output_file(self : CodemlControl) -> String {
  self.output_file
}

///|
pub fn CodemlControl::tree_file(self : CodemlControl) -> String {
  self.tree_file
}

///|
pub fn CodemlControl::options(self : CodemlControl) -> Array[CodemlOption] {
  codeml_copy_options(self.options)
}

///|
pub fn CodemlControl::option(self : CodemlControl, name : String) -> String? {
  for option in self.options {
    if option.name == name {
      return Some(option.value)
    }
  }
  None
}

///|
pub fn CodemlControl::ns_sites(self : CodemlControl) -> Array[Int] {
  let result : Array[Int] = []
  match self.option("NSsites") {
    Some(value) =>
      for word in codeml_words(value) {
        let number = codeml_parse_int(word, "CODEML NSsites class") catch {
          _ => 0
        }
        result.push(number)
      }
    None => ()
  }
  result
}

///|
pub fn CodemlOption::name(self : CodemlOption) -> String {
  self.name
}

///|
pub fn CodemlOption::value(self : CodemlOption) -> String {
  self.value
}

///|
fn codeml_site_model_number(name : String) -> Int {
  match name {
    "" | "one-ratio" => 0
    "NearlyNeutral" => 1
    "PositiveSelection" => 2
    "discrete" => 3
    "beta" => 7
    "beta&w>1" => 8
    "M2a_rel" => 22
    _ => 0
  }
}

///|
fn codeml_model_header(line : String) -> (Int, String)? {
  let value = codeml_trim(line)
  if !codeml_starts_with(value, "Model ") {
    return None
  }
  let colon = codeml_index_of_char(value, ':'.to_int(), 6)
  if colon < 0 {
    return None
  }
  let number_text = codeml_trim(value[6:colon].to_owned())
  if number_text.length() == 0 {
    return None
  }
  let number = codeml_parse_int(number_text, "CODEML model number") catch {
    _ => return None
  }
  Some((number, codeml_trim(value[colon + 1:value.length()].to_owned())))
}

///|
fn codeml_gene_header(line : String) -> Int? {
  let value = codeml_trim(line)
  if !codeml_starts_with(value, "Gene ") {
    return None
  }
  let words = codeml_words(value)
  if words.length() < 2 {
    return None
  }
  let raw = if codeml_ends_with(words[1], ":") {
    words[1][0:words[1].length() - 1].to_owned()
  } else {
    words[1]
  }
  let number = codeml_parse_int(raw, "CODEML gene number") catch {
    _ => return None
  }
  Some(number)
}

///|
fn codeml_set_site_proportions(
  values : Array[Double],
) -> Array[CodemlSiteClass] {
  let classes : Array[CodemlSiteClass] = []
  for index in 0..<values.length() {
    classes.push(CodemlSiteClass::{
      index,
      proportion: values[index],
      omega: None,
      branch_types: [],
    })
  }
  classes
}

///|
fn codeml_set_site_omegas(
  classes : Array[CodemlSiteClass],
  values : Array[Double],
) -> Array[CodemlSiteClass] {
  let updated : Array[CodemlSiteClass] = []
  for index in 0..<classes.length() {
    updated.push(CodemlSiteClass::{
      index: classes[index].index,
      proportion: classes[index].proportion,
      omega: if index < values.length() {
        Some(values[index])
      } else {
        classes[index].omega
      },
      branch_types: codeml_copy_named_values(classes[index].branch_types),
    })
  }
  updated
}

///|
fn codeml_set_site_branch_type(
  classes : Array[CodemlSiteClass],
  name : String,
  values : Array[Double],
) -> Array[CodemlSiteClass] {
  let updated : Array[CodemlSiteClass] = []
  for index in 0..<classes.length() {
    let branch_types = codeml_copy_named_values(classes[index].branch_types)
    if index < values.length() {
      branch_types.push(CodemlNamedValue::{ name, value: values[index] })
    }
    updated.push(CodemlSiteClass::{
      index: classes[index].index,
      proportion: classes[index].proportion,
      omega: classes[index].omega,
      branch_types,
    })
  }
  updated
}

///|
fn codeml_parse_named_parameters(line : String) -> Array[CodemlNamedValue] {
  let output : Array[CodemlNamedValue] = []
  let names = ["p0", "p1", "p", "q", "w"]
  for name in names {
    let mut from = 0
    while from < line.length() {
      let position = codeml_index_of(line, name, from)
      if position < 0 {
        break
      }
      let before_ok = position == 0 ||
        codeml_is_space(line.unsafe_get(position - 1).to_int()) ||
        line.unsafe_get(position - 1).to_int() == '('.to_int()
      let mut equals = position + name.length()
      while equals < line.length() &&
            codeml_is_space(line.unsafe_get(equals).to_int()) {
        equals = equals + 1
      }
      if before_ok &&
        equals < line.length() &&
        line.unsafe_get(equals).to_int() == '='.to_int() {
        let numbers = codeml_scan_numbers(
          line[equals + 1:line.length()].to_owned(),
        )
        if numbers.length() > 0 {
          output.push(CodemlNamedValue::{ name, value: numbers[0] })
          break
        }
      }
      from = position + name.length()
    }
  }
  output
}

///|
fn codeml_parse_positive_site(
  line : String,
  analysis_method : String,
) -> CodemlPositiveSite? {
  if analysis_method.length() == 0 {
    return None
  }
  let words = codeml_words(codeml_trim(line))
  if words.length() < 3 {
    return None
  }
  let position = codeml_parse_int(words[0], "CODEML site position") catch {
    _ => return None
  }
  if words[1].length() != 1 {
    return None
  }
  let mut probability_text = words[2]
  let mut significance = ""
  while codeml_ends_with(probability_text, "*") {
    significance = significance + "*"
    probability_text = probability_text[0:probability_text.length() - 1].to_owned()
  }
  let probability = codeml_parse_double(
    probability_text, "CODEML site probability",
  ) catch {
    _ => return None
  }
  let numbers = codeml_scan_numbers(line)
  let posterior = if numbers.length() >= 3 { Some(numbers[2]) } else { None }
  let standard_error = if numbers.length() >= 4 {
    Some(numbers[3])
  } else {
    None
  }
  Some(CodemlPositiveSite::{
    analysis_method,
    position,
    amino_acid: words[1],
    probability,
    significance,
    posterior_mean_omega: posterior,
    standard_error,
  })
}

///|
fn codeml_parse_model_segment(
  lines : Array[String],
  start : Int,
  end : Int,
  number : Int,
  gene : Int,
  description : String,
) -> CodemlModelResult? raise CodemlError {
  let mut ln_likelihood : Double? = None
  let mut parameter_count = 0
  let mut parameter_list = ""
  let mut standard_errors = ""
  let mut tree_length : Double? = None
  let mut tree = ""
  let mut dn_tree = ""
  let mut ds_tree = ""
  let mut omega_tree = ""
  let mut kappa : Double? = None
  let omegas : Array[Double] = []
  let mut dn_tree_length : Double? = None
  let mut ds_tree_length : Double? = None
  let rates : Array[Double] = []
  let gene_parameters : Array[CodemlGeneParameter] = []
  let named_parameters : Array[CodemlNamedValue] = []
  let mut site_classes : Array[CodemlSiteClass] = []
  let branches : Array[CodemlBranchResult] = []
  let positive_sites : Array[CodemlPositiveSite] = []
  let mut expect_parameter_list = false
  let mut expect_standard_errors = false
  let mut expect_dn_tree = false
  let mut expect_ds_tree = false
  let mut expect_omega_tree = false
  let mut positive_method = ""
  for line_index in start..<end {
    let line = lines[line_index]
    let trimmed = codeml_trim(line)
    if codeml_contains(trimmed, "lnL(ntime:") {
      let np = codeml_index_of(trimmed, "np:", 0)
      if np >= 0 {
        let close = codeml_index_of_char(trimmed, ')'.to_int(), np)
        if close > np {
          parameter_count = codeml_parse_int(
            trimmed[np + 3:close].to_owned(),
            "CODEML parameter count",
          )
        }
      }
      let colon = codeml_last_index_of_char(
        if codeml_contains(trimmed, "):") {
          trimmed[0:codeml_index_of(trimmed, "+", 0).max(0)].to_owned()
        } else {
          trimmed
        },
        ':'.to_int(),
      )
      let numbers = if colon >= 0 {
        codeml_scan_numbers(trimmed[colon + 1:trimmed.length()].to_owned())
      } else {
        []
      }
      if numbers.length() == 0 {
        codeml_fail(
          "CODEML lnL line " + (line_index + 1).to_string() + " is malformed",
        )
      }
      ln_likelihood = Some(numbers[0])
      expect_parameter_list = parameter_count > 0
      continue
    }
    if expect_parameter_list && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() == parameter_count {
        parameter_list = trimmed
        expect_parameter_list = false
        continue
      }
    }
    if codeml_starts_with(trimmed, "SEs for parameters:") {
      expect_standard_errors = true
      continue
    }
    if expect_standard_errors && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() == parameter_count {
        standard_errors = trimmed
        expect_standard_errors = false
        continue
      }
    }
    if codeml_starts_with(trimmed, "tree length for dN:") {
      let values = codeml_scan_numbers(
        trimmed[codeml_index_of_char(trimmed, ':'.to_int(), 0) + 1:trimmed.length()].to_owned(),
      )
      if values.length() > 0 {
        dn_tree_length = Some(values[0])
      }
      continue
    }
    if codeml_starts_with(trimmed, "tree length for dS:") {
      let values = codeml_scan_numbers(
        trimmed[codeml_index_of_char(trimmed, ':'.to_int(), 0) + 1:trimmed.length()].to_owned(),
      )
      if values.length() > 0 {
        ds_tree_length = Some(values[0])
      }
      continue
    }
    if codeml_starts_with(trimmed, "tree length =") {
      let values = codeml_scan_numbers(
        trimmed[codeml_index_of_char(trimmed, '='.to_int(), 0) + 1:trimmed.length()].to_owned(),
      )
      if values.length() > 0 {
        tree_length = Some(values[0])
      }
      continue
    }
    if trimmed == "dS tree:" {
      expect_ds_tree = true
      continue
    }
    if trimmed == "dN tree:" {
      expect_dn_tree = true
      continue
    }
    if codeml_starts_with(trimmed, "w ratios as labels for TreeView:") {
      expect_omega_tree = true
      continue
    }
    if codeml_starts_with(trimmed, "(") &&
      codeml_ends_with(trimmed, ";") &&
      (codeml_contains(trimmed, ":") || codeml_contains(trimmed, "#")) {
      if expect_ds_tree {
        ds_tree = trimmed
        expect_ds_tree = false
      } else if expect_dn_tree {
        dn_tree = trimmed
        expect_dn_tree = false
      } else if expect_omega_tree {
        omega_tree = trimmed
        expect_omega_tree = false
      } else {
        tree = trimmed
      }
      continue
    }
    if codeml_contains(trimmed, "kappa (ts/tv)") {
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      let values = if equals >= 0 {
        codeml_scan_numbers(trimmed[equals + 1:trimmed.length()].to_owned())
      } else {
        []
      }
      if values.length() > 0 {
        kappa = Some(values[0])
      }
      continue
    }
    if codeml_starts_with(trimmed, "omega (dN/dS)") ||
      codeml_starts_with(trimmed, "w (dN/dS)") {
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), 0)
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      let split = if colon >= 0 { colon } else { equals }
      let values = if split >= 0 {
        codeml_scan_numbers(trimmed[split + 1:trimmed.length()].to_owned())
      } else {
        []
      }
      for value in values {
        omegas.push(value)
      }
      continue
    }
    if codeml_starts_with(trimmed, "rates for") {
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), 0)
      if colon >= 0 {
        let values = codeml_scan_numbers(
          trimmed[colon + 1:trimmed.length()].to_owned(),
        )
        if values.length() > 0 && values[0] != 1.0 {
          rates.push(1.0)
        }
        for value in values {
          rates.push(value)
        }
      }
      continue
    }
    if codeml_starts_with(trimmed, "gene #") {
      let hash = codeml_index_of_char(trimmed, '#'.to_int(), 0)
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), hash)
      if hash >= 0 && colon > hash {
        let gene_number = codeml_parse_int(
          trimmed[hash + 1:colon].to_owned(),
          "CODEML gene parameter number",
        )
        let values = codeml_scan_numbers(
          trimmed[colon + 1:trimmed.length()].to_owned(),
        )
        if values.length() >= 2 {
          gene_parameters.push(CodemlGeneParameter::{
            gene: gene_number,
            kappa: values[0],
            omega: values[1],
          })
        }
      }
      continue
    }
    if codeml_starts_with(trimmed, "p:") ||
      codeml_starts_with(trimmed, "proportion") {
      let values = if codeml_starts_with(trimmed, "p:") {
        codeml_scan_numbers(trimmed[2:trimmed.length()].to_owned())
      } else {
        codeml_scan_numbers(trimmed[10:trimmed.length()].to_owned())
      }
      site_classes = codeml_set_site_proportions(values)
      continue
    }
    if codeml_starts_with(trimmed, "w:") {
      let mut values = codeml_scan_fixed_five_numbers(
        trimmed[2:trimmed.length()].to_owned(),
      )
      if values.length() == 0 {
        values = codeml_scan_numbers(trimmed[2:trimmed.length()].to_owned())
      }
      site_classes = codeml_set_site_omegas(site_classes, values)
      continue
    }
    if codeml_starts_with(trimmed, "background w") ||
      codeml_starts_with(trimmed, "foreground w") {
      let foreground = codeml_starts_with(trimmed, "foreground")
      let values = codeml_scan_numbers(trimmed[12:trimmed.length()].to_owned())
      site_classes = codeml_set_site_branch_type(
        site_classes,
        if foreground {
          "foreground"
        } else {
          "background"
        },
        values,
      )
      continue
    }
    if codeml_starts_with(trimmed, "branch type ") {
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), 0)
      if colon >= 0 {
        let type_text = codeml_trim(trimmed[12:colon].to_owned())
        let values = codeml_scan_numbers(
          trimmed[colon + 1:trimmed.length()].to_owned(),
        )
        site_classes = codeml_set_site_branch_type(
          site_classes, type_text, values,
        )
      }
      continue
    }
    let words = codeml_words(trimmed)
    if words.length() == 9 && codeml_contains(words[0], "..") {
      branches.push(CodemlBranchResult::{
        branch: words[0],
        time: codeml_parse_optional_double(words[1]),
        nonsynonymous_sites: codeml_parse_optional_double(words[2]),
        synonymous_sites: codeml_parse_optional_double(words[3]),
        omega: codeml_parse_optional_double(words[4]),
        dn: codeml_parse_optional_double(words[5]),
        ds: codeml_parse_optional_double(words[6]),
        n_times_dn: codeml_parse_optional_double(words[7]),
        s_times_ds: codeml_parse_optional_double(words[8]),
      })
      continue
    }
    if codeml_starts_with(trimmed, "Naive Empirical Bayes") {
      positive_method = "NEB"
      continue
    }
    if codeml_starts_with(trimmed, "Bayes Empirical Bayes") {
      positive_method = "BEB"
      continue
    }
    match codeml_parse_positive_site(trimmed, positive_method) {
      Some(site) => {
        positive_sites.push(site)
        continue
      }
      None => ()
    }
    let parsed_named = codeml_parse_named_parameters(trimmed)
    for value in parsed_named {
      let mut seen = false
      for existing in named_parameters {
        if existing.name == value.name {
          seen = true
          break
        }
      }
      if !seen {
        named_parameters.push(value)
      }
    }
  }
  match ln_likelihood {
    Some(value) =>
      Some(CodemlModelResult::{
        number,
        gene,
        description,
        ln_likelihood: value,
        parameter_count,
        parameter_list,
        standard_errors,
        tree_length,
        tree,
        dn_tree,
        ds_tree,
        omega_tree,
        kappa,
        omegas,
        dn_tree_length,
        ds_tree_length,
        rates,
        gene_parameters,
        named_parameters,
        site_classes,
        branches,
        positive_sites,
      })
    None => None
  }
}

///|
fn codeml_parse_pairwise(
  lines : Array[String],
) -> Array[CodemlPairwiseResult] raise CodemlError {
  let output : Array[CodemlPairwiseResult] = []
  let mut pending : CodemlPendingPair? = None
  for line_index in 0..<lines.length() {
    let trimmed = codeml_trim(lines[line_index])
    let first_open = codeml_index_of_char(trimmed, '('.to_int(), 0)
    let first_close = if first_open >= 0 {
      codeml_index_of_char(trimmed, ')'.to_int(), first_open + 1)
    } else {
      -1
    }
    let dots = codeml_index_of(trimmed, "...", 0)
    let second_open = if dots >= 0 {
      codeml_index_of_char(trimmed, '('.to_int(), dots + 3)
    } else {
      -1
    }
    let second_close = if second_open >= 0 {
      codeml_index_of_char(trimmed, ')'.to_int(), second_open + 1)
    } else {
      -1
    }
    if first_open > 0 &&
      first_close > first_open &&
      dots > first_close &&
      second_open > dots &&
      second_close > second_open {
      pending = Some(CodemlPendingPair::{
        first: trimmed[first_open + 1:first_close].to_owned(),
        second: trimmed[second_open + 1:second_close].to_owned(),
        ln_likelihood: None,
      })
      continue
    }
    match pending {
      Some(pair) =>
        if codeml_starts_with(trimmed, "lnL =") {
          let values = codeml_scan_numbers(
            trimmed[5:trimmed.length()].to_owned(),
          )
          if values.length() != 1 {
            codeml_fail(
              "Malformed CODEML pairwise lnL at line " +
              (line_index + 1).to_string(),
            )
          }
          pair.ln_likelihood = Some(values[0])
        } else if codeml_starts_with(trimmed, "t=") {
          let values = codeml_scan_numbers(trimmed)
          if values.length() != 6 || pair.ln_likelihood is None {
            codeml_fail(
              "Malformed CODEML pairwise estimate at line " +
              (line_index + 1).to_string(),
            )
          }
          output.push(CodemlPairwiseResult::{
            first: pair.first,
            second: pair.second,
            ln_likelihood: pair.ln_likelihood.unwrap(),
            time: values[0],
            synonymous_sites: values[1],
            nonsynonymous_sites: values[2],
            omega: values[3],
            dn: values[4],
            ds: values[5],
          })
          pending = None
        }
      None => ()
    }
  }
  output
}

///|
fn codeml_parse_distances(lines : Array[String]) -> Array[CodemlDistance] {
  let output : Array[CodemlDistance] = []
  let mut kind = ""
  let names : Array[String] = []
  let mut started = false
  for line in lines {
    let trimmed = codeml_trim(line)
    if codeml_contains(trimmed, "AA distances") {
      kind = "raw"
      names.clear()
      started = false
      continue
    }
    if codeml_starts_with(trimmed, "ML distances of aa seqs.") {
      kind = "ml"
      names.clear()
      started = false
      continue
    }
    if kind.length() == 0 {
      continue
    }
    if trimmed.length() == 0 {
      if started {
        kind = ""
      }
      continue
    }
    let words = codeml_words(trimmed)
    if words.length() == 0 {
      continue
    }
    let row_name = words[0]
    let values : Array[Double] = []
    let mut valid = true
    for index in 1..<words.length() {
      if codeml_is_number(words[index]) {
        match parse_double(words[index]) {
          Some(value) => values.push(value)
          None => valid = false
        }
      } else {
        valid = false
      }
    }
    if !valid || values.length() > names.length() {
      kind = ""
      continue
    }
    for index in 0..<values.length() {
      output.push(CodemlDistance::{
        kind,
        first: row_name,
        second: names[index],
        value: values[index],
      })
    }
    names.push(row_name)
    started = true
  }
  output
}

///|
fn codeml_parse_counts(line : String) -> (Int, Int)? {
  let ns = codeml_index_of(line, "ns =", 0)
  let ls = codeml_index_of(line, "ls =", 0)
  if ns < 0 || ls < 0 || ls <= ns {
    return None
  }
  let left = codeml_scan_numbers(line[ns + 4:ls].to_owned())
  let right = codeml_scan_numbers(line[ls + 4:line.length()].to_owned())
  if left.length() == 0 || right.length() == 0 {
    None
  } else {
    Some((left[0].to_int(), right[0].to_int()))
  }
}

///|
/// Parse CODONML or AAML result text.
pub fn codeml_parse_results(text : String) -> CodemlResults raise CodemlError {
  if codeml_trim(text).length() == 0 {
    codeml_fail("CODEML results are empty")
  }
  let lines = codeml_lines(text)
  let mut program = ""
  let mut version = ""
  let mut model = ""
  let mut codon_frequency_model = ""
  let mut site_class_model = ""
  let mut sequence_count = 0
  let mut site_count = 0
  let mut ln_l_max : Double? = None
  for line in lines {
    let trimmed = codeml_trim(line)
    if program.length() == 0 &&
      (
        codeml_starts_with(trimmed, "CODONML ") ||
        codeml_starts_with(trimmed, "AAML ")
      ) {
      let words = codeml_words(trimmed)
      program = words[0]
      let marker = codeml_index_of(trimmed, "(in paml version ", 0)
      if marker >= 0 {
        let start = marker + 17
        let comma = codeml_index_of_char(trimmed, ','.to_int(), start)
        let close = codeml_index_of_char(trimmed, ')'.to_int(), start)
        let end = if comma >= 0 { comma } else { close }
        if end > start {
          version = trimmed[start:end].to_owned()
        }
      }
    } else if codeml_starts_with(trimmed, "Model:") {
      model = codeml_trim(trimmed[6:trimmed.length()].to_owned())
      match codeml_parse_counts(trimmed) {
        Some((sequences, sites)) => {
          sequence_count = sequences
          site_count = sites
        }
        None => ()
      }
    } else if codeml_starts_with(trimmed, "Codon frequency model:") {
      codon_frequency_model = codeml_trim(
        trimmed[22:trimmed.length()].to_owned(),
      )
    } else if codeml_starts_with(trimmed, "Codon frequencies:") {
      codon_frequency_model = codeml_trim(
        trimmed[18:trimmed.length()].to_owned(),
      )
    } else if codeml_starts_with(trimmed, "Site-class models:") {
      site_class_model = codeml_trim(trimmed[18:trimmed.length()].to_owned())
    } else if codeml_contains(trimmed, "ns =") &&
      codeml_contains(trimmed, "ls =") {
      match codeml_parse_counts(trimmed) {
        Some((sequences, sites)) => {
          sequence_count = sequences
          site_count = sites
        }
        None => ()
      }
    } else if codeml_contains(trimmed, "ln Lmax") {
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      if equals >= 0 {
        let values = codeml_scan_numbers(
          trimmed[equals + 1:trimmed.length()].to_owned(),
        )
        if values.length() > 0 {
          ln_l_max = Some(values[0])
        }
      }
    }
  }
  if program.length() == 0 || version.length() == 0 {
    codeml_fail("Invalid CODEML results header")
  }
  let models : Array[CodemlModelResult] = []
  let model_starts : Array[Int] = []
  let model_numbers : Array[Int] = []
  let model_descriptions : Array[String] = []
  for index in 0..<lines.length() {
    match codeml_model_header(lines[index]) {
      Some((number, description)) => {
        model_starts.push(index)
        model_numbers.push(number)
        model_descriptions.push(description)
      }
      None => ()
    }
  }
  if model_starts.length() > 0 {
    for index in 0..<model_starts.length() {
      let end = if index + 1 < model_starts.length() {
        model_starts[index + 1]
      } else {
        lines.length()
      }
      match
        codeml_parse_model_segment(
          lines,
          model_starts[index],
          end,
          model_numbers[index],
          0,
          model_descriptions[index],
        ) {
        Some(result) => models.push(result)
        None =>
          codeml_fail(
            "CODEML model " +
            model_numbers[index].to_string() +
            " has no likelihood",
          )
      }
    }
  } else {
    let gene_starts : Array[Int] = []
    let gene_numbers : Array[Int] = []
    for index in 0..<lines.length() {
      match codeml_gene_header(lines[index]) {
        Some(gene) => {
          gene_starts.push(index)
          gene_numbers.push(gene)
        }
        None => ()
      }
    }
    if gene_starts.length() > 0 {
      for index in 0..<gene_starts.length() {
        let end = if index + 1 < gene_starts.length() {
          gene_starts[index + 1]
        } else {
          lines.length()
        }
        match
          codeml_parse_model_segment(
            lines,
            gene_starts[index],
            end,
            codeml_site_model_number(site_class_model),
            gene_numbers[index],
            if site_class_model.length() > 0 {
              site_class_model
            } else {
              "one-ratio"
            },
          ) {
          Some(result) => models.push(result)
          None => ()
        }
      }
    } else {
      match
        codeml_parse_model_segment(
          lines,
          0,
          lines.length(),
          codeml_site_model_number(site_class_model),
          0,
          if site_class_model.length() > 0 {
            site_class_model
          } else {
            "one-ratio"
          },
        ) {
        Some(result) => models.push(result)
        None => ()
      }
    }
  }
  let pairwise = codeml_parse_pairwise(lines)
  let distances = codeml_parse_distances(lines)
  if models.length() == 0 &&
    pairwise.length() == 0 &&
    distances.length() == 0 &&
    ln_l_max is None {
    codeml_fail("CODEML results contain no recognized estimates")
  }
  CodemlResults::{
    program,
    version,
    model,
    codon_frequency_model,
    site_class_model,
    sequence_count,
    site_count,
    ln_l_max,
    models,
    pairwise,
    distances,
  }
}

///|
pub fn CodemlResults::program(self : CodemlResults) -> String {
  self.program
}

///|
pub fn CodemlResults::version(self : CodemlResults) -> String {
  self.version
}

///|
pub fn CodemlResults::model(self : CodemlResults) -> String {
  self.model
}

///|
pub fn CodemlResults::codon_frequency_model(self : CodemlResults) -> String {
  self.codon_frequency_model
}

///|
pub fn CodemlResults::site_class_model(self : CodemlResults) -> String {
  self.site_class_model
}

///|
pub fn CodemlResults::sequence_count(self : CodemlResults) -> Int {
  self.sequence_count
}

///|
pub fn CodemlResults::site_count(self : CodemlResults) -> Int {
  self.site_count
}

///|
pub fn CodemlResults::ln_l_max(self : CodemlResults) -> Double? {
  self.ln_l_max
}

///|
pub fn CodemlResults::models(self : CodemlResults) -> Array[CodemlModelResult] {
  codeml_copy_models(self.models)
}

///|
pub fn CodemlResults::pairwise(
  self : CodemlResults,
) -> Array[CodemlPairwiseResult] {
  codeml_copy_pairwise(self.pairwise)
}

///|
pub fn CodemlResults::distances(self : CodemlResults) -> Array[CodemlDistance] {
  codeml_copy_distances(self.distances)
}

///|
pub fn CodemlResults::model_by_number(
  self : CodemlResults,
  number : Int,
) -> CodemlModelResult? {
  for model in self.models {
    if model.number == number {
      return Some(codeml_clone_model(model))
    }
  }
  None
}

///|
pub fn CodemlModelResult::number(self : CodemlModelResult) -> Int {
  self.number
}

///|
pub fn CodemlModelResult::gene(self : CodemlModelResult) -> Int {
  self.gene
}

///|
pub fn CodemlModelResult::description(self : CodemlModelResult) -> String {
  self.description
}

///|
pub fn CodemlModelResult::ln_likelihood(self : CodemlModelResult) -> Double {
  self.ln_likelihood
}

///|
pub fn CodemlModelResult::parameter_count(self : CodemlModelResult) -> Int {
  self.parameter_count
}

///|
pub fn CodemlModelResult::parameter_list(self : CodemlModelResult) -> String {
  self.parameter_list
}

///|
pub fn CodemlModelResult::standard_errors(self : CodemlModelResult) -> String {
  self.standard_errors
}

///|
pub fn CodemlModelResult::tree_length(self : CodemlModelResult) -> Double? {
  self.tree_length
}

///|
pub fn CodemlModelResult::tree(self : CodemlModelResult) -> String {
  self.tree
}

///|
pub fn CodemlModelResult::dn_tree(self : CodemlModelResult) -> String {
  self.dn_tree
}

///|
pub fn CodemlModelResult::ds_tree(self : CodemlModelResult) -> String {
  self.ds_tree
}

///|
pub fn CodemlModelResult::omega_tree(self : CodemlModelResult) -> String {
  self.omega_tree
}

///|
pub fn CodemlModelResult::kappa(self : CodemlModelResult) -> Double? {
  self.kappa
}

///|
pub fn CodemlModelResult::omegas(self : CodemlModelResult) -> Array[Double] {
  codeml_copy_doubles(self.omegas)
}

///|
pub fn CodemlModelResult::dn_tree_length(self : CodemlModelResult) -> Double? {
  self.dn_tree_length
}

///|
pub fn CodemlModelResult::ds_tree_length(self : CodemlModelResult) -> Double? {
  self.ds_tree_length
}

///|
pub fn CodemlModelResult::rates(self : CodemlModelResult) -> Array[Double] {
  codeml_copy_doubles(self.rates)
}

///|
pub fn CodemlModelResult::gene_parameters(
  self : CodemlModelResult,
) -> Array[CodemlGeneParameter] {
  codeml_copy_gene_parameters(self.gene_parameters)
}

///|
pub fn CodemlModelResult::named_parameters(
  self : CodemlModelResult,
) -> Array[CodemlNamedValue] {
  codeml_copy_named_values(self.named_parameters)
}

///|
pub fn CodemlModelResult::site_classes(
  self : CodemlModelResult,
) -> Array[CodemlSiteClass] {
  codeml_copy_site_classes(self.site_classes)
}

///|
pub fn CodemlModelResult::branches(
  self : CodemlModelResult,
) -> Array[CodemlBranchResult] {
  codeml_copy_branches(self.branches)
}

///|
pub fn CodemlModelResult::positive_sites(
  self : CodemlModelResult,
) -> Array[CodemlPositiveSite] {
  codeml_copy_positive_sites(self.positive_sites)
}

///|
pub fn CodemlNamedValue::name(self : CodemlNamedValue) -> String {
  self.name
}

///|
pub fn CodemlNamedValue::value(self : CodemlNamedValue) -> Double {
  self.value
}

///|
pub fn CodemlGeneParameter::gene(self : CodemlGeneParameter) -> Int {
  self.gene
}

///|
pub fn CodemlGeneParameter::kappa(self : CodemlGeneParameter) -> Double {
  self.kappa
}

///|
pub fn CodemlGeneParameter::omega(self : CodemlGeneParameter) -> Double {
  self.omega
}

///|
pub fn CodemlSiteClass::index(self : CodemlSiteClass) -> Int {
  self.index
}

///|
pub fn CodemlSiteClass::proportion(self : CodemlSiteClass) -> Double {
  self.proportion
}

///|
pub fn CodemlSiteClass::omega(self : CodemlSiteClass) -> Double? {
  self.omega
}

///|
pub fn CodemlSiteClass::branch_types(
  self : CodemlSiteClass,
) -> Array[CodemlNamedValue] {
  codeml_copy_named_values(self.branch_types)
}

///|
pub fn CodemlBranchResult::branch(self : CodemlBranchResult) -> String {
  self.branch
}

///|
pub fn CodemlBranchResult::time(self : CodemlBranchResult) -> Double? {
  self.time
}

///|
pub fn CodemlBranchResult::nonsynonymous_sites(
  self : CodemlBranchResult,
) -> Double? {
  self.nonsynonymous_sites
}

///|
pub fn CodemlBranchResult::synonymous_sites(
  self : CodemlBranchResult,
) -> Double? {
  self.synonymous_sites
}

///|
pub fn CodemlBranchResult::omega(self : CodemlBranchResult) -> Double? {
  self.omega
}

///|
pub fn CodemlBranchResult::dn(self : CodemlBranchResult) -> Double? {
  self.dn
}

///|
pub fn CodemlBranchResult::ds(self : CodemlBranchResult) -> Double? {
  self.ds
}

///|
pub fn CodemlBranchResult::n_times_dn(self : CodemlBranchResult) -> Double? {
  self.n_times_dn
}

///|
pub fn CodemlBranchResult::s_times_ds(self : CodemlBranchResult) -> Double? {
  self.s_times_ds
}

///|
pub fn CodemlPositiveSite::analysis_method(self : CodemlPositiveSite) -> String {
  self.analysis_method
}

///|
pub fn CodemlPositiveSite::position(self : CodemlPositiveSite) -> Int {
  self.position
}

///|
pub fn CodemlPositiveSite::amino_acid(self : CodemlPositiveSite) -> String {
  self.amino_acid
}

///|
pub fn CodemlPositiveSite::probability(self : CodemlPositiveSite) -> Double {
  self.probability
}

///|
pub fn CodemlPositiveSite::significance(self : CodemlPositiveSite) -> String {
  self.significance
}

///|
pub fn CodemlPositiveSite::posterior_mean_omega(
  self : CodemlPositiveSite,
) -> Double? {
  self.posterior_mean_omega
}

///|
pub fn CodemlPositiveSite::standard_error(self : CodemlPositiveSite) -> Double? {
  self.standard_error
}

///|
pub fn CodemlPairwiseResult::first(self : CodemlPairwiseResult) -> String {
  self.first
}

///|
pub fn CodemlPairwiseResult::second(self : CodemlPairwiseResult) -> String {
  self.second
}

///|
pub fn CodemlPairwiseResult::ln_likelihood(
  self : CodemlPairwiseResult,
) -> Double {
  self.ln_likelihood
}

///|
pub fn CodemlPairwiseResult::time(self : CodemlPairwiseResult) -> Double {
  self.time
}

///|
pub fn CodemlPairwiseResult::synonymous_sites(
  self : CodemlPairwiseResult,
) -> Double {
  self.synonymous_sites
}

///|
pub fn CodemlPairwiseResult::nonsynonymous_sites(
  self : CodemlPairwiseResult,
) -> Double {
  self.nonsynonymous_sites
}

///|
pub fn CodemlPairwiseResult::omega(self : CodemlPairwiseResult) -> Double {
  self.omega
}

///|
pub fn CodemlPairwiseResult::dn(self : CodemlPairwiseResult) -> Double {
  self.dn
}

///|
pub fn CodemlPairwiseResult::ds(self : CodemlPairwiseResult) -> Double {
  self.ds
}

///|
pub fn CodemlDistance::kind(self : CodemlDistance) -> String {
  self.kind
}

///|
pub fn CodemlDistance::first(self : CodemlDistance) -> String {
  self.first
}

///|
pub fn CodemlDistance::second(self : CodemlDistance) -> String {
  self.second
}

///|
pub fn CodemlDistance::value(self : CodemlDistance) -> Double {
  self.value
}

///|
/// Akaike information criterion for a fitted model.
pub fn CodemlModelResult::aic(self : CodemlModelResult) -> Double {
  2.0 * self.parameter_count.to_double() - 2.0 * self.ln_likelihood
}

///|
/// Bayesian information criterion for a fitted model.
pub fn CodemlModelResult::bic(
  self : CodemlModelResult,
  observations : Int,
) -> Double raise CodemlError {
  if observations <= 0 {
    codeml_fail("CODEML BIC observations must be positive")
  }
  self.parameter_count.to_double() * @math.ln(observations.to_double()) -
  2.0 * self.ln_likelihood
}

///|
fn codeml_log_gamma(value : Double) -> Double {
  let coefficients = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  let shifted = value - 1.0
  let mut sum = coefficients[0]
  for index in 1..<coefficients.length() {
    sum = sum + coefficients[index] / (shifted + index.to_double())
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(sum)
}

///|
fn codeml_clamp(value : Double, lower : Double, upper : Double) -> Double {
  if value < lower {
    lower
  } else if value > upper {
    upper
  } else {
    value
  }
}

///|
fn codeml_gamma_q(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 1.0
  }
  if x < a + 1.0 {
    let mut term = 1.0 / a
    let mut sum = term
    let mut index = 1
    while index < 300 {
      term = term * x / (a + index.to_double())
      sum = sum + term
      if term.abs() < sum.abs() * 1.0e-12 {
        break
      }
      index = index + 1
    }
    return codeml_clamp(
      1.0 - sum * @math.exp(-x + a * @math.ln(x) - codeml_log_gamma(a)),
      0.0,
      1.0,
    )
  }
  let tiny = 1.0e-300
  let mut b = x + 1.0 - a
  let mut c = 1.0 / tiny
  let mut d = 1.0 / b
  let mut value = d
  let mut index = 1
  while index < 300 {
    let an = -index.to_double() * (index.to_double() - a)
    b = b + 2.0
    d = an * d + b
    if d.abs() < tiny {
      d = tiny
    }
    c = b + an / c
    if c.abs() < tiny {
      c = tiny
    }
    d = 1.0 / d
    let delta = d * c
    value = value * delta
    if (delta - 1.0).abs() < 1.0e-12 {
      break
    }
    index = index + 1
  }
  codeml_clamp(
    value * @math.exp(-x + a * @math.ln(x) - codeml_log_gamma(a)),
    0.0,
    1.0,
  )
}

///|
/// Compare nested CODEML models with a chi-square likelihood-ratio test.
pub fn codeml_likelihood_ratio(
  null_model : CodemlModelResult,
  alternative_model : CodemlModelResult,
  alpha? : Double = 0.05,
) -> CodemlLikelihoodRatio raise CodemlError {
  if alpha <= 0.0 || alpha >= 1.0 || alpha.is_nan() {
    codeml_fail("CODEML LRT alpha must be between zero and one")
  }
  let degrees = alternative_model.parameter_count - null_model.parameter_count
  if degrees <= 0 {
    codeml_fail(
      "CODEML alternative model must have more parameters than the null model",
    )
  }
  if alternative_model.ln_likelihood + 1.0e-10 < null_model.ln_likelihood {
    codeml_fail(
      "CODEML alternative model has a lower likelihood than the null model",
    )
  }
  let statistic = 2.0 *
    (alternative_model.ln_likelihood - null_model.ln_likelihood)
  let p_value = codeml_gamma_q(0.5 * degrees.to_double(), 0.5 * statistic)
  CodemlLikelihoodRatio::{
    statistic,
    degrees_of_freedom: degrees,
    p_value,
    alpha,
    significant: p_value < alpha,
  }
}

///|
pub fn CodemlLikelihoodRatio::statistic(self : CodemlLikelihoodRatio) -> Double {
  self.statistic
}

///|
pub fn CodemlLikelihoodRatio::degrees_of_freedom(
  self : CodemlLikelihoodRatio,
) -> Int {
  self.degrees_of_freedom
}

///|
pub fn CodemlLikelihoodRatio::p_value(self : CodemlLikelihoodRatio) -> Double {
  self.p_value
}

///|
pub fn CodemlLikelihoodRatio::alpha(self : CodemlLikelihoodRatio) -> Double {
  self.alpha
}

///|
pub fn CodemlLikelihoodRatio::significant(self : CodemlLikelihoodRatio) -> Bool {
  self.significant
}

///|
/// Select the model with the smallest AIC.
pub fn codeml_best_aic(
  models : Array[CodemlModelResult],
) -> CodemlModelResult raise CodemlError {
  if models.length() == 0 {
    codeml_fail("CODEML model collection must not be empty")
  }
  let mut best = models[0]
  let mut best_aic = best.aic()
  for index in 1..<models.length() {
    let candidate = models[index].aic()
    if candidate < best_aic {
      best = models[index]
      best_aic = candidate
    }
  }
  codeml_clone_model(best)
}

///|
/// Compact official-style control fixture for examples and tests.
pub fn codeml_example_control_text() -> String {
  "seqfile = alignment.phylip\n" +
  "outfile = codeml.out\n" +
  "treefile = species.tree\n" +
  "verbose = 0\n" +
  "runmode = 0\n" +
  "seqtype = 1\n" +
  "CodonFreq = 2\n" +
  "model = 0\n" +
  "NSsites = 0 1 2\n" +
  "fix_kappa = 0\n" +
  "kappa = 4.54006\n" +
  "fix_omega = 0\n" +
  "omega = 1\n"
}

///|
/// Compact result fixture adapted from Biopython's PAML 4.7 test output.
pub fn codeml_example_results_text() -> String {
  "CODONML (in paml version 4.7, January 2013)  alignment.phylip\n" +
  "Model: One dN/dS ratio for branches,\n" +
  "Codon frequency model: F3x4\n" +
  "Site-class models:\n" +
  "ns =   5  ls =  74\n" +
  "\n" +
  "Model 0: one-ratio\n" +
  "lnL(ntime:  7  np:  9):   -308.032865      +0.000000\n" +
  "0.000004 0.000004 0.000005 0.013437 0.013431 0.000004 0.027837 1.519192 0.000100\n" +
  "tree length =   0.05472\n" +
  "(((Homo_sapie: 0.00000, Pan_troglo: 0.01344): 0.00000, Gorilla_go: 0.01343): 0.00000, Pongo_pygm: 0.00000, Macaca_mul: 0.02784);\n" +
  "kappa (ts/tv) =  1.51919\n" +
  "omega (dN/dS) =  0.00010\n" +
  " branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS\n" +
  " 6..7 0.000 167.7 54.3 0.0001 0.0000 0.0000 0.0 0.0\n" +
  " 8..2 0.013 167.7 54.3 0.0001 0.0000 0.0183 0.0 1.0\n" +
  "tree length for dN: 0.0000\n" +
  "tree length for dS: 0.0746\n" +
  "\n" +
  "Model 2: PositiveSelection (3 categories)\n" +
  "lnL(ntime:  7  np: 12):   -302.031579      +0.000000\n" +
  "0.000004 0.000004 0.000004 0.013438 0.013431 0.000004 0.027839 1.518255 1.000000 0.000000 0.000001 1.000000\n" +
  "tree length =   0.05472\n" +
  "(((Homo_sapie: 0.00000, Pan_troglo: 0.01344): 0.00000, Gorilla_go: 0.01343): 0.00000, Pongo_pygm: 0.00000, Macaca_mul: 0.02784);\n" +
  "kappa (ts/tv) =  1.51826\n" +
  "p:   0.90000  0.08000  0.02000\n" +
  "w:   0.10000  1.00000  5.00000\n" +
  " 6..7 0.000 167.7 54.3 0.1000 0.0000 0.0000 0.0 0.0\n" +
  "Bayes Empirical Bayes (BEB) analysis\n" +
  " 42 K 0.997** 5.200 +- 1.100\n"
}
