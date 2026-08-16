// Biopython-compatible Bio.Phylo.PAML.baseml support.
//
// This module parses BASEML control files and result text. It deliberately
// does not launch the external PAML executable.

///|
/// Error raised for malformed BASEML configuration, output, or comparisons.
pub suberror BasemlError {
  BasemlError(String)
}

///|
/// One validated BASEML option retained in canonical text form.
pub struct BasemlOption {
  name : String
  value : String
} derive(Eq, Debug)

///|
/// A BASEML control file.
pub struct BasemlControl {
  sequence_file : String
  output_file : String
  tree_file : String
  options : Array[BasemlOption]
} derive(Eq, Debug)

///|
/// Nucleotide frequencies in PAML's T, C, A, G order.
pub struct BasemlBaseFrequencies {
  thymine : Double
  cytosine : Double
  adenine : Double
  guanine : Double
} derive(Eq, Debug)

///|
/// A nucleotide rate matrix and its reported average transition/transversion
/// ratio.
pub struct BasemlQMatrix {
  average_ts_tv : Double?
  rows : Array[Array[Double]]
} derive(Eq, Debug)

///|
/// Branch-specific quantities emitted by a non-homogeneous kappa model.
pub struct BasemlBranchParameter {
  branch : String
  time : Double
  kappa : Double
  transitions : Double
  transversions : Double
} derive(Eq, Debug)

///|
/// Frequency parameters and optional realized frequencies for one node.
pub struct BasemlNodeFrequencies {
  node : Int
  root : Bool
  frequency_parameters : Array[Double]
  base_frequencies : BasemlBaseFrequencies?
} derive(Eq, Debug)

///|
/// Parsed BASEML maximum-likelihood result.
pub struct BasemlResults {
  version : String
  model_description : String
  ln_l_max : Double?
  ln_likelihood : Double
  parameter_count : Int
  parameter_list : Array[Double]
  standard_errors : Array[Double]
  tree_length : Double
  tree : String
  kappas : Array[Double]
  branch_parameters : Array[BasemlBranchParameter]
  rate_parameters : Array[Double]
  rates : Array[Double]
  rate_frequencies : Array[Double]
  q_matrix : BasemlQMatrix?
  alpha : Double?
  rho : Double?
  transition_probabilities : Array[Array[Double]]
  base_frequencies : BasemlBaseFrequencies?
  nodes : Array[BasemlNodeFrequencies]
} derive(Eq, Debug)

///|
/// Likelihood-ratio comparison between nested BASEML models.
pub struct BasemlLikelihoodRatio {
  statistic : Double
  degrees_of_freedom : Int
  p_value : Double
  alpha : Double
  significant : Bool
} derive(Eq, Debug)

///|
fn baseml_fail(message : String) -> Unit raise BasemlError {
  raise BasemlError(message)
}

///|
fn baseml_parse_int(text : String, label : String) -> Int raise BasemlError {
  codeml_parse_int(text, label) catch {
    CodemlError(_) => {
      baseml_fail(label + " must be an integer")
      0
    }
  }
}

///|
fn baseml_parse_double(
  text : String,
  label : String,
) -> Double raise BasemlError {
  codeml_parse_double(text, label) catch {
    CodemlError(_) => {
      baseml_fail(label + " must be a finite number")
      0.0
    }
  }
}

///|
fn baseml_copy_doubles(values : Array[Double]) -> Array[Double] {
  let copied : Array[Double] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn baseml_copy_matrix(values : Array[Array[Double]]) -> Array[Array[Double]] {
  let copied : Array[Array[Double]] = []
  for row in values {
    copied.push(baseml_copy_doubles(row))
  }
  copied
}

///|
fn baseml_copy_options(values : Array[BasemlOption]) -> Array[BasemlOption] {
  let copied : Array[BasemlOption] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn baseml_copy_branches(
  values : Array[BasemlBranchParameter],
) -> Array[BasemlBranchParameter] {
  let copied : Array[BasemlBranchParameter] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn baseml_copy_base_frequencies(
  value : BasemlBaseFrequencies?,
) -> BasemlBaseFrequencies? {
  match value {
    Some(frequencies) =>
      Some(BasemlBaseFrequencies::{
        thymine: frequencies.thymine,
        cytosine: frequencies.cytosine,
        adenine: frequencies.adenine,
        guanine: frequencies.guanine,
      })
    None => None
  }
}

///|
fn baseml_copy_q_matrix(value : BasemlQMatrix?) -> BasemlQMatrix? {
  match value {
    Some(matrix) =>
      Some(BasemlQMatrix::{
        average_ts_tv: matrix.average_ts_tv,
        rows: baseml_copy_matrix(matrix.rows),
      })
    None => None
  }
}

///|
fn baseml_copy_nodes(
  values : Array[BasemlNodeFrequencies],
) -> Array[BasemlNodeFrequencies] {
  let copied : Array[BasemlNodeFrequencies] = []
  for value in values {
    copied.push(BasemlNodeFrequencies::{
      node: value.node,
      root: value.root,
      frequency_parameters: baseml_copy_doubles(value.frequency_parameters),
      base_frequencies: baseml_copy_base_frequencies(value.base_frequencies),
    })
  }
  copied
}

///|
fn baseml_is_allowed_option(name : String) -> Bool {
  match name {
    "noisy"
    | "verbose"
    | "runmode"
    | "model"
    | "Mgene"
    | "ndata"
    | "clock"
    | "fix_kappa"
    | "kappa"
    | "fix_alpha"
    | "alpha"
    | "Malpha"
    | "ncatG"
    | "fix_rho"
    | "rho"
    | "nparK"
    | "nhomo"
    | "getSE"
    | "RateAncestor"
    | "Small_Diff"
    | "cleandata"
    | "icode"
    | "fix_blength"
    | "method" => true
    _ => false
  }
}

///|
fn baseml_is_double_option(name : String) -> Bool {
  match name {
    "kappa" | "alpha" | "rho" | "Small_Diff" => true
    _ => false
  }
}

///|
fn baseml_validate_bounded_option(
  name : String,
  value : Int,
) -> Unit raise BasemlError {
  match name {
    "model" =>
      if value < 0 || value > 10 {
        baseml_fail("BASEML model must be between 0 and 10")
      }
    "nparK" =>
      if value < 0 || value > 4 {
        baseml_fail("BASEML nparK must be between 0 and 4")
      }
    "nhomo" =>
      if value < 0 || value > 4 {
        baseml_fail("BASEML nhomo must be between 0 and 4")
      }
    "ncatG" => if value <= 0 { baseml_fail("BASEML ncatG must be positive") }
    "ndata" => if value <= 0 { baseml_fail("BASEML ndata must be positive") }
    "fix_kappa"
    | "fix_alpha"
    | "fix_rho"
    | "Malpha"
    | "getSE"
    | "cleandata"
    | "method" =>
      if value < 0 || value > 1 {
        baseml_fail("BASEML " + name + " must be 0 or 1")
      }
    _ => ()
  }
}

///|
fn baseml_normalize_model(value : String) -> String raise BasemlError {
  let words = codeml_words(value)
  if words.length() == 0 {
    baseml_fail("BASEML model must not be empty")
  }
  let model = baseml_parse_int(words[0], "BASEML model")
  baseml_validate_bounded_option("model", model)
  if words.length() > 1 && model != 9 && model != 10 {
    baseml_fail("BASEML model options are only valid for model 9 or 10")
  }
  let mut normalized = model.to_string()
  if words.length() > 1 {
    let options : Array[String] = []
    for index in 1..<words.length() {
      options.push(words[index])
    }
    normalized = normalized + " " + options.join(" ")
  }
  normalized
}

///|
fn baseml_normalize_option(
  name : String,
  value : String,
) -> String raise BasemlError {
  if name == "model" {
    return baseml_normalize_model(value)
  }
  if baseml_is_double_option(name) {
    return baseml_parse_double(value, "BASEML " + name).to_string()
  }
  let parsed = baseml_parse_int(value, "BASEML " + name)
  baseml_validate_bounded_option(name, parsed)
  parsed.to_string()
}

///|
/// Parse and strictly validate BASEML control text.
pub fn baseml_parse_control(text : String) -> BasemlControl raise BasemlError {
  let mut sequence_file = ""
  let mut output_file = ""
  let mut tree_file = ""
  let options : Array[BasemlOption] = []
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
      baseml_fail(
        "Malformed BASEML control line " + (line_index + 1).to_string(),
      )
    }
    let name = codeml_trim(uncommented[0:equals].to_owned())
    let value = codeml_trim(
      uncommented[equals + 1:uncommented.length()].to_owned(),
    )
    if name.length() == 0 || value.length() == 0 {
      baseml_fail(
        "Malformed BASEML control line " + (line_index + 1).to_string(),
      )
    }
    for previous in seen {
      if previous == name {
        baseml_fail("Duplicate BASEML control key: " + name)
      }
    }
    seen.push(name)
    match name {
      "seqfile" => sequence_file = value
      "outfile" => output_file = value
      "treefile" => tree_file = value
      _ => {
        if !baseml_is_allowed_option(name) {
          baseml_fail("Unknown BASEML option: " + name)
        }
        options.push(BasemlOption::{
          name,
          value: baseml_normalize_option(name, value),
        })
      }
    }
  }
  if sequence_file.length() == 0 {
    baseml_fail("BASEML control is missing seqfile")
  }
  if output_file.length() == 0 {
    baseml_fail("BASEML control is missing outfile")
  }
  if tree_file.length() == 0 {
    baseml_fail("BASEML control is missing treefile")
  }
  BasemlControl::{ sequence_file, output_file, tree_file, options }
}

///|
/// Write canonical BASEML control text.
pub fn baseml_write_control(control : BasemlControl) -> String {
  let output : Array[String] = [
    "seqfile = " + control.sequence_file,
    "outfile = " + control.output_file,
    "treefile = " + control.tree_file,
  ]
  for option in control.options {
    output.push(option.name + " = " + option.value)
  }
  output.join("\n") + "\n"
}

///|
pub fn BasemlOption::name(self : BasemlOption) -> String {
  self.name
}

///|
pub fn BasemlOption::value(self : BasemlOption) -> String {
  self.value
}

///|
pub fn BasemlControl::sequence_file(self : BasemlControl) -> String {
  self.sequence_file
}

///|
pub fn BasemlControl::output_file(self : BasemlControl) -> String {
  self.output_file
}

///|
pub fn BasemlControl::tree_file(self : BasemlControl) -> String {
  self.tree_file
}

///|
pub fn BasemlControl::options(self : BasemlControl) -> Array[BasemlOption] {
  baseml_copy_options(self.options)
}

///|
pub fn BasemlControl::option(self : BasemlControl, name : String) -> String? {
  for option in self.options {
    if option.name == name {
      return Some(option.value)
    }
  }
  None
}

///|
pub fn BasemlControl::model_number(self : BasemlControl) -> Int? {
  match self.option("model") {
    Some(value) => {
      let words = codeml_words(value)
      if words.length() == 0 {
        None
      } else {
        Some(baseml_parse_int(words[0], "BASEML model") catch { _ => 0 })
      }
    }
    None => None
  }
}

///|
pub fn BasemlControl::model_options(self : BasemlControl) -> String {
  match self.option("model") {
    Some(value) => {
      let words = codeml_words(value)
      let options : Array[String] = []
      for index in 1..<words.length() {
        options.push(words[index])
      }
      options.join(" ")
    }
    None => ""
  }
}

///|
/// Return PAML's name for a BASEML model number.
pub fn baseml_model_name(model : Int) -> String raise BasemlError {
  match model {
    0 => "JC69"
    1 => "K80"
    2 => "F81"
    3 => "F84"
    4 => "HKY85"
    5 => "T92"
    6 => "TN93"
    7 => "REV"
    8 => "UNREST"
    9 => "REVu"
    10 => "UNRESTu"
    _ => {
      baseml_fail("BASEML model must be between 0 and 10")
      ""
    }
  }
}

///|
fn baseml_frequencies(values : Array[Double]) -> BasemlBaseFrequencies {
  BasemlBaseFrequencies::{
    thymine: values[0],
    cytosine: values[1],
    adenine: values[2],
    guanine: values[3],
  }
}

///|
fn baseml_after_delimiter(line : String, delimiter : Int) -> String {
  let index = codeml_index_of_char(line, delimiter, 0)
  if index < 0 {
    ""
  } else {
    line[index + 1:line.length()].to_owned()
  }
}

///|
fn baseml_header(line : String) -> (String, String)? {
  let marker = "BASEML (in paml version "
  let start = codeml_index_of(line, marker, 0)
  if start < 0 {
    return None
  }
  let version_start = start + marker.length()
  let comma = codeml_index_of_char(line, ','.to_int(), version_start)
  if comma < 0 {
    return None
  }
  let version = codeml_trim(line[version_start:comma].to_owned())
  let close = codeml_index_of_char(line, ')'.to_int(), comma)
  if close < 0 || version.length() == 0 {
    return None
  }
  let remainder = codeml_trim(line[close + 1:line.length()].to_owned())
  let words = codeml_words(remainder)
  let description : Array[String] = []
  for index in 1..<words.length() {
    description.push(words[index])
  }
  Some((version, description.join(" ")))
}

///|
fn baseml_branch_from_line(line : String) -> BasemlBranchParameter? {
  let words = codeml_words(line)
  if words.length() == 0 || !codeml_contains(words[0], "..") {
    return None
  }
  let start = codeml_index_of(line, words[0], 0)
  let values = codeml_scan_numbers(
    line[start + words[0].length():line.length()].to_owned(),
  )
  if values.length() < 4 {
    return None
  }
  Some(BasemlBranchParameter::{
    branch: words[0],
    time: values[0],
    kappa: values[1],
    transitions: values[2],
    transversions: values[3],
  })
}

///|
fn baseml_node_from_line(line : String) -> BasemlNodeFrequencies? {
  let trimmed = codeml_trim(line)
  if !codeml_starts_with(trimmed, "Node #") {
    return None
  }
  let open = codeml_index_of_char(trimmed, '('.to_int(), 0)
  if open < 0 {
    return None
  }
  let node_text = codeml_trim(trimmed[6:open].to_owned())
  let node = baseml_parse_int(node_text, "BASEML node number") catch {
    _ => return None
  }
  let values = codeml_scan_numbers(
    trimmed[open + 1:trimmed.length()].to_owned(),
  )
  if values.length() != 4 && values.length() != 8 {
    return None
  }
  Some(BasemlNodeFrequencies::{
    node,
    root: false,
    frequency_parameters: values[0:4].to_owned(),
    base_frequencies: if values.length() == 8 {
      Some(baseml_frequencies(values[4:8].to_owned()))
    } else {
      None
    },
  })
}

///|
fn baseml_set_root(
  nodes : Array[BasemlNodeFrequencies],
  root : Int,
) -> Array[BasemlNodeFrequencies] {
  let updated : Array[BasemlNodeFrequencies] = []
  for node in nodes {
    updated.push(BasemlNodeFrequencies::{
      node: node.node,
      root: node.node == root,
      frequency_parameters: baseml_copy_doubles(node.frequency_parameters),
      base_frequencies: baseml_copy_base_frequencies(node.base_frequencies),
    })
  }
  updated
}

///|
/// Parse BASEML 4.x result text.
pub fn baseml_parse_results(text : String) -> BasemlResults raise BasemlError {
  if codeml_trim(text).length() == 0 {
    baseml_fail("BASEML results are empty")
  }
  let lines = codeml_lines(text)
  let mut version = ""
  let mut model_description = ""
  let mut ln_l_max : Double? = None
  let mut ln_likelihood : Double? = None
  let mut parameter_count = -1
  let mut parameter_list : Array[Double] = []
  let mut standard_errors : Array[Double] = []
  let mut tree_length : Double? = None
  let mut tree = ""
  let mut kappas : Array[Double] = []
  let branch_parameters : Array[BasemlBranchParameter] = []
  let mut rate_parameters : Array[Double] = []
  let mut rates : Array[Double] = []
  let mut rate_frequencies : Array[Double] = []
  let mut q_seen = false
  let mut q_average : Double? = None
  let q_rows : Array[Array[Double]] = []
  let mut alpha : Double? = None
  let mut rho : Double? = None
  let transition_probabilities : Array[Array[Double]] = []
  let mut base_frequencies : BasemlBaseFrequencies? = None
  let mut nodes : Array[BasemlNodeFrequencies] = []
  let mut seek_parameter_list = false
  let mut expect_standard_errors = false
  let mut expect_kappa = false
  let mut expect_base_frequencies = false
  let mut expect_q_rows = false
  let mut expect_transition_rows = false
  for line_index in 0..<lines.length() {
    let line = lines[line_index]
    let trimmed = codeml_trim(line)
    match baseml_header(line) {
      Some((parsed_version, description)) => {
        version = parsed_version
        model_description = description
        continue
      }
      None => ()
    }
    if codeml_contains(trimmed, "ln Lmax") {
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      let values = codeml_scan_numbers(
        if equals >= 0 {
          trimmed[equals + 1:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      if values.length() != 1 {
        baseml_fail(
          "BASEML ln Lmax line " +
          (line_index + 1).to_string() +
          " is malformed",
        )
      }
      ln_l_max = Some(values[0])
      continue
    }
    if codeml_contains(trimmed, "lnL(ntime:") {
      let np = codeml_index_of(trimmed, "np:", 0)
      let close = codeml_index_of_char(trimmed, ')'.to_int(), np)
      if np < 0 || close < 0 {
        baseml_fail(
          "BASEML lnL line " + (line_index + 1).to_string() + " is malformed",
        )
      }
      parameter_count = baseml_parse_int(
        trimmed[np + 3:close].to_owned(),
        "BASEML parameter count",
      )
      if parameter_count < 0 {
        baseml_fail("BASEML parameter count must not be negative")
      }
      let suffix = codeml_index_of(trimmed, "):", 0)
      let values = codeml_scan_numbers(
        if suffix >= 0 {
          trimmed[suffix + 2:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      if values.length() == 0 {
        baseml_fail(
          "BASEML lnL line " + (line_index + 1).to_string() + " is malformed",
        )
      }
      ln_likelihood = Some(values[0])
      seek_parameter_list = parameter_count > 0
      continue
    }
    if codeml_starts_with(trimmed, "SEs for parameters:") {
      expect_standard_errors = true
      continue
    }
    if expect_standard_errors && trimmed.length() > 0 {
      standard_errors = codeml_scan_numbers(trimmed)
      expect_standard_errors = false
      continue
    }
    if seek_parameter_list && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() == parameter_count {
        parameter_list = values
        seek_parameter_list = false
        continue
      }
    }
    if codeml_contains(trimmed, "tree length") {
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      let values = codeml_scan_numbers(
        if equals >= 0 {
          trimmed[equals + 1:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      if values.length() != 1 {
        baseml_fail(
          "BASEML tree length line " +
          (line_index + 1).to_string() +
          " is malformed",
        )
      }
      tree_length = Some(values[0])
      continue
    }
    if codeml_starts_with(trimmed, "(") && codeml_contains(trimmed, ":") {
      tree = trimmed
      continue
    }
    if codeml_contains(trimmed, "Parameters (kappa)") {
      expect_kappa = true
      continue
    }
    if codeml_contains(trimmed, "kappa under") {
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), 0)
      kappas = codeml_scan_numbers(
        if colon >= 0 {
          trimmed[colon + 1:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      continue
    }
    if expect_kappa && trimmed.length() > 0 {
      match baseml_branch_from_line(trimmed) {
        Some(branch) => {
          branch_parameters.push(branch)
          continue
        }
        None => {
          let values = codeml_scan_numbers(trimmed)
          if values.length() > 0 {
            kappas = values
            expect_kappa = false
            continue
          }
        }
      }
    }
    if codeml_contains(trimmed, "Rate parameters:") {
      rate_parameters = codeml_scan_numbers(
        baseml_after_delimiter(trimmed, ':'.to_int()),
      )
      continue
    }
    let lower = trimmed.to_lower()
    if codeml_contains(lower, "matrix q") {
      q_seen = true
      let equals = codeml_index_of_char(trimmed, '='.to_int(), 0)
      if equals >= 0 {
        let values = codeml_scan_numbers(
          trimmed[equals + 1:trimmed.length()].to_owned(),
        )
        if values.length() > 0 {
          q_average = Some(values[0])
        }
      }
      expect_q_rows = true
      continue
    }
    if expect_q_rows && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() == 4 {
        q_rows.push(values)
        if q_rows.length() == 4 {
          expect_q_rows = false
        }
        continue
      }
    }
    if codeml_starts_with(trimmed, "alpha") {
      let equals = codeml_last_index_of_char(trimmed, '='.to_int())
      let values = codeml_scan_numbers(
        if equals >= 0 {
          trimmed[equals + 1:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      if values.length() > 0 {
        alpha = Some(values[0])
      }
      continue
    }
    if codeml_contains(trimmed, "rho for") {
      let colon = codeml_index_of_char(trimmed, ':'.to_int(), 0)
      let values = codeml_scan_numbers(
        if colon >= 0 {
          trimmed[colon + 1:trimmed.length()].to_owned()
        } else {
          ""
        },
      )
      if values.length() > 0 {
        rho = Some(values[0])
      }
      continue
    }
    if codeml_starts_with(trimmed, "rate:") {
      rates = codeml_scan_numbers(baseml_after_delimiter(trimmed, ':'.to_int()))
      continue
    }
    if codeml_starts_with(trimmed, "freq:") {
      rate_frequencies = codeml_scan_numbers(
        baseml_after_delimiter(trimmed, ':'.to_int()),
      )
      continue
    }
    if codeml_contains(trimmed, "transition probabilities") {
      expect_transition_rows = true
      continue
    }
    if expect_transition_rows && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() > 0 {
        transition_probabilities.push(values)
        if rates.length() > 0 &&
          transition_probabilities.length() == rates.length() {
          expect_transition_rows = false
        }
        continue
      }
    }
    if codeml_contains(trimmed, "base frequency parameters") ||
      (
        codeml_starts_with(trimmed, "Base frequencies") &&
        codeml_scan_numbers(trimmed).length() == 0
      ) {
      expect_base_frequencies = true
      continue
    }
    if codeml_starts_with(trimmed, "Base frequencies:") {
      let values = codeml_scan_numbers(
        baseml_after_delimiter(trimmed, ':'.to_int()),
      )
      if values.length() == 4 {
        base_frequencies = Some(baseml_frequencies(values))
      }
      continue
    }
    if expect_base_frequencies && trimmed.length() > 0 {
      let values = codeml_scan_numbers(trimmed)
      if values.length() != 4 {
        baseml_fail(
          "BASEML base frequencies line " +
          (line_index + 1).to_string() +
          " is malformed",
        )
      }
      base_frequencies = Some(baseml_frequencies(values))
      expect_base_frequencies = false
      continue
    }
    match baseml_node_from_line(trimmed) {
      Some(node) => {
        for previous in nodes {
          if previous.node == node.node {
            baseml_fail(
              "Duplicate BASEML node frequency entry: " + node.node.to_string(),
            )
          }
        }
        nodes.push(node)
        continue
      }
      None => ()
    }
    if codeml_starts_with(trimmed, "Note: node ") &&
      codeml_contains(trimmed, " is root") {
      let values = codeml_scan_numbers(trimmed)
      if values.length() != 1 {
        baseml_fail(
          "BASEML root node line " +
          (line_index + 1).to_string() +
          " is malformed",
        )
      }
      nodes = baseml_set_root(nodes, values[0].to_int())
    }
  }
  if version.length() == 0 {
    baseml_fail("BASEML results do not contain a valid BASEML header")
  }
  let fitted_likelihood = match ln_likelihood {
    Some(value) => value
    None => {
      baseml_fail("BASEML results do not contain a likelihood estimate")
      0.0
    }
  }
  if parameter_count < 0 {
    baseml_fail("BASEML results do not contain a parameter count")
  }
  if parameter_list.length() != parameter_count {
    baseml_fail(
      "BASEML parameter vector length does not match the reported count",
    )
  }
  if standard_errors.length() > 0 && standard_errors.length() != parameter_count {
    baseml_fail(
      "BASEML standard-error vector length does not match the parameter count",
    )
  }
  let fitted_tree_length = match tree_length {
    Some(value) => value
    None => {
      baseml_fail("BASEML results do not contain a tree length")
      0.0
    }
  }
  if tree.length() == 0 {
    baseml_fail("BASEML results do not contain a branch-length tree")
  }
  if q_seen && q_rows.length() != 4 {
    baseml_fail("BASEML Q matrix must contain four rows")
  }
  if transition_probabilities.length() > 0 {
    if rates.length() == 0 ||
      transition_probabilities.length() != rates.length() {
      baseml_fail(
        "BASEML transition matrix dimensions do not match rate categories",
      )
    }
    for row in transition_probabilities {
      if row.length() != rates.length() {
        baseml_fail(
          "BASEML transition matrix dimensions do not match rate categories",
        )
      }
    }
  }
  BasemlResults::{
    version,
    model_description,
    ln_l_max,
    ln_likelihood: fitted_likelihood,
    parameter_count,
    parameter_list,
    standard_errors,
    tree_length: fitted_tree_length,
    tree,
    kappas,
    branch_parameters,
    rate_parameters,
    rates,
    rate_frequencies,
    q_matrix: if q_seen {
      Some(BasemlQMatrix::{ average_ts_tv: q_average, rows: q_rows })
    } else {
      None
    },
    alpha,
    rho,
    transition_probabilities,
    base_frequencies,
    nodes,
  }
}

///|
pub fn BasemlBaseFrequencies::thymine(self : BasemlBaseFrequencies) -> Double {
  self.thymine
}

///|
pub fn BasemlBaseFrequencies::cytosine(self : BasemlBaseFrequencies) -> Double {
  self.cytosine
}

///|
pub fn BasemlBaseFrequencies::adenine(self : BasemlBaseFrequencies) -> Double {
  self.adenine
}

///|
pub fn BasemlBaseFrequencies::guanine(self : BasemlBaseFrequencies) -> Double {
  self.guanine
}

///|
pub fn BasemlQMatrix::average_ts_tv(self : BasemlQMatrix) -> Double? {
  self.average_ts_tv
}

///|
pub fn BasemlQMatrix::rows(self : BasemlQMatrix) -> Array[Array[Double]] {
  baseml_copy_matrix(self.rows)
}

///|
pub fn BasemlBranchParameter::branch(self : BasemlBranchParameter) -> String {
  self.branch
}

///|
pub fn BasemlBranchParameter::time(self : BasemlBranchParameter) -> Double {
  self.time
}

///|
pub fn BasemlBranchParameter::kappa(self : BasemlBranchParameter) -> Double {
  self.kappa
}

///|
pub fn BasemlBranchParameter::transitions(
  self : BasemlBranchParameter,
) -> Double {
  self.transitions
}

///|
pub fn BasemlBranchParameter::transversions(
  self : BasemlBranchParameter,
) -> Double {
  self.transversions
}

///|
pub fn BasemlNodeFrequencies::node(self : BasemlNodeFrequencies) -> Int {
  self.node
}

///|
pub fn BasemlNodeFrequencies::is_root(self : BasemlNodeFrequencies) -> Bool {
  self.root
}

///|
pub fn BasemlNodeFrequencies::frequency_parameters(
  self : BasemlNodeFrequencies,
) -> Array[Double] {
  baseml_copy_doubles(self.frequency_parameters)
}

///|
pub fn BasemlNodeFrequencies::base_frequencies(
  self : BasemlNodeFrequencies,
) -> BasemlBaseFrequencies? {
  baseml_copy_base_frequencies(self.base_frequencies)
}

///|
pub fn BasemlResults::version(self : BasemlResults) -> String {
  self.version
}

///|
pub fn BasemlResults::model_description(self : BasemlResults) -> String {
  self.model_description
}

///|
pub fn BasemlResults::ln_l_max(self : BasemlResults) -> Double? {
  self.ln_l_max
}

///|
pub fn BasemlResults::ln_likelihood(self : BasemlResults) -> Double {
  self.ln_likelihood
}

///|
pub fn BasemlResults::parameter_count(self : BasemlResults) -> Int {
  self.parameter_count
}

///|
pub fn BasemlResults::parameter_list(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.parameter_list)
}

///|
pub fn BasemlResults::standard_errors(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.standard_errors)
}

///|
pub fn BasemlResults::tree_length(self : BasemlResults) -> Double {
  self.tree_length
}

///|
pub fn BasemlResults::tree(self : BasemlResults) -> String {
  self.tree
}

///|
pub fn BasemlResults::kappas(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.kappas)
}

///|
pub fn BasemlResults::branch_parameters(
  self : BasemlResults,
) -> Array[BasemlBranchParameter] {
  baseml_copy_branches(self.branch_parameters)
}

///|
pub fn BasemlResults::rate_parameters(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.rate_parameters)
}

///|
pub fn BasemlResults::rates(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.rates)
}

///|
pub fn BasemlResults::rate_frequencies(self : BasemlResults) -> Array[Double] {
  baseml_copy_doubles(self.rate_frequencies)
}

///|
pub fn BasemlResults::q_matrix(self : BasemlResults) -> BasemlQMatrix? {
  baseml_copy_q_matrix(self.q_matrix)
}

///|
pub fn BasemlResults::alpha(self : BasemlResults) -> Double? {
  self.alpha
}

///|
pub fn BasemlResults::rho(self : BasemlResults) -> Double? {
  self.rho
}

///|
pub fn BasemlResults::transition_probabilities(
  self : BasemlResults,
) -> Array[Array[Double]] {
  baseml_copy_matrix(self.transition_probabilities)
}

///|
pub fn BasemlResults::base_frequencies(
  self : BasemlResults,
) -> BasemlBaseFrequencies? {
  baseml_copy_base_frequencies(self.base_frequencies)
}

///|
pub fn BasemlResults::nodes(
  self : BasemlResults,
) -> Array[BasemlNodeFrequencies] {
  baseml_copy_nodes(self.nodes)
}

///|
/// Akaike information criterion.
pub fn BasemlResults::aic(self : BasemlResults) -> Double {
  2.0 * self.parameter_count.to_double() - 2.0 * self.ln_likelihood
}

///|
/// Bayesian information criterion.
pub fn BasemlResults::bic(
  self : BasemlResults,
  observations : Int,
) -> Double raise BasemlError {
  if observations <= 0 {
    baseml_fail("BASEML BIC observations must be positive")
  }
  self.parameter_count.to_double() * @math.ln(observations.to_double()) -
  2.0 * self.ln_likelihood
}

///|
/// Compare two nested BASEML models with a chi-square likelihood-ratio test.
pub fn baseml_likelihood_ratio(
  null_model : BasemlResults,
  alternative_model : BasemlResults,
  alpha? : Double = 0.05,
) -> BasemlLikelihoodRatio raise BasemlError {
  if alpha <= 0.0 || alpha >= 1.0 || alpha.is_nan() {
    baseml_fail("BASEML LRT alpha must be between 0 and 1")
  }
  let degrees = alternative_model.parameter_count - null_model.parameter_count
  if degrees <= 0 {
    baseml_fail(
      "BASEML alternative model must have more parameters than the null model",
    )
  }
  if alternative_model.ln_likelihood + 1.0e-10 < null_model.ln_likelihood {
    baseml_fail(
      "BASEML alternative model has a lower likelihood than the null model",
    )
  }
  let statistic = 2.0 *
    (alternative_model.ln_likelihood - null_model.ln_likelihood)
  let p_value = codeml_gamma_q(degrees.to_double() / 2.0, statistic / 2.0)
  BasemlLikelihoodRatio::{
    statistic,
    degrees_of_freedom: degrees,
    p_value,
    alpha,
    significant: p_value < alpha,
  }
}

///|
pub fn BasemlLikelihoodRatio::statistic(self : BasemlLikelihoodRatio) -> Double {
  self.statistic
}

///|
pub fn BasemlLikelihoodRatio::degrees_of_freedom(
  self : BasemlLikelihoodRatio,
) -> Int {
  self.degrees_of_freedom
}

///|
pub fn BasemlLikelihoodRatio::p_value(self : BasemlLikelihoodRatio) -> Double {
  self.p_value
}

///|
pub fn BasemlLikelihoodRatio::alpha(self : BasemlLikelihoodRatio) -> Double {
  self.alpha
}

///|
pub fn BasemlLikelihoodRatio::significant(self : BasemlLikelihoodRatio) -> Bool {
  self.significant
}

///|
/// Compact official-style control fixture used by tests and examples.
pub fn baseml_example_control_text() -> String {
  "seqfile = alignment.phylip\n" +
  "treefile = species.tree\n" +
  "outfile = baseml.out\n" +
  "noisy = 0\n" +
  "verbose = 0\n" +
  "runmode = 0\n" +
  "model = 7\n" +
  "Mgene = 1\n" +
  "clock = 0\n" +
  "fix_kappa = 0\n" +
  "kappa = 5\n" +
  "fix_alpha = 0\n" +
  "alpha = 0.5\n" +
  "Malpha = 0\n" +
  "ncatG = 5\n" +
  "fix_rho = 1\n" +
  "rho = 0\n" +
  "nparK = 0\n" +
  "nhomo = 0\n" +
  "getSE = 1\n" +
  "RateAncestor = 0\n" +
  "Small_Diff = 7e-6\n" +
  "cleandata = 1\n" +
  "method = 0\n"
}

///|
/// Compact result fixture adapted from Biopython's PAML 4.7 BASEML data.
pub fn baseml_example_results_text() -> String {
  "BASEML (in paml version 4.7, January 2013) alignment.phylip REV dGamma (ncatG=5)\n" +
  "ln Lmax (unconstrained) = -316.049385\n" +
  "lnL(ntime: 7 np: 13): -319.774741 +0.000000\n" +
  "0.000004 0.000004 0.000004 0.004533 0.004529 0.000004 0.009174 998.999980 130.949080 0.000010 0.000010 0.000010 200.957530\n" +
  "SEs for parameters:\n" +
  "0.001 0.001 0.001 0.002 0.002 0.001 0.003 2.0 1.0 0.1 0.1 0.1 4.0\n" +
  "tree length = 0.01826\n" +
  "(((Homo_sapie:0.000004,Pan_troglo:0.004533):0.000004,Gorilla_go:0.004529):0.000004,Pongo_pygm:0.000004,Macaca_mul:0.009174);\n" +
  "Rate parameters: 998.99998 130.94908 0.00001 0.00001 0.00001\n" +
  "Base frequencies: 0.20090 0.16306 0.37027 0.26577\n" +
  "Rate matrix Q, Average Ts/Tv = 3.3698\n" +
  "-2.483179 1.865730 0.617449 0.000000\n" +
  "2.298662 -2.298662 0.000000 0.000000\n" +
  "0.335015 0.000000 -0.338059 0.003044\n" +
  "0.000000 0.000000 0.004241 -0.004241\n" +
  "alpha (gamma, K=5) = 200.95753\n" +
  "rate: 0.90327 0.96137 0.99838 1.03635 1.10064\n" +
  "freq: 0.20000 0.20000 0.20000 0.20000 0.20000\n"
}
