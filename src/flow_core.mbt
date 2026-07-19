///|
/// flowCore - Flow Cytometry FCS file processing and analysis.
/// Inspired by Bioconductor flowCore package.
///
/// Provides FCS file parsing, data transformation, compensation,
/// and gating (rectangle, polygon, ellipsoid, quadrant) for flow cytometry data.

///|
/// Square root via Newton's method.
pub fn flow_core_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 1.0
  }
  let mut guess = x / 2.0
  if guess < 1.0 {
    guess = 1.0
  }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = next - guess
    if diff < 0.0 {
      if -diff < 0.00000000000001 {
        return next
      }
    } else if diff < 0.00000000000001 {
      return next
    }
    guess = next
    i = i + 1
  }
  guess
}

///|
/// Natural logarithm via series expansion.
pub fn flow_core_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -10000000000.0
  }
  if x == 1.0 {
    return 0.0
  }
  let mut v = x
  let mut offset = 0.0
  while v > 2.0 {
    v = v / 2.718281828459045
    offset = offset + 1.0
  }
  while v < 0.5 {
    v = v * 2.718281828459045
    offset = offset - 1.0
  }
  let y = (v - 1.0) / (v + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i <= 50 {
    term = term * y2
    let delta = term / (2 * i + 1).to_double()
    sum = sum + delta
    if delta < 0.0 {
      if -delta < 0.0000000000000001 {
        return 2.0 * sum + offset
      }
    } else if delta < 0.0000000000000001 {
      return 2.0 * sum + offset
    }
    i = i + 1
  }
  2.0 * sum + offset
}

///|
/// Exponential function via Taylor series.
pub fn flow_core_exp(x : Double) -> Double {
  if x == 0.0 {
    return 1.0
  }
  if x < 0.0 {
    return 1.0 / flow_core_exp(-x)
  }
  let mut k = 0
  let mut r = x
  while r >= 1.0 {
    r = r - 1.0
    k = k + 1
  }
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 30 {
    term = term * r / i.to_double()
    sum = sum + term
    if term < 0.0 {
      if -term < 0.0000000000000001 {
        break
      }
    } else if term < 0.0000000000000001 {
      break
    }
    i = i + 1
  }
  let mut e_pow_k = 1.0
  let mut j = 0
  while j < k {
    e_pow_k = e_pow_k * 2.718281828459045
    j = j + 1
  }
  e_pow_k * sum
}

///|
/// Base-10 logarithm.
pub fn flow_core_log10(x : Double) -> Double {
  flow_core_ln(x) / flow_core_ln(10.0)
}

///|
/// Inverse hyperbolic sine (asinh).
pub fn flow_core_asinh(x : Double) -> Double {
  flow_core_ln(x + flow_core_sqrt(x * x + 1.0))
}

///|
/// Gate type enumeration.
pub enum GateType {
  RectangleGate
  PolygonGate
  EllipsoidGate
  QuadrantGate
  BooleanGate
} derive(Eq)

///|
/// Transformation type enumeration.
pub enum TransformType {
  LogTransform
  AsinhTransform
  LogicleTransform
  LinearTransform
} derive(Eq)

///|
/// String representation of GateType.
pub fn GateType::to_string(self : GateType) -> String {
  match self {
    RectangleGate => "RectangleGate"
    PolygonGate => "PolygonGate"
    EllipsoidGate => "EllipsoidGate"
    QuadrantGate => "QuadrantGate"
    BooleanGate => "BooleanGate"
  }
}

///|
/// String representation of TransformType.
pub fn TransformType::to_string(self : TransformType) -> String {
  match self {
    LogTransform => "LogTransform"
    AsinhTransform => "AsinhTransform"
    LogicleTransform => "LogicleTransform"
    LinearTransform => "LinearTransform"
  }
}

///|
/// FCS file representation with text, data, and analysis segments.
pub struct FCSFile {
  text_segment : Map[String, String]
  data_segment : Array[Array[Double]]
  analysis_segment : String
  n_events : Int
  n_parameters : Int
  parameter_names : Array[String]
  parameter_ranges : Array[(Double, Double)]
}

///|
/// Flow frame containing expression matrix and metadata.
pub struct FlowFrame {
  exprs : Array[Array[Double]]
  params : Array[ParameterDescription]
  description : Map[String, String]
  n_events : Int
  n_params : Int
}

///|
/// Parameter description for a flow frame column.
pub struct ParameterDescription {
  name : String
  range : Double
  min_range : Double
  max_range : Double
  compensation : Array[Double]
}

///|
/// Flow set containing multiple flow frames.
pub struct FlowSet {
  frames : Array[FlowFrame]
  sample_names : Array[String]
  pheno_data : Map[String, String]
}

///|
/// Gate for filtering flow cytometry events.
pub struct Gate {
  gate_type : GateType
  parameters : Array[String]
  coordinates : Array[(Double, Double)]
  filter_id : String
}

///|
/// Transformation specification for flow data.
pub struct Transformation {
  transform_type : TransformType
  parameters : Array[Double]
}

///|
/// Compensation matrix for fluorescence spillover correction.
pub struct CompensationMatrix {
  matrix : Array[Array[Double]]
  fluorochromes : Array[String]
  spillover : Array[Double]
}

///|
/// Result of gating operation.
pub struct GatingResult {
  filtered_events : Array[Array[Double]]
  n_filtered : Int
  n_total : Int
  percentages : Array[Double]
}

///|
/// Construct a new ParameterDescription.
pub fn ParameterDescription::new(
  name : String,
  range : Double,
  min_range : Double,
  max_range : Double,
) -> ParameterDescription {
  ParameterDescription::{
    name,
    range,
    min_range,
    max_range,
    compensation: Array::new(),
  }
}

///|
/// String representation of ParameterDescription.
pub fn ParameterDescription::to_string(self : ParameterDescription) -> String {
  let mut result = "Parameter: " + self.name + "\n"
  result = result + "  range: " + self.range.to_string() + "\n"
  result = result + "  min_range: " + self.min_range.to_string() + "\n"
  result = result + "  max_range: " + self.max_range.to_string()
  result
}

///|
/// Construct a new FCSFile with defaults.
pub fn FCSFile::new() -> FCSFile {
  FCSFile::{
    text_segment: Map([], capacity=32),
    data_segment: Array::new(),
    analysis_segment: "",
    n_events: 0,
    n_parameters: 0,
    parameter_names: Array::new(),
    parameter_ranges: Array::new(),
  }
}

///|
/// String representation of FCSFile.
pub fn FCSFile::to_string(self : FCSFile) -> String {
  let mut result = "FCSFile:\n"
  result = result + "  Events: " + self.n_events.to_string() + "\n"
  result = result + "  Parameters: " + self.n_parameters.to_string() + "\n"
  result = result + "  Parameter names:"
  let mut i = 0
  while i < self.parameter_names.length() {
    result = result + " " + self.parameter_names[i]
    i = i + 1
  }
  result = result +
    "\n  Text segment keys: " +
    self.text_segment.length().to_string()
  result
}

///|
/// Construct a new FlowFrame from expression matrix and parameters.
pub fn FlowFrame::new(
  exprs : Array[Array[Double]],
  params : Array[ParameterDescription],
) -> FlowFrame {
  let n_events = exprs.length()
  let n_params = params.length()
  FlowFrame::{
    exprs,
    params,
    description: Map([], capacity=16),
    n_events,
    n_params,
  }
}

///|
/// Create a FlowFrame from an FCSFile.
pub fn FlowFrame::from_fcs(fcs : FCSFile) -> FlowFrame {
  let params : Array[ParameterDescription] = Array::new()
  let mut i = 0
  while i < fcs.parameter_names.length() {
    let name = fcs.parameter_names[i]
    let (min_val, max_val) = if i < fcs.parameter_ranges.length() {
      fcs.parameter_ranges[i]
    } else {
      (0.0, 0.0)
    }
    params.push(ParameterDescription::{
      name,
      range: max_val - min_val,
      min_range: min_val,
      max_range: max_val,
      compensation: Array::new(),
    })
    i = i + 1
  }
  let n_events = fcs.data_segment.length()
  let n_params = fcs.n_parameters
  FlowFrame::{
    exprs: fcs.data_segment,
    params,
    description: fcs.text_segment,
    n_events,
    n_params,
  }
}

///|
/// String representation of FlowFrame.
pub fn FlowFrame::to_string(self : FlowFrame) -> String {
  let mut result = "FlowFrame:\n"
  result = result + "  Events: " + self.n_events.to_string() + "\n"
  result = result + "  Parameters: " + self.n_params.to_string() + "\n"
  result = result + "  Parameters:"
  let mut i = 0
  while i < self.params.length() {
    result = result + " " + self.params[i].name
    i = i + 1
  }
  result = result +
    "\n  Description keys: " +
    self.description.length().to_string()
  result
}

///|
/// Construct a new FlowSet.
pub fn FlowSet::new(
  frames : Array[FlowFrame],
  sample_names : Array[String],
) -> FlowSet {
  FlowSet::{ frames, sample_names, pheno_data: Map([], capacity=8) }
}

///|
/// String representation of FlowSet.
pub fn FlowSet::to_string(self : FlowSet) -> String {
  let mut result = "FlowSet:\n"
  result = result + "  Frames: " + self.frames.length().to_string() + "\n"
  result = result + "  Sample names:"
  let mut i = 0
  while i < self.sample_names.length() {
    result = result + " " + self.sample_names[i]
    i = i + 1
  }
  result = result +
    "\n  Pheno data keys: " +
    self.pheno_data.length().to_string()
  result
}

///|
/// Construct a new Gate.
pub fn Gate::new(
  gate_type : GateType,
  parameters : Array[String],
  coordinates : Array[(Double, Double)],
  filter_id : String,
) -> Gate {
  Gate::{ gate_type, parameters, coordinates, filter_id }
}

///|
/// String representation of Gate.
pub fn Gate::to_string(self : Gate) -> String {
  let mut result = "Gate:\n"
  result = result + "  Type: " + self.gate_type.to_string() + "\n"
  result = result + "  Filter ID: " + self.filter_id + "\n"
  result = result + "  Parameters:"
  let mut i = 0
  while i < self.parameters.length() {
    result = result + " " + self.parameters[i]
    i = i + 1
  }
  result = result + "\n  Coordinates: " + self.coordinates.length().to_string()
  result
}

///|
/// Construct a new Transformation.
pub fn Transformation::new(
  transform_type : TransformType,
  parameters : Array[Double],
) -> Transformation {
  Transformation::{ transform_type, parameters }
}

///|
/// String representation of Transformation.
pub fn Transformation::to_string(self : Transformation) -> String {
  let mut result = "Transformation:\n"
  result = result + "  Type: " + self.transform_type.to_string() + "\n"
  result = result + "  Parameters:"
  let mut i = 0
  while i < self.parameters.length() {
    result = result + " " + self.parameters[i].to_string()
    i = i + 1
  }
  result
}

///|
/// Construct a new CompensationMatrix.
pub fn CompensationMatrix::new(
  matrix : Array[Array[Double]],
  fluorochromes : Array[String],
) -> CompensationMatrix {
  let spillover : Array[Double] = Array::new()
  let mut i = 0
  while i < matrix.length() {
    let mut j = 0
    while j < matrix[i].length() {
      spillover.push(matrix[i][j])
      j = j + 1
    }
    i = i + 1
  }
  CompensationMatrix::{ matrix, fluorochromes, spillover }
}

///|
/// String representation of CompensationMatrix.
pub fn CompensationMatrix::to_string(self : CompensationMatrix) -> String {
  let mut result = "CompensationMatrix:\n"
  result = result + "  Fluorochromes:"
  let mut i = 0
  while i < self.fluorochromes.length() {
    result = result + " " + self.fluorochromes[i]
    i = i + 1
  }
  result = result + "\n  Matrix (" + self.matrix.length().to_string() + "x):"
  i = 0
  while i < self.matrix.length() {
    result = result + "\n    ["
    let mut j = 0
    while j < self.matrix[i].length() {
      if j > 0 {
        result = result + ", "
      }
      result = result + self.matrix[i][j].to_string()
      j = j + 1
    }
    result = result + "]"
    i = i + 1
  }
  result
}

///|
/// Construct a new GatingResult.
pub fn GatingResult::new(
  filtered_events : Array[Array[Double]],
  n_filtered : Int,
  n_total : Int,
  percentages : Array[Double],
) -> GatingResult {
  GatingResult::{ filtered_events, n_filtered, n_total, percentages }
}

///|
/// String representation of GatingResult.
pub fn GatingResult::to_string(self : GatingResult) -> String {
  let mut result = "GatingResult:\n"
  result = result +
    "  Filtered: " +
    self.n_filtered.to_string() +
    "/" +
    self.n_total.to_string()
  if self.n_total > 0 {
    result = result +
      " (" +
      (self.n_filtered.to_double() * 100.0 / self.n_total.to_double()).to_string() +
      "%)"
  }
  result = result + "\n  Percentages:"
  let mut i = 0
  while i < self.percentages.length() {
    result = result + " " + self.percentages[i].to_string()
    i = i + 1
  }
  result
}

///|
/// Parse an integer from a string, returning 0 on failure.
fn flow_core_parse_int(s : String) -> Int {
  match parse_double(s) {
    Some(v) => v.to_int()
    None => 0
  }
}

///|
/// Trim whitespace from both ends of a string.
fn flow_core_trim(s : String) -> String {
  let n = s.length()
  if n == 0 {
    return ""
  }
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let nl : UInt16 = '\n'.to_int().to_uint16()
  let mut start = 0
  let mut end = n
  while start < end {
    let c = s.unsafe_get(start)
    if c == sp || c == tab || c == cr || c == nl {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == sp || c == tab || c == cr || c == nl {
      end = end - 1
    } else {
      break
    }
  }
  if start >= end {
    ""
  } else {
    s[start:end].to_owned()
  }
}

///|
/// Parse a line of numbers separated by spaces, tabs, or commas.
fn flow_core_parse_numbers(s : String) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let parts = split_by_char(s, ' '.to_int())
  let mut i = 0
  while i < parts.length() {
    let part = parts[i]
    if part.length() > 0 {
      let tab_parts = split_by_char(part, '\t'.to_int())
      let mut j = 0
      while j < tab_parts.length() {
        let tp = tab_parts[j]
        if tp.length() > 0 {
          let comma_parts = split_by_char(tp, ','.to_int())
          let mut k = 0
          while k < comma_parts.length() {
            let val = comma_parts[k]
            if val.length() > 0 {
              match parse_double(val) {
                Some(v) => result.push(v)
                None => ()
              }
            }
            k = k + 1
          }
        }
        j = j + 1
      }
    }
    i = i + 1
  }
  result
}

///|
/// Check if content starts with "FCS".
fn flow_core_is_fcs_header(content : String) -> Bool {
  content.length() >= 3 &&
  content.unsafe_get(0) == 'F'.to_int().to_uint16() &&
  content.unsafe_get(1) == 'C'.to_int().to_uint16() &&
  content.unsafe_get(2) == 'S'.to_int().to_uint16()
}

///|
/// Parse FCS text segment using delimiter (first character).
fn flow_core_parse_text_segment(
  text : String,
  map : Map[String, String],
) -> Unit {
  let n = text.length()
  if n < 2 {
    return
  }
  let delim = text.unsafe_get(0)
  let rest = text[1:n].to_owned()
  let parts = split_by_char(rest, delim.to_int())
  let mut i = 0
  while i + 1 < parts.length() {
    let key = parts[i]
    let value = parts[i + 1]
    if key.length() > 0 {
      map[key] = value
    }
    i = i + 2
  }
}

///|
/// Parse simplified text format: lines starting with $ are key/value, others are data.
fn flow_core_parse_simplified(
  content : String,
  text_map : Map[String, String],
  data : Array[Array[Double]],
) -> Unit {
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let raw = lines[i]
    let line = flow_core_trim(raw)
    if line.length() == 0 {
      i = i + 1
      continue
    }
    let dollar : UInt16 = '$'.to_int().to_uint16()
    if line.unsafe_get(0) == dollar {
      let sp_parts = split_by_char(line, ' '.to_int())
      if sp_parts.length() >= 2 {
        let key = sp_parts[0]
        let value = sp_parts[1]
        text_map[key] = value
      } else {
        let eq_parts = split_by_char(line, '='.to_int())
        if eq_parts.length() >= 2 {
          let key = flow_core_trim(eq_parts[0])
          let value = flow_core_trim(eq_parts[1])
          text_map[key] = value
        }
      }
    } else if line != "---DATA---" && line != "[DATA]" {
      let nums = flow_core_parse_numbers(line)
      if nums.length() > 0 {
        data.push(nums)
      }
    }
    i = i + 1
  }
}

///|
/// Parse FCS file content into FCSFile.
pub fn FCSFile::parse(content : String) -> FCSFile {
  let text_segment : Map[String, String] = Map([], capacity=32)
  let data_segment : Array[Array[Double]] = Array::new()
  let analysis_segment = ""
  if flow_core_is_fcs_header(content) && content.length() >= 58 {
    let text_start_str = flow_core_trim(content[10:18].to_owned())
    let text_end_str = flow_core_trim(content[18:26].to_owned())
    let data_start_str = flow_core_trim(content[26:34].to_owned())
    let data_end_str = flow_core_trim(content[34:42].to_owned())
    let text_start = flow_core_parse_int(text_start_str)
    let text_end = flow_core_parse_int(text_end_str)
    let data_start = flow_core_parse_int(data_start_str)
    let data_end = flow_core_parse_int(data_end_str)
    if text_start > 0 && text_end > text_start && text_end <= content.length() {
      let text = content[text_start:text_end + 1].to_owned()
      flow_core_parse_text_segment(text, text_segment)
    } else if content.length() >= 10 {
      let text = content[10:content.length()].to_owned()
      flow_core_parse_text_segment(text, text_segment)
    }
    if data_start > 0 && data_end > data_start && data_end <= content.length() {
      let data_str = content[data_start:data_end + 1].to_owned()
      let lines = split_by_char(data_str, '\n'.to_int())
      let mut i = 0
      while i < lines.length() {
        let nums = flow_core_parse_numbers(lines[i])
        if nums.length() > 0 {
          data_segment.push(nums)
        }
        i = i + 1
      }
    }
  } else {
    flow_core_parse_simplified(content, text_segment, data_segment)
  }
  let n_events = if text_segment.contains("$TOT") {
    flow_core_parse_int(text_segment["$TOT"])
  } else {
    data_segment.length()
  }
  let n_parameters = if text_segment.contains("$PAR") {
    flow_core_parse_int(text_segment["$PAR"])
  } else if data_segment.length() > 0 {
    data_segment[0].length()
  } else {
    0
  }
  let parameter_names : Array[String] = Array::new()
  let parameter_ranges : Array[(Double, Double)] = Array::new()
  let mut i = 1
  while i <= n_parameters {
    let name_key = "$P" + i.to_string() + "N"
    let range_key = "$P" + i.to_string() + "R"
    let name = if text_segment.contains(name_key) {
      text_segment[name_key]
    } else {
      "P" + i.to_string()
    }
    parameter_names.push(name)
    let range_val = if text_segment.contains(range_key) {
      flow_core_parse_int(text_segment[range_key]).to_double()
    } else {
      0.0
    }
    parameter_ranges.push((0.0, range_val))
    i = i + 1
  }
  FCSFile::{
    text_segment,
    data_segment,
    analysis_segment,
    n_events,
    n_parameters,
    parameter_names,
    parameter_ranges,
  }
}

///|
/// Read an FCS file. Note: treats file_path as FCS content since this library does not perform file I/O.
pub fn flow_core_read_fcs(file_path : String) -> FlowFrame {
  let fcs = FCSFile::parse(file_path)
  FlowFrame::from_fcs(fcs)
}

///|
/// Write a FlowFrame to FCS format string.
pub fn flow_core_write_fcs(frame : FlowFrame, file_path : String) -> String {
  let _ = file_path
  let buf = StringBuilder::new()
  buf.write_string("FCS3.1\n")
  buf.write_string("$FCSVERSION\t3.1\n")
  buf.write_string("$TOT\t" + frame.n_events.to_string() + "\n")
  buf.write_string("$PAR\t" + frame.n_params.to_string() + "\n")
  let mut i = 0
  while i < frame.params.length() {
    let p = frame.params[i]
    buf.write_string("$P" + (i + 1).to_string() + "N\t" + p.name + "\n")
    buf.write_string(
      "$P" + (i + 1).to_string() + "R\t" + p.range.to_string() + "\n",
    )
    i = i + 1
  }
  let desc_keys = frame.description.keys().collect()
  i = 0
  while i < desc_keys.length() {
    let key = desc_keys[i]
    buf.write_string(key + "\t" + frame.description[key] + "\n")
    i = i + 1
  }
  buf.write_string("---DATA---\n")
  i = 0
  while i < frame.exprs.length() {
    let row = frame.exprs[i]
    let mut j = 0
    while j < row.length() {
      if j > 0 {
        buf.write_string("\t")
      }
      buf.write_string(row[j].to_string())
      j = j + 1
    }
    buf.write_string("\n")
    i = i + 1
  }
  buf.to_string()
}

///|
/// Get marker (parameter) names from a FlowFrame.
pub fn flow_core_get_markers(frame : FlowFrame) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < frame.params.length() {
    result.push(frame.params[i].name)
    i = i + 1
  }
  result
}

///|
/// Log transform: a * ln(x) + c (b is base, used if > 1).
pub fn flow_core_log_transform(
  x : Double,
  a : Double,
  b : Double,
  c : Double,
) -> Double {
  if x <= 0.0 {
    return c
  }
  if b > 1.0 {
    a * flow_core_ln(x) / flow_core_ln(b) + c
  } else {
    a * flow_core_ln(x) + c
  }
}

///|
/// Asinh transform: a * asinh(x / b) + c.
pub fn flow_core_asinh_transform(
  x : Double,
  a : Double,
  b : Double,
  c : Double,
) -> Double {
  if b == 0.0 {
    return c
  }
  a * flow_core_asinh(x / b) + c
}

///|
/// Logicle transform with T (top), W (width), M (max), A (additional) parameters.
pub fn flow_core_logicle_transform(
  x : Double,
  t : Double,
  w : Double,
  m : Double,
  a : Double,
) -> Double {
  if t <= 0.0 || m <= 0.0 {
    return 0.0
  }
  let xn = x / t
  if xn >= 1.0 {
    return m
  }
  let total = m + a
  if total <= 0.0 {
    return 0.0
  }
  let a_coeff = flow_core_exp(w * 2.302585093) - 1.0
  let denom = flow_core_ln(1.0 + a_coeff)
  if denom == 0.0 {
    return 0.0
  }
  if xn <= 0.0 {
    return xn * m * a_coeff / denom
  }
  m * flow_core_ln(1.0 + a_coeff * xn) / denom
}

///|
/// Apply a transformation to all values in a FlowFrame.
pub fn flow_core_transform(
  frame : FlowFrame,
  transform : Transformation,
) -> FlowFrame {
  let new_exprs : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < frame.exprs.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < frame.exprs[i].length() {
      let x = frame.exprs[i][j]
      let transformed = match transform.transform_type {
        LogTransform => {
          let a = if transform.parameters.length() > 0 {
            transform.parameters[0]
          } else {
            1.0
          }
          let b = if transform.parameters.length() > 1 {
            transform.parameters[1]
          } else {
            10.0
          }
          let c = if transform.parameters.length() > 2 {
            transform.parameters[2]
          } else {
            0.0
          }
          flow_core_log_transform(x, a, b, c)
        }
        AsinhTransform => {
          let a = if transform.parameters.length() > 0 {
            transform.parameters[0]
          } else {
            1.0
          }
          let b = if transform.parameters.length() > 1 {
            transform.parameters[1]
          } else {
            1.0
          }
          let c = if transform.parameters.length() > 2 {
            transform.parameters[2]
          } else {
            0.0
          }
          flow_core_asinh_transform(x, a, b, c)
        }
        LogicleTransform => {
          let t = if transform.parameters.length() > 0 {
            transform.parameters[0]
          } else {
            262144.0
          }
          let w = if transform.parameters.length() > 1 {
            transform.parameters[1]
          } else {
            0.5
          }
          let m = if transform.parameters.length() > 2 {
            transform.parameters[2]
          } else {
            4.5
          }
          let a = if transform.parameters.length() > 3 {
            transform.parameters[3]
          } else {
            0.0
          }
          flow_core_logicle_transform(x, t, w, m, a)
        }
        LinearTransform => {
          let a = if transform.parameters.length() > 0 {
            transform.parameters[0]
          } else {
            1.0
          }
          let b = if transform.parameters.length() > 1 {
            transform.parameters[1]
          } else {
            0.0
          }
          a * x + b
        }
      }
      row.push(transformed)
      j = j + 1
    }
    new_exprs.push(row)
    i = i + 1
  }
  FlowFrame::{
    exprs: new_exprs,
    params: frame.params,
    description: frame.description,
    n_events: frame.n_events,
    n_params: frame.n_params,
  }
}

///|
/// Compute matrix inverse using Gaussian elimination with partial pivoting.
fn flow_core_matrix_inverse(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = matrix.length()
  if n == 0 {
    return Array::new()
  }
  let aug : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(matrix[i][j])
      j = j + 1
    }
    j = 0
    while j < n {
      if i == j {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
      j = j + 1
    }
    aug.push(row)
    i = i + 1
  }
  let mut col = 0
  while col < n {
    let mut max_val = aug[col][col].abs()
    let mut pivot_row = col
    let mut r = col + 1
    while r < n {
      let av = aug[r][col].abs()
      if av > max_val {
        max_val = av
        pivot_row = r
      }
      r = r + 1
    }
    if max_val != 0.0 {
      if pivot_row != col {
        let tmp = aug[col]
        aug[col] = aug[pivot_row]
        aug[pivot_row] = tmp
      }
      let pivot = aug[col][col]
      if pivot != 0.0 {
        let mut j = 0
        while j < 2 * n {
          aug[col][j] = aug[col][j] / pivot
          j = j + 1
        }
        let mut r2 = 0
        while r2 < n {
          if r2 != col {
            let factor = aug[r2][col]
            let mut j = 0
            while j < 2 * n {
              aug[r2][j] = aug[r2][j] - factor * aug[col][j]
              j = j + 1
            }
          }
          r2 = r2 + 1
        }
      }
    }
    col = col + 1
  }
  let result : Array[Array[Double]] = Array::new()
  i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = n
    while j < 2 * n {
      row.push(aug[i][j])
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Apply fluorescence compensation to a FlowFrame.
pub fn flow_core_compensate(
  frame : FlowFrame,
  comp_matrix : CompensationMatrix,
) -> FlowFrame {
  let n_channels = comp_matrix.fluorochromes.length()
  if n_channels == 0 || comp_matrix.matrix.length() == 0 {
    return frame
  }
  let inv = flow_core_matrix_inverse(comp_matrix.matrix)
  if inv.length() == 0 {
    return frame
  }
  let channel_indices : Array[Int] = Array::new()
  let mut i = 0
  while i < comp_matrix.fluorochromes.length() {
    let fluor = comp_matrix.fluorochromes[i]
    let mut idx = -1
    let mut j = 0
    while j < frame.params.length() {
      if frame.params[j].name == fluor {
        idx = j
        break
      }
      j = j + 1
    }
    channel_indices.push(idx)
    i = i + 1
  }
  let new_exprs : Array[Array[Double]] = Array::new()
  i = 0
  while i < frame.exprs.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < frame.exprs[i].length() {
      row.push(frame.exprs[i][j])
      j = j + 1
    }
    let mut c = 0
    while c < n_channels {
      let idx = channel_indices[c]
      if idx >= 0 && idx < row.length() {
        let mut sum = 0.0
        let mut k = 0
        while k < n_channels {
          let kidx = channel_indices[k]
          if kidx >= 0 && kidx < frame.exprs[i].length() {
            sum = sum + inv[c][k] * frame.exprs[i][kidx]
          }
          k = k + 1
        }
        row[idx] = sum
      }
      c = c + 1
    }
    new_exprs.push(row)
    i = i + 1
  }
  FlowFrame::{
    exprs: new_exprs,
    params: frame.params,
    description: frame.description,
    n_events: frame.n_events,
    n_params: frame.n_params,
  }
}

///|
/// Test if a point is inside a polygon using ray casting algorithm.
pub fn flow_core_point_in_polygon(
  x : Double,
  y : Double,
  polygon : Array[(Double, Double)],
) -> Bool {
  let n = polygon.length()
  if n < 3 {
    return false
  }
  let mut inside = false
  let mut j = n - 1
  let mut i = 0
  while i < n {
    let (xi, yi) = polygon[i]
    let (xj, yj) = polygon[j]
    let cond1 = (yi > y) != (yj > y)
    if cond1 {
      let x_intersect = (xj - xi) * (y - yi) / (yj - yi) + xi
      if x < x_intersect {
        inside = not(inside)
      }
    }
    j = i
    i = i + 1
  }
  inside
}

///|
/// Test if a point is inside an ellipse.
pub fn flow_core_point_in_ellipse(
  x : Double,
  y : Double,
  center : (Double, Double),
  radii : (Double, Double),
) -> Bool {
  let (cx, cy) = center
  let (rx, ry) = radii
  if rx == 0.0 || ry == 0.0 {
    return false
  }
  let dx = (x - cx) / rx
  let dy = (y - cy) / ry
  dx * dx + dy * dy <= 1.0
}

///|
/// Find the index of a parameter in a FlowFrame by name.
fn flow_core_find_param(frame : FlowFrame, name : String) -> Int {
  let mut i = 0
  while i < frame.params.length() {
    if frame.params[i].name == name {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Apply a rectangle gate to a FlowFrame.
pub fn flow_core_rectangle_gate(frame : FlowFrame, gate : Gate) -> GatingResult {
  let filtered : Array[Array[Double]] = Array::new()
  let n_total = frame.exprs.length()
  let n_params = gate.parameters.length()
  let param_indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n_params {
    let idx = flow_core_find_param(frame, gate.parameters[i])
    param_indices.push(idx)
    i = i + 1
  }
  i = 0
  while i < frame.exprs.length() {
    let row = frame.exprs[i]
    let mut keep = true
    let mut j = 0
    while j < n_params {
      let idx = param_indices[j]
      if idx >= 0 && idx < row.length() && j < gate.coordinates.length() {
        let (min_val, max_val) = gate.coordinates[j]
        let val = row[idx]
        if val < min_val || val > max_val {
          keep = false
          break
        }
      }
      j = j + 1
    }
    if keep {
      filtered.push(row)
    }
    i = i + 1
  }
  let n_filtered = filtered.length()
  let pct = if n_total > 0 {
    n_filtered.to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  let fail_pct = if n_total > 0 {
    (n_total - n_filtered).to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  GatingResult::{
    filtered_events: filtered,
    n_filtered,
    n_total,
    percentages: [pct, fail_pct],
  }
}

///|
/// Apply a polygon gate to a FlowFrame (2D gate).
pub fn flow_core_polygon_gate(frame : FlowFrame, gate : Gate) -> GatingResult {
  let filtered : Array[Array[Double]] = Array::new()
  let n_total = frame.exprs.length()
  if gate.parameters.length() < 2 || gate.coordinates.length() < 3 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 100.0],
    }
  }
  let idx_x = flow_core_find_param(frame, gate.parameters[0])
  let idx_y = flow_core_find_param(frame, gate.parameters[1])
  if idx_x < 0 || idx_y < 0 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 100.0],
    }
  }
  let polygon = gate.coordinates
  let mut i = 0
  while i < frame.exprs.length() {
    let row = frame.exprs[i]
    if idx_x < row.length() && idx_y < row.length() {
      let x = row[idx_x]
      let y = row[idx_y]
      if flow_core_point_in_polygon(x, y, polygon) {
        filtered.push(row)
      }
    }
    i = i + 1
  }
  let n_filtered = filtered.length()
  let pct = if n_total > 0 {
    n_filtered.to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  let fail_pct = if n_total > 0 {
    (n_total - n_filtered).to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  GatingResult::{
    filtered_events: filtered,
    n_filtered,
    n_total,
    percentages: [pct, fail_pct],
  }
}

///|
/// Apply an ellipsoid gate to a FlowFrame (2D gate).
pub fn flow_core_ellipsoid_gate(frame : FlowFrame, gate : Gate) -> GatingResult {
  let filtered : Array[Array[Double]] = Array::new()
  let n_total = frame.exprs.length()
  if gate.parameters.length() < 2 || gate.coordinates.length() < 2 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 100.0],
    }
  }
  let idx_x = flow_core_find_param(frame, gate.parameters[0])
  let idx_y = flow_core_find_param(frame, gate.parameters[1])
  if idx_x < 0 || idx_y < 0 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 100.0],
    }
  }
  let center = gate.coordinates[0]
  let radii = gate.coordinates[1]
  let mut i = 0
  while i < frame.exprs.length() {
    let row = frame.exprs[i]
    if idx_x < row.length() && idx_y < row.length() {
      let x = row[idx_x]
      let y = row[idx_y]
      if flow_core_point_in_ellipse(x, y, center, radii) {
        filtered.push(row)
      }
    }
    i = i + 1
  }
  let n_filtered = filtered.length()
  let pct = if n_total > 0 {
    n_filtered.to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  let fail_pct = if n_total > 0 {
    (n_total - n_filtered).to_double() * 100.0 / n_total.to_double()
  } else {
    0.0
  }
  GatingResult::{
    filtered_events: filtered,
    n_filtered,
    n_total,
    percentages: [pct, fail_pct],
  }
}

///|
/// Apply a quadrant gate to a FlowFrame (2D gate, returns Q2 upper-right).
pub fn flow_core_quadrant_gate(frame : FlowFrame, gate : Gate) -> GatingResult {
  let filtered : Array[Array[Double]] = Array::new()
  let n_total = frame.exprs.length()
  if gate.parameters.length() < 2 || gate.coordinates.length() < 1 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 0.0, 0.0, 0.0],
    }
  }
  let idx_x = flow_core_find_param(frame, gate.parameters[0])
  let idx_y = flow_core_find_param(frame, gate.parameters[1])
  if idx_x < 0 || idx_y < 0 {
    return GatingResult::{
      filtered_events: filtered,
      n_filtered: 0,
      n_total,
      percentages: [0.0, 0.0, 0.0, 0.0],
    }
  }
  let (x_thresh, y_thresh) = gate.coordinates[0]
  let mut q1 = 0
  let mut q2 = 0
  let mut q3 = 0
  let mut q4 = 0
  let mut i = 0
  while i < frame.exprs.length() {
    let row = frame.exprs[i]
    if idx_x < row.length() && idx_y < row.length() {
      let x = row[idx_x]
      let y = row[idx_y]
      let x_pos = x >= x_thresh
      let y_pos = y >= y_thresh
      if x_pos && y_pos {
        q2 = q2 + 1
        filtered.push(row)
      } else if x_pos && not(y_pos) {
        q4 = q4 + 1
      } else if not(x_pos) && y_pos {
        q1 = q1 + 1
      } else {
        q3 = q3 + 1
      }
    }
    i = i + 1
  }
  let n_filtered = q2
  let denom = if n_total > 0 { n_total.to_double() } else { 1.0 }
  GatingResult::{
    filtered_events: filtered,
    n_filtered,
    n_total,
    percentages: [
      q1.to_double() * 100.0 / denom,
      q2.to_double() * 100.0 / denom,
      q3.to_double() * 100.0 / denom,
      q4.to_double() * 100.0 / denom,
    ],
  }
}

///|
/// Main gating entry point: dispatch by gate type.
pub fn flow_core_filter_events(frame : FlowFrame, gate : Gate) -> GatingResult {
  match gate.gate_type {
    RectangleGate => flow_core_rectangle_gate(frame, gate)
    PolygonGate => flow_core_polygon_gate(frame, gate)
    EllipsoidGate => flow_core_ellipsoid_gate(frame, gate)
    QuadrantGate => flow_core_quadrant_gate(frame, gate)
    BooleanGate => flow_core_rectangle_gate(frame, gate)
  }
}

///|
/// Compute spillover (flattened matrix) from a 2D matrix.
pub fn flow_core_compute_spillover(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < matrix.length() {
    let mut j = 0
    while j < matrix[i].length() {
      result.push(matrix[i][j])
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Create a compensation matrix from fluorochromes and values.
pub fn flow_core_create_compensation_matrix(
  fluorochromes : Array[String],
  values : Array[Array[Double]],
) -> CompensationMatrix {
  let n = fluorochromes.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      if i < values.length() && j < values[i].length() {
        row.push(values[i][j])
      } else {
        row.push(0.0)
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  CompensationMatrix::new(matrix, fluorochromes)
}

///|
/// Normalize a FlowFrame using the specified method.
pub fn flow_core_normalize(
  frame : FlowFrame,
  norm_method : String,
) -> FlowFrame {
  let n_events = frame.exprs.length()
  if n_events == 0 {
    return frame
  }
  let n_params = frame.exprs[0].length()
  match norm_method {
    "scale" => {
      let mins : Array[Double] = Array::new()
      let maxs : Array[Double] = Array::new()
      let mut j = 0
      while j < n_params {
        mins.push(0.0)
        maxs.push(0.0)
        j = j + 1
      }
      j = 0
      while j < n_params {
        let mut mn = 999999999999999999999999999.0
        let mut mx = -999999999999999999999999999.0
        let mut i = 0
        while i < n_events {
          if j < frame.exprs[i].length() {
            let v = frame.exprs[i][j]
            if v < mn {
              mn = v
            }
            if v > mx {
              mx = v
            }
          }
          i = i + 1
        }
        mins[j] = mn
        maxs[j] = mx
        j = j + 1
      }
      let new_exprs : Array[Array[Double]] = Array::new()
      let mut i = 0
      while i < n_events {
        let row : Array[Double] = Array::new()
        let mut j = 0
        while j < n_params {
          let v = if j < frame.exprs[i].length() {
            frame.exprs[i][j]
          } else {
            0.0
          }
          let range = maxs[j] - mins[j]
          let norm = if range > 0.0 { (v - mins[j]) / range } else { 0.0 }
          row.push(norm)
          j = j + 1
        }
        new_exprs.push(row)
        i = i + 1
      }
      FlowFrame::{
        exprs: new_exprs,
        params: frame.params,
        description: frame.description,
        n_events: frame.n_events,
        n_params: frame.n_params,
      }
    }
    "log" => {
      let new_exprs : Array[Array[Double]] = Array::new()
      let mut i = 0
      while i < n_events {
        let row : Array[Double] = Array::new()
        let mut j = 0
        while j < n_params {
          let v = if j < frame.exprs[i].length() {
            frame.exprs[i][j]
          } else {
            0.0
          }
          if v > 0.0 {
            row.push(flow_core_ln(v))
          } else {
            row.push(0.0)
          }
          j = j + 1
        }
        new_exprs.push(row)
        i = i + 1
      }
      FlowFrame::{
        exprs: new_exprs,
        params: frame.params,
        description: frame.description,
        n_events: frame.n_events,
        n_params: frame.n_params,
      }
    }
    _ => frame
  }
}

///|
/// Generate a summary string of a FlowFrame.
pub fn flow_core_summary(frame : FlowFrame) -> String {
  let mut result = "FlowFrame Summary\n"
  result = result + "================\n"
  result = result + "Events: " + frame.n_events.to_string() + "\n"
  result = result + "Parameters: " + frame.n_params.to_string() + "\n"
  result = result + "\nParameter Statistics:\n"
  result = result + "Name\tMin\tMax\tMean\tRange\n"
  let n_events = frame.exprs.length()
  let mut i = 0
  while i < frame.params.length() {
    let p = frame.params[i]
    let mut mn = 999999999999999999999999999.0
    let mut mx = -999999999999999999999999999.0
    let mut sum = 0.0
    let mut count = 0
    let mut j = 0
    while j < n_events {
      if i < frame.exprs[j].length() {
        let v = frame.exprs[j][i]
        if v < mn {
          mn = v
        }
        if v > mx {
          mx = v
        }
        sum = sum + v
        count = count + 1
      }
      j = j + 1
    }
    let mean = if count > 0 { sum / count.to_double() } else { 0.0 }
    let range = mx - mn
    result = result +
      p.name +
      "\t" +
      mn.to_string() +
      "\t" +
      mx.to_string() +
      "\t" +
      mean.to_string() +
      "\t" +
      range.to_string() +
      "\n"
    i = i + 1
  }
  result = result +
    "\nDescription keys: " +
    frame.description.length().to_string() +
    "\n"
  let keys = frame.description.keys().collect()
  i = 0
  while i < keys.length() {
    result = result + "  " + keys[i] + " = " + frame.description[keys[i]] + "\n"
    i = i + 1
  }
  result
}

///|
/// Merge multiple FlowFrames into one.
pub fn flow_core_merge_frames(frames : Array[FlowFrame]) -> FlowFrame {
  let merged_exprs : Array[Array[Double]] = Array::new()
  let merged_params : Array[ParameterDescription] = Array::new()
  let merged_desc : Map[String, String] = Map([], capacity=16)
  let mut n_events = 0
  let mut n_params = 0
  let mut first = true
  let mut i = 0
  while i < frames.length() {
    let f = frames[i]
    let mut j = 0
    while j < f.exprs.length() {
      merged_exprs.push(f.exprs[j])
      j = j + 1
    }
    if first {
      let mut k = 0
      while k < f.params.length() {
        merged_params.push(f.params[k])
        k = k + 1
      }
      let keys = f.description.keys().collect()
      k = 0
      while k < keys.length() {
        merged_desc[keys[k]] = f.description[keys[k]]
        k = k + 1
      }
      n_params = f.n_params
      first = false
    }
    n_events = n_events + f.n_events
    i = i + 1
  }
  FlowFrame::{
    exprs: merged_exprs,
    params: merged_params,
    description: merged_desc,
    n_events,
    n_params,
  }
}

///|
/// Create an example FCS file with 5 events and 4 parameters.
pub fn flow_core_create_example_fcs() -> FCSFile {
  let text_segment : Map[String, String] = Map([], capacity=16)
  text_segment["$FCSVERSION"] = "3.1"
  text_segment["$TOT"] = "5"
  text_segment["$PAR"] = "4"
  text_segment["$P1N"] = "FSC-A"
  text_segment["$P2N"] = "SSC-A"
  text_segment["$P3N"] = "FITC-A"
  text_segment["$P4N"] = "PE-A"
  text_segment["$DATATYPE"] = "D"
  let data_segment : Array[Array[Double]] = [
    [100.0, 200.0, 10.0, 5.0],
    [150.0, 250.0, 20.0, 15.0],
    [200.0, 300.0, 30.0, 25.0],
    [250.0, 350.0, 40.0, 35.0],
    [300.0, 400.0, 50.0, 45.0],
  ]
  let parameter_names = ["FSC-A", "SSC-A", "FITC-A", "PE-A"]
  let parameter_ranges = [
    (0.0, 1000.0),
    (0.0, 1000.0),
    (0.0, 100.0),
    (0.0, 100.0),
  ]
  FCSFile::{
    text_segment,
    data_segment,
    analysis_segment: "",
    n_events: 5,
    n_parameters: 4,
    parameter_names,
    parameter_ranges,
  }
}

///|
/// Create an example FlowFrame.
pub fn flow_core_create_example_frame() -> FlowFrame {
  let fcs = flow_core_create_example_fcs()
  FlowFrame::from_fcs(fcs)
}

///|
/// Create an example FlowSet with 2 frames.
pub fn flow_core_create_example_set() -> FlowSet {
  let frame1 = flow_core_create_example_frame()
  let data2 : Array[Array[Double]] = [
    [110.0, 210.0, 11.0, 6.0],
    [160.0, 260.0, 21.0, 16.0],
    [210.0, 310.0, 31.0, 26.0],
  ]
  let params : Array[ParameterDescription] = Array::new()
  let mut i = 0
  while i < frame1.params.length() {
    let p = frame1.params[i]
    params.push(ParameterDescription::{
      name: p.name,
      range: p.range,
      min_range: p.min_range,
      max_range: p.max_range,
      compensation: Array::new(),
    })
    i = i + 1
  }
  let frame2 = FlowFrame::new(data2, params)
  let frames = [frame1, frame2]
  let sample_names = ["sample1", "sample2"]
  FlowSet::new(frames, sample_names)
}

///|
/// Create an example rectangle gate on FITC-A and PE-A.
pub fn flow_core_create_example_gate() -> Gate {
  Gate::new(
    RectangleGate,
    ["FITC-A", "PE-A"],
    [(10.0, 50.0), (5.0, 45.0)],
    "example_rect_gate",
  )
}

///|
/// Create an example compensation matrix for 3 fluorochromes.
pub fn flow_core_create_example_compensation() -> CompensationMatrix {
  let fluorochromes = ["FITC-A", "PE-A", "PerCP-A"]
  let values = [[1.0, 0.1, 0.05], [0.08, 1.0, 0.12], [0.03, 0.15, 1.0]]
  flow_core_create_compensation_matrix(fluorochromes, values)
}
