///|
/// Smart - SMART (Simple Modular Architecture Research Tool) domain database format parser.
/// Inspired by Biopython Bio.Smart.
///
/// SMART is a protein domain database that identifies and annotates protein domains.
/// The SMART output format is tab-delimited and contains domain annotations for
/// protein sequences, including domain boundaries, evalues, and functional descriptions.
///
/// Reference: http://smart.embl-heidelberg.de/
/// Schultz et al. (2000) Nucl. Acids Res. 28:231-234.

// ============================================================================
// SMART Data Structures
// ============================================================================

///|
/// Represents a single SMART domain annotation entry.
pub struct SmartDomain {
  domain_id : String
  domain_name : String
  domain_type : String
  start : Int
  end : Int
  evalue : Double
  sequence : String
  description : String
  pfam_id : String
  interpro_id : String
  go_annotations : Array[String]
}

///|
/// Construct a SmartDomain.
pub fn SmartDomain::new(
  domain_id : String,
  domain_name : String,
  start : Int,
  end : Int,
  evalue? : Double = 0.0
) -> SmartDomain {
  SmartDomain::{
    domain_id,
    domain_name,
    domain_type: "",
    start,
    end,
    evalue,
    sequence: "",
    description: "",
    pfam_id: "",
    interpro_id: "",
    go_annotations: [],
  }
}

///|
/// Represents a complete SMART annotation result for a protein sequence.
pub struct SmartResult {
  mut sequence_id : String
  mut sequence_length : Int
  domains : Array[SmartDomain]
  pfam_domains : Array[SmartDomain]
  signal_peptides : Array[SmartDomain]
  transmembrane : Array[SmartDomain]
}

///|
/// Construct a SmartResult.
pub fn SmartResult::new(
  sequence_id : String,
  sequence_length? : Int = 0
) -> SmartResult {
  SmartResult::{
    sequence_id,
    sequence_length,
    domains: [],
    pfam_domains: [],
    signal_peptides: [],
    transmembrane: [],
  }
}

// ============================================================================
// SMART Parser (Tab-delimited text to SmartResult)
// ============================================================================

///|
/// Parse SMART tab-delimited output into a SmartResult.
///
/// The SMART output format is tab-delimited with the following columns:
/// - #domain_id: SMART domain identifier (e.g., "SM00001")
/// - domain_name: Short domain name (e.g., "ABC_membrane")
/// - domain_type: Type of domain (e.g., "SMART", "Pfam", "SignalP")
/// - start: Start position in protein sequence
/// - end: End position in protein sequence
/// - evalue: E-value for the domain hit
/// - sequence: Matched sequence region
/// - description: Functional description of the domain
/// - pfam_id: Corresponding Pfam accession
/// - interpro_id: Corresponding InterPro accession
/// - go_annotations: GO terms separated by semicolons
pub fn smart_parse(content : String) -> SmartResult {
  let result = SmartResult::new("")
  let lines = smart_split_lines(content)
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    let trimmed = line.trim().to_string()
    if trimmed.length() == 0 {
      i = i + 1
      continue
    }
    if trimmed[0:1].to_owned() == "#" {
      if trimmed.contains("sequence_id") || trimmed.contains("Sequence") {
        result.sequence_id = smart_extract_value(trimmed)
      } else if trimmed.contains("sequence_length") {
        result.sequence_length = smart_parse_int(smart_extract_value(trimmed))
      }
      i = i + 1
      continue
    }
    let cols = smart_split_tab(trimmed)
    if cols.length() >= 5 {
      let domain = smart_parse_domain_line(cols)
      match domain.domain_type {
        "SignalP" | "signal_peptide" => result.signal_peptides.push(domain)
        "TM" | "transmembrane" => result.transmembrane.push(domain)
        "Pfam" => result.pfam_domains.push(domain)
        _ => result.domains.push(domain)
      }
    }
    i = i + 1
  }
  result
}

///|
/// Parse a single domain line from column array.
fn smart_parse_domain_line(cols : Array[String]) -> SmartDomain {
  let domain_id = if cols.length() > 0 { cols[0] } else { "" }
  let domain_name = if cols.length() > 1 { cols[1] } else { "" }
  let domain_type = if cols.length() > 2 { cols[2] } else { "SMART" }
  let start = if cols.length() > 3 { smart_parse_int(cols[3]) } else { 0 }
  let end = if cols.length() > 4 { smart_parse_int(cols[4]) } else { 0 }
  let evalue = if cols.length() > 5 { smart_parse_double(cols[5]) } else { 0.0 }
  let sequence = if cols.length() > 6 { cols[6] } else { "" }
  let description = if cols.length() > 7 { cols[7] } else { "" }
  let pfam_id = if cols.length() > 8 { cols[8] } else { "" }
  let interpro_id = if cols.length() > 9 { cols[9] } else { "" }
  let go_str = if cols.length() > 10 { cols[10] } else { "" }
  let go_annotations = smart_split_semicolon(go_str)

  SmartDomain::{
    domain_id,
    domain_name,
    domain_type,
    start,
    end,
    evalue,
    sequence,
    description,
    pfam_id,
    interpro_id,
    go_annotations,
  }
}

// ============================================================================
// SMART Search / Query Functions
// ============================================================================

///|
/// Find all domains matching a given domain name (substring match).
pub fn smart_find_domains(result : SmartResult, name : String) -> Array[SmartDomain] {
  let matches : Array[SmartDomain] = Array::new()
  let mut i = 0
  while i < result.domains.length() {
    if result.domains[i].domain_name.contains(name) {
      matches.push(result.domains[i])
    }
    i = i + 1
  }
  let mut j = 0
  while j < result.pfam_domains.length() {
    if result.pfam_domains[j].domain_name.contains(name) {
      matches.push(result.pfam_domains[j])
    }
    j = j + 1
  }
  matches
}

///|
/// Find all domains overlapping with a given region [start, end].
pub fn smart_find_overlapping_domains(result : SmartResult, start : Int, end : Int) -> Array[SmartDomain] {
  let matches : Array[SmartDomain] = Array::new()
  let mut i = 0
  while i < result.domains.length() {
    let d = result.domains[i]
    if d.start <= end && d.end >= start {
      matches.push(d)
    }
    i = i + 1
  }
  let mut j = 0
  while j < result.pfam_domains.length() {
    let d = result.pfam_domains[j]
    if d.start <= end && d.end >= start {
      matches.push(d)
    }
    j = j + 1
  }
  matches
}

///|
/// Get the lowest E-value domain in the result.
pub fn smart_best_domain(result : SmartResult) -> SmartDomain? {
  let mut best : SmartDomain? = None
  let mut i = 0
  while i < result.domains.length() {
    match best {
      Some(cur) => {
        if result.domains[i].evalue < cur.evalue {
          best = Some(result.domains[i])
        }
      }
      None => best = Some(result.domains[i])
    }
    i = i + 1
  }
  let mut j = 0
  while j < result.pfam_domains.length() {
    match best {
      Some(cur) => {
        if result.pfam_domains[j].evalue < cur.evalue {
          best = Some(result.pfam_domains[j])
        }
      }
      None => best = Some(result.pfam_domains[j])
    }
    j = j + 1
  }
  best
}

///|
/// Get total number of domains (SMART + Pfam).
pub fn smart_total_domains(result : SmartResult) -> Int {
  result.domains.length() + result.pfam_domains.length()
}

///|
/// Filter domains by E-value threshold.
pub fn smart_filter_evalue(result : SmartResult, max_evalue : Double) -> Array[SmartDomain] {
  let filtered : Array[SmartDomain] = Array::new()
  let mut i = 0
  while i < result.domains.length() {
    if result.domains[i].evalue <= max_evalue {
      filtered.push(result.domains[i])
    }
    i = i + 1
  }
  let mut j = 0
  while j < result.pfam_domains.length() {
    if result.pfam_domains[j].evalue <= max_evalue {
      filtered.push(result.pfam_domains[j])
    }
    j = j + 1
  }
  filtered
}

///|
/// Generate a sample SMART output string for testing.
pub fn smart_sample() -> String {
  let lines = [
    "#sequence_id\tProtein123",
    "#sequence_length\t350",
    "SM00001\tABC_membrane\tSMART\t10\t100\t1.5e-30\tMFWVLVVVGGVLALAVLVLAVL\tABC transporter membrane domain\tPF00001\tIPR00001\tGO:0005215;GO:0016020",
    "SM00002\tABC_subdomain\tSMART\t110\t250\t2.3e-25\tGCGKSTEQRLLLDELPGRRL\tABC transporter nucleotide binding domain\tPF00002\tIPR00002\tGO:0005524",
    "SM00003\tABC_membrane\tSMART\t260\t340\t5.1e-28\tVVLLVAVVGCGVALALVAVLL\tABC transporter membrane domain\tPF00001\tIPR00001\tGO:0005215;GO:0016020",
    "PF00001\tABC_membrane\tPfam\t10\t100\t1.5e-30\tMFWVLVVVGGVLALAVLVLAVL\tABC transporter membrane domain\tPF00001\tIPR00001\tGO:0005215;GO:0016020",
    "SignalP\tSignalP\tSignalP\t1\t22\t0.001\tMFWVLVVVGGVLALAVLVLA\tSignal peptide\t\t\t",
    "TM1\tTM\ttransmembrane\t240\t262\t0.001\tVVLLVAVVGCGVALALVAVLL\tTransmembrane region\t\t\t",
  ]
  let mut result = ""
  let mut i = 0
  while i < lines.length() {
    if i > 0 {
      result = result + "\n"
    }
    result = result + lines[i]
    i = i + 1
  }
  result
}

///|
/// Convert a SmartResult to a summary string.
pub fn smart_summary(result : SmartResult) -> String {
  let mut out = ""
  out = out + "SMART Analysis Summary\n"
  out = out + "Sequence: " + result.sequence_id + "\n"
  out = out + "Length: " + result.sequence_length.to_string() + "\n"
  out = out + "SMART domains: " + result.domains.length().to_string() + "\n"
  out = out + "Pfam domains: " + result.pfam_domains.length().to_string() + "\n"
  out = out + "Signal peptides: " + result.signal_peptides.length().to_string() + "\n"
  out = out + "Transmembrane: " + result.transmembrane.length().to_string() + "\n"
  let mut i = 0
  while i < result.domains.length() {
    let d = result.domains[i]
    out = out + "  " + d.domain_id + "\t" + d.domain_name + "\t" + d.start.to_string() + "-" + d.end.to_string() + "\tE=" + d.evalue.to_string() + "\n"
    i = i + 1
  }
  out
}

// ============================================================================
// Private utility functions
// ============================================================================

///|
/// Split a string into lines.
fn smart_split_lines(content : String) -> Array[String] {
  let lines : Array[String] = Array::new()
  let mut start = 0
  let mut i = 0
  while i < content.length() {
    let ch = content[i:i + 1].to_owned()
    if ch == "\n" || ch == "\r" {
      if start < i {
        lines.push(content[start:i].trim().to_string())
      }
      if ch == "\r" && i + 1 < content.length() && content[i + 1:i + 2].to_owned() == "\n" {
        i = i + 1
      }
      start = i + 1
    }
    i = i + 1
  }
  if start < content.length() {
    lines.push(content[start:content.length()].trim().to_string())
  }
  lines
}

///|
/// Split a string by tab delimiter.
fn smart_split_tab(line : String) -> Array[String] {
  let cols : Array[String] = Array::new()
  let mut start = 0
  let mut i = 0
  while i < line.length() {
    if line[i:i + 1].to_owned() == "\t" {
      cols.push(line[start:i].to_string())
      start = i + 1
    }
    i = i + 1
  }
  cols.push(line[start:line.length()].to_string())
  cols
}

///|
/// Split a string by semicolon delimiter.
fn smart_split_semicolon(s : String) -> Array[String] {
  let items : Array[String] = Array::new()
  if s.length() == 0 { return items }
  let mut start = 0
  let mut i = 0
  while i < s.length() {
    if s[i:i + 1].to_owned() == ";" {
      let item = s[start:i].trim().to_string()
      if item.length() > 0 {
        items.push(item)
      }
      start = i + 1
    }
    i = i + 1
  }
  let last = s[start:s.length()].trim().to_string()
  if last.length() > 0 {
    items.push(last)
  }
  items
}

///|
/// Extract value after a colon or tab from a comment line.
fn smart_extract_value(line : String) -> String {
  let colon_pos = line.find(":")
  match colon_pos {
    Some(pos) => line[pos + 1:line.length()].trim().to_string()
    None => {
      let tab_pos = line.find("\t")
      match tab_pos {
        Some(tp) => line[tp + 1:line.length()].trim().to_string()
        None => line.trim().to_string()
      }
    }
  }
}

///|
/// Parse an integer from a string, returning 0 on failure.
fn smart_parse_int(s : String) -> Int {
  let sv = s.to_string().trim()
  let n = sv.length()
  if n == 0 { return 0 }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 { break }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Parse a double from a string, returning 0.0 on failure.
fn smart_parse_double(s : String) -> Double {
  let sv = s.to_string().trim()
  if sv.length() == 0 { return 0.0 }
  let n = sv.length()
  let mut result = 0.0
  let mut decimal = 0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  if sv.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 { break }
    result = result * 10.0 + (c - zero).to_double()
    i = i + 1
  }
  if i < n && sv.unsafe_get(i) == '.' {
    i = i + 1
    while i < n {
      let c = sv.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 { break }
      result = result * 10.0 + (c - zero).to_double()
      decimal = decimal + 1
      i = i + 1
    }
  }
  while decimal > 0 {
    result = result / 10.0
    decimal = decimal - 1
  }
  // Handle scientific notation (e.g., 1.5e-30)
  let mut exponent = 0
  let mut exp_sign = 1
  if i < n && (sv.unsafe_get(i) == 'e' || sv.unsafe_get(i) == 'E') {
    i = i + 1
    if i < n && sv.unsafe_get(i) == '-' {
      exp_sign = -1
      i = i + 1
    } else if i < n && sv.unsafe_get(i) == '+' {
      i = i + 1
    }
    while i < n {
      let c = sv.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 { break }
      exponent = exponent * 10 + (c - zero)
      i = i + 1
    }
    exponent = exponent * exp_sign
  }
  result * sign * @math.pow(10.0, exponent.to_double())
}
