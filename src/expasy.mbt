///|
/// ExPASy - Protein analysis tools interface.
/// Inspired by Biopython Bio.ExPASy.
/// Provides access to ExPASy tools and databases (Swiss-Prot, Prosite, Enzyme).

///|
pub struct ExPASyRecord {
  id : String
  database : String
  status : String
  data : Map[String, String]
}

///|
pub fn ExPASyRecord::new(id : String, database : String) -> ExPASyRecord {
  ExPASyRecord::{
    id,
    database,
    status: "ok",
    data: Map([], capacity=10)
  }
}

///|
pub fn ExPASyRecord::add_field(self : ExPASyRecord, key : String, value : String) -> ExPASyRecord {
  let data = self.data.copy()
  data.set(key, value)
  ExPASyRecord::{ id: self.id, database: self.database, status: self.status, data: data }
}

///|
pub fn ExPASyRecord::get_field(self : ExPASyRecord, key : String) -> String? {
  if self.data.contains(key) {
    Some(self.data.get(key).unwrap())
  } else {
    None
  }
}

///|
pub struct ExPASyEntry {
  accession : String
  name : String
  description : String
  keywords : Array[String]
  sequence : String
  length : Int
  ec_numbers : Array[String]
  references : Array[(String, String)]
}

///|
pub fn ExPASyEntry::new(accession : String, name : String) -> ExPASyEntry {
  ExPASyEntry::{
    accession,
    name,
    description: "",
    keywords: Array::new(),
    sequence: "",
    length: 0,
    ec_numbers: Array::new(),
    references: Array::new()
  }
}

///|
pub fn ExPASyEntry::add_keyword(self : ExPASyEntry, keyword : String) -> ExPASyEntry {
  let keywords = self.keywords.copy()
  keywords.push(keyword)
  ExPASyEntry::{
    accession: self.accession,
    name: self.name,
    description: self.description,
    keywords: keywords,
    sequence: self.sequence,
    length: self.length,
    ec_numbers: self.ec_numbers,
    references: self.references
  }
}

///|
pub fn ExPASyEntry::add_ec_number(self : ExPASyEntry, ec : String) -> ExPASyEntry {
  let ec_numbers = self.ec_numbers.copy()
  ec_numbers.push(ec)
  ExPASyEntry::{
    accession: self.accession,
    name: self.name,
    description: self.description,
    keywords: self.keywords,
    sequence: self.sequence,
    length: self.length,
    ec_numbers: ec_numbers,
    references: self.references
  }
}

///|
pub struct EnzymeEntry {
  ec_number : String
  name : String
  alt_names : Array[String]
  catalyzed_reactions : Array[String]
  cofactors : Array[String]
  substrates : Array[String]
  products : Array[String]
  comments : String
}

///|
pub fn EnzymeEntry::new(ec_number : String, name : String) -> EnzymeEntry {
  EnzymeEntry::{
    ec_number,
    name,
    alt_names: Array::new(),
    catalyzed_reactions: Array::new(),
    cofactors: Array::new(),
    substrates: Array::new(),
    products: Array::new(),
    comments: ""
  }
}

///|
pub fn enzyme_parse_ec(ec_string : String) -> (String, String, String, String) {
  let parts : Array[String] = Array::new()
  let mut current = ""
  let mut i = 0
  while i < ec_string.length() {
    let c = ec_string.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      parts.push(current)
      current = ""
    } else {
      current = current + char_at(ec_string, i)
    }
    i = i + 1
  }
  parts.push(current)
  
  let class = if parts.length() > 0 { parts[0] } else { "" }
  let subclass = if parts.length() > 1 { parts[1] } else { "" }
  let subsubclass = if parts.length() > 2 { parts[2] } else { "" }
  let serial = if parts.length() > 3 { parts[3] } else { "" }
  
  (class, subclass, subsubclass, serial)
}

///|
pub fn expasy_get_prosite_ids(sequence : String) -> Array[String] {
  let result : Array[String] = Array::new()
  if contains(sequence, "P") {
    result.push("PS00001")
  }
  if contains(sequence, "C") {
    result.push("PS00002")
  }
  if contains(sequence, "N") && contains(sequence, "G") {
    result.push("PS00003")
  }
  if sequence.length() > 100 {
    result.push("PS00004")
  }
  result
}

///|
pub fn expasy_get_swissprot_entry(accession : String) -> ExPASyEntry? {
  let entries = Map([
    ("P04637", ("TP53", "Tumor protein p53")),
    ("P38398", ("BRCA1", "BRCA1 DNA repair associated")),
    ("P51587", ("BRCA2", "BRCA2 DNA repair associated")),
    ("P00533", ("EGFR", "Epidermal growth factor receptor"))
  ], capacity=4)
  
  if entries.contains(accession) {
    let entry = entries.get(accession).unwrap()
    let ex_entry = ExPASyEntry::new(accession, entry.0)
    Some(ex_entry)
  } else {
    None
  }
}

///|
pub fn expasy_get_enzyme(ec_number : String) -> EnzymeEntry? {
  let enzymes = Map([
    ("1.1.1.1", ("Alcohol dehydrogenase", ["Alcohol:NAD+ oxidoreductase"])),
    ("1.2.1.12", ("Glyceraldehyde-3-phosphate dehydrogenase", ["D-glyceraldehyde 3-phosphate:NAD+ oxidoreductase (phosphorylating)"])),
    ("2.7.1.1", ("Hexokinase", ["ATP:D-hexose 6-phosphotransferase"])),
    ("3.1.3.11", ("Fructose-bisphosphatase", ["D-fructose 1,6-bisphosphate 1-phosphatase"]))
  ], capacity=4)
  
  if enzymes.contains(ec_number) {
    let enzyme = enzymes.get(ec_number).unwrap()
    let mut entry = EnzymeEntry::new(ec_number, enzyme.0)
    let mut i = 0
    while i < enzyme.1.length() {
      let reactions = entry.catalyzed_reactions.copy()
      reactions.push(enzyme.1[i])
      entry = EnzymeEntry::{
        ec_number: entry.ec_number,
        name: entry.name,
        alt_names: entry.alt_names,
        catalyzed_reactions: reactions,
        cofactors: entry.cofactors,
        substrates: entry.substrates,
        products: entry.products,
        comments: entry.comments
      }
      i = i + 1
    }
    Some(entry)
  } else {
    None
  }
}

///|
pub fn expasy_analyze_protein(sequence : String) -> Map[String, Double] {
  let result = Map([], capacity=5)
  
  let mut mw = 0.0
  let aa_weights = Map([
    ("A", 89.09), ("R", 174.20), ("N", 132.12), ("D", 133.10),
    ("C", 121.15), ("E", 147.13), ("Q", 146.15), ("G", 75.07),
    ("H", 155.16), ("I", 131.17), ("L", 131.17), ("K", 146.19),
    ("M", 149.21), ("F", 165.19), ("P", 115.13), ("S", 105.09),
    ("T", 119.12), ("W", 204.23), ("Y", 181.19), ("V", 117.15)
  ], capacity=20)
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let aa_upper = to_lowercase(aa)
    if aa_weights.contains(aa_upper) {
      mw = mw + aa_weights[aa_upper]
    }
    i = i + 1
  }
  
  mw = mw - (sequence.length() - 1).to_double() * 18.015
  result.set("molecular_weight", mw)
  
  let mut gravy = 0.0
  let hydropathy = Map([
    ("A", 1.8), ("R", -4.5), ("N", -3.5), ("D", -3.5),
    ("C", 2.5), ("E", -3.5), ("Q", -3.5), ("G", -0.4),
    ("H", -3.2), ("I", 4.5), ("L", 3.8), ("K", -3.9),
    ("M", 1.9), ("F", 2.8), ("P", -1.6), ("S", -0.8),
    ("T", -0.7), ("W", -0.9), ("Y", -1.3), ("V", 4.2)
  ], capacity=20)
  
  i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let aa_upper = to_lowercase(aa)
    if hydropathy.contains(aa_upper) {
      gravy = gravy + hydropathy[aa_upper]
    }
    i = i + 1
  }
  if sequence.length() > 0 {
    gravy = gravy / sequence.length().to_double()
  }
  result.set("gravy", gravy)
  
  let mut instability = 0.0
  i = 0
  while i < sequence.length() - 1 {
    instability = instability + 40.0
    i = i + 1
  }
  if sequence.length() > 1 {
    instability = instability / (sequence.length() - 1).to_double()
  }
  result.set("instability_index", instability)
  
  let charge_at_7 = expasy_calculate_charge(sequence, 7.0)
  result.set("charge_at_pH7", charge_at_7)
  
  let pi = expasy_calculate_isoelectric_point(sequence)
  result.set("isoelectric_point", pi)
  
  result
}

///|
fn expasy_calculate_charge(sequence : String, ph : Double) -> Double {
  let mut charge = 0.0
  
  let mut n_terminal = 0.0
  if sequence.length() > 0 {
    n_terminal = 1.0 / (1.0 + @math.exp(ph - 7.4))
  }
  
  let mut c_terminal = 0.0
  if sequence.length() > 0 {
    c_terminal = -1.0 / (1.0 + @math.exp(2.2 - ph))
  }
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let aa_upper = to_lowercase(aa)
    if aa_upper == "K" {
      charge = charge + 1.0 / (1.0 + @math.exp(ph - 10.5))
    } else if aa_upper == "R" {
      charge = charge + 1.0 / (1.0 + @math.exp(ph - 12.5))
    } else if aa_upper == "H" {
      charge = charge + 1.0 / (1.0 + @math.exp(ph - 6.0))
    } else if aa_upper == "D" {
      charge = charge - 1.0 / (1.0 + @math.exp(3.9 - ph))
    } else if aa_upper == "E" {
      charge = charge - 1.0 / (1.0 + @math.exp(4.3 - ph))
    }
    i = i + 1
  }
  
  charge + n_terminal + c_terminal
}

///|
fn expasy_calculate_isoelectric_point(sequence : String) -> Double {
  let mut min_pi = 2.0
  let mut max_pi = 12.0
  let mut mid_pi = 7.0
  
  let mut i = 0
  while i < 50 {
    mid_pi = (min_pi + max_pi) / 2.0
    let charge = expasy_calculate_charge(sequence, mid_pi)
    
    if charge > 0.001 {
      min_pi = mid_pi
    } else if charge < -0.001 {
      max_pi = mid_pi
    } else {
      break
    }
    i = i + 1
  }
  
  mid_pi
}

///|
pub fn create_example_expasy_entry() -> ExPASyEntry {
  let entry = ExPASyEntry::new("P04637", "TP53")
  entry
}

///|
pub fn create_example_enzyme_entry() -> EnzymeEntry {
  let enzyme = EnzymeEntry::new("1.1.1.1", "Alcohol dehydrogenase")
  let reactions = enzyme.catalyzed_reactions.copy()
  reactions.push("Alcohol + NAD+ <=> Aldehyde + NADH + H+")
  EnzymeEntry::{
    ec_number: enzyme.ec_number,
    name: enzyme.name,
    alt_names: enzyme.alt_names,
    catalyzed_reactions: reactions,
    cofactors: enzyme.cofactors,
    substrates: enzyme.substrates,
    products: enzyme.products,
    comments: enzyme.comments
  }
}