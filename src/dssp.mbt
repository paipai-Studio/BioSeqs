///|
/// Bio.PDB.DSSP module, inspired by Biopython Bio.PDB.DSSP.
/// Provides functions for secondary structure analysis using DSSP algorithm.

///|
pub struct DSSPRecord {
  chain_id : String
  residue_number : Int32
  insertion_code : String
  amino_acid : String
  secondary_structure : String
  solvent_accessibility : Int32
  phi : Double
  psi : Double
}

///|
pub fn DSSPRecord::new(
  chain_id : String,
  residue_number : Int32,
  insertion_code : String,
  amino_acid : String,
  secondary_structure : String,
  solvent_accessibility : Int32,
  phi : Double,
  psi : Double,
) -> DSSPRecord {
  DSSPRecord::{
    chain_id,
    residue_number,
    insertion_code,
    amino_acid,
    secondary_structure,
    solvent_accessibility,
    phi,
    psi,
  }
}

///|
pub struct DSSPResult {
  records : Array[DSSPRecord]
}

///|
pub fn DSSPResult::new(records : Array[DSSPRecord]) -> DSSPResult {
  DSSPResult::{ records, }
}

///|
pub fn predict_secondary_structure(sequence : String) -> String {
  let mut result = ""

  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let ss = classify_secondary_structure(aa)
    result = result + ss
    i = i + 1
  }

  result
}

///|
fn classify_secondary_structure(amino_acid : String) -> String {
  let helix_formers = ["E", "D", "L", "M", "K", "R", "H"]
  let sheet_formers = ["V", "I", "Y", "F", "W", "T"]
  let coil_formers = ["G", "A", "S", "P", "N", "C", "Q"]

  let mut i = 0
  while i < helix_formers.length() {
    if helix_formers[i] == amino_acid {
      return "H"
    }
    i = i + 1
  }

  i = 0
  while i < sheet_formers.length() {
    if sheet_formers[i] == amino_acid {
      return "E"
    }
    i = i + 1
  }

  i = 0
  while i < coil_formers.length() {
    if coil_formers[i] == amino_acid {
      return "-"
    }
    i = i + 1
  }

  "-"
}

///|
pub fn calculate_solvent_accessibility(
  sequence : String,
  structure : String,
) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)

  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let ss = if i < structure.length() { char_at(structure, i) } else { "-" }

    let accessibility = estimate_accessibility(aa, ss, i, sequence.length())
    result.push(accessibility)

    i = i + 1
  }

  result
}

///|
fn estimate_accessibility(
  amino_acid : String,
  secondary_structure : String,
  position : Int32,
  total_length : Int32,
) -> Int32 {
  let max_accessibility = Map(
    [
      ("A", 129),
      ("R", 274),
      ("N", 195),
      ("D", 193),
      ("C", 167),
      ("Q", 225),
      ("E", 223),
      ("G", 104),
      ("H", 224),
      ("I", 197),
      ("L", 201),
      ("K", 236),
      ("M", 224),
      ("F", 240),
      ("P", 159),
      ("S", 155),
      ("T", 172),
      ("W", 285),
      ("Y", 263),
      ("V", 174),
    ],
    capacity=20,
  )

  let max_acc = max_accessibility.get(amino_acid).unwrap_or(150)

  let mut accessibility = max_acc

  if secondary_structure == "H" {
    accessibility = (max_acc.to_double() * 0.3).to_int()
  } else if secondary_structure == "E" {
    accessibility = (max_acc.to_double() * 0.4).to_int()
  } else {
    accessibility = (max_acc.to_double() * 0.6).to_int()
  }

  if position < 5 || position > total_length - 5 {
    accessibility = (accessibility.to_double() * 0.5).to_int()
  }

  accessibility
}

///|
pub fn calculate_ramachandran_plot(
  records : Array[DSSPRecord],
) -> Array[(Double, Double)] {
  let result : Array[(Double, Double)] = Array(capacity=0)

  let mut i = 0
  while i < records.length() {
    result.push((records[i].phi, records[i].psi))
    i = i + 1
  }

  result
}

///|
pub fn analyze_structure_composition(
  records : Array[DSSPRecord],
) -> Map[String, Int32] {
  let composition : Map[String, Int32] = Map(
    [
      ("H", 0),
      ("B", 0),
      ("E", 0),
      ("G", 0),
      ("I", 0),
      ("T", 0),
      ("S", 0),
      ("-", 0),
    ],
    capacity=8,
  )

  let mut i = 0
  while i < records.length() {
    let ss = records[i].secondary_structure
    let count = composition.get(ss).unwrap_or(0)
    composition.set(ss, count + 1)
    i = i + 1
  }

  composition
}

///|
pub fn create_example_dssp_data() -> DSSPResult {
  let records : Array[DSSPRecord] = Array(capacity=0)

  records.push(DSSPRecord::new("A", 1, "", "M", "-", 100, -60.0, -40.0))
  records.push(DSSPRecord::new("A", 2, "", "S", "H", 50, -57.0, -47.0))
  records.push(DSSPRecord::new("A", 3, "", "G", "H", 40, -59.0, -47.0))
  records.push(DSSPRecord::new("A", 4, "", "L", "H", 35, -57.0, -47.0))
  records.push(DSSPRecord::new("A", 5, "", "P", "-", 80, -60.0, 120.0))
  records.push(DSSPRecord::new("A", 6, "", "V", "E", 30, -139.0, 135.0))
  records.push(DSSPRecord::new("A", 7, "", "A", "E", 25, -139.0, 135.0))
  records.push(DSSPRecord::new("A", 8, "", "L", "-", 60, -60.0, -40.0))

  DSSPResult::new(records)
}
