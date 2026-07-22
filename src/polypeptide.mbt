///|
/// Bio.PDB.Polypeptide module, inspired by Biopython Bio.PDB.Polypeptide.
/// Provides functions for polypeptide analysis, amino acid composition, and hydrophobicity.

///|
pub struct PolypeptideComposition {
  amino_acids : Map[String, Int]
  total_residues : Int
}

///|
pub fn PolypeptideComposition::new(amino_acids : Map[String, Int], total_residues : Int) -> PolypeptideComposition {
  PolypeptideComposition::{ amino_acids, total_residues }
}

///|
pub struct HydrophobicityProfile {
  positions : Array[Int]
  values : Array[Double]
}

///|
pub fn HydrophobicityProfile::new(positions : Array[Int], values : Array[Double]) -> HydrophobicityProfile {
  HydrophobicityProfile::{ positions, values }
}

///|
pub struct TransmembraneRegion {
  start : Int
  end : Int
  score : Double
}

///|
pub fn TransmembraneRegion::new(start : Int, end : Int, score : Double) -> TransmembraneRegion {
  TransmembraneRegion::{ start, end, score }
}

///|
pub fn calculate_composition(sequence : String) -> PolypeptideComposition {
  let amino_acids : Map[String, Int] = Map([
    ("A", 0), ("R", 0), ("N", 0), ("D", 0), ("C", 0),
    ("Q", 0), ("E", 0), ("G", 0), ("H", 0), ("I", 0),
    ("L", 0), ("K", 0), ("M", 0), ("F", 0), ("P", 0),
    ("S", 0), ("T", 0), ("W", 0), ("Y", 0), ("V", 0),
  ], capacity=20)
  
  let mut total_residues = 0
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let count = amino_acids.get(aa).unwrap_or(0)
    amino_acids.set(aa, count + 1)
    total_residues = total_residues + 1
    i = i + 1
  }
  
  PolypeptideComposition::new(amino_acids, total_residues)
}

///|
pub fn calculate_percent_composition(composition : PolypeptideComposition) -> Map[String, Double] {
  let percent : Map[String, Double] = Map([], capacity=20)
  
  if composition.total_residues == 0 {
    return percent
  }
  
  let amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  
  let mut i = 0
  while i < amino_acids.length() {
    let aa = amino_acids[i]
    let count = composition.amino_acids.get(aa).unwrap_or(0)
    let percentage = count.to_double() / composition.total_residues.to_double() * 100.0
    percent.set(aa, percentage)
    i = i + 1
  }
  
  percent
}

///|
pub fn calculate_hydrophobicity(sequence : String, scale : String) -> HydrophobicityProfile {
  let positions : Array[Int] = Array::new()
  let values : Array[Double] = Array::new()
  
  let hydrophobicity_map = if scale == "kyte-doolittle" {
    kyte_doolittle_scale()
  } else if scale == "hopp-woods" {
    hopp_woods_scale()
  } else {
    kyte_doolittle_scale()
  }
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    let score = hydrophobicity_map.get(aa).unwrap_or(0.0)
    
    positions.push(i + 1)
    values.push(score)
    
    i = i + 1
  }
  
  HydrophobicityProfile::new(positions, values)
}

///|
fn kyte_doolittle_scale() -> Map[String, Double] {
  Map([
    ("A", 1.8), ("R", -4.5), ("N", -3.5), ("D", -3.5), ("C", 2.5),
    ("Q", -3.5), ("E", -3.5), ("G", -0.4), ("H", -3.2), ("I", 4.5),
    ("L", 3.8), ("K", -3.9), ("M", 1.9), ("F", 2.8), ("P", -1.6),
    ("S", -0.8), ("T", -0.7), ("W", -0.9), ("Y", -1.3), ("V", 4.2),
  ], capacity=20)
}

///|
fn hopp_woods_scale() -> Map[String, Double] {
  Map([
    ("A", -0.5), ("R", 3.0), ("N", 0.2), ("D", 3.0), ("C", -1.0),
    ("Q", 0.2), ("E", 3.0), ("G", 0.0), ("H", -0.5), ("I", -1.8),
    ("L", -1.8), ("K", 3.0), ("M", -1.3), ("F", -2.5), ("P", 0.0),
    ("S", 0.3), ("T", -0.4), ("W", -3.4), ("Y", -2.3), ("V", -1.5),
  ], capacity=20)
}

///|
pub fn calculate_hydrophobicity_window(
  sequence : String,
  window_size : Int,
  scale : String,
) -> HydrophobicityProfile {
  let positions : Array[Int] = Array::new()
  let values : Array[Double] = Array::new()
  
  let hydrophobicity_map = if scale == "kyte-doolittle" {
    kyte_doolittle_scale()
  } else if scale == "hopp-woods" {
    hopp_woods_scale()
  } else {
    kyte_doolittle_scale()
  }
  
  let half_window = window_size / 2
  
  let mut i = 0
  while i < sequence.length() {
    let mut sum = 0.0
    let mut count = 0
    
    let start = if i < half_window { 0 } else { i - half_window }
    let end = if i + half_window >= sequence.length() { sequence.length() } else { i + half_window + 1 }
    
    let mut j = start
    while j < end {
      let aa = char_at(sequence, j)
      sum = sum + hydrophobicity_map.get(aa).unwrap_or(0.0)
      count = count + 1
      j = j + 1
    }
    
    let avg = sum / count.to_double()
    
    positions.push(i + 1)
    values.push(avg)
    
    i = i + 1
  }
  
  HydrophobicityProfile::new(positions, values)
}

///|
pub fn predict_transmembrane_regions(
  sequence : String,
  window_size : Int,
  threshold : Double,
) -> Array[TransmembraneRegion] {
  let regions : Array[TransmembraneRegion] = Array::new()
  
  let profile = calculate_hydrophobicity_window(sequence, window_size, "kyte-doolittle")
  
  let mut in_region = false
  let mut region_start = 0
  let mut region_scores : Array[Double] = Array::new()
  
  let mut i = 0
  while i < profile.values.length() {
    let score = profile.values[i]
    
    if score >= threshold && !in_region {
      in_region = true
      region_start = profile.positions[i]
      region_scores = Array::new()
      region_scores.push(score)
    } else if score >= threshold && in_region {
      region_scores.push(score)
    } else if score < threshold && in_region {
      in_region = false
      
      if region_scores.length() >= 5 {
        let mut avg_score = 0.0
        let mut j = 0
        while j < region_scores.length() {
          avg_score = avg_score + region_scores[j]
          j = j + 1
        }
        avg_score = avg_score / region_scores.length().to_double()
        
        regions.push(TransmembraneRegion::new(region_start, profile.positions[i - 1], avg_score))
      }
    }
    
    i = i + 1
  }
  
  if in_region && region_scores.length() >= 5 {
    let mut avg_score = 0.0
    let mut j = 0
    while j < region_scores.length() {
      avg_score = avg_score + region_scores[j]
      j = j + 1
    }
    avg_score = avg_score / region_scores.length().to_double()
    
    regions.push(TransmembraneRegion::new(region_start, profile.positions[profile.positions.length() - 1], avg_score))
  }
  
  regions
}

///|
pub fn calculate_isoelectric_point(sequence : String) -> Double {
  let mut low = 0.0
  let mut high = 14.0
  let mut mid = 7.0
  
  let mut i = 0
  while i < 50 {
    mid = (low + high) / 2.0
    let charge = calculate_charge(sequence, mid)
    
    if charge > 0.0 {
      low = mid
    } else {
      high = mid
    }
    
    i = i + 1
  }
  
  mid
}

///|
fn calculate_charge(sequence : String, ph : Double) -> Double {
  let mut charge = 0.0
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    charge = charge + amino_acid_charge(aa, ph)
    i = i + 1
  }
  
  charge
}

///|
fn amino_acid_charge(amino_acid : String, ph : Double) -> Double {
  let pka_values = Map([
    ("A", (2.34, 9.69, 0.0)), ("R", (2.17, 9.04, 12.48)), ("N", (2.02, 8.80, 0.0)),
    ("D", (1.88, 9.60, 3.65)), ("C", (1.96, 10.28, 8.18)), ("Q", (2.17, 9.13, 0.0)),
    ("E", (2.19, 9.67, 4.25)), ("G", (2.34, 9.60, 0.0)), ("H", (1.82, 9.17, 6.00)),
    ("I", (2.36, 9.60, 0.0)), ("L", (2.36, 9.60, 0.0)), ("K", (2.18, 8.95, 10.53)),
    ("M", (2.28, 9.21, 0.0)), ("F", (1.83, 9.13, 0.0)), ("P", (1.99, 10.60, 0.0)),
    ("S", (2.21, 9.15, 0.0)), ("T", (2.11, 9.62, 0.0)), ("W", (2.83, 9.39, 0.0)),
    ("Y", (2.32, 9.62, 10.07)), ("V", (2.32, 9.62, 0.0)),
  ], capacity=20)
  
  let pka_tuple = pka_values.get(amino_acid).unwrap_or((2.0, 9.0, 0.0))
  let (pka_n, pka_c, pka_side) = pka_tuple
  
  let mut charge = 0.0
  
  charge = charge + @math.pow(10.0, pka_n) / (@math.pow(10.0, pka_n) + @math.pow(10.0, ph))
  charge = charge - @math.pow(10.0, ph) / (@math.pow(10.0, pka_c) + @math.pow(10.0, ph))
  
  if pka_side > 0.0 {
    if amino_acid == "D" || amino_acid == "E" {
      charge = charge - @math.pow(10.0, ph) / (@math.pow(10.0, pka_side) + @math.pow(10.0, ph))
    } else if amino_acid == "R" || amino_acid == "K" {
      charge = charge + @math.pow(10.0, pka_side) / (@math.pow(10.0, pka_side) + @math.pow(10.0, ph))
    } else if amino_acid == "H" {
      charge = charge + @math.pow(10.0, pka_side) / (@math.pow(10.0, pka_side) + @math.pow(10.0, ph))
    } else if amino_acid == "C" {
      charge = charge - @math.pow(10.0, ph) / (@math.pow(10.0, pka_side) + @math.pow(10.0, ph))
    } else if amino_acid == "Y" {
      charge = charge - @math.pow(10.0, ph) / (@math.pow(10.0, pka_side) + @math.pow(10.0, ph))
    }
  }
  
  charge
}

///|
pub fn calculate_aromaticity(sequence : String) -> Double {
  let aromatic = ["F", "Y", "W"]
  
  let mut count = 0
  let mut total = 0
  
  let mut i = 0
  while i < sequence.length() {
    let aa = char_at(sequence, i)
    total = total + 1
    
    let mut j = 0
    while j < aromatic.length() {
      if aromatic[j] == aa {
        count = count + 1
        break
      }
      j = j + 1
    }
    
    i = i + 1
  }
  
  if total == 0 { 0.0 } else { count.to_double() / total.to_double() }
}

///|
pub fn create_example_polypeptide_data() -> (String, PolypeptideComposition, HydrophobicityProfile) {
  let sequence = "MSGLVLLLFLLVSLTGVASSTPGGGEETLAALGRALEEELARHPYFYAPELLYYANKYNGVFQEAAELRKTVVLDVTSEQLRRLDTHLETYGKHEEHAQEAGERGPLRESLRALYNVPQLIESKGGVVPRRYKPVNDSGFFYVEELKDIKEKVAADDDDMVLKAVREAAKEGADFTMVTEIIEELKARLEKFNGEEDKARLAKDALKESDDWPLVKVLKNKDDYVHGFEEA";
  
  let composition = calculate_composition(sequence)
  let profile = calculate_hydrophobicity(sequence, "kyte-doolittle")
  
  (sequence, composition, profile)
}