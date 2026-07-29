///|
/// Proteomics module - implements proteolytic cleavage agents, peptide mass calculation,
/// and isotope patterns. Inspired by Biopython Bio.SeqUtils.Proteomics.

///|
pub struct CleavageSite {
  position : Int
  residue : String
  probability : Double
}

///|
pub struct CleavageAgent {
  name : String
  cleavage_sites : Array[CleavageSite]
}

///|
pub struct PeptideMassResult {
  peptide_sequence : String
  monoisotopic_mass : Double
  average_mass : Double
  start_position : Int
  end_position : Int
}

///|
pub struct IsotopePattern {
  masses : Array[Double]
  abundances : Array[Double]
  charge : Int
}

///|
pub struct DigestResult {
  peptides : Array[PeptideMassResult]
  protein_mass : Double
  num_peptides : Int
}

///|
let average_weights : Map[String, Double] = Map([
  ("A", 89.0935),
  ("R", 174.2017),
  ("N", 132.1184),
  ("D", 133.1032),
  ("C", 121.1590),
  ("E", 147.1299),
  ("Q", 146.1451),
  ("G", 75.0669),
  ("H", 155.1552),
  ("I", 131.1736),
  ("L", 131.1736),
  ("K", 146.1882),
  ("M", 149.2124),
  ("F", 165.1900),
  ("P", 115.1310),
  ("S", 105.0930),
  ("T", 119.1197),
  ("W", 204.2262),
  ("Y", 181.1894),
  ("V", 117.1469),
], capacity=20)

///|
let monoisotopic_weights : Map[String, Double] = Map([
  ("A", 71.03711),
  ("R", 156.10111),
  ("N", 114.04293),
  ("D", 115.02694),
  ("C", 103.00919),
  ("E", 129.04259),
  ("Q", 128.05858),
  ("G", 57.02146),
  ("H", 137.05891),
  ("I", 113.08406),
  ("L", 113.08406),
  ("K", 128.09496),
  ("M", 131.04049),
  ("F", 147.06841),
  ("P", 97.05276),
  ("S", 87.03203),
  ("T", 101.04768),
  ("W", 186.07931),
  ("Y", 163.06333),
  ("V", 99.06841),
], capacity=20)

///|
let water_mono : Double = 18.010565
let water_avg : Double = 18.01524

///|
let proton_mass : Double = 1.007276

///|
let enzyme_data : Map[String, (String, Bool, String)] = Map([
  ("trypsin", ("KR", false, "P")),
  ("chymotrypsin", ("FYWML", false, "P")),
  ("chymotrypsin_high_specificity", ("FYW", false, "P")),
  ("chymotrypsin_low_specificity", ("FYWML", false, "P")),
  ("pepsin_ph2", ("FL", false, "")),
  ("pepsin_ph1.3", ("FL", false, "")),
  ("lysc", ("K", false, "P")),
  ("argc", ("R", false, "P")),
  ("cnbr", ("M", false, "")),
  ("gluc", ("DE", false, "")),
  ("aspn", ("D", true, "")),
  ("aspn_nterm", ("D", true, "")),
], capacity=15)

///|
pub fn proteomics_select_enzyme(name : String) -> CleavageAgent? {
  let key = name.to_lower()
  match enzyme_data.get(key) {
    Some((residues, _, _)) => {
      let sites : Array[CleavageSite] = Array::new()
      let mut i = 0
      while i < residues.length() {
        let r = residues[i:i + 1].to_owned()
        sites.push(CleavageSite::{
          position: 0,
          residue: r,
          probability: 1.0,
        })
        i = i + 1
      }
      Some(CleavageAgent::{
        name: key,
        cleavage_sites: sites,
      })
    }
    None => None
  }
}

///|
pub fn proteomics_available_enzymes() -> Array[String] {
  enzyme_data.keys().collect()
}

///|
pub fn proteomics_cleavage_sites(
  protein_sequence : String,
  enzyme : String,
) -> Array[Int] {
  let sites : Array[Int] = Array::new()
  let seq = protein_sequence.to_upper()
  let key = enzyme.to_lower()

  match enzyme_data.get(key) {
    Some((residues, before, except)) => {
      let mut i = 0
      while i < seq.length() {
        let c = seq[i:i + 1].to_owned()
        let is_cleavage_residue = proteomics_contains_char(residues, c)
        if is_cleavage_residue {
          if before {
            if i > 0 {
              let prev = seq[i - 1:i].to_owned()
              let is_except = proteomics_contains_char(except, prev)
              if !is_except {
                sites.push(i)
              }
            }
          } else {
            if i < seq.length() - 1 {
              let next = seq[i + 1:i + 2].to_owned()
              let is_except = proteomics_contains_char(except, next)
              if !is_except {
                sites.push(i + 1)
              }
            }
          }
        }
        i = i + 1
      }
    }
    None => ()
  }

  sites
}

///|
pub fn proteomics_digest(protein_sequence : String, enzyme : String) -> DigestResult {
  let seq = protein_sequence.to_upper()
  let sites = proteomics_cleavage_sites(seq, enzyme)

  let sorted_sites = proteomics_sort_ints(sites)

  let peptides : Array[PeptideMassResult] = Array::new()

  let starts : Array[Int] = Array::new()
  starts.push(0)
  for site in sorted_sites {
    starts.push(site)
  }

  let ends : Array[Int] = Array::new()
  for site in sorted_sites {
    ends.push(site)
  }
  ends.push(seq.length())

  let total_mass = proteomics_molecular_weight(seq)

  let mut i = 0
  while i < starts.length() {
    let s = starts[i]
    let e = ends[i]
    if s < e && e <= seq.length() {
      let pep_seq = seq[s:e].to_owned()
      let mono = proteomics_calculate_mass(pep_seq, true)
      let avg = proteomics_calculate_mass(pep_seq, false)
      peptides.push(PeptideMassResult::{
        peptide_sequence: pep_seq,
        monoisotopic_mass: mono,
        average_mass: avg,
        start_position: s,
        end_position: e,
      })
    }
    i = i + 1
  }

  DigestResult::{
    peptides,
    protein_mass: total_mass,
    num_peptides: peptides.length(),
  }
}

///|
pub fn proteomics_trypsin_digest(
  protein : String,
  missed_cleavages : Int,
) -> DigestResult {
  if missed_cleavages <= 0 {
    return proteomics_digest(protein, "trypsin")
  }

  let seq = protein.to_upper()
  let sites = proteomics_cleavage_sites(seq, "trypsin")
  let sorted_sites = proteomics_sort_ints(sites)

  let peptides : Array[PeptideMassResult] = Array::new()

  let starts : Array[Int] = Array::new()
  starts.push(0)
  for site in sorted_sites {
    starts.push(site)
  }

  let ends : Array[Int] = Array::new()
  for site in sorted_sites {
    ends.push(site)
  }
  ends.push(seq.length())

  let total_mass = proteomics_molecular_weight(seq)

  let mut i = 0
  while i < starts.length() {
    let s = starts[i]
    let e = ends[i]
    if s < e && e <= seq.length() {
      let pep_seq = seq[s:e].to_owned()
      let mono = proteomics_calculate_mass(pep_seq, true)
      let avg = proteomics_calculate_mass(pep_seq, false)
      peptides.push(PeptideMassResult::{
        peptide_sequence: pep_seq,
        monoisotopic_mass: mono,
        average_mass: avg,
        start_position: s,
        end_position: e,
      })
    }
    i = i + 1
  }

  let n_peptides = starts.length()
  let max_missed = if missed_cleavages < n_peptides - 1 {
    missed_cleavages
  } else {
    if n_peptides > 1 { n_peptides - 1 } else { 0 }
  }

  let mut mc = 1
  while mc <= max_missed {
    let mut idx = 0
    while idx < starts.length() - mc {
      let s = starts[idx]
      let e = ends[idx + mc]
      if s < e && e <= seq.length() {
        let pep_seq = seq[s:e].to_owned()
        let mono = proteomics_calculate_mass(pep_seq, true)
        let avg = proteomics_calculate_mass(pep_seq, false)
        peptides.push(PeptideMassResult::{
          peptide_sequence: pep_seq,
          monoisotopic_mass: mono,
          average_mass: avg,
          start_position: s,
          end_position: e,
        })
      }
      idx = idx + 1
    }
    mc = mc + 1
  }

  DigestResult::{
    peptides,
    protein_mass: total_mass,
    num_peptides: peptides.length(),
  }
}

///|
pub fn proteomics_calculate_mass(peptide : String, is_monoisotopic : Bool) -> Double {
  let seq = peptide.to_upper()
  if seq.length() == 0 {
    return 0.0
  }

  let weights = if is_monoisotopic {
    monoisotopic_weights
  } else {
    average_weights
  }

  let water = if is_monoisotopic { water_mono } else { water_avg }

  let mut mass = 0.0
  let mut i = 0
  while i < seq.length() {
    let c = seq[i:i + 1].to_owned()
    match weights.get(c) {
      Some(w) => mass = mass + w
      None => ()
    }
    i = i + 1
  }

  mass + water
}

///|
pub fn proteomics_molecular_weight(peptide : String) -> Double {
  proteomics_calculate_mass(peptide, false)
}

///|
pub fn proteomics_isotope_pattern(peptide : String, charge : Int) -> IsotopePattern {
  let seq = peptide.to_upper()
  let n = seq.length()

  let num_peaks = if n < 5 { 3 } else { if n < 15 { 4 } else { if n < 30 { 5 } else { 6 } } }

  let mono_mass = proteomics_calculate_mass(seq, true)
  let avg_mass = proteomics_calculate_mass(seq, false)

  let masses : Array[Double] = Array::new()
  let abundances : Array[Double] = Array::new()

  let delta_mass = if n > 0 { (avg_mass - mono_mass) / n.to_double() } else { 1.0 }

  let mut i = 0
  while i < num_peaks {
    let mass = mono_mass + (delta_mass * i.to_double())
    let m = if charge > 0 {
      (mass / charge.to_double()) + proton_mass
    } else {
      mass
    }
    masses.push(m)
    let abundance = proteomics_isotope_abundance(i, n)
    abundances.push(abundance)
    i = i + 1
  }

  let max_ab = abundances.get(0).unwrap_or(1.0)
  let normalized : Array[Double] = Array::new()
  let mut j = 0
  while j < abundances.length() {
    normalized.push(abundances[j] / max_ab)
    j = j + 1
  }

  IsotopePattern::{
    masses,
    abundances: normalized,
    charge,
  }
}

///|
fn proteomics_isotope_abundance(peak_index : Int, peptide_length : Int) -> Double {
  if peptide_length <= 0 {
    return if peak_index == 0 { 1.0 } else { 0.0 }
  }

  let n = peptide_length.to_double()

  let lambda = n * 0.0107
  let log_p = -lambda + peak_index.to_double() * proteomics_log_approx(lambda) - proteomics_log_factorial(peak_index)
  if log_p < -50.0 {
    return 0.0
  }
  proteomics_exp_approx(log_p)
}

///|
fn proteomics_log_factorial(n : Int) -> Double {
  if n <= 1 { return 0.0 }
  let mut result = 0.0
  let mut i = 2
  while i <= n {
    result = result + proteomics_log_approx(i.to_double())
    i = i + 1
  }
  result
}

///|
fn proteomics_log_approx(x : Double) -> Double {
  if x <= 0.0 {
    return -1000000000000000000000000000000000.0
  }
  let y = (x - 1.0) / (x + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i < 50 {
    term = term * y2
    sum = sum + term / (2 * i + 1).to_double()
    i = i + 1
  }
  2.0 * sum
}

///|
fn proteomics_exp_approx(x : Double) -> Double {
  let term = 1.0 + x / 25.0
  proteomics_pow_approx(term, 25)
}

///|
fn proteomics_pow_approx(base : Double, exp : Int) -> Double {
  if exp <= 0 { return 1.0 }
  let mut result = 1.0
  let mut b = base
  let mut e = exp
  while e > 0 {
    if (e & 1) != 0 {
      result = result * b
    }
    e = e >> 1
    b = b * b
  }
  result
}

///|
pub fn proteomics_fragment_ions(
  peptide : String,
  charge : Int,
) -> Array[PeptideMassResult] {
  let seq = peptide.to_upper()
  let frags : Array[PeptideMassResult] = Array::new()

  let n = seq.length()
  if n == 0 {
    return frags
  }

  let mut i = 1
  while i < n {
    let b_seq = seq[0:i].to_owned()
    let y_seq = seq[i:n].to_owned()

    let b_mono = proteomics_calculate_mass(b_seq, true)
    let b_avg = proteomics_calculate_mass(b_seq, false)

    let y_mono = proteomics_calculate_mass(y_seq, true)
    let y_avg = proteomics_calculate_mass(y_seq, false)

    let b_z = if charge > 0 { charge } else { 1 }
    let b_ion_mono = (b_mono + (b_z.to_double() - 1.0) * proton_mass) / b_z.to_double()
    let b_ion_avg = (b_avg + (b_z.to_double() - 1.0) * proton_mass) / b_z.to_double()

    frags.push(PeptideMassResult::{
      peptide_sequence: "b" + i.to_string(),
      monoisotopic_mass: b_ion_mono,
      average_mass: b_ion_avg,
      start_position: 0,
      end_position: i,
    })

    let y_ion_mono = (y_mono + (b_z.to_double() - 1.0) * proton_mass) / b_z.to_double()
    let y_ion_avg = (y_avg + (b_z.to_double() - 1.0) * proton_mass) / b_z.to_double()

    frags.push(PeptideMassResult::{
      peptide_sequence: "y" + (n - i).to_string(),
      monoisotopic_mass: y_ion_mono,
      average_mass: y_ion_avg,
      start_position: i,
      end_position: n,
    })

    i = i + 1
  }

  frags
}

///|
fn proteomics_contains_char(s : String, c : String) -> Bool {
  let mut i = 0
  while i < s.length() {
    if s[i:i + 1] == c {
      return true
    }
    i = i + 1
  }
  false
}

///|
fn proteomics_sort_ints(arr : Array[Int]) -> Array[Int] {
  let result = arr
  ignore(proteomics_quick_sort(result, 0, result.length() - 1))
  result
}

///|
fn proteomics_quick_sort(arr : Array[Int], lo : Int, hi : Int) -> Int {
  if lo < hi {
    let p = proteomics_partition(arr, lo, hi)
    ignore(proteomics_quick_sort(arr, lo, p - 1))
    ignore(proteomics_quick_sort(arr, p + 1, hi))
  }
  0
}

///|
fn proteomics_partition(arr : Array[Int], lo : Int, hi : Int) -> Int {
  let pivot = arr[hi]
  let mut i = lo - 1
  let mut j = lo
  while j < hi {
    if arr[j] <= pivot {
      i = i + 1
      let tmp = arr[i]
      arr[i] = arr[j]
      arr[j] = tmp
    }
    j = j + 1
  }
  let tmp = arr[i + 1]
  arr[i + 1] = arr[hi]
  arr[hi] = tmp
  i + 1
}