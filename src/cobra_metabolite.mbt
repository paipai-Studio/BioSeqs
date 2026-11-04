/// cobra_metabolite.mbt — Enhanced Metabolite (COBRApy core.Metabolite)
/// Metabolite with formula, charge, compartment, molecular weight.

///|
/// Enhanced metabolite with chemical properties.
pub struct CobraMetabolite {
  id : String
  name : String
  formula : String
  charge : Double
  compartment : String
}

///|
/// Parse a chemical formula string into element->count map.
/// e.g. "C6H12O6" -> {C: 6, H: 12, O: 6}
pub fn cobra_parse_formula(formula : String) -> Map[String, Int32] {
  let elements : Map[String, Int32] = Map([], capacity=0)
  let mut i = 0
  let s = formula
  while i < s.length() {
    let ch = s[i]
    // Element starts with uppercase letter
    if ch >= 'A' && ch <= 'Z' {
      let mut elem = ch.to_string()
      i = i + 1
      // Optional lowercase letter
      if i < s.length() && s[i] >= 'a' && s[i] <= 'z' {
        elem = elem + s[i].to_string()
        i = i + 1
      }
      // Optional count (digits)
      let mut count = 0
      let mut has_count = false
      while i < s.length() && s[i] >= '0' && s[i] <= '9' {
        count = count * 10 + (s[i].to_int() - '0'.to_int())
        has_count = true
        i = i + 1
      }
      if !has_count {
        count = 1
      }
      let prev = elements.get(elem).unwrap_or(0)
      elements.set(elem, prev + count)
    } else {
      i = i + 1
    }
  }
  elements
}

///|
/// Get element dictionary from a metabolite.
pub fn cobra_metabolite_elements(met : CobraMetabolite) -> Map[String, Int32] {
  cobra_parse_formula(met.formula)
}

///|
/// Compute molecular weight from formula.
/// Uses approximate atomic weights for common elements.
pub fn cobra_metabolite_formula_weight(met : CobraMetabolite) -> Double {
  let elements = cobra_parse_formula(met.formula)
  // Atomic weights (approximate)
  let weights : Map[String, Double] = Map(
    [
      ("H", 1.008),
      ("C", 12.011),
      ("N", 14.007),
      ("O", 15.999),
      ("P", 30.974),
      ("S", 32.065),
      ("Fe", 55.845),
      ("Mg", 24.305),
      ("Ca", 40.078),
      ("Cl", 35.453),
      ("K", 39.098),
      ("Na", 22.990),
      ("Zn", 65.38),
      ("Cu", 63.546),
      ("Mn", 54.938),
      ("Co", 58.933),
      ("Mo", 95.96),
      ("Ni", 58.693),
      ("Se", 78.971),
      ("F", 18.998),
      ("Br", 79.904),
      ("I", 126.904),
      ("B", 10.81),
      ("Si", 28.085),
      ("Ar", 39.948),
      ("He", 4.003),
      ("Li", 6.94),
      ("Be", 9.012),
      ("Al", 26.982),
      ("Ti", 47.867),
      ("V", 50.942),
      ("Cr", 51.996),
      ("Rb", 85.468),
      ("Sr", 87.62),
      ("Ag", 107.87),
      ("Cd", 112.41),
      ("Sn", 118.71),
      ("Ba", 137.33),
      ("W", 183.84),
      ("Pt", 195.08),
      ("Au", 196.97),
      ("Hg", 200.59),
      ("Pb", 207.2),
      ("U", 238.03),
    ],
    capacity=0,
  )
  let mut weight = 0.0
  for elem in elements.keys() {
    let count = elements.get(elem).unwrap_or(0)
    let w = weights.get(elem).unwrap_or(0.0)
    weight = weight + count.to_double() * w
  }
  weight
}

///|
/// Create a metabolite.
pub fn CobraMetabolite::new(
  id : String,
  formula~ : String = "",
  name~ : String = "",
  charge~ : Double = 0.0,
  compartment~ : String = "",
) -> CobraMetabolite {
  { id, name, formula, charge, compartment, }
}

///|
/// Check if two metabolites have the same formula and charge.
pub fn cobra_metabolite_compatible(
  a : CobraMetabolite,
  b : CobraMetabolite,
) -> Bool {
  a.formula == b.formula && a.charge == b.charge
}
