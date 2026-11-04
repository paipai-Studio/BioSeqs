/// PAINS (Pan Assay Interference Compounds) filters.
///
/// RDKit equivalent: `rdkit.Chem.FilterCatalog` with PAINS catalog.
///
/// PAINS compounds are frequent hitters in biochemical assays due to
/// reactive or promiscuous substructures. This module checks molecules
/// against a curated set of alert substructures using SMARTS matching.

///|
/// A single PAINS alert
pub struct RdkitPainsAlert {
  name : String
  smarts : String
  matched : Bool
}

///|
pub fn RdkitPainsAlert::name(self : RdkitPainsAlert) -> String {
  self.name
}

///|
pub fn RdkitPainsAlert::smarts(self : RdkitPainsAlert) -> String {
  self.smarts
}

///|
pub fn RdkitPainsAlert::matched(self : RdkitPainsAlert) -> Bool {
  self.matched
}

///|
/// PAINS alert substructures (simplified set of common offenders).
/// Each entry: (name, SMARTS pattern)
let rdkit_pains_alerts : Array[(String, String)] = [
  // Rhodanine — known PAINS
  ("Rhodanine", "C1(=S)NC(=S)CC1"),
  // Enone (Michael acceptor)
  ("Enone_Michael_acceptor", "C=CC(=O)"),
  // Catechol — oxidizes easily
  ("Catechol", "c1cc(O)c(O)cc1"),
  // Hydroquinone
  ("Hydroquinone", "c1cc(O)cc(O)c1"),
  // Quinone
  ("Quinone", "C1(=O)C=CC(=O)C=C1"),
  // Azo group
  ("Azo", "N=N"),
  // Alkyl halide — reactive
  ("Alkyl_halide", "C-[F,Cl,Br,I]"),
  // Thiol — forms disulfides
  ("Thiol", "[SH]"),
  // Aldehyde — reactive
  ("Aldehyde", "[CX1]=O"),
  // Anhydride
  ("Anhydride", "C(=O)OC(=O)"),
  // Beta-lactam
  ("Beta_lactam", "C1(=O)NC1"),
  // Furan — unstable
  ("Furan", "c1ccoc1"),
  // Thiophene S-oxide
  ("Thiophene", "c1ccsc1"),
  // Peroxide
  ("Peroxide", "COO"),
  // Nitroso
  ("Nitroso", "N=O"),
  // Iminium ion
  ("Iminium", "C=[N+]"),
]

///|
/// Check a molecule against all PAINS alerts.
/// Returns list of alerts that matched.
pub fn rdkit_pains_check(mol : Molecule) -> Array[RdkitPainsAlert] {
  let results : Array[RdkitPainsAlert] = Array(capacity=0)
  for alert in rdkit_pains_alerts {
    let matched = smarts_match_str(mol, alert.1)
    if matched {
      results.push({ name: alert.0, smarts: alert.1, matched: true, })
    }
  }
  results
}

///|
/// Quick check: does molecule contain any PAINS substructure?
pub fn rdkit_has_pains(mol : Molecule) -> Bool {
  for alert in rdkit_pains_alerts {
    if smarts_match_str(mol, alert.1) {
      return true
    }
  }
  false
}

///|
/// Count number of PAINS alerts triggered.
pub fn rdkit_pains_count(mol : Molecule) -> Int32 {
  let mut count = 0
  for alert in rdkit_pains_alerts {
    if smarts_match_str(mol, alert.1) {
      count = count + 1
    }
  }
  count
}

///|
/// Get names of all triggered PAINS alerts.
pub fn rdkit_pains_names(mol : Molecule) -> Array[String] {
  let names : Array[String] = Array(capacity=0)
  for alert in rdkit_pains_alerts {
    if smarts_match_str(mol, alert.1) {
      names.push(alert.0)
    }
  }
  names
}
