/// rdkit_props_predict.mbt — toxicity alerts, solubility estimation, similarity map
/// Items 15-17

///|
/// Toxicity alert result
pub struct RdkitToxAlert {
  alert_name : String
  category : String // "mutagen", "hERG", "hepatotox", "skin_sens"
  severity : String // "high", "medium", "low"
  matched : Bool
  message : String
}

///|
/// Check for potential Ames mutagenicity using structural alerts.
pub fn rdkit_ames_alerts(smiles : String) -> Array[RdkitToxAlert] {
  let alerts : Array[RdkitToxAlert] = Array(capacity=0)
  // Aromatic amine: potentially mutagenic
  if rdkit_has_substruct(smiles, "cN") {
    alerts.push(RdkitToxAlert::{
      alert_name: "aromatic_amine",
      category: "mutagen",
      severity: "high",
      matched: true,
      message: "Aromatic amine detected — potential Ames positive",
    })
  }
  // Nitro group
  if rdkit_tox_has_nitro(smiles) {
    alerts.push(RdkitToxAlert::{
      alert_name: "nitro_group",
      category: "mutagen",
      severity: "high",
      matched: true,
      message: "Nitro group detected — known mutagenic alert",
    })
  }
  // Alkyl halide: potential electrophile
  if rdkit_tox_has_alkyl_halide(smiles) {
    alerts.push(RdkitToxAlert::{
      alert_name: "alkyl_halide",
      category: "mutagen",
      severity: "medium",
      matched: true,
      message: "Alkyl halide — potential electrophilic mutagen",
    })
  }
  // Aromatic nitro
  if rdkit_tox_has_aromatic_nitro(smiles) {
    alerts.push(RdkitToxAlert::{
      alert_name: "aromatic_nitro",
      category: "mutagen",
      severity: "high",
      matched: true,
      message: "Aromatic nitro — strong mutagenic alert",
    })
  }
  // If no alerts, add a "clean" result
  if alerts.length() == 0 {
    alerts.push(RdkitToxAlert::{
      alert_name: "none",
      category: "mutagen",
      severity: "low",
      matched: false,
      message: "No known mutagenic structural alerts",
    })
  }
  alerts
}

///|
fn rdkit_tox_has_nitro(smiles : String) -> Bool {
  smiles.contains("N(=O)=O") || smiles.contains("[N+](=O)[O-]")
}

///|
fn rdkit_tox_has_alkyl_halide(smiles : String) -> Bool {
  smiles.contains("CF") ||
  smiles.contains("CCl") ||
  smiles.contains("CBr") ||
  smiles.contains("CI")
}

///|
fn rdkit_tox_has_aromatic_nitro(smiles : String) -> Bool {
  smiles.contains("c[N+](=O)[O-]") || smiles.contains("cN(=O)=O")
}

///|
/// Check for potential hERG channel blockade.
pub fn rdkit_herg_alerts(smiles : String) -> Array[RdkitToxAlert] {
  let alerts : Array[RdkitToxAlert] = Array(capacity=0)
  let mol = parse_smiles(smiles)
  let _n = mol.num_atoms()
  let mw = rdkit_tox_mol_wt(mol)
  // hERG liability: MW > 400 and contains basic nitrogen
  let has_basic_n = rdkit_tox_has_basic_n(smiles)
  if mw > 400.0 && has_basic_n {
    alerts.push(RdkitToxAlert::{
      alert_name: "herg_basic_n_heavy",
      category: "hERG",
      severity: "high",
      matched: true,
      message: "High MW + basic N — high hERG liability risk",
    })
  } else if mw > 300.0 && has_basic_n {
    alerts.push(RdkitToxAlert::{
      alert_name: "herg_basic_n_moderate",
      category: "hERG",
      severity: "medium",
      matched: true,
      message: "Moderate MW + basic N — moderate hERG risk",
    })
  } else {
    alerts.push(RdkitToxAlert::{
      alert_name: "none",
      category: "hERG",
      severity: "low",
      matched: false,
      message: "Low hERG liability based on simple filters",
    })
  }
  alerts
}

///|
fn rdkit_tox_has_basic_n(smiles : String) -> Bool {
  // Simplified: contains N not in amide
  smiles.contains("N") && !smiles.contains("C(=O)N")
}

///|
fn rdkit_tox_mol_wt(mol : Molecule) -> Double {
  let symbols = mol.atom_symbols()
  let mut total = 0.0
  for i = 0; i < symbols.length(); i = i + 1 {
    total = total + rdkit_tox_atomic_wt(symbols[i])
  }
  total
}

///|
fn rdkit_tox_atomic_wt(sym : String) -> Double {
  match sym {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "P" => 30.974
    "S" => 32.06
    "Cl" => 35.45
    "Br" => 79.904
    "I" => 126.904
    _ => 0.0
  }
}

///|
/// Check for potential hepatotoxicity.
pub fn rdkit_hepatotox_alerts(smiles : String) -> Array[RdkitToxAlert] {
  let alerts : Array[RdkitToxAlert] = Array(capacity=0)
  // Bromine-containing compounds
  if smiles.contains("Br") {
    alerts.push(RdkitToxAlert::{
      alert_name: "bromine",
      category: "hepatotox",
      severity: "medium",
      matched: true,
      message: "Bromine-containing — potential hepatotoxicity",
    })
  }
  // Aniline
  if rdkit_has_substruct(smiles, "cN") {
    alerts.push(RdkitToxAlert::{
      alert_name: "aniline",
      category: "hepatotox",
      severity: "medium",
      matched: true,
      message: "Aniline moiety — potential hepatotoxicity",
    })
  }
  // Thiophene
  if smiles.contains("c1ccsc1") || smiles.contains("c1sccc1") {
    alerts.push(RdkitToxAlert::{
      alert_name: "thiophene",
      category: "hepatotox",
      severity: "medium",
      matched: true,
      message: "Thiophene ring — potential hepatotoxicity",
    })
  }
  if alerts.length() == 0 {
    alerts.push(RdkitToxAlert::{
      alert_name: "none",
      category: "hepatotox",
      severity: "low",
      matched: false,
      message: "No known hepatotoxicity alerts",
    })
  }
  alerts
}

///|
/// Combined toxicity assessment.
pub fn rdkit_tox_assessment(smiles : String) -> Array[RdkitToxAlert] {
  let all : Array[RdkitToxAlert] = Array(capacity=0)
  let ames = rdkit_ames_alerts(smiles)
  let herg = rdkit_herg_alerts(smiles)
  let hep = rdkit_hepatotox_alerts(smiles)
  for i = 0; i < ames.length(); i = i + 1 {
    all.push(ames[i])
  }
  for i = 0; i < herg.length(); i = i + 1 {
    all.push(herg[i])
  }
  for i = 0; i < hep.length(); i = i + 1 {
    all.push(hep[i])
  }
  all
}

///|
/// Solubility result
pub struct RdkitSolubility {
  esol : Double // Estimated SOLubility (log S)
  solubility_mol_per_l : Double // 10^logS
  class_ : String // "high", "moderate", "low", "very low"
}

///|
/// Estimate aqueous solubility using ESOL (Estimated SOLubility).
/// ESOL = 0.808 - 0.873*logP - 0.0198*MW + 0.0052*TPSA - 0.0056*N_rot
pub fn rdkit_esol(smiles : String) -> RdkitSolubility {
  let mol = parse_smiles(smiles)
  let mw = rdkit_tox_mol_wt(mol)
  let logp = rdkit_sol_logp(mol)
  let tpsa = rdkit_sol_tpsa(mol)
  let n_rot = rdkit_sol_rotatable(mol)
  let log_s = 0.808 -
    0.873 * logp -
    0.0198 * mw +
    0.0052 * tpsa -
    0.0056 * n_rot.to_double()
  let sol_mol = @math.pow(10.0, log_s)
  let class_ = if log_s > -2.0 {
    "high"
  } else if log_s > -4.0 {
    "moderate"
  } else if log_s > -6.0 {
    "low"
  } else {
    "very low"
  }
  RdkitSolubility::{ esol: log_s, solubility_mol_per_l: sol_mol, class_, }
}

///|
fn rdkit_sol_logp(mol : Molecule) -> Double {
  // Use existing descriptor if available
  mol_logp(mol)
}

///|
fn rdkit_sol_tpsa(mol : Molecule) -> Double {
  mol_tpsa(mol)
}

///|
fn rdkit_sol_rotatable(mol : Molecule) -> Int32 {
  mol_rotatable_bonds(mol)
}

///|
/// Similarity map: compute per-atom contribution to Tanimoto similarity.
pub fn rdkit_similarity_map(
  smiles1 : String,
  smiles2 : String,
  radius~ : Int32 = 2,
) -> Array[(Int32, Double)] {
  let mol1 = parse_smiles(smiles1)
  let n1 = mol1.num_atoms()
  let fp1 = rdkit_ecfp_fp(smiles1, radius~, n_bits=2048)
  let fp2 = rdkit_ecfp_fp(smiles2, radius~, n_bits=2048)
  let base_sim = rdkit_fp_tanimoto(fp1, fp2)
  // For each atom, remove its contribution and measure similarity change
  let contributions : Array[(Int32, Double)] = Array(capacity=n1)
  let symbols1 = mol1.atom_symbols()
  let adj1 = rdkit_build_adj_matrix(mol1)
  for i = 0; i < n1; i = i + 1 {
    // Recompute FP without atom i
    let fp_without = rdkit_simmap_fp_without(
      n1, symbols1, adj1, i, radius, 2048,
    )
    let sim_without = rdkit_fp_tanimoto(fp_without, fp2)
    let contribution = base_sim - sim_without
    contributions.push((i, contribution))
  }
  contributions
}

///|
fn rdkit_simmap_fp_without(
  n : Int32,
  symbols : Array[String],
  adj : Array[Array[Int32]],
  exclude : Int32,
  radius : Int32,
  n_bits : Int32,
) -> Array[Int32] {
  let fp = Array(capacity=n_bits)
  for i = 0; i < n_bits; i = i + 1 {
    fp.push(0)
  }
  // Compute ECFP excluding one atom
  let env_hashes : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    if i == exclude {
      env_hashes.push(0)
    } else {
      let mut h = 0
      let s = symbols[i]
      for c = 0; c < s.length(); c = c + 1 {
        h = h * 31 + s[c].to_int()
      }
      env_hashes.push(h)
    }
  }
  for r = 0; r < radius; r = r + 1 {
    let new_hashes : Array[Int32] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      if i == exclude {
        new_hashes.push(0)
      } else {
        let mut h = env_hashes[i] * 31
        for j = 0; j < n; j = j + 1 {
          if j != exclude && adj[i][j] > 0 {
            h = h + env_hashes[j] * adj[i][j]
          }
        }
        new_hashes.push(h & 0x7FFFFFFF)
      }
    }
    for i = 0; i < n; i = i + 1 {
      env_hashes[i] = new_hashes[i]
    }
  }
  for i = 0; i < n; i = i + 1 {
    if i != exclude {
      let bit = env_hashes[i] % n_bits
      fp[bit] = 1
    }
  }
  fp
}
