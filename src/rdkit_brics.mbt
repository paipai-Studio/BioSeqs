/// BRICS fragmentation.
///
/// RDKit equivalent: `rdkit.Chem.BRICS.BRICSDecompose`
///
/// BRICS (Breaking Retrosynthetically Interesting Chemical Substructures)
/// identifies bonds that can be broken for retrosynthetic analysis.
/// It uses a set of 16 SMARTS-defined bond types that are synthetically
/// meaningful to cleave.

///|
/// BRICS bond types: (type_name, atom1_smarts, atom2_smarts)
/// Simplified set of the most common BRICS bond types.
let rdkit_brics_bond_types : Array[(String, String, String)] = [
  ("L1", "[C;D3;R0]", "[N;D2;R0]"), // amine C-N
  ("L2", "[C;D3;R0]", "[O;D1;R0]"), // ester/ether C-O
  ("L3", "[C;D3;R0]", "[S;D1;R0]"), // thioether C-S
  ("L4", "[C;D3;R0]", "[C;D3;R0]"), // C-C single bond
  ("L5", "[C;D2;R0]", "[N;D2;R0]"), // amide C-N
  ("L6", "[C;D2;R0]", "[O;D1;R0]"), // carbonyl C-O
  ("L7", "[N;D2;R0]", "[C;D2;R0]"), // N-C (reverse amide)
  ("L8", "[C;D3;R0]", "[N;D3;R0]"), // C-N tertiary
  ("L9", "[C;D3;R0]", "[C;D2;R0]"), // C-C (sp3-sp2)
  ("L10", "[C;D2;R0]", "[C;D2;R0]"), // C-C (sp2-sp2)
  ("L11", "[C;D3;R0]", "[F,Cl,Br,I]"), // C-halogen
  ("L12", "[O;D2;R0]", "[C;D3;R0]"), // O-C ether
  ("L13", "[N;D2;R0]", "[N;D2;R0]"), // N-N
  ("L14", "[C;D3;R0]", "[P;D2;R0]"), // C-P
  ("L15", "[C;D3;R0]", "[B;D3;R0]"), // C-B
  ("L16", "[S;D2;R0]", "[C;D3;R0]"), // S-C
]

///|
/// A BRICS fragment result
pub struct RdkitBricsFragment {
  smiles : String
  bond_type : String
  atom_indices : Array[Int32]
}

///|
pub fn RdkitBricsFragment::smiles(self : RdkitBricsFragment) -> String {
  self.smiles
}

///|
pub fn RdkitBricsFragment::bond_type(self : RdkitBricsFragment) -> String {
  self.bond_type
}

///|
pub fn RdkitBricsFragment::atom_indices(
  self : RdkitBricsFragment,
) -> Array[Int32] {
  self.atom_indices
}

///|
/// Check if a bond matches a BRICS bond type.
/// Uses simple heuristic: check atom symbols and degrees.
fn brics_bond_matches(
  mol : Molecule,
  bond_idx : Int32,
  type_name : String,
) -> Bool {
  let bond = mol.bonds[bond_idx]
  let a1 = mol.atoms[bond.from]
  let a2 = mol.atoms[bond.to]
  let d1 = mol_degree(mol, bond.from)
  let d2 = mol_degree(mol, bond.to)
  // Don't break ring bonds
  if mol_is_in_ring(mol, bond.from) && mol_is_in_ring(mol, bond.to) {
    return false
  }
  match type_name {
    "L1" =>
      (a1.symbol == "C" && a2.symbol == "N") ||
      (a1.symbol == "N" && a2.symbol == "C")
    "L2" =>
      (a1.symbol == "C" && a2.symbol == "O") ||
      (a1.symbol == "O" && a2.symbol == "C")
    "L3" =>
      (a1.symbol == "C" && a2.symbol == "S") ||
      (a1.symbol == "S" && a2.symbol == "C")
    "L4" => a1.symbol == "C" && a2.symbol == "C" && bond.order == 1
    "L5" | "L7" =>
      (
        (a1.symbol == "C" && a2.symbol == "N") ||
        (a1.symbol == "N" && a2.symbol == "C")
      ) &&
      bond.order == 1
    "L6" => (a1.symbol == "C" && a2.symbol == "O") && bond.order == 2
    "L8" =>
      (a1.symbol == "C" && a2.symbol == "N" && d2 >= 3) ||
      (a2.symbol == "C" && a1.symbol == "N" && d1 >= 3)
    "L9" | "L10" => a1.symbol == "C" && a2.symbol == "C"
    "L11" =>
      (
        a1.symbol == "C" &&
        (
          a2.symbol == "F" ||
          a2.symbol == "Cl" ||
          a2.symbol == "Br" ||
          a2.symbol == "I"
        )
      ) ||
      (
        a2.symbol == "C" &&
        (
          a1.symbol == "F" ||
          a1.symbol == "Cl" ||
          a1.symbol == "Br" ||
          a1.symbol == "I"
        )
      )
    "L12" => a1.symbol == "O" && a2.symbol == "C" && d1 == 2
    "L13" => a1.symbol == "N" && a2.symbol == "N"
    "L14" =>
      (a1.symbol == "C" && a2.symbol == "P") ||
      (a1.symbol == "P" && a2.symbol == "C")
    "L15" =>
      (a1.symbol == "C" && a2.symbol == "B") ||
      (a1.symbol == "B" && a2.symbol == "C")
    "L16" => a1.symbol == "S" && a2.symbol == "C" && d1 == 2
    _ => false
  }
}

///|
/// Decompose a molecule into BRICS fragments.
/// Returns all cleavable bonds and their types.
pub fn rdkit_brics_decompose(mol : Molecule) -> Array[RdkitBricsFragment] {
  let fragments : Array[RdkitBricsFragment] = Array(capacity=0)
  for bi = 0; bi < mol.bonds.length(); bi = bi + 1 {
    for bt in rdkit_brics_bond_types {
      if brics_bond_matches(mol, bi, bt.0) {
        let bond = mol.bonds[bi]
        fragments.push({
          smiles: mol_to_smiles(mol),
          bond_type: bt.0,
          atom_indices: [bond.from, bond.to],
        })
        break
      }
    }
  }
  fragments
}

///|
/// Get unique BRICS bond types found in a molecule.
pub fn rdkit_brics_bond_types_found(mol : Molecule) -> Array[String] {
  let decomposed = rdkit_brics_decompose(mol)
  let types : Array[String] = Array(capacity=0)
  for f in decomposed {
    let mut found = false
    for t in types {
      if t == f.bond_type {
        found = true
        break
      }
    }
    if !found {
      types.push(f.bond_type)
    }
  }
  types
}

///|
/// Count number of BRICS-cleavable bonds.
pub fn rdkit_brics_count(mol : Molecule) -> Int32 {
  rdkit_brics_decompose(mol).length()
}
