/// rdkit_mmp.mbt — matched molecular pair cutting (RDKit gap fill, A4)
///
/// Single and double cuts on non-ring single bonds. Each cut yields a
/// context (the larger remainder, multi-component contexts joined with
/// ".") and the detached fragment(s) as canonical SMILES.

///|
/// A matched molecular pair cut: context SMILES + detached fragments.
pub struct RdkitMmpCut {
  context : String
  fragments : Array[String]
  cut_bonds : Array[(Int32, Int32)]
}

///|
fn mmp_candidate_bonds(mol : Molecule) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    let b = mol.bonds[i]
    if b.order == 1 && !mol_is_in_ring(mol, b.from) && !mol_is_in_ring(mol, b.to) {
      out.push(i)
    }
  }
  out
}

///|
/// Extract the component containing `seed` while skipping the given bond
/// indices. Returns member atom indices.
fn mmp_component(
  mol : Molecule,
  seed : Int32,
  cut_bonds : Set[Int32],
) -> Array[Int32] {
  let n = mol.num_atoms()
  let visited : Array[Bool] = Array::make(n, false)
  let out : Array[Int32] = Array(capacity=0)
  let queue : Array[Int32] = Array(capacity=0)
  queue.push(seed)
  visited[seed] = true
  let mut head = 0
  while head < queue.length() {
    let cur = queue[head]
    head = head + 1
    out.push(cur)
    for b in mol.bonds {
      let other = if b.from == cur {
        b.to
      } else if b.to == cur {
        b.from
      } else {
        -1
      }
      if other >= 0 && !visited[other] {
        let mut is_cut = false
        for i = 0; i < mol.bonds.length(); i = i + 1 {
          if cut_bonds.contains(i) {
            let bb = mol.bonds[i]
            if (bb.from == cur && bb.to == other) || (bb.from == other && bb.to == cur) {
              is_cut = true
            }
          }
        }
        if is_cut {
          continue
        }
        visited[other] = true
        queue.push(other)
      }
    }
  }
  out
}

///|
fn mmp_sub_mol(mol : Molecule, atoms : Array[Int32]) -> Molecule {
  let in_set : Array[Bool] = Array::make(mol.num_atoms(), false)
  for a in atoms {
    in_set[a] = true
  }
  let index_map : Map[Int32, Int32] = Map([], capacity=0)
  let new_atoms : Array[MolAtom] = Array(capacity=0)
  for a in atoms {
    index_map[a] = new_atoms.length()
    new_atoms.push(MolAtom::{ symbol: mol.atoms[a].symbol, charge: mol.atoms[a].charge, aromatic: mol.atoms[a].aromatic, })
  }
  let new_bonds : Array[MolBond] = Array(capacity=0)
  for b in mol.bonds {
    if in_set[b.from] && in_set[b.to] {
      new_bonds.push(MolBond::{ from: index_map[b.from], to: index_map[b.to], order: b.order, })
    }
  }
  Molecule::{ atoms: new_atoms, bonds: new_bonds, }
}

///|
fn mmp_component_smiles(mol : Molecule, atoms : Array[Int32]) -> String {
  mol_to_smiles(mol_perceive_aromaticity(mmp_sub_mol(mol, atoms)))
}

///|
/// Split the molecule into connected components after cutting the given
/// bond index set. Components are returned sorted by size (small first).
fn mmp_components_sorted(
  mol : Molecule,
  cuts : Set[Int32],
) -> Array[Array[Int32]] {
  let n = mol.num_atoms()
  let visited : Array[Bool] = Array::make(n, false)
  let comps : Array[Array[Int32]] = Array(capacity=0)
  for seed = 0; seed < n; seed = seed + 1 {
    if !visited[seed] {
      let comp = mmp_component(mol, seed, cuts)
      for a in comp {
        visited[a] = true
      }
      comps.push(comp)
    }
  }
  // insertion sort by size ascending
  for i = 1; i < comps.length(); i = i + 1 {
    let key = comps[i]
    let mut j = i - 1
    while j >= 0 && comps[j].length() > key.length() {
      comps[j + 1] = comps[j]
      j = j - 1
    }
    comps[j + 1] = key
  }
  comps
}

///|
/// Enumerate matched molecular pair cuts with 1 or 2 bond cuts.
/// Candidate non-ring single bonds are capped at 10 to bound the
/// combinatorics. Fragments are the smaller detached components;
/// the context is the remainder (joined with "." when several pieces).
pub fn rdkit_mmp_cuts(
  mol : Molecule,
  max_cuts~ : Int32 = 1,
) -> Array[RdkitMmpCut] {
  let out : Array[RdkitMmpCut] = Array(capacity=0)
  let cands = mmp_candidate_bonds(mol)
  let limit = if cands.length() > 10 { 10 } else { cands.length() }
  let total_atoms = mol.num_atoms()
  for i = 0; i < limit; i = i + 1 {
    // Single cut
    let cut1 : Set[Int32] = Set([cands[i]])
    let comps1 = mmp_components_sorted(mol, cut1)
    if comps1.length() == 2 {
      let frag = mmp_component_smiles(mol, comps1[0])
      let ctx = mmp_component_smiles(mol, comps1[1])
      out.push(RdkitMmpCut::{
          context: ctx,
          fragments: [frag],
          cut_bonds: [(mol.bonds[cands[i]].from, mol.bonds[cands[i]].to)],
        })
    }
    if max_cuts >= 2 {
      for j = i + 1; j < limit; j = j + 1 {
        let cut2 : Set[Int32] = Set([cands[i], cands[j]])
        let comps2 = mmp_components_sorted(mol, cut2)
        if comps2.length() >= 2 {
          // fragments = components that are <= half of the molecule
          let frags : Array[String] = Array(capacity=0)
          let rest : Array[String] = Array(capacity=0)
          for k = 0; k < comps2.length(); k = k + 1 {
            let comp = comps2[k]
            if comp.length() * 2 <= total_atoms && k < comps2.length() - 1 {
              frags.push(mmp_component_smiles(mol, comp))
            } else {
              rest.push(mmp_component_smiles(mol, comp))
            }
          }
          if frags.length() > 0 && rest.length() > 0 {
            let sb = StringBuilder()
            for k = 0; k < rest.length(); k = k + 1 {
              if k > 0 {
                let _ = sb.write(".")
              }
              let _ = sb.write(rest[k])
            }
            out.push(RdkitMmpCut::{
                context: sb.to_string(),
                fragments: frags,
                cut_bonds: [
                  (mol.bonds[cands[i]].from, mol.bonds[cands[i]].to),
                  (mol.bonds[cands[j]].from, mol.bonds[cands[j]].to),
                ],
              })
          }
        }
      }
    }
  }
  out
}
