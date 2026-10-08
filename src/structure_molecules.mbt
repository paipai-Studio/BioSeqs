///|
/// Connected molecule detection — find separate molecules in a structure via
/// bond graph connected components.
/// Ported from biotite.structure.molecules.

///|
/// A molecule identified in a structure: atom indices forming a connected component.
pub struct ConnectedMolecule {
  atom_indices : Array[Int]   // indices into the original atom array
  center : Vector3            // geometric center
  n_atoms : Int
}

///|
/// Detect connected molecules using a bond-distance criterion.
/// Two atoms are bonded if their distance ≤ (r_cov_a + r_cov_b + tolerance).
/// Uses union-find for O(n·α(n)) connected components.
///
/// Parameters:
///   `atoms` — flat atom list
///   `tolerance` — extra tolerance on covalent radii sum (Å), default 0.4
pub fn detect_molecules(
  atoms : Array[Atom],
  tolerance~ : Double = 0.4,
) -> Array[ConnectedMolecule] {
  let n = atoms.length()
  // Union-find parent array
  let parent : Array[Int] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    parent.push(i)
  }
  // Find with path compression
  fn find(parent : Array[Int], x : Int) -> Int {
    let mut root = x
    while parent[root] != root {
      root = parent[root]
    }
    // Path compression
    let mut curr = x
    while parent[curr] != root {
      let next = parent[curr]
      parent[curr] = root
      curr = next
    }
    root
  }
  // Union by checking bond connectivity
  for i = 0; i < n; i = i + 1 {
    let ri = covalent_radius(atoms[i].element)
    for j = i + 1; j < n; j = j + 1 {
      let rj = covalent_radius(atoms[j].element)
      let bond_threshold = ri + rj + tolerance
      let dist = atoms[i].coord.distance(atoms[j].coord)
      if dist <= bond_threshold && dist > 0.3 {
        let ri_root = find(parent, i)
        let rj_root = find(parent, j)
        if ri_root != rj_root {
          parent[ri_root] = rj_root
        }
      }
    }
  }
  // Group atoms by root
  let groups : Map[Int, Array[Int]] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    let root = find(parent, i)
    match groups.get(root) {
      Some(arr) => arr.push(i)
      None => {
        let arr : Array[Int] = Array(capacity=1)
        arr.push(i)
        groups[root] = arr
      }
    }
  }
  // Build Molecule structs
  let molecules : Array[ConnectedMolecule] = Array(capacity=groups.length())
  for kv in groups {
    let indices = kv.1
    let mut cx = 0.0
    let mut cy = 0.0
    let mut cz = 0.0
    for idx in indices {
      cx = cx + atoms[idx].coord.x
      cy = cy + atoms[idx].coord.y
      cz = cz + atoms[idx].coord.z
    }
    let m = indices.length().to_double()
    molecules.push(ConnectedMolecule::{
      atom_indices: indices,
      center: Vector3::new(cx / m, cy / m, cz / m),
      n_atoms: indices.length(),
    })
  }
  // Sort by size descending
  for i = 1; i < molecules.length(); i = i + 1 {
    let key = molecules[i]
    let mut j = i - 1
    while j >= 0 && molecules[j].n_atoms < key.n_atoms {
      molecules[j + 1] = molecules[j]
      j = j - 1
    }
    molecules[j + 1] = key
  }
  molecules
}
