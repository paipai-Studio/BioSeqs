/// rdkit_rings_depict.mbt — SSSR ring perception, aromatic perception, 2D depiction (SVG)
/// Items 18-20

///|
/// SSSR (Smallest Set of Smallest Rings) result
pub struct RdkitSSSR {
  rings : Array[Array[Int32]] // each ring is a list of atom indices
  count : Int32
}

///|
/// Compute SSSR using a simplified approach: find all cycles, select smallest.
pub fn rdkit_sssr(smiles : String) -> RdkitSSSR {
  let mol = parse_smiles(smiles)
  let n = mol.num_atoms()
  let adj = rdkit_build_adj_matrix(mol)
  // Build neighbor lists
  let neighbors : Array[Array[Int32]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let nb : Array[Int32] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        nb.push(j)
      }
    }
    neighbors.push(nb)
  }
  // Find all cycles via DFS
  let all_rings : Array[Array[Int32]] = Array(capacity=0)
  for start = 0; start < n; start = start + 1 {
    let path : Array[Int32] = Array(capacity=0)
    path.push(start)
    rdkit_sssr_dfs(start, start, path, neighbors, all_rings, n)
  }
  // Deduplicate rings (same ring can be found from different starting points)
  let unique_rings = rdkit_sssr_dedup(all_rings)
  // Select SSSR: smallest rings first
  let sorted_rings = rdkit_sssr_sort(unique_rings)
  // Number of rings in SSSR = bonds - atoms + components (Euler formula)
  let mut n_bonds = 0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        n_bonds = n_bonds + 1
      }
    }
  }
  // Simplified: assume 1 component
  let sssr_count = if n_bonds - n + 1 > 0 { n_bonds - n + 1 } else { 0 }
  // Take the smallest sssr_count rings
  let selected : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < sorted_rings.length(); i = i + 1 {
    if selected.length() >= sssr_count {
      break
    }
    selected.push(sorted_rings[i])
  }
  RdkitSSSR::{ rings: selected, count: selected.length(), }
}

///|
fn rdkit_sssr_dfs(
  start : Int32,
  current : Int32,
  path : Array[Int32],
  neighbors : Array[Array[Int32]],
  rings : Array[Array[Int32]],
  max_depth : Int32,
) -> Unit {
  if path.length() > 12 { // limit ring size
    return
  }
  let nb = neighbors[current]
  for i = 0; i < nb.length(); i = i + 1 {
    let next = nb[i]
    if next == start && path.length() >= 3 {
      // Found a ring
      let ring : Array[Int32] = Array(capacity=path.length())
      for j = 0; j < path.length(); j = j + 1 {
        ring.push(path[j])
      }
      rings.push(ring)
    } else if !rdkit_sssr_in_path(path, next) && next > start {
      path.push(next)
      rdkit_sssr_dfs(start, next, path, neighbors, rings, max_depth)
      let _ = path.remove(path.length() - 1)
    }
  }
}

///|
fn rdkit_sssr_in_path(path : Array[Int32], atom : Int32) -> Bool {
  for i = 0; i < path.length(); i = i + 1 {
    if path[i] == atom {
      return true
    }
  }
  false
}

///|
fn rdkit_sssr_dedup(rings : Array[Array[Int32]]) -> Array[Array[Int32]] {
  let unique : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < rings.length(); i = i + 1 {
    let mut is_dup = false
    for j = 0; j < unique.length(); j = j + 1 {
      if rdkit_sssr_same_ring(rings[i], unique[j]) {
        is_dup = true
      }
    }
    if !is_dup {
      unique.push(rings[i])
    }
  }
  unique
}

///|
fn rdkit_sssr_same_ring(r1 : Array[Int32], r2 : Array[Int32]) -> Bool {
  if r1.length() != r2.length() {
    return false
  }
  let s1 : Set[Int32] = Set([])
  for i = 0; i < r1.length(); i = i + 1 {
    s1.add(r1[i])
  }
  for i = 0; i < r2.length(); i = i + 1 {
    if !s1.contains(r2[i]) {
      return false
    }
  }
  true
}

///|
fn rdkit_sssr_sort(rings : Array[Array[Int32]]) -> Array[Array[Int32]] {
  let sorted = Array(capacity=rings.length())
  for i = 0; i < rings.length(); i = i + 1 {
    sorted.push(rings[i])
  }
  sorted.sort_by(fn(a, b) { a.length() - b.length() })
  sorted
}

///|
/// Aromatic perception result
pub struct RdkitAromaticInfo {
  aromatic_atoms : Array[Int32]
  aromatic_rings : Array[Array[Int32]]
  n_aromatic_atoms : Int32
  n_aromatic_rings : Int32
  n_hetero_aromatic_rings : Int32
}

///|
/// Perceive aromaticity in detail: count aromatic atoms and rings.
pub fn rdkit_perceive_aromaticity(smiles : String) -> RdkitAromaticInfo {
  let mol = parse_smiles(smiles)
  let _n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let sssr = rdkit_sssr(smiles)
  let aromatic_atoms : Array[Int32] = Array(capacity=0)
  let aromatic_rings : Array[Array[Int32]] = Array(capacity=0)
  let mut n_hetero = 0
  // Check each ring for aromaticity (simplified: 4n+2 rule)
  for i = 0; i < sssr.rings.length(); i = i + 1 {
    let ring = sssr.rings[i]
    let ring_size = ring.length()
    // Hückel: 4n+2 = ring_size → n = (ring_size - 2) / 4
    let n_huckel = ring_size - 2
    if n_huckel % 4 == 0 && n_huckel >= 0 {
      // Potentially aromatic
      let mut has_hetero = false
      for j = 0; j < ring_size; j = j + 1 {
        let atom_idx = ring[j]
        let sym = symbols[atom_idx]
        if sym == "N" || sym == "O" || sym == "S" {
          has_hetero = true
        }
        // Add to aromatic atoms if not already there
        let mut found = false
        for k = 0; k < aromatic_atoms.length(); k = k + 1 {
          if aromatic_atoms[k] == atom_idx {
            found = true
          }
        }
        if !found {
          aromatic_atoms.push(atom_idx)
        }
      }
      aromatic_rings.push(ring)
      if has_hetero {
        n_hetero = n_hetero + 1
      }
    }
  }
  RdkitAromaticInfo::{
    aromatic_atoms,
    aromatic_rings,
    n_aromatic_atoms: aromatic_atoms.length(),
    n_aromatic_rings: aromatic_rings.length(),
    n_hetero_aromatic_rings: n_hetero,
  }
}

///|
/// Count aromatic atoms in a molecule.
pub fn rdkit_count_aromatic_atoms(smiles : String) -> Int32 {
  rdkit_perceive_aromaticity(smiles).n_aromatic_atoms
}

///|
/// Count aromatic rings in a molecule.
pub fn rdkit_count_aromatic_rings(smiles : String) -> Int32 {
  rdkit_perceive_aromaticity(smiles).n_aromatic_rings
}

///|
/// Check if a ring is heteroaromatic.
pub fn rdkit_is_heteroaromatic(smiles : String, ring_idx : Int32) -> Bool {
  let info = rdkit_perceive_aromaticity(smiles)
  if ring_idx < 0 || ring_idx >= info.aromatic_rings.length() {
    return false
  }
  let mol = parse_smiles(smiles)
  let symbols = mol.atom_symbols()
  let ring = info.aromatic_rings[ring_idx]
  for i = 0; i < ring.length(); i = i + 1 {
    let s = symbols[ring[i]]
    if s == "N" || s == "O" || s == "S" {
      return true
    }
  }
  false
}

///|
/// Generate SVG depiction of a molecule (2D text-based rendering).
pub fn rdkit_depict_svg(
  smiles : String,
  width~ : Int32 = 300,
  height~ : Int32 = 300,
) -> String {
  let mol = parse_smiles(smiles)
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let adj = rdkit_build_adj_matrix(mol)
  // Generate 2D coordinates
  let coords = rdkit_depict_layout(n, adj, width, height)
  let sb = StringBuilder()
  let _ = sb.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"")
  let _ = sb.write(width.to_string())
  let _ = sb.write("\" height=\"")
  let _ = sb.write(height.to_string())
  let _ = sb.write("\">\n")
  // Background
  let _ = sb.write("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n")
  // Draw bonds
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        let (x1, y1) = coords[i]
        let (x2, y2) = coords[j]
        let _ = sb.write("<line x1=\"")
        let _ = sb.write(x1.to_string())
        let _ = sb.write("\" y1=\"")
        let _ = sb.write(y1.to_string())
        let _ = sb.write("\" x2=\"")
        let _ = sb.write(x2.to_string())
        let _ = sb.write("\" y2=\"")
        let _ = sb.write(y2.to_string())
        let _ = sb.write("\" stroke=\"black\" stroke-width=\"1.5\"/>\n")
      }
    }
  }
  // Draw atoms
  for i = 0; i < n; i = i + 1 {
    let (x, y) = coords[i]
    let _ = sb.write("<text x=\"")
    let _ = sb.write(x.to_string())
    let _ = sb.write("\" y=\"")
    let _ = sb.write(y.to_string())
    let _ = sb.write("\" text-anchor=\"middle\" font-size=\"14\" fill=\"")
    let color = rdkit_depict_atom_color(symbols[i])
    let _ = sb.write(color)
    let _ = sb.write("\">")
    let _ = sb.write(symbols[i])
    let _ = sb.write("</text>\n")
  }
  let _ = sb.write("</svg>")
  sb.to_string()
}

///|
fn rdkit_depict_layout(
  n : Int32,
  _adj : Array[Array[Int32]],
  width : Int32,
  height : Int32,
) -> Array[(Double, Double)] {
  let coords : Array[(Double, Double)] = Array(capacity=n)
  let cx = width.to_double() / 2.0
  let cy = height.to_double() / 2.0
  if n == 0 {
    return coords
  }
  if n == 1 {
    coords.push((cx, cy))
    return coords
  }
  // Place atoms in a circle
  let min_dim = if width < height { width } else { height }
  let radius = min_dim.to_double() * 0.3
  for i = 0; i < n; i = i + 1 {
    let angle = i.to_double() * 2.0 * 3.141592653589793 / n.to_double()
    let x = cx + radius * @math.cos(angle)
    let y = cy + radius * @math.sin(angle)
    coords.push((x, y))
  }
  coords
}

///|
fn rdkit_depict_atom_color(sym : String) -> String {
  match sym {
    "C" => "black"
    "N" => "blue"
    "O" => "red"
    "S" => "yellow"
    "F" => "green"
    "Cl" => "green"
    "Br" => "brown"
    "I" => "purple"
    "P" => "orange"
    "H" => "gray"
    _ => "black"
  }
}

///|
/// Generate a simplified SVG with just atom labels (no bonds).
pub fn rdkit_depict_atoms_only(smiles : String) -> String {
  let mol = parse_smiles(smiles)
  let symbols = mol.atom_symbols()
  let n = mol.num_atoms()
  let sb = StringBuilder()
  let _ = sb.write(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"100\">\n",
  )
  for i = 0; i < n; i = i + 1 {
    let x = 20 + i * 30
    let _ = sb.write("<text x=\"")
    let _ = sb.write(x.to_string())
    let _ = sb.write("\" y=\"50\" font-size=\"14\">")
    let _ = sb.write(symbols[i])
    let _ = sb.write("</text>\n")
  }
  let _ = sb.write("</svg>")
  sb.to_string()
}
