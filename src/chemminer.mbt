// Bioconductor ChemmineR - Cheminformatics toolkit for small-molecule analysis.
//
// Inspired by the Bioconductor ChemmineR package (Backman et al. 2011,
// Bioinformatics 27:2400-1), which provides infrastructure for analyzing
// small molecules. This pure-MoonBit module implements:
// - SDF (Structure-Data File) parsing and serialization
// - Molecular descriptors (MW, formula, atom/bond counts, ring detection)
// - Atom-pair fingerprint generation
// - Tanimoto and Dice similarity coefficients
// - Simplified substructure search
//
// Reference:
//   https://bioconductor.org/packages/ChemmineR/

// ============================================================================
// Element data - atomic weights and symbols
// ============================================================================

///|
/// Atomic weights for common organic-chemistry elements.
/// Values are monoisotopic-average hybrid weights sufficient for
/// descriptor-level calculations.
fn chem_atomic_weight(symbol : String) -> Double {
  let up = symbol.to_upper()
  match up {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "P" => 30.974
    "S" => 32.06
    "CL" => 35.45
    "BR" => 79.904
    "I" => 126.904
    "B" => 10.81
    "SI" => 28.085
    "SE" => 78.971
    "NA" => 22.990
    "K" => 39.098
    "MG" => 24.305
    "CA" => 40.078
    "FE" => 55.845
    "ZN" => 65.38
    "MN" => 54.938
    "CU" => 63.546
    "CO" => 58.933
    "NI" => 58.693
    "AL" => 26.982
    "AS" => 74.922
    "LI" => 6.94
    _ => 0.0
  }
}

///|
/// Whether a symbol is a recognized element.
fn chem_is_element(symbol : String) -> Bool {
  chem_atomic_weight(symbol) > 0.0
}

// ============================================================================
// SDF atom and bond
// ============================================================================

///|
/// A single atom in an SDF molecule: element symbol and 3D coordinates.
pub struct SdfAtom {
  symbol : String
  x : Double
  y : Double
  z : Double
  charge : Int // formal charge (0 = neutral)
}

///|
/// Construct an SdfAtom.
pub fn SdfAtom::new(
  symbol : String,
  x : Double,
  y : Double,
  z : Double,
  charge : Int
) -> SdfAtom {
  SdfAtom::{ symbol, x, y, z, charge }
}

///|
pub fn SdfAtom::symbol(self : SdfAtom) -> String {
  self.symbol
}

///|
pub fn SdfAtom::x(self : SdfAtom) -> Double {
  self.x
}

///|
pub fn SdfAtom::y(self : SdfAtom) -> Double {
  self.y
}

///|
pub fn SdfAtom::z(self : SdfAtom) -> Double {
  self.z
}

///|
pub fn SdfAtom::charge(self : SdfAtom) -> Int {
  self.charge
}

///|
/// Euclidean distance between two atoms.
pub fn SdfAtom::distance(self : SdfAtom, other : SdfAtom) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// A single bond in an SDF molecule: atom indices (1-based in SDF, stored
/// 0-based internally) and bond type.
pub struct SdfBond {
  atom1 : Int // 0-based index into the molecule's atom array
  atom2 : Int
  bond_type : Int // 1=single, 2=double, 3=triple, 4=aromatic
}

///|
/// Construct an SdfBond.
pub fn SdfBond::new(atom1 : Int, atom2 : Int, bond_type : Int) -> SdfBond {
  SdfBond::{ atom1, atom2, bond_type }
}

///|
pub fn SdfBond::atom1(self : SdfBond) -> Int {
  self.atom1
}

///|
pub fn SdfBond::atom2(self : SdfBond) -> Int {
  self.atom2
}

///|
pub fn SdfBond::bond_type(self : SdfBond) -> Int {
  self.bond_type
}

///|
/// Human-readable bond type label.
pub fn SdfBond::bond_label(self : SdfBond) -> String {
  match self.bond_type {
    1 => "single"
    2 => "double"
    3 => "triple"
    4 => "aromatic"
    _ => "unknown"
  }
}

// ============================================================================
// SDF molecule
// ============================================================================

///|
/// A parsed SDF molecule: name, atoms, bonds, and properties.
pub struct SdfMolecule {
  name : String
  atoms : Array[SdfAtom]
  bonds : Array[SdfBond]
  properties : Map[String, String]
}

///|
/// Construct an SdfMolecule from parts.
pub fn SdfMolecule::new(
  name : String,
  atoms : Array[SdfAtom],
  bonds : Array[SdfBond]
) -> SdfMolecule {
  SdfMolecule::{ name, atoms, bonds, properties: Map([], capacity=8) }
}

///|
pub fn SdfMolecule::name(self : SdfMolecule) -> String {
  self.name
}

///|
pub fn SdfMolecule::atoms(self : SdfMolecule) -> Array[SdfAtom] {
  self.atoms
}

///|
pub fn SdfMolecule::bonds(self : SdfMolecule) -> Array[SdfBond] {
  self.bonds
}

///|
pub fn SdfMolecule::properties(self : SdfMolecule) -> Map[String, String] {
  self.properties
}

///|
/// Number of atoms.
pub fn SdfMolecule::n_atoms(self : SdfMolecule) -> Int {
  self.atoms.length()
}

///|
/// Number of bonds.
pub fn SdfMolecule::n_bonds(self : SdfMolecule) -> Int {
  self.bonds.length()
}

///|
/// Get a property value by key. Returns empty string if missing.
pub fn SdfMolecule::get_property(
  self : SdfMolecule,
  key : String
) -> String {
  self.properties.get_or_default(key, "")
}

///|
/// Set a property value.
pub fn SdfMolecule::set_property(
  self : SdfMolecule,
  key : String,
  value : String
) -> Unit {
  self.properties[key] = value
}

///|
/// Whether the molecule has a property.
pub fn SdfMolecule::has_property(
  self : SdfMolecule,
  key : String
) -> Bool {
  self.properties.contains(key)
}

// ============================================================================
// Molecular descriptors
// ============================================================================

///|
/// Molecular weight: sum of atomic weights of all atoms.
pub fn chem_sdf_molecular_weight(mol : SdfMolecule) -> Double {
  let mut mw = 0.0
  for a in mol.atoms {
    mw = mw + chem_atomic_weight(a.symbol)
  }
  mw
}

///|
/// Molecular formula as a Hill-order string (C first, H second, then
/// remaining elements alphabetically).
pub fn chem_molecular_formula(mol : SdfMolecule) -> String {
  let counts : Map[String, Int] = Map([], capacity=16)
  for a in mol.atoms {
    let up = a.symbol.to_upper()
    counts[up] = counts.get_or_default(up, 0) + 1
  }
  // Hill order: C, H, then alphabetical.
  let sb = StringBuilder::new()
  // Carbon
  if counts.contains("C") {
    let c = counts["C"]
    sb.write_string("C")
    if c > 1 {
      sb.write_string(c.to_string())
    }
    counts.remove("C")
  }
  // Hydrogen
  if counts.contains("H") {
    let h = counts["H"]
    sb.write_string("H")
    if h > 1 {
      sb.write_string(h.to_string())
    }
    counts.remove("H")
  }
  // Remaining elements in alphabetical order.
  let remaining : Array[String] = Array::new()
  for k in counts.keys() {
    remaining.push(k)
  }
  remaining.sort_by(fn(a : String, b : String) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  for elem in remaining {
    let c = counts[elem]
    sb.write_string(elem)
    if c > 1 {
      sb.write_string(c.to_string())
    }
  }
  sb.to_string()
}

///|
/// Count of heavy atoms (non-hydrogen).
pub fn chem_heavy_atom_count(mol : SdfMolecule) -> Int {
  let mut n = 0
  for a in mol.atoms {
    if a.symbol.to_upper() != "H" {
      n = n + 1
    }
  }
  n
}

///|
/// Count atoms of a specific element (case-insensitive).
pub fn chem_count_element(mol : SdfMolecule, element : String) -> Int {
  let up = element.to_upper()
  let mut n = 0
  for a in mol.atoms {
    if a.symbol.to_upper() == up {
      n = n + 1
    }
  }
  n
}

///|
/// Number of rings detected using a simple approach: for a connected
/// graph, the number of independent cycles (cyclomatic number) equals
/// E - V + C where E = edges, V = vertices, C = connected components.
pub fn chem_ring_count(mol : SdfMolecule) -> Int {
  let n_atoms = mol.atoms.length()
  if n_atoms == 0 {
    return 0
  }
  let n_bonds = mol.bonds.length()
  // Count connected components via union-find.
  let parent : Array[Int] = Array::make(n_atoms, 0)
  for i in 0..<n_atoms {
    parent[i] = i
  }
  let find = fn(parent : Array[Int], x : Int) -> Int {
    let mut root = x
    while parent[root] != root {
      root = parent[root]
    }
    // Path compression.
    let mut cur = x
    while parent[cur] != root {
      let next = parent[cur]
      parent[cur] = root
      cur = next
    }
    root
  }
  for b in mol.bonds {
    let ra = find(parent, b.atom1)
    let rb = find(parent, b.atom2)
    if ra != rb {
      parent[ra] = rb
    }
  }
  let mut n_components = 0
  for i in 0..<n_atoms {
    if find(parent, i) == i {
      n_components = n_components + 1
    }
  }
  // Cyclomatic number.
  let rings = n_bonds - n_atoms + n_components
  if rings < 0 { 0 } else { rings }
}

// ============================================================================
// Atom-pair fingerprints
// ============================================================================

///|
/// Encode an atom's "type" as an integer for fingerprinting.
/// Type = element index * 10 + degree.
/// Element indices: H=1, C=2, N=3, O=4, F=5, P=6, S=7, Cl=8, Br=9, I=10,
/// others=0.
fn chem_atom_type(atom : SdfAtom, degree : Int) -> Int {
  let idx = match atom.symbol.to_upper() {
    "H" => 1
    "C" => 2
    "N" => 3
    "O" => 4
    "F" => 5
    "P" => 6
    "S" => 7
    "CL" => 8
    "BR" => 9
    "I" => 10
    _ => 0
  }
  idx * 10 + degree
}

///|
/// Compute the degree (number of bonds) for each atom.
fn chem_atom_degrees(mol : SdfMolecule) -> Array[Int] {
  let degrees = Array::make(mol.n_atoms(), 0)
  for b in mol.bonds {
    degrees[b.atom1] = degrees[b.atom1] + 1
    degrees[b.atom2] = degrees[b.atom2] + 1
  }
  degrees
}

///|
/// BFS shortest-path distance between two atom indices in the molecular
/// graph. Returns -1 if unreachable.
fn chem_bfs_distance(
  mol : SdfMolecule,
  start : Int,
  end : Int
) -> Int {
  if start == end {
    return 0
  }
  let n = mol.n_atoms()
  let visited = Array::make(n, false)
  let dist = Array::make(n, -1)
  // Build adjacency list.
  let adj : Array[Array[Int]] = Array::new()
  for _ in 0..<n {
    adj.push(Array::new())
  }
  for b in mol.bonds {
    adj[b.atom1].push(b.atom2)
    adj[b.atom2].push(b.atom1)
  }
  // BFS.
  let queue : Array[Int] = Array::new()
  queue.push(start)
  visited[start] = true
  dist[start] = 0
  let mut head = 0
  while head < queue.length() {
    let cur = queue[head]
    head = head + 1
    for nb in adj[cur] {
      if !visited[nb] {
        visited[nb] = true
        dist[nb] = dist[cur] + 1
        if nb == end {
          return dist[nb]
        }
        queue.push(nb)
      }
    }
  }
  -1
}

///|
/// Generate the atom-pair fingerprint of a molecule as a sorted set of
/// integer hash codes. Each atom pair (i, j) with i < j encodes:
///   hash = type_i * 10000 + dist * 100 + type_j
/// where type = element_index * 10 + degree and dist is the topological
/// (graph) distance. The set of unique hashes is the fingerprint.
pub fn chem_atom_pair_fingerprint(mol : SdfMolecule) -> Array[Int] {
  let n = mol.n_atoms()
  if n < 2 {
    return []
  }
  let degrees = chem_atom_degrees(mol)
  let types = Array::make(n, 0)
  for i in 0..<n {
    types[i] = chem_atom_type(mol.atoms[i], degrees[i])
  }
  let fp_set : Map[Int, Bool] = Map([], capacity=64)
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d = chem_bfs_distance(mol, i, j)
      if d < 0 {
        continue
      }
      let ti = types[i]
      let tj = types[j]
      // Canonical: smaller type first.
      let hash = if ti <= tj {
        ti * 10000 + d * 100 + tj
      } else {
        tj * 10000 + d * 100 + ti
      }
      fp_set[hash] = true
    }
  }
  let result : Array[Int] = Array::new()
  for k in fp_set.keys() {
    result.push(k)
  }
  result.sort_by(fn(a : Int, b : Int) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  result
}

// ============================================================================
// Similarity coefficients
// ============================================================================

///|
/// Intersection size of two sorted integer arrays.
fn chem_sorted_intersection(a : Array[Int], b : Array[Int]) -> Int {
  let mut i = 0
  let mut j = 0
  let mut count = 0
  while i < a.length() && j < b.length() {
    if a[i] == b[j] {
      count = count + 1
      i = i + 1
      j = j + 1
    } else if a[i] < b[j] {
      i = i + 1
    } else {
      j = j + 1
    }
  }
  count
}

///|
/// Tanimoto similarity between two fingerprints (sorted int arrays):
///   T = |A ∩ B| / |A ∪ B|
/// Returns 0.0 if both are empty.
pub fn chem_tanimoto(fp1 : Array[Int], fp2 : Array[Int]) -> Double {
  let inter = chem_sorted_intersection(fp1, fp2)
  let union_size = fp1.length() + fp2.length() - inter
  if union_size == 0 {
    return 0.0
  }
  inter.to_double() / union_size.to_double()
}

///|
/// Dice similarity between two fingerprints:
///   D = 2|A ∩ B| / (|A| + |B|)
/// Returns 0.0 if both are empty.
pub fn chem_dice(fp1 : Array[Int], fp2 : Array[Int]) -> Double {
  let inter = chem_sorted_intersection(fp1, fp2)
  let total = fp1.length() + fp2.length()
  if total == 0 {
    return 0.0
  }
  2.0 * inter.to_double() / total.to_double()
}

// ============================================================================
// Substructure search
// ============================================================================

///|
/// Recursive backtracking helper for `chem_substructure_match`. MoonBit does
/// not support recursive local functions, so this is a top-level function
/// taking all needed parameters.
fn chem_substructure_try_match(
  qn : Int,
  query : SdfMolecule,
  target : SdfMolecule,
  q_adj : Array[Array[(Int, Int)]],
  t_adj : Array[Array[(Int, Int)]],
  mapping : Array[Int],
  used : Array[Bool],
  i : Int
) -> Bool {
  if i == qn {
    return true
  }
  let tn = target.n_atoms()
  let qsym = query.atoms[i].symbol.to_upper()
  let mut t = 0
  while t < tn {
    if !used[t] && target.atoms[t].symbol.to_upper() == qsym {
      // Check bond consistency with already-mapped query atoms.
      let mut ok = true
      let mut k = 0
      while k < q_adj[i].length() {
        let (qnb, qbt) = q_adj[i][k]
        if qnb < i && mapping[qnb] >= 0 {
          // Hydrogen atoms are often implicit in chemistry; skip bond
          // checks involving H so that explicit-H queries still match
          // targets whose heavy-atom skeleton is compatible.
          let i_is_h = qsym == "H"
          let qnb_sym = query.atoms[qnb].symbol.to_upper()
          if !i_is_h && qnb_sym != "H" {
            let mapped = mapping[qnb]
            // Find the bond in target.
            let mut bond_ok = false
            let mut m = 0
            while m < t_adj[t].length() {
              let (tnb, tbt) = t_adj[t][m]
              if tnb == mapped && tbt >= qbt {
                bond_ok = true
                break
              }
              m = m + 1
            }
            if !bond_ok {
              ok = false
              break
            }
          }
        }
        k = k + 1
      }
      if ok {
        mapping[i] = t
        used[t] = true
        if chem_substructure_try_match(
          qn, query, target, q_adj, t_adj, mapping, used, i + 1
        ) {
          return true
        } else {
          mapping[i] = -1
          used[t] = false
        }
      }
    }
    t = t + 1
  }
  false
}

///|
/// Check whether `query` is a substructure of `target`. A substructure match
/// requires that every atom in the query maps to a distinct atom in the
/// target with the same element symbol, and every bond in the query maps
/// to a bond in the target with the same or higher bond type. Uses a
/// greedy backtracking search.
pub fn chem_substructure_match(
  query : SdfMolecule,
  target : SdfMolecule
) -> Bool {
  let qn = query.n_atoms()
  let tn = target.n_atoms()
  if qn == 0 {
    return true
  }
  if qn > tn {
    return false
  }
  // Build target adjacency: for each atom, list of (neighbor, bond_type).
  let t_adj : Array[Array[(Int, Int)]] = Array::new()
  for _ in 0..<tn {
    t_adj.push(Array::new())
  }
  for b in target.bonds {
    t_adj[b.atom1].push((b.atom2, b.bond_type))
    t_adj[b.atom2].push((b.atom1, b.bond_type))
  }
  // Build query adjacency.
  let q_adj : Array[Array[(Int, Int)]] = Array::new()
  for _ in 0..<qn {
    q_adj.push(Array::new())
  }
  for b in query.bonds {
    q_adj[b.atom1].push((b.atom2, b.bond_type))
    q_adj[b.atom2].push((b.atom1, b.bond_type))
  }
  // Backtracking assignment: query atom i -> target atom mapping[i].
  let mapping = Array::make(qn, -1)
  let used = Array::make(tn, false)
  chem_substructure_try_match(qn, query, target, q_adj, t_adj, mapping, used, 0)
}

// ============================================================================
// SDF parsing
// ============================================================================

///|
/// Trim leading/trailing whitespace from a string.
fn chem_trim(s : StringView) -> String {
  let n = s.length()
  let mut start = 0
  while start < n {
    let ch = s[start : start + 1]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let ch = s[end - 1 : end]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_owned()
}

///|
/// Parse a double from a substring field of an SDF line.
fn chem_parse_double_field(s : StringView) -> Double {
  let trimmed = chem_trim(s)
  match parse_double(trimmed) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Parse an integer from a substring field of an SDF line.
fn chem_parse_int_field(s : StringView) -> Int {
  let trimmed = chem_trim(s)
  match parse_double(trimmed) {
    Some(v) => v.to_int()
    None => 0
  }
}

///|
/// Parse SDF text into an array of molecules.
///
/// The SDF format:
/// ```
/// molecule_name
///   timestamp
/// comment
/// counts_line  (aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv)
/// atom_block   (aaa lines)
/// bond_block   (bbb lines)
/// > <property_name>
/// property_value
/// (blank line)
/// $$$$
/// ```
pub fn chem_parse_sdf(content : String) -> Array[SdfMolecule] {
  let molecules = Array::new()
  let lines = content.split("\n").to_array()
  let n_lines = lines.length()
  let mut idx = 0
  while idx < n_lines {
    // Skip blank lines between molecules.
    while idx < n_lines && chem_trim(lines[idx]).length() == 0 {
      idx = idx + 1
    }
    if idx >= n_lines {
      break
    }
    // Line 1: molecule name.
    let name = chem_trim(lines[idx])
    idx = idx + 1
    // Line 2: timestamp (skip).
    idx = idx + 1
    // Line 3: comment (skip).
    idx = idx + 1
    if idx >= n_lines {
      break
    }
    // Line 4: counts line.
    let counts_line = lines[idx]
    idx = idx + 1
    let n_atoms = chem_parse_int_field(counts_line[0 : 3])
    let n_bonds = chem_parse_int_field(counts_line[3 : 6])
    // Parse atom block.
    let atoms : Array[SdfAtom] = Array::new()
    for _ in 0..<n_atoms {
      if idx >= n_lines {
        break
      }
      let line = lines[idx]
      idx = idx + 1
      // SDF atom line: xxxx.xxxxyyyy.yyyyzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
      let x = chem_parse_double_field(line[0 : 10])
      let y = chem_parse_double_field(line[10 : 20])
      let z = chem_parse_double_field(line[20 : 30])
      let symbol = chem_trim(line[31 : 34])
      // Charge field is at columns 36-38 (0-based from 36).
      let charge = if line.length() >= 39 {
        let c = chem_parse_int_field(line[36 : 39])
        // SDF charge encoding: 0=0, 1=+3, 2=+2, 3=+1, 4=0, 5=-1, 6=-2, 7=-3
        match c {
          1 => 3
          2 => 2
          3 => 1
          5 => -1
          6 => -2
          7 => -3
          _ => 0
        }
      } else {
        0
      }
      atoms.push(SdfAtom::new(symbol, x, y, z, charge))
    }
    // Parse bond block.
    let bonds : Array[SdfBond] = Array::new()
    for _ in 0..<n_bonds {
      if idx >= n_lines {
        break
      }
      let line = lines[idx]
      idx = idx + 1
      let a1 = chem_parse_int_field(line[0 : 3]) - 1 // convert 1-based to 0-based
      let a2 = chem_parse_int_field(line[3 : 6]) - 1
      let bt = chem_parse_int_field(line[6 : 9])
      bonds.push(SdfBond::new(a1, a2, bt))
    }
    let mol = SdfMolecule::new(name, atoms, bonds)
    // Parse properties until $$$$
    while idx < n_lines {
      let line = lines[idx]
      let trimmed = chem_trim(line)
      if trimmed == "$$$$" {
        idx = idx + 1
        break
      }
      // Property line: "> <name>" (or ">  <name>" with extra spaces).
      if trimmed.length() > 3 && trimmed[0 : 1] == ">" {
        // Find the first '<' and the last '>'.
        let lt_pos = trimmed.find("<")
        // Find the last '>' by scanning from the end.
        let mut gt_idx = -1
        let tn = trimmed.length()
        let mut k = tn - 1
        while k >= 0 {
          if trimmed[k : k + 1] == ">" {
            gt_idx = k
            break
          }
          k = k - 1
        }
        let prop_name = match lt_pos {
          Some(lt) => {
            if gt_idx > lt + 1 {
              trimmed[lt + 1 : gt_idx].to_owned()
            } else {
              ""
            }
          }
          None => ""
        }
        idx = idx + 1
        if idx < n_lines && prop_name.length() > 0 {
          let value = chem_trim(lines[idx])
          mol.set_property(prop_name, value)
          idx = idx + 1
        }
      } else {
        idx = idx + 1
      }
    }
    molecules.push(mol)
  }
  molecules
}

// ============================================================================
// SDF serialization
// ============================================================================

///|
/// Format a double into a fixed-width 10-character field (right-aligned).
fn chem_format_coord(v : Double) -> String {
  let s = v.to_string()
  // Pad to 10 characters, right-aligned.
  if s.length() >= 10 {
    s
  } else {
    let pad = 10 - s.length()
    let mut result = ""
    for _ in 0..<pad {
      result = result + " "
    }
    result + s
  }
}

///|
/// Format an integer into a fixed-width 3-character field.
fn chem_format_int3(v : Int) -> String {
  let s = v.to_string()
  if s.length() >= 3 {
    s
  } else {
    let pad = 3 - s.length()
    let mut result = ""
    for _ in 0..<pad {
      result = result + " "
    }
    result + s
  }
}

///|
/// Serialize a single molecule to SDF text.
pub fn chem_molecule_to_sdf(mol : SdfMolecule) -> String {
  let sb = StringBuilder::new()
  // Header: name, timestamp, comment.
  sb.write_string(mol.name)
  sb.write_char('\n')
  sb.write_string("  -ChemmineR-MoonBit\n")
  sb.write_string("\n")
  // Counts line.
  sb.write_string(chem_format_int3(mol.n_atoms()))
  sb.write_string(chem_format_int3(mol.n_bonds()))
  sb.write_string("  0  0  0  0  0  0  0  0  0 V2000\n")
  // Atom block.
  for a in mol.atoms {
    sb.write_string(chem_format_coord(a.x))
    sb.write_string(chem_format_coord(a.y))
    sb.write_string(chem_format_coord(a.z))
    sb.write_string(" ")
    let sym = a.symbol
    // Pad symbol to 3 chars.
    if sym.length() >= 3 {
      sb.write_string(sym)
    } else {
      sb.write_string(sym)
      let pad = 3 - sym.length()
      for _ in 0..<pad {
        sb.write_char(' ')
      }
    }
    sb.write_string(" 0  0  0  0  0  0  0  0  0  0  0  0\n")
  }
  // Bond block.
  for b in mol.bonds {
    sb.write_string(chem_format_int3(b.atom1 + 1))
    sb.write_string(chem_format_int3(b.atom2 + 1))
    sb.write_string(chem_format_int3(b.bond_type))
    sb.write_string("  0  0  0  0\n")
  }
  // Properties.
  for key in mol.properties.keys() {
    sb.write_string("> <")
    sb.write_string(key)
    sb.write_string(">\n")
    sb.write_string(mol.properties[key])
    sb.write_char('\n')
    sb.write_char('\n')
  }
  sb.write_string("M  END\n")
  sb.write_string("$$$$\n")
  sb.to_string()
}

///|
/// Serialize an array of molecules to SDF text.
pub fn chem_molecules_to_sdf(molecules : Array[SdfMolecule]) -> String {
  let sb = StringBuilder::new()
  for mol in molecules {
    sb.write_string(chem_molecule_to_sdf(mol))
  }
  sb.to_string()
}

// ============================================================================
// Summary and convenience
// ============================================================================

///|
/// A text summary of a molecule: name, formula, MW, atom/bond counts.
pub fn chem_molecule_summary(mol : SdfMolecule) -> String {
  let formula = chem_molecular_formula(mol)
  let mw = chem_sdf_molecular_weight(mol)
  let rings = chem_ring_count(mol)
  "Molecule: " +
  mol.name +
  "\n  Formula: " +
  formula +
  "\n  MW: " +
  mw.to_string() +
  "\n  Atoms: " +
  mol.n_atoms().to_string() +
  "\n  Bonds: " +
  mol.n_bonds().to_string() +
  "\n  Heavy atoms: " +
  chem_heavy_atom_count(mol).to_string() +
  "\n  Rings: " +
  rings.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Build a sample water molecule (H2O).
pub fn chem_sample_water() -> SdfMolecule {
  let atoms = [
    SdfAtom::new("O", 0.0, 0.0, 0.0, 0),
    SdfAtom::new("H", 0.96, 0.0, 0.0, 0),
    SdfAtom::new("H", -0.24, 0.93, 0.0, 0),
  ]
  let bonds = [SdfBond::new(0, 1, 1), SdfBond::new(0, 2, 1)]
  SdfMolecule::new("water", atoms, bonds)
}

///|
/// Build a sample methane molecule (CH4).
pub fn chem_sample_methane() -> SdfMolecule {
  let atoms = [
    SdfAtom::new("C", 0.0, 0.0, 0.0, 0),
    SdfAtom::new("H", 0.63, 0.63, 0.63, 0),
    SdfAtom::new("H", -0.63, -0.63, 0.63, 0),
    SdfAtom::new("H", -0.63, 0.63, -0.63, 0),
    SdfAtom::new("H", 0.63, -0.63, -0.63, 0),
  ]
  let bonds = [
    SdfBond::new(0, 1, 1),
    SdfBond::new(0, 2, 1),
    SdfBond::new(0, 3, 1),
    SdfBond::new(0, 4, 1),
  ]
  SdfMolecule::new("methane", atoms, bonds)
}

///|
/// Build a sample ethanol molecule (C2H6O).
pub fn chem_sample_ethanol() -> SdfMolecule {
  let atoms = [
    SdfAtom::new("C", 1.16, -0.09, 0.0, 0), // C1 (methyl)
    SdfAtom::new("C", 0.0, 0.59, 0.0, 0), // C2 (methylene)
    SdfAtom::new("O", -1.13, -0.12, 0.0, 0), // O
    SdfAtom::new("H", 2.09, 0.45, 0.0, 0), // H on C1
    SdfAtom::new("H", 1.18, -1.18, 0.0, 0), // H on C1
    SdfAtom::new("H", 2.09, -0.67, 0.89, 0), // H on C1
    SdfAtom::new("H", 0.02, 1.68, 0.0, 0), // H on C2
    SdfAtom::new("H", 0.02, 0.01, 0.89, 0), // H on C2
    SdfAtom::new("H", -1.98, 0.37, 0.0, 0), // H on O
  ]
  let bonds = [
    SdfBond::new(0, 1, 1), // C1-C2
    SdfBond::new(1, 2, 1), // C2-O
    SdfBond::new(0, 3, 1), // C1-H
    SdfBond::new(0, 4, 1), // C1-H
    SdfBond::new(0, 5, 1), // C1-H
    SdfBond::new(1, 6, 1), // C2-H
    SdfBond::new(1, 7, 1), // C2-H
    SdfBond::new(2, 8, 1), // O-H
  ]
  let mol = SdfMolecule::new("ethanol", atoms, bonds)
  mol.set_property("MW", "46.069")
  mol
}

///|
/// Build a sample benzene molecule (C6H6) for ring detection.
pub fn chem_sample_benzene() -> SdfMolecule {
  let atoms = [
    SdfAtom::new("C", 1.40, 0.0, 0.0, 0),
    SdfAtom::new("C", 0.70, 1.21, 0.0, 0),
    SdfAtom::new("C", -0.70, 1.21, 0.0, 0),
    SdfAtom::new("C", -1.40, 0.0, 0.0, 0),
    SdfAtom::new("C", -0.70, -1.21, 0.0, 0),
    SdfAtom::new("C", 0.70, -1.21, 0.0, 0),
    SdfAtom::new("H", 2.49, 0.0, 0.0, 0),
    SdfAtom::new("H", 1.25, 2.16, 0.0, 0),
    SdfAtom::new("H", -1.25, 2.16, 0.0, 0),
    SdfAtom::new("H", -2.49, 0.0, 0.0, 0),
    SdfAtom::new("H", -1.25, -2.16, 0.0, 0),
    SdfAtom::new("H", 1.25, -2.16, 0.0, 0),
  ]
  let bonds = [
    SdfBond::new(0, 1, 4), // aromatic
    SdfBond::new(1, 2, 4),
    SdfBond::new(2, 3, 4),
    SdfBond::new(3, 4, 4),
    SdfBond::new(4, 5, 4),
    SdfBond::new(5, 0, 4),
    SdfBond::new(0, 6, 1),
    SdfBond::new(1, 7, 1),
    SdfBond::new(2, 8, 1),
    SdfBond::new(3, 9, 1),
    SdfBond::new(4, 10, 1),
    SdfBond::new(5, 11, 1),
  ]
  SdfMolecule::new("benzene", atoms, bonds)
}

///|
/// Build a sample aspirin molecule (C9H8O4, acetylsalicylic acid).
pub fn chem_sample_aspirin() -> SdfMolecule {
  // Simplified 2D coordinates for aspirin.
  // Benzene ring: C1-C6, with C1 bearing COOH and C2 bearing O-CO-CH3.
  let atoms = [
    SdfAtom::new("C", 1.40, 0.0, 0.0, 0), // C1 (ring, attached to COOH)
    SdfAtom::new("C", 0.70, 1.21, 0.0, 0), // C2 (ring, attached to ester O)
    SdfAtom::new("C", -0.70, 1.21, 0.0, 0), // C3
    SdfAtom::new("C", -1.40, 0.0, 0.0, 0), // C4
    SdfAtom::new("C", -0.70, -1.21, 0.0, 0), // C5
    SdfAtom::new("C", 0.70, -1.21, 0.0, 0), // C6
    SdfAtom::new("C", 2.80, 0.0, 0.0, 0), // C7 (carboxyl C)
    SdfAtom::new("O", 3.50, 1.21, 0.0, 0), // O (carboxyl =O)
    SdfAtom::new("O", 3.50, -1.21, 0.0, 0), // O (carboxyl -OH)
    SdfAtom::new("H", 4.49, -1.21, 0.0, 0), // H (carboxyl OH)
    SdfAtom::new("O", 1.40, 2.49, 0.0, 0), // O (ester)
    SdfAtom::new("C", 2.80, 2.49, 0.0, 0), // C (acetyl C)
    SdfAtom::new("O", 3.50, 3.70, 0.0, 0), // O (acetyl =O)
    SdfAtom::new("C", 2.10, 3.70, 0.0, 0), // C (methyl)
    SdfAtom::new("H", -2.49, 0.0, 0.0, 0), // H on C4
    SdfAtom::new("H", -1.25, -2.16, 0.0, 0), // H on C5
    SdfAtom::new("H", 1.25, -2.16, 0.0, 0), // H on C6
    SdfAtom::new("H", -1.25, 2.16, 0.0, 0), // H on C3
    SdfAtom::new("H", 2.10, 4.79, 0.0, 0), // H on methyl
    SdfAtom::new("H", 1.01, 3.70, 0.0, 0), // H on methyl
    SdfAtom::new("H", 2.10, 4.79, 0.89, 0), // H on methyl
  ]
  let bonds = [
    SdfBond::new(0, 1, 4), // C1-C2 aromatic
    SdfBond::new(1, 2, 4), // C2-C3
    SdfBond::new(2, 3, 4), // C3-C4
    SdfBond::new(3, 4, 4), // C4-C5
    SdfBond::new(4, 5, 4), // C5-C6
    SdfBond::new(5, 0, 4), // C6-C1
    SdfBond::new(0, 6, 1), // C1-C7 (carboxyl)
    SdfBond::new(6, 7, 2), // C7=O
    SdfBond::new(6, 8, 1), // C7-OH
    SdfBond::new(8, 9, 1), // O-H
    SdfBond::new(1, 10, 1), // C2-O (ester)
    SdfBond::new(10, 11, 1), // O-C (acetyl)
    SdfBond::new(11, 12, 2), // C=O (acetyl)
    SdfBond::new(11, 13, 1), // C-CH3
    SdfBond::new(3, 14, 1), // C4-H
    SdfBond::new(4, 15, 1), // C5-H
    SdfBond::new(5, 16, 1), // C6-H
    SdfBond::new(2, 17, 1), // C3-H
    SdfBond::new(13, 18, 1), // CH3-H
    SdfBond::new(13, 19, 1), // CH3-H
    SdfBond::new(13, 20, 1), // CH3-H
  ]
  let mol = SdfMolecule::new("aspirin", atoms, bonds)
  mol.set_property("MW", "180.158")
  mol.set_property("DRUG", "true")
  mol
}

///|
/// A sample SDF text with two molecules (ethanol and methane) for round-trip
/// testing.
pub fn chem_sample_sdf_text() -> String {
  let mols = [chem_sample_ethanol(), chem_sample_methane()]
  chem_molecules_to_sdf(mols)
}
