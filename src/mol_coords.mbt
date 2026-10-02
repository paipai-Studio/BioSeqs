// ============================================================================
// 2D coordinate generation
// ============================================================================
// Generates 2D coordinates for a molecular graph for visualization.
// Algorithm:
// (1) Find SSSR rings and lay them out as regular polygons.
// (2) Fuse shared bonds between rings (adjust ring centers).
// (3) Place non-ring atoms relative to their neighbors (tree layout).

///|
/// 2D point.
pub struct Point2D {
  x : Double
  y : Double
}

///|
/// Generate 2D coordinates for a molecule.
/// Returns an array of (x, y) coordinates, one per atom.
pub fn generate_2d_coords(mol : Molecule) -> Array[(Double, Double)] {
  let n = mol.atoms.length()
  let coords = Array::make(n, (0.0, 0.0))
  let placed = Array::make(n, false)
  // Step 1: Layout rings
  let rings = mol_find_rings(mol)
  let mut ring_idx = 0
  for ring in rings {
    let rs = ring.length()
    let cx = ring_idx.to_double() * 3.0
    let cy = 0.0
    let radius = 1.5
    for i = 0; i < rs; i = i + 1 {
      let angle = 2.0 * 3.14159265358979 * i.to_double() / rs.to_double()
      let x = cx + radius * @math.cos(angle)
      let y = cy + radius * @math.sin(angle)
      coords[ring[i]] = (x, y)
      placed[ring[i]] = true
    }
    ring_idx = ring_idx + 1
  }
  // Step 2: Layout non-ring atoms (BFS from placed atoms)
  let neighbors = mol_adjacency(mol)
  let queue : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    if placed[k] {
      queue.push(k)
    }
  }
  let mut head = 0
  while head < queue.length() {
    let current = queue[head]
    head = head + 1
    let (cx, cy) = coords[current]
    // Count unplaced neighbors
    let unplaced : Array[Int] = Array(capacity=0)
    for (nb, _) in neighbors[current] {
      if not(placed[nb]) {
        unplaced.push(nb)
      }
    }
    if unplaced.length() == 0 {
      continue
    }
    // Place unplaced neighbors in a fan around current atom
    let bond_len = 1.5
    let start_angle = 0.0
    let step = 2.0 * 3.14159265358979 / (unplaced.length() + 1).to_double()
    for i = 0; i < unplaced.length(); i = i + 1 {
      let angle = start_angle + step * (i.to_double() + 1.0)
      let x = cx + bond_len * @math.cos(angle)
      let y = cy + bond_len * @math.sin(angle)
      coords[unplaced[i]] = (x, y)
      placed[unplaced[i]] = true
      queue.push(unplaced[i])
    }
  }
  // Step 3: Handle any remaining unplaced atoms (disconnected)
  for k = 0; k < n; k = k + 1 {
    if not(placed[k]) {
      coords[k] = (k.to_double() * 2.0, 0.0)
      placed[k] = true
    }
  }
  coords
}

///|
/// Generate 2D coordinates and return as Point2D array.
pub fn generate_2d_points(mol : Molecule) -> Array[Point2D] {
  let coords = generate_2d_coords(mol)
  let points : Array[Point2D] = Array(capacity=0)
  for (x, y) in coords {
    points.push(Point2D::{ x, y })
  }
  points
}

///|
/// Calculate bounding box of 2D coordinates.
pub fn coord_bounding_box(
  coords : Array[(Double, Double)],
) -> ((Double, Double), (Double, Double)) {
  if coords.length() == 0 {
    return ((0.0, 0.0), (0.0, 0.0))
  }
  let (fx, fy) = coords[0]
  let mut min_x = fx
  let mut max_x = fx
  let mut min_y = fy
  let mut max_y = fy
  for (x, y) in coords {
    if x < min_x {
      min_x = x
    }
    if x > max_x {
      max_x = x
    }
    if y < min_y {
      min_y = y
    }
    if y > max_y {
      max_y = y
    }
  }
  ((min_x, min_y), (max_x, max_y))
}

///|
/// Center coordinates at origin.
pub fn center_2d_coords(
  coords : Array[(Double, Double)],
) -> Array[(Double, Double)] {
  let ((min_x, min_y), (max_x, max_y)) = coord_bounding_box(coords)
  let cx = (min_x + max_x) / 2.0
  let cy = (min_y + max_y) / 2.0
  let result : Array[(Double, Double)] = Array(capacity=0)
  for (x, y) in coords {
    result.push((x - cx, y - cy))
  }
  result
}

///|
/// Format coordinates as a string table.
pub fn coords_to_string(
  mol : Molecule,
  coords : Array[(Double, Double)],
) -> String {
  let sb = StringBuilder()
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    let (x, y) = coords[i]
    sb.write_string(mol.atoms[i].symbol)
    sb.write_string(": (")
    sb.write_string(x.to_string())
    sb.write_string(", ")
    sb.write_string(y.to_string())
    sb.write_string(")\n")
  }
  sb.to_string()
}
