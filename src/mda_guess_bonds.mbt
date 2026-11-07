///|
/// Covalent bond guessing from pure geometry, mirroring MDAnalysis
/// `topology.guessers.guess_bonds`: for every atom pair `i < j` a bond is
/// created when the distance between the atoms is at most the sum of their
/// covalent radii plus a tolerance, while being larger than a small lower
/// cutoff (0.4 A) that removes spurious overlaps. An optional orthorhombic
/// box enables the minimum image convention so that bonds across periodic
/// boundaries are detected.

///|
/// Covalent radius (in Angstrom) for common elements, used by
/// `mda_guess_bonds`. Unknown or unsupported elements fall back to a
/// generic carbon-like radius of 0.77 A.
pub fn mda_cov_radius(element : String) -> Double {
  match element {
    "H" => 0.31
    "C" => 0.76
    "N" => 0.71
    "O" => 0.66
    "F" => 0.57
    "P" => 1.07
    "S" => 1.05
    "Cl" => 1.02
    "Na" => 1.66
    "Mg" => 1.41
    "K" => 2.03
    "Ca" => 1.76
    "Fe" => 1.32
    "Zn" => 1.22
    "Br" => 1.20
    "I" => 1.39
    _ => 0.77
  }
}

///|
/// Displacement vector `a - b` folded into the primary cell by the minimum
/// image convention for an orthorhombic box with edge lengths
/// `box = [Lx, Ly, Lz]`. Each component is wrapped as
/// `d - L * round(d / L)`, yielding the shortest periodic image.
fn mda_guess_min_image(
  a : Array[Double],
  b : Array[Double],
  box : Array[Double],
) -> Array[Double] {
  let mut dx = a[0] - b[0]
  let mut dy = a[1] - b[1]
  let mut dz = a[2] - b[2]
  dx = dx - box[0] * (dx / box[0]).round()
  dy = dy - box[1] * (dy / box[1]).round()
  dz = dz - box[2] * (dz / box[2]).round()
  [dx, dy, dz]
}

///|
/// Guess covalent bonds from atom coordinates and element symbols
/// (MDAnalysis `guess_bonds`). Two atoms `i < j` are bonded when
/// `dist(i, j) <= cov_radius(elem_i) + cov_radius(elem_j) + tolerance`
/// and `dist(i, j) > 0.4` A. When `box` is given (orthorhombic edge
/// lengths `[Lx, Ly, Lz]`) distances use the minimum image convention so
/// that bonds across periodic boundaries are found.
///
/// coords: atom coordinates `[atom][xyz]`
/// elements: element symbol per atom (same length as `coords`)
/// tolerance: extra distance allowed on top of the summed covalent radii
/// box: optional periodic box edge lengths
///
/// Returns the bond list as `(i, j)` index pairs with `i < j`, ordered by
/// the O(n^2) pair scan.
pub fn mda_guess_bonds(
  coords : Array[Array[Double]],
  elements : Array[String],
  tolerance~ : Double = 0.45,
  box~ : Array[Double]? = None,
) -> Array[(Int32, Int32)] {
  let n = coords.length()
  let bonds : Array[(Int32, Int32)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let cutoff = mda_cov_radius(elements[i]) +
        mda_cov_radius(elements[j]) +
        tolerance
      let d = match box {
        Some(bx) => {
          let v = mda_guess_min_image(coords[i], coords[j], bx)
          (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
        }
        None => mda_euclid(coords[i], coords[j])
      }
      if d > 0.4 && d <= cutoff {
        bonds.push((i, j))
      }
    }
  }
  bonds
}

///|
/// Build an adjacency list from a bond list: entry `i` holds the sorted-by-
/// insertion order indices of all atoms bonded to atom `i`. Isolated atoms
/// map to empty lists; bonds referencing out-of-range atoms are ignored.
pub fn mda_bond_list_adjacency(
  n_atoms : Int32,
  bonds : Array[(Int32, Int32)],
) -> Array[Array[Int32]] {
  let adj : Array[Array[Int32]] = Array(capacity=n_atoms)
  for i = 0; i < n_atoms; i = i + 1 {
    adj.push([])
  }
  for bond in bonds {
    let (i, j) = bond
    if i >= 0 && i < n_atoms && j >= 0 && j < n_atoms {
      adj[i].push(j)
      adj[j].push(i)
    }
  }
  adj
}
