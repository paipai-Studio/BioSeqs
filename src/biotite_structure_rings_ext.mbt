///|
/// biotite_structure_rings_ext.mbt — ring orders and histograms on top of
/// the SSSR/aromaticity detector in structure_rings.mbt.

///|
/// For every atom, the size of the smallest ring containing it (0 = acyclic).
pub fn biotite_get_ring_orders(
  atoms : Array[Atom],
  max_size~ : Int32 = 7,
) -> Array[Int32] {
  let order : Array[Int32] = Array::make(atoms.length(), 0)
  let rings = detect_rings(atoms, max_size~)
  for ring in rings {
    for idx in ring.atom_indices {
      if idx >= 0 && idx < order.length() {
        if order[idx] == 0 || ring.size < order[idx] {
          order[idx] = ring.size
        }
      }
    }
  }
  order
}

///|
/// Histogram indexed by ring size (index 0 unused).
pub fn biotite_ring_histogram(
  atoms : Array[Atom],
  max_size~ : Int32 = 7,
) -> Array[Int32] {
  let histogram : Array[Int32] = Array::make(max_size + 1, 0)
  let rings = detect_rings(atoms, max_size~)
  for ring in rings {
    if ring.size >= 0 && ring.size <= max_size {
      histogram[ring.size] = histogram[ring.size] + 1
    }
  }
  histogram
}

///|
/// Number of rings classified as aromatic by the underlying detector.
pub fn biotite_count_aromatic_rings(
  atoms : Array[Atom],
  max_size~ : Int32 = 7,
) -> Int32 {
  let rings = detect_rings(atoms, max_size~)
  let mut count = 0
  for ring in rings {
    if ring.is_aromatic {
      count = count + 1
    }
  }
  count
}

///|
/// Atom indices belonging to at least one detected ring.
pub fn biotite_ring_atoms(
  atoms : Array[Atom],
  max_size~ : Int32 = 7,
) -> Array[Int32] {
  let orders = biotite_get_ring_orders(atoms, max_size~)
  let result : Array[Int32] = []
  for i in 0..<orders.length() {
    if orders[i] > 0 {
      result.push(i)
    }
  }
  result
}
