///|
/// NeighborSearch example: KD-tree based spatial search for protein structures.
///
/// Run: moon run examples/neighbor_search_demo/main.mbt
fn main {
  println("=== Neighbor Search Demo (KD-Tree) ===")
  println("")

  // === Create sample atoms ===
  println("--- 1. Create Sample Atoms ---")

  let atoms = [
    @bio.Atom::new(name="CA", coord=@bio.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
    @bio.Atom::new(name="CB", coord=@bio.Vector3::new(1.5, 0.0, 0.0), resname="ALA"),
    @bio.Atom::new(name="CG", coord=@bio.Vector3::new(2.5, 0.5, 0.0), resname="VAL"),
    @bio.Atom::new(name="CA", coord=@bio.Vector3::new(5.0, 0.0, 0.0), resname="GLY"),
    @bio.Atom::new(name="CB", coord=@bio.Vector3::new(6.0, 0.5, 0.0), resname="GLY"),
    @bio.Atom::new(name="CA", coord=@bio.Vector3::new(10.0, 0.0, 0.0), resname="SER"),
    @bio.Atom::new(name="OG", coord=@bio.Vector3::new(10.5, 1.0, 0.0), resname="SER"),
    @bio.Atom::new(name="CA", coord=@bio.Vector3::new(0.0, 5.0, 0.0), resname="ASP"),
    @bio.Atom::new(name="CA", coord=@bio.Vector3::new(0.0, 0.0, 5.0), resname="LYS"),
  ]

  println("  Created \{atoms.length()} atoms")
  for i = 0; i < atoms.length(); i = i + 1 {
    let a = atoms[i]
    let c = a.coord
    println("    Atom \{i}: name=\{a.name}, res=\{a.resname}, pos=(\{c.x}, \{c.y}, \{c.z})")
  }
  println("")

  // === Build KD-tree ===
  println("--- 2. Build Neighbor Search (KD-tree) ---")

  let ns = @bio.NeighborSearch::new(atoms)
  println("  KD-tree built with \{ns.size()} atoms")
  println("")

  // === Radius search ===
  println("--- 3. Radius Search ---")

  let center1 = @bio.Vector3::new(0.0, 0.0, 0.0)
  let radius = 3.0
  let neighbors1 = ns.search(center1, radius)

  println("  Center: (\{center1.x}, \{center1.y}, \{center1.z})")
  println("  Radius: \{radius} A")
  println("  Found \{neighbors1.length()} neighbors:")
  for a in neighbors1 {
    let d = a.coord.distance(center1)
    println("    - \{a.name} (\{a.resname}): \{d} A")
  }
  println("")

  // === Another radius search ===
  let center2 = @bio.Vector3::new(5.5, 0.0, 0.0)
  let radius2 = 2.0
  let neighbors2 = ns.search(center2, radius2)

  println("  Center: (\{center2.x}, \{center2.y}, \{center2.z})")
  println("  Radius: \{radius2} A")
  println("  Found \{neighbors2.length()} neighbors")
  println("")

  // === Nearest neighbor ===
  println("--- 4. Nearest Neighbor Search ---")

  let query = @bio.Vector3::new(3.0, 0.0, 0.0)
  let nearest = ns.nearest(query)

  println("  Query point: (\{query.x}, \{query.y}, \{query.z})")
  match nearest {
    Some(a) => {
      let d = a.coord.distance(query)
      println("  Nearest: \{a.name} (\{a.resname}), distance: \{d} A")
    }
    None => println("  No atoms found")
  }
  println("")

  // === Atom-based radius search ===
  println("--- 5. Atom-Based Radius Search ---")

  let atom0 = atoms[0]
  let atom_radius = 4.0
  let atom_neighbors = ns.search_atom(atom0, atom_radius)

  println("  Query atom: \{atom0.name} (\{atom0.resname})")
  println("  Radius: \{atom_radius} A")
  println("  Found \{atom_neighbors.length()} neighboring atoms")
  println("")

  // === Search all pairs within distance ===
  println("--- 6. All Atom Pairs Within Distance ---")

  let pair_radius = 2.5
  let pairs = ns.search_all(pair_radius)

  println("  Distance cutoff: \{pair_radius} A")
  println("  Found \{pairs.length()} atom pairs:")
  for i = 0; i < pairs.length(); i = i + 1 {
    let (a, b) = pairs[i]
    let d = a.distance(b)
    println("    \{i+1}. \{a.name}-\{b.name}: \{d} A")
  }
  println("")

  println("=== Demo Complete ===")
}
