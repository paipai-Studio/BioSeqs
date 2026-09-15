///|
/// ic_rebuild demo - Rebuild Cartesian from internal coordinates (NeRF).
/// Inspired by Biopython Bio.PDB.ic_rebuild.

fn main {
  println("=== ic_rebuild Demo: NeRF Algorithm ===")
  println("")

  // Step 1: Seed atoms
  println("--- Step 1: Seed Atoms ---")
  let p0 = @src.IcPoint::new(0.0, 0.0, 0.0)
  let p1 = @src.IcPoint::new(1.5, 0.0, 0.0)
  let p2 = @src.IcPoint::new(2.5, 1.0, 0.0)
  println("P0: (" + p0.x.to_string() + ", " + p0.y.to_string() + ", " + p0.z.to_string() + ")")
  println("P1: (" + p1.x.to_string() + ", " + p1.y.to_string() + ", " + p1.z.to_string() + ")")
  println("P2: (" + p2.x.to_string() + ", " + p2.y.to_string() + ", " + p2.z.to_string() + ")")
  println("")

  // Step 2: Internal coords
  println("--- Step 2: Internal Coordinates ---")
  let ics : Array[@src.IcInternalCoord] = [
    @src.IcInternalCoord::{ atom_name: "CA", bond_length: 1.46, bond_angle: 1.91, dihedral: -0.52 },
    @src.IcInternalCoord::{ atom_name: "CB", bond_length: 1.53, bond_angle: 1.91, dihedral: 1.57 },
    @src.IcInternalCoord::{ atom_name: "CG", bond_length: 1.50, bond_angle: 1.91, dihedral: -1.05 },
  ]
  let mut i = 0
  while i < ics.length() {
    let ic = ics[i]
    println(ic.atom_name + ": bl=" + ic.bond_length.to_string() + " ba=" + ic.bond_angle.to_string() + " dh=" + ic.dihedral.to_string())
    i = i + 1
  }
  println("")

  // Step 3: Rebuild
  println("--- Step 3: Rebuilt Cartesian Coordinates ---")
  let seeds : Array[@src.IcPoint] = [p0, p1, p2]
  let rebuilt = @src.ic_rebuild_chain(seeds, ics)
  let mut j = 0
  while j < rebuilt.length() {
    let atom = rebuilt[j]
    println(atom.name + ": (" + atom.pos.x.to_string() + ", " + atom.pos.y.to_string() + ", " + atom.pos.z.to_string() + ")")
    j = j + 1
  }
  println("")

  // Step 4: Round-trip verification
  println("--- Step 4: Round-trip Verification ---")
  let positions : Array[@src.IcPoint] = [
    @src.IcPoint::new(0.0, 0.0, 0.0),
    @src.IcPoint::new(1.5, 0.0, 0.0),
    @src.IcPoint::new(2.5, 1.0, 0.0),
    @src.IcPoint::new(3.0, 2.0, 0.5),
    @src.IcPoint::new(2.5, 3.0, 1.0),
    @src.IcPoint::new(1.5, 3.5, 0.8),
  ]
  let dev = @src.ic_verify_roundtrip(positions)
  println("Max deviation: " + dev.to_string())
  println("")

  println("=== ic_rebuild Demo Complete ===")
}