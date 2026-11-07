///|
/// Temporary debug output (to be deleted before commit).

///|
test "zz_debug2" {
  let pts = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ]
  let ca : Array[@src.AtomCoordinate] = []
  let cb : Array[@src.AtomCoordinate] = []
  for p in pts {
    let _ = ca.push(@src.AtomCoordinate::new(p.0, p.1, p.2))
    let _ = cb.push(@src.AtomCoordinate::new(p.0, p.1, p.2))
  }
  let sup = @src.SVDSuperimposer::set(ca, cb)
  println("svd_same_rmsd=" + sup.get_rmsd().to_string())
  let atoms = ["C", "C", "C", "C"]
  let planar = @src.RdkitConformer::from_coords(atoms, pts)
  println("pbf_planar=" + @src.rdkit_pbf(planar).to_string())
  let buckled = @src.RdkitConformer::from_coords(
    atoms,
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.6), (1.0, 1.0, 0.0), (0.0, 1.0, -0.6)],
  )
  println("pbf_buckled=" + @src.rdkit_pbf(buckled).to_string())
}
