///|
/// Hydrogen bond detection — geometric donor-acceptor criteria.
/// Ported from biotite.structure.hbond.
///
/// Criteria (Baker & Hubbard 1984, as used in biotite):
///   D-A distance ≤ max_distance (default 3.5 Å)
///   H-A distance ≤ max_h_distance (default 2.5 Å)
///   D-H···A angle ≥ min_angle (default 120°)

///|
/// A detected hydrogen bond.
pub struct HBond {
  donor_resname : String
  donor_resseq : Int32
  donor_chainid : Char
  donor_atom : String // e.g. "N", "OG"
  acceptor_resname : String
  acceptor_resseq : Int32
  acceptor_chainid : Char
  acceptor_atom : String // e.g. "O", "OD1"
  distance : Double // D-A distance (Å)
  angle : Double // D-H-A angle (degrees)
}
