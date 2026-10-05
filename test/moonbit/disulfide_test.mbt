///|
/// Helper: create a cysteine residue with backbone + SG atoms.
fn make_cys(
  chainid : Char,
  resseq : Int,
  n_pos : @src.Vector3,
  ca_pos : @src.Vector3,
  cb_pos : @src.Vector3,
  sg_pos : @src.Vector3,
) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=n_pos, resname="CYS", chainid=chainid, resseq=resseq),
    @src.Atom::new(name="CA", coord=ca_pos, resname="CYS", chainid=chainid, resseq=resseq),
    @src.Atom::new(name="CB", coord=cb_pos, resname="CYS", chainid=chainid, resseq=resseq),
    @src.Atom::new(name="SG", coord=sg_pos, resname="CYS", chainid=chainid, resseq=resseq),
  ]
  @src.Residue::new(resname="CYS", chainid=chainid, resseq=resseq, atoms~)
}

///|
/// Helper: create a non-cysteine residue.
fn make_gly(
  chainid : Char,
  resseq : Int,
  ca_pos : @src.Vector3,
) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(ca_pos.x - 1.5, ca_pos.y, ca_pos.z), resname="GLY", chainid=chainid, resseq=resseq),
    @src.Atom::new(name="CA", coord=ca_pos, resname="GLY", chainid=chainid, resseq=resseq),
    @src.Atom::new(name="C", coord=@src.Vector3::new(ca_pos.x + 1.5, ca_pos.y, ca_pos.z), resname="GLY", chainid=chainid, resseq=resseq),
  ]
  @src.Residue::new(resname="GLY", chainid=chainid, resseq=resseq, atoms~)
}

///|
test "disulfide: find disulfide bond" {
  // Two cysteines with SG atoms 2.05 Å apart
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(6.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(4.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5 + 2.05, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  // SG-SG distance is 2.05, should be found with default cutoff 2.5
  assert_eq(bonds.length(), 1)
}

///|
test "disulfide: no disulfide (distant Cys)" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(20.0, 0.0, 0.0),
    @src.Vector3::new(18.5, 0.0, 0.0),
    @src.Vector3::new(17.0, 1.0, 0.0),
    @src.Vector3::new(16.0, 1.5, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 0)
}

///|
test "disulfide: no cysteines" {
  let gly1 = make_gly('A', 1, @src.Vector3::new(0.0, 0.0, 0.0))
  let gly2 = make_gly('A', 2, @src.Vector3::new(3.8, 0.0, 0.0))
  let chain = @src.Chain::new(id='A', residues=[gly1, gly2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 0)
}

///|
test "disulfide: single cysteine" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 0)
}

///|
test "disulfide: oxidation state" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(6.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(4.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5 + 2.05, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 1)
  assert_eq(bonds[0].oxidation_state, "oxidized")
}

///|
test "disulfide: sg_sg_distance" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(10.0, 0.0, 0.0),
    @src.Vector3::new(8.5, 0.0, 0.0),
    @src.Vector3::new(7.0, 1.0, 0.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 1)
  assert_true((bonds[0].sg_sg_distance - 2.0).abs() < 0.01)
}

///|
test "disulfide: distance cutoff filtering" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(10.0, 0.0, 0.0),
    @src.Vector3::new(8.5, 0.0, 0.0),
    @src.Vector3::new(7.0, 1.0, 0.0),
    @src.Vector3::new(2.3, 0.0, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  // SG-SG = 2.3, found with cutoff 2.5
  let bonds1 = @src.find_disulfide_bonds(chain)
  assert_eq(bonds1.length(), 1)
  // Not found with cutoff 2.0
  let bonds2 = @src.find_disulfide_bonds(chain, distance_cutoff=2.0)
  assert_eq(bonds2.length(), 0)
}

///|
test "disulfide: summary output" {
  let bonds : Array[@src.DisulfideBond] = []
  let summary = @src.disulfide_bond_summary(bonds)
  assert_true(summary.contains("No disulfide"))
}

///|
test "disulfide: chi angles present" {
  let cys1 = make_cys(
    'A', 1,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let cys2 = make_cys(
    'A', 2,
    @src.Vector3::new(6.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(4.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5 + 2.05, 0.0),
  )
  let chain = @src.Chain::new(id='A', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 1)
  // Chi angles should be finite numbers (not NaN)
  let bond = bonds[0]
  assert_true(bond.chi1.abs() < 180.1)
  assert_true(bond.chi2.abs() < 180.1)
  assert_true(bond.chi_ss.abs() < 180.1)
}

///|
test "disulfide: resseq and chainid" {
  let cys1 = make_cys(
    'B', 10,
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5, 0.0),
  )
  let cys2 = make_cys(
    'B', 20,
    @src.Vector3::new(6.0, 0.0, 0.0),
    @src.Vector3::new(5.0, 0.0, 0.0),
    @src.Vector3::new(4.0, 1.0, 0.0),
    @src.Vector3::new(3.0, 1.5 + 2.05, 0.0),
  )
  let chain = @src.Chain::new(id='B', residues=[cys1, cys2])
  let bonds = @src.find_disulfide_bonds(chain)
  assert_eq(bonds.length(), 1)
  assert_eq(bonds[0].chainid1, 'B')
  assert_eq(bonds[0].resseq1, 10)
  assert_eq(bonds[0].chainid2, 'B')
  assert_eq(bonds[0].resseq2, 20)
}

///|
test "disulfide: format_disulfide_bond" {
  let bond = @src.DisulfideBond::{
    chainid1: 'A',
    resseq1: 5,
    chainid2: 'A',
    resseq2: 15,
    sg_sg_distance: 2.05,
    chi1: -60.0,
    chi2: 120.0,
    chi_ss: -90.0,
    chi2p: -120.0,
    chi1p: 60.0,
    oxidation_state: "oxidized",
  }
  let s = @src.format_disulfide_bond(bond)
  assert_true(s.contains("SS-bond"))
  assert_true(s.contains("oxidized"))
}
