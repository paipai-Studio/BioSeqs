///|
/// Tests for protein_interactions module (InSty).

///|
/// Build a minimal ARG residue with charged atoms.
fn make_arg(resseq : Int, nz_x : Double, nz_y : Double, nz_z : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ARG", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="ARG", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(2.0, 1.0, 0.0), resname="ARG", resseq=resseq),
    @src.Atom::new(name="NE", coord=@src.Vector3::new(nz_x, nz_y, nz_z), resname="ARG", resseq=resseq),
    @src.Atom::new(name="NH1", coord=@src.Vector3::new(nz_x + 0.5, nz_y, nz_z), resname="ARG", resseq=resseq),
    @src.Atom::new(name="NH2", coord=@src.Vector3::new(nz_x - 0.5, nz_y, nz_z), resname="ARG", resseq=resseq),
  ]
  @src.Residue::new(resname="ARG", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a minimal ASP residue with charged atoms.
fn make_asp(resseq : Int, od_x : Double, od_y : Double, od_z : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ASP", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="ASP", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(2.0, 1.0, 0.0), resname="ASP", resseq=resseq),
    @src.Atom::new(name="OD1", coord=@src.Vector3::new(od_x, od_y, od_z), resname="ASP", resseq=resseq),
    @src.Atom::new(name="OD2", coord=@src.Vector3::new(od_x + 0.5, od_y, od_z), resname="ASP", resseq=resseq),
  ]
  @src.Residue::new(resname="ASP", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a minimal LYS residue.
fn make_lys(resseq : Int, nz_x : Double, nz_y : Double, nz_z : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="LYS", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="LYS", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(2.0, 1.0, 0.0), resname="LYS", resseq=resseq),
    @src.Atom::new(name="NZ", coord=@src.Vector3::new(nz_x, nz_y, nz_z), resname="LYS", resseq=resseq),
  ]
  @src.Residue::new(resname="LYS", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a minimal PHE residue with ring atoms in xy-plane.
fn make_phe(resseq : Int, cx : Double, cy : Double, cz : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(cx, cy, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(cx + 1.5, cy, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(cx + 2.0, cy + 1.0, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CG", coord=@src.Vector3::new(cx + 2.0, cy + 2.0, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CD1", coord=@src.Vector3::new(cx + 2.7, cy + 2.7, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CD2", coord=@src.Vector3::new(cx + 1.3, cy + 2.7, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CE1", coord=@src.Vector3::new(cx + 2.7, cy + 4.0, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CE2", coord=@src.Vector3::new(cx + 1.3, cy + 4.0, cz), resname="PHE", resseq=resseq),
    @src.Atom::new(name="CZ", coord=@src.Vector3::new(cx + 2.0, cy + 4.5, cz), resname="PHE", resseq=resseq),
  ]
  @src.Residue::new(resname="PHE", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a minimal VAL residue.
fn make_val(resseq : Int, cx : Double, cy : Double, cz : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(cx, cy, cz), resname="VAL", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(cx + 1.5, cy, cz), resname="VAL", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(cx + 2.0, cy + 1.0, cz), resname="VAL", resseq=resseq),
    @src.Atom::new(name="CG1", coord=@src.Vector3::new(cx + 2.5, cy + 2.0, cz), resname="VAL", resseq=resseq),
    @src.Atom::new(name="CG2", coord=@src.Vector3::new(cx + 1.5, cy + 2.0, cz), resname="VAL", resseq=resseq),
  ]
  @src.Residue::new(resname="VAL", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a minimal LEU residue.
fn make_leu(resseq : Int, cx : Double, cy : Double, cz : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(cx, cy, cz), resname="LEU", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(cx + 1.5, cy, cz), resname="LEU", resseq=resseq),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(cx + 2.0, cy + 1.0, cz), resname="LEU", resseq=resseq),
    @src.Atom::new(name="CG", coord=@src.Vector3::new(cx + 2.5, cy + 2.0, cz), resname="LEU", resseq=resseq),
    @src.Atom::new(name="CD1", coord=@src.Vector3::new(cx + 3.0, cy + 3.0, cz), resname="LEU", resseq=resseq),
    @src.Atom::new(name="CD2", coord=@src.Vector3::new(cx + 2.0, cy + 3.0, cz), resname="LEU", resseq=resseq),
  ]
  @src.Residue::new(resname="LEU", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a water molecule.
fn make_water(resseq : Int, x : Double, y : Double, z : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="O", coord=@src.Vector3::new(x, y, z), resname="HOH", resseq=resseq),
  ]
  @src.Residue::new(resname="HOH", chainid='A', resseq=resseq, atoms=atoms)
}

///|
/// Build a GLY residue with N and O atoms.
fn make_gly_full(resseq : Int, n_x : Double, n_y : Double, n_z : Double, o_x : Double, o_y : Double, o_z : Double) -> @src.Residue {
  let atoms = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(n_x, n_y, n_z), resname="GLY", resseq=resseq),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(n_x + 1.5, n_y, n_z), resname="GLY", resseq=resseq),
    @src.Atom::new(name="C", coord=@src.Vector3::new(n_x + 2.5, n_y, n_z), resname="GLY", resseq=resseq),
    @src.Atom::new(name="O", coord=@src.Vector3::new(o_x, o_y, o_z), resname="GLY", resseq=resseq),
  ]
  @src.Residue::new(resname="GLY", chainid='A', resseq=resseq, atoms=atoms)
}

// ─── Salt Bridge Tests ─────────────────────────────────────

///|
test "insty: salt bridge Arg-Asp close" {
  let arg = make_arg(1, 5.0, 0.0, 0.0)
  let asp = make_asp(2, 7.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let bridges = @src.find_salt_bridges(chain)
  assert_true(bridges.length() > 0)
  assert_true(bridges[0].resname1 == "ARG" || bridges[0].resname1 == "ASP")
}

///|
test "insty: salt bridge Arg-Asp far" {
  let arg = make_arg(1, 0.0, 0.0, 0.0)
  let asp = make_asp(2, 20.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let bridges = @src.find_salt_bridges(chain)
  assert_eq(bridges.length(), 0)
}

///|
test "insty: salt bridge cutoff" {
  let arg = make_arg(1, 0.0, 0.0, 0.0)
  let asp = make_asp(2, 5.5, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let bridges_default = @src.find_salt_bridges(chain)
  assert_eq(bridges_default.length(), 0)
  let bridges_wide = @src.find_salt_bridges(chain, distance_cutoff=7.0)
  assert_true(bridges_wide.length() > 0)
}

// ─── π-π Stacking Tests ────────────────────────────────────

///|
test "insty: pi-stacking parallel Phe-Phe" {
  let phe1 = make_phe(1, 0.0, 0.0, 0.0)
  let phe2 = make_phe(2, 0.0, 0.0, 4.0)
  let chain = @src.Chain::new(id='A', residues=[phe1, phe2])
  let stackings = @src.find_pi_stacking(chain)
  assert_true(stackings.length() > 0)
  assert_true(stackings[0].centroid_distance < 6.0)
}

///|
test "insty: pi-stacking far apart" {
  let phe1 = make_phe(1, 0.0, 0.0, 0.0)
  let phe2 = make_phe(2, 0.0, 0.0, 20.0)
  let chain = @src.Chain::new(id='A', residues=[phe1, phe2])
  let stackings = @src.find_pi_stacking(chain)
  assert_eq(stackings.length(), 0)
}

// ─── π-Cation Tests ────────────────────────────────────────

///|
test "insty: pi-cation Lys-Phe" {
  let phe = make_phe(1, 0.0, 0.0, 0.0)
  let lys = make_lys(2, 2.0, 3.3, 3.0)
  let chain = @src.Chain::new(id='A', residues=[phe, lys])
  let cations = @src.find_pi_cation(chain)
  assert_true(cations.length() > 0)
  assert_true(cations[0].cation_resname == "LYS")
  assert_true(cations[0].aromatic_resname == "PHE")
}

///|
test "insty: pi-cation far" {
  let phe = make_phe(1, 0.0, 0.0, 0.0)
  let lys = make_lys(2, 50.0, 50.0, 50.0)
  let chain = @src.Chain::new(id='A', residues=[phe, lys])
  let cations = @src.find_pi_cation(chain)
  assert_eq(cations.length(), 0)
}

// ─── Hydrophobic Tests ─────────────────────────────────────

///|
test "insty: hydrophobic Val-Leu close" {
  let val = make_val(1, 0.0, 0.0, 0.0)
  let leu = make_leu(2, 3.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[val, leu])
  let interactions = @src.find_hydrophobic_interactions(chain)
  assert_true(interactions.length() > 0)
}

///|
test "insty: hydrophobic far" {
  let val = make_val(1, 0.0, 0.0, 0.0)
  let leu = make_leu(2, 20.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[val, leu])
  let interactions = @src.find_hydrophobic_interactions(chain)
  assert_eq(interactions.length(), 0)
}

///|
test "insty: hydrophobic non-hydrophobic excluded" {
  let arg = make_arg(1, 0.0, 0.0, 0.0)
  let asp = make_asp(2, 2.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let interactions = @src.find_hydrophobic_interactions(chain)
  assert_eq(interactions.length(), 0)
}

// ─── Repulsive Ionic Tests ─────────────────────────────────

///|
test "insty: repulsive Lys-Lys close" {
  let lys1 = make_lys(1, 0.0, 0.0, 0.0)
  let lys2 = make_lys(2, 3.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[lys1, lys2])
  let repulsions = @src.find_repulsive_ionic(chain)
  assert_true(repulsions.length() > 0)
}

///|
test "insty: repulsive Asp-Asp close" {
  let asp1 = make_asp(1, 0.0, 0.0, 0.0)
  let asp2 = make_asp(2, 3.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[asp1, asp2])
  let repulsions = @src.find_repulsive_ionic(chain)
  assert_true(repulsions.length() > 0)
}

///|
test "insty: repulsive far" {
  let lys1 = make_lys(1, 0.0, 0.0, 0.0)
  let lys2 = make_lys(2, 20.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[lys1, lys2])
  let repulsions = @src.find_repulsive_ionic(chain)
  assert_eq(repulsions.length(), 0)
}

// ─── Water Bridge Tests ────────────────────────────────────

///|
test "insty: water bridge" {
  // Two GLY residues with O atoms near a water
  let gly1 = make_gly_full(1, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0)
  let gly2 = make_gly_full(2, 6.0, 0.0, 0.0, 6.0, 0.0, 0.0)
  let water = make_water(100, 4.5, 0.0, 0.0)  // O at 4.5
  let chain = @src.Chain::new(id='A', residues=[gly1, gly2, water])
  let bridges = @src.find_water_bridges(chain)
  // Water O at 4.5, gly1 O at 3.0 → d=1.5; gly2 O at 6.0 → d=1.5
  assert_true(bridges.length() > 0)
}

///|
test "insty: water bridge none" {
  let gly1 = make_gly_full(1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  let gly2 = make_gly_full(2, 20.0, 0.0, 0.0, 20.0, 0.0, 0.0)
  let water = make_water(100, 50.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[gly1, gly2, water])
  let bridges = @src.find_water_bridges(chain)
  assert_eq(bridges.length(), 0)
}

// ─── Analyze All Tests ─────────────────────────────────────

///|
test "insty: analyze_interactions" {
  let arg = make_arg(1, 5.0, 0.0, 0.0)
  let asp = make_asp(2, 7.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let result = @src.analyze_interactions(chain)
  assert_true(result.salt_bridges.length() > 0)
}

///|
test "insty: empty chain" {
  let chain = @src.Chain::new(id='A', residues=[])
  let result = @src.analyze_interactions(chain)
  assert_eq(result.salt_bridges.length(), 0)
  assert_eq(result.pi_stackings.length(), 0)
  assert_eq(result.pi_cations.length(), 0)
  assert_eq(result.hydrophobic.length(), 0)
  assert_eq(result.repulsive_ionic.length(), 0)
  assert_eq(result.water_bridges.length(), 0)
}

///|
test "insty: summary output" {
  let arg = make_arg(1, 5.0, 0.0, 0.0)
  let asp = make_asp(2, 7.0, 0.0, 0.0)
  let chain = @src.Chain::new(id='A', residues=[arg, asp])
  let result = @src.analyze_interactions(chain)
  let summary = @src.interaction_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Salt Bridges"))
}
