// ===========================================================================
// component_separation_test.mbt — Tests for ligand/water/ion separation
// ===========================================================================

test "is_water_resname recognizes common water names" {
  assert_true(@src.is_water_resname("HOH"))
  assert_true(@src.is_water_resname("WAT"))
  assert_true(@src.is_water_resname("TIP"))
  assert_true(@src.is_water_resname("TIP3"))
  assert_true(@src.is_water_resname("TIP4"))
  assert_true(@src.is_water_resname("SOL"))
  assert_true(@src.is_water_resname("DOD"))
  assert_true(@src.is_water_resname("H2O"))
  assert_false(@src.is_water_resname("ALA"))
  assert_false(@src.is_water_resname("NA"))
  assert_false(@src.is_water_resname("HEM"))
}

test "is_ion_resname recognizes metal cations" {
  assert_true(@src.is_ion_resname("NA"))
  assert_true(@src.is_ion_resname("K"))
  assert_true(@src.is_ion_resname("MG"))
  assert_true(@src.is_ion_resname("CA"))
  assert_true(@src.is_ion_resname("ZN"))
  assert_true(@src.is_ion_resname("FE"))
  assert_true(@src.is_ion_resname("CU"))
  assert_true(@src.is_ion_resname("MN"))
  assert_false(@src.is_ion_resname("ALA"))
  assert_false(@src.is_ion_resname("HOH"))
}

test "is_ion_resname recognizes halide anions" {
  assert_true(@src.is_ion_resname("CL"))
  assert_true(@src.is_ion_resname("BR"))
  assert_true(@src.is_ion_resname("I"))
  assert_true(@src.is_ion_resname("F"))
  assert_false(@src.is_ion_resname("CLN"))
}

test "is_ion_resname handles charge notation" {
  assert_true(@src.is_ion_resname("NA+"))
  assert_true(@src.is_ion_resname("K+"))
  assert_true(@src.is_ion_resname("MG2"))
  assert_true(@src.is_ion_resname("CA2"))
  assert_true(@src.is_ion_resname("FE3"))
  assert_true(@src.is_ion_resname("CL-"))
}

test "classify_residue identifies water" {
  let atom = @src.Atom::new(name="O", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HOH")
  let res = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[atom], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Water)
}

test "classify_residue identifies ion by resname" {
  let atom = @src.Atom::new(name="NA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="NA", element="NA")
  let res = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[atom], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ion)
}

test "classify_residue identifies zinc ion" {
  let atom = @src.Atom::new(name="ZN", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ZN", element="ZN")
  let res = @src.Residue::new(resname="ZN", chainid='A', resseq=300, atoms=[atom], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ion)
}

test "classify_residue identifies chloride ion" {
  let atom = @src.Atom::new(name="CL", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="CL", element="CL")
  let res = @src.Residue::new(resname="CL", chainid='A', resseq=301, atoms=[atom], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ion)
}

test "classify_residue identifies ion by element fallback" {
  let atom = @src.Atom::new(name="MG", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="XXX", element="MG")
  let res = @src.Residue::new(resname="XXX", chainid='A', resseq=200, atoms=[atom], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ion)
}

test "classify_residue identifies protein" {
  let n = @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA")
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="ALA")
  let c = @src.Atom::new(name="C", coord=@src.Vector3::new(2.5, 0.0, 0.0), resname="ALA")
  let o = @src.Atom::new(name="O", coord=@src.Vector3::new(3.0, 0.0, 0.0), resname="ALA")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[n, ca, c, o])
  assert_true(@src.classify_residue(res) is @src.Protein)
}

test "classify_residue identifies all 20 standard AAs" {
  let aa_names = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
  for name in aa_names {
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname=name)
    let res = @src.Residue::new(resname=name, chainid='A', resseq=1, atoms=[ca])
    assert_true(@src.classify_residue(res) is @src.Protein)
  }
}

test "classify_residue identifies nucleic acid" {
  let p = @src.Atom::new(name="P", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="DA")
  let res = @src.Residue::new(resname="DA", chainid='B', resseq=1, atoms=[p])
  assert_true(@src.classify_residue(res) is @src.Nucleic)
}

test "classify_residue identifies ligand HEM" {
  let fe = @src.Atom::new(name="FE", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HEM")
  let c1 = @src.Atom::new(name="C1A", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="HEM")
  let c2 = @src.Atom::new(name="C2A", coord=@src.Vector3::new(2.0, 0.0, 0.0), resname="HEM")
  let res = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[fe, c1, c2], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ligand)
}

test "classify_residue identifies ATP ligand" {
  let p1 = @src.Atom::new(name="PG", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ATP")
  let p2 = @src.Atom::new(name="PB", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="ATP")
  let res = @src.Residue::new(resname="ATP", chainid='A', resseq=600, atoms=[p1, p2], hetfield="HETATM")
  assert_true(@src.classify_residue(res) is @src.Ligand)
}

test "extract_waters from chain" {
  let protein_res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let water1 = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(5.0, 5.0, 5.0), resname="HOH"),
  ], hetfield="HETATM")
  let water2 = @src.Residue::new(resname="WAT", chainid='A', resseq=101, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(6.0, 5.0, 5.0), resname="WAT"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[protein_res, water1, water2])
  let waters = @src.extract_waters(chain)
  assert_eq(waters.length(), 2)
}

test "extract_ions from chain" {
  let protein_res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let na = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let zn = @src.Residue::new(resname="ZN", chainid='A', resseq=201, atoms=[
    @src.Atom::new(name="ZN", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="ZN", element="ZN"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[protein_res, na, zn])
  let ions = @src.extract_ions(chain)
  assert_eq(ions.length(), 2)
}

test "extract_ligands from chain" {
  let protein_res = @src.Residue::new(resname="GLY", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="GLY"),
  ])
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[protein_res, hem])
  let ligands = @src.extract_ligands(chain)
  assert_eq(ligands.length(), 1)
  assert_eq(ligands[0].resname, "HEM")
}

test "extract_protein from chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let gly = @src.Residue::new(resname="GLY", chainid='A', resseq=2, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(3.8, 0.0, 0.0), resname="GLY"),
  ])
  let water = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(5.0, 5.0, 5.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, gly, water])
  let protein = @src.extract_protein(chain)
  assert_eq(protein.length(), 2)
}

test "separate_chain full separation" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let hoh = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let na = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(7.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hem, hoh, na])
  let result = @src.separate_chain(chain)
  assert_eq(result.protein.length(), 1)
  assert_eq(result.ligands.length(), 1)
  assert_eq(result.waters.length(), 1)
  assert_eq(result.ions.length(), 1)
  assert_eq(result.total_residues(), 4)
}

test "separate_chain with nucleic acid" {
  let da = @src.Residue::new(resname="DA", chainid='B', resseq=1, atoms=[
    @src.Atom::new(name="P", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="DA"),
  ])
  let chain = @src.Chain::new(id='B', residues=[da])
  let result = @src.separate_chain(chain)
  assert_eq(result.nucleic.length(), 1)
  assert_eq(result.protein.length(), 0)
}

test "SeparationResult summary" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
    @src.Atom::new(name="N", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="ALA"),
  ])
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C2A", coord=@src.Vector3::new(7.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let hoh = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hem, hoh])
  let result = @src.separate_chain(chain)
  let s = result.summary()
  assert_eq(s.num_protein, 1)
  assert_eq(s.num_ligands, 1)
  assert_eq(s.num_waters, 1)
  assert_eq(s.num_ligand_atoms, 3)
  assert_eq(s.num_water_atoms, 1)
}

test "water_chain builds water-only chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hoh = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hoh])
  let wc = @src.water_chain(chain)
  assert_eq(wc.get_id(), 'A')
  assert_eq(wc.get_num_residues(), 1)
}

test "protein_chain builds protein-only chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hoh = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hoh])
  let pc = @src.protein_chain(chain)
  assert_eq(pc.get_num_residues(), 1)
}

test "separate_model across multiple chains" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hoh_a = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain_a = @src.Chain::new(id='A', residues=[ala, hoh_a])
  let gly = @src.Residue::new(resname="GLY", chainid='B', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="GLY"),
  ])
  let na = @src.Residue::new(resname="NA", chainid='B', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let chain_b = @src.Chain::new(id='B', residues=[gly, na])
  let model = @src.Model::new(id=1, chains=[chain_a, chain_b])
  let result = @src.separate_model(model)
  assert_eq(result.protein.length(), 2)
  assert_eq(result.waters.length(), 1)
  assert_eq(result.ions.length(), 1)
}

test "group_ligands groups multi-residue ligand" {
  // Two consecutive residues of same ligand (e.g., a large cofactor)
  let r1 = @src.Residue::new(resname="COA", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="COA"),
  ], hetfield="HETATM")
  let r2 = @src.Residue::new(resname="COA", chainid='A', resseq=501, atoms=[
    @src.Atom::new(name="C2", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="COA"),
  ], hetfield="HETATM")
  let r3 = @src.Residue::new(resname="HEM", chainid='A', resseq=600, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([r1, r2, r3])
  assert_eq(groups.length(), 2)
  // COA should have 2 residues, HEM should have 1
  let mut coa_found = false
  let mut hem_found = false
  for g in groups {
    if g.resname == "COA" {
      assert_eq(g.residues.length(), 2)
      assert_eq(g.num_atoms, 2)
      coa_found = true
    }
    if g.resname == "HEM" {
      assert_eq(g.residues.length(), 1)
      hem_found = true
    }
  }
  assert_true(coa_found)
  assert_true(hem_found)
}

test "group_ligands empty input" {
  let groups = @src.group_ligands([])
  assert_eq(groups.length(), 0)
}

test "group_ligands single ligand" {
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([hem])
  assert_eq(groups.length(), 1)
  assert_eq(groups[0].num_atoms, 2)
}

test "LigandMolecule center_of_mass" {
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(2.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([hem])
  let com = groups[0].center_of_mass()
  assert_eq(com.x, 1.0)
  assert_eq(com.y, 0.0)
  assert_eq(com.z, 0.0)
}

test "LigandMolecule get_atoms" {
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C2A", coord=@src.Vector3::new(2.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([hem])
  let atoms = groups[0].get_atoms()
  assert_eq(atoms.length(), 3)
}

test "contacting_ligands finds nearby ligands" {
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let far_lig = @src.Residue::new(resname="ATP", chainid='A', resseq=600, atoms=[
    @src.Atom::new(name="PG", coord=@src.Vector3::new(50.0, 50.0, 50.0), resname="ATP"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([hem, far_lig])
  let protein_atoms = [
    @src.Atom::new(name="CA", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="ALA"),
  ]
  let contacts = @src.contacting_ligands(groups, protein_atoms, 5.0)
  assert_eq(contacts.length(), 1)
  assert_eq(contacts[0], 0)
}

test "contacting_ligands finds none far away" {
  let far_lig = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(50.0, 50.0, 50.0), resname="HEM"),
  ], hetfield="HETATM")
  let groups = @src.group_ligands([far_lig])
  let protein_atoms = [
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ]
  let contacts = @src.contacting_ligands(groups, protein_atoms, 5.0)
  assert_eq(contacts.length(), 0)
}

test "ion_coordination finds coordinating atoms" {
  let na = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let protein_atoms = [
    @src.Atom::new(name="OD1", coord=@src.Vector3::new(5.5, 0.0, 0.0), resname="ASP"),
    @src.Atom::new(name="OD2", coord=@src.Vector3::new(5.0, 0.5, 0.0), resname="ASP"),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(20.0, 0.0, 0.0), resname="ALA"),
  ]
  let coord = @src.ion_coordination(na, protein_atoms, 3.0)
  assert_eq(coord.length(), 2)
}

test "ion_coordination no contacts" {
  let na = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(50.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let protein_atoms = [
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ]
  let coord = @src.ion_coordination(na, protein_atoms, 3.0)
  assert_eq(coord.length(), 0)
}

test "ion_chain builds ion-only chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let na = @src.Residue::new(resname="NA", chainid='A', resseq=200, atoms=[
    @src.Atom::new(name="NA", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="NA", element="NA"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, na])
  let ic = @src.ion_chain(chain)
  assert_eq(ic.get_num_residues(), 1)
}

test "ligand_chain builds ligand-only chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hem = @src.Residue::new(resname="HEM", chainid='A', resseq=500, atoms=[
    @src.Atom::new(name="FE", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="HEM"),
    @src.Atom::new(name="C1A", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="HEM"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hem])
  let lc = @src.ligand_chain(chain)
  assert_eq(lc.get_num_residues(), 1)
}

test "nucleic_chain builds nucleic-only chain" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let da = @src.Residue::new(resname="DA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="P", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="DA"),
  ])
  let chain = @src.Chain::new(id='A', residues=[ala, da])
  let nc = @src.nucleic_chain(chain)
  assert_eq(nc.get_num_residues(), 1)
}

test "SeparationSummary to_string produces output" {
  let ala = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA"),
  ])
  let hoh = @src.Residue::new(resname="HOH", chainid='A', resseq=100, atoms=[
    @src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 10.0), resname="HOH"),
  ], hetfield="HETATM")
  let chain = @src.Chain::new(id='A', residues=[ala, hoh])
  let result = @src.separate_chain(chain)
  let s = result.summary()
  let str = s.to_string()
  assert_true(str.contains("Protein"))
  assert_true(str.contains("Water"))
}
