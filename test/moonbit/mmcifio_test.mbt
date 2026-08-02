///|
/// Tests for Bio.PDB.mmcifio module (mmCIF writer).

///|
test "mmcifio_format_coord_positive" {
  let s = @src.mmcif_format_coord(1.0)
  assert_eq(s, "1.000")
  let s2 = @src.mmcif_format_coord(1.5)
  assert_eq(s2, "1.500")
  let s3 = @src.mmcif_format_coord(123.456789)
  assert_eq(s3, "123.456")
}

///|
test "mmcifio_format_coord_zero" {
  let s = @src.mmcif_format_coord(0.0)
  assert_eq(s, "0.000")
}

///|
test "mmcifio_format_coord_negative" {
  let s = @src.mmcif_format_coord(-1.5)
  assert_eq(s, "-1.500")
  let s2 = @src.mmcif_format_coord(-0.0)
  assert_eq(s2, "0.000")
}

///|
test "mmcifio_format_bfactor" {
  let s = @src.mmcif_format_bfactor(0.0)
  assert_eq(s, "0.00")
  let s2 = @src.mmcif_format_bfactor(20.0)
  assert_eq(s2, "20.00")
  let s3 = @src.mmcif_format_bfactor(99.999)
  assert_eq(s3, "99.99")
  let s4 = @src.mmcif_format_bfactor(1.5)
  assert_eq(s4, "1.50")
}

///|
test "mmcifio_escape_value_empty" {
  let s = @src.mmcif_escape_value("")
  assert_eq(s, "?")
}

///|
test "mmcifio_escape_value_plain" {
  let s = @src.mmcif_escape_value("ALA")
  assert_eq(s, "ALA")
  let s2 = @src.mmcif_escape_value("CA")
  assert_eq(s2, "CA")
  let s3 = @src.mmcif_escape_value("1ABC")
  assert_eq(s3, "1ABC")
}

///|
test "mmcifio_escape_value_with_space" {
  let s = @src.mmcif_escape_value("hello world")
  assert_eq(s, "'hello world'")
  let s2 = @src.mmcif_escape_value(" alpha")
  assert_eq(s2, "' alpha'")
}

///|
test "mmcifio_escape_value_with_quote" {
  let s = @src.mmcif_escape_value("it's")
  assert_eq(s, "'it''s'")
  let s2 = @src.mmcif_escape_value("'quoted'")
  assert_eq(s2, "'''quoted'''")
}

///|
test "mmcifio_format_atom_basic" {
  let atom = @src.Atom::new(
    name="N",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    occupancy=1.0,
    bfactor=0.0,
    element="N",
  )
  let line = @src.mmcif_format_atom(atom, 1)
  let parts = line.split("   ").to_array()
  // Expected columns (filtered for non-empty):
  // ATOM | 1 | N | N | . | ALA | A | 1 | 1 | ? | 1.000 | 2.000 | 3.000 | 1.00 | 0.00 | ? | 1 | ALA | A | N
  let cols : Array[String] = []
  for p in parts {
    let t = p.to_owned()
    if t.length() > 0 {
      cols.push(t)
    }
  }
  assert_eq(cols[0], "ATOM")
  assert_eq(cols[1], "1")
  assert_eq(cols[2], "N")
  assert_eq(cols[3], "N")
  assert_eq(cols[4], ".")
  assert_eq(cols[5], "ALA")
  assert_eq(cols[6], "A")
  assert_eq(cols[7], "1")
  assert_eq(cols[8], "1")
  assert_eq(cols[9], "?")
  assert_eq(cols[10], "1.000")
  assert_eq(cols[11], "2.000")
  assert_eq(cols[12], "3.000")
  assert_eq(cols[13], "1.00")
  assert_eq(cols[14], "0.00")
  assert_eq(cols[15], "?")
  assert_eq(cols[16], "1")
  assert_eq(cols[17], "ALA")
  assert_eq(cols[18], "A")
  assert_eq(cols[19], "N")
}

///|
test "mmcifio_format_atom_with_group_hetatm" {
  let atom = @src.Atom::new(
    name="FE",
    coord=@src.Vector3::new(0.0, 0.0, 0.0),
    resname="FE",
    chainid='A',
    resseq=1,
    element="FE",
    charge=2,
  )
  let line = @src.mmcif_format_atom_with_group(atom, 5, "HETATM")
  assert_eq(line.contains("HETATM"), true)
  assert_eq(line.contains("FE"), true)
  // Charge of +2 should appear as "2+"
  assert_eq(line.contains("2+"), true)
}

///|
test "mmcifio_format_atom_negative_charge" {
  let atom = @src.Atom::new(
    name="CL",
    coord=@src.Vector3::new(0.0, 0.0, 0.0),
    resname="CL",
    chainid='A',
    resseq=1,
    element="CL",
    charge=-1,
  )
  let line = @src.mmcif_format_atom(atom, 1)
  // Charge of -1 should appear as "1-"
  assert_eq(line.contains("1-"), true)
}

///|
test "mmcifio_header_basic" {
  let structure = @src.Structure::new(id="1ABC")
  let header = @src.write_mmcif_header(structure)
  assert_eq(header.contains("_entry.id   1ABC"), true)
  assert_eq(header.contains("_entity.id   1"), true)
  assert_eq(header.contains("_entity.type   polymer"), true)
  assert_eq(header.contains("_entity.src_method   man"), true)
}

///|
test "mmcifio_header_empty_id" {
  let structure = @src.Structure::new(id="")
  let header = @src.write_mmcif_header(structure)
  // Empty id should fall back to XXXX
  assert_eq(header.contains("_entry.id   XXXX"), true)
}

///|
test "mmcifio_atom_site_empty" {
  let structure = @src.Structure::new(id="1ABC")
  let site = @src.write_mmcif_atom_site(structure)
  assert_eq(site.contains("loop_"), true)
  assert_eq(site.contains("_atom_site.group_PDB"), true)
  assert_eq(site.contains("_atom_site.id"), true)
  assert_eq(site.contains("_atom_site.Cartn_x"), true)
  assert_eq(site.contains("_atom_site.auth_atom_id"), true)
}

///|
test "mmcifio_atom_site_single_atom" {
  let atom = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(1.5, 2.5, 3.5),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
  )
  let residue = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    atom,
  ])
  let chain = @src.Chain::new(id='A', residues=[residue])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1XXX", models=[model])
  let site = @src.write_mmcif_atom_site(structure)
  assert_eq(site.contains("ATOM"), true)
  assert_eq(site.contains("1.500"), true)
  assert_eq(site.contains("2.500"), true)
  assert_eq(site.contains("3.500"), true)
  assert_eq(site.contains("ALA"), true)
}

///|
test "mmcifio_write_empty_structure" {
  let structure = @src.Structure::new(id="1ABC")
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("data_1ABC"), true)
  assert_eq(output.contains("_entry.id   1ABC"), true)
  assert_eq(output.contains("loop_"), true)
  assert_eq(output.contains("_atom_site.group_PDB"), true)
}

///|
test "mmcifio_write_single_atom" {
  let atom = @src.Atom::new(
    name="N",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="N",
  )
  let residue = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    atom,
  ])
  let chain = @src.Chain::new(id='A', residues=[residue])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1ABC", models=[model])
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("data_1ABC"), true)
  assert_eq(output.contains("ATOM"), true)
  assert_eq(output.contains("1.000"), true)
  assert_eq(output.contains("2.000"), true)
  assert_eq(output.contains("3.000"), true)
}

///|
test "mmcifio_write_single_residue" {
  let atoms = [
    @src.Atom::new(
      name="N",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="N",
    ),
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(1.5, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
  ]
  let residue = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms~)
  let chain = @src.Chain::new(id='A', residues=[residue])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1ALA", models=[model])
  let output = @src.write_mmcif(structure)
  // Two atoms -> two ATOM lines in the atom_site loop
  let lines = output.split("\n").to_array()
  let mut atom_count = 0
  for line in lines {
    if line.contains("ATOM") && line.length() > 4 {
      atom_count = atom_count + 1
    }
  }
  assert_eq(atom_count, 2)
}

///|
test "mmcifio_write_single_chain" {
  let r1 = @src.Residue::new(
    resname="ALA",
    chainid='A',
    resseq=1,
    atoms=[
      @src.Atom::new(
        name="N",
        coord=@src.Vector3::new(0.0, 0.0, 0.0),
        resname="ALA",
        chainid='A',
        resseq=1,
        element="N",
      ),
    ],
  )
  let r2 = @src.Residue::new(
    resname="GLY",
    chainid='A',
    resseq=2,
    atoms=[
      @src.Atom::new(
        name="N",
        coord=@src.Vector3::new(3.0, 0.0, 0.0),
        resname="GLY",
        chainid='A',
        resseq=2,
        element="N",
      ),
    ],
  )
  let chain = @src.Chain::new(id='A', residues=[r1, r2])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1CG", models=[model])
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("ALA"), true)
  assert_eq(output.contains("GLY"), true)
}

///|
test "mmcifio_write_multiple_chains" {
  let r1 = @src.Residue::new(
    resname="ALA",
    chainid='A',
    resseq=1,
    atoms=[
      @src.Atom::new(
        name="N",
        coord=@src.Vector3::new(0.0, 0.0, 0.0),
        resname="ALA",
        chainid='A',
        resseq=1,
        element="N",
      ),
    ],
  )
  let r2 = @src.Residue::new(
    resname="GLY",
    chainid='B',
    resseq=1,
    atoms=[
      @src.Atom::new(
        name="N",
        coord=@src.Vector3::new(10.0, 0.0, 0.0),
        resname="GLY",
        chainid='B',
        resseq=1,
        element="N",
      ),
    ],
  )
  let chain_a = @src.Chain::new(id='A', residues=[r1])
  let chain_b = @src.Chain::new(id='B', residues=[r2])
  let model = @src.Model::new(id=1, chains=[chain_a, chain_b])
  let structure = @src.Structure::new(id="2CHN", models=[model])
  let output = @src.write_mmcif(structure)
  // Both chains should be present
  assert_eq(output.contains("ALA"), true)
  assert_eq(output.contains("GLY"), true)
  // Count atom rows: should be 2 (one per chain)
  let lines = output.split("\n").to_array()
  let mut atom_rows = 0
  for line in lines {
    if line.contains("ATOM") && line.contains("10.000") {
      atom_rows = atom_rows + 1
    }
  }
  // The second atom has x=10.000
  assert_eq(atom_rows, 1)
}

///|
test "mmcifio_write_multiple_models" {
  let atom1 = @src.Atom::new(
    name="N",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="N",
  )
  let residue1 = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    atom1,
  ])
  let chain1 = @src.Chain::new(id='A', residues=[residue1])
  let model1 = @src.Model::new(id=1, chains=[chain1])
  let atom2 = @src.Atom::new(
    name="N",
    coord=@src.Vector3::new(4.0, 5.0, 6.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="N",
  )
  let residue2 = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    atom2,
  ])
  let chain2 = @src.Chain::new(id='A', residues=[residue2])
  let model2 = @src.Model::new(id=2, chains=[chain2])
  let structure = @src.Structure::new(id="2MDL", models=[model1, model2])
  let output = @src.write_mmcif(structure)
  // Both atoms from both models should be present (sequential atom ids 1 and 2)
  assert_eq(output.contains("1.000"), true)
  assert_eq(output.contains("4.000"), true)
}

///|
test "mmcifio_write_hetatm" {
  let atom = @src.Atom::new(
    name="FE",
    coord=@src.Vector3::new(5.0, 5.0, 5.0),
    resname="FE",
    chainid='A',
    resseq=1,
    element="FE",
  )
  // hetfield="H" makes this a hetero-residue
  let residue = @src.Residue::new(
    resname="FE",
    chainid='A',
    resseq=1,
    atoms=[atom],
    hetfield="H",
  )
  let chain = @src.Chain::new(id='A', residues=[residue])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1HET", models=[model])
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("HETATM"), true)
  assert_eq(output.contains("5.000"), true)
}

///|
test "mmcifio_sample_structure" {
  let structure = @src.mmcif_sample_structure()
  assert_eq(structure.get_id(), "1ABC")
  assert_eq(structure.get_num_models(), 1)
  assert_eq(structure.get_num_chains(), 1)
  assert_eq(structure.get_num_residues(), 2)
  assert_eq(structure.get_num_atoms(), 8)
}

///|
test "mmcifio_sample_structure_output" {
  let structure = @src.mmcif_sample_structure()
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("data_1ABC"), true)
  assert_eq(output.contains("_entry.id   1ABC"), true)
  assert_eq(output.contains("ALA"), true)
  assert_eq(output.contains("GLY"), true)
  // First atom (N of ALA at 1.0, 2.0, 3.0)
  assert_eq(output.contains("1.000"), true)
  assert_eq(output.contains("2.000"), true)
  assert_eq(output.contains("3.000"), true)
}

///|
test "mmcifio_roundtrip_parse" {
  let structure = @src.mmcif_sample_structure()
  let output = @src.write_mmcif(structure)
  // Parse the written output back using the existing mmcif parser
  let block = @src.parse_mmcif(output)
  assert_eq(block.name, "1ABC")
  let atoms = block.to_atoms()
  // All 8 atoms should be parsed back
  assert_eq(atoms.length(), 8)
  // First atom should be N of ALA at (1.0, 2.0, 3.0)
  let first = atoms[0]
  assert_eq(first.name, "N")
  assert_eq(first.resname, "ALA")
  assert_eq(first.coord.x, 1.0)
  assert_eq(first.coord.y, 2.0)
  assert_eq(first.coord.z, 3.0)
}

///|
test "mmcifio_roundtrip_preserves_residues" {
  let structure = @src.mmcif_sample_structure()
  let output = @src.write_mmcif(structure)
  let block = @src.parse_mmcif(output)
  let atoms = block.to_atoms()
  // Count unique residues by (resname, resseq) pairs
  let mut ala_count = 0
  let mut gly_count = 0
  for atom in atoms {
    if atom.resname == "ALA" {
      ala_count = ala_count + 1
    }
    if atom.resname == "GLY" {
      gly_count = gly_count + 1
    }
  }
  assert_eq(ala_count, 4)
  assert_eq(gly_count, 4)
}

///|
test "mmcifio_sequential_atom_ids" {
  let structure = @src.mmcif_sample_structure()
  let output = @src.write_mmcif(structure)
  // The 8 atoms should have ids 1..8
  let lines = output.split("\n").to_array()
  let ids : Array[Int] = []
  for line in lines {
    if line.contains("ATOM") {
      // Split by "   " (3-space separator) and grab the second column (id)
      let parts = line.split("   ").to_array()
      let cols : Array[String] = []
      for p in parts {
        let t = p.to_owned()
        if t.length() > 0 {
          cols.push(t)
        }
      }
      if cols.length() >= 2 {
        let id = @src.parse_int(cols[1])
        ids.push(id)
      }
    }
  }
  assert_eq(ids.length(), 8)
  assert_eq(ids[0], 1)
  assert_eq(ids[7], 8)
}

///|
test "mmcifio_negative_coordinates" {
  let atom = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(-1.5, -2.5, -3.5),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
  )
  let residue = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    atom,
  ])
  let chain = @src.Chain::new(id='A', residues=[residue])
  let model = @src.Model::new(id=1, chains=[chain])
  let structure = @src.Structure::new(id="1NEG", models=[model])
  let output = @src.write_mmcif(structure)
  assert_eq(output.contains("-1.500"), true)
  assert_eq(output.contains("-2.500"), true)
  assert_eq(output.contains("-3.500"), true)
}

///|
test "mmcifio_altloc_and_icode" {
  // Build an atom with explicit altloc and icode via direct construction.
  let atom = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
  )
  let line = @src.mmcif_format_atom(atom, 1)
  // Default altloc=' ' -> label_alt_id is "."
  assert_eq(line.contains("   .   "), true)
  // Default icode=' ' -> pdbx_PDB_ins_code is "?"
  assert_eq(line.contains("   ?   "), true)
}

///|
test "mmcifio_atom_site_has_all_columns" {
  let structure = @src.mmcif_sample_structure()
  let site = @src.write_mmcif_atom_site(structure)
  // Verify all 20 column headers are present
  let columns = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_formal_charge",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
  ]
  for col in columns {
    assert_eq(site.contains(col), true)
  }
}

///|
test "mmcifio_full_output_structure" {
  let structure = @src.mmcif_sample_structure()
  let output = @src.write_mmcif(structure)
  let lines = output.split("\n").to_array()
  // First line should be data_1ABC
  let first_non_empty = lines[0].to_owned()
  assert_eq(first_non_empty, "data_1ABC")
  // Should contain header section
  let mut has_entry = false
  let mut has_entity = false
  let mut has_loop = false
  for line in lines {
    if line.contains("_entry.id") {
      has_entry = true
    }
    if line.contains("_entity.id") {
      has_entity = true
    }
    if line.to_owned() == "loop_" {
      has_loop = true
    }
  }
  assert_eq(has_entry, true)
  assert_eq(has_entity, true)
  assert_eq(has_loop, true)
}

///|
test "mmcifio_write_empty_id_fallback" {
  let structure = @src.Structure::new(id="")
  let output = @src.write_mmcif(structure)
  // Empty id should fall back to XXXX in the data_ line
  assert_eq(output.contains("data_XXXX"), true)
}

///|
test "mmcifio_format_atom_id_increments" {
  let atom = @src.Atom::new(
    name="N",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="N",
  )
  let line1 = @src.mmcif_format_atom(atom, 1)
  let line2 = @src.mmcif_format_atom(atom, 42)
  // The id column should reflect the passed atom_id
  let cols1 = line1.split("   ").to_array()
  let cols2 = line2.split("   ").to_array()
  let id1 = cols1[1].to_owned()
  let id2 = cols2[1].to_owned()
  assert_eq(id1, "1")
  assert_eq(id2, "42")
}

///|
test "mmcifio_large_structure_values" {
  // Edge case: large coordinate values
  let atom = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(999.999, -999.999, 0.001),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
    bfactor=99.99,
  )
  let line = @src.mmcif_format_atom(atom, 1)
  assert_eq(line.contains("999.999"), true)
  assert_eq(line.contains("-999.999"), true)
  assert_eq(line.contains("0.001"), true)
  assert_eq(line.contains("99.99"), true)
}
