///|
test "pdb_header_minimal_creation" {
  let header = @src.PDBHeader::new()
  assert_eq(header.pdb_id, "")
  assert_eq(header.title, "")
  assert_eq(header.resolution, 0.0)
  assert_eq(header.rfactor, 0.0)
  assert_eq(header.keywords.length(), 0)
  assert_eq(header.authors.length(), 0)
  assert_eq(header.chain_ids.length(), 0)
}

///|
test "pdb_header_with_values" {
  let header = @src.PDBHeader::new(
    pdb_id="1ABC",
    title="Test Structure",
    resolution=2.5,
    rfactor=0.18,
  )
  assert_eq(header.pdb_id, "1ABC")
  assert_eq(header.title, "Test Structure")
  assert_eq(header.resolution, 2.5)
  assert_eq(header.rfactor, 0.18)
}

///|
test "compound_info_creation" {
  let ci = @src.CompoundInfo::new(
    molecule_id="1",
    chain_ids=["A"],
    compound_name="Test Protein",
    formula="C500H800N100O200",
    ec_number="1.2.3.4",
  )
  assert_eq(ci.molecule_id, "1")
  assert_eq(ci.chain_ids.length(), 1)
  assert_eq(ci.compound_name, "Test Protein")
  assert_eq(ci.formula, "C500H800N100O200")
  assert_eq(ci.ec_number, "1.2.3.4")
  assert_false(ci.engineered)
  assert_false(ci.mutation)
}

///|
test "source_info_creation" {
  let si = @src.SourceInfo::new(
    organism_scientific="Homo sapiens",
    organism_common="Human",
    tax_id="9606",
  )
  assert_eq(si.organism_scientific, "Homo sapiens")
  assert_eq(si.organism_common, "Human")
  assert_eq(si.tax_id, "9606")
  assert_eq(si.strain, "")
}

///|
test "dbref_entry_creation" {
  let entry = @src.DBRefEntry::new(
    db_name="PDB",
    db_accession="1ABC",
    db_id="1",
  )
  assert_eq(entry.db_name, "PDB")
  assert_eq(entry.db_accession, "1ABC")
  assert_eq(entry.db_id, "1")
}

///|
test "parse_title_record_single_line" {
  let lines = [
    "TITLE     Crystal structure of a protein",
  ]
  let title = @src.parse_title_record(lines, 0)
  assert_true(title.contains("Crystal structure"))
}

///|
test "parse_title_record_multi_line" {
  let lines = [
    "TITLE     Crystal structure of a protein",
    "TITLE     in complex with ligand",
  ]
  let title = @src.parse_title_record(lines, 0)
  assert_true(title.contains("Crystal structure"))
  assert_true(title.contains("in complex with ligand"))
}

///|
test "parse_title_record_missing" {
  let lines = [
    "HEADER                                              01-JAN-20   1ABC",
  ]
  let title = @src.parse_title_record(lines, 0)
  assert_eq(title, "")
}

///|
test "parse_resolution_record" {
  let line = "REMARK   2 RESOLUTION. 2.50 ANGSTROMS"
  let resolution = @src.parse_resolution_record(line)
  assert_true(resolution > 0.0)
}

///|
test "parse_resolution_record_no_keyword" {
  let line = "REMARK   1   SOME OTHER TEXT"
  let resolution = @src.parse_resolution_record(line)
  assert_eq(resolution, 0.0)
}

///|
test "parse_rfactor_record" {
  let line = "REMARK   3 R FACTOR. 0.18"
  let rfactor = @src.parse_rfactor_record(line)
  assert_true(rfactor > 0.0)
}

///|
test "parse_rfactor_record_no_keyword" {
  let line = "REMARK   2 RESOLUTION. 2.50"
  let rfactor = @src.parse_rfactor_record(line)
  assert_eq(rfactor, 0.0)
}

///|
test "extract_chain_ids" {
  let compound = @src.CompoundInfo::new(
    chain_ids=["A", "B", "C"],
  )
  let header = @src.PDBHeader::new(compound=compound)
  let chains = @src.extract_chain_ids(header)
  assert_eq(chains.length(), 3)
  assert_eq(chains[0], "A")
  assert_eq(chains[1], "B")
  assert_eq(chains[2], "C")
}

///|
test "extract_chain_ids_empty" {
  let header = @src.PDBHeader::new()
  let chains = @src.extract_chain_ids(header)
  assert_eq(chains.length(), 0)
}

///|
test "parse_pdb_header_empty" {
  let header = @src.parse_pdb_header("")
  assert_eq(header.pdb_id, "")
  assert_eq(header.title, "")
  assert_eq(header.resolution, 0.0)
  assert_eq(header.rfactor, 0.0)
}

///|
test "parse_pdb_header_missing_records" {
  let pdb_text =
    "HEADER   VIRUS 01-JAN-20                            1ABC\n"
  let header = @src.parse_pdb_header(pdb_text)
  assert_eq(header.title, "")
  assert_eq(header.compound.molecule_id, "")
  assert_eq(header.source.organism_scientific, "")
  assert_eq(header.keywords.length(), 0)
  assert_eq(header.resolution, 0.0)
  assert_eq(header.rfactor, 0.0)
  assert_eq(header.authors.length(), 0)
  assert_eq(header.journal, "")
  assert_eq(header.dbrefs.length(), 0)
}

///|
test "parse_pdb_header_lines" {
  let lines = [
    "HEADER   VIRUS 01-JAN-20                            1ABC",
    "TITLE     Crystal structure",
  ]
  let header = @src.parse_pdb_header_lines(lines)
  assert_true(header.title.length() > 0)
}

///|
test "parse_resolution_record_no_resolution" {
  let line = "REMARK   2   No resolution data"
  let resolution = @src.parse_resolution_record(line)
  assert_eq(resolution, 0.0)
}

///|
test "parse_rfactor_record_no_rfactor" {
  let line = "REMARK   3   No r-factor data"
  let rfactor = @src.parse_rfactor_record(line)
  assert_eq(rfactor, 0.0)
}

///|
test "compound_info_default_creation" {
  let ci = @src.CompoundInfo::new()
  assert_eq(ci.molecule_id, "")
  assert_eq(ci.chain_ids.length(), 0)
  assert_eq(ci.compound_name, "")
  assert_eq(ci.formula, "")
  assert_false(ci.engineered)
  assert_false(ci.mutation)
}

///|
test "source_info_default_creation" {
  let si = @src.SourceInfo::new()
  assert_eq(si.organism_scientific, "")
  assert_eq(si.organism_common, "")
  assert_eq(si.tax_id, "")
  assert_eq(si.strain, "")
}

///|
test "dbref_entry_default_creation" {
  let entry = @src.DBRefEntry::new()
  assert_eq(entry.db_name, "")
  assert_eq(entry.db_accession, "")
  assert_eq(entry.db_id, "")
}
