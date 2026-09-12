///|
/// Tests for the Bio.SCOP.Raf (ASTRAL RAF sequence maps) and Bio.SCOP.Dom
/// ports. All data is synthetic and mirrors the fixed-width RAF format:
/// header of 38 columns followed by 7-column residue fields.

///|
/// Build a synthetic RAF record line for the given residues.
/// Each residue is (resid, atom-letter, seqres-letter).
fn raf_line(
  pdbid : String,
  chainid : String,
  version : String,
  flags : String,
  residues : Array[(String, String, String)],
) -> String {
  let mut line = pdbid +
    chainid +
    " " +
    version +
    "    " +
    "26JUN9" +
    " " +
    flags
  // pad the header to 38 columns
  while line.length() < 38 {
    line = line + " "
  }
  for r in residues {
    let (resid, atom, seqres) = r
    // resid right-justified in 5 columns
    let mut field = resid
    while field.length() < 5 {
      field = " " + field
    }
    line = line + field + atom + seqres
  }
  line
}

///|
test "scop_raf_seqmap_parse_fields" {
  let line = raf_line("1sct", "A", "0.01", "150212", [
    ("1", "A", "A"),
    ("2", "G", "G"),
    ("10A", "K", "K"),
    ("11", "R", "R"),
  ])
  let sm = @src.scop_raf_seqmap_parse(line)
  assert_eq(sm.pdbid, "1sct")
  assert_eq(sm.pdb_datestamp, "26JUN9")
  assert_eq(sm.version, "0.01")
  assert_eq(sm.flags, "150212")
  assert_eq(sm.length(), 4)
  let r0 = sm.get(0)
  assert_eq(r0.chainid, "A")
  assert_eq(r0.resid, "1")
  assert_eq(r0.atom, "A")
  assert_eq(r0.seqres, "A")
  let r2 = sm.get(2)
  assert_eq(r2.resid, "10A")
  assert_eq(r2.atom, "K")
  assert_eq(r2.seqres, "K")
}

///|
test "scop_raf_seqmap_normalize_unknown" {
  // "." means unknown in RAF and is normalized to X; letters are uppercased
  let line = raf_line("1xyz", "_", "0.02", "      ", [
    ("1", ".", "a"),
    ("2", "x", "X"),
  ])
  let sm = @src.scop_raf_seqmap_parse(line)
  assert_eq(sm.version, "0.02")
  assert_eq(sm.get(0).atom, "X")
  assert_eq(sm.get(0).seqres, "A")
  assert_eq(sm.get(1).atom, "X")
  // blank chain id keeps the underscore
  assert_eq(sm.get(0).chainid, "_")
}

///|
test "scop_raf_seqmap_version_check" {
  let r = try {
    let _ = @src.scop_raf_seqmap_parse(
      raf_line("1sct", "A", "0.03", "", [("1", "A", "A")]),
    )
    false
  } catch {
    ScopRafError(msg) => msg == "Incompatible RAF version: 0.03"
  }
  assert_true(r)
}

///|
test "scop_raf_seqmap_index_and_slice" {
  let line = raf_line("1sct", "A", "0.01", "150212", [
    ("1", "A", "A"),
    ("2", "G", "G"),
    ("10A", "K", "K"),
    ("11", "R", "R"),
  ])
  let sm = @src.scop_raf_seqmap_parse(line)
  // Biopython's SeqMap.index defaults to the blank "_" chain, so the
  // chain id must be given explicitly for chain "A" residues.
  assert_eq(sm.index("1", chainid="A"), 0)
  assert_eq(sm.index("10A", chainid="A"), 2)
  assert_eq(sm.index("11", chainid="A"), 3)
  let r = try {
    let _ = sm.index("99", chainid="A")
    false
  } catch {
    ScopRafError(msg) => msg == "No such residue A99"
  }
  assert_true(r)
  let frag = sm.slice(1, 3)
  assert_eq(frag.length(), 2)
  assert_eq(frag.get(0).resid, "2")
  assert_eq(frag.get(1).resid, "10A")
}

///|
test "scop_raf_seqmap_extend_and_flags" {
  let a = @src.scop_raf_seqmap_parse(
    raf_line("1sct", "A", "0.01", "150212", [("1", "A", "A")]),
  )
  let b = @src.scop_raf_seqmap_parse(
    raf_line("1sct", "B", "0.01", "150212", [("1", "V", "V")]),
  )
  // matching everything: flags kept
  a.extend(b)
  assert_eq(a.length(), 2)
  assert_eq(a.get(1).chainid, "B")
  assert_eq(a.flags, "150212")
  // inconsistent flags are erased
  let c = @src.scop_raf_seqmap_parse(
    raf_line("1sct", "C", "0.01", "      ", [("2", "L", "L")]),
  )
  a.extend(c)
  assert_eq(a.length(), 3)
  assert_eq(a.flags, "")
  // different proteins are rejected
  let d = @src.scop_raf_seqmap_parse(
    raf_line("1xyz", "A", "0.01", "      ", [("1", "A", "A")]),
  )
  let r = try {
    a.extend(d)
    false
  } catch {
    ScopRafError(msg) => msg == "Cannot add fragments from different proteins"
  }
  assert_true(r)
}

///|
test "scop_raf_index_and_get_seqmap" {
  let raf_text = raf_line("1sct", "A", "0.01", "150212", [
      ("1", "A", "A"),
      ("2", "G", "G"),
      ("10A", "K", "K"),
      ("11", "R", "R"),
    ]) +
    "\n" +
    raf_line("1flp", "_", "0.01", "150212", [("5", "D", "D"), ("6", "E", "E")])
  let index = @src.scop_raf_index(raf_text)
  // fragment cut: residues 2 to 10A of chain A
  let residues = @src.scop_parse_residues("1sct A:2-10A")
  let sm = @src.scop_raf_get_seqmap(index, residues)
  assert_eq(sm.length(), 2)
  assert_eq(sm.get(0).resid, "2")
  assert_eq(sm.get(1).resid, "10A")
  // whole chain
  let whole = @src.scop_raf_get_seqmap(
    index,
    @src.scop_parse_residues("1sct A:"),
  )
  assert_eq(whole.length(), 4)
  // unnamed chain uses the underscore key
  let unnamed = @src.scop_raf_get_seqmap(
    index,
    @src.scop_parse_residues("1flp -"),
  )
  assert_eq(unnamed.length(), 2)
  assert_eq(unnamed.get(0).resid, "5")
  // unknown key raises
  let r = try {
    let _ = @src.scop_raf_get_seqmap(index, @src.scop_parse_residues("9xyz A:"))
    false
  } catch {
    ScopRafError(msg) => msg == "RAF record not found for 9xyzA"
  }
  assert_true(r)
}

///|
test "scop_raf_get_atoms" {
  let raf_text = raf_line("1sct", "A", "0.01", "150212", [
    ("1", "A", "A"),
    ("2", "G", "G"),
    ("10A", "K", "K"),
  ])
  let index = @src.scop_raf_index(raf_text)
  let sm = @src.scop_raf_get_seqmap(index, @src.scop_parse_residues("1sct A:"))
  let pdb = "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 12.79           N\nATOM      2  CA  GLY A   2      10.2     6.9    -6.0    1.00 13.00           C\nATOM      3  N   LYS A  10A     12.0     7.0    -5.0    1.00 14.00           N\nHETATM    4  O   HOH A  99      30.0    30.0    30.0    1.00 20.00           O\nEND\n"
  let out = @src.scop_raf_get_atoms(sm, pdb)
  let lines = @src.split_by_char(out, '\n'.to_int())
  // the water residue is not part of the sequence map and is skipped;
  // the trailing newline yields one empty final element
  assert_eq(lines.length(), 4)
  assert_eq(lines[3], "")
  assert_eq(
    lines[0],
    "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 12.79           N",
  )
  assert_eq(
    lines[1],
    "ATOM      2  CA  GLY A   2      10.2     6.9    -6.0    1.00 13.00           C",
  )
  assert_eq(
    lines[2],
    "ATOM      3  N   LYS A  10A     12.0     7.0    -5.0    1.00 14.00           N",
  )
}

///|
test "scop_raf_get_atoms_missing_residue" {
  let raf_text = raf_line("1sct", "A", "0.01", "150212", [
    ("1", "A", "A"),
    ("2", "G", "G"),
  ])
  let index = @src.scop_raf_index(raf_text)
  let sm = @src.scop_raf_get_seqmap(index, @src.scop_parse_residues("1sct A:"))
  // the GLY 2 ATOM record is missing
  let pdb = "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 12.79           N\nEND\n"
  let r = try {
    let _ = @src.scop_raf_get_atoms(sm, pdb)
    false
  } catch {
    ScopRafError(msg) =>
      msg.contains("Could not find at least one ATOM or HETATM")
  }
  assert_true(r)
}

///|
test "scop_raf_parse_multi_line" {
  let text = raf_line("1sct", "A", "0.01", "150212", [("1", "A", "A")]) +
    "\n\n" +
    raf_line("1flp", "_", "0.01", "150212", [("5", "D", "D")]) +
    "\n"
  let maps = @src.scop_raf_parse(text)
  assert_eq(maps.length(), 2)
  assert_eq(maps[0].pdbid, "1sct")
  assert_eq(maps[1].pdbid, "1flp")
}

///|
test "scop_dom_parse_and_roundtrip" {
  let text = "# comment line\n\n\td1sctg_\t1sct\tg:\t1.001.001.001.001.001\n\td1flp__\t1flp\t-\t1.001.001.001.001.002\n"
  // NOTE: the first data line intentionally starts with an empty (sid)
  // column in the raw string above; rebuild without that leading tab:
  let text2 = "# comment line\nd1sctg_\t1sct\tg:\t1.001.001.001.001.001\nd1flp__\t1flp\t-\t1.001.001.001.001.002\n"
  let _ = text
  let records = @src.scop_dom_parse(text2)
  assert_eq(records.length(), 2)
  let rec = records[0]
  assert_eq(rec.sid, "d1sctg_")
  assert_eq(rec.residues.pdbid(), "1sct")
  assert_eq(rec.residues.fragments().length(), 1)
  assert_eq(rec.residues.fragments()[0].chain, "g")
  assert_eq(rec.hierarchy, "1.001.001.001.001.001")
  // round-trip: to_string reproduces the input line (tabs, trailing newline)
  assert_eq(rec.to_string(), "d1sctg_\t1sct\tg:\t1.001.001.001.001.001\n")
  // "-" residues keep their dash
  assert_eq(records[1].to_string(), "d1flp__\t1flp\t-\t1.001.001.001.001.002\n")
}

///|
test "scop_dom_invalid_format" {
  let r = try {
    let _ = @src.scop_dom_record_parse("only-three\tcolumns\there")
    false
  } catch {
    ScopDomError(msg) => msg.contains("I don't understand the format of")
  }
  assert_true(r)
}
