///|
/// Bio.SCOP.Raf / Bio.SCOP.Dom Demo.
///
/// Demonstrates the ASTRAL RAF (Rapid Access Format) sequence maps and the
/// (deprecated) SCOP DOMain file: parsing RAF records, indexing, cutting a
/// domain's sequence map, extracting the matching PDB ATOM records, and
/// round-tripping DOM records.

///|
/// Build a synthetic RAF record line for the given residues
/// (resid, ATOM-letter, SEQRES-letter); header is padded to 38 columns.
fn raf_line(
  pdbid : String,
  chainid : String,
  version : String,
  flags : String,
  residues : Array[(String, String, String)],
) -> String {
  let mut line = pdbid + chainid + " " + version + "    " + "26JUN9" + " " + flags
  while line.length() < 38 {
    line = line + " "
  }
  for residue in residues {
    let (resid, atom, seqres) = residue
    let mut field = resid
    while field.length() < 5 {
      field = " " + field
    }
    line = line + field + atom + seqres
  }
  line
}

///|
fn main raise {
  println("=== 1. Parse a RAF sequence map ===")
  let line = raf_line(
    "1sct",
    "A",
    "0.01",
    "150212",
    [("1", "A", "A"), ("2", "G", "G"), ("10A", "K", "K"), ("11", "R", "R")],
  )
  let seqmap = @src.scop_raf_seqmap_parse(line)
  println("pdbid:       \{seqmap.pdbid}")
  println("datestamp:   \{seqmap.pdb_datestamp}")
  println("version:     \{seqmap.version}")
  println("flags:       \{seqmap.flags}")
  println("residues:    \{seqmap.length()}")
  for i in 0..<seqmap.length() {
    let res = seqmap.get(i)
    println(
      "  res[\{i}] chain=\{res.chainid} resid=\{res.resid} atom=\{res.atom} seqres=\{res.seqres}",
    )
  }

  println("\n=== 2. RAF index and domain sequence map ===")
  let raf_text = line + "\n" + raf_line(
    "1flp",
    "_",
    "0.01",
    "150212",
    [("5", "D", "D"), ("6", "E", "E")],
  )
  let index = @src.scop_raf_index(raf_text)
  println("index keys: \{index.keys().collect().to_string()}")
  // SCOP-style residue definition: residues 2 to 10A of chain A
  let residues = @src.scop_parse_residues("1sct A:2-10A")
  let domain = @src.scop_raf_get_seqmap(index, residues)
  println("domain length: \{domain.length()}")
  for i in 0..<domain.length() {
    println("  resid=\{domain.get(i).resid}")
  }

  println("\n=== 3. Extract the matching PDB ATOM records ===")
  let pdb = "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 12.79           N\nATOM      2  CA  GLY A   2      10.200   6.900  -6.000  1.00 13.00           C\nATOM      3  N   LYS A  10A     12.000   7.000  -5.000  1.00 14.00           N\nHETATM    4  O   HOH A  99      30.000  30.000  30.000  1.00 20.00           O\nEND\n"
  let atoms = @src.scop_raf_get_atoms(domain, pdb)
  println(atoms)

  println("\n=== 4. SCOP DOM records ===")
  let dom_text = "# comment lines start with '#'\nd1sctg_\t1sct\tg:\t1.001.001.001.001.001\nd1flp__\t1flp\t-\t1.001.001.001.001.002\n"
  let records = @src.scop_dom_parse(dom_text)
  for rec in records {
    println("sid=\{rec.sid}  hierarchy=\{rec.hierarchy}  fragments=\{rec.residues.fragments().length()}")
  }
  println("--- round trip ---")
  for rec in records {
    println(rec.to_string())
  }

  println("\nDone.")
}
