///|
/// PDB SeqIO Demo.
///
/// Demonstrates Bio.SeqIO.PdbIO module: extracting amino acid sequences
/// from PDB structures via SEQRES and ATOM records.
/// Inspired by BioPython's Bio.SeqIO.PdbIO module.

fn main {
  println("=== PDB SeqIO Demo ===")
  println("Bio.SeqIO.PdbIO: Sequence extraction from PDB structures")
  println("Inspired by BioPython's Bio.SeqIO.PdbIO module\n")

  // === 1. Construct PDB text with SEQRES and ATOM records ===
  println("=== 1. Construct PDB Text ===")
  // SEQRES has 9 residues; ATOM has 8 (missing PRO at position 5,
  // simulating an unresolved residue common in real structures).
  let pdb_text = "HEADER    PROTEIN" +
    " ".repeat(33) +
    "01-JAN-24   1ABC\n" +
    "SEQRES   1 A   9  ALA GLY VAL LYS PRO SER THR MET PHE\n" +
    "ATOM      1  CA  ALA A   1      11.104  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      2  CA  GLY A   2      12.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      3  CA  VAL A   3      13.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      4  CA  LYS A   4      14.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      5  CA  SER A   6      15.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      6  CA  THR A   7      16.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      7  CA  MET A   8      17.000  6.134  10.298  1.00 20.00           C\n" +
    "ATOM      8  CA  PHE A   9      18.000  6.134  10.298  1.00 20.00           C\n" +
    "TER       9      PHE A   9\n" +
    "END\n"
  println("  Constructed PDB text with:")
  println("    - HEADER (PDB ID: 1ABC)")
  println("    - SEQRES: 9 residues (ALA GLY VAL LYS PRO SER THR MET PHE)")
  println("    - ATOM:   8 residues (PRO at position 5 is unresolved)")

  // === 2. Extract sequence from SEQRES records ===
  println("\n=== 2. Extract Sequence from SEQRES ===")
  let seqres_records = @bio.pdb_seqres_parser(pdb_text)
  let mut ri = 0
  while ri < seqres_records.length() {
    let r = seqres_records[ri]
    println("  Record: id=\{r.id}, name=\{r.name}")
    println("  Sequence: \{r.seq.to_string()}")
    println("  Description: \{r.description}")
    ri = ri + 1
  }

  // === 3. Extract sequence from ATOM records ===
  println("\n=== 3. Extract Sequence from ATOM ===")
  let atom_records = @bio.pdb_atom_parser(pdb_text)
  ri = 0
  while ri < atom_records.length() {
    let r = atom_records[ri]
    println("  Record: id=\{r.id}, name=\{r.name}")
    println("  Sequence: \{r.seq.to_string()}")
    println("  Description: \{r.description}")
    ri = ri + 1
  }

  // === 4. Compare SEQRES vs ATOM sequences ===
  println("\n=== 4. Compare SEQRES vs ATOM Sequences ===")
  let seqres_seq = if seqres_records.length() > 0 {
    seqres_records[0].seq.to_string()
  } else {
    ""
  }
  let atom_seq = if atom_records.length() > 0 {
    atom_records[0].seq.to_string()
  } else {
    ""
  }
  println("  SEQRES sequence (\{seqres_seq.length()} residues): \{seqres_seq}")
  println("  ATOM   sequence (\{atom_seq.length()} residues): \{atom_seq}")
  if seqres_seq == atom_seq {
    println("  Result: sequences are identical")
  } else {
    println("  Result: sequences differ (ATOM missing unresolved residues)")
  }

  // === 5. three_to_one_aa conversion ===
  println("\n=== 5. Three-to-One Amino Acid Conversion ===")
  let three_letter_codes = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "XXX",
  ]
  println("  Three-letter -> One-letter:")
  let mut ci = 0
  while ci < three_letter_codes.length() {
    let three = three_letter_codes[ci]
    let one = @bio.three_to_one_aa(three)
    println("    \{three} -> \{one}")
    ci = ci + 1
  }

  // === 6. extract_sequence_from_pdb_lines ===
  println("\n=== 6. extract_sequence_from_pdb_lines (chain A) ===")
  let split_result = pdb_text.split("\n").to_array()
  let pdb_lines : Array[String] = Array::new()
  let mut li = 0
  while li < split_result.length() {
    pdb_lines.push(split_result[li].to_owned())
    li = li + 1
  }
  let chain_a_seq = @bio.extract_sequence_from_pdb_lines(pdb_lines, "A")
  println("  Chain A sequence from ATOM lines: \{chain_a_seq}")
  println("  Length: \{chain_a_seq.length()} residues")

  // === 7. Write SEQRES records back ===
  println("\n=== 7. Write PDB SEQRES Records ===")
  let rewritten_seqres = @bio.write_pdb_seqrecords(seqres_records)
  println("  Rewritten SEQRES records:")
  println(rewritten_seqres)

  // === 8. pdb_to_seqrecords unified interface ===
  println("=== 8. pdb_to_seqrecords Unified Interface ===")
  let atom_via_unified = @bio.pdb_to_seqrecords(pdb_text, source="atom")
  let seqres_via_unified = @bio.pdb_to_seqrecords(pdb_text, source="seqres")
  println("  source=\"atom\":")
  ri = 0
  while ri < atom_via_unified.length() {
    let r = atom_via_unified[ri]
    println("    id=\{r.id}, seq=\{r.seq.to_string()}")
    ri = ri + 1
  }
  println("  source=\"seqres\":")
  ri = 0
  while ri < seqres_via_unified.length() {
    let r = seqres_via_unified[ri]
    println("    id=\{r.id}, seq=\{r.seq.to_string()}")
    ri = ri + 1
  }

  println("\n=== PDB SeqIO Demo Complete ===")
}
