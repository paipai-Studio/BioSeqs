///|
/// FragmentMapper demo example.
///
/// Demonstrates protein structure fragment mapping and classification
/// using the Bio.PDB.FragmentMapper module.
///
/// Run: moon run examples/fragment_mapper_demo/main.mbt

fn main {
  // === 1. Classify DSSP secondary structure ===
  println("=== 1. DSSP Secondary Structure Classification ===")
  let ss_string = "HHHHHHLLLLLLEEEEELLLTTSSCC"
  let classified = @bio.fm_classify_secondary_structure(ss_string)
  println("  Input:  " + ss_string)
  let mut class_str = ""
  for ft in classified {
    class_str = class_str + ft.to_char()
  }
  println("  Output: " + class_str)

  // Count types
  let mut helix_count = 0
  let mut sheet_count = 0
  let mut loop_count = 0
  let mut turn_count = 0
  for ft in classified {
    if ft == @bio.ft_helix() {
      helix_count = helix_count + 1
    } else if ft == @bio.ft_sheet() {
      sheet_count = sheet_count + 1
    } else if ft == @bio.ft_loop() {
      loop_count = loop_count + 1
    } else if ft == @bio.ft_turn() {
      turn_count = turn_count + 1
    }
  }
  println("  Helix: " + helix_count.to_string() + ", Sheet: " + sheet_count.to_string() + ", Loop: " + loop_count.to_string() + ", Turn: " + turn_count.to_string())

  // === 2. Assign fragments from sequence + SS ===
  println("\n=== 2. Fragment Assignment ===")
  let sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLR"
  let ss = "HHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHCCCEEEEEEEECCCHHHHHHHHHHHHHHHHHHCCCCCCEEEEEECCCCCHHHHHH"
  let result = @bio.fm_assign_fragments(sequence, ss)
  println("  Sequence length: " + result.sequence_length.to_string())
  println("  Number of fragments: " + result.num_fragments().to_string())
  println("  Fragments:")
  for frag in result.fragments {
    let (start, end) = @bio.fm_fragment_boundaries(frag)
    let frag_seq = @bio.fm_fragment_sequence(frag)
    let frag_type = @bio.fm_fragment_type(frag)
    println(
      "    " + frag.id + ": type=" + frag_type.to_string() + ", residues " + start.to_string() + "-" + end.to_string() + ", length=" + frag.length().to_string() + ", seq=" + frag_seq,
    )
  }

  // === 3. Merge adjacent fragments ===
  println("\n=== 3. Merge Adjacent Fragments ===")
  let merged = @bio.fm_merge_fragments(result.fragments)
  println("  Before merge: " + result.num_fragments().to_string() + " fragments")
  println("  After merge:  " + merged.length().to_string() + " fragments")

  // === 4. Filter by length ===
  println("\n=== 4. Filter by Length (min 5 residues) ===")
  let filtered = @bio.fm_filter_by_length(merged, 5)
  println("  After filtering (>= 5 residues): " + filtered.length().to_string() + " fragments")
  for frag in filtered {
    println(
      "    " + frag.id + ": type=" + frag.fragment_type.to_string() + ", length=" + frag.length().to_string(),
    )
  }

  // === 5. Get fragment by ID ===
  println("\n=== 5. Get Fragment by ID ===")
  match @bio.fm_get_fragment(result, "fragment_0") {
    Some(frag) => {
      println(
        "  Found fragment_0: residues " + frag.start_residue.to_string() + "-" + frag.end_residue.to_string() + ", type=" + frag.fragment_type.to_string(),
      )
    }
    None => println("  fragment_0 not found")
  }

  // === 6. Parse PDB content ===
  println("\n=== 6. Parse PDB Content ===")
  let pdb_content = "ATOM      1  N   ALA A   1      11.104  6.134  10.298  1.00 20.00           N\n" +
    "ATOM      2  CA  ALA A   1      11.639  6.168  10.426  1.00 20.00           C\n" +
    "ATOM      3  C   ALA A   1      12.221  6.201  10.542  1.00 20.00           C\n" +
    "ATOM      4  O   ALA A   1      12.803  6.234  10.658  1.00 20.00           O\n" +
    "ATOM      5  N   GLY A   2      13.385  6.267  10.774  1.00 20.00           N\n" +
    "ATOM      6  CA  GLY A   2      13.968  6.300  10.891  1.00 20.00           C\n" +
    "ATOM      7  C   GLY A   2      14.549  6.333  11.007  1.00 20.00           C\n" +
    "ATOM      8  O   GLY A   2      15.131  6.366  11.123  1.00 20.00           O\n" +
    "ATOM      9  N   SER A   3      15.712  6.399  11.239  1.00 20.00           N\n" +
    "ATOM     10  CA  SER A   3      16.294  6.432  11.355  1.00 20.00           C\n" +
    "END\n"
  let pdb_result = @bio.fm_parse_fragments(pdb_content)
  println("  Sequence: " + pdb_result.protein_sequence)
  println("  Sequence length: " + pdb_result.sequence_length.to_string())
  println("  Fragments: " + pdb_result.num_fragments().to_string())
  for frag in pdb_result.fragments {
    println(
      "    " + frag.id + ": chain=" + frag.chain_id + ", residues " + frag.start_residue.to_string() + "-" + frag.end_residue.to_string() + ", atom_count=" + frag.residues.length().to_string(),
    )
  }

  // === 7. Coverage analysis ===
  println("\n=== 7. Coverage Analysis ===")
  let cov = @bio.fm_coverage(result)
  let mut covered = 0
  let mut total = 0
  for _, val in cov {
    if val > 0.0 {
      covered = covered + 1
    }
    total = total + 1
  }
  println("  Total positions: " + total.to_string())
  println("  Covered positions: " + covered.to_string())
  if total > 0 {
    let pct = (covered.to_double() / total.to_double()) * 100.0
    println("  Coverage: " + pct.to_string() + "%")
  }

  println("\n=== Done ===")
}