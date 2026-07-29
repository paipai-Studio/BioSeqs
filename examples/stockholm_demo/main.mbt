///|
/// Stockholm format parsing example.
///
/// Run: moon run examples/stockholm_demo/main.mbt
fn main {
  println("=== Stockholm Format Parser Demo ===")
  println("")

  let stockholm_content = "# STOCKHOLM 1.0\n" +
    "#=GF ID   TestFamily\n" +
    "#=GF AC   PF00001.20\n" +
    "#=GF DE   A test protein family with secondary structure\n" +
    "#=GC SS_cons    ........((((...)))).....\n" +
    "seq1/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    "seq2/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    "seq3/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    "seq4/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    "#=GR seq1    ........((((...)))).....\n" +
    "//\n"

  println("1. Parsing Stockholm alignment...")
  let ali = @bio.stockholm_parse(stockholm_content)
  println("   Version: \{ali.version}")
  println("   Blocks: \{ali.blocks.length()}")
  println("   GF annotations: \{ali.gf_annotations.keys().length()}")

  println("")
  println("2. GF Annotations:")
  let gf_keys = ali.gf_annotations.keys().collect()
  let mut gi = 0
  while gi < gf_keys.length() {
    println("   \{gf_keys[gi]}: \{ali.gf_annotations[gf_keys[gi]]}")
    gi = gi + 1
  }

  println("")
  println("3. Sequences in first block:")
  if ali.blocks.length() > 0 {
    let block = ali.blocks[0]
    let mut si = 0
    while si < block.sequences.length() {
      let seq = block.sequences[si]
      println("   \{seq.name}/\{seq.start}-\{seq.end} (\{seq.aligned_seq.length()} chars)")
      si = si + 1
    }

    println("")
    println("4. Secondary structure: \{block.secondary_structure}")

    println("")
    println("5. Percent identity: \{@bio.stockholm_percent_identity(block)}")

    println("")
    println("6. Per-column conservation (first 10 columns):")
    let cons = @bio.stockholm_conservation(block)
    let mut ci = 0
    while ci < 10 && ci < cons.length() {
      println("   Column \{ci}: \{cons[ci]}")
      ci = ci + 1
    }
  }

  println("")
  println("7. Converting to FASTA...")
  let fasta = @bio.stockholm_to_fasta(ali)
  println(fasta)

  println("8. Writing Stockholm format...")
  let output = @bio.stockholm_write(ali)
  println(output)

  println("9. Merging blocks...")
  let merged = @bio.stockholm_merge_blocks(ali)
  println("   Merged blocks: \{merged.blocks.length()}")
  if merged.blocks.length() > 0 {
    println("   Sequences in merged block: \{merged.blocks[0].sequences.length()}")
  }

  println("")
  println("=== Done ===")
}