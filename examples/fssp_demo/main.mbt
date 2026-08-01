///|
/// FSSP module demo - Fold classification based on Structure-Structure
/// alignment of Proteins database parser.
///
/// Run: moon run examples/fssp_demo
fn main {
  println("=== FSSP Module Demo ===")
  println("Structure-Structure alignment of Proteins database parser")
  println("Reference: Holm L, Sander C (1996) Science 273:595-603")
  println("")

  // ------------------------------------------------------------------
  // 1. Generate and parse sample FSSP file
  // ------------------------------------------------------------------
  println("=== 1. Generate and parse sample FSSP file ===")
  let text = @bio.fssp_sample_text()
  let data = @bio.fssp_parse(text)
  println("Parsed FSSP file with " + data.n_alignments().to_string() + " alignments")
  println("")

  // ------------------------------------------------------------------
  // 2. Display header information
  // ------------------------------------------------------------------
  println("=== 2. Header information ===")
  let h = data.header()
  println("Reference PDB ID: " + h.pdbid())
  println("Date: " + h.pdbid_date())
  println("Title: " + h.title())
  println("Compound: " + h.compound())
  println("Source: " + h.source())
  println("Author: " + h.author())
  println("Sequence length: " + h.seq_length().to_string())
  println("Number of alignments: " + h.n_alignments().to_string())
  println("Z-score threshold: " + h.threshold())
  println("")

  // ------------------------------------------------------------------
  // 3. Display reference alignment
  // ------------------------------------------------------------------
  println("=== 3. Reference alignment ===")
  let ref_aln = data.reference().unwrap()
  println("PDB ID: " + ref_aln.pdbid())
  println("Start position: " + ref_aln.start_position().to_string())
  println("Aligned sequence: " + ref_aln.aligned_sequence())
  println("Z-score: " + ref_aln.z_score().to_string())
  println("RMSD: " + ref_aln.rmsd().to_string())
  println("Alignment length: " + ref_aln.lali().to_string())
  println("Percent identity: " + ref_aln.pid().to_string() + "%")
  println("")

  // ------------------------------------------------------------------
  // 4. Display all alignments
  // ------------------------------------------------------------------
  println("=== 4. All alignments ===")
  let alignments = data.alignments()
  for i in 0..<alignments.length() {
    let a = alignments[i]
    let marker = if i == 0 { " (reference)" } else { "" }
    println("  " + (i + 1).to_string() + ". " + a.pdbid() + marker)
    println("     Z-score=" + a.z_score().to_string() + " RMSD=" + a.rmsd().to_string() + " PID=" + a.pid().to_string() + "%")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Filter by Z-score
  // ------------------------------------------------------------------
  println("=== 5. Filter alignments by Z-score >= 10.0 ===")
  let high_z = data.filter_by_zscore(10.0)
  println("Found " + high_z.length().to_string() + " alignment(s) with Z-score >= 10.0")
  for a in high_z {
    println("  " + a.pdbid() + ": Z=" + a.z_score().to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Full summary
  // ------------------------------------------------------------------
  println("=== 6. Full FSSP summary ===")
  println(data.to_string())

  println("=== Demo complete ===")
}
