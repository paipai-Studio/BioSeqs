///|
/// dnashape module demo - DNA shape prediction from sequence.
///
/// Demonstrates prediction of four DNA shape features (MGW, Roll, ProT, HelT)
/// from nucleotide sequences using dinucleotide lookup tables.
///
/// Reference: Chiu TP et al. (2016) DNAshapeR. Bioinformatics 32:3415-3417.
/// Run: moon run examples/dnashape_demo
fn main {
  println("=== dnashape Module Demo ===")
  println("DNA shape prediction from sequence")
  println("Reference: Bioconductor DNAshapeR")
  println("")

  // ------------------------------------------------------------------
  // 1. Predict all shape features for a sample sequence
  // ------------------------------------------------------------------
  println("=== 1. Predict All Shape Features ===")
  let seq = @bio.dnashape_sample_sequence()
  println("Sequence: " + seq)
  println("Length: " + seq.length().to_string() + " bp")
  println("")

  let pred = @bio.dnashape_predict_all(seq)
  println("MGW  mean: " + pred.mgw.mean().to_string() + " (range " + pred.mgw.min().to_string() + " to " + pred.mgw.max().to_string() + ")")
  println("Roll mean: " + pred.roll.mean().to_string() + " (range " + pred.roll.min().to_string() + " to " + pred.roll.max().to_string() + ")")
  println("ProT mean: " + pred.prot.mean().to_string() + " (range " + pred.prot.min().to_string() + " to " + pred.prot.max().to_string() + ")")
  println("HelT mean: " + pred.helt.mean().to_string() + " (range " + pred.helt.min().to_string() + " to " + pred.helt.max().to_string() + ")")
  println("")

  // ------------------------------------------------------------------
  // 2. Compare A-tract vs GC-rich sequences
  // ------------------------------------------------------------------
  println("=== 2. A-tract vs GC-rich Comparison ===")
  let atract = @bio.dnashape_sample_atract()
  let gc_rich = @bio.dnashape_sample_gc_rich()
  let mgw_feat = @bio.dnashape_make_feature("MGW")

  let mgw_atract = @bio.dnashape_predict(atract, mgw_feat)
  let mgw_gc = @bio.dnashape_predict(gc_rich, mgw_feat)

  println("A-tract MGW mean: " + mgw_atract.mean().to_string() + " (narrow groove)")
  println("GC-rich MGW mean: " + mgw_gc.mean().to_string() + " (wide groove)")
  println("A-tracts have narrower minor groove than GC-rich sequences.")
  println("")

  // ------------------------------------------------------------------
  // 3. Individual feature prediction
  // ------------------------------------------------------------------
  println("=== 3. Individual Feature Prediction ===")
  let test_seq = "ATCGATCGATCGATCGATCG"
  println("Test sequence: " + test_seq)

  let roll_feat = @bio.dnashape_make_feature("Roll")
  let roll_result = @bio.dnashape_predict(test_seq, roll_feat)
  println("Roll values: " + @bio.dnashape_result_to_string(roll_result))
  println("")

  // ------------------------------------------------------------------
  // 4. Dinucleotide composition
  // ------------------------------------------------------------------
  println("=== 4. Dinucleotide Composition ===")
  let comp_seq = "ATCGATCG"
  let counts = @bio.dnashape_dinuc_composition(comp_seq)
  println("Sequence: " + comp_seq)
  let bases = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
  let mut i = 0
  while i < 16 {
    if counts[i] > 0 {
      println("  " + bases[i] + ": " + counts[i].to_string())
    }
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. GC content
  // ------------------------------------------------------------------
  println("=== 5. GC Content ===")
  println("A-tract GC content: " + @bio.dnashape_gc_content(atract).to_string())
  println("GC-rich GC content: " + @bio.dnashape_gc_content(gc_rich).to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. Heatmap visualization
  // ------------------------------------------------------------------
  println("=== 6. Heatmap Visualization ===")
  let heatmap_seq = "ATCGATCGATCGATCGATCG"
  let helt_feat = @bio.dnashape_make_feature("HelT")
  let helt_result = @bio.dnashape_predict(heatmap_seq, helt_feat)
  println("HelT heatmap for '" + heatmap_seq + "':")
  println("  " + @bio.dnashape_heatmap(helt_result))
  println("  (low=space, high=@)")
  println("")

  // ------------------------------------------------------------------
  // 7. Reverse complement
  // ------------------------------------------------------------------
  println("=== 7. Reverse Complement ===")
  let sample_seq = "ATCGATCG"
  let rc = @bio.dnashape_revcomp(sample_seq)
  println("Sequence: " + sample_seq)
  println("RevComp:  " + rc)
  println("")

  println("=== Demo complete ===")
}
