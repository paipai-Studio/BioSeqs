///|
/// Demo for orfik module: ORFik-style ORF detection.

fn main {
  println("=== orfik Demo: ORFik-style ORF Detection ===")
  println("")

  // Example 1: Find all ORFs in six frames
  println("--- Step 1: Find All ORFs ---")
  let seq = "ATGAAATAAATGCCCTAGATGGGGGTGA"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  println("Sequence: " + seq)
  println("ORFs found: " + orfs.length().to_string())
  for orf in orfs {
    println("  start=" + orf.start.to_string() + " stop=" + orf.stop.to_string() +
            " frame=" + @src.orff_frame_label(orf.frame) +
            " protein=" + orf.protein +
            " kozak=" + orf.kozak_score.to_string())
  }
  println("")

  // Example 2: Kozak score
  println("--- Step 2: Kozak Score ---")
  let kozak_seq = "GCCGCCAUGG"
  let score = @src.orff_kozak_score(kozak_seq, 6)
  println("Sequence: " + kozak_seq)
  println("Kozak score at ATG (pos 6): " + score.to_string())
  println("")

  // Example 3: Classify ORFs relative to CDS
  println("--- Step 3: Classify ORFs ---")
  let classified = @src.orff_classify_all(orfs, 0, 9)
  for orf in classified {
    println("  start=" + orf.start.to_string() + " stop=" + orf.stop.to_string() +
            " class=" + @src.orff_class_label(orf.orf_class))
  }
  println("")

  // Example 4: Find uORFs
  println("--- Step 4: Find uORFs ---")
  let uorfs = @src.orff_find_uorfs(seq, 9, min_length=9)
  println("uORFs found: " + uorfs.length().to_string())
  for orf in uorfs {
    println("  start=" + orf.start.to_string() + " stop=" + orf.stop.to_string() + " protein=" + orf.protein)
  }
  println("")

  // Example 5: Count by frame and class
  println("--- Step 5: Summary Statistics ---")
  let by_frame = @src.orff_count_by_frame(orfs)
  println("By frame: F0=" + by_frame.get("F0").unwrap_or(0).to_string() +
          " F1=" + by_frame.get("F1").unwrap_or(0).to_string() +
          " F2=" + by_frame.get("F2").unwrap_or(0).to_string())
  let by_class = @src.orff_count_by_class(classified)
  println("By class: CDS=" + by_class.get("CDS").unwrap_or(0).to_string() +
          " uORF=" + by_class.get("uORF").unwrap_or(0).to_string())
  println("")

  // Example 6: Export ORFs
  println("--- Step 6: Export ORFs ---")
  let exported = @src.orff_export_orfs(orfs)
  println(exported)
}
