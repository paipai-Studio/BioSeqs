// Bio.codonalign-inspired advanced codon alignment and selection pressure
// analysis workflow.
//
// Exercises the full selection-pressure pipeline on synthetic codon
// alignments:
//   1. Z-test for selection (positive / purifying / neutrality).
//   2. Fisher's exact test for neutrality (2×2 contingency table).
//   3. Codon alignment builder from protein alignment + coding sequences.
//   4. Sliding-window dN/dS scan for selection hotspots.
//   5. Benjamini–Hochberg FDR correction of multiple p-values.
//   6. Pairwise Ka/Ks table with Z-test and FDR across a small panel.

///|
fn caa_demo_round(value : Double, digits : Int) -> Double {
  let scale = @math.pow(10.0, digits.to_double())
  (value * scale).round() / scale
}

///|
fn caa_demo_print_d(label : String, value : Double) -> Unit {
  println("  " + label + ": " + caa_demo_round(value, 6).to_string())
}

///|
fn caa_demo_print_s(label : String, value : String) -> Unit {
  println("  " + label + ": " + value)
}

///|
fn main {
  println("=== Bio.codonalign Advanced Demo ===")

  // -----------------------------------------------------------------------
  // 1. Z-test for selection
  // -----------------------------------------------------------------------
  println("\n1. Z-test for selection (Nei–Gojobori approximate variance)")
  let pos_seqs = @src.create_positive_selection_alignment()
  let pur_seqs = @src.create_purifying_selection_alignment()
  let s1_pos = pos_seqs[0]
  let s2_pos = pos_seqs[1]
  let s1_pur = pur_seqs[0]
  let s2_pur = pur_seqs[1]
  let pos_result = @src.codon_test_selection(s1_pos, s2_pos, test_type="positive") catch {
    CodonAlignAdvancedError(msg) => abort("Z-test positive failed: " + msg)
  }
  println("  Positive-selection pair:")
  caa_demo_print_d("dN", pos_result.dn())
  caa_demo_print_d("dS", pos_result.ds())
  caa_demo_print_d("dN/dS", pos_result.dnds())
  caa_demo_print_d("Z-score", pos_result.z_score())
  caa_demo_print_d("p-value", pos_result.p_value())
  caa_demo_print_s("conclusion", pos_result.conclusion())
  let pur_result = @src.codon_test_selection(s1_pur, s2_pur, test_type="purifying") catch {
    CodonAlignAdvancedError(msg) => abort("Z-test purifying failed: " + msg)
  }
  println("  Purifying-selection pair:")
  caa_demo_print_d("dN", pur_result.dn())
  caa_demo_print_d("dS", pur_result.ds())
  caa_demo_print_d("dN/dS", pur_result.dnds())
  caa_demo_print_d("Z-score", pur_result.z_score())
  caa_demo_print_d("p-value", pur_result.p_value())
  caa_demo_print_s("conclusion", pur_result.conclusion())

  // -----------------------------------------------------------------------
  // 2. Fisher's exact test for neutrality
  // -----------------------------------------------------------------------
  println("\n2. Fisher's exact test for neutrality")
  let fisher_pos = @src.codon_test_neutrality(s1_pos, s2_pos, test_type="greater") catch {
    CodonAlignAdvancedError(msg) => abort("Fisher test failed: " + msg)
  }
  println("  Positive-selection pair (greater):")
  caa_demo_print_d("p-value", fisher_pos.p_value())
  caa_demo_print_d("odds ratio", fisher_pos.odds_ratio())
  println(
    "  Nd=" +
    fisher_pos.n_diff().to_string() +
    ", Sd=" +
    fisher_pos.s_diff().to_string(),
  )
  caa_demo_print_s("conclusion", fisher_pos.conclusion())

  // -----------------------------------------------------------------------
  // 3. Codon alignment builder
  // -----------------------------------------------------------------------
  println("\n3. Codon alignment builder (protein alignment + CDS)")
  let (protein_aln, coding_seqs) = @src.create_demo_protein_alignment()
  let codon_aln = @src.build_codon_alignment(
    protein_aln,
    coding_seqs,
    names=["geneA", "geneB", "geneC"],
  ) catch {
    CodonAlignAdvancedError(msg) => abort("Codon alignment builder failed: " + msg)
  }
  println("  n_codons (columns): " + codon_aln.n_codons().to_string())
  let aln_seqs = codon_aln.sequences()
  let aln_names = codon_aln.names()
  for i in 0..<aln_seqs.length() {
    println("  " + aln_names[i] + ": " + aln_seqs[i])
  }

  // -----------------------------------------------------------------------
  // 4. Sliding-window dN/dS scan
  // -----------------------------------------------------------------------
  println("\n4. Sliding-window dN/dS scan")
  // 21-codon (63 nt) sequences.
  let sw_s1 = "ATGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGG"
  let sw_s2 = "ATGGGTCTGGAGAGACACCCCAGCATGGGGTCTGGAGAGACACCCCAGCATGGGGTCTGGAGA"
  let sw_result = @src.sliding_window_dnds(
    sw_s1,
    sw_s2,
    window_size=9,
    step_size=6,
  ) catch {
    CodonAlignAdvancedError(msg) => abort("Sliding window failed: " + msg)
  }
  let windows = sw_result.windows()
  println("  " + windows.length().to_string() + " windows (size=9, step=6):")
  for w in windows {
    println(
      "  codons " +
      w.start_codon().to_string() +
      "-" +
      w.end_codon().to_string() +
      " (n=" +
      w.n_codons().to_string() +
      "): dN=" +
      caa_demo_round(w.dn(), 6).to_string() +
      ", dS=" +
      caa_demo_round(w.ds(), 6).to_string() +
      ", dN/dS=" +
      caa_demo_round(w.dnds(), 4).to_string(),
    )
  }

  // -----------------------------------------------------------------------
  // 5. Benjamini–Hochberg FDR correction
  // -----------------------------------------------------------------------
  println("\n5. Benjamini–Hochberg FDR correction")
  let raw_pvalues = [0.001, 0.04, 0.03, 0.5, 0.01, 0.2]
  let adjusted = @src.codon_align_advanced_bh_fdr(raw_pvalues)
  println("  raw → adjusted:")
  for i in 0..<raw_pvalues.length() {
    println(
      "  " +
      caa_demo_round(raw_pvalues[i], 4).to_string() +
      " → " +
      caa_demo_round(adjusted[i], 4).to_string(),
    )
  }

  // -----------------------------------------------------------------------
  // 6. Pairwise Ka/Ks table with Z-test and FDR
  // -----------------------------------------------------------------------
  println("\n6. Pairwise Ka/Ks table with Z-test and FDR")
  let full_panel = @src.create_positive_selection_alignment()
  // Use a 4-sequence subset (6 pairwise comparisons) for a concise table.
  let panel : Array[String] = [
    full_panel[0], full_panel[1], full_panel[2], full_panel[3],
  ]
  let kaks_table = @src.pairwise_kaks_table(
    panel,
    names=["pos01", "pos02", "pos03", "pos04"],
    fdr_threshold=0.05,
  ) catch {
    CodonAlignAdvancedError(msg) => abort("Pairwise Ka/Ks table failed: " + msg)
  }
  println("  " + kaks_table.length().to_string() + " pairwise comparisons:")
  for row in kaks_table {
    println(
      "  " +
      row.seq1_name() +
      " vs " +
      row.seq2_name() +
      ": dN/dS=" +
      caa_demo_round(row.dnds(), 4).to_string() +
      ", p=" +
      caa_demo_round(row.p_value(), 4).to_string() +
      ", FDR=" +
      caa_demo_round(row.fdr(), 4).to_string() +
      " — " +
      row.conclusion(),
    )
  }

  println("\n=== Demo complete ===")
}
