///|
/// Wise module demo - Demonstrates parsing and analyzing GeneWise/GenomeWise output.
///
/// Run: moon run examples/wise_demo --target wasm-gc
fn main {
  println("=== Wise Module Demo ===")
  println("Parser for Wise2 (GeneWise/GenomeWise) output")
  println("Reference: Birney E et al. (2004) Genome Research 14:988-995")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample GeneWise output
  // ------------------------------------------------------------------
  println("=== 1. Load sample GeneWise output ===")
  let text = @bio.wise_sample_output()
  let result = @bio.wise_parse(text)
  println("Parsed GeneWise output successfully")
  println("")

  // ------------------------------------------------------------------
  // 2. Print result info
  // ------------------------------------------------------------------
  println("=== 2. Result information ===")
  println("Protein ID: \{result.protein_id}")
  println("Gene ID: \{result.gene_id}")
  println("Score: \{result.score}")
  println("Bits score: \{result.bits_score}")
  println("Number of exons: \{result.get_num_exons()}")
  println("Number of introns: \{result.get_num_introns()}")
  println("Total exon length: \{result.get_total_exon_length()} nt")
  println("")

  // ------------------------------------------------------------------
  // 3. Show exon details
  // ------------------------------------------------------------------
  println("=== 3. Exon details ===")
  let exons = result.get_exons()
  for i in 0..<exons.length() {
    let e = exons[i]
    println("Exon \{i + 1}:")
    println("  start: \{e.start}")
    println("  end: \{e.end}")
    println("  length: \{e.length()} nt")
    println("  phase: \{e.phase}")
    println("  end_phase: \{e.end_phase}")
    println("  score: \{e.score}")
    println("  protein_start: \{e.protein_start}")
    println("  protein_end: \{e.protein_end}")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Show intron details
  // ------------------------------------------------------------------
  println("=== 4. Intron details ===")
  let introns = result.get_introns()
  for i in 0..<introns.length() {
    let intron = introns[i]
    println("Intron \{i + 1}:")
    println("  start: \{intron.start}")
    println("  end: \{intron.end}")
    println("  length: \{intron.length} nt")
    println("  donor score (5' splice): \{intron.donor_score}")
    println("  acceptor score (3' splice): \{intron.acceptor_score}")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Demonstrate gene structure
  // ------------------------------------------------------------------
  println("=== 5. Gene structure ===")
  let structure = @bio.wise_gene_structure(result)
  println(structure)
  println("")

  // ------------------------------------------------------------------
  // 6. Demonstrate translation
  // ------------------------------------------------------------------
  println("=== 6. Translation ===")
  // Build a synthetic gene sequence matching the exon coordinates.
  // Exon 1: 0-99 (100 nt), Exon 2: 200-299 (100 nt), Exon 3: 400-479 (80 nt)
  // Total gene length needed: 480 nt.
  // Fill exon regions with codons that translate to known amino acids,
  // and intron regions with arbitrary sequence.
  let gene_seq = wise_build_gene_sequence()
  println("Gene sequence length: \{gene_seq.length()} nt")
  let translated = @bio.wise_translate_gene(result, gene_seq)
  println("Translated protein (first 20 aa): \{translated.substring(start=0, end=20).to_string()}")
  println("Total translated length: \{translated.length()} aa")
  // Build a matching protein to compute identity.
  // Exon 1 (100 nt) -> 33 M's, Exon 2 (100 nt) -> 33 A's, Exon 3 (80 nt) -> 26 G's
  let expected_protein = wise_build_expected_protein()
  let pct = @bio.wise_percent_identity(result, expected_protein, gene_seq)
  println("Percent identity vs expected protein: \{pct}%")
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize back to text
  // ------------------------------------------------------------------
  println("=== 7. Serialize back to text ===")
  let rewritten = @bio.wise_write(result)
  println(rewritten)
  println("")

  // ------------------------------------------------------------------
  // 8. Print summary
  // ------------------------------------------------------------------
  println("=== 8. Summary ===")
  let summary = @bio.wise_summary(result)
  println(summary)
  println("")

  // ------------------------------------------------------------------
  // 9. Demonstrate GenomeWise result
  // ------------------------------------------------------------------
  println("=== 9. GenomeWise result (multi-segment) ===")
  let gw = @bio.GenomeWiseResult::new()
  let seg1 = @bio.GenomeWiseSegment::new(
    segment_id="segment_1",
    start=0,
    end=500,
    score=45.5,
  )
  seg1.wise_result.set_protein_id("P12345")
  seg1.wise_result.set_gene_id("G67890")
  seg1.wise_result.add_exon(@bio.WiseExon::new(start=0, end=99, phase=0, score=45.5))
  let seg2 = @bio.GenomeWiseSegment::new(
    segment_id="segment_2",
    start=1000,
    end=1500,
    score=52.3,
  )
  seg2.wise_result.set_protein_id("P12345")
  seg2.wise_result.set_gene_id("G67890")
  seg2.wise_result.add_exon(@bio.WiseExon::new(start=0, end=99, phase=0, score=52.3))
  gw.add_segment(seg1)
  gw.add_segment(seg2)
  println("GenomeWise segments: \{gw.get_num_segments()}")
  let segs = gw.get_segments()
  for i in 0..<segs.length() {
    let s = segs[i]
    println("Segment \{i + 1}: id=\{s.segment_id}, range=\{s.start}-\{s.end}, score=\{s.score}, exons=\{s.wise_result.get_num_exons()}")
  }
  println("")
  println("=== Wise Module Demo Complete ===")
}

///|
/// Build a synthetic gene sequence matching the sample exon coordinates.
/// Exon 1: 0-99 (100 nt) of "ATG" -> translates to Met (M)
/// Intron 1: 100-199 (100 nt) of arbitrary sequence
/// Exon 2: 200-299 (100 nt) of "GCC" -> translates to Ala (A)
/// Intron 2: 300-399 (100 nt) of arbitrary sequence
/// Exon 3: 400-479 (80 nt) of "GGT" -> translates to Gly (G)
/// Total: 480 nt
fn wise_build_gene_sequence() -> String {
  let exon1 = wise_repeat_pattern("ATG", 100)
  let intron1 = wise_repeat_pattern("GTAAC", 100)
  let exon2 = wise_repeat_pattern("GCC", 100)
  let intron2 = wise_repeat_pattern("GTAAC", 100)
  let exon3 = wise_repeat_pattern("GGT", 80)
  exon1 + intron1 + exon2 + intron2 + exon3
}

///|
/// Build the expected protein from the sample exons.
/// Exon 1 (100 nt) -> 33 M's, Exon 2 (100 nt) -> 33 A's, Exon 3 (80 nt) -> 26 G's
fn wise_build_expected_protein() -> String {
  wise_repeat_char('M', 33) + wise_repeat_char('A', 33) + wise_repeat_char('G', 26)
}

///|
/// Repeat a pattern string to fill an exact length.
fn wise_repeat_pattern(pattern : String, length : Int) -> String {
  let sb = StringBuilder::new()
  let chars = pattern.to_array()
  let plen = chars.length()
  for i in 0..<length {
    sb.write_char(chars[i % plen])
  }
  sb.to_string()
}

///|
/// Repeat a character to fill an exact length.
fn wise_repeat_char(c : Char, length : Int) -> String {
  let sb = StringBuilder::new()
  for _ in 0..<length {
    sb.write_char(c)
  }
  sb.to_string()
}
