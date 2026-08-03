///|
/// AlignACE module demo - AlignACE motif discovery output format parser and analysis.
///
/// Run: moon run examples/alignace_demo --target wasm-gc
fn main {
  println("=== AlignACE Module Demo ===")
  println("AlignACE motif discovery output format")
  println("Reference: BioPython Bio.Motifs.AlignAce")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample AlignACE record
  // ------------------------------------------------------------------
  println("=== 1. Load sample AlignACE record ===")
  let record = @bio.alignace_sample()
  println("Loaded sample AlignACE record")
  println("")

  // ------------------------------------------------------------------
  // 2. Print record info
  // ------------------------------------------------------------------
  println("=== 2. Record info ===")
  println("Version: \{record.version}")
  println("Command: \{record.command}")
  println("Number of motifs: \{@bio.alignace_num_motifs(record)}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show motif details
  // ------------------------------------------------------------------
  println("=== 3. Motif details ===")
  for i in 0..<@bio.alignace_num_motifs(record) {
    let m_opt = @bio.alignace_get_motif(record, i)
    match m_opt {
      Some(m) => {
        println("Motif \{i + 1}:")
        println("  width: \{m.width}")
        println("  num_sites: \{m.num_sites}")
        println("  consensus: \{@bio.alignace_motif_consensus(m)}")
        println("  total IC: \{@bio.alignace_total_ic(m)} bits")
        // Show count matrix
        let bases = ["A", "C", "G", "T"]
        for row in 0..<4 {
          let sb = StringBuilder::new()
          sb.write_string("  " + bases[row] + ":")
          for col in 0..<m.width {
            sb.write_string(" " + m.count_matrix[row][col].to_string())
          }
          println(sb.to_string())
        }
        // Show sites
        println("  sites: \{m.sites.length()}")
        for j in 0..<m.sites.length() {
          let s = m.sites[j]
          println("    site \{j}: \{s.sequence_id} @ \{s.position} [\{s.strand}] \{s.sequence}")
        }
        println("")
      }
      None => println("Motif \{i + 1}: not found")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Compute PWM (position weight matrix)
  // ------------------------------------------------------------------
  println("=== 4. Position Weight Matrix (PWM) ===")
  let m0_opt = @bio.alignace_get_motif(record, 0)
  match m0_opt {
    Some(m0) => {
      let pwm = @bio.alignace_to_pwm(m0)
      let bases = ["A", "C", "G", "T"]
      println("PWM for motif 1 (frequencies):")
      for row in 0..<4 {
        let sb = StringBuilder::new()
        sb.write_string("  " + bases[row] + ":")
        for col in 0..<m0.width {
          sb.write_string(" " + pwm[row][col].to_string())
        }
        println(sb.to_string())
      }
      println("")
    }
    None => println("Motif 1 not found")
  }

  // ------------------------------------------------------------------
  // 5. Information content per position
  // ------------------------------------------------------------------
  println("=== 5. Information content per position ===")
  match @bio.alignace_get_motif(record, 0) {
    Some(m0) => {
      let ic = @bio.alignace_information_content(m0)
      let sb = StringBuilder::new()
      sb.write_string("  IC (bits):")
      for col in 0..<ic.length() {
        sb.write_string(" " + ic[col].to_string())
      }
      println(sb.to_string())
      let total = @bio.alignace_total_ic(m0)
      println("  Total IC: \{total} bits")
      println("")
    }
    None => println("Motif 1 not found")
  }

  // ------------------------------------------------------------------
  // 6. Serialize to text and back
  // ------------------------------------------------------------------
  println("=== 6. Serialize to text and back ===")
  let text = @bio.alignace_write(record)
  println("Serialized to \{text.length()} chars")
  let preview_len = if text.length() < 80 { text.length() } else { 80 }
  println("  first \{preview_len} chars: \{text[0:preview_len]}")
  let parsed = @bio.alignace_parse(text)
  println("Deserialized from text:")
  println("  version: \{parsed.version}")
  println("  command: \{parsed.command}")
  println("  motifs: \{@bio.alignace_num_motifs(parsed)}")
  let orig_count = @bio.alignace_num_motifs(record)
  let parsed_count = @bio.alignace_num_motifs(parsed)
  println("Round-trip verify (motif count): \{orig_count == parsed_count}")
  println("")

  // ------------------------------------------------------------------
  // 7. Build a custom AlignACE record
  // ------------------------------------------------------------------
  println("=== 7. Build a custom AlignACE record ===")
  let custom = @bio.AlignAceRecord::new()
  custom.set_version("4.0")
  custom.set_command("alignace -i=custom.fa -numcols=6")
  let m = @bio.AlignAceMotif::new(
    [
      [8, 0, 0, 8, 0, 0],
      [0, 0, 8, 0, 0, 0],
      [0, 8, 0, 0, 0, 8],
      [0, 0, 0, 0, 8, 0],
    ],
  )
  m.add_site(
    @bio.AlignAceSite::new(
      sequence_id="custom1",
      position=5,
      strand="+",
      sequence="ACGAGT",
    ),
  )
  custom.add_motif(m)
  println("Custom record built:")
  println("  version: \{custom.version}")
  println("  command: \{custom.command}")
  println("  motifs: \{@bio.alignace_num_motifs(custom)}")
  match @bio.alignace_get_motif(custom, 0) {
    Some(cm) => {
      println("  motif 1 consensus: \{@bio.alignace_motif_consensus(cm)}")
      println("  motif 1 width: \{cm.width}")
      println("  motif 1 total IC: \{@bio.alignace_total_ic(cm)} bits")
    }
    None => println("  motif 1 not found")
  }
  let custom_text = @bio.alignace_write(custom)
  let custom_parsed = @bio.alignace_parse(custom_text)
  println("Round-trip verify (motif count): \{@bio.alignace_num_motifs(custom) == @bio.alignace_num_motifs(custom_parsed)}")
  println("")

  // ------------------------------------------------------------------
  // 8. Build a motif from scratch and compute IC
  // ------------------------------------------------------------------
  println("=== 8. Information content of a pure motif ===")
  // A pure motif (only one base per column) has IC = 2 bits per position
  let pure = @bio.AlignAceMotif::new(
    [
      [10, 0],
      [0, 0],
      [0, 10],
      [0, 0],
    ],
  )
  let pure_ic = @bio.alignace_information_content(pure)
  println("Pure motif (A then G):")
  println("  consensus: \{@bio.alignace_motif_consensus(pure)}")
  println("  IC[0]: \{pure_ic[0]} bits (should be ~2.0)")
  println("  IC[1]: \{pure_ic[1]} bits (should be ~2.0)")
  println("  Total IC: \{@bio.alignace_total_ic(pure)} bits")
  println("")

  // ------------------------------------------------------------------
  // 9. Summary
  // ------------------------------------------------------------------
  println("=== 9. Summary ===")
  println(@bio.alignace_summary(record))
  println("=== Demo complete ===")
}
