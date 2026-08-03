///|
/// seqLogo demo - Sequence logo computation and rendering.
///
/// Demonstrates PWM construction, information content, logo computation,
/// consensus extraction, and ASCII/table rendering.
///
/// Reference: Bembom O (2007) seqLogo: Sequence logos for DNA sequence
/// alignments. R package.
/// Run: moon run examples/seqlogo_demo
fn main {
  println("=== seqLogo Demo ===")
  println("Sequence logo computation and rendering")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample PWM
  // ------------------------------------------------------------------
  println("=== 1. Load sample PWM ===")
  let pwm = @bio.seqlogo_sample_pwm()
  println("Loaded sample PWM (4 bases x \{pwm.width()} positions)")
  println("")

  // ------------------------------------------------------------------
  // 2. Print PWM info (width, to_string first few positions)
  // ------------------------------------------------------------------
  println("=== 2. PWM info ===")
  println("Width: \{pwm.width()}")
  println("Alphabet size: \{pwm.alphabet_size()}")
  let alphabet = pwm.alphabet()
  println("Alphabet: \{alphabet[0]}\{alphabet[1]}\{alphabet[2]}\{alphabet[3]}")
  println("First 3 positions (A, C, G, T):")
  let pos_limit = if pwm.width() < 3 { pwm.width() } else { 3 }
  for pos in 0..<pos_limit {
    println(
      "  pos \{pos}: A=\{pwm.get(0, pos)}, C=\{pwm.get(1, pos)}, G=\{pwm.get(2, pos)}, T=\{pwm.get(3, pos)}",
    )
  }
  println("")
  println("Full PWM:")
  println(pwm.to_string())

  // ------------------------------------------------------------------
  // 3. Compute information content per position, show IC values
  // ------------------------------------------------------------------
  println("=== 3. Per-position information content ===")
  let ic = @bio.seqlogo_information_content(pwm)
  for pos in 0..<ic.length() {
    println("  pos \{pos}: IC=\{ic[pos]} bits")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Compute total IC, show max possible IC
  // ------------------------------------------------------------------
  println("=== 4. Total IC and max possible IC ===")
  let total_ic = @bio.seqlogo_total_information_content(pwm)
  let max_ic_per_pos = @bio.seqlogo_max_information_content(pwm.alphabet_size())
  let max_total = max_ic_per_pos * pwm.width().to_double()
  println("Total IC: \{total_ic} bits")
  println("Max IC per position: \{max_ic_per_pos} bits")
  println("Max possible total IC: \{max_total} bits")
  println("")

  // ------------------------------------------------------------------
  // 5. Compute logo, show letter heights per column
  // ------------------------------------------------------------------
  println("=== 5. Logo letter heights per column ===")
  let logo = @bio.seqlogo_compute_logo(pwm)
  let columns = logo.columns()
  for col in columns {
    let buf = StringBuilder::new()
    buf.write_string("  pos ")
    buf.write_string(col.position().to_string())
    buf.write_string(" (IC=")
    buf.write_string(col.info_content().to_string())
    buf.write_string("): ")
    let letters = col.letters()
    for i in 0..<letters.length() {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_char(letters[i].letter())
      buf.write_string("=")
      buf.write_string(letters[i].height().to_string())
    }
    println(buf.to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Show consensus sequence
  // ------------------------------------------------------------------
  println("=== 6. Consensus sequence ===")
  let consensus = @bio.seqlogo_consensus_sequence(pwm)
  println("Consensus: \{consensus}")
  println("")

  // ------------------------------------------------------------------
  // 7. Render ASCII art logo
  // ------------------------------------------------------------------
  println("=== 7. ASCII art logo ===")
  let ascii = @bio.seqlogo_to_ascii(logo, 10)
  println(ascii)

  // ------------------------------------------------------------------
  // 8. Show text table (position, base, height, IC)
  // ------------------------------------------------------------------
  println("=== 8. Text table ===")
  let table = @bio.seqlogo_to_text_table(logo)
  println(table)

  // ------------------------------------------------------------------
  // 9. Build PWM from sample sequences, compare with sample PWM
  // ------------------------------------------------------------------
  println("=== 9. Build PWM from sample sequences ===")
  let seqs = @bio.seqlogo_sample_sequences()
  println("Loaded \{seqs.length()} aligned sequences of length \{seqs[0].length()}")
  let built_pwm = @bio.SeqLogoPwm::from_sequences(seqs)
  let built_consensus = @bio.seqlogo_consensus_sequence(built_pwm)
  let built_total_ic = @bio.seqlogo_total_information_content(built_pwm)
  println("Built PWM width: \{built_pwm.width()}")
  println("Built PWM total IC: \{built_total_ic} bits")
  println("Built consensus:  \{built_consensus}")
  println("Sample consensus: \{consensus}")
  // Compare consensus bases position-by-position between the two PWMs.
  let mut n_match = 0
  for pos in 0..<built_pwm.width() {
    let mut built_top = 0
    let mut built_max = -1.0
    for base in 0..<built_pwm.alphabet_size() {
      let f = built_pwm.get(base, pos)
      if f > built_max {
        built_max = f
        built_top = base
      }
    }
    let mut sample_top = 0
    let mut sample_max = -1.0
    for base in 0..<pwm.alphabet_size() {
      let f = pwm.get(base, pos)
      if f > sample_max {
        sample_max = f
        sample_top = base
      }
    }
    if built_top == sample_top {
      n_match = n_match + 1
    }
  }
  println(
    "Consensus positions matching: \{n_match} / \{built_pwm.width()}",
  )
  println("")

  // ------------------------------------------------------------------
  // 10. Show logo summary
  // ------------------------------------------------------------------
  println("=== 10. Logo summary ===")
  println(@bio.seqlogo_logo_summary(logo))

  println("=== seqLogo Demo complete ===")
}
