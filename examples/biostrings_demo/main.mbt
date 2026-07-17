///|
/// Biostrings sequence analysis example.
///
/// Demonstrates IUPAC handling, k-mer frequency, RSCU, complexity,
/// melting temperature, and sequence identity functions.
///
/// Run: moon run examples/biostrings_demo/main.mbt
fn main {
  println("=== Biostrings Demo ===")
  println("")

  let seq = "ACGTMRWSYKVHDBN"
  println("Sequence: \{seq}")

  let freq = @bio.letter_frequency_dna(seq)
  println(
    "Letter Frequency: A=\{freq.a}, C=\{freq.c}, G=\{freq.g}, T=\{freq.t}",
  )

  let di_freq = @bio.oligonucleotide_frequency(seq, 2, false)
  println("Di-nucleotide Frequency: AC=\{di_freq.get("AC").unwrap_or(0.0)}")

  let gc_content = @bio.gc_content_by_position(seq, 1)
  println("GC Content by position: \{gc_content[0]}")

  let entropy = @bio.shannon_entropy(seq)
  println("Shannon Entropy: \{entropy}")

  let rc = @bio.iupac_reverse_complement(seq)
  println("IUPAC Reverse Complement: \{rc}")

  let tm = @bio.tm_nn("AGCTAGCT")
  println("Nearest-neighbor Tm: \{tm}")

  let pid = @bio.percent_identity("AGCT", "AGTT")
  println("Percent Identity: \{pid}%")

  let edits = @bio.nedit_at("AGCT", "AGTT", true)
  println("Number of edits: \{edits}")

  println("")
  println("=== RSCU Example ===")
  let coding_seq = "ATGTGCTGAATGAA"
  let rscu = @bio.rscu(coding_seq)
  println("Coding sequence: \{coding_seq}")
  println("RSCU table size: \{rscu.length()} entries")

  println("")
  println("=== DUST Complexity Example ===")
  let simple_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  let dust_score = @bio.dust_complexity(simple_seq)
  println("Simple sequence DUST score: \{dust_score}")
}
