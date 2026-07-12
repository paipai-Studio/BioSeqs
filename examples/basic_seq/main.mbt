///|
/// Basic sequence operations example.
///
/// Demonstrates DNA/RNA/Protein creation, complement, transcription,
/// translation, and common sequence utilities.
///
/// Run: moon run examples/basic_seq/main.mbt
fn main {
  // === Create sequences ===
  let dna = @bio.DNA::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  println("DNA: \{dna.to_string()} (len=\{dna.length()})")

  // === Complement & reverse complement ===
  let seq = @bio.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  println("Complement:         \{seq.complement().to_string()}")
  println("Reverse complement: \{seq.reverse_complement().to_string()}")

  // === Transcription (DNA -> RNA) ===
  let rna = seq.transcribe()
  println("Transcribed RNA: \{rna.to_string()}")

  // === Back transcription (RNA -> DNA) ===
  let back_dna = rna.back_transcribe()
  println("Back-transcribed: \{back_dna.to_string()}")

  // === Translation (RNA -> Protein) ===
  try {
    let protein = rna.translate()
    println("Translated protein: \{protein.to_string()}")

    let protein_stop = rna.translate(stop_symbol='@')
    println("With stop symbol:   \{protein_stop.to_string()}")

    let protein_trunc = rna.translate(to_stop=true)
    println("Translate to stop:  \{protein_trunc.to_string()}")
  } catch {
    _ => println("Translation error occurred")
  }

  // === Sequence search & count ===
  let s = @bio.Seq::new("ACGTACGTACGTACGTAAATTT")
  let pattern = "ACGT"
  println("\nSearch in: \{s.to_string()}")
  println("  count ACGT:         \{s.count(pattern)}")
  println("  count_overlap ACGT:  \{s.count_overlap(pattern)}")
  println("  find ACGT:          \{s.find(pattern)}")
  println("  rfind ACGT:         \{s.rfind(pattern)}")

  let pat_start = "ACGT"
  let pat_end = "TTT"
  let pat_contains = "AAA"
  println("  startswith ACGT:    \{s.startswith(pat_start)}")
  println("  endswith TTT:       \{s.endswith(pat_end)}")
  println("  contains AAA:       \{s.contains(pat_contains)}")

  // === Sequence manipulation ===
  let upper_seq = @bio.Seq::new("acgtACGT").upper()
  let lower_seq = @bio.Seq::new("ACGTacgt").lower()
  let rev_seq = @bio.Seq::new("ATGC").reverse()
  let sliced = s.slice(start=0, end=8)
  println("\nManipulation:")
  println("  upper:   \{upper_seq.to_string()}")
  println("  lower:   \{lower_seq.to_string()}")
  println("  reverse: \{rev_seq.to_string()}")
  println("  slice:   \{sliced.to_string()}")

  // === GC content & molecular weight ===
  let dna_str = dna.to_string()
  let gc_content = @bio.gc(dna_str)
  println("\nGC content: \{gc_content}")

  // === Melting temperature ===
  let tm_w = @bio.tm_wallace(dna_str)
  let tm_g = @bio.tm_gc(dna_str)
  println("Tm (Wallace): \{tm_w}")
  println("Tm (GC):      \{tm_g}")

  // === Checksums ===
  let crc = @bio.crc32(dna_str)
  let gcg_val = @bio.gcg(dna_str)
  let seg = @bio.seguid(dna_str)
  println("\nCRC32:  \{crc}")
  println("GCG:    \{gcg_val}")
  println("SEGUID: \{seg}")
}
