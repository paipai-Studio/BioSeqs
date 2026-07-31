///|
/// Biostrings sequence analysis example.
///
/// Demonstrates IUPAC handling, k-mer frequency, RSCU, complexity,
/// melting temperature, and sequence identity functions.
/// Also demonstrates advanced pattern matching with mismatches and indels.
/// Includes translation, reverse complement, and consensus sequence.
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
  println("=== Translation Example ===")
  let coding_seq = "ATGGCTTAACCAGT"
  println("Coding sequence: \{coding_seq}")
  let protein = @bio.biostrings_translate(coding_seq)
  println("Translated protein: \{protein}")

  println("")
  println("=== Reverse Complement Example ===")
  let dna_seq = "ATGC"
  println("DNA sequence: \{dna_seq}")
  let dna_rc = @bio.biostrings_reverse_complement(dna_seq)
  println("Reverse complement: \{dna_rc}")

  let iupac_seq = "AR"
  println("IUPAC sequence: \{iupac_seq}")
  let iupac_rc = @bio.biostrings_reverse_complement(iupac_seq)
  println("IUPAC reverse complement: \{iupac_rc}")

  println("")
  println("=== Letter Frequency Matrix Example ===")
  let matrix_seq = "ACGTACGT"
  println("Sequence: \{matrix_seq}")
  let matrix = @bio.letter_frequency_matrix(matrix_seq)
  let a_counts = matrix.get_or_default("A", Array::make(8, 0))
  println("A counts per position: \{a_counts}")

  println("")
  println("=== Consensus Sequence Example ===")
  let alignment = ["ACGT", "ACGT", "TCGA"]
  println("Alignment: \{alignment}")
  let consensus = @bio.biostrings_consensus_sequence(alignment)
  println("Consensus: \{consensus}")

  println("")
  println("=== Expected Matches Example ===")
  let expected = @bio.expected_matches("ACGT", 12, 0.5)
  println("Expected matches of ACGT in 12bp sequence (50% GC): \{expected}")

  println("")
  println("=== Sequence Complexity Example ===")
  let complex_seq = "ACGTACGTACGT"
  let complexity = @bio.sequence_complexity(complex_seq, 2)
  println("Sequence complexity (word_size=2): \{complexity}")

  println("")
  println("=== RSCU Example ===")
  let coding_seq2 = "ATGTGCTGAATGAA"
  let rscu = @bio.rscu(coding_seq2)
  println("Coding sequence: \{coding_seq2}")
  println("RSCU table size: \{rscu.length()} entries")

  println("")
  println("=== DUST Complexity Example ===")
  let simple_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  let dust_score = @bio.dust_complexity(simple_seq)
  println("Simple sequence DUST score: \{dust_score}")

  println("")
  println("=== Pattern Matching with Mismatches ===")
  let subject = "ACGTACGTACGTACGT"
  let pattern = "CGT"
  println("Subject: \{subject}")
  println("Pattern: \{pattern}")

  // Exact matching
  let result_exact = @bio.match_pattern(pattern=pattern, subject=subject)
  println("Exact match hits: \{result_exact.hits.length()}")
  let mut i = 0
  while i < result_exact.hits.length() {
    let hit = result_exact.hits[i]
    println("  Hit \{i}: start=\{hit.start}, end=\{hit.end}, mismatches=\{hit.mismatches}")
    i = i + 1
  }

  // With 1 mismatch
  let result_mm = @bio.match_pattern(
    pattern=pattern,
    subject=subject,
    max_mismatches=1,
  )
  println("")
  println("With max_mismatches=1 hits: \{result_mm.hits.length()}")

  // With indels
  let result_indel = @bio.match_pattern(
    pattern="ACG",
    subject="ACGTACGT",
    max_mismatches=1,
    with_indels=true,
  )
  println("")
  println("With indels hits: \{result_indel.hits.length()}")

  println("")
  println("=== Finding Palindromes ===")
  let palindrome_seq = "AAAAAAAAAACGTACGTGGGGGGGGGG"
  println("Sequence: \{palindrome_seq}")
  let palindromes = @bio.find_palindromes(seq=palindrome_seq, min_length=4)
  println("Palindromes found: \{palindromes.length()}")

  println("")
  println("=== Vectorized Pattern Matching ===")
  let subjects = ["ACGTACGT", "TGCATGCA", "AAAAAAAA"]
  println("Subjects array: \{subjects.length()} entries")
  let vmatch_results = @bio.vmatch_pattern(pattern="CGT", subjects=subjects)
  println("Vectorized match results: \{vmatch_results.length()} entries")
}
