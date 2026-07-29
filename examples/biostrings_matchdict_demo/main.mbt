///|
/// Biostrings PDict dictionary matching demo.
///
/// Demonstrates creating a pattern dictionary (PDict) with restriction enzyme
/// recognition patterns, matching patterns against DNA sequences, vectorized
/// counting across multiple subjects, and using which/which_index functions.
///
/// Run: moon run examples/biostrings_matchdict_demo/main.mbt
fn main {
  println("=== PDict Dictionary Matching Demo ===")
  println("")

  println("1. Creating a PDict with restriction enzyme patterns")
  // Restriction enzymes recognize specific short DNA sequences
  let patterns = ["GAATTC", "AAGCTT", "GGATCC", "GTCGAC"]
  let enzyme_names = ["EcoRI", "HindIII", "BamHI", "SalI"]
  let pdict = @bio.bmd_create_pdict(patterns=patterns)
  println("   Patterns loaded: \{pdict.patterns.length()}")
  let mut p = 0
  while p < patterns.length() {
    println("     \{enzyme_names[p]}: \{patterns[p]}")
    p = p + 1
  }

  println("")
  println("2. Matching patterns against a DNA sequence")
  let subject = "ATGCGAATTCTTAAGCTTGGATCCGTCGACATGC"
  println("   Subject: \{subject}")
  let match_result = @bio.bmd_match_pdict(pdict=pdict, subject=subject)
  println("   Total hits found: \{match_result.hits.length()}")
  let mut h = 0
  while h < match_result.hits.length() {
    let hit = match_result.hits[h]
    println(
      "     Hit \{h}: pattern=\"\{hit.pattern}\" (\{enzyme_names[hit.pattern_idx]}), start=\{hit.start}, end=\{hit.end}, mismatches=\{hit.mismatches}",
    )
    h = h + 1
  }

  println("")
  println("3. Matching with allowed mismatches")
  // Create a PDict allowing 1 mismatch
  let pdict_mm = @bio.bmd_create_pdict(patterns=patterns, max_mismatches=1)
  let mm_result = @bio.bmd_match_pdict(
    pdict=pdict_mm,
    subject="ATGCGAATTXCTTAAGCTTGGATCC",
  )
  println("   With max_mismatches=1, hits found: \{mm_result.hits.length()}")
  h = 0
  while h < mm_result.hits.length() {
    let hit = mm_result.hits[h]
    println(
      "     Hit \{h}: pattern=\"\{hit.pattern}\", start=\{hit.start}, mismatches=\{hit.mismatches}",
    )
    h = h + 1
  }

  println("")
  println("4. Vectorized countPattern across multiple sequences")
  let subjects = [
    "GAATTCGAATTCGAATTC",
    "AAGCTTGCATGCAAGCTT",
    "TTTTTTTTTTTTTTTTTT",
    "GAATTCAAGCTTGGATCC",
  ]
  println("   Subjects (\{subjects.length()} sequences):")
  let mut s = 0
  while s < subjects.length() {
    println("     \{s}: \{subjects[s]}")
    s = s + 1
  }

  let vcount_result = @bio.bmd_vcount_pattern(
    pattern="GAATTC",
    subjects=subjects,
  )
  println("")
  println("   vcountPattern for \"GAATTC\":")
  println("     Total occurrences: \{vcount_result.total}")
  s = 0
  while s < vcount_result.counts.length() {
    println("     Subject \{s}: \{vcount_result.counts[s]} occurrences")
    s = s + 1
  }

  println("")
  println("5. Vectorized matching of multiple patterns")
  let multi_patterns = ["GAATTC", "AAGCTT"]
  let multi_results = @bio.bmd_vmatch_pattern(
    patterns=multi_patterns,
    subjects=subjects,
  )
  let mut mp = 0
  while mp < multi_results.length() {
    let r = multi_results[mp]
    println("     Pattern \"\{r.pattern}\": total=\{r.total}, counts=\{r.counts}")
    mp = mp + 1
  }

  println("")
  println("6. which/which_index: finding subjects with matches")
  // Create a PDict with a single pattern for which/which_index
  let pdict_single = @bio.bmd_create_pdict(patterns=["GAATTC"])
  let which_result = @bio.bmd_which(pdict=pdict_single, subjects=subjects)
  println("   which() results (GAATTC in each subject):")
  s = 0
  while s < which_result.length() {
    println("     Subject \{s}: \{which_result[s]}")
    s = s + 1
  }

  let which_indices = @bio.bmd_which_index(pdict=pdict_single, subjects=subjects)
  println("   which_index() returns matching indices: \{which_indices}")

  println("")
  println("7. Counting total pattern occurrences")
  let total = @bio.bmd_count_occurrences(pdict=pdict, subject=subject)
  println("   Total pattern occurrences in subject: \{total}")

  println("")
  println("8. Finding best match (fewest mismatches)")
  let best = @bio.bmd_find_best_match(
    pdict=pdict_mm,
    subject="ATGCGAATTXCTTAAGCTTGGATCC",
  )
  match best {
    Some(hit) =>
      println(
        "   Best match: pattern=\"\{hit.pattern}\", mismatches=\{hit.mismatches}, start=\{hit.start}",
      )
    None => println("   No matches found")
  }

  println("")
  println("=== Demo Complete ===")
}