///|
/// Bio.Seq.Approximate module example - Approximate string matching.
///
/// Run: moon run examples/seq_approx_demo/main.mbt
fn main {
  println("=== Bio.Seq.Approximate Module Demo ===")
  println("")

  // Count mismatches
  println("=== Count Mismatches ===")
  println("ACGT vs ACGT: \{@bio.count_mismatches("ACGT", "ACGT")}")
  println("ACGT vs ACGA: \{@bio.count_mismatches("ACGT", "ACGA")}")
  println("ACGT vs TTTT: \{@bio.count_mismatches("ACGT", "TTTT")}")
  println("")

  // Exact search
  println("=== Approximate Search (Exact) ===")
  let m1 = @bio.approx_search("ACGT", "AACGTG", 0)
  println("Search 'ACGT' in 'AACGTG' (max 0 mismatches):")
  println("  Found: \{m1.is_match}")
  println("  Start: \{m1.start}, End: \{m1.end}")
  println("  Mismatches: \{m1.mismatches}")
  println("")

  // Search with mismatches
  println("=== Approximate Search (With Mismatches) ===")
  let m2 = @bio.approx_search("ACGT", "AACGAG", 1)
  println("Search 'ACGT' in 'AACGAG' (max 1 mismatch):")
  println("  Found: \{m2.is_match}")
  println("  Mismatches: \{m2.mismatches}")
  println("")

  // No match
  println("=== Approximate Search (No Match) ===")
  let m3 = @bio.approx_search("ACGT", "AAAAAA", 0)
  println("Search 'ACGT' in 'AAAAAA' (max 0 mismatches):")
  println("  Found: \{m3.is_match}")
  println("")

  // Levenshtein distance
  println("=== Levenshtein Distance ===")
  println("ACGT vs ACGT: \{@bio.levenshtein_distance("ACGT", "ACGT")}")
  println("ACGT vs ACGA: \{@bio.levenshtein_distance("ACGT", "ACGA")}")
  println("ACGT vs XXGT: \{@bio.levenshtein_distance("ACGT", "XXGT")}")
  println("ACGT vs ACGAT: \{@bio.levenshtein_distance("ACGT", "ACGAT")}")
  println("ACGT vs AGT: \{@bio.levenshtein_distance("ACGT", "AGT")}")
  println("ACGT vs '' (empty): \{@bio.levenshtein_distance("ACGT", "")}")
  println("'' vs '' (both empty): \{@bio.levenshtein_distance("", "")}")
  println("")

  // Approximate search with indels
  println("=== Approximate Search With Indels ===")
  let m4 = @bio.approx_search_with_indels("ACGT", "ACXGT", 1)
  println("Search 'ACGT' in 'ACXGT' (max 1 error):")
  println("  Found: \{m4.is_match}")
  println("  Errors: \{m4.mismatches}")
  println("")

  let m5 = @bio.approx_search_with_indels("ACGT", "ACGXT", 1)
  println("Search 'ACGT' in 'ACGXT' (max 1 error):")
  println("  Found: \{m5.is_match}")
  println("  Errors: \{m5.mismatches}")
  println("")

  let m6 = @bio.approx_search_with_indels("ACGT", "AGT", 1)
  println("Search 'ACGT' in 'AGT' (max 1 error, deletion):")
  println("  Found: \{m6.is_match}")
  println("  Errors: \{m6.mismatches}")
  println("")

  let m7 = @bio.approx_search_with_indels("AAAA", "TCGTCG", 1)
  println("Search 'AAAA' in 'TCGTCG' (max 1 error):")
  println("  Found: \{m7.is_match}")
  println("")

  // Find all approximate matches
  println("=== Find All Approximate Matches ===")
  let results = @bio.approx_find_all("ACG", "AACGACGTCG", 0)
  println("Find all 'ACG' in 'AACGACGTCG' (exact): \{results.length()} matches")
  let results2 = @bio.approx_find_all("ACG", "AACGACGTCG", 1)
  println("Find all 'ACG' in 'AACGACGTCG' (max 1 mismatch): \{results2.length()} matches")
  println("")

  // Approx word search
  println("=== Approx Word Search ===")
  let w1 = @bio.approx_word_search("ACGT", "AACGTG", 0)
  println("Word search 'ACGT' in 'AACGTG' (max 0 errors):")
  println("  Found: \{w1.is_found}")
  println("  Edit distance: \{w1.edit_distance}")
  println("  Match start: \{w1.match_start}")
  println("")

  let w2 = @bio.approx_word_search("ZZZZ", "ACGTACGT", 0)
  println("Word search 'ZZZZ' in 'ACGTACGT' (max 0 errors):")
  println("  Found: \{w2.is_found}")
  println("")

  // Approx best match
  println("=== Approx Best Match ===")
  let bm = @bio.approx_best_match("ACGT", "ACXGTACGT", 0.25)
  println("Best match 'ACGT' in 'ACXGTACGT' (max 25%% error rate):")
  println("  Found: \{bm.is_match}")
  println("  Mismatches: \{bm.mismatches}")
  println("")

  // Match to string
  println("=== Match String Representation ===")
  let mstr = @bio.approx_match_new("ACGT", "ACGT", 0, 4, 0, 4, true)
  println("Match: \{mstr.to_string()}")
  println("")

  println("=== Done ===")
}