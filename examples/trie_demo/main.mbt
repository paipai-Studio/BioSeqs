///|
/// Trie (prefix tree) demo - efficient prefix matching for bioinformatics.
///
/// Demonstrates the trie data structure for:
///   - Restriction enzyme site detection
///   - Primer matching in sequences
///   - Prefix-based search
///
/// Inspired by Biopython's Bio.trie and Bio.triefind modules.
/// Run: moon run examples/trie_demo
fn main {
  println("=== Trie Module Demo ===")
  println("Prefix tree for efficient string matching")
  println("Reference: Biopython Bio.trie / Bio.triefind")
  println("")

  // ------------------------------------------------------------------
  // 1. Restriction enzyme site detection
  // ------------------------------------------------------------------
  println("=== 1. Restriction Enzyme Site Detection ===")
  let trie = @bio.trie_sample_restriction_sites()
  println("Loaded restriction enzyme sites:")
  let items = @bio.trie_items(trie)
  let mut i = 0
  while i < items.length() {
    println("  " + items[i].0 + " -> " + items[i].1)
    i = i + 1
  }
  println("")

  // Scan a sequence for restriction sites
  let sequence = "ATGCGAATTCGGATCCGCGGCCGCTTAAAGAATTCTTT"
  println("Scanning sequence: " + sequence)
  let matches = @bio.triefind_find(trie, sequence)
  println("Found " + matches.length().to_string() + " restriction sites:")
  let mut m = 0
  while m < matches.length() {
    let (key, enzyme, start, end) = matches[m]
    println("  " + enzyme + " (" + key + ") at position " + start.to_string() + "-" + end.to_string())
    m = m + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Primer matching
  // ------------------------------------------------------------------
  println("=== 2. Primer Matching ===")
  let primer_trie = @bio.trie_sample_primers()
  println("Loaded primers:")
  let primers = @bio.trie_items(primer_trie)
  let mut p = 0
  while p < primers.length() {
    println("  " + primers[p].0 + " -> " + primers[p].1)
    p = p + 1
  }
  println("")

  // Find longest matching primer at sequence start
  let seq5 = "ATCGATCGATCGGGGCC"
  println("Matching primers at 5' end of: " + seq5)
  let result = @bio.triefind_match(primer_trie, seq5)
  match result {
    Some((key, val)) => println("  Best match: " + val + " (" + key + ")")
    None => println("  No match found")
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Prefix search
  // ------------------------------------------------------------------
  println("=== 3. Prefix-Based Search ===")
  let trie3 = @bio.Trie::new()
  @bio.trie_insert(trie3, "ATCGATCGAT", "Primer_Fwd_1")
  @bio.trie_insert(trie3, "ATCGATCGAC", "Primer_Fwd_2")
  @bio.trie_insert(trie3, "ATCGATCGTT", "Primer_Fwd_3")
  @bio.trie_insert(trie3, "GCTAGCTAGC", "Primer_Rev_1")
  println("Searching for keys with prefix 'ATCGATCG':")
  let prefix_results = @bio.trie_with_prefix(trie3, "ATCGATCG")
  let mut pr = 0
  while pr < prefix_results.length() {
    println("  " + prefix_results[pr].0 + " -> " + prefix_results[pr].1)
    pr = pr + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Trie operations
  // ------------------------------------------------------------------
  println("=== 4. Trie Operations ===")
  let trie4 = @bio.Trie::new()
  println("Inserting keys: AAA, AAB, AAC, BBB")
  @bio.trie_insert(trie4, "AAA", "val1")
  @bio.trie_insert(trie4, "AAB", "val2")
  @bio.trie_insert(trie4, "AAC", "val3")
  @bio.trie_insert(trie4, "BBB", "val4")
  println("Trie size: " + @bio.trie_size(trie4).to_string())
  println("Has key 'AAA': " + @bio.trie_has_key(trie4, "AAA").to_string())
  println("Has key 'XXX': " + @bio.trie_has_key(trie4, "XXX").to_string())
  println("Has prefix 'AA': " + @bio.trie_has_prefix(trie4, "AA").to_string())
  println("Has prefix 'XX': " + @bio.trie_has_prefix(trie4, "XX").to_string())
  println("")

  println("Removing key 'AAB'...")
  let removed = @bio.trie_remove(trie4, "AAB")
  println("Removed: " + removed.to_string())
  println("Trie size: " + @bio.trie_size(trie4).to_string())
  println("Has key 'AAB': " + @bio.trie_has_key(trie4, "AAB").to_string())
  println("")

  // ------------------------------------------------------------------
  // 5. Word matching
  // ------------------------------------------------------------------
  println("=== 5. Word Boundary Matching ===")
  let word_trie = @bio.Trie::new()
  @bio.trie_insert(word_trie, "EcoRI", "enzyme1")
  @bio.trie_insert(word_trie, "BamHI", "enzyme2")
  @bio.trie_insert(word_trie, "HindIII", "enzyme3")
  let text = "Use EcoRI and BamHI for DNA digestion, not HindIII today."
  println("Text: " + text)
  println("Word matches:")
  let words = @bio.triefind_find_words(word_trie, text)
  let mut w = 0
  while w < words.length() {
    let (key, val, start, end) = words[w]
    println("  " + val + " (" + key + ") at position " + start.to_string() + "-" + end.to_string())
    w = w + 1
  }
  println("")
  println("=== Trie Demo Complete ===")
}
