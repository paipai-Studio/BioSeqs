///|
/// CircSeq Demo - Demonstrates circular DNA sequence operations.
/// Inspired by Biopython's Bio.SeqUtils.CircSeq module.
///
/// Run with: moon run examples/circ_seq_demo/main.mbt
fn main {
  println("=== CircSeq Demo ===\n")

  // 1. Create a circular sequence
  println("1. Creating circular sequence:")
  let circ = @bio.CircSeq::new(sequence="ATCGATCGATCGATCG")
  println("   Length: \{circ.circ_length()}")
  println("   Name: \{circ.circ_name()}")
  println("   GC Content: \{circ.circ_gc_content()}")

  // 2. Extract subsequence (circular coordinates)
  println("\n2. Circular subsequence extraction:")
  let subseq = circ.circ_subseq(14, 4)
  println("   Subseq(14, 4): \{subseq}")

  let slice = circ.circ_slice(14, 6)
  println("   Slice(14, 6): \{slice}")

  // 3. Rotate sequence
  println("\n3. Rotating circular sequence:")
  let rotated = circ.circ_rotate(4)
  println("   Rotated (origin=4): \{rotated.circ_sequence()}")

  // 4. Reverse complement
  println("\n4. Reverse complement:")
  let rc = circ.circ_reverse_complement()
  println("   Original: \{circ.circ_sequence()}")
  println("   RC: \{rc.circ_sequence()}")

  // 5. Find restriction sites
  println("\n5. Finding restriction sites:")
  let sites = circ.circ_find_sites("GATC")
  let mut s = 0
  while s < sites.length() {
    let (pos, site) = sites[s]
    println("   Site at \{pos}: \{site}")
    s = s + 1
  }

  // 6. Simulate restriction digest
  println("\n6. Restriction digest:")
  let fragments = circ.circ_digest(["GATC", "ATCG"])
  let mut f = 0
  while f < fragments.length() {
    println("   Fragment: \{fragments[f]}")
    f = f + 1
  }

  // 7. Circular distance
  println("\n7. Circular distance calculations:")
  let d1 = @bio.circular_distance(14, 4, 16)
  println("   Distance(14->4, len=16): \{d1}")

  // 8. Insert and delete
  println("\n8. Sequence manipulation:")
  let inserted = circ.circ_insert(2, "GGGG")
  println("   After insert at 2: \{inserted.circ_sequence()}")

  let deleted = circ.circ_delete(4, 7)
  println("   After delete(4,7): \{deleted.circ_sequence()}")

  // 9. Random circular sequence
  println("\n9. Random sequence generation:")
  let random_circ = @bio.circ_random_sequence(50, 0.5)
  println("   Random seq length: \{random_circ.circ_length()}")
  println("   GC content: \{random_circ.circ_gc_content()}")

  // 10. Common restriction sites
  println("\n10. Common restriction enzymes:")
  let common = @bio.circ_common_restriction_sites()
  let mut c = 0
  while c < common.length() && c < 5 {
    let (name, site) = common[c]
    println("   \{name}: \{site}")
    c = c + 1
  }

  println("\n=== Demo Complete ===")
}