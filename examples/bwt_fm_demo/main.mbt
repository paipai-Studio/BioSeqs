///|
/// Burrows-Wheeler Transform and FM-index demo.
///
/// Demonstrates:
/// 1. BWT forward and inverse transform
/// 2. Suffix array construction
/// 3. FM-index exact pattern matching (count + locate)
/// 4. LF mapping
///
/// Run: moon run examples/bwt_fm_demo/main.mbt
fn main {
  println("=== Burrows-Wheeler Transform & FM-index ===")
  println("")

  // --- 1. BWT on banana (classic example) ---
  println("1. BWT of 'banana':")
  let s = "banana"
  let bwt_str = @src.bwt(s)
  println("   Original: \{s}")
  println("   BWT:      \{bwt_str}")
  println("")

  // --- 2. Inverse BWT ---
  println("2. Inverse BWT:")
  let reconstructed = @src.bwt_inverse(bwt_str)
  println("   BWT:        \{bwt_str}")
  println("   Reconstructed: \{reconstructed}")
  println("   Match: \{reconstructed == s}")
  println("")

  // --- 3. Suffix Array ---
  println("3. Suffix Array of 'banana':")
  let sa = @src.suffix_array("banana")
  println("   Index  SA[i]  Suffix")
  for i = 0; i < sa.length(); i = i + 1 {
    let suffix = if sa[i] < s.length() { s[sa[i]:].to_owned() } else { "$" }
    println("   \{i}      \{sa[i]}      \{suffix}")
  }
  println("")

  // --- 4. FM-index Count ---
  println("4. FM-index Pattern Counting:")
  let fm = @src.FMIndex::new("banana")
  let patterns = ["a", "na", "ana", "banana", "xyz", "ban"]
  for i = 0; i < patterns.length(); i = i + 1 {
    let count = fm.count(patterns[i])
    println("   '\{patterns[i]}': \{count} occurrence(s)")
  }
  println("")

  // --- 5. FM-index Locate ---
  println("5. FM-index Pattern Locating:")
  let pattern = "ana"
  let positions = fm.locate(pattern)
  println("   Pattern: '\{pattern}'")
  println("   Found at positions:")
  for i = 0; i < positions.length(); i = i + 1 {
    let pos = positions[i]
    let context = if pos + pattern.length() <= s.length() {
      s[pos:pos + pattern.length()].to_owned()
    } else {
      s[pos:].to_owned()
    }
    println(
      "     position \{pos}: ...\{s[0:pos].to_owned()}[\{context}]\{s[pos + pattern.length():].to_owned()}...",
    )
  }
  println("")

  // --- 6. DNA Sequence Example ---
  println("6. DNA Sequence Alignment:")
  let dna = "ACGTACGTACGTACGT"
  let fm_dna = @src.FMIndex::new(dna)
  println("   Sequence: \{dna}")
  println("   Length: \{dna.length()}")
  println("   BWT length: \{fm_dna.len()}")
  let dna_patterns = ["ACGT", "CGTA", "TTTT", "AC", "GTAC"]
  for i = 0; i < dna_patterns.length(); i = i + 1 {
    let pat = dna_patterns[i]
    let cnt = fm_dna.count(pat)
    let pos = fm_dna.locate(pat)
    let pos_str = if pos.length() == 0 {
      "none"
    } else {
      let buf = StringBuilder::new()
      for j = 0; j < pos.length(); j = j + 1 {
        if j > 0 {
          buf.write_string(", ")
        }
        buf.write_string(pos[j].to_string())
      }
      buf.to_string()
    }
    println("   '\{pat}': \{cnt} times at [\{pos_str}]")
  }
  println("")

  // --- 7. LF Mapping ---
  println("7. LF Mapping Demonstration:")
  let fm2 = @src.FMIndex::new("banana")
  println("   BWT: \{fm2.bwt()}")
  println("   LF chain from position 3:")
  let mut cur = 3
  for step = 0; step < fm2.len(); step = step + 1 {
    let c = fm2.bwt().unsafe_get(cur).unsafe_to_char()
    println("     step \{step}: position \{cur}, char = '\{c}'")
    cur = fm2.lf(cur)
  }
  println("")

  // --- 8. Length preservation ---
  println("8. BWT Length Preservation:")
  let test_strings = ["", "A", "ACGT", "AAAAA", "GATTACA"]
  for i = 0; i < test_strings.length(); i = i + 1 {
    let t = test_strings[i]
    let b = @src.bwt(t)
    let inv = @src.bwt_inverse(b)
    let ok = if inv == t { "✓" } else { "✗" }
    println("   '\{t}': BWT len=\{b.length()}, inverse match: \{ok}")
  }
  println("")

  // --- 9. Overlapping patterns ---
  println("9. Overlapping Patterns (AAAAA):")
  let fm3 = @src.FMIndex::new("AAAAA")
  println("   Sequence: AAAAA")
  for k = 1; k <= 5; k = k + 1 {
    let pat = "AAAAA"[0:k].to_owned()
    let cnt = fm3.count(pat)
    println("   '\{pat}' (k=\{k}): \{cnt} occurrence(s)")
  }
  println("")

  println("=== BWT & FM-index Demo Complete ===")
}
