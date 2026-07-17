///|
/// Suffix Array and Suffix Tree Demo
fn main {
  println("=== Suffix Array and Suffix Tree Demo ===")

  println("\n1. Suffix Array - Basic Operations:")
  println("===================================")
  let sa = @src.SuffixArray::new("banana")
  println("Text: " + sa.text())
  println("Suffix array length: " + sa.len().to_string())
  println("\nSorted suffixes:")
  for i = 0; i < sa.len(); i = i + 1 {
    println("  [" + i.to_string() + "] " + sa.get_sorted_suffix(i))
  }

  println("\n2. Suffix Array - Pattern Matching:")
  println("===================================")
  println("Count of 'ana': " + sa.count("ana").to_string())
  println("Count of 'an': " + sa.count("an").to_string())
  println("Count of 'b': " + sa.count("b").to_string())
  println("Count of 'xyz': " + sa.count("xyz").to_string())

  let pos = sa.locate("ana")
  let pos_str = if pos.length() == 0 {
    "none"
  } else {
    let buf = StringBuilder::new()
    for i = 0; i < pos.length(); i = i + 1 {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(pos[i].to_string())
    }
    buf.to_string()
  }
  println("Positions of 'ana': " + pos_str)

  println("\n3. LCP Array (Longest Common Prefix):")
  println("=====================================")
  let lcp = @src.LCPArray::new(sa)
  println("LCP values:")
  let lcp_values = lcp.lcp()
  for i = 0; i < lcp_values.length(); i = i + 1 {
    println("  [" + i.to_string() + "] " + lcp_values[i].to_string())
  }
  println("Max LCP: " + lcp.max_lcp().to_string())

  println("\n4. Suffix Tree - Basic Operations:")
  println("===================================")
  let st = @src.SuffixTree::new("banana")
  println("Text: " + st.text())
  println("Tree size (nodes): " + st.size().to_string())
  println("Contains 'ana': " + st.contains("ana").to_string())
  println("Contains 'ban': " + st.contains("ban").to_string())
  println("Contains 'xyz': " + st.contains("xyz").to_string())

  println("\n5. Suffix Tree - Pattern Matching:")
  println("==================================")
  println("Count of 'ana': " + st.count("ana").to_string())
  println("Count of 'an': " + st.count("an").to_string())

  let st_pos = st.locate("ana")
  let st_pos_str = if st_pos.length() == 0 {
    "none"
  } else {
    let buf = StringBuilder::new()
    for i = 0; i < st_pos.length(); i = i + 1 {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(st_pos[i].to_string())
    }
    buf.to_string()
  }
  println("Positions of 'ana': " + st_pos_str)

  println("\n6. Suffix Tree - Longest Repeated Substring:")
  println("=============================================")
  println("LRS of 'banana': " + st.longest_repeated_substring())

  let st_dna = @src.SuffixTree::new("ACGTACGT")
  println("LRS of 'ACGTACGT': " + st_dna.longest_repeated_substring())

  println("\n7. DNA Sequence Example:")
  println("=========================")
  let dna_sa = @src.SuffixArray::new("ACGTACGT")
  let dna_st = @src.SuffixTree::new("ACGTACGT")
  println("Text: ACGTACGT")
  println("SA count of 'ACGT': " + dna_sa.count("ACGT").to_string())
  println("ST count of 'ACGT': " + dna_st.count("ACGT").to_string())

  let dna_pos = dna_sa.locate("ACGT")
  let dna_pos_str = if dna_pos.length() == 0 {
    "none"
  } else {
    let buf = StringBuilder::new()
    for i = 0; i < dna_pos.length(); i = i + 1 {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(dna_pos[i].to_string())
    }
    buf.to_string()
  }
  println("SA positions of 'ACGT': " + dna_pos_str)

  println("\n8. Suffix Tree Structure Output:")
  println("================================")
  let st_abc = @src.SuffixTree::new("abc")
  println(st_abc.to_string())

  println("\n=== Demo Complete ===")
}
