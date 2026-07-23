///|
/// MeltingTemp demo.
///
/// Run: moon run examples/melting_temp_demo/main.mbt
fn main {
  println("=== MeltingTemp Demo ===\n")

  let seq1 = "AGCT"
  let seq2 = "GGCC"
  let seq3 = "AGCTGATCGATCGATCGATCGATCGATCGATCGATCG"

  println("1. Wallace rule (short sequences):")
  let result1 = @bio.mt_wallace(seq1)
  println("   Seq: \{seq1}")
  println("   Tm: \{result1.tm} C")
  println("   Method: \{result1.method}")

  let result2 = @bio.mt_wallace(seq2)
  println("\n   Seq: \{seq2}")
  println("   Tm: \{result2.tm} C")
  println("   Method: \{result2.method}")

  println("\n2. GC content method (longer sequences):")
  let result3 = @bio.mt_gc_content(seq3)
  println("   Seq: \{seq3}")
  println("   Length: \{seq3.length()}")
  println("   Tm: \{result3.tm} C")
  println("   Method: \{result3.method}")

  println("\n3. Recommended method:")
  let method1 = @bio.recommend_tm_method(seq1)
  println("   Seq length \{seq1.length()}: recommended = \{method1}")
  let method2 = @bio.recommend_tm_method(seq3)
  println("   Seq length \{seq3.length()}: recommended = \{method2}")

  println("\n4. Example DNA:")
  let example_seq = @bio.create_example_dna()
  println("   Example DNA: \{example_seq}")

  println("\n=== Demo Complete ===")
}