///|
fn main {
  println("=== CAPS (Cleaved Amplified Polymorphic Sequences) Demo ===")
  println("Finding differential restriction enzyme cut sites in alignments.")
  println("")

  // 1. Load sample sequences and enzymes
  println("1. Load Sample Data:")
  println("====================")
  let seqs = @bio.caps_sample_sequences()
  let enzymes = @bio.caps_sample_enzymes()
  println("Sequences: " + seqs.length().to_string())
  for i = 0; i < seqs.length(); i = i + 1 {
    println("  seq" + i.to_string() + " (" + seqs[i].length().to_string() + " bp): " + seqs[i])
  }
  println("Enzymes: " + enzymes.length().to_string())
  for e in enzymes {
    println("  " + e.name + " (site=" + e.recognition_site + ", cut=" + e.cutting_position.to_string() + ")")
  }

  // 2. Build CAPS map
  println("\n2. Build CAPS Map:")
  println("==================")
  let m = @bio.caps_map(seqs, enzymes)
  println(m.to_string())

  // 3. View differential cutsites
  println("\n3. Differential Cutsites Found:")
  println("===============================")
  let dcuts = m.dcuts()
  println("Total: " + dcuts.length().to_string())
  for dc in dcuts {
    println("  " + dc.to_string())
  }

  // 4. Query by enzyme
  println("\n4. Query by Enzyme:")
  println("====================")
  let ecori_dcuts = m.get_dcuts_by_enzyme("EcoRI")
  println("EcoRI differential cutsites: " + ecori_dcuts.length().to_string())
  for dc in ecori_dcuts {
    println("  Position " + dc.start().to_string() + ": cuts seq " + dc.cuts_in()[0].to_string() + ", blocked in seq " + dc.blocked_in()[0].to_string())
  }

  // 5. Check specific enzyme
  println("\n5. Check Enzyme Presence:")
  println("=========================")
  println("Has EcoRI cutsites: " + m.has_dcuts_for_enzyme("EcoRI").to_string())
  println("Has HindIII cutsites: " + m.has_dcuts_for_enzyme("HindIII").to_string())

  // 6. Full report
  println("\n6. Full Report:")
  println("===============")
  println(m.report())

  // 7. Custom alignment example
  println("7. Custom Alignment Example:")
  println("============================")
  // Two sequences: one with HindIII site, one without
  let custom_seqs = [
    "GGATCCAAGCTTGGATCC", // has BamHI (GGATCC) and HindIII (AAGCTT)
    "GGATCCAAGCTTGGATCC", // identical - no differential
  ]
  let bamHI = @bio.restriction_enzyme("BamHI", "GGATCC", 1, "GATC", "")
  let custom_enzymes = [bamHI]
  let m2 = @bio.caps_map(custom_seqs, custom_enzymes)
  println("Custom alignment: " + m2.to_string())
  println("Differential cutsites: " + m2.dcut_count().to_string())
  println("(Expected: 0 - sequences are identical)")

  // 8. Three-way comparison
  println("\n8. Three-Way Comparison:")
  println("========================")
  let three_seqs = [
    "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", // has EcoRI site
    "AAAAAAAAGAGAATTTAAAAAAAAAAAAAA", // EcoRI blocked
    "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", // has EcoRI site
  ]
  let m3 = @bio.caps_map(three_seqs, enzymes)
  let dc3 = m3.get_dcuts_by_enzyme("EcoRI")
  if dc3.length() > 0 {
    let dc = dc3[0]
    println("EcoRI cutsite at position " + dc.start().to_string())
    println("  Cuts in: " + array_to_string(dc.cuts_in()))
    println("  Blocked in: " + array_to_string(dc.blocked_in()))
  }

  println("\n=== Demo Complete ===")
}

fn array_to_string(arr : Array[Int]) -> String {
  let mut s = "["
  let mut first = true
  for i in arr {
    if !first {
      s = s + ", "
    }
    first = false
    s = s + i.to_string()
  }
  s = s + "]"
  s
}
