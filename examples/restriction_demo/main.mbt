///|
/// Restriction enzyme analysis example inspired by Biopython Bio.Restriction.
///
/// Run: moon run examples/restriction_demo/main.mbt
fn main {
  println("=== Bio.Restriction Demo ===")

  println("\n1. List available enzymes:")
  let enzymes = @bio.list_enzymes()
  println("   Available enzymes: \{enzymes.length()}")
  let mut first5 = ""
  for i = 0; i < 5 && i < enzymes.length(); i = i + 1 {
    if i > 0 {
      first5 = first5 + ", "
    }
    first5 = first5 + enzymes[i]
  }
  println("   First 5: \{first5}")

  println("\n2. Get enzyme by name:")
  match @bio.get_enzyme_by_name("EcoRI") {
    Some(e) =>
      println(
        "   EcoRI: recognition_site=\{e.recognition_site}, cutting_position=\{e.cutting_position}",
      )
    None => println("   EcoRI not found")
  }

  let seq = "GAATTCGAAGCTTGAATTC"
  println("\n3. Sequence: \{seq}")

  println("\n4. Find cut sites with EcoRI:")
  let ecoRI = @bio.restriction_enzyme("EcoRI", "GAATTC", 1, "AATT", "")
  let sites = @bio.find_cut_sites(seq, ecoRI)
  println("   Cut sites: \{sites}")

  println("\n5. Digest sequence with EcoRI:")
  let fragments = @bio.digest(seq, ecoRI)
  println("   Fragments: \{fragments}")

  println("\n6. Fragment sizes:")
  let sizes = @bio.fragment_sizes(fragments)
  println("   Sizes: \{sizes}")

  println("\n7. Digest with multiple enzymes:")
  let hindIII = @bio.restriction_enzyme("HindIII", "AAGCTT", 1, "AGCT", "")
  let multi_fragments = @bio.digest_with_multiple(seq, [ecoRI, hindIII])
  println("   Fragments: \{multi_fragments}")

  println("\n8. Enzymes with 6-base recognition sites:")
  let six_base = @bio.enzymes_by_site_length(6, 6)
  println("   Count: \{six_base.length()}")

  println("\n9. Compatible overhangs:")
  let bamHI = @bio.restriction_enzyme("BamHI", "GGATCC", 1, "GATC", "")
  let bglII = @bio.restriction_enzyme("BglII", "AGATCT", 1, "GATC", "")
  println(
    "   BamHI & BglII compatible: \{@bio.is_compatible_overhang(bamHI, bglII)}",
  )
}
