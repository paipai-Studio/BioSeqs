///|
/// SeqCode demo: Biopython Bio.SeqUtils helpers for amino-acid code
/// conversion, GC fraction / codon-position GC, IUPAC-aware nucleotide
/// search and six-frame translation.
///
/// Run: moon run examples/seqcode_demo
fn main {
  println("=== SeqCode Demo (Bio.SeqUtils advanced) ===\n")

  // 1. One-letter -> three-letter amino-acid codes (seq3)
  println("1. seq3: one-letter to three-letter codes")
  let protein = "MAIVMGRWKGAR*"
  println("   seq3(MAIVMGRWKGAR*) = \{@bio.seq3(protein)}")
  let star3 = @bio.seq3(protein, custom_map=Map([("*", "***")]))
  println("   seq3 with custom_map star->*** = \{star3}")
  let extended3 = @bio.seq3("BXZJUO*")
  println("   seq3(BXZJUO*) (IUPAC extended) = \{extended3}")

  // 2. Three-letter -> one-letter codes (seq1), case insensitive
  println("\n2. seq1: three-letter to one-letter codes")
  let three = "MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer"
  println("   seq1(MetAla..ArgTer) = \{@bio.seq1(three)}")
  let mixed_case = @bio.seq1("METalaIlEValMetGLYArgTer")
  println("   seq1 is case insensitive: \{mixed_case}")
  let extended1 = @bio.seq1("AsxGlxXleSecPylTer")
  println("   seq1(AsxGlxXleSecPylTer) = \{extended1}")

  // 3. GC fraction with ambiguous base handling (gc_fraction)
  println("\n3. gc_fraction: GC content between 0 and 1")
  let plain = @bio.gc_fraction("ACTG") catch { _ => 0.0 }
  println("   gc_fraction(ACTG) = \{plain}")
  let ign = @bio.gc_fraction("ACTGN", ambiguous="ignore") catch { _ => 0.0 }
  println("   gc_fraction(ACTGN, ignore) = \{ign}")
  let wgt = @bio.gc_fraction("ACTGN", ambiguous="weighted") catch { _ => 0.0 }
  println("   gc_fraction(ACTGN, weighted) = \{wgt}")
  let rem = @bio.gc_fraction("GDVV", ambiguous="remove") catch { _ => 0.0 }
  println("   gc_fraction(GDVV, remove) = \{rem}")

  // 4. GC content per codon position (GC123)
  println("\n4. gc123: total and per-codon-position GC percentages")
  let coding = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let r = @bio.gc123(coding)
  println("   sequence (39 nt), GC total = \{r.gc_all}%")
  println("   GC pos1 = \{r.gc1}%  GC pos2 = \{r.gc2}%  GC pos3 = \{r.gc3}%")

  // 5. IUPAC ambiguous nucleotide search (nt_search)
  println("\n5. nt_search: IUPAC-aware forward-strand search")
  let hits = @bio.nt_search("GATCGATCGA", "NGA") catch {
    _ => abort("valid IUPAC query should never fail")
  }
  println("   pattern for NGA = \{hits.pattern}")
  println("   positions in GATCGATCGA = \{join_positions(hits.positions)}")

  // 6. Six-frame translation with GC content
  println("\n6. six_frame_translations")
  let sft = @bio.six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA") catch {
    _ => abort("standard genetic code should never fail")
  }
  println(sft)

  println("=== Demo Complete ===")
}

///|
/// Render a list of integers as "1, 2, 3" for display.
fn join_positions(positions : Array[Int]) -> String {
  let buf = StringBuilder()
  for i = 0; i < positions.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(positions[i].to_string())
  }
  buf.to_string()
}
