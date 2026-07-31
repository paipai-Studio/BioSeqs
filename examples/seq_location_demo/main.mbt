///|
/// seq_location module example - Sequence location types for feature annotations.
///
/// Run: moon run examples/seq_location_demo
fn main {
  println("=== Bio.SeqLocation Module Demo ===")
  println("Sequence Location Types for Feature Annotations")
  println("")

  // Position types
  println("=== Position Types ===")
  let ep = @bio.exact_position(42)
  println("ExactPosition: " + ep.to_string() + " (exact=" + ep.is_exact().to_string() + ")")
  let bp = @bio.before_position(50)
  println("BeforePosition: " + bp.to_string() + " (exact=" + bp.is_exact().to_string() + ")")
  let ap = @bio.after_position(100)
  println("AfterPosition: " + ap.to_string() + " (exact=" + ap.is_exact().to_string() + ")")
  println("")

  // OneOf position
  println("=== OneOf Position ===")
  let opts : Array[@bio.Pos] = Array::new()
  opts.push(@bio.exact_position(5))
  opts.push(@bio.exact_position(7))
  opts.push(@bio.exact_position(9))
  let oop = @bio.one_of_position(7, opts)
  println("OneOfPosition: " + oop.to_string() + " (get=" + oop.get().to_string() + ")")
  println("")

  // Within position
  println("=== Within Position ===")
  let wp = @bio.within_position(50, @bio.exact_position(40), @bio.exact_position(60))
  println("WithinPosition: " + wp.to_string() + " (get=" + wp.get().to_string() + ")")
  println("")

  // Simple location
  println("=== Simple Location ===")
  let loc1 = @bio.simple_location(@bio.exact_position(10), @bio.exact_position(50), strand="+")
  println("Forward strand: " + loc1.to_string() + " (len=" + loc1.len().to_string() + ")")
  println("  Start: " + loc1.start().to_string() + ", End: " + loc1.end().to_string() + ", Strand: " + loc1.strand())
  let loc2 = @bio.simple_location(@bio.exact_position(100), @bio.exact_position(200), strand="-")
  println("Reverse strand: " + loc2.to_string() + " (len=" + loc2.len().to_string() + ")")
  println("")

  // Location operations
  println("=== Location Operations ===")
  let loc3 = @bio.simple_location(@bio.exact_position(30), @bio.exact_position(70))
  println("Overlaps check:")
  println("  loc1 (10-50) overlaps loc3 (30-70): " + loc1.overlaps(loc3).to_string())
  println("  loc1 (10-50) overlaps loc2 (100-200): " + loc1.overlaps(loc2).to_string())
  println("Contains check:")
  println("  loc1 contains 25: " + loc1.contains(25).to_string())
  println("  loc1 contains 100: " + loc1.contains(100).to_string())
  println("")

  // Compound location (exon structure)
  println("=== Compound Location (Gene with 3 exons) ===")
  let exon1 = @bio.simple_location(@bio.exact_position(0), @bio.exact_position(100), strand="+")
  let exon2 = @bio.simple_location(@bio.exact_position(200), @bio.exact_position(350), strand="+")
  let exon3 = @bio.simple_location(@bio.exact_position(500), @bio.exact_position(600), strand="+")
  let exon_locs : Array[@bio.Loc] = Array::new()
  exon_locs.push(exon1)
  exon_locs.push(exon2)
  exon_locs.push(exon3)
  let gene = @bio.compound_location(exon_locs, strand="+")
  println("Gene compound location: " + gene.to_string())
  println("  Is compound: " + gene.is_compound().to_string())
  println("  Number of parts: " + gene.parts().length().to_string())
  println("  Total length: " + gene.len().to_string())
  println("  GenBank format: " + gene.to_one_based())
  println("")

  // Compound from simple arrays
  println("=== Compound from Arrays ===")
  let starts : Array[Int] = [100, 300, 500, 700]
  let ends : Array[Int] = [200, 400, 600, 800]
  let cl = @bio.compound_from_simple(starts, ends, strand="-")
  println("Compound location from arrays:")
  println("  Parts: " + cl.parts().length().to_string())
  println("  Total length: " + cl.len().to_string())
  println("  Is compound: " + cl.is_compound().to_string())
  println("  GenBank: " + cl.to_one_based())
  println("")

  // GenBank conversion
  println("=== GenBank Coordinate Conversion ===")
  let gb_loc = @bio.from_genbank_coords(10, 50, strand="+")
  println("From GenBank (10, 50): 0-based start=" + gb_loc.start().to_string() + ", end=" + gb_loc.end().to_string())
  println("  GenBank output: " + gb_loc.to_one_based())
  println("")

  // Sample data
  println("=== Sample Locations ===")
  let samples = @bio.seq_location_sample()
  println("Number of sample locations: " + samples.length().to_string())
  let mut i = 0
  while i < samples.length() {
    println("  Location " + (i + 1).to_string() + ": " + samples[i].to_string() + " (compound=" + samples[i].is_compound().to_string() + ")")
    i = i + 1
  }
}
