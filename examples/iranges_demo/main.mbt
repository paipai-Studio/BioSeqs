///|
/// IRanges advanced operations demo.
///
/// Demonstrates integer range construction, interval operations,
/// overlap detection, distance calculations, coverage, and nearest
/// neighbor search for genomic intervals.
///
/// Run: moon run examples/iranges_demo/main.mbt
fn main {
  println("=== IRanges Demo ===\n")

  println("1. Creating IRanges...")
  let ir = @src.iranges([1, 5, 10], [10, 15, 20])
  println(
    "   Starts: \{@src.ir_start(ir).map(fn(x) { x.to_string() }).join(", ")}",
  )
  println("   Ends: \{@src.ir_end(ir).map(fn(x) { x.to_string() }).join(", ")}")
  println(
    "   Widths: \{@src.ir_width(ir).map(fn(x) { x.to_string() }).join(", ")}",
  )

  println("\n2. Shifting ranges by +3...")
  let shifted = @src.ir_shift(ir, 3)
  println(
    "   New starts: \{shifted.starts.map(fn(x) { x.to_string() }).join(", ")}",
  )

  println("\n3. Resizing to width 5...")
  let resized = @src.ir_resize(ir, 5, "start")
  println(
    "   New widths: \{resized.widths.map(fn(x) { x.to_string() }).join(", ")}",
  )

  println("\n4. Reducing overlapping ranges...")
  let ir2 = @src.iranges([1, 3, 8], [5, 6, 10])
  let reduced = @src.ir_reduce(ir2, 0)
  println("   Reduced to \{reduced.starts.length()} ranges")

  println("\n5. Intersecting two IRanges...")
  let ir3 = @src.iranges([1, 8], [10, 15])
  let ir4 = @src.iranges([5, 12], [12, 20])
  let intersect = @src.ir_intersect(ir3, ir4)
  println("   Intersections: \{intersect.starts.length()}")

  println("\n6. Union of two IRanges...")
  let union = @src.ir_union(ir3, ir4)
  println("   Union ranges: \{union.starts.length()}")

  println("\n7. Finding overlaps...")
  let ir5 = @src.iranges([1, 20], [10, 30])
  let ir6 = @src.iranges([5, 25], [15, 35])
  let overlaps = @src.ir_overlaps(ir5, ir6)
  println("   Overlapping pairs: \{overlaps.length()}")

  println("\n8. Advanced overlap types...")
  let ir_query = @src.iranges([1, 10], [5, 15])
  let ir_subject = @src.iranges([3, 8, 12], [8, 14, 18])

  let any_overlaps = @src.find_overlaps_advanced(query=ir_query, subject=ir_subject, type_="any")
  println("   Any overlap count: \{any_overlaps.length()}")

  let within_overlaps = @src.find_overlaps_advanced(query=ir_query, subject=ir_subject, type_="within")
  println("   Within overlap count: \{within_overlaps.length()}")

  println("\n9. IRanges nearest...")
  let ir_n1 = @src.iranges([1, 20, 40], [5, 25, 45])
  let ir_n2 = @src.iranges([10, 30], [15, 35])
  let nearest_indices = @src.ir_nearest(ir_n1, ir_n2)
  println("   Nearest indices: \{nearest_indices}")
  let mut i = 0
  while i < nearest_indices.length() {
    println("     Query \{i} nearest to Subject \{nearest_indices[i]}")
    i = i + 1
  }

  println("\n10. IRanges coverage...")
  let ir_cov = @src.iranges([1, 5, 10], [8, 15, 20])
  let coverage = @src.ir_coverage(ir_cov, 25)
  println("   Coverage array (length 25):")
  println("   \{coverage}")

  println("\n11. IRanges distance matrix...")
  let ir_d1 = @src.iranges([1, 15], [5, 20])
  let ir_d2 = @src.iranges([10, 30], [12, 35])
  let distance_matrix = @src.ir_distance(ir_d1, ir_d2)
  println("   Distance matrix (2x2):")
  let mut row = 0
  while row < distance_matrix.length() {
    println("     Row \{row}: \{distance_matrix[row]}")
    row = row + 1
  }

  println("\n12. IRanges distance to nearest...")
  let ir_dn = @src.iranges([1, 10, 25], [5, 15, 30])
  let distances_to_nearest = @src.ir_distance_to_nearest(ir_dn)
  println("   Distance to nearest for each range:")
  i = 0
  while i < ir_dn.starts.length() {
    println("     Range \{i} (\{ir_dn.starts[i]}-\{ir_dn.ends[i]}): distance = \{distances_to_nearest[i]}")
    i = i + 1
  }

  println("\n=== IRanges Demo Complete ===")
}
