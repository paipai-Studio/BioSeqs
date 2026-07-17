///|
/// IRanges Demo.
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

  println("\n=== IRanges Demo Complete ===")
}
