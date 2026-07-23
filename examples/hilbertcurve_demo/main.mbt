///|
/// HilbertCurve Demo - Demonstrates Hilbert curve coordinate mapping.

fn main {
  println("=== HilbertCurve Demo ===")
  
  let hc = @bio.HilbertCurve::new(4, 2)
  println("Levels: \{hc.levels}")
  println("Dimension: \{hc.dimension}")
  println("Max coordinate: \{hc.max_coordinate}")
  
  println("")
  println("Encode/decode roundtrip:")
  let coords = [3, 5]
  let encoded = @bio.hilbert_encode(hc, coords)
  let decoded = @bio.hilbert_decode(hc, encoded)
  println("  Coords [3,5] -> encoded: \{encoded} -> decoded: [\{decoded[0]},\{decoded[1]}]")
  
  println("")
  println("Distance between [1,2] and [3,4]:")
  let dist = @bio.hilbert_distance(hc, [1, 2], [3, 4])
  println("  Distance: \{dist}")
  
  println("")
  println("Linearize genome (10 positions):")
  let segments = @bio.hilbert_linearize_genome(10, 3)
  let mut i = 0
  while i < segments.length() && i < 5 {
    println("  Position \{segments[i].0} -> Hilbert index \{segments[i].1}")
    i = i + 1
  }
  
  println("")
  println("Map values to grid:")
  let values = [1.0, 2.0, 3.0, 4.0]
  let grid = @bio.hilbert_map_to_grid(@bio.HilbertCurve::new(2, 2), values)
  println("  Grid size: \{grid.length()}x\{grid.length()}")
}
