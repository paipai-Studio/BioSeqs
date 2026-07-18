///|
pub fn main() {
  println("=== DropletUtils Demo ===")
  
  let (counts, barcodes) = @src.create_example_droplet_data()
  
  println("\n1. Calculate droplet statistics:")
  let stats = @src.calculate_droplet_stats(counts, barcodes)
  println("  Number of barcodes: " + stats.total_counts.length().to_string())
  println("  Mean total counts: " + (stats.total_counts[0..10].reduce((a, b) => a + b, 0.0) / 10.0).to_string())
  
  println("\n2. Barcode ranking:")
  let ranked = @src.barcode_ranking(stats.total_counts)
  println("  Top 5 barcodes by UMI count:")
  let mut i = 0
  while i < 5 && i < ranked.length() {
    println("    " + barcodes[ranked[i].1] + ": " + ranked[i].0.to_string())
    i = i + 1
  }
  
  println("\n3. Find knee point:")
  let knee = @src.find_knee_point(ranked)
  println("  Knee point at index: " + knee.to_string())
  println("  Estimated cell number: " + @src.estimate_cell_number(stats.total_counts).to_string())
  
  println("\n4. Empty droplet detection (emptyDrops):")
  let result = @src.empty_drops(counts, barcodes, 10, 100)
  println("  Tested barcodes: " + result.barcode.length().to_string())
  
  let mut n_cells = 0
  i = 0
  while i < result.is_cell.length() {
    if result.is_cell[i] {
      n_cells = n_cells + 1
    }
    i = i + 1
  }
  println("  Detected cells: " + n_cells.to_string())
  
  println("\n5. Filter cells by emptyDrops results:")
  let (filtered_counts, filtered_barcodes) = @src.filter_cells_by_empty_drops(counts, barcodes, result)
  println("  Filtered barcodes: " + filtered_barcodes.length().to_string())
  
  println("\nDropletUtils demo completed!")
}