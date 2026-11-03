///|
/// S4Vectors module demo.

///|
pub fn main() {
  println("=== S4Vectors Demo ===")
  
  println("\n1. Rle (Run-Length Encoding)")
  let vec = ["A", "A", "A", "B", "B", "C", "C", "C", "C", "A", "A"]
  let rle = @src.Rle::from_vector(vec)
  println("Original vector length: {}", vec.length())
  println("Rle compressed values: {}", rle.values.length())
  println("Rle lengths: {}", rle.lengths.length())
  println("Rle total length: {}", rle.length())
  
  println("\n2. S4DataFrame")
  let col1 = @src.S4DataFrameColumn::new("gene_id", ["gene1", "gene2", "gene3", "gene4", "gene5"])
  let col2 = @src.S4DataFrameColumn::new("expression", ["10.5", "25.3", "5.8", "100.2", "15.7"])
  let col3 = @src.S4DataFrameColumn::new("p_value", ["0.001", "0.05", "0.8", "0.0001", "0.1"])
  
  let df = @src.S4DataFrame::new([col1, col2, col3], ["row1", "row2", "row3", "row4", "row5"])
  println("DataFrame rows: {}", df.nrow())
  println("DataFrame cols: {}", df.ncol())
  println("Column names: {}", df.colnames().length())
  
  println("\n3. Hits")
  let hits = @src.Hits::new([0, 0, 1, 2, 2], [1, 2, 0, 1, 2], 3, 3)
  println("Number of hits: {}", hits.length())
  
  println("\n=== Demo complete ===")
}