///|
test "s4vectors_rle_from_vector" {
  let vec = ["A", "A", "A", "B", "B", "C", "C", "C", "C", "A", "A"]
  let rle = @src.Rle::from_vector(vec)
  
  assert_eq(rle.values.length(), 4)
  assert_eq(rle.lengths.length(), 4)
}

///|
test "s4vectors_rle_length" {
  let vec = ["A", "A", "B", "B", "B"]
  let rle = @src.Rle::from_vector(vec)
  
  assert_eq(rle.length(), 5)
}

///|
test "s4vectors_rle_get" {
  let vec = ["A", "A", "B", "C", "C", "C"]
  let rle = @src.Rle::from_vector(vec)
  
  assert_eq(rle.get(0), "A")
  assert_eq(rle.get(2), "B")
  assert_eq(rle.get(5), "C")
}

///|
test "s4vectors_dataframe_basic" {
  let col1 = @src.S4DataFrameColumn::new("gene_id", ["gene1", "gene2", "gene3"])
  let col2 = @src.S4DataFrameColumn::new("expression", ["10.5", "25.3", "5.8"])
  
  let df = @src.S4DataFrame::new([col1, col2], ["row1", "row2", "row3"])
  
  assert_eq(df.nrow(), 3)
  assert_eq(df.ncol(), 2)
}

///|
test "s4vectors_dataframe_colnames" {
  let col1 = @src.S4DataFrameColumn::new("gene_id", ["gene1", "gene2"])
  let col2 = @src.S4DataFrameColumn::new("expression", ["10.5", "25.3"])
  
  let df = @src.S4DataFrame::new([col1, col2], [])
  let names = df.colnames()
  
  assert_eq(names.length(), 2)
}

///|
test "s4vectors_hits_basic" {
  let hits = @src.Hits::new([0, 0, 1, 2], [1, 2, 0, 1], 3, 3)
  
  assert_eq(hits.n_hits(), 4)
}

///|
test "s4vectors_hits_count_query_hits" {
  let hits = @src.Hits::new([0, 0, 1, 2], [1, 2, 0, 1], 3, 3)
  
  let counts = hits.count_query_hits()
  
  assert_eq(counts.length(), 3)
}