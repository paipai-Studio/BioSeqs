///|
/// Beachmat module demo - Bioconductor API for matrix access.
///
/// Run: moon run examples/beachmat_demo
fn main {
  println("=== Beachmat Module Demo ===")
  println("Bioconductor API for matrix access and block processing")
  println("")

  // ------------------------------------------------------------------
  // 1. Create sample data and Bmat
  // ------------------------------------------------------------------
  println("=== 1. Create sample data and Bmat ===")
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let param = @bio.BmatParam::column_param(2)
  let bmat = @bio.Bmat::new(data, param)
  println("Created Bmat with column block size = \{param.col_block_size()}")
  println(@bio.Bmat::to_string(bmat))
  println("")

  // ------------------------------------------------------------------
  // 2. Basic accessors
  // ------------------------------------------------------------------
  println("=== 2. Basic accessors ===")
  println("dim_rows: \{@bio.Bmat::dim_rows(bmat)}")
  println("dim_cols: \{@bio.Bmat::dim_cols(bmat)}")
  println("get(0, 0): \{@bio.Bmat::get(bmat, 0, 0)}")
  println("get(1, 2): \{@bio.Bmat::get(bmat, 1, 2)}")
  println("get(2, 1): \{@bio.Bmat::get(bmat, 2, 1)}")
  println("get_row(1): \{@bio.Bmat::get_row(bmat, 1)}")
  println("get_col(2): \{@bio.Bmat::get_col(bmat, 2)}")
  println("")

  // ------------------------------------------------------------------
  // 3. Column block processing
  // ------------------------------------------------------------------
  println("=== 3. Column block processing (block_size=2) ===")
  @bio.bmat_apply_col_blocks(bmat, 2, (block) => {
    println(
      "  Col block: cols \{@bio.BmatBlock::col_start(block)}-\{@bio.BmatBlock::col_end(block)} (\{@bio.BmatBlock::block_ncol(block)} cols x \{@bio.BmatBlock::block_nrow(block)} rows)",
    )
    let block_data = @bio.BmatBlock::get_data(block)
    println("    data: \{block_data}")
  })
  println("")

  // ------------------------------------------------------------------
  // 4. Row block processing
  // ------------------------------------------------------------------
  println("=== 4. Row block processing (block_size=2) ===")
  @bio.bmat_apply_row_blocks(bmat, 2, (block) => {
    println(
      "  Row block: rows \{@bio.BmatBlock::row_start(block)}-\{@bio.BmatBlock::row_end(block)} (\{@bio.BmatBlock::block_nrow(block)} rows x \{@bio.BmatBlock::block_ncol(block)} cols)",
    )
    let block_data = @bio.BmatBlock::get_data(block)
    println("    data: \{block_data}")
  })
  println("")

  // ------------------------------------------------------------------
  // 5. BmatIterator traversal
  // ------------------------------------------------------------------
  println("=== 5. BmatIterator traversal ===")
  let iter = @bio.BmatIterator::new(bmat)
  let mut count = 0
  while @bio.BmatIterator::has_next(iter) {
    let val = @bio.BmatIterator::next(iter)
    if val != None {
      let v = val.unwrap()
      let r = @bio.BmatIterator::cur_row(iter)
      let c = @bio.BmatIterator::cur_col(iter)
      println(
        "  [\{r}, \{c}] = \{v}",
      )
      count = count + 1
    }
  }
  println("Total elements traversed: \{count}")
  println("")

  // ------------------------------------------------------------------
  // 6. Transpose
  // ------------------------------------------------------------------
  println("=== 6. Transpose ===")
  let transposed = @bio.bmat_transpose(bmat)
  println("Original: \{@bio.Bmat::dim_rows(bmat)} x \{@bio.Bmat::dim_cols(bmat)}")
  println("Transposed: \{@bio.Bmat::dim_rows(transposed)} x \{@bio.Bmat::dim_cols(transposed)}")
  println(@bio.Bmat::to_string(transposed))
  println("")

  // ------------------------------------------------------------------
  // 7. Bind columns
  // ------------------------------------------------------------------
  println("=== 7. Bind columns ===")
  let extra = @bio.bmat_from_array([[10.0], [20.0], [30.0]])
  let bound_cols = @bio.bmat_bind_cols(bmat, extra)
  println("Left: \{@bio.Bmat::dim_rows(bmat)} x \{@bio.Bmat::dim_cols(bmat)}")
  println("Right: \{@bio.Bmat::dim_rows(extra)} x \{@bio.Bmat::dim_cols(extra)}")
  println("Bound: \{@bio.Bmat::dim_rows(bound_cols)} x \{@bio.Bmat::dim_cols(bound_cols)}")
  println(@bio.Bmat::to_string(bound_cols))
  println("")

  // ------------------------------------------------------------------
  // 8. Bind rows
  // ------------------------------------------------------------------
  println("=== 8. Bind rows ===")
  let extra_row = @bio.bmat_from_array([[100.0, 200.0, 300.0]])
  let bound_rows = @bio.bmat_bind_rows(bmat, extra_row)
  println("Top: \{@bio.Bmat::dim_rows(bmat)} x \{@bio.Bmat::dim_cols(bmat)}")
  println("Bottom: \{@bio.Bmat::dim_rows(extra_row)} x \{@bio.Bmat::dim_cols(extra_row)}")
  println("Bound: \{@bio.Bmat::dim_rows(bound_rows)} x \{@bio.Bmat::dim_cols(bound_rows)}")
  println(@bio.Bmat::to_string(bound_rows))
  println("")

  // ------------------------------------------------------------------
  // 9. Subset
  // ------------------------------------------------------------------
  println("=== 9. Subset (rows 0-2, cols 1-3) ===")
  let sub = @bio.bmat_subset(bmat, 0, 2, 1, 3)
  println("Subset: \{@bio.Bmat::dim_rows(sub)} x \{@bio.Bmat::dim_cols(sub)}")
  println(@bio.Bmat::to_string(sub))
  println("")

  // ------------------------------------------------------------------
  // 10. Element-wise operations
  // ------------------------------------------------------------------
  println("=== 10. Element-wise operations ===")
  let doubled = @bio.bmat_apply_elementwise(bmat, (x) => { x * 2.0 })
  println("Doubled:")
  println(@bio.Bmat::to_string(doubled))
  let squared = @bio.bmat_apply_elementwise(bmat, (x) => { x * x })
  println("Squared:")
  println(@bio.Bmat::to_string(squared))
  println("")

  println("=== Demo complete ===")
}