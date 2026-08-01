///|
/// Tests for beachmat module.

// ---------------------------------------------------------------------------
// BmatParam creation and accessors
// ---------------------------------------------------------------------------

test "bmat_param_creation" {
  let p = @src.BmatParam::new()
  assert_eq(p.row_block_size(), 0)
  assert_eq(p.col_block_size(), 0)
}

test "bmat_param_with_blocks" {
  let p = @src.BmatParam::with_blocks(10, 20)
  assert_eq(p.row_block_size(), 10)
  assert_eq(p.col_block_size(), 20)
}

test "bmat_param_column_param" {
  let p = @src.BmatParam::column_param(15)
  assert_eq(p.row_block_size(), 0)
  assert_eq(p.col_block_size(), 15)
}

test "bmat_param_row_param" {
  let p = @src.BmatParam::row_param(25)
  assert_eq(p.row_block_size(), 25)
  assert_eq(p.col_block_size(), 0)
}

test "bmat_param_accessors" {
  let p = @src.BmatParam::with_blocks(5, 10)
  assert_eq(p.row_block_size(), 5)
  assert_eq(p.col_block_size(), 10)
}

// ---------------------------------------------------------------------------
// Bmat creation and basic accessors
// ---------------------------------------------------------------------------

test "bmat_creation_from_array" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 3)
  assert_eq(bmat.dim_cols(), 3)
}

test "bmat_creation_zeros" {
  let bmat = @src.Bmat::zeros(4, 5)
  assert_eq(bmat.dim_rows(), 4)
  assert_eq(bmat.dim_cols(), 5)
  assert_eq(bmat.get(0, 0), 0.0)
  assert_eq(bmat.get(3, 4), 0.0)
}

test "bmat_dim_rows" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 3)
}

test "bmat_dim_cols" {
  let data = [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_cols(), 4)
}

test "bmat_get" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.get(0, 0), 1.0)
  assert_eq(bmat.get(0, 1), 2.0)
  assert_eq(bmat.get(0, 2), 3.0)
  assert_eq(bmat.get(1, 0), 4.0)
  assert_eq(bmat.get(1, 1), 5.0)
  assert_eq(bmat.get(1, 2), 6.0)
  assert_eq(bmat.get(2, 0), 7.0)
  assert_eq(bmat.get(2, 1), 8.0)
  assert_eq(bmat.get(2, 2), 9.0)
}

test "bmat_set" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  bmat.set(0, 0, 10.0)
  bmat.set(1, 1, 20.0)
  assert_eq(bmat.get(0, 0), 10.0)
  assert_eq(bmat.get(0, 1), 2.0)
  assert_eq(bmat.get(1, 0), 3.0)
  assert_eq(bmat.get(1, 1), 20.0)
}

test "bmat_get_col" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col0 = bmat.get_col(0)
  assert_eq(col0.length(), 3)
  assert_eq(col0[0], 1.0)
  assert_eq(col0[1], 4.0)
  assert_eq(col0[2], 7.0)
  let col1 = bmat.get_col(1)
  assert_eq(col1[0], 2.0)
  assert_eq(col1[1], 5.0)
  assert_eq(col1[2], 8.0)
  let col2 = bmat.get_col(2)
  assert_eq(col2[0], 3.0)
  assert_eq(col2[1], 6.0)
  assert_eq(col2[2], 9.0)
}

test "bmat_get_row" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row0 = bmat.get_row(0)
  assert_eq(row0.length(), 3)
  assert_eq(row0[0], 1.0)
  assert_eq(row0[1], 2.0)
  assert_eq(row0[2], 3.0)
  let row1 = bmat.get_row(1)
  assert_eq(row1[0], 4.0)
  assert_eq(row1[1], 5.0)
  assert_eq(row1[2], 6.0)
  let row2 = bmat.get_row(2)
  assert_eq(row2[0], 7.0)
  assert_eq(row2[1], 8.0)
  assert_eq(row2[2], 9.0)
}

// ---------------------------------------------------------------------------
// BmatBlock creation and properties
// ---------------------------------------------------------------------------

test "bmat_block_properties" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_starts : Array[Int] = []
  let row_ends : Array[Int] = []
  let col_starts : Array[Int] = []
  let col_ends : Array[Int] = []
  let block_nrows : Array[Int] = []
  let block_ncols : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 2, fn(block) {
    row_starts.push(block.row_start())
    row_ends.push(block.row_end())
    col_starts.push(block.col_start())
    col_ends.push(block.col_end())
    block_nrows.push(block.block_nrow())
    block_ncols.push(block.block_ncol())
  })
  assert_eq(row_starts.length(), 2)
  assert_eq(row_starts[0], 0)
  assert_eq(row_ends[0], 3)
  assert_eq(block_nrows[0], 3)
  assert_eq(col_starts[0], 0)
  assert_eq(col_ends[0], 2)
  assert_eq(block_ncols[0], 2)
  assert_eq(col_starts[1], 2)
  assert_eq(col_ends[1], 3)
  assert_eq(block_ncols[1], 1)
}

test "bmat_block_get_data" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let block0_vals : Array[Double] = []
  let block1_vals : Array[Double] = []
  let mut block_idx = 0
  @src.bmat_apply_col_blocks(bmat, 2, fn(block) {
    let bd = block.get_data()
    if block_idx == 0 {
      for i in 0..<bd.length() {
        for j in 0..<bd[i].length() {
          block0_vals.push(bd[i][j])
        }
      }
    } else {
      for i in 0..<bd.length() {
        for j in 0..<bd[i].length() {
          block1_vals.push(bd[i][j])
        }
      }
    }
    block_idx = block_idx + 1
  })
  assert_eq(block0_vals.length(), 6)
  assert_eq(block0_vals[0], 1.0)
  assert_eq(block0_vals[1], 2.0)
  assert_eq(block0_vals[2], 4.0)
  assert_eq(block0_vals[3], 5.0)
  assert_eq(block0_vals[4], 7.0)
  assert_eq(block0_vals[5], 8.0)
  assert_eq(block1_vals.length(), 3)
  assert_eq(block1_vals[0], 3.0)
  assert_eq(block1_vals[1], 6.0)
  assert_eq(block1_vals[2], 9.0)
}

// ---------------------------------------------------------------------------
// Block processing
// ---------------------------------------------------------------------------

test "bmat_apply_col_blocks_basic" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col_starts : Array[Int] = []
  let col_ends : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 2, fn(block) {
    col_starts.push(block.col_start())
    col_ends.push(block.col_end())
  })
  assert_eq(col_starts.length(), 2)
  assert_eq(col_starts[0], 0)
  assert_eq(col_ends[0], 2)
  assert_eq(col_starts[1], 2)
  assert_eq(col_ends[1], 3)
}

test "bmat_apply_col_blocks_three_blocks" {
  let data = [[1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col_starts : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 2, fn(block) {
    col_starts.push(block.col_start())
  })
  assert_eq(col_starts.length(), 3)
  assert_eq(col_starts[0], 0)
  assert_eq(col_starts[1], 2)
  assert_eq(col_starts[2], 4)
}

test "bmat_apply_col_blocks_verification" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let all_vals : Array[Double] = []
  @src.bmat_apply_col_blocks(bmat, 2, fn(block) {
    let bd = block.get_data()
    for i in 0..<bd.length() {
      for j in 0..<bd[i].length() {
        all_vals.push(bd[i][j])
      }
    }
  })
  assert_eq(all_vals.length(), 9)
  assert_eq(all_vals[0], 1.0)
  assert_eq(all_vals[1], 2.0)
  assert_eq(all_vals[2], 4.0)
  assert_eq(all_vals[3], 5.0)
  assert_eq(all_vals[4], 7.0)
  assert_eq(all_vals[5], 8.0)
  assert_eq(all_vals[6], 3.0)
  assert_eq(all_vals[7], 6.0)
  assert_eq(all_vals[8], 9.0)
}

test "bmat_apply_row_blocks_basic" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_starts : Array[Int] = []
  let row_ends : Array[Int] = []
  @src.bmat_apply_row_blocks(bmat, 2, fn(block) {
    row_starts.push(block.row_start())
    row_ends.push(block.row_end())
  })
  assert_eq(row_starts.length(), 2)
  assert_eq(row_starts[0], 0)
  assert_eq(row_ends[0], 2)
  assert_eq(row_starts[1], 2)
  assert_eq(row_ends[1], 3)
}

test "bmat_apply_row_blocks_verification" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let all_vals : Array[Double] = []
  @src.bmat_apply_row_blocks(bmat, 2, fn(block) {
    let bd = block.get_data()
    for i in 0..<bd.length() {
      for j in 0..<bd[i].length() {
        all_vals.push(bd[i][j])
      }
    }
  })
  assert_eq(all_vals.length(), 9)
  assert_eq(all_vals[0], 1.0)
  assert_eq(all_vals[1], 2.0)
  assert_eq(all_vals[2], 3.0)
  assert_eq(all_vals[3], 4.0)
  assert_eq(all_vals[4], 5.0)
  assert_eq(all_vals[5], 6.0)
  assert_eq(all_vals[6], 7.0)
  assert_eq(all_vals[7], 8.0)
  assert_eq(all_vals[8], 9.0)
}

test "bmat_apply_row_blocks_three_blocks" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 10.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_starts : Array[Int] = []
  @src.bmat_apply_row_blocks(bmat, 2, fn(block) {
    row_starts.push(block.row_start())
  })
  assert_eq(row_starts.length(), 3)
  assert_eq(row_starts[0], 0)
  assert_eq(row_starts[1], 2)
  assert_eq(row_starts[2], 4)
}

test "bmat_foreach" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let values : Array[Double] = []
  @src.bmat_foreach(bmat, fn(row, col, val) {
    values.push(val)
  })
  assert_eq(values.length(), 6)
  assert_eq(values[0], 1.0)
  assert_eq(values[1], 2.0)
  assert_eq(values[2], 3.0)
  assert_eq(values[3], 4.0)
  assert_eq(values[4], 5.0)
  assert_eq(values[5], 6.0)
}

test "bmat_foreach_row_col_verification" {
  let data = [[10.0, 20.0], [30.0, 40.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let rows : Array[Int] = []
  let cols : Array[Int] = []
  let vals : Array[Double] = []
  @src.bmat_foreach(bmat, fn(row, col, val) {
    rows.push(row)
    cols.push(col)
    vals.push(val)
  })
  assert_eq(rows.length(), 4)
  assert_eq(rows[0], 0)
  assert_eq(rows[1], 0)
  assert_eq(rows[2], 1)
  assert_eq(rows[3], 1)
  assert_eq(cols[0], 0)
  assert_eq(cols[1], 1)
  assert_eq(cols[2], 0)
  assert_eq(cols[3], 1)
  assert_eq(vals[0], 10.0)
  assert_eq(vals[1], 20.0)
  assert_eq(vals[2], 30.0)
  assert_eq(vals[3], 40.0)
}

// ---------------------------------------------------------------------------
// BmatIterator
// ---------------------------------------------------------------------------

test "bmat_iterator_basic" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  let v0 = iter.next()
  assert_true(v0 is Some(_))
  assert_eq(v0.unwrap(), 1.0)
  let v1 = iter.next()
  assert_eq(v1.unwrap(), 2.0)
  let v2 = iter.next()
  assert_eq(v2.unwrap(), 3.0)
  let v3 = iter.next()
  assert_eq(v3.unwrap(), 4.0)
  let v4 = iter.next()
  assert_true(v4 is None)
}

test "bmat_iterator_has_next" {
  let data = [[1.0, 2.0, 3.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  assert_true(iter.has_next())
  let _ = iter.next()
  assert_true(iter.has_next())
  let _ = iter.next()
  assert_true(iter.has_next())
  let _ = iter.next()
  assert_false(iter.has_next())
}

test "bmat_iterator_exhausted_next" {
  let data = [[1.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  assert_true(iter.has_next())
  let v = iter.next()
  assert_eq(v.unwrap(), 1.0)
  assert_false(iter.has_next())
  let v2 = iter.next()
  assert_true(v2 is None)
}

test "bmat_iterator_cur_row_col" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  assert_eq(iter.cur_row(), 0)
  assert_eq(iter.cur_col(), 0)
  let _ = iter.next()
  assert_eq(iter.cur_row(), 0)
  assert_eq(iter.cur_col(), 1)
  let _ = iter.next()
  assert_eq(iter.cur_row(), 0)
  assert_eq(iter.cur_col(), 2)
  let _ = iter.next()
  assert_eq(iter.cur_row(), 1)
  assert_eq(iter.cur_col(), 0)
  let _ = iter.next()
  assert_eq(iter.cur_row(), 1)
  assert_eq(iter.cur_col(), 1)
  let _ = iter.next()
  assert_eq(iter.cur_row(), 1)
  assert_eq(iter.cur_col(), 2)
}

test "bmat_iterator_full_traversal" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  let values : Array[Double] = []
  while iter.has_next() {
    let v = iter.next()
    values.push(v.unwrap())
  }
  assert_eq(values.length(), 9)
  assert_eq(values[0], 1.0)
  assert_eq(values[1], 2.0)
  assert_eq(values[2], 3.0)
  assert_eq(values[3], 4.0)
  assert_eq(values[4], 5.0)
  assert_eq(values[5], 6.0)
  assert_eq(values[6], 7.0)
  assert_eq(values[7], 8.0)
  assert_eq(values[8], 9.0)
}

// ---------------------------------------------------------------------------
// Conversion
// ---------------------------------------------------------------------------

test "bmat_from_array" {
  let data = [[5.0, 10.0], [15.0, 20.0]]
  let bmat = @src.bmat_from_array(data)
  assert_eq(bmat.dim_rows(), 2)
  assert_eq(bmat.dim_cols(), 2)
  assert_eq(bmat.get(0, 0), 5.0)
  assert_eq(bmat.get(0, 1), 10.0)
  assert_eq(bmat.get(1, 0), 15.0)
  assert_eq(bmat.get(1, 1), 20.0)
}

test "bmat_to_array" {
  let data = [[2.0, 4.0, 6.0], [8.0, 10.0, 12.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let arr = @src.bmat_to_array(bmat)
  assert_eq(arr.length(), 2)
  assert_eq(arr[0].length(), 3)
  assert_eq(arr[0][0], 2.0)
  assert_eq(arr[0][1], 4.0)
  assert_eq(arr[0][2], 6.0)
  assert_eq(arr[1][0], 8.0)
  assert_eq(arr[1][1], 10.0)
  assert_eq(arr[1][2], 12.0)
}

test "bmat_roundtrip" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.bmat_from_array(data)
  let arr = @src.bmat_to_array(bmat)
  assert_eq(arr.length(), 3)
  for i in 0..<3 {
    for j in 0..<3 {
      assert_eq(arr[i][j], data[i][j])
    }
  }
}

// ---------------------------------------------------------------------------
// Subset, transpose, bind
// ---------------------------------------------------------------------------

test "bmat_subset" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let sub = @src.bmat_subset(bmat, 0, 2, 0, 2)
  assert_eq(sub.dim_rows(), 2)
  assert_eq(sub.dim_cols(), 2)
  assert_eq(sub.get(0, 0), 1.0)
  assert_eq(sub.get(0, 1), 2.0)
  assert_eq(sub.get(1, 0), 4.0)
  assert_eq(sub.get(1, 1), 5.0)
}

test "bmat_subset_bottom_right" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let sub = @src.bmat_subset(bmat, 1, 3, 1, 3)
  assert_eq(sub.dim_rows(), 2)
  assert_eq(sub.dim_cols(), 2)
  assert_eq(sub.get(0, 0), 5.0)
  assert_eq(sub.get(0, 1), 6.0)
  assert_eq(sub.get(1, 0), 8.0)
  assert_eq(sub.get(1, 1), 9.0)
}

test "bmat_transpose" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let t = @src.bmat_transpose(bmat)
  assert_eq(t.dim_rows(), 3)
  assert_eq(t.dim_cols(), 2)
  assert_eq(t.get(0, 0), 1.0)
  assert_eq(t.get(0, 1), 4.0)
  assert_eq(t.get(1, 0), 2.0)
  assert_eq(t.get(1, 1), 5.0)
  assert_eq(t.get(2, 0), 3.0)
  assert_eq(t.get(2, 1), 6.0)
}

test "bmat_transpose_square" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let t = @src.bmat_transpose(bmat)
  assert_eq(t.get(0, 0), 1.0)
  assert_eq(t.get(0, 1), 3.0)
  assert_eq(t.get(1, 0), 2.0)
  assert_eq(t.get(1, 1), 4.0)
}

test "bmat_bind_cols" {
  let left = @src.Bmat::new([[1.0, 2.0], [3.0, 4.0]], @src.BmatParam::new())
  let right = @src.Bmat::new([[5.0], [6.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_cols(left, right)
  assert_eq(result.dim_rows(), 2)
  assert_eq(result.dim_cols(), 3)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(0, 1), 2.0)
  assert_eq(result.get(0, 2), 5.0)
  assert_eq(result.get(1, 0), 3.0)
  assert_eq(result.get(1, 1), 4.0)
  assert_eq(result.get(1, 2), 6.0)
}

test "bmat_bind_cols_equal_width" {
  let left = @src.Bmat::new([[1.0, 2.0], [3.0, 4.0]], @src.BmatParam::new())
  let right = @src.Bmat::new([[5.0, 6.0], [7.0, 8.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_cols(left, right)
  assert_eq(result.dim_rows(), 2)
  assert_eq(result.dim_cols(), 4)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(0, 1), 2.0)
  assert_eq(result.get(0, 2), 5.0)
  assert_eq(result.get(0, 3), 6.0)
  assert_eq(result.get(1, 0), 3.0)
  assert_eq(result.get(1, 1), 4.0)
  assert_eq(result.get(1, 2), 7.0)
  assert_eq(result.get(1, 3), 8.0)
}

test "bmat_bind_rows" {
  let top = @src.Bmat::new([[1.0, 2.0], [3.0, 4.0]], @src.BmatParam::new())
  let bottom = @src.Bmat::new([[5.0, 6.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_rows(top, bottom)
  assert_eq(result.dim_rows(), 3)
  assert_eq(result.dim_cols(), 2)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(0, 1), 2.0)
  assert_eq(result.get(1, 0), 3.0)
  assert_eq(result.get(1, 1), 4.0)
  assert_eq(result.get(2, 0), 5.0)
  assert_eq(result.get(2, 1), 6.0)
}

test "bmat_bind_rows_equal_height" {
  let top = @src.Bmat::new([[1.0, 2.0], [3.0, 4.0]], @src.BmatParam::new())
  let bottom = @src.Bmat::new([[5.0, 6.0], [7.0, 8.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_rows(top, bottom)
  assert_eq(result.dim_rows(), 4)
  assert_eq(result.dim_cols(), 2)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(1, 0), 3.0)
  assert_eq(result.get(2, 0), 5.0)
  assert_eq(result.get(3, 0), 7.0)
  assert_eq(result.get(0, 1), 2.0)
  assert_eq(result.get(1, 1), 4.0)
  assert_eq(result.get(2, 1), 6.0)
  assert_eq(result.get(3, 1), 8.0)
}

// ---------------------------------------------------------------------------
// Element-wise operations
// ---------------------------------------------------------------------------

test "bmat_apply_elementwise_double" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let result = @src.bmat_apply_elementwise(bmat, fn(x) { x * 2.0 })
  assert_eq(result.dim_rows(), 2)
  assert_eq(result.dim_cols(), 3)
  assert_eq(result.get(0, 0), 2.0)
  assert_eq(result.get(0, 1), 4.0)
  assert_eq(result.get(0, 2), 6.0)
  assert_eq(result.get(1, 0), 8.0)
  assert_eq(result.get(1, 1), 10.0)
  assert_eq(result.get(1, 2), 12.0)
}

test "bmat_apply_elementwise_square" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let result = @src.bmat_apply_elementwise(bmat, fn(x) { x * x })
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(0, 1), 4.0)
  assert_eq(result.get(1, 0), 9.0)
  assert_eq(result.get(1, 1), 16.0)
}

test "bmat_apply_elementwise_negate" {
  let data = [[1.0, -2.0], [3.0, -4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let result = @src.bmat_apply_elementwise(bmat, fn(x) { -x })
  assert_eq(result.get(0, 0), -1.0)
  assert_eq(result.get(0, 1), 2.0)
  assert_eq(result.get(1, 0), -3.0)
  assert_eq(result.get(1, 1), 4.0)
}

test "bmat_apply_elementwise_no_change" {
  let data = [[5.0, 10.0], [15.0, 20.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let result = @src.bmat_apply_elementwise(bmat, fn(x) { x })
  assert_eq(result.get(0, 0), 5.0)
  assert_eq(result.get(0, 1), 10.0)
  assert_eq(result.get(1, 0), 15.0)
  assert_eq(result.get(1, 1), 20.0)
}

// ---------------------------------------------------------------------------
// to_string
// ---------------------------------------------------------------------------

test "bmat_to_string" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let s = bmat.to_string()
  assert_true(s.contains("2 x 2"))
  assert_true(s.contains("1"))
  assert_true(s.contains("2"))
  assert_true(s.contains("3"))
  assert_true(s.contains("4"))
}

test "bmat_to_string_single" {
  let data = [[7.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let s = bmat.to_string()
  assert_true(s.contains("1 x 1"))
  assert_true(s.contains("7"))
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "bmat_empty_matrix" {
  let data : Array[Array[Double]] = []
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 0)
  assert_eq(bmat.dim_cols(), 0)
}

test "bmat_empty_matrix_zeros" {
  let bmat = @src.Bmat::zeros(0, 0)
  assert_eq(bmat.dim_rows(), 0)
  assert_eq(bmat.dim_cols(), 0)
}

test "bmat_empty_matrix_iterator" {
  let data : Array[Array[Double]] = []
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  assert_false(iter.has_next())
  let v = iter.next()
  assert_true(v is None)
  assert_eq(iter.cur_row(), 0)
  assert_eq(iter.cur_col(), 0)
}

test "bmat_empty_matrix_roundtrip" {
  let data : Array[Array[Double]] = []
  let bmat = @src.bmat_from_array(data)
  assert_eq(bmat.dim_rows(), 0)
  assert_eq(bmat.dim_cols(), 0)
  let arr = @src.bmat_to_array(bmat)
  assert_eq(arr.length(), 0)
}

test "bmat_single_row" {
  let data = [[1.0, 2.0, 3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 1)
  assert_eq(bmat.dim_cols(), 4)
  assert_eq(bmat.get(0, 0), 1.0)
  assert_eq(bmat.get(0, 3), 4.0)
  let row = bmat.get_row(0)
  assert_eq(row.length(), 4)
  assert_eq(row[2], 3.0)
}

test "bmat_single_col" {
  let data = [[1.0], [2.0], [3.0], [4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 4)
  assert_eq(bmat.dim_cols(), 1)
  assert_eq(bmat.get(0, 0), 1.0)
  assert_eq(bmat.get(3, 0), 4.0)
  let col = bmat.get_col(0)
  assert_eq(col.length(), 4)
  assert_eq(col[2], 3.0)
}

test "bmat_single_element" {
  let data = [[42.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  assert_eq(bmat.dim_rows(), 1)
  assert_eq(bmat.dim_cols(), 1)
  assert_eq(bmat.get(0, 0), 42.0)
  bmat.set(0, 0, 99.0)
  assert_eq(bmat.get(0, 0), 99.0)
}

test "bmat_apply_col_blocks_single_block" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col_starts : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 10, fn(block) {
    col_starts.push(block.col_start())
  })
  assert_eq(col_starts.length(), 1)
  assert_eq(col_starts[0], 0)
}

test "bmat_apply_row_blocks_single_block" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_starts : Array[Int] = []
  @src.bmat_apply_row_blocks(bmat, 10, fn(block) {
    row_starts.push(block.row_start())
  })
  assert_eq(row_starts.length(), 1)
  assert_eq(row_starts[0], 0)
}

test "bmat_apply_col_blocks_exact_fit" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col_starts : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 3, fn(block) {
    col_starts.push(block.col_start())
  })
  assert_eq(col_starts.length(), 1)
  assert_eq(col_starts[0], 0)
}

test "bmat_apply_row_blocks_exact_fit" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_starts : Array[Int] = []
  @src.bmat_apply_row_blocks(bmat, 3, fn(block) {
    row_starts.push(block.row_start())
  })
  assert_eq(row_starts.length(), 1)
  assert_eq(row_starts[0], 0)
}

test "bmat_iterator_single_element" {
  let data = [[42.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let iter = @src.BmatIterator::new(bmat)
  assert_true(iter.has_next())
  assert_eq(iter.cur_row(), 0)
  assert_eq(iter.cur_col(), 0)
  let v = iter.next()
  assert_eq(v.unwrap(), 42.0)
  assert_false(iter.has_next())
}

test "bmat_subset_single_element" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let sub = @src.bmat_subset(bmat, 1, 2, 1, 2)
  assert_eq(sub.dim_rows(), 1)
  assert_eq(sub.dim_cols(), 1)
  assert_eq(sub.get(0, 0), 4.0)
}

test "bmat_transpose_empty" {
  let data : Array[Array[Double]] = []
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let t = @src.bmat_transpose(bmat)
  assert_eq(t.dim_rows(), 0)
  assert_eq(t.dim_cols(), 0)
}

test "bmat_foreach_empty" {
  let data : Array[Array[Double]] = []
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let values : Array[Double] = []
  @src.bmat_foreach(bmat, fn(row, col, val) {
    values.push(val)
  })
  assert_eq(values.length(), 0)
}

test "bmat_zeros_set_get" {
  let bmat = @src.Bmat::zeros(3, 3)
  assert_eq(bmat.get(1, 1), 0.0)
  bmat.set(1, 1, 5.5)
  assert_eq(bmat.get(1, 1), 5.5)
  assert_eq(bmat.get(0, 0), 0.0)
  assert_eq(bmat.get(2, 2), 0.0)
}

test "bmat_apply_elementwise_empty" {
  let data : Array[Array[Double]] = []
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let result = @src.bmat_apply_elementwise(bmat, fn(x) { x * 2.0 })
  assert_eq(result.dim_rows(), 0)
  assert_eq(result.dim_cols(), 0)
}

test "bmat_bind_cols_empty" {
  let left = @src.Bmat::zeros(3, 0)
  let right = @src.Bmat::new([[1.0], [2.0], [3.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_cols(left, right)
  assert_eq(result.dim_rows(), 3)
  assert_eq(result.dim_cols(), 1)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(2, 0), 3.0)
}

test "bmat_bind_rows_empty" {
  let top = @src.Bmat::zeros(0, 2)
  let bottom = @src.Bmat::new([[1.0, 2.0]], @src.BmatParam::new())
  let result = @src.bmat_bind_rows(top, bottom)
  assert_eq(result.dim_rows(), 1)
  assert_eq(result.dim_cols(), 2)
  assert_eq(result.get(0, 0), 1.0)
  assert_eq(result.get(0, 1), 2.0)
}

test "bmat_get_set_multiple" {
  let bmat = @src.Bmat::zeros(3, 3)
  bmat.set(0, 0, 1.0)
  bmat.set(0, 1, 2.0)
  bmat.set(0, 2, 3.0)
  bmat.set(1, 0, 4.0)
  bmat.set(1, 1, 5.0)
  bmat.set(1, 2, 6.0)
  bmat.set(2, 0, 7.0)
  bmat.set(2, 1, 8.0)
  bmat.set(2, 2, 9.0)
  assert_eq(bmat.get(0, 0), 1.0)
  assert_eq(bmat.get(0, 2), 3.0)
  assert_eq(bmat.get(1, 1), 5.0)
  assert_eq(bmat.get(2, 0), 7.0)
  assert_eq(bmat.get(2, 2), 9.0)
  let row1 = bmat.get_row(1)
  assert_eq(row1[0], 4.0)
  assert_eq(row1[1], 5.0)
  assert_eq(row1[2], 6.0)
  let col2 = bmat.get_col(2)
  assert_eq(col2[0], 3.0)
  assert_eq(col2[1], 6.0)
  assert_eq(col2[2], 9.0)
}

test "bmat_apply_col_blocks_last_block_partial" {
  let data = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let col_ends : Array[Int] = []
  @src.bmat_apply_col_blocks(bmat, 3, fn(block) {
    col_ends.push(block.col_end())
  })
  assert_eq(col_ends.length(), 3)
  assert_eq(col_ends[0], 3)
  assert_eq(col_ends[1], 6)
  assert_eq(col_ends[2], 7)
}

test "bmat_apply_row_blocks_last_block_partial" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 10.0]]
  let bmat = @src.Bmat::new(data, @src.BmatParam::new())
  let row_ends : Array[Int] = []
  @src.bmat_apply_row_blocks(bmat, 2, fn(block) {
    row_ends.push(block.row_end())
  })
  assert_eq(row_ends.length(), 3)
  assert_eq(row_ends[0], 2)
  assert_eq(row_ends[1], 4)
  assert_eq(row_ends[2], 5)
}