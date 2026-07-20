///|
/// Batchelor tests - Single-cell batch correction tests.

import "mooncake/test"

///|
test "batchelor_create_example_data" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 1.0)
  
  assert_eq(batches.length(), 2)
  assert_eq(batch_names.length(), 2)
  assert_eq(batches[0].length(), 50)
  assert_eq(batches[1].length(), 50)
  assert_eq(batches[0][0].length(), 100)
  assert_eq(batches[1][0].length(), 100)
  
  assert_eq(batch_names[0], "Batch_0")
  assert_eq(batch_names[1], "Batch_1")
}

///|
test "rescale_batches" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 1.0)
  
  let result = @bio.rescale_batches(batches, batch_names)
  
  assert_eq(result.corrected.length(), 100)
  assert_eq(result.batch_indices.length(), 100)
  assert_eq(result.corrected[0].length(), 100)
  
  let mut batch0_count = 0
  let mut batch1_count = 0
  for idx in result.batch_indices {
    if idx == 0 {
      batch0_count = batch0_count + 1
    } else {
      batch1_count = batch1_count + 1
    }
  }
  
  assert_eq(batch0_count, 50)
  assert_eq(batch1_count, 50)
}

///|
test "find_mutual_nearest_neighbors" {
  let (batches, _) = @bio.batchelor_create_example_data(2, 20, 50, 0.5)
  
  let mnn_pairs = @bio.find_mutual_nearest_neighbors(batches[0], batches[1], 3)
  
  assert_true(mnn_pairs.length() >= 0)
  
  for pair in mnn_pairs {
    assert_true(pair.0 >= 0 && pair.0 < 20)
    assert_true(pair.1 >= 0 && pair.1 < 20)
  }
}

///|
test "fast_mnn" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 1.0)
  
  let result = @bio.fast_mnn(batches, batch_names, 5, 20)
  
  assert_eq(result.corrected.length(), 100)
  assert_eq(result.batch_indices.length(), 100)
  assert_eq(result.corrected[0].length(), 20)
  assert_eq(result.var_explained.length(), 20)
  
  let mut sum_var = 0.0
  for var in result.var_explained {
    assert_true(var >= 0.0 && var <= 1.0)
    sum_var = sum_var + var
  }
  assert_true(sum_var > 0.5)
}

///|
test "compute_batch_mixing_score" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 1.0)
  
  let result = @bio.fast_mnn(batches, batch_names, 5, 20)
  let score = @bio.compute_batch_mixing_score(result.corrected, result.batch_indices)
  
  assert_true(score >= 0.0 && score <= 1.0)
  
  let raw_data : Array[Array[Double]] = Array::new()
  for cell in batches[0] {
    raw_data.push(cell)
  }
  for cell in batches[1] {
    raw_data.push(cell)
  }
  
  let raw_batch_indices : Array[Int] = Array::new()
  let mut idx = 0
  while idx < 50 {
    raw_batch_indices.push(0)
    idx = idx + 1
  }
  idx = 0
  while idx < 50 {
    raw_batch_indices.push(1)
    idx = idx + 1
  }
  
  let raw_score = @bio.compute_batch_mixing_score(raw_data, raw_batch_indices)
  
  assert_true(score > raw_score * 0.8)
}

///|
test "rescale_batches_correction_effect" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 2.0)
  
  let raw_data : Array[Array[Double]] = Array::new()
  for cell in batches[0] {
    raw_data.push(cell)
  }
  for cell in batches[1] {
    raw_data.push(cell)
  }
  
  let raw_batch_indices : Array[Int] = Array::new()
  let mut idx = 0
  while idx < 50 {
    raw_batch_indices.push(0)
    idx = idx + 1
  }
  idx = 0
  while idx < 50 {
    raw_batch_indices.push(1)
    idx = idx + 1
  }
  
  let raw_score = @bio.compute_batch_mixing_score(raw_data, raw_batch_indices)
  
  let result = @bio.rescale_batches(batches, batch_names)
  let corrected_score = @bio.compute_batch_mixing_score(result.corrected, result.batch_indices)
  
  assert_true(corrected_score > raw_score)
}

///|
test "fast_mnn_correction_effect" {
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 100, 2.0)
  
  let raw_data : Array[Array[Double]] = Array::new()
  for cell in batches[0] {
    raw_data.push(cell)
  }
  for cell in batches[1] {
    raw_data.push(cell)
  }
  
  let raw_batch_indices : Array[Int] = Array::new()
  let mut idx = 0
  while idx < 50 {
    raw_batch_indices.push(0)
    idx = idx + 1
  }
  idx = 0
  while idx < 50 {
    raw_batch_indices.push(1)
    idx = idx + 1
  }
  
  let raw_score = @bio.compute_batch_mixing_score(raw_data, raw_batch_indices)
  
  let result = @bio.fast_mnn(batches, batch_names, 5, 20)
  let corrected_score = @bio.compute_batch_mixing_score(result.corrected, result.batch_indices)
  
  assert_true(corrected_score > raw_score)
}
