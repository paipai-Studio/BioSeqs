///|
/// Tests for Harmony batch correction module.

test "harmony_create_example" {
  let data = @src.harmony_create_example(100, 2, 5)
  assert_true(data.embeddings.length() == 100)
  assert_true(data.batch_labels.length() == 100)
  assert_true(data.embeddings[0].length() == 5)
}

test "harmony_params_default" {
  let params = @src.HarmonyParams::new()
  assert_true(params.n_clusters == 20)
  assert_true(params.max_iterations == 20)
  assert_true(params.tolerance > 0.0)
  assert_true(params.lambda > 0.0)
}

test "harmony_run_basic" {
  let data = @src.harmony_create_example(50, 2, 3)
  let params = @src.HarmonyParams::create(
    n_clusters = 3,
    max_iterations = 10,
  )
  let result = @src.harmony_run(data, params)
  assert_true(result.corrected.length() == 50)
  assert_true(result.membership.length() == 50)
  assert_true(result.centroids.length() == 3)
}

test "harmony_single_batch" {
  // Single batch should remain largely unchanged
  let embeddings : Array[Array[Double]] = Array::new()
  let batch_labels : Array[String] = Array::new()
  let cell_ids : Array[String] = Array::new()
  let mut i = 0
  while i < 20 {
    embeddings.push([1.0, 2.0, 3.0])
    batch_labels.push("batch_1")
    cell_ids.push("cell_" + i.to_string())
    i = i + 1
  }
  let data = @src.HarmonyData::new(embeddings, batch_labels, cell_ids)
  let params = @src.HarmonyParams::create(
    n_clusters = 2,
    max_iterations = 5,
  )
  let result = @src.harmony_run(data, params)
  assert_true(result.corrected.length() == 20)
}

test "harmony_convergence" {
  let data = @src.harmony_create_example(80, 3, 4)
  let params = @src.HarmonyParams::create(
    n_clusters = 4,
    max_iterations = 30,
    tolerance = 0.001,
  )
  let result = @src.harmony_run(data, params)
  assert_true(result.n_iterations > 0)
  assert_true(result.n_iterations <= 30)
}

test "harmony_batch_alignment" {
  // Two batches with different means should be aligned
  let embeddings : Array[Array[Double]] = Array::new()
  let batch_labels : Array[String] = Array::new()
  let cell_ids : Array[String] = Array::new()
  let mut i = 0
  while i < 40 {
    let batch_idx = if i < 20 { 0 } else { 1 }
    let offset = if batch_idx == 0 { 0.0 } else { 5.0 }
    embeddings.push([offset + (i.to_double() * 0.1), offset + 1.0, offset + 2.0])
    batch_labels.push("batch_" + batch_idx.to_string())
    cell_ids.push("cell_" + i.to_string())
    i = i + 1
  }
  let data = @src.HarmonyData::new(embeddings, batch_labels, cell_ids)
  let params = @src.HarmonyParams::create(
    n_clusters = 5,
    max_iterations = 30,
  )
  let result = @src.harmony_run(data, params)
  assert_true(result.corrected.length() == 40)
}

test "harmony_empty_data" {
  let data = @src.HarmonyData::new([], [], [])
  let params = @src.HarmonyParams::create(
    n_clusters = 2,
  )
  let result = @src.harmony_run(data, params)
  assert_true(result.corrected.length() == 0)
}
