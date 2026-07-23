///|
/// Tests for BiocParallel module.

test "BPPARAM creation" {
  let params = @src.create_example_bpparam()
  assert_eq(params.workers, 4)
  assert_eq(params.progressbar, true)
}

test "BPJob creation" {
  let job = @src.create_example_bpjob()
  assert_eq(job.name, "example_job")
  assert_eq(job.tasks.length(), 4)
}

test "bp_sum" {
  let data : Array[Array[Double]] = Array::new()
  let chunk1 : Array[Double] = Array::new()
  chunk1.push(1.0)
  chunk1.push(2.0)
  chunk1.push(3.0)
  data.push(chunk1)
  
  let chunk2 : Array[Double] = Array::new()
  chunk2.push(4.0)
  chunk2.push(5.0)
  chunk2.push(6.0)
  data.push(chunk2)
  
  let params = @src.BPPARAM::new(2, false, 120)
  let result = @src.bp_sum(data, params)
  assert_eq(result, 21.0)
}

test "bp_mean" {
  let data : Array[Double] = Array::new()
  data.push(1.0)
  data.push(2.0)
  data.push(3.0)
  data.push(4.0)
  data.push(5.0)
  
  let result = @src.bp_mean_simple(data, 2)
  assert_eq(result, 3.0)
}