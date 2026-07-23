///|
/// BiocParallel demo.

pub fn main() {
  let params = create_example_bpparam()
  println("BPPARAM created with")
  println("- workers: " + params.workers.to_string())
  println("- progressbar: " + params.progressbar.to_string())
  println("- timeout: " + params.timeout.to_string())
  
  let job = create_example_bpjob()
  println("\nBPJob: " + job.name)
  println("- tasks: " + job.tasks.length().to_string())
  
  let data : Array[Array[Double]] = Array::new()
  let mut chunk1 : Array[Double] = Array::new()
  chunk1.push(1.0)
  chunk1.push(2.0)
  chunk1.push(3.0)
  chunk1.push(4.0)
  chunk1.push(5.0)
  data.push(chunk1)
  
  let mut chunk2 : Array[Double] = Array::new()
  chunk2.push(6.0)
  chunk2.push(7.0)
  chunk2.push(8.0)
  chunk2.push(9.0)
  chunk2.push(10.0)
  data.push(chunk2)
  
  let sum_result = bp_sum(data, params)
  println("\nParallel sum: " + sum_result.to_string())
  
  let mean_data : Array[Double] = Array::new()
  mean_data.push(1.0)
  mean_data.push(2.0)
  mean_data.push(3.0)
  mean_data.push(4.0)
  mean_data.push(5.0)
  
  let mean_result = bp_mean(mean_data, 2)
  println("Parallel mean: " + mean_result.to_string())
}