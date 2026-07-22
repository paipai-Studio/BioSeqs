///|
/// BiocParallel - Parallel computation framework.

///|
pub struct BPPARAM {
  workers : Int
  progressbar : Bool
  timeout : Int
  log_file : String
}

///|
pub fn BPPARAM::new(workers : Int, progressbar : Bool, timeout : Int) -> BPPARAM {
  BPPARAM::{ workers, progressbar, timeout, log_file: "" }
}

///|
pub fn BPPARAM::new_with_log(workers : Int, progressbar : Bool, timeout : Int, log_file : String) -> BPPARAM {
  BPPARAM::{ workers, progressbar, timeout, log_file }
}

///|
pub struct BPResult {
  success : Bool
  result : Double
  error : String
  worker_id : Int
}

///|
pub fn BPResult::success(result : Double, worker_id : Int) -> BPResult {
  BPResult::{ success: true, result, error: "", worker_id }
}

///|
pub fn BPResult::error(error : String, worker_id : Int) -> BPResult {
  BPResult::{ success: false, result: 0.0, error, worker_id }
}

///|
pub struct Task {
  id : String
  data : Array[Double]
  priority : Int
}

///|
pub fn Task::new(id : String, data : Array[Double], priority : Int) -> Task {
  Task::{ id, data, priority }
}

///|
pub struct BPJob {
  name : String
  tasks : Array[Task]
  params : BPPARAM
}

///|
pub fn BPJob::new(name : String, params : BPPARAM) -> BPJob {
  BPJob::{ name, tasks: Array::new(), params }
}

///|
pub fn bp_add_task(job : BPJob, task : Task) -> BPJob {
  job.tasks.push(task)
  job
}

///|
pub fn bp_sum(data : Array[Array[Double]], params : BPPARAM) -> Double {
  let mut total = 0.0
  
  let chunk_size = if data.length() % params.workers == 0 {
    data.length() / params.workers
  } else {
    data.length() / params.workers + 1
  }
  
  let mut start = 0
  while start < data.length() {
    let end = if start + chunk_size < data.length() {
      start + chunk_size
    } else {
      data.length()
    }
    
    let mut i = start
    while i < end {
      let mut j = 0
      while j < data[i].length() {
        total = total + data[i][j]
        j = j + 1
      }
      i = i + 1
    }
    
    start = end
  }
  
  total
}

///|
pub fn bp_mean(data : Array[Array[Double]], params : BPPARAM) -> Double {
  let total = bp_sum(data, params)
  let mut count = 0.0
  
  let mut i = 0
  while i < data.length() {
    count = count + data[i].length().to_double()
    i = i + 1
  }
  
  if count > 0.0 { total / count } else { 0.0 }
}

///|
pub fn bp_mean_simple(data : Array[Double], workers : Int) -> Double {
  let chunks = bp_parallelize(data, workers)
  let params = BPPARAM::new(workers, false, 120)
  bp_mean(chunks, params)
}

///|
pub fn bp_colsum(data : Array[Array[Double]], params : BPPARAM) -> Array[Double] {
  if data.length() == 0 {
    return Array::new()
  }
  
  let col_count = data[0].length()
  let results : Array[Double] = Array::new()
  
  let mut j = 0
  while j < col_count {
    let mut sum = 0.0
    let mut i = 0
    while i < data.length() {
      if j < data[i].length() {
        sum = sum + data[i][j]
      }
      i = i + 1
    }
    results.push(sum)
    j = j + 1
  }
  
  results
}

///|
pub fn bp_rowsum(data : Array[Array[Double]], params : BPPARAM) -> Array[Double] {
  let results : Array[Double] = Array::new()
  
  let mut i = 0
  while i < data.length() {
    let mut sum = 0.0
    let mut j = 0
    while j < data[i].length() {
      sum = sum + data[i][j]
      j = j + 1
    }
    results.push(sum)
    i = i + 1
  }
  
  results
}

///|
pub fn bp_parallelize(data : Array[Double], workers : Int) -> Array[Array[Double]] {
  let chunks : Array[Array[Double]] = Array::new()
  
  let chunk_size = if data.length() % workers == 0 {
    data.length() / workers
  } else {
    data.length() / workers + 1
  }
  
  let mut start = 0
  while start < data.length() {
    let chunk : Array[Double] = Array::new()
    let end = if start + chunk_size < data.length() {
      start + chunk_size
    } else {
      data.length()
    }
    
    let mut i = start
    while i < end {
      chunk.push(data[i])
      i = i + 1
    }
    
    chunks.push(chunk)
    start = end
  }
  
  chunks
}

///|
pub fn bp_run(job : BPJob) -> Array[BPResult] {
  let results : Array[BPResult] = Array::new()
  
  let mut i = 0
  while i < job.tasks.length() {
    let task = job.tasks[i]
    let mut sum = 0.0
    let mut j = 0
    while j < task.data.length() {
      sum = sum + task.data[j]
      j = j + 1
    }
    results.push(BPResult::success(sum, i % job.params.workers))
    i = i + 1
  }
  
  results
}

///|
pub fn bp_map_sum(data : Array[Array[Double]], params : BPPARAM) -> Array[BPResult] {
  let results : Array[BPResult] = Array::new()
  
  let chunk_size = if data.length() % params.workers == 0 {
    data.length() / params.workers
  } else {
    data.length() / params.workers + 1
  }
  
  let mut worker_id = 0
  let mut start = 0
  
  while start < data.length() {
    let end = if start + chunk_size < data.length() {
      start + chunk_size
    } else {
      data.length()
    }
    
    let mut i = start
    while i < end {
      let mut sum = 0.0
      let mut j = 0
      while j < data[i].length() {
        sum = sum + data[i][j]
        j = j + 1
      }
      results.push(BPResult::success(sum, worker_id))
      i = i + 1
    }
    
    worker_id = worker_id + 1
    start = end
  }
  
  results
}

///|
pub fn bp_map_mean(data : Array[Array[Double]], params : BPPARAM) -> Array[BPResult] {
  let results : Array[BPResult] = Array::new()
  
  let chunk_size = if data.length() % params.workers == 0 {
    data.length() / params.workers
  } else {
    data.length() / params.workers + 1
  }
  
  let mut worker_id = 0
  let mut start = 0
  
  while start < data.length() {
    let end = if start + chunk_size < data.length() {
      start + chunk_size
    } else {
      data.length()
    }
    
    let mut i = start
    while i < end {
      let mut sum = 0.0
      let mut j = 0
      while j < data[i].length() {
        sum = sum + data[i][j]
        j = j + 1
      }
      let mean = if data[i].length() > 0 { sum / data[i].length().to_double() } else { 0.0 }
      results.push(BPResult::success(mean, worker_id))
      i = i + 1
    }
    
    worker_id = worker_id + 1
    start = end
  }
  
  results
}

///|
pub fn create_example_bpparam() -> BPPARAM {
  BPPARAM::new_with_log(4, true, 300, "bp.log")
}

///|
pub fn create_example_bpjob() -> BPJob {
  let params = BPPARAM::new(2, false, 120)
  let job = BPJob::new("example_job", params)
  
  let task1_data : Array[Double] = Array::new()
  task1_data.push(1.0)
  task1_data.push(2.0)
  task1_data.push(3.0)
  task1_data.push(4.0)
  task1_data.push(5.0)
  job.tasks.push(Task::new("task1", task1_data, 1))
  
  let task2_data : Array[Double] = Array::new()
  task2_data.push(6.0)
  task2_data.push(7.0)
  task2_data.push(8.0)
  task2_data.push(9.0)
  task2_data.push(10.0)
  job.tasks.push(Task::new("task2", task2_data, 1))
  
  let task3_data : Array[Double] = Array::new()
  task3_data.push(11.0)
  task3_data.push(12.0)
  task3_data.push(13.0)
  task3_data.push(14.0)
  task3_data.push(15.0)
  job.tasks.push(Task::new("task3", task3_data, 2))
  
  let task4_data : Array[Double] = Array::new()
  task4_data.push(16.0)
  task4_data.push(17.0)
  task4_data.push(18.0)
  task4_data.push(19.0)
  task4_data.push(20.0)
  job.tasks.push(Task::new("task4", task4_data, 2))
  
  job
}