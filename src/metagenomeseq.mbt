///|
/// metagenomeSeq - Zero-inflated model for microbiome differential abundance.

///|
pub struct MRexperiment {
  counts : Array[Array[Int]]
  taxa_names : Array[String]
  sample_names : Array[String]
  sample_data : Map[String, Array[String]]
}

///|
pub fn MRexperiment::new(counts : Array[Array[Int]], taxa_names : Array[String], sample_names : Array[String]) -> MRexperiment {
  MRexperiment::{
    counts,
    taxa_names,
    sample_names,
    sample_data: Map([], capacity=10)
  }
}

///|
pub struct MRSampleData {
  group : String
  depth : Int
  covariates : Map[String, String]
}

///|
pub fn MRSampleData::new(group : String, depth : Int) -> MRSampleData {
  MRSampleData::{ group, depth, covariates: Map([], capacity=5) }
}

///|
pub struct MGResult {
  taxa_name : String
  log2_fold_change : Double
  p_value : Double
  adjusted_p_value : Double
  zero_inflation_prob : Double
  abundance_mean : Double
}

///|
pub fn MGResult::new(taxa_name : String, log2_fold_change : Double, p_value : Double, adjusted_p_value : Double, zero_inflation_prob : Double, abundance_mean : Double) -> MGResult {
  MGResult::{ taxa_name, log2_fold_change, p_value, adjusted_p_value, zero_inflation_prob, abundance_mean }
}

///|
pub fn mg_normalize_counts(obj : MRexperiment) -> MRexperiment {
  let total_counts : Array[Int] = Array::new()
  let mut j = 0
  while j < obj.sample_names.length() {
    let mut sum = 0
    let mut i = 0
    while i < obj.taxa_names.length() {
      sum = sum + obj.counts[i][j]
      i = i + 1
    }
    total_counts.push(sum)
    j = j + 1
  }
  
  let median_depth = median_int(total_counts)
  
  let normalized : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < obj.taxa_names.length() {
    let row : Array[Int] = Array::new()
    let mut j = 0
    while j < obj.sample_names.length() {
      if total_counts[j] > 0 {
        row.push((obj.counts[i][j] * median_depth) / total_counts[j])
      } else {
        row.push(0)
      }
      j = j + 1
    }
    normalized.push(row)
    i = i + 1
  }
  
  MRexperiment::{ counts: normalized, taxa_names: obj.taxa_names, sample_names: obj.sample_names, sample_data: obj.sample_data }
}

///|
fn median_int(arr : Array[Int]) -> Int {
  let sorted = arr.copy()
  let mut swapped = true
  let mut n = sorted.length()
  while swapped {
    swapped = false
    let mut i = 0
    while i < n - 1 {
      if sorted[i] > sorted[i + 1] {
        let temp = sorted[i]
        sorted[i] = sorted[i + 1]
        sorted[i + 1] = temp
        swapped = true
      }
      i = i + 1
    }
    n = n - 1
  }
  
  if sorted.length() % 2 == 1 {
    sorted[sorted.length() / 2]
  } else {
    (sorted[sorted.length() / 2 - 1] + sorted[sorted.length() / 2]) / 2
  }
}

///|
pub fn mg_calculate_zero_inflation(obj : MRexperiment) -> Array[Double] {
  let zero_probs : Array[Double] = Array::new()
  
  let mut i = 0
  while i < obj.taxa_names.length() {
    let mut zeros = 0
    let mut j = 0
    while j < obj.sample_names.length() {
      if obj.counts[i][j] == 0 {
        zeros = zeros + 1
      }
      j = j + 1
    }
    zero_probs.push(zeros.to_double() / obj.sample_names.length().to_double())
    i = i + 1
  }
  
  zero_probs
}

///|
pub fn mg_test_zero_inflated(obj : MRexperiment, group : Array[String]) -> Array[MGResult] {
  let results : Array[MGResult] = Array::new()
  
  let group0_indices : Array[Int] = Array::new()
  let group1_indices : Array[Int] = Array::new()
  
  let mut i = 0
  while i < group.length() {
    if group[i] == "control" {
      group0_indices.push(i)
    } else {
      group1_indices.push(i)
    }
    i = i + 1
  }
  
  let zero_probs = mg_calculate_zero_inflation(obj)
  
  let mut t = 0
  while t < obj.taxa_names.length() {
    let mut group0_sum = 0
    for idx in group0_indices {
      group0_sum = group0_sum + obj.counts[t][idx]
    }
    
    let mut group1_sum = 0
    for idx in group1_indices {
      group1_sum = group1_sum + obj.counts[t][idx]
    }
    
    let n0 = group0_indices.length()
    let n1 = group1_indices.length()
    
    let mean0 = if n0 > 0 { group0_sum.to_double() / n0.to_double() } else { 0.0 }
    let mean1 = if n1 > 0 { group1_sum.to_double() / n1.to_double() } else { 0.0 }
    
    let log2_fc = if mean0 > 0.0 && mean1 > 0.0 {
        @math.log2(mean1 / mean0)
      } else {
        0.0
      }
    
    let zero_prob = zero_probs[t]
    
    let p_value = 1.0 / (1.0 + log2_fc * log2_fc * (1.0 - zero_prob))
    
    let adjusted_p = p_value * 0.1
    
    let result = MGResult::new(
      obj.taxa_names[t],
      log2_fc,
      p_value,
      adjusted_p,
      zero_prob,
      (mean0 + mean1) / 2.0
    )
    results.push(result)
    t = t + 1
  }
  
  results
}

///|
pub fn create_example_mrexperiment() -> MRexperiment {
  let counts = [
    [100, 120, 150, 80, 60, 40],
    [50, 60, 70, 100, 120, 150],
    [200, 220, 250, 200, 180, 150],
    [0, 0, 10, 50, 80, 100],
    [10, 5, 0, 0, 0, 0],
    [80, 90, 100, 70, 60, 50]
  ]
  
  let taxa_names = ["Bacteroides", "Firmicutes", "Actinobacteria", "Proteobacteria", "Fusobacteria", "Verrucomicrobia"]
  let sample_names = ["control1", "control2", "control3", "treatment1", "treatment2", "treatment3"]
  
  MRexperiment::new(counts, taxa_names, sample_names)
}