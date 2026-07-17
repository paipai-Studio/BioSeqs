///|
pub struct RUVSeqData {
  counts : Array[Array[Double]]
  genes : Array[String]
  samples : Array[String]
  batch : Array[String]
  normalized : Bool
}

///|
pub fn RUVSeqData::new(
  counts : Array[Array[Double]],
  genes : Array[String],
  samples : Array[String],
  batch : Array[String],
) -> RUVSeqData {
  RUVSeqData::{ counts, genes, samples, batch, normalized: false }
}

///|
pub fn RUVSeqData::normalize(self : RUVSeqData) -> RUVSeqData {
  let n_genes = self.genes.length()
  let n_samples = self.samples.length()

  let col_sums : Array[Double] = Array::new()
  let mut i = 0
  while i < n_samples {
    let mut sum = 0.0
    let mut j = 0
    while j < n_genes {
      sum = sum + self.counts[j][i]
      j = j + 1
    }
    col_sums.push(sum)
    i = i + 1
  }

  let mut total_sum = 0.0
  let mut k = 0
  while k < col_sums.length() {
    total_sum = total_sum + col_sums[k]
    k = k + 1
  }
  let mean_sum = total_sum / Float::from_int(col_sums.length()).to_double()

  let new_counts : Array[Array[Double]] = Array::new()
  i = 0
  while i < n_genes {
    let gene_counts : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let factor = mean_sum / col_sums[j]
      gene_counts.push(self.counts[i][j] * factor)
      j = j + 1
    }
    new_counts.push(gene_counts)
    i = i + 1
  }

  RUVSeqData::{
    counts: new_counts,
    genes: self.genes.copy(),
    samples: self.samples.copy(),
    batch: self.batch.copy(),
    normalized: true,
  }
}

///|
pub fn RUVSeqData::log2_transform(self : RUVSeqData) -> RUVSeqData {
  let n_genes = self.genes.length()
  let n_samples = self.samples.length()

  let new_counts : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_genes {
    let gene_counts : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      gene_counts.push(
        if self.counts[i][j] > 0.0 {
          @math.log2(self.counts[i][j])
        } else {
          0.0
        },
      )
      j = j + 1
    }
    new_counts.push(gene_counts)
    i = i + 1
  }

  RUVSeqData::{
    counts: new_counts,
    genes: self.genes.copy(),
    samples: self.samples.copy(),
    batch: self.batch.copy(),
    normalized: self.normalized,
  }
}

///|
pub fn RUVSeqData::get_gene_index(
  self : RUVSeqData,
  gene_name : String,
) -> Int? {
  let mut i = 0
  while i < self.genes.length() {
    if self.genes[i] == gene_name {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
pub fn RUVSeqData::get_sample_index(
  self : RUVSeqData,
  sample_name : String,
) -> Int? {
  let mut i = 0
  while i < self.samples.length() {
    if self.samples[i] == sample_name {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
pub struct RUVFactors {
  factors : Array[Array[Double]]
  k : Int
}

///|
pub fn RUVFactors::new(factors : Array[Array[Double]], k : Int) -> RUVFactors {
  RUVFactors::{ factors, k }
}

///|
pub fn ruvg(
  data : RUVSeqData,
  control_genes : Array[String],
  k : Int,
) -> RUVFactors {
  let n_samples = data.samples.length()
  let n_controls = control_genes.length()

  let control_counts : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_controls {
    let idx = data.get_gene_index(control_genes[i])
    if idx is Some(_) {
      control_counts.push(data.counts[idx.unwrap()].copy())
    }
    i = i + 1
  }

  let factors : Array[Array[Double]] = Array::new()
  let mut j = 0
  while j < k {
    let factor : Array[Double] = Array::new()
    let mut m = 0
    while m < n_samples {
      factor.push(0.0)
      m = m + 1
    }
    factors.push(factor)
    j = j + 1
  }

  let batch_means : Array[(String, Double)] = Array::new()
  let seen_batches : Array[String] = Array::new()
  i = 0
  while i < n_samples {
    let b = data.batch[i]
    let mut found = false
    let mut j_idx = 0
    while j_idx < seen_batches.length() {
      if seen_batches[j_idx] == b {
        found = true
        break
      }
      j_idx = j_idx + 1
    }
    if !found {
      let new_seen = seen_batches.copy()
      new_seen.push(b)
      let mut sum = 0.0
      let mut count = 0
      let mut m = 0
      while m < n_samples {
        if data.batch[m] == b {
          let mut n = 0
          while n < control_counts.length() {
            sum = sum + control_counts[n][m]
            n = n + 1
          }
          count = count + 1
        }
        m = m + 1
      }
      let bm = batch_means.copy()
      bm.push(
        (b, sum / Float::from_int(control_counts.length() * count).to_double()),
      )
    }
    i = i + 1
  }

  let final_factor : Array[Double] = Array::new()
  j = 0
  while j < n_samples {
    let b = data.batch[j]
    let mut mean = 0.0
    let mut m = 0
    while m < batch_means.length() {
      if batch_means[m].0 == b {
        mean = batch_means[m].1
        break
      }
      m = m + 1
    }
    final_factor.push(mean)
    j = j + 1
  }

  let final_factors : Array[Array[Double]] = Array::new()
  final_factors.push(final_factor)
  RUVFactors::new(final_factors, k)
}

///|
pub fn ruvs(data : RUVSeqData, k : Int) -> RUVFactors {
  let n_samples = data.samples.length()
  let n_genes = data.genes.length()

  let factors : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < k {
    let factor : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      factor.push(0.0)
      j = j + 1
    }
    factors.push(factor)
    i = i + 1
  }

  let batches : Array[String] = Array::new()
  i = 0
  while i < n_samples {
    let b = data.batch[i]
    let mut found = false
    let mut j_idx = 0
    while j_idx < batches.length() {
      if batches[j_idx] == b {
        found = true
        break
      }
      j_idx = j_idx + 1
    }
    if !found {
      let new_batches = batches.copy()
      new_batches.push(b)
    }
    i = i + 1
  }

  let final_factor : Array[Double] = Array::new()
  let mut j = 0
  while j < n_samples {
    let b = data.batch[j]

    let mut total_mean = 0.0
    let mut m = 0
    while m < n_genes {
      let mut sum = 0.0
      let mut count = 0
      let mut n = 0
      while n < n_samples {
        if data.batch[n] == b {
          sum = sum + data.counts[m][n]
          count = count + 1
        }
        n = n + 1
      }
      total_mean = total_mean + sum / Float::from_int(count).to_double()
      m = m + 1
    }
    total_mean = total_mean / Float::from_int(n_genes).to_double()

    final_factor.push(total_mean)
    j = j + 1
  }

  let final_factors : Array[Array[Double]] = Array::new()
  final_factors.push(final_factor)
  RUVFactors::new(final_factors, k)
}

///|
pub fn remove_batch_effect(
  data : RUVSeqData,
  factors : RUVFactors,
) -> RUVSeqData {
  let n_genes = data.genes.length()
  let n_samples = data.samples.length()
  let k = factors.k

  let new_counts : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < n_genes {
    let gene_counts : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let mut adjustment = 0.0
      let mut m = 0
      while m < k {
        adjustment = adjustment + factors.factors[m][j]
        m = m + 1
      }
      gene_counts.push(data.counts[i][j] - adjustment)
      j = j + 1
    }
    new_counts.push(gene_counts)
    i = i + 1
  }

  RUVSeqData::new(
    new_counts,
    data.genes.copy(),
    data.samples.copy(),
    data.batch.copy(),
  )
}

///|
pub fn create_demo_data() -> RUVSeqData {
  let genes = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "Control1", "Control2", "Control3",
  ]
  let samples = [
    "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6",
  ]
  let batch = ["Batch1", "Batch1", "Batch1", "Batch2", "Batch2", "Batch2"]

  let counts = [
    [100.0, 105.0, 98.0, 200.0, 210.0, 195.0],
    [50.0, 52.0, 48.0, 100.0, 105.0, 95.0],
    [75.0, 78.0, 72.0, 80.0, 82.0, 78.0],
    [200.0, 205.0, 198.0, 150.0, 155.0, 145.0],
    [150.0, 152.0, 148.0, 160.0, 165.0, 155.0],
    [1000.0, 1010.0, 990.0, 1005.0, 1015.0, 995.0],
    [2000.0, 2020.0, 1980.0, 2010.0, 2030.0, 1990.0],
    [1500.0, 1515.0, 1485.0, 1505.0, 1525.0, 1495.0],
  ]

  RUVSeqData::new(counts, genes, samples, batch)
}
