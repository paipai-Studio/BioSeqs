///|
test "ruvseq_create" {
  let counts : Array[Array[Double]] = Array::new()
  let gene_counts : Array[Double] = Array::new()
  gene_counts.push(10.0)
  gene_counts.push(20.0)
  counts.push(gene_counts)

  let genes : Array[String] = Array::new()
  genes.push("gene1")

  let samples : Array[String] = Array::new()
  samples.push("sample1")
  samples.push("sample2")

  let batch : Array[String] = Array::new()
  batch.push("batch1")
  batch.push("batch2")

  let data = @src.RUVSeqData::new(counts, genes, samples, batch)
  assert_eq(data.genes.length(), 1)
}

///|
test "ruvseq_normalize" {
  let counts : Array[Array[Double]] = Array::new()
  let gene_counts : Array[Double] = Array::new()
  gene_counts.push(10.0)
  gene_counts.push(20.0)
  counts.push(gene_counts)

  let genes : Array[String] = Array::new()
  genes.push("gene1")

  let samples : Array[String] = Array::new()
  samples.push("sample1")
  samples.push("sample2")

  let batch : Array[String] = Array::new()
  batch.push("batch1")
  batch.push("batch2")

  let data = @src.RUVSeqData::new(counts, genes, samples, batch)
  let normalized = data.normalize()
  assert_true(normalized.normalized)
}

///|
test "ruvseq_log2_transform" {
  let counts : Array[Array[Double]] = Array::new()
  let gene_counts : Array[Double] = Array::new()
  gene_counts.push(8.0)
  gene_counts.push(16.0)
  counts.push(gene_counts)

  let genes : Array[String] = Array::new()
  genes.push("gene1")

  let samples : Array[String] = Array::new()
  samples.push("sample1")
  samples.push("sample2")

  let batch : Array[String] = Array::new()
  batch.push("batch1")
  batch.push("batch2")

  let data = @src.RUVSeqData::new(counts, genes, samples, batch)
  let log_data = data.log2_transform()
  assert_eq(log_data.counts.length(), 1)
}

///|
test "ruvseq_remove_batch_effect" {
  let counts : Array[Array[Double]] = Array::new()
  let gene_counts : Array[Double] = Array::new()
  gene_counts.push(10.0)
  gene_counts.push(20.0)
  counts.push(gene_counts)

  let genes : Array[String] = Array::new()
  genes.push("gene1")

  let samples : Array[String] = Array::new()
  samples.push("sample1")
  samples.push("sample2")

  let batch : Array[String] = Array::new()
  batch.push("batch1")
  batch.push("batch2")

  let data = @src.RUVSeqData::new(counts, genes, samples, batch)

  let factors_arr : Array[Array[Double]] = Array::new()
  let factor : Array[Double] = Array::new()
  factor.push(1.0)
  factor.push(1.0)
  factors_arr.push(factor)

  let factors = @src.RUVFactors::new(factors_arr, 1)
  let corrected = @src.remove_batch_effect(data, factors)
  assert_eq(corrected.counts.length(), 1)
}
