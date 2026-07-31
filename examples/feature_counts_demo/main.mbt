///|
/// featureCounts module example - Read counting over genomic features.
///
/// Run: moon run examples/feature_counts_demo
fn main {
  println("=== featureCounts Module Demo ===")
  println("Read Counting over Genomic Features")
  println("")

  // Create feature annotation
  println("=== Feature Annotation ===")
  let features : Array[@bio.FeatureAnnotation] = Array::new()
  features.push(@bio.FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+",
    gene_id="BRCA1", feature_id="exon_1",
  ))
  features.push(@bio.FeatureAnnotation::new(
    chr="chr1", start=3000, end_=5000, strand="+",
    gene_id="BRCA1", feature_id="exon_2",
  ))
  features.push(@bio.FeatureAnnotation::new(
    chr="chr17", start=43044295, end_=43045000, strand="+",
    gene_id="TP53", feature_id="exon_1",
  ))
  println("Features created: " + features.length().to_string())
  let mut i = 0
  while i < features.length() {
    let f = features[i]
    println("  " + f.gene_id + ": " + f.chr + ":" + f.start.to_string() + "-" + f.end_.to_string() + " (len=" + f.length().to_string() + ")")
    i = i + 1
  }
  println("")

  // Create reads
  println("=== Read Alignments ===")
  let reads : Array[@bio.ReadAlignment] = Array::new()
  reads.push(@bio.ReadAlignment::new(
    read_id="read_1", chr="chr1", start=1500, end_=1600, strand="+", mapq=60,
  ))
  reads.push(@bio.ReadAlignment::new(
    read_id="read_2", chr="chr1", start=3500, end_=3600, strand="+", mapq=60,
  ))
  reads.push(@bio.ReadAlignment::new(
    read_id="read_3", chr="chr1", start=1500, end_=1700, strand="+", mapq=60,
  ))
  reads.push(@bio.ReadAlignment::new(
    read_id="read_4", chr="chr17", start=43044500, end_=43044600, strand="+", mapq=60,
  ))
  reads.push(@bio.ReadAlignment::new(
    read_id="read_5", chr="chr1", start=1500, end_=1600, strand="+", mapq=20,
  ))
  reads.push(@bio.ReadAlignment::new(
    read_id="read_6", chr="chr5", start=1000, end_=1100, strand="+", mapq=60,
  ))
  println("Reads created: " + reads.length().to_string())
  println("")

  // Run counting
  println("=== Feature Counting ===")
  let config = @bio.FeatureCountsConfig::new()
  println("Config: min_overlap=" + config.min_overlap.to_string() + ", min_mapq=" + config.min_mapq.to_string())
  let result = @bio.feature_counts_count(features, reads, sample_names=["tumor"], config=config)
  println(result.summary())
  println("")

  // Show counts
  println("=== Per-Feature Counts ===")
  let mut j = 0
  while j < features.length() {
    let fid = features[j].feature_id
    let count = result.get_count(feature_id=fid, sample_idx=0)
    let cpm = result.cpm(feature_id=fid, sample_idx=0)
    println("  " + fid + " (" + features[j].gene_id + "): count=" + count.to_string() + ", CPM=" + cpm.to_string())
    j = j + 1
  }
  println("")

  // Configuration options
  println("=== Configuration Options ===")
  let cfg_stranded = @bio.FeatureCountsConfig::new().set_strand_mode(mode=@bio.strand_stranded())
  let result_stranded = @bio.feature_counts_count(features, reads, sample_names=["tumor"], config=cfg_stranded)
  println("Stranded counting assigned: " + result_stranded.assigned_reads[0].to_string())

  let cfg_reversed = @bio.FeatureCountsConfig::new().set_strand_mode(mode=@bio.strand_reversed())
  let result_reversed = @bio.feature_counts_count(features, reads, sample_names=["tumor"], config=cfg_reversed)
  println("Reversed stranded assigned: " + result_reversed.assigned_reads[0].to_string())

  let cfg_mapq = @bio.FeatureCountsConfig::new().set_min_mapq(val=30)
  let result_mapq = @bio.feature_counts_count(features, reads, sample_names=["tumor"], config=cfg_mapq)
  println("MAPQ>=30 assigned: " + result_mapq.assigned_reads[0].to_string())
  println("")

  // Paired-end reads
  println("=== Paired-End Read Counting ===")
  let pe_reads : Array[@bio.ReadAlignment] = Array::new()
  pe_reads.push(@bio.ReadAlignment::new_paired(
    read_id="pe_read_1", chr="chr1", start=1500, end_=1600,
    strand="+", mate_start=1800,
  ))
  let pe_result = @bio.feature_counts_count(features, pe_reads, sample_names=["tumor"])
  println("PE reads assigned: " + pe_result.assigned_reads[0].to_string())
  println("")

  // Multi-sample counting
  println("=== Multi-Sample Counting ===")
  let multi_reads : Array[@bio.ReadAlignment] = Array::new()
  multi_reads.push(@bio.ReadAlignment::new(read_id="r1", chr="chr1", start=1500, end_=1600))
  multi_reads.push(@bio.ReadAlignment::new(read_id="r2", chr="chr1", start=3500, end_=3600))
  let multi_result = @bio.feature_counts_count(
    features, multi_reads,
    sample_names=["tumor", "normal"],
  )
  println("Samples: " + multi_result.n_samples().to_string())
  println("Sample 0 assigned: " + multi_result.assigned_reads[0].to_string())
  println("Sample 1 assigned: " + multi_result.assigned_reads[1].to_string())
  println("")

  // Sample data
  println("=== Sample Data ===")
  let (sf, sr, ss) = @bio.feature_counts_sample_data()
  let sample_result = @bio.feature_counts_count(sf, sr, sample_names=ss)
  println(sample_result.to_string())
}
