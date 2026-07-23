///|
/// CNVkit Demo - Copy Number Variation Detection.
/// Inspired by Bioconductor CNVkit package.

///|
fn main {
  println("=== CNVkit Copy Number Variation Demo ===")
  println("Inspired by Bioconductor CNVkit package")
  println("")
  
  println("1. Creating example CNV dataset...")
  let dataset = @src.create_example_cnv_dataset()
  println("   Total probes: \{dataset.count_probes().to_string()}")
  println("   Chromosomes: \{dataset.chromosomes.length().to_string()}")
  println("")
  
  println("2. CBS segmentation analysis...")
  let segments : Array[@src.CNVSegment] = Array::new()
  let probes = dataset.probes
  let n = probes.length()
  
  // Perform CBS on each chromosome
  let chromosomes = dataset.chromosomes
  let mut i = 0
  while i < chromosomes.length() {
    let chrom = chromosomes[i]
    let chrom_probes : Array[@src.CNVProbe] = Array::new()
    
    let mut j = 0
    while j < n {
      if probes[j].chromosome == chrom {
        chrom_probes.push(probes[j])
      }
      j = j + 1
    }
    
    if chrom_probes.length() >= 3 {
      let result = @src.cbs_segment(chrom_probes, 0.05)
      let segs = result.segments
      
      let mut k = 0
      while k < segs.length() {
        segments.push(segs[k])
        k = k + 1
      }
    }
    i = i + 1
  }
  
  println("   Segments found: \{segments.length().to_string()}")
  println("")
  
  println("3. Segment details:")
  let mut seg_i = 0
  while seg_i < segments.length() {
    let seg = segments[seg_i]
    println("   Segment \{seg_i.to_string()}: \{seg.chromosome}:\{seg.start.to_string()}-\{seg.end.to_string()} (\{seg.probes.to_string()} probes, mean log2 = \{seg.mean_log2.to_string()}, state = \{seg.state})")
    seg_i = seg_i + 1
  }
  println("")
  
  println("4. Copy number state calling...")
  let calls = @src.call_copy_numbers(segments, 2)
  let mut del_count = 0
  let mut amp_count = 0
  let mut neu_count = 0
  let mut c = 0
  while c < calls.length() {
    if calls[c].state == "deletion" {
      del_count = del_count + 1
    } else if calls[c].state == "amplification" {
      amp_count = amp_count + 1
    } else {
      neu_count = neu_count + 1
    }
    c = c + 1
  }
  println("   Deletions: \{del_count.to_string()}")
  println("   Amplifications: \{amp_count.to_string()}")
  println("   Neutral: \{neu_count.to_string()}")
  println("")
  
  println("5. Log2 ratio smoothing...")
  let smoothed = @src.smooth_log2_ratios(probes, 5)
  println("   Smoothed probes: \{smoothed.length().to_string()}")
  println("")
  
  println("6. Breakpoint detection...")
  let breakpoints = @src.detect_breakpoints(segments, 0.3)
  println("   Breakpoints detected: \{breakpoints.length().to_string()}")
  println("")
  
  println("7. Dataset summary:")
  let summary_str = @src.cbs_summarize_dataset(dataset)
  println("   \{summary_str}")
  
  println("=== CNVkit demo completed ===")
}
