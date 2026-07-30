///|
/// NOISeq module example - Noise-robust differential expression for RNA-seq.
///
/// Run: moon run examples/noiseq_demo
fn main {
  println("=== NOISeq Module Demo ===")
  println("Noise-Robust Differential Expression Analysis")
  println("")

  // Generate sample data
  println("=== Sample Data ===")
  let (ctrl, trt) = @bio.noiseq_sample_data()
  println("Control samples: \{ctrl.length()}")
  println("Treatment samples: \{trt.length()}")
  let mut i = 0
  while i < ctrl.length() {
    println("  Control \{ctrl[i].sample_id}: lib_size=\{ctrl[i].library_size()}, expressed=\{ctrl[i].n_expressed()}")
    i = i + 1
  }
  println("")

  // Normalization
  println("=== Normalization (TMM) ===")
  let norm_ctrl = @bio.noiseq_normalize(ctrl, norm_method=@bio.noiseq_norm_tmm())
  let norm_trt = @bio.noiseq_normalize(trt, norm_method=@bio.noiseq_norm_tmm())
  println("Normalized control samples: \{norm_ctrl.length()}")
  println("Normalized treatment samples: \{norm_trt.length()}")
  println("")

  // Run NOISeq analysis
  println("=== NOISeqBio Analysis ===")
  let results = @bio.noiseq_run(
    ctrl, trt,
    ds_method=@bio.noiseq_method_bio(),
    norm=@bio.noiseq_norm_tmm(),
    prob_threshold=0.3,
  )
  println(results.summary())
  println("")

  // Top genes
  println("=== Top 5 Genes ===")
  let top = results.get_top_genes(5)
  let mut j = 0
  while j < top.length() {
    let g = top[j]
    println("  \{g.gene_id}: log2FC=\{g.log2fc}, prob=\{g.prob}, rank=\{g.ranking}")
    j = j + 1
  }
  println("")

  // Significant genes
  println("=== Significant Genes (prob >= 0.3) ===")
  let sig = results.get_significant()
  println("Total: \{sig.length()}")
  let mut k = 0
  while k < sig.length() {
    let g = sig[k]
    let direction = if g.log2fc > 0.0 { "UP" } else { "DOWN" }
    println("  \{g.gene_id}: \{direction} log2FC=\{g.log2fc}, prob=\{g.prob}")
    k = k + 1
  }
  println("")

  // Up and down regulated
  println("=== Direction of Change ===")
  let up = results.get_up_regulated()
  let down = results.get_down_regulated()
  println("Up-regulated: \{up.length()}")
  println("Down-regulated: \{down.length()}")
  println("")

  // QC diagnostics
  println("=== QC Diagnostics ===")
  let all : Array[@bio.NOISeqSample] = Array::new()
  let mut ci = 0
  while ci < ctrl.length() {
    all.push(ctrl[ci])
    ci = ci + 1
  }
  let mut ti = 0
  while ti < trt.length() {
    all.push(trt[ti])
    ti = ti + 1
  }
  let qc = @bio.noiseq_qc(all)
  println("GC content: \{qc.gc_content}%")
  println("Bias detected: \{qc.bias_detected}")
  println("")

  // Normalization comparison
  println("=== Normalization Methods ===")
  let norm_none = @bio.noiseq_normalize(ctrl, norm_method=@bio.noiseq_norm_none())
  let norm_tmm = @bio.noiseq_normalize(ctrl, norm_method=@bio.noiseq_norm_tmm())
  let norm_uq = @bio.noiseq_normalize(ctrl, norm_method=@bio.noiseq_norm_upper())
  println("None: Gene1 = \{norm_none[0].get_count("Gene1")}")
  println("TMM: Gene1 = \{norm_tmm[0].get_count("Gene1")}")
  println("UpperQuartile: Gene1 = \{norm_uq[0].get_count("Gene1")}")
}
