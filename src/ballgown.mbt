///|
/// Ballgown - Transcript-level differential expression analysis.
/// Inspired by Bioconductor ballgown package.

///|
/// Transcript structure for ballgown.
pub struct BgTranscript {
  id : String
  gene_id : String
  gene_name : String
  chr : String
  strand : String
  num_exons : Int
  fpkm : Array[Double]
  cov : Array[Double]
}

///|
/// Exon structure for ballgown.
pub struct BgExon {
  id : String
  chr : String
  start : Int
  end : Int
  strand : String
  transcript_ids : Array[String]
}

///|
/// Intron structure for ballgown.
pub struct BgIntron {
  id : String
  chr : String
  start : Int
  end : Int
  strand : String
  transcript_ids : Array[String]
}

///|
/// Create a new transcript.
pub fn BgTranscript::new(
  id : String,
  gene_id : String,
  gene_name : String,
  chr : String,
  strand : String,
  num_exons : Int,
  fpkm : Array[Double],
  cov : Array[Double],
) -> BgTranscript {
  BgTranscript::{
    id,
    gene_id,
    gene_name,
    chr,
    strand,
    num_exons,
    fpkm,
    cov,
  }
}

///|
/// Create a new exon.
pub fn BgExon::new(
  id : String,
  chr : String,
  start : Int,
  end : Int,
  strand : String,
  transcript_ids : Array[String],
) -> BgExon {
  BgExon::{ id, chr, start, end, strand, transcript_ids }
}

///|
/// Create a new intron.
pub fn BgIntron::new(
  id : String,
  chr : String,
  start : Int,
  end : Int,
  strand : String,
  transcript_ids : Array[String],
) -> BgIntron {
  BgIntron::{ id, chr, start, end, strand, transcript_ids }
}

///|
/// Ballgown main structure.
pub struct Ballgown {
  sample_names : Array[String]
  transcripts : Array[BgTranscript]
  exons : Array[BgExon]
  introns : Array[BgIntron]
}

///|
/// Create a new Ballgown object.
pub fn Ballgown::new(
  sample_names : Array[String],
  transcripts : Array[BgTranscript],
  exons : Array[BgExon],
  introns : Array[BgIntron],
) -> Ballgown {
  Ballgown::{ sample_names, transcripts, exons, introns }
}

///|
/// Get FPKM matrix for all transcripts.
pub fn Ballgown::get_fpkm_matrix(self : Ballgown) -> Array[Array[Double]] {
  let mat : Array[Array[Double]] = []
  for t in self.transcripts {
    mat.push(t.fpkm)
  }
  mat
}

///|
/// Get coverage matrix for all transcripts.
pub fn Ballgown::get_cov_matrix(self : Ballgown) -> Array[Array[Double]] {
  let mat : Array[Array[Double]] = []
  for t in self.transcripts {
    mat.push(t.cov)
  }
  mat
}

///|
/// Get transcript by ID.
pub fn Ballgown::get_transcript(self : Ballgown, id : String) -> BgTranscript? {
  for t in self.transcripts {
    if t.id == id {
      return Some(t)
    }
  }
  None
}

///|
/// Get transcripts by gene ID.
pub fn Ballgown::get_gene_transcripts(self : Ballgown, gene_id : String) -> Array[BgTranscript] {
  let result : Array[BgTranscript] = []
  for t in self.transcripts {
    if t.gene_id == gene_id {
      result.push(t)
    }
  }
  result
}

///|
/// Get unique gene IDs.
pub fn Ballgown::get_gene_ids(self : Ballgown) -> Array[String] {
  let result : Array[String] = []
  for t in self.transcripts {
    let found = {
      let mut f = false
      for r in result {
        if r == t.gene_id {
          f = true
          break
        }
      }
      f
    }
    if !found {
      result.push(t.gene_id)
    }
  }
  result
}

///|
/// Get unique gene names.
pub fn Ballgown::get_gene_names(self : Ballgown) -> Array[String] {
  let result : Array[String] = []
  for t in self.transcripts {
    let found = {
      let mut f = false
      for r in result {
        if r == t.gene_name {
          f = true
          break
        }
      }
      f
    }
    if !found && t.gene_name.length() > 0 {
      result.push(t.gene_name)
    }
  }
  result
}

///|
/// Differential expression result.
pub struct DeResult {
  feature_id : String
  gene_id : String
  gene_name : String
  log2_fc : Double
  p_value : Double
  q_value : Double
  mean_fpkm : Double
}

///|
/// Create a new DE result.
pub fn DeResult::new(
  feature_id : String,
  gene_id : String,
  gene_name : String,
  log2_fc : Double,
  p_value : Double,
  q_value : Double,
  mean_fpkm : Double,
) -> DeResult {
  DeResult::{ feature_id, gene_id, gene_name, log2_fc, p_value, q_value, mean_fpkm }
}

///|
/// Welch's t-test p-value approximation (simplified).
fn t_test_pvalue(x : Array[Double], y : Array[Double]) -> Double {
  let n1 = x.length()
  let n2 = y.length()
  if n1 < 2 || n2 < 2 {
    return 1.0
  }

  let mean1 = {
    let mut sum = 0.0
    for v in x { sum = sum + v }
    sum / n1.to_double()
  }
  let mean2 = {
    let mut sum = 0.0
    for v in y { sum = sum + v }
    sum / n2.to_double()
  }

  let var1 = {
    let mut sum = 0.0
    for v in x {
      let d = v - mean1
      sum = sum + d * d
    }
    sum / (n1 - 1).to_double()
  }
  let var2 = {
    let mut sum = 0.0
    for v in y {
      let d = v - mean2
      sum = sum + d * d
    }
    sum / (n2 - 1).to_double()
  }

  let se = (var1 / n1.to_double() + var2 / n2.to_double()).sqrt()
  if se == 0.0 {
    return 1.0
  }
  let t_stat = (mean1 - mean2) / se

  let t_abs = t_stat.abs()
  let approx = 2.0 * (1.0 - normal_cdf(t_abs))
  if approx < 0.0 { 0.0 } else if approx > 1.0 { 1.0 } else { approx }
}

///|
/// Normal CDF approximation (Abramowitz and Stegun).
fn normal_cdf(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.2316419 * x.abs())
  let d = 0.3989423 * @math.exp(-x * x / 2.0)
  let prob = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x > 0.0 { 1.0 - prob } else { prob }
}

///|
/// Benjamini-Hochberg FDR correction.
fn fdr_correct(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return []
  }

  let indexed : Array[(Int, Double)] = []
  for i = 0; i < n; i = i + 1 {
    indexed.push((i, pvalues[i]))
  }

  let sorted = {
    let arr = indexed
    for i = 0; i < arr.length() - 1; i = i + 1 {
      for j = i + 1; j < arr.length(); j = j + 1 {
        if arr[i].1 > arr[j].1 {
          let tmp = arr[i]
          arr[i] = arr[j]
          arr[j] = tmp
        }
      }
    }
    arr
  }

  let qvalues : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    qvalues.push(0.0)
  }

  let mut prev_q = 1.0
  let n_d = n.to_double()
  for i = n - 1; i >= 0; i = i - 1 {
    let (idx, p) = sorted[i]
    let q = p * n_d / (i + 1).to_double()
    let adj_q = if q < prev_q { q } else { prev_q }
    let clamped = if adj_q > 1.0 { 1.0 } else { adj_q }
    qvalues[idx] = clamped
    prev_q = clamped
  }

  qvalues
}

///|
/// Perform statistical test for differential expression.
/// feature can be "transcript", "gene", "exon", or "intron".
pub fn stattest(
  bg : Ballgown,
  group : Array[Int],
  feature? : String = "transcript",
) -> Array[DeResult] {
  let results : Array[DeResult] = []
  let pvalues : Array[Double] = []
  let n_samples = bg.sample_names.length()

  let group1_idx : Array[Int] = []
  let group2_idx : Array[Int] = []
  for i = 0; i < n_samples; i = i + 1 {
    if group[i] == 0 {
      group1_idx.push(i)
    } else {
      group2_idx.push(i)
    }
  }
  let n1 = group1_idx.length()
  let n2 = group2_idx.length()

  match feature {
    "transcript" => {
      for t in bg.transcripts {
        let group1 : Array[Double] = []
        let group2 : Array[Double] = []
        for idx in group1_idx { group1.push(t.fpkm[idx] + 0.001) }
        for idx in group2_idx { group2.push(t.fpkm[idx] + 0.001) }

        let mean1 = {
          let mut sum = 0.0
          for v in group1 { sum = sum + v }
          sum / n1.to_double()
        }
        let mean2 = {
          let mut sum = 0.0
          for v in group2 { sum = sum + v }
          sum / n2.to_double()
        }

        let log2_fc = if mean2 > 0.0 && mean1 > 0.0 {
          @math.log2(mean2) - @math.log2(mean1)
        } else {
          0.0
        }

        let p_value = t_test_pvalue(group1, group2)
        pvalues.push(p_value)

        let mean_fpkm = (mean1 + mean2) / 2.0
        results.push(DeResult::new(
          t.id,
          t.gene_id,
          t.gene_name,
          log2_fc,
          p_value,
          0.0,
          mean_fpkm,
        ))
      }
    }
    "gene" => {
      let gene_ids = bg.get_gene_ids()
      for gene_id in gene_ids {
        let gene_transcripts = bg.get_gene_transcripts(gene_id)
        if gene_transcripts.length() > 0 {
          let gene_name = gene_transcripts[0].gene_name
          let gene_fpkm : Array[Double] = []
          for i = 0; i < n_samples; i = i + 1 {
            let mut total = 0.0
            for t in gene_transcripts {
              total = total + t.fpkm[i]
            }
            gene_fpkm.push(total)
          }

          let group1 : Array[Double] = []
          let group2 : Array[Double] = []
          for idx in group1_idx { group1.push(gene_fpkm[idx] + 0.001) }
          for idx in group2_idx { group2.push(gene_fpkm[idx] + 0.001) }

          let mean1 = {
            let mut sum = 0.0
            for v in group1 { sum = sum + v }
            sum / n1.to_double()
          }
          let mean2 = {
            let mut sum = 0.0
            for v in group2 { sum = sum + v }
            sum / n2.to_double()
          }

          let log2_fc = if mean2 > 0.0 && mean1 > 0.0 {
            @math.log2(mean2) - @math.log2(mean1)
          } else {
            0.0
          }

          let p_value = t_test_pvalue(group1, group2)
          pvalues.push(p_value)

          let mean_fpkm = (mean1 + mean2) / 2.0
          results.push(DeResult::new(
            gene_id,
            gene_id,
            gene_name,
            log2_fc,
            p_value,
            0.0,
            mean_fpkm,
          ))
        }
      }
    }
    _ => {
      let mut i = 0
      while i < bg.transcripts.length() {
        let t = bg.transcripts[i]
        let group1 : Array[Double] = []
        let group2 : Array[Double] = []
        for idx in group1_idx { group1.push(t.fpkm[idx] + 0.001) }
        for idx in group2_idx { group2.push(t.fpkm[idx] + 0.001) }

        let mean1 = {
          let mut sum = 0.0
          for v in group1 { sum = sum + v }
          sum / n1.to_double()
        }
        let mean2 = {
          let mut sum = 0.0
          for v in group2 { sum = sum + v }
          sum / n2.to_double()
        }

        let log2_fc = if mean2 > 0.0 && mean1 > 0.0 {
          @math.log2(mean2) - @math.log2(mean1)
        } else {
          0.0
        }

        let p_value = t_test_pvalue(group1, group2)
        pvalues.push(p_value)

        let mean_fpkm = (mean1 + mean2) / 2.0
        results.push(DeResult::new(
          t.id,
          t.gene_id,
          t.gene_name,
          log2_fc,
          p_value,
          0.0,
          mean_fpkm,
        ))
        i = i + 1
      }
    }
  }

  let qvalues = fdr_correct(pvalues)
  let final_results : Array[DeResult] = []
  for i = 0; i < results.length(); i = i + 1 {
    let r = results[i]
    final_results.push(DeResult::new(
      r.feature_id,
      r.gene_id,
      r.gene_name,
      r.log2_fc,
      r.p_value,
      qvalues[i],
      r.mean_fpkm,
    ))
  }

  final_results
}

///|
/// Get gene-level expression (sum of transcript FPKM).
pub fn gene_expression(bg : Ballgown, gene_id : String) -> Array[Double] {
  let transcripts = bg.get_gene_transcripts(gene_id)
  if transcripts.length() == 0 {
    return []
  }
  let n_samples = bg.sample_names.length()
  let result : Array[Double] = []
  for i = 0; i < n_samples; i = i + 1 {
    let mut total = 0.0
    for t in transcripts {
      total = total + t.fpkm[i]
    }
    result.push(total)
  }
  result
}

///|
/// Filter transcripts by expression threshold.
pub fn filter_low_expression(
  bg : Ballgown,
  min_fpkm? : Double = 1.0,
) -> Ballgown {
  let filtered_transcripts : Array[BgTranscript] = []
  for t in bg.transcripts {
    let mut keep = false
    for f in t.fpkm {
      if f > min_fpkm {
        keep = true
        break
      }
    }
    if keep {
      filtered_transcripts.push(t)
    }
  }

  Ballgown::new(
    bg.sample_names,
    filtered_transcripts,
    bg.exons,
    bg.introns,
  )
}

///|
/// Create sample ballgown data for testing.
pub fn sample_ballgown() -> Ballgown {
  let sample_names = ["s1", "s2", "s3", "s4", "s5", "s6"]
  let group = [0, 0, 0, 1, 1, 1]

  let transcripts : Array[BgTranscript] = []
  let exons : Array[BgExon] = []
  let introns : Array[BgIntron] = []

  let gene_names = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
  let chrs = ["chr1", "chr1", "chr2", "chr2", "chr3"]
  let strands = ["+", "-", "+", "-", "+"]
  let tx_counts = [2, 2, 3, 1, 2]
  let tx_exon_counts = [
    [3, 2],
    [4, 2],
    [2, 3, 2],
    [5],
    [2, 3],
  ]

  let mut exon_counter = 1
  let mut intron_counter = 1
  let mut tx_start = 1000

  for gene_idx = 0; gene_idx < 5; gene_idx = gene_idx + 1 {
    let gene_name = gene_names[gene_idx]
    let chr = chrs[gene_idx]
    let strand = strands[gene_idx]
    let gene_id = "ENSG" + (gene_idx + 1).to_string()
    let mut current_pos = tx_start
    let n_tx = tx_counts[gene_idx]

    for tx_idx = 0; tx_idx < n_tx; tx_idx = tx_idx + 1 {
      let num_exons = tx_exon_counts[gene_idx][tx_idx]
      let tx_id = "ENST" + (gene_idx * 10 + tx_idx + 1).to_string()
      let fpkm : Array[Double] = []
      let cov : Array[Double] = []

      let base_expr = 5.0 + gene_idx.to_double() * 3.0
      let de_effect = if gene_idx % 2 == 0 { 2.0 } else { 0.5 }

      for i = 0; i < 6; i = i + 1 {
        let noise = (gene_idx * 7 + i * 13 + tx_idx * 3).to_double() % 10.0 / 10.0 - 0.5
        let expr = if group[i] == 1 {
          base_expr * de_effect + noise
        } else {
          base_expr + noise
        }
        fpkm.push(if expr > 0.0 { expr } else { 0.1 })
        cov.push(if expr > 0.0 { expr * 10.0 } else { 1.0 })
      }

      let transcript_ids_for_features = [tx_id]
      let tx_exon_starts : Array[Int] = []
      let tx_exon_ends : Array[Int] = []

      for e = 0; e < num_exons; e = e + 1 {
        let exon_start = current_pos + e * 200
        let exon_end = exon_start + 100
        tx_exon_starts.push(exon_start)
        tx_exon_ends.push(exon_end)

        let exon_id = "exon_" + exon_counter.to_string()
        exon_counter = exon_counter + 1
        exons.push(BgExon::new(
          exon_id,
          chr,
          exon_start,
          exon_end,
          strand,
          transcript_ids_for_features,
        ))

        if e > 0 {
          let intron_start = tx_exon_ends[e - 1] + 1
          let intron_end = exon_start - 1
          let intron_id = "intron_" + intron_counter.to_string()
          intron_counter = intron_counter + 1
          introns.push(BgIntron::new(
            intron_id,
            chr,
            intron_start,
            intron_end,
            strand,
            transcript_ids_for_features,
          ))
        }
      }

      transcripts.push(BgTranscript::new(
        tx_id,
        gene_id,
        gene_name,
        chr,
        strand,
        num_exons,
        fpkm,
        cov,
      ))

      current_pos = current_pos + num_exons * 200 + 500
    }

    tx_start = tx_start + 10000
  }

  Ballgown::new(sample_names, transcripts, exons, introns)
}
