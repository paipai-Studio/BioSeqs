///|
pub struct PeakAnnotation {
  peak_id : String
  chr : String
  start : Int
  end : Int
  annotation : String
  gene_name : String
  distance_to_tss : Int
  feature_type : String
  strand : String
}

///|
pub fn PeakAnnotation::new(
  peak_id : String,
  chr : String,
  start : Int,
  end : Int,
  annotation : String,
  gene_name : String,
  distance_to_tss : Int,
  feature_type : String,
  strand : String,
) -> PeakAnnotation {
  PeakAnnotation::{
    peak_id,
    chr,
    start,
    end,
    annotation,
    gene_name,
    distance_to_tss,
    feature_type,
    strand,
  }
}

///|
pub struct AnnotationResult {
  annotations : Array[PeakAnnotation]
  summary : Map[String, Int]
}

///|
pub struct Peak {
  chr : String
  start : Int
  end : Int
  peak_id : String
}

///|
pub fn Peak::new(chr : String, start : Int, end : Int, peak_id : String) -> Peak {
  Peak::{ chr, start, end, peak_id }
}

///|
pub fn annotate_peaks(
  peaks : Array[Peak],
  txdb : TxDb,
) -> AnnotationResult {
  let annotations : Array[PeakAnnotation] = Array::new()
  let summary = Map([], capacity=16)

  let mut i = 0
  while i < peaks.length() {
    let peak = peaks[i]
    let annotation = annotate_single_peak(peak, txdb)
    annotations.push(annotation)

    let feature_type = annotation.feature_type
    let current = summary.get(feature_type)
    summary[feature_type] = if current.is_some() { current.unwrap() + 1 } else { 1 }

    i = i + 1
  }

  AnnotationResult::{ annotations, summary }
}

///|
fn annotate_single_peak(peak : Peak, txdb : TxDb) -> PeakAnnotation {
  let peak_center = (peak.start + peak.end) / 2

  let (nearest_gene, distance, strand) = find_nearest_gene(peak.chr, peak_center, txdb)

  let (annotation, feature_type) = classify_genomic_region(distance)

  PeakAnnotation::new(
    peak.peak_id,
    peak.chr,
    peak.start,
    peak.end,
    annotation,
    nearest_gene,
    distance,
    feature_type,
    strand,
  )
}

///|
fn find_nearest_gene(chr : String, pos : Int, txdb : TxDb) -> (String, Int, String) {
  let mut min_distance = 100000000
  let mut nearest_gene = ""
  let mut gene_strand = "."

  let mut i = 0
  while i < txdb.genes.length() {
    let gene = txdb.genes[i]
    if gene.seqname != chr {
      i = i + 1
      continue
    }

    let gene_start = gene.start
    let gene_end = gene.end

    let distance = if pos < gene_start {
      gene_start - pos
    } else if pos > gene_end {
      pos - gene_end
    } else {
      0
    }

    if distance < min_distance {
      min_distance = distance
      nearest_gene = gene.gene_id
      gene_strand = match gene.strand {
        Strand::Plus => "+"
        Strand::Minus => "-"
        _ => "."
      }
    }
    i = i + 1
  }

  if nearest_gene == "" {
    ("Intergenic", 100000000, ".")
  } else {
    (nearest_gene, min_distance, gene_strand)
  }
}

///|
fn classify_genomic_region(distance : Int) -> (String, String) {
  if distance == 0 {
    ("Promoter (<=0bp)", "Promoter")
  } else if distance <= 2000 {
    ("Promoter (0-2kb)", "Promoter")
  } else if distance <= 5000 {
    ("Promoter (2-5kb)", "Promoter")
  } else if distance <= 10000 {
    ("Distal Intergenic (5-10kb)", "Distal")
  } else if distance <= 50000 {
    ("Distal Intergenic (10-50kb)", "Distal")
  } else {
    ("Intergenic (>50kb)", "Intergenic")
  }
}

///|
pub fn plot_annotation_bar(summary : Map[String, Int]) -> String {
  let mut result = "Annotation Distribution Bar Plot\n"
  result = result + "=============================\n\n"

  let keys = summary.keys().collect()
  let mut max_count = 0
  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    let count = summary[key]
    if count > max_count {
      max_count = count
    }
    i = i + 1
  }

  if max_count == 0 {
    return result + "No annotations found\n"
  }

  let sorted_keys : Array[String] = Array::make(keys.length(), "")
  let mut j = 0
  while j < keys.length() {
    sorted_keys[j] = keys[j]
    j = j + 1
  }

  let mut k = 0
  while k < sorted_keys.length() - 1 {
    let mut l = k + 1
    while l < sorted_keys.length() {
      if summary[sorted_keys[k]] < summary[sorted_keys[l]] {
        let temp = sorted_keys[k]
        sorted_keys[k] = sorted_keys[l]
        sorted_keys[l] = temp
      }
      l = l + 1
    }
    k = k + 1
  }

  let mut m = 0
  while m < sorted_keys.length() {
    let key = sorted_keys[m]
    let count = summary[key]
    let bar_length = (count.to_double() / max_count.to_double() * 40.0).to_int()
    let mut bar = ""
    let mut n = 0
    while n < bar_length {
      bar = bar + "="
      n = n + 1
    }
    result = result + key + ": " + count.to_string() + " [" + bar + "]\n"
    m = m + 1
  }

  result
}

///|
pub fn plot_dist_to_tss(annotations : Array[PeakAnnotation]) -> String {
  let mut result = "Distance to TSS Distribution\n"
  result = result + "============================\n\n"

  let bins = [
    ("0-100bp", 0, 100),
    ("100-500bp", 100, 500),
    ("500-1000bp", 500, 1000),
    ("1-2kb", 1000, 2000),
    ("2-5kb", 2000, 5000),
    ("5-10kb", 5000, 10000),
    (">10kb", 10000, 100000000),
  ]

  let mut max_count = 0
  let counts : Array[Int] = Array::make(bins.length(), 0)

  let mut i = 0
  while i < annotations.length() {
    let dist = annotations[i].distance_to_tss
    let mut j = 0
    while j < bins.length() {
      let (_, low, high) = bins[j]
      if dist >= low && dist < high {
        counts[j] = counts[j] + 1
        if counts[j] > max_count {
          max_count = counts[j]
        }
        break
      }
      j = j + 1
    }
    i = i + 1
  }

  let mut k = 0
  while k < bins.length() {
    let (label, _, _) = bins[k]
    let count = counts[k]
    let bar_length = if max_count > 0 {
      (count.to_double() / max_count.to_double() * 30.0).to_int()
    } else {
      0
    }
    let mut bar = ""
    let mut l = 0
    while l < bar_length {
      bar = bar + "*"
      l = l + 1
    }
    result = result + label + ": " + count.to_string() + " " + bar + "\n"
    k = k + 1
  }

  result
}

///|
pub fn get_annotation_statistics(result : AnnotationResult) -> Map[String, Double] {
  let stats = Map([], capacity=16)
  let total = result.annotations.length().to_double()

  for key in result.summary.keys() {
    let value = result.summary.get(key)
    if value.is_some() {
      let percentage = value.unwrap().to_double() / total * 100.0
      stats[key] = percentage
    }
  }

  stats
}

///|
pub fn get_promoter_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let promoters : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Promoter" {
      promoters.push(annotation)
    }
    i = i + 1
  }

  promoters
}

///|
pub fn get_distal_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let distal : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Distal" {
      distal.push(annotation)
    }
    i = i + 1
  }

  distal
}

///|
pub fn get_intergenic_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let intergenic : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Intergenic" {
      intergenic.push(annotation)
    }
    i = i + 1
  }

  intergenic
}

///|
pub fn annotation_summary(result : AnnotationResult) -> String {
  let mut result_str = "Peak Annotation Summary\n"
  result_str = result_str + "=======================\n"
  result_str = result_str + "Total peaks: " + result.annotations.length().to_string() + "\n\n"

  let stats = get_annotation_statistics(result)
  for key in stats.keys() {
    let value = stats.get(key)
    if value.is_some() {
      result_str = result_str + key + ": " + value.unwrap().to_string() + "%\n"
    }
  }

  let promoters = get_promoter_annotations(result)
  let distal = get_distal_annotations(result)
  let intergenic = get_intergenic_annotations(result)

  result_str = result_str + "\nPromoter peaks: " + promoters.length().to_string() + "\n"
  result_str = result_str + "Distal peaks: " + distal.length().to_string() + "\n"
  result_str = result_str + "Intergenic peaks: " + intergenic.length().to_string() + "\n"

  result_str
}

///|
pub fn cs_create_example_peaks() -> Array[Peak] {
  let peaks : Array[Peak] = Array::new()

  let mut i = 0
  while i < 50 {
    let chr = if i < 15 { "chr1" } else if i < 30 { "chr2" } else if i < 40 { "chr3" } else { "chr4" }
    let start = i * 1000 + 10000
    let end = start + 500
    peaks.push(Peak::new(chr, start, end, "peak_" + i.to_string()))
    i = i + 1
  }

  peaks
}

///|
pub fn cs_create_example_txdb() -> TxDb {
  let genes : Array[Gene] = [
    Gene::new("TP53", "chr1", 10000, 25000, Strand::Plus),
    Gene::new("BRCA1", "chr1", 50000, 70000, Strand::Minus),
    Gene::new("BRCA2", "chr2", 30000, 55000, Strand::Plus),
    Gene::new("ATM", "chr2", 80000, 100000, Strand::Minus),
    Gene::new("CDK4", "chr3", 20000, 35000, Strand::Plus),
    Gene::new("MYC", "chr3", 60000, 85000, Strand::Minus),
    Gene::new("EGFR", "chr4", 15000, 40000, Strand::Plus),
    Gene::new("GAPDH", "chr4", 90000, 110000, Strand::Minus),
  ]

  let transcripts : Array[Transcript] = [
    Transcript::new("TP53-001", "TP53", "chr1", 10000, 25000, Strand::Plus),
    Transcript::new("BRCA1-001", "BRCA1", "chr1", 50000, 70000, Strand::Minus),
    Transcript::new("BRCA2-001", "BRCA2", "chr2", 30000, 55000, Strand::Plus),
    Transcript::new("ATM-001", "ATM", "chr2", 80000, 100000, Strand::Minus),
    Transcript::new("CDK4-001", "CDK4", "chr3", 20000, 35000, Strand::Plus),
    Transcript::new("MYC-001", "MYC", "chr3", 60000, 85000, Strand::Minus),
    Transcript::new("EGFR-001", "EGFR", "chr4", 15000, 40000, Strand::Plus),
    Transcript::new("GAPDH-001", "GAPDH", "chr4", 90000, 110000, Strand::Minus),
  ]

  let exons : Array[Exon] = []
  let cds : Array[CDS] = []

  TxDb::new(genes, transcripts, exons, cds, "hg38", "Homo sapiens")
}
