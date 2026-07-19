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
  transcript_id : String
  exon_id : String
  utr_type : String
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
  transcript_id : String,
  exon_id : String,
  utr_type : String,
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
    transcript_id,
    exon_id,
    utr_type,
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

  let (is_utr, utr_tx_id, utr_gene_id, utr_type) = is_inside_utr(peak, txdb)
  if is_utr {
    return PeakAnnotation::new(
      peak.peak_id,
      peak.chr,
      peak.start,
      peak.end,
      utr_type,
      utr_gene_id,
      0,
      utr_type,
      strand,
      utr_tx_id,
      "",
      utr_type,
    )
  }

  let (is_exon, exon_id, exon_gene_id) = is_inside_exon(peak, txdb)
  if is_exon {
    return PeakAnnotation::new(
      peak.peak_id,
      peak.chr,
      peak.start,
      peak.end,
      "Exon",
      exon_gene_id,
      0,
      "Exon",
      strand,
      "",
      exon_id,
      "",
    )
  }

  let (is_intron, intron_tx_id, intron_gene_id) = is_inside_intron(peak, txdb)
  if is_intron {
    return PeakAnnotation::new(
      peak.peak_id,
      peak.chr,
      peak.start,
      peak.end,
      "Intron",
      intron_gene_id,
      0,
      "Intron",
      strand,
      intron_tx_id,
      "",
      "",
    )
  }

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
    "",
    "",
    "",
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
fn is_inside_exon(peak : Peak, txdb : TxDb) -> (Bool, String, String) {
  let peak_center = (peak.start + peak.end) / 2
  
  let mut i = 0
  while i < txdb.exons.length() {
    let exon = txdb.exons[i]
    if exon.seqname == peak.chr && peak_center >= exon.start && peak_center <= exon.end {
      return (true, exon.exon_id, exon.gene_id)
    }
    i = i + 1
  }
  
  (false, "", "")
}

///|
fn is_inside_intron(peak : Peak, txdb : TxDb) -> (Bool, String, String) {
  let peak_center = (peak.start + peak.end) / 2
  
  let mut i = 0
  while i < txdb.transcripts.length() {
    let tx = txdb.transcripts[i]
    if tx.seqname == peak.chr && peak_center >= tx.start && peak_center <= tx.end {
      let mut inside_exon = false
      let mut j = 0
      while j < txdb.exons.length() {
        let exon = txdb.exons[j]
        if exon.transcript_id == tx.transcript_id && 
           peak_center >= exon.start && peak_center <= exon.end {
          inside_exon = true
          break
        }
        j = j + 1
      }
      if !inside_exon {
        return (true, tx.transcript_id, tx.gene_id)
      }
    }
    i = i + 1
  }
  
  (false, "", "")
}

///|
fn is_inside_utr(peak : Peak, txdb : TxDb) -> (Bool, String, String, String) {
  let peak_center = (peak.start + peak.end) / 2
  let utr_length = 500
  
  let mut i = 0
  while i < txdb.transcripts.length() {
    let tx = txdb.transcripts[i]
    if tx.seqname != peak.chr {
      i = i + 1
      continue
    }
    
    match tx.strand {
      Strand::Plus => {
        let utr5_start = tx.start
        let utr5_end = tx.start + utr_length
        if peak_center >= utr5_start && peak_center <= utr5_end {
          return (true, tx.transcript_id, tx.gene_id, "5'UTR")
        }
        let utr3_start = tx.end - utr_length
        let utr3_end = tx.end
        if peak_center >= utr3_start && peak_center <= utr3_end {
          return (true, tx.transcript_id, tx.gene_id, "3'UTR")
        }
      }
      Strand::Minus => {
        let utr5_start = tx.end - utr_length
        let utr5_end = tx.end
        if peak_center >= utr5_start && peak_center <= utr5_end {
          return (true, tx.transcript_id, tx.gene_id, "5'UTR")
        }
        let utr3_start = tx.start
        let utr3_end = tx.start + utr_length
        if peak_center >= utr3_start && peak_center <= utr3_end {
          return (true, tx.transcript_id, tx.gene_id, "3'UTR")
        }
      }
      _ => ()
    }
    i = i + 1
  }
  
  (false, "", "", "")
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
  let exons = get_exon_annotations(result)
  let introns = get_intron_annotations(result)
  let utrs = get_utr_annotations(result)

  result_str = result_str + "\nPromoter peaks: " + promoters.length().to_string() + "\n"
  result_str = result_str + "Exon peaks: " + exons.length().to_string() + "\n"
  result_str = result_str + "Intron peaks: " + introns.length().to_string() + "\n"
  result_str = result_str + "UTR peaks: " + utrs.length().to_string() + "\n"
  result_str = result_str + "Distal peaks: " + distal.length().to_string() + "\n"
  result_str = result_str + "Intergenic peaks: " + intergenic.length().to_string() + "\n"

  result_str
}

///|
pub fn read_peak_file(file_content : String) -> Array[Peak] {
  let peaks : Array[Peak] = Array::new()
  let lines = file_content.split("\n")
  let mut peak_id = 0

  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 || trimmed.has_prefix("#") {
      continue
    }

    let fields_iter = trimmed.split("\t")
    let fields : Array[String] = Array::new()
    for f in fields_iter {
      fields.push(f.to_string())
    }

    if fields.length() >= 3 {
      let chr = fields[0]
      let start = parse_int(fields[1])
      let end = parse_int(fields[2])
      let id = if fields.length() >= 4 { fields[3] } else { "peak_" + peak_id.to_string() }
      peaks.push(Peak::new(chr, start, end, id))
      peak_id = peak_id + 1
    }
  }

  peaks
}

///|
pub struct PeakGeneAssociation {
  peak_id : String
  gene_id : String
  distance : Int
  feature_type : String
  strand : String
}

///|
pub fn PeakGeneAssociation::new(
  peak_id : String,
  gene_id : String,
  distance : Int,
  feature_type : String,
  strand : String,
) -> PeakGeneAssociation {
  PeakGeneAssociation::{ peak_id, gene_id, distance, feature_type, strand }
}

///|
pub fn peak2gene(result : AnnotationResult) -> Array[PeakGeneAssociation] {
  let associations : Array[PeakGeneAssociation] = Array::new()

  for annotation in result.annotations {
    associations.push(PeakGeneAssociation::new(
      annotation.peak_id,
      annotation.gene_name,
      annotation.distance_to_tss,
      annotation.feature_type,
      annotation.strand,
    ))
  }

  associations
}

///|
pub fn get_exon_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let exons : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Exon" {
      exons.push(annotation)
    }
    i = i + 1
  }

  exons
}

///|
pub fn get_intron_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let introns : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Intron" {
      introns.push(annotation)
    }
    i = i + 1
  }

  introns
}

///|
pub fn get_utr_annotations(result : AnnotationResult) -> Array[PeakAnnotation] {
  let utrs : Array[PeakAnnotation] = Array::new()

  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.utr_type != "" {
      utrs.push(annotation)
    }
    i = i + 1
  }

  utrs
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

  let exons : Array[Exon] = [
    Exon::new("TP53-exon1", "TP53-001", "TP53", "chr1", 10000, 11000, Strand::Plus, 1),
    Exon::new("TP53-exon2", "TP53-001", "TP53", "chr1", 13000, 14000, Strand::Plus, 2),
    Exon::new("TP53-exon3", "TP53-001", "TP53", "chr1", 16000, 17000, Strand::Plus, 3),
    Exon::new("TP53-exon4", "TP53-001", "TP53", "chr1", 20000, 25000, Strand::Plus, 4),
    Exon::new("BRCA1-exon1", "BRCA1-001", "BRCA1", "chr1", 50000, 51000, Strand::Minus, 1),
    Exon::new("BRCA1-exon2", "BRCA1-001", "BRCA1", "chr1", 53000, 54000, Strand::Minus, 2),
    Exon::new("BRCA1-exon3", "BRCA1-001", "BRCA1", "chr1", 65000, 70000, Strand::Minus, 3),
    Exon::new("BRCA2-exon1", "BRCA2-001", "BRCA2", "chr2", 30000, 31000, Strand::Plus, 1),
    Exon::new("BRCA2-exon2", "BRCA2-001", "BRCA2", "chr2", 35000, 36000, Strand::Plus, 2),
    Exon::new("BRCA2-exon3", "BRCA2-001", "BRCA2", "chr2", 50000, 55000, Strand::Plus, 3),
    Exon::new("ATM-exon1", "ATM-001", "ATM", "chr2", 80000, 81000, Strand::Minus, 1),
    Exon::new("ATM-exon2", "ATM-001", "ATM", "chr2", 85000, 86000, Strand::Minus, 2),
    Exon::new("ATM-exon3", "ATM-001", "ATM", "chr2", 95000, 100000, Strand::Minus, 3),
    Exon::new("CDK4-exon1", "CDK4-001", "CDK4", "chr3", 20000, 21000, Strand::Plus, 1),
    Exon::new("CDK4-exon2", "CDK4-001", "CDK4", "chr3", 25000, 26000, Strand::Plus, 2),
    Exon::new("CDK4-exon3", "CDK4-001", "CDK4", "chr3", 30000, 35000, Strand::Plus, 3),
    Exon::new("MYC-exon1", "MYC-001", "MYC", "chr3", 60000, 61000, Strand::Minus, 1),
    Exon::new("MYC-exon2", "MYC-001", "MYC", "chr3", 70000, 71000, Strand::Minus, 2),
    Exon::new("MYC-exon3", "MYC-001", "MYC", "chr3", 80000, 85000, Strand::Minus, 3),
    Exon::new("EGFR-exon1", "EGFR-001", "EGFR", "chr4", 15000, 16000, Strand::Plus, 1),
    Exon::new("EGFR-exon2", "EGFR-001", "EGFR", "chr4", 20000, 21000, Strand::Plus, 2),
    Exon::new("EGFR-exon3", "EGFR-001", "EGFR", "chr4", 35000, 40000, Strand::Plus, 3),
    Exon::new("GAPDH-exon1", "GAPDH-001", "GAPDH", "chr4", 90000, 91000, Strand::Minus, 1),
    Exon::new("GAPDH-exon2", "GAPDH-001", "GAPDH", "chr4", 95000, 96000, Strand::Minus, 2),
    Exon::new("GAPDH-exon3", "GAPDH-001", "GAPDH", "chr4", 105000, 110000, Strand::Minus, 3),
  ]

  let cds : Array[CDS] = [
    CDS::new("TP53-cds1", "TP53-001", "TP53", "chr1", 10200, 11000, Strand::Plus, 1),
    CDS::new("TP53-cds2", "TP53-001", "TP53", "chr1", 13000, 14000, Strand::Plus, 2),
    CDS::new("TP53-cds3", "TP53-001", "TP53", "chr1", 16000, 17000, Strand::Plus, 3),
    CDS::new("TP53-cds4", "TP53-001", "TP53", "chr1", 20000, 24800, Strand::Plus, 4),
  ]

  TxDb::new(genes, transcripts, exons, cds, "hg38", "Homo sapiens")
}
