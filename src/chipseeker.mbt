///|
/// ChIPseeker - ChIP-seq peak annotation module.
///
/// Inspired by the Bioconductor ChIPseeker package. This implementation
/// focuses on peak annotation: mapping ChIP-seq peaks to nearby genes and
/// genomic features (Promoter, 5'UTR, 3'UTR, Exon, Intron, Downstream,
/// Distal Intergenic).
///
/// Feature assignment uses the ChIPseeker priority order:
///   Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Distal Intergenic

///| Length of the 5'UTR / 3'UTR approximation (in bp).
const CS_UTR_LENGTH : Int = 500

///| Size of each implicit "exon" segment inside the gene body (in bp).
const CS_EXON_SIZE : Int = 1000

///| Size of the downstream region past the gene body (in bp).
const CS_DOWNSTREAM_WINDOW : Int = 3000

///| A ChIP-seq peak (genomic interval).
pub struct Peak {
  chr : String
  start : Int
  mut end : Int
  summit : Int
  fold_enrichment : Double
  p_value : Double
}

///| Construct a Peak. `chr` and `start` are required; all other fields
///| are optional and default to zero (summit=0 means "use peak center").
pub fn Peak::new(
  chr~ : String,
  start~ : Int,
  end? : Int = 0,
  summit? : Int = 0,
  fold_enrichment? : Double = 0.0,
  p_value? : Double = 0.0,
) -> Peak {
  Peak::{ chr, start, end, summit, fold_enrichment, p_value }
}

///| Get the chromosome of the peak.
pub fn Peak::get_chr(self : Peak) -> String {
  self.chr
}

///| Get the start position of the peak.
pub fn Peak::get_start(self : Peak) -> Int {
  self.start
}

///| Get the end position of the peak.
pub fn Peak::get_end(self : Peak) -> Int {
  self.end
}

///| Set the end position of the peak.
pub fn Peak::set_end(self : Peak, new_end : Int) -> Unit {
  self.end = new_end
}

///| Get the summit position of the peak.
pub fn Peak::get_summit(self : Peak) -> Int {
  self.summit
}

///| Get the fold enrichment of the peak.
pub fn Peak::get_fold_enrichment(self : Peak) -> Double {
  self.fold_enrichment
}

///| Get the p-value of the peak.
pub fn Peak::get_p_value(self : Peak) -> Double {
  self.p_value
}

///| A gene annotation record used for peak annotation.
pub struct GeneAnnotation {
  gene_id : String
  gene_name : String
  chr : String
  gene_start : Int
  gene_end : Int
  strand : String // "+" or "-"
  tss : Int // transcription start site
}

///| Construct a GeneAnnotation. The TSS is auto-computed from the strand:
///|   "+" strand -> tss = gene_start
///|   "-" strand -> tss = gene_end
pub fn GeneAnnotation::new(
  gene_id~ : String,
  gene_name~ : String,
  chr~ : String,
  gene_start~ : Int,
  gene_end~ : Int,
  strand~ : String,
) -> GeneAnnotation {
  let tss = if strand == "-" { gene_end } else { gene_start }
  GeneAnnotation::{ gene_id, gene_name, chr, gene_start, gene_end, strand, tss }
}

///| Get the gene id.
pub fn GeneAnnotation::get_gene_id(self : GeneAnnotation) -> String {
  self.gene_id
}

///| Get the gene name.
pub fn GeneAnnotation::get_gene_name(self : GeneAnnotation) -> String {
  self.gene_name
}

///| Get the chromosome.
pub fn GeneAnnotation::get_chr(self : GeneAnnotation) -> String {
  self.chr
}

///| Get the gene start.
pub fn GeneAnnotation::get_gene_start(self : GeneAnnotation) -> Int {
  self.gene_start
}

///| Get the gene end.
pub fn GeneAnnotation::get_gene_end(self : GeneAnnotation) -> Int {
  self.gene_end
}

///| Get the strand.
pub fn GeneAnnotation::get_strand(self : GeneAnnotation) -> String {
  self.strand
}

///| Get the transcription start site (TSS).
pub fn GeneAnnotation::get_tss(self : GeneAnnotation) -> Int {
  self.tss
}

///| The result of annotating a single peak.
pub struct PeakAnnotation {
  peak : Peak
  gene_id : String
  gene_name : String
  chr : String
  tss : Int
  distance_to_tss : Int
  feature : String // "Promoter", "5'UTR", etc.
  annotation : String // detailed annotation string
}

///| Construct a PeakAnnotation.
pub fn PeakAnnotation::new(
  peak~ : Peak,
  gene_id~ : String,
  gene_name~ : String,
  chr~ : String,
  tss~ : Int,
  distance_to_tss~ : Int,
  feature~ : String,
  annotation~ : String,
) -> PeakAnnotation {
  PeakAnnotation::{
    peak,
    gene_id,
    gene_name,
    chr,
    tss,
    distance_to_tss,
    feature,
    annotation,
  }
}

///| Get the peak.
pub fn PeakAnnotation::get_peak(self : PeakAnnotation) -> Peak {
  self.peak
}

///| Get the gene id.
pub fn PeakAnnotation::get_gene_id(self : PeakAnnotation) -> String {
  self.gene_id
}

///| Get the gene name.
pub fn PeakAnnotation::get_gene_name(self : PeakAnnotation) -> String {
  self.gene_name
}

///| Get the chromosome.
pub fn PeakAnnotation::get_chr(self : PeakAnnotation) -> String {
  self.chr
}

///| Get the TSS.
pub fn PeakAnnotation::get_tss(self : PeakAnnotation) -> Int {
  self.tss
}

///| Get the signed distance to TSS.
pub fn PeakAnnotation::get_distance_to_tss(self : PeakAnnotation) -> Int {
  self.distance_to_tss
}

///| Get the assigned feature string.
pub fn PeakAnnotation::get_feature(self : PeakAnnotation) -> String {
  self.feature
}

///| Get the detailed annotation string.
pub fn PeakAnnotation::get_annotation(self : PeakAnnotation) -> String {
  self.annotation
}

///| Summary counts of peaks assigned to each feature type.
pub struct AnnotationSummary {
  total_peaks : Int
  promoter : Int
  utr5 : Int
  utr3 : Int
  exon : Int
  intron : Int
  downstream : Int
  distal_intergenic : Int
}

///| Construct an AnnotationSummary.
pub fn AnnotationSummary::new(
  total_peaks~ : Int,
  promoter~ : Int,
  utr5~ : Int,
  utr3~ : Int,
  exon~ : Int,
  intron~ : Int,
  downstream~ : Int,
  distal_intergenic~ : Int,
) -> AnnotationSummary {
  AnnotationSummary::{
    total_peaks,
    promoter,
    utr5,
    utr3,
    exon,
    intron,
    downstream,
    distal_intergenic,
  }
}

///| Get total peaks.
pub fn AnnotationSummary::get_total_peaks(self : AnnotationSummary) -> Int {
  self.total_peaks
}

///| Get promoter count.
pub fn AnnotationSummary::get_promoter(self : AnnotationSummary) -> Int {
  self.promoter
}

///| Get 5'UTR count.
pub fn AnnotationSummary::get_utr5(self : AnnotationSummary) -> Int {
  self.utr5
}

///| Get 3'UTR count.
pub fn AnnotationSummary::get_utr3(self : AnnotationSummary) -> Int {
  self.utr3
}

///| Get exon count.
pub fn AnnotationSummary::get_exon(self : AnnotationSummary) -> Int {
  self.exon
}

///| Get intron count.
pub fn AnnotationSummary::get_intron(self : AnnotationSummary) -> Int {
  self.intron
}

///| Get downstream count.
pub fn AnnotationSummary::get_downstream(self : AnnotationSummary) -> Int {
  self.downstream
}

///| Get distal intergenic count.
pub fn AnnotationSummary::get_distal_intergenic(self : AnnotationSummary) -> Int {
  self.distal_intergenic
}

///| Compute the effective single position used for a peak. Uses the summit
///| if it is non-zero, otherwise the midpoint of [start, end].
pub fn cs_peak_position(peak : Peak) -> Int {
  if peak.summit != 0 {
    peak.summit
  } else {
    (peak.start + peak.end) / 2
  }
}

///| Build a sentinel "no gene" annotation. Used when no gene is found on
///| the same chromosome as the peak.
fn cs_na_gene(chr : String) -> GeneAnnotation {
  GeneAnnotation::{
    gene_id: "",
    gene_name: "",
    chr,
    gene_start: 0,
    gene_end: 0,
    strand: "+",
    tss: 0,
  }
}

///| Compute the signed distance from a peak to a gene's TSS.
///|
///| Positive distance means the peak is downstream of the TSS (in the
///| transcription direction); negative means upstream. For "-" strand
///| genes the orientation is flipped so the sign convention still holds.
pub fn compute_tss_distance(peak : Peak, gene : GeneAnnotation) -> Int {
  if gene.gene_id == "" {
    return 0
  }
  let peak_pos = cs_peak_position(peak)
  if gene.strand == "-" {
    gene.tss - peak_pos
  } else {
    peak_pos - gene.tss
  }
}

///| Find the gene whose TSS is closest to the peak. Only genes on the same
///| chromosome as the peak are considered. If no such gene exists, a
///| sentinel annotation with empty `gene_id` is returned.
pub fn find_nearest_gene(
  peak : Peak,
  genes : Array[GeneAnnotation],
) -> GeneAnnotation {
  let peak_pos = cs_peak_position(peak)
  let mut nearest : GeneAnnotation = cs_na_gene(peak.chr)
  let mut min_dist : Int = 0
  let mut found = false
  for gene in genes {
    if gene.chr != peak.chr {
      continue
    }
    let dist = (gene.tss - peak_pos).abs()
    if !found || dist < min_dist {
      min_dist = dist
      nearest = gene
      found = true
    }
  }
  nearest
}

///| Assign a genomic feature to a peak relative to a gene, following the
///| ChIPseeker priority:
///|   Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Distal Intergenic
pub fn assign_feature(
  peak : Peak,
  gene : GeneAnnotation,
  promoter_window : Int,
) -> String {
  // No real gene -> everything is distal intergenic.
  if gene.gene_id == "" || gene.chr != peak.chr {
    return "Distal Intergenic"
  }
  let peak_pos = cs_peak_position(peak)
  let distance = compute_tss_distance(peak, gene)
  // 1. Promoter: within +/- promoter_window of the TSS.
  if distance >= -promoter_window && distance <= promoter_window {
    return "Promoter"
  }
  // Compute strand-aware region boundaries. Each gene is approximated as
  // having a 5'UTR and 3'UTR at the two ends of the gene body, two implicit
  // "exon" segments just inside the UTRs, and a single "intron" segment in
  // the middle of the gene body.
  let utr_len = CS_UTR_LENGTH
  let exon_size = CS_EXON_SIZE
  let (utr5_start, utr5_end, utr3_start, utr3_end, exon1_start, exon1_end, exon2_start, exon2_end, intron_start, intron_end, down_start, down_end) = if gene.strand == "-" {
    // "-" strand: TSS = gene_end; transcription proceeds toward gene_start.
    let utr5_s = gene.gene_end - utr_len
    let utr5_e = gene.gene_end
    let utr3_s = gene.gene_start
    let utr3_e = gene.gene_start + utr_len
    let exon1_s = gene.gene_end - utr_len - exon_size
    let exon1_e = gene.gene_end - utr_len
    let exon2_s = gene.gene_start + utr_len
    let exon2_e = gene.gene_start + utr_len + exon_size
    let intron_s = gene.gene_start + utr_len + exon_size
    let intron_e = gene.gene_end - utr_len - exon_size
    let down_s = gene.gene_start - CS_DOWNSTREAM_WINDOW
    let down_e = gene.gene_start
    (
      utr5_s,
      utr5_e,
      utr3_s,
      utr3_e,
      exon1_s,
      exon1_e,
      exon2_s,
      exon2_e,
      intron_s,
      intron_e,
      down_s,
      down_e,
    )
  } else {
    // "+" strand: TSS = gene_start; transcription proceeds toward gene_end.
    let utr5_s = gene.gene_start
    let utr5_e = gene.gene_start + utr_len
    let utr3_s = gene.gene_end - utr_len
    let utr3_e = gene.gene_end
    let exon1_s = gene.gene_start + utr_len
    let exon1_e = gene.gene_start + utr_len + exon_size
    let exon2_s = gene.gene_end - utr_len - exon_size
    let exon2_e = gene.gene_end - utr_len
    let intron_s = gene.gene_start + utr_len + exon_size
    let intron_e = gene.gene_end - utr_len - exon_size
    let down_s = gene.gene_end
    let down_e = gene.gene_end + CS_DOWNSTREAM_WINDOW
    (
      utr5_s,
      utr5_e,
      utr3_s,
      utr3_e,
      exon1_s,
      exon1_e,
      exon2_s,
      exon2_e,
      intron_s,
      intron_e,
      down_s,
      down_e,
    )
  }
  // 2. 5'UTR
  if peak_pos >= utr5_start && peak_pos <= utr5_end {
    return "5'UTR"
  }
  // 3. 3'UTR
  if peak_pos >= utr3_start && peak_pos <= utr3_end {
    return "3'UTR"
  }
  // 4. Exon (either of the two implicit exon segments)
  if (peak_pos >= exon1_start && peak_pos <= exon1_end) ||
     (peak_pos >= exon2_start && peak_pos <= exon2_end) {
    return "Exon"
  }
  // 5. Intron
  if peak_pos >= intron_start && peak_pos <= intron_end {
    return "Intron"
  }
  // 6. Downstream
  if peak_pos >= down_start && peak_pos <= down_end {
    return "Downstream"
  }
  // 7. Distal Intergenic
  "Distal Intergenic"
}

///| Format the detailed annotation string for a peak/gene/feature triple.
fn cs_format_annotation(
  gene : GeneAnnotation,
  distance : Int,
  feature : String,
  promoter_window : Int,
) -> String {
  match feature {
    "Promoter" =>
      "Promoter ([-\{promoter_window}, +\{promoter_window}] relative to TSS of \{gene.gene_name})"
    "5'UTR" => "5' UTR of \{gene.gene_name}"
    "3'UTR" => "3' UTR of \{gene.gene_name}"
    "Exon" => "Exon of \{gene.gene_name}"
    "Intron" => "Intron of \{gene.gene_name}"
    "Downstream" =>
      "Downstream (<=\{CS_DOWNSTREAM_WINDOW}bp) of \{gene.gene_name}"
    _ => {
      if gene.gene_id == "" {
        "Distal Intergenic (no nearby gene on same chromosome)"
      } else {
        "Distal Intergenic (nearest gene: \{gene.gene_name}, distance to TSS: \{distance})"
      }
    }
  }
}

///| Annotate a single peak: find the nearest gene, compute the signed
///| distance to its TSS, assign a genomic feature, and build the detailed
///| annotation string.
pub fn annotate_peak(
  peak : Peak,
  genes : Array[GeneAnnotation],
  promoter_window : Int,
) -> PeakAnnotation {
  let gene = find_nearest_gene(peak, genes)
  let distance = compute_tss_distance(peak, gene)
  let feature = assign_feature(peak, gene, promoter_window)
  let annotation = cs_format_annotation(gene, distance, feature, promoter_window)
  PeakAnnotation::new(
    peak=peak,
    gene_id=gene.gene_id,
    gene_name=gene.gene_name,
    chr=peak.chr,
    tss=gene.tss,
    distance_to_tss=distance,
    feature=feature,
    annotation=annotation,
  )
}

///| Annotate a collection of peaks against a set of gene annotations.
pub fn annotate_peaks(
  peaks : Array[Peak],
  genes : Array[GeneAnnotation],
  promoter_window : Int,
) -> Array[PeakAnnotation] {
  let result : Array[PeakAnnotation] = Array::new()
  for peak in peaks {
    result.push(annotate_peak(peak, genes, promoter_window))
  }
  result
}

///| Tally the feature assignments of a list of peak annotations.
pub fn annotation_summary(
  annotations : Array[PeakAnnotation],
) -> AnnotationSummary {
  let mut promoter = 0
  let mut utr5 = 0
  let mut utr3 = 0
  let mut exon = 0
  let mut intron = 0
  let mut downstream = 0
  let mut distal_intergenic = 0
  for ann in annotations {
    match ann.feature {
      "Promoter" => promoter = promoter + 1
      "5'UTR" => utr5 = utr5 + 1
      "3'UTR" => utr3 = utr3 + 1
      "Exon" => exon = exon + 1
      "Intron" => intron = intron + 1
      "Downstream" => downstream = downstream + 1
      "Distal Intergenic" => distal_intergenic = distal_intergenic + 1
      _ => ()
    }
  }
  AnnotationSummary::new(
    total_peaks=annotations.length(),
    promoter=promoter,
    utr5=utr5,
    utr3=utr3,
    exon=exon,
    intron=intron,
    downstream=downstream,
    distal_intergenic=distal_intergenic,
  )
}

///| Return only those annotations whose feature matches `feature`.
pub fn filter_by_feature(
  annotations : Array[PeakAnnotation],
  feature : String,
) -> Array[PeakAnnotation] {
  let result : Array[PeakAnnotation] = Array::new()
  for ann in annotations {
    if ann.feature == feature {
      result.push(ann)
    }
  }
  result
}

///| Return only those annotations whose absolute distance to TSS is at
///| most `max_distance`.
pub fn filter_by_distance(
  annotations : Array[PeakAnnotation],
  max_distance : Int,
) -> Array[PeakAnnotation] {
  let result : Array[PeakAnnotation] = Array::new()
  for ann in annotations {
    if ann.distance_to_tss.abs() <= max_distance {
      result.push(ann)
    }
  }
  result
}

///| Format an AnnotationSummary as a human-readable string.
pub fn chipseeker_summary_string(summary : AnnotationSummary) -> String {
  let mut s = "ChIPseeker Annotation Summary\n"
  s = s + "=============================\n"
  s = s + "Total peaks: \{summary.total_peaks}\n"
  s = s + "Promoter: \{summary.promoter}\n"
  s = s + "5'UTR: \{summary.utr5}\n"
  s = s + "3'UTR: \{summary.utr3}\n"
  s = s + "Exon: \{summary.exon}\n"
  s = s + "Intron: \{summary.intron}\n"
  s = s + "Downstream: \{summary.downstream}\n"
  s = s + "Distal Intergenic: \{summary.distal_intergenic}\n"
  s
}

///| Sample ChIP-seq peaks for testing and demonstration. The peaks are
///| laid out so that, when annotated against `chipseeker_sample_genes()`
///| with the default promoter window (3000), they exercise every feature
///| type except 5'UTR (which is only reachable with a smaller window
///| because the 5'UTR sits inside the default promoter window).
pub fn chipseeker_sample_peaks() -> Array[Peak] {
  [
    // chr1, GeneA (+ strand, TSS=10000, gene_end=25000)
    Peak::new(chr="chr1", start=7900, end=8100, summit=8000), // Promoter (dist=-2000)
    Peak::new(chr="chr1", start=11900, end=12100, summit=12000), // Promoter (dist=+2000)
    Peak::new(chr="chr1", start=13900, end=14100, summit=14000), // Intron (dist=+4000)
    Peak::new(chr="chr1", start=23900, end=24100, summit=24000), // Exon (dist=+14000)
    Peak::new(chr="chr1", start=24600, end=24800, summit=24700), // 3'UTR (dist=+14700)
    Peak::new(chr="chr1", start=25900, end=26100, summit=26000), // Downstream (dist=+16000)
    Peak::new(chr="chr1", start=38900, end=39100, summit=39000), // Distal Intergenic
    // chr1, GeneB (- strand, TSS=70000, gene_start=50000)
    Peak::new(chr="chr1", start=50900, end=51100, summit=51000), // Exon (dist=+19000)
    Peak::new(chr="chr1", start=59900, end=60100, summit=60000), // Intron (dist=+10000)
    Peak::new(chr="chr1", start=68900, end=69100, summit=69000), // Promoter (dist=+1000)
    Peak::new(chr="chr1", start=48900, end=49100, summit=49000), // Downstream (dist=+21000)
    // chr2, GeneC (+ strand, TSS=20000, gene_end=40000)
    Peak::new(chr="chr2", start=19900, end=20100, summit=20000), // Promoter (dist=0, at TSS)
    Peak::new(chr="chr2", start=29900, end=30100, summit=30000), // Intron (dist=+10000)
    // chr3, no gene -> Distal Intergenic with empty gene_id
    Peak::new(chr="chr3", start=10000, end=10200, summit=10100),
  ]
}

///| Sample gene annotations for testing and demonstration.
pub fn chipseeker_sample_genes() -> Array[GeneAnnotation] {
  [
    GeneAnnotation::new(
      gene_id="GENE1",
      gene_name="GeneA",
      chr="chr1",
      gene_start=10000,
      gene_end=25000,
      strand="+",
    ),
    GeneAnnotation::new(
      gene_id="GENE2",
      gene_name="GeneB",
      chr="chr1",
      gene_start=50000,
      gene_end=70000,
      strand="-",
    ),
    GeneAnnotation::new(
      gene_id="GENE3",
      gene_name="GeneC",
      chr="chr2",
      gene_start=20000,
      gene_end=40000,
      strand="+",
    ),
  ]
}
