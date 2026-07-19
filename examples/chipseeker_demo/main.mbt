///|
fn main {
  let peaks = @bio.cs_create_example_peaks()
  let txdb = @bio.cs_create_example_txdb()

  println("ChIPseeker Peak Annotation Demo")
  println("==============================")
  println("")

  let result = @bio.annotate_peaks(peaks, txdb)
  println("Total peaks annotated: " + result.annotations.length().to_string())
  println("")

  let summary = @bio.annotation_summary(result)
  println(summary)

  let bar_plot = @bio.plot_annotation_bar(result.summary)
  println("")
  println(bar_plot)

  let dist_plot = @bio.plot_dist_to_tss(result.annotations)
  println(dist_plot)

  let promoters = @bio.get_promoter_annotations(result)
  println("Promoter peaks: " + promoters.length().to_string())

  let exons = @bio.get_exon_annotations(result)
  println("Exon peaks: " + exons.length().to_string())

  let introns = @bio.get_intron_annotations(result)
  println("Intron peaks: " + introns.length().to_string())

  let utrs = @bio.get_utr_annotations(result)
  println("UTR peaks: " + utrs.length().to_string())

  let distal = @bio.get_distal_annotations(result)
  println("Distal peaks: " + distal.length().to_string())

  let intergenic = @bio.get_intergenic_annotations(result)
  println("Intergenic peaks: " + intergenic.length().to_string())

  println("")
  println("=== Advanced Features ===")
  println("")

  println("1. read_peak_file - BED format parsing")
  let bed_content = "chr1\t10000\t11000\tTP53_promoter\nchr1\t50000\t51000\tBRCA1_peak\nchr2\t35000\t36000\tBRCA2_exon"
  let parsed_peaks = @bio.read_peak_file(bed_content)
  println("   Peaks parsed from BED: " + parsed_peaks.length().to_string())
  println("")

  println("2. peak2gene - Peak to gene associations")
  let associations = @bio.peak2gene(result)
  println("   Total associations: " + associations.length().to_string())
  if associations.length() > 0 {
    println("   First 3 associations:")
    let mut i = 0
    while i < 3 && i < associations.length() {
      let assoc = associations[i]
      println("     Peak: " + assoc.peak_id + " -> Gene: " + assoc.gene_id + " (" + assoc.feature_type + ")")
      i = i + 1
    }
  }
  println("")

  println("=== Demo completed! ===")
}
