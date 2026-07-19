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

  let distal = @bio.get_distal_annotations(result)
  println("Distal peaks: " + distal.length().to_string())

  let intergenic = @bio.get_intergenic_annotations(result)
  println("Intergenic peaks: " + intergenic.length().to_string())
}
