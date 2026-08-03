///|
/// StructuralVariantAnnotation module demo - SV annotation from VCF.
///
/// Run: moon run examples/sv_demo
fn main {
  println("=== StructuralVariantAnnotation Module Demo ===")
  println("Structural variant parsing, partner pairing, and gene annotation")
  println("Reference: Bioconductor StructuralVariantAnnotation package")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample VCF
  // ------------------------------------------------------------------
  println("=== 1. Parse Sample VCF ===")
  let vcf = @bio.sv_sample_vcf()
  let records = @bio.sv_parse_vcf(vcf)
  println("Parsed " + records.length().to_string() + " SV records")
  println("")

  // ------------------------------------------------------------------
  // 2. Count by type
  // ------------------------------------------------------------------
  println("=== 2. SV Type Distribution ===")
  let counts = @bio.sv_count_by_type(records)
  let type_labels = ["DEL", "DUP", "INV", "INS", "TRA", "BND"]
  for label in type_labels {
    let count = counts.get_or_default(label, 0)
    println("  " + label + ": " + count.to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Show all records
  // ------------------------------------------------------------------
  println("=== 3. All SV Records ===")
  for i in 0..<records.length() {
    let rec = records[i]
    let summary = @bio.sv_record_summary(rec)
    println(summary)
  }

  // ------------------------------------------------------------------
  // 4. Partner pairing
  // ------------------------------------------------------------------
  println("=== 4. Breakend Partner Pairing ===")
  let pairs = @bio.sv_find_partners(records)
  println("Found " + pairs.length().to_string() + " partner pair(s)")
  for pair in pairs {
    let (a, b) = pair
    println("  " + a + " <-> " + b)
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Gene annotation
  // ------------------------------------------------------------------
  println("=== 5. Gene Overlap Annotation ===")
  let genes = @bio.sv_sample_genes()
  println("Gene regions: " + genes.length().to_string())
  for g in genes {
    println("  " + g.gene() + " (" + g.chrom() + ":" + g.start().to_string() + "-" + g.end().to_string() + ")")
  }
  println("")
  let annotations = @bio.sv_annotate(records, genes)
  for ann in annotations {
    let rec = ann.record()
    let gene_str = if ann.genes().length() > 0 {
      let mut s = ""
      for i in 0..<ann.genes().length() {
        if i > 0 {
          s = s + ", "
        }
        s = s + ann.genes()[i]
      }
      s
    } else {
      "(none)"
    }
    println(
      "  " +
        rec.id() +
        " (" +
        @bio.sv_type_label(rec.sv_type()) +
        "): " +
        gene_str,
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Filtering
  // ------------------------------------------------------------------
  println("=== 6. Filtering ===")
  let high_qual = @bio.sv_filter_by_quality(records, 400.0)
  println("Quality >= 400: " + high_qual.length().to_string() + " records")
  let large_sv = @bio.sv_filter_by_size(records, 1000)
  println("Size >= 1000: " + large_sv.length().to_string() + " records")
  let bnds = @bio.sv_filter_by_type(records, @bio.sv_type_bnd())
  println("BND records: " + bnds.length().to_string())
  println("")

  println("=== Demo complete ===")
}
