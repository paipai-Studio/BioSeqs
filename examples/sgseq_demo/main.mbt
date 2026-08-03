///|
/// sgseq module demo - Splice variant detection from RNA-seq data.
///
/// Demonstrates splice graph construction, junction parsing, and
/// detection of five alternative splicing event types (SE, A5SS,
/// A3SS, MXE, RI) from RNA-seq junction reads.
///
/// Reference: Goldstein LD et al. (2016) SGSeq. Bioinformatics 32:1188-1195.
/// Run: moon run examples/sgseq_demo
fn main {
  println("=== sgseq Module Demo ===")
  println("Splice variant detection from RNA-seq data")
  println("Reference: Bioconductor SGSeq")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample junction data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Junction Data ===")
  let junctions = @bio.sgseq_sample_junctions()
  println("Total junctions: " + junctions.length().to_string())
  for j in junctions {
    println(
      "  " +
        j.chr +
        ":" +
        j.start.to_string() +
        "-" +
        j.end.to_string() +
        " (" +
        j.strand +
        ") count=" +
        j.count.to_string() +
        " intron_len=" +
        j.intron_length().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Load sample exons
  // ------------------------------------------------------------------
  println("=== 2. Sample Exons ===")
  let exons = @bio.sgseq_sample_exons()
  println("Total exons: " + exons.length().to_string())
  for e in exons {
    println(
      "  " +
        e.chr +
        ":" +
        e.start.to_string() +
        "-" +
        e.end.to_string() +
        " (" +
        e.strand +
        ") length=" +
        e.length().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Build splice graph
  // ------------------------------------------------------------------
  println("=== 3. Build Splice Graph ===")
  let graph_exons = @bio.sgseq_build_splice_graph(junctions)
  println("Splice graph exons: " + graph_exons.length().to_string())
  for e in graph_exons {
    println(
      "  " +
        e.chr +
        ":" +
        e.start.to_string() +
        "-" +
        e.end.to_string() +
        " length=" +
        e.length().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Convert to features
  // ------------------------------------------------------------------
  println("=== 4. Splice Graph Features ===")
  let features = @bio.sgseq_to_features(junctions, graph_exons)
  println("Total features: " + features.length().to_string())
  let mut n_exon = 0
  let mut n_junc = 0
  for f in features {
    if f.type_ == "E" {
      n_exon = n_exon + 1
    } else {
      n_junc = n_junc + 1
    }
  }
  println("  Exon features: " + n_exon.to_string())
  println("  Junction features: " + n_junc.to_string())
  println("")

  // ------------------------------------------------------------------
  // 5. Detect all splice variants
  // ------------------------------------------------------------------
  println("=== 5. Detect Splice Variants ===")
  let variants = @bio.sgseq_detect_variants(junctions, exons)
  println("Total variants detected: " + variants.length().to_string())
  for v in variants {
    println("  " + @bio.sgseq_variant_to_string(v))
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Variant summary by type
  // ------------------------------------------------------------------
  println("=== 6. Variant Summary ===")
  let summary = @bio.sgseq_summary(variants)
  for key in summary.keys() {
    println("  " + key + ": " + summary[key].to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Parse junctions from text
  // ------------------------------------------------------------------
  println("=== 7. Parse Junctions from Text ===")
  let text = "chr1\t100\t200\t+\t50\nchr1\t300\t400\t+\t30\nchr1\t100\t350\t+\t10"
  let parsed = @bio.sgseq_parse_junctions(text)
  println("Parsed " + parsed.length().to_string() + " junctions from text:")
  for j in parsed {
    println(
      "  " +
        j.chr +
        ":" +
        j.start.to_string() +
        "-" +
        j.end.to_string() +
        " count=" +
        j.count.to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Filter variants by read support
  // ------------------------------------------------------------------
  println("=== 8. Filter Variants ===")
  let filtered = @bio.sgseq_filter_variants(variants, min_reads=10)
  println("Variants with >= 10 reads: " + filtered.length().to_string() + " (out of " + variants.length().to_string() + ")")
  println("")

  // ------------------------------------------------------------------
  // 9. Specific variant type detection
  // ------------------------------------------------------------------
  println("=== 9. Individual Variant Type Detection ===")
  let se_variants = @bio.sgseq_detect_se(junctions, exons)
  println("SE variants: " + se_variants.length().to_string())

  let a5ss_variants = @bio.sgseq_detect_a5ss(junctions)
  println("A5SS variants: " + a5ss_variants.length().to_string())

  let a3ss_variants = @bio.sgseq_detect_a3ss(junctions)
  println("A3SS variants: " + a3ss_variants.length().to_string())

  let mxe_variants = @bio.sgseq_detect_mxe(junctions, exons)
  println("MXE variants: " + mxe_variants.length().to_string())

  let ri_variants = @bio.sgseq_detect_ri(junctions, exons)
  println("RI variants: " + ri_variants.length().to_string())
  println("")

  println("=== Demo complete ===")
}
