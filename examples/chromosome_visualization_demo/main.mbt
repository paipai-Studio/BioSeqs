///|
/// Chromosome Visualization Demo - SVG-based chromosome rendering.
/// Inspired by Biopython's Bio.Graphics.Chromosome module.

///|
fn main {
  println("=== Chromosome Visualization Demo ===")
  println("Inspired by Biopython's Bio.Graphics.Chromosome module")
  println("")

  println("1. Creating Chromosome Features")
  let gene1 = @bio.ChrFeature::new(
    name="BRCA1",
    start=1000.0,
    end=5000.0,
  )
  let gene2 = @bio.ChrFeature::new(
    name="TP53",
    start=10000.0,
    end=15000.0,
  )
  let gene3 = @bio.ChrFeature::new(
    name="EGFR",
    start=20000.0,
    end=25000.0,
  )
  println("   Created features:")
  println("     - BRCA1")
  println("     - TP53")
  println("     - EGFR")
  println("")

  println("2. Creating a Chromosome")
  let chr1 = @bio.Chromosome::new(
    name="chr1",
    length=30000.0,
    centromere_pos=15000.0,
  )
  println("   Chromosome: chr1")
  println("   Length: 30000 bp")
  println("   Centromere position: 15000")
  println("")

  println("3. Adding Features to Chromosome")
  let chr1_with_features = chr1
    .add_feature(gene1)
    .add_feature(gene2)
    .add_feature(gene3)
  println("   Added 3 features to chromosome")
  println("   Total features: " + chr1_with_features.features.length().to_string())
  println("")

  println("4. Adding Regions and Bands")
  let region1 = @bio.ChrRegion::new(
    name="Region1",
    start=0.0,
    end=8000.0,
  )
  let region2 = @bio.ChrRegion::new(
    name="Region2",
    start=18000.0,
    end=28000.0,
  )
  let chr_full = chr1_with_features
    .add_region(region1)
    .add_region(region2)
    .add_band("p36", 0.0, 2000.0, "gpos75")
    .add_band("p35", 2000.0, 5000.0, "gpos50")
    .add_band("cen", 14000.0, 16000.0, "acen")
    .add_band("q21", 20000.0, 24000.0, "gpos25")
  println("   Added regions: 2")
  println("   Added bands: 4")
  println("")

  println("5. Feature Query Operations")
  let all_features = chr_full.features
  let features_in_region = chr_full.features_in_region(5000.0, 12000.0)
  println("   All features: " + all_features.length().to_string())
  println("   Features in region 5000-12000: " + features_in_region.length().to_string())
  println("")

  println("6. Band Color Calculation")
  let stains = ["gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"]
  for stain in stains {
    let color = @bio.chr_band_color(stain)
    println("   Stain: " + stain + " -> Color: " + color)
  }
  println("")

  println("7. Creating Layout Configuration")
  let layout = @bio.ChrLayout::new(
    width=900.0,
    height=500.0,
    show_labels=true,
    show_centromere=true,
    show_bands=true,
    chromosome_thickness=12.0,
  )
  println("   Layout: 900 x 500")
  println("   Show labels: true")
  println("   Show centromere: true")
  println("   Show bands: true")
  println("")

  println("8. Generating SVG Output")
  let svg_elements = @bio.chr_chromosome_to_svg(chr_full, layout, 50.0, 100.0, 0.02)
  println("   SVG elements generated: " + svg_elements.length().to_string())
  println("   First element type: " + svg_elements[0].element_type)
  println("")

  println("9. Creating a Bacterial Chromosome (Circular)")
  let bacterial_chr = @bio.Chromosome::new(
    name="E.coli",
    length=4600000.0,
    centromere_pos=0.0,
    is_circular=true,
  )
  let bacterial_svg = @bio.chr_chromosome_to_svg(bacterial_chr, layout, 50.0, 100.0, 0.001)
  println("   Circular chromosome SVG elements: " + bacterial_svg.length().to_string())
  println("")

  println("=== Demo completed successfully! ===")
  println("")
  println("Key concepts:")
  println("- Chromosome visualization uses SVG for rendering")
  println("- Features represent genes, markers, and annotations")
  println("- Bands represent staining patterns (G-banding)")
  println("- Linear chromosomes have centromere positions")
  println("- Circular chromosomes are rendered differently")
}
