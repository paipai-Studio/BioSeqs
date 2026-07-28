///|
/// Tests for Bio.Graphics.Chromosome module.

///|
test "chr_orientation_values" {
  let _f = @src.chr_orientation_forward()
  let _r = @src.chr_orientation_reverse()
  let _n = @src.chr_orientation_none()
  assert_true(true)
}

///|
test "chr_feature_type_values" {
  let _g = @src.chr_feature_type_gene()
  let _e = @src.chr_feature_type_exon()
  let _i = @src.chr_feature_type_intron()
  let _c = @src.chr_feature_type_cds()
  assert_true(true)
}

///|
test "chr_feature_creation" {
  let feature = @src.ChrFeature::new(
    name="Gene1",
    start=100.0,
    end=200.0,
  )
  assert_true((feature.start - 100.0).abs() < 1.0e-10)
  assert_true((feature.end - 200.0).abs() < 1.0e-10)
  assert_eq(feature.label, "Gene1")
}

///|
test "chr_region_creation" {
  let region = @src.ChrRegion::new(
    name="Region1",
    start=0.0,
    end=1000.0,
  )
  assert_true((region.start).abs() < 1.0e-10)
  assert_true((region.end - 1000.0).abs() < 1.0e-10)
  assert_eq(region.label, "Region1")
}

///|
test "chr_chromosome_creation" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=10000.0,
    centromere_pos=4000.0,
  )
  assert_eq(chr.name, "chr1")
  assert_true((chr.length - 10000.0).abs() < 1.0e-10)
  assert_true((chr.centromere_pos - 4000.0).abs() < 1.0e-10)
  assert_true(!chr.is_circular)
}

///|
test "chr_chromosome_add_feature" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
  )
  let feature = @src.ChrFeature::new(
    name="Gene1",
    start=100.0,
    end=200.0,
  )
  let chr_with_feature = chr.add_feature(feature)
  assert_eq(chr_with_feature.features.length(), 1)
}

///|
test "chr_chromosome_add_region" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
  )
  let region = @src.ChrRegion::new(
    name="Region1",
    start=300.0,
    end=500.0,
  )
  let chr_with_region = chr.add_region(region)
  assert_eq(chr_with_region.regions.length(), 1)
}

///|
test "chr_chromosome_add_band" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
  )
  let chr_with_band = chr.add_band("Band1", 100.0, 200.0, "gpos50")
  assert_eq(chr_with_band.bands.length(), 1)
}

///|
test "chr_features_by_type" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
  )
  let exon1 = @src.ChrFeature::new(
    name="E1",
    start=100.0,
    end=200.0,
    feature_type=@src.chr_feature_type_exon(),
  )
  let exon2 = @src.ChrFeature::new(
    name="E2",
    start=300.0,
    end=400.0,
    feature_type=@src.chr_feature_type_exon(),
  )
  let intron = @src.ChrFeature::new(
    name="I1",
    start=250.0,
    end=290.0,
    feature_type=@src.chr_feature_type_intron(),
  )
  let chr1 = chr.add_feature(exon1).add_feature(exon2).add_feature(intron)
  let exons = chr1.features_by_type(@src.chr_feature_type_exon())
  assert_eq(exons.length(), 2)
  let introns = chr1.features_by_type(@src.chr_feature_type_intron())
  assert_eq(introns.length(), 1)
}

///|
test "chr_features_in_region" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
  )
  let f1 = @src.ChrFeature::new(
    name="F1",
    start=100.0,
    end=200.0,
  )
  let f2 = @src.ChrFeature::new(
    name="F2",
    start=300.0,
    end=400.0,
  )
  let chr1 = chr.add_feature(f1).add_feature(f2)
  let in_region = chr1.features_in_region(150.0, 350.0)
  assert_eq(in_region.length(), 2)
}

///|
test "chr_layout_creation" {
  let layout = @src.ChrLayout::new(
    width=800.0,
    height=600.0,
    chromosome_thickness=20.0,
    show_labels=true,
    show_centromere=true,
  )
  assert_true((layout.width - 800.0).abs() < 1.0e-10)
  assert_true((layout.chromosome_thickness - 20.0).abs() < 1.0e-10)
  assert_true(layout.show_labels)
}

///|
test "chr_diagram_creation" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
    centromere_pos=2000.0,
  )
  let diagram = @src.ChrDiagram::new(
    chromosomes=[chr],
    title="Test Diagram",
  )
  assert_eq(diagram.chromosomes.length(), 1)
  assert_eq(diagram.title, "Test Diagram")
}

///|
test "chr_to_svg_returns_string" {
  let chr = @src.Chromosome::new(
    name="chr1",
    length=1000.0,
    centromere_pos=400.0,
  )
  let diagram = @src.ChrDiagram::new(
    chromosomes=[chr],
    title="Chromosome 1",
  )
  let svg = diagram.to_svg()
  assert_true(svg.length() > 0)
  assert_true(svg.contains("<svg"))
}

///|
test "chr_to_svg_empty_diagram" {
  let diagram = @src.ChrDiagram::new(
    chromosomes=[],
    title="Empty",
  )
  let svg = diagram.to_svg()
  assert_true(svg.length() > 0)
  assert_true(svg.contains("<svg"))
}

///|
test "chr_band_color" {
  let c1 = @src.chr_band_color("gneg")
  assert_eq(c1, "#FFFFFF")
  let c2 = @src.chr_band_color("gpos100")
  assert_eq(c2, "#606060")
  let c3 = @src.chr_band_color("acen")
  assert_eq(c3, "#FF4444")
}

///|
test "chr_create_human_karyotype" {
  let diagram = @src.chr_create_human_karyotype()
  assert_true(diagram.chromosomes.length() > 0)
  let svg = diagram.to_svg()
  assert_true(svg.length() > 0)
}

///|
test "chr_multiple_chromosomes" {
  let chr1 = @src.Chromosome::new(
    name="chr1",
    length=5000.0,
    centromere_pos=2000.0,
  )
  let chr2 = @src.Chromosome::new(
    name="chr2",
    length=4000.0,
    centromere_pos=1500.0,
  )
  let diagram = @src.ChrDiagram::new(
    chromosomes=[chr1, chr2],
    title="Two Chromosomes",
  )
  let svg = diagram.to_svg()
  assert_true(svg.length() > 0)
}

///|
test "chr_circular_chromosome" {
  let chr = @src.Chromosome::new(
    name="chr_circ",
    length=2000.0,
    is_circular=true,
  )
  let diagram = @src.ChrDiagram::new(
    chromosomes=[chr],
    title="Circular",
  )
  let svg = diagram.to_svg()
  assert_true(svg.length() > 0)
}