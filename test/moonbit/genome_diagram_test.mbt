///|
/// Tests for GenomeDiagram module.

test "gd_create_diagram - basic creation" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  assert_eq(d.name(), "test")
  assert_eq(d.start(), 0)
  assert_eq(d.end(), 1000)
  assert_eq(@src.gd_diagram_length(d), 1000)
}

test "gd_create_diagram - negative coordinates" {
  let d = @src.gd_create_diagram("neg", -100, 200)
  assert_eq(d.start(), -100)
  assert_eq(d.end(), 200)
  assert_eq(@src.gd_diagram_length(d), 300)
}

test "gd_add_track - single track" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let d2 = @src.gd_add_track(d, "genes")
  assert_eq(@src.gd_track_count(d2), 1)
  let tracks = d2.tracks()
  assert_eq(tracks[0].name(), "genes")
}

test "gd_add_track - multiple tracks" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_track(d2, "exons")
  let d4 = @src.gd_add_track(d3, "regulatory")
  assert_eq(@src.gd_track_count(d4), 3)
}

test "gd_add_feature - basic feature" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 300, "+", "GeneA", @src.gd_rectangle_shape())
  assert_eq(@src.gd_feature_count(d2), 1)
  let feats = d2.features()
  assert_eq(feats[0].start(), 100)
  assert_eq(feats[0].end(), 300)
  assert_eq(feats[0].strand(), "+")
  assert_eq(feats[0].label(), "GeneA")
}

test "gd_add_feature - multiple features" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 50, 150, "+", "Gene1", @src.gd_arrow_shape())
  let d3 = @src.gd_add_feature(d2, 200, 400, "-", "Gene2", @src.gd_diamond_shape())
  let d4 = @src.gd_add_feature(d3, 500, 800, "+", "Gene3", @src.gd_rectangle_shape())
  assert_eq(@src.gd_feature_count(d4), 3)
}

test "gd_add_track_feature - basic track feature" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_track_feature(d2, "genes", 100, 300, "+", "GeneA")
  let feats = @src.gd_get_features(d3, "genes")
  assert_eq(feats.length(), 1)
  assert_eq(feats[0].start(), 100)
  assert_eq(feats[0].end(), 300)
  assert_eq(feats[0].label(), "GeneA")
}

test "gd_add_track_feature - multiple features in same track" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_track_feature(d2, "genes", 100, 200, "+", "G1")
  let d4 = @src.gd_add_track_feature(d3, "genes", 300, 400, "+", "G2")
  let d5 = @src.gd_add_track_feature(d4, "genes", 500, 600, "-", "G3")
  let feats = @src.gd_get_features(d5, "genes")
  assert_eq(feats.length(), 3)
}

test "gd_add_track_feature - multiple tracks" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_track(d2, "exons")
  let d4 = @src.gd_add_track_feature(d3, "genes", 100, 300, "+", "GeneA")
  let d5 = @src.gd_add_track_feature(d4, "exons", 150, 200, "+", "Exon1")
  let d6 = @src.gd_add_track_feature(d5, "genes", 500, 700, "-", "GeneB")
  assert_eq(@src.gd_feature_count(d6), 3)
  let gene_feats = @src.gd_get_features(d6, "genes")
  assert_eq(gene_feats.length(), 2)
  let exon_feats = @src.gd_get_features(d6, "exons")
  assert_eq(exon_feats.length(), 1)
}

test "gd_get_features - nonexistent track returns empty" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let feats = @src.gd_get_features(d, "nonexistent")
  assert_eq(feats.length(), 0)
}

test "gd_find_overlapping_features - basic overlap" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 300, "+", "GeneA", @src.gd_rectangle_shape())
  let d3 = @src.gd_add_feature(d2, 400, 600, "+", "GeneB", @src.gd_rectangle_shape())
  let d4 = @src.gd_add_feature(d3, 250, 500, "+", "GeneC", @src.gd_rectangle_shape())
  let overlaps = @src.gd_find_overlapping_features(d4, 200, 400)
  assert_eq(overlaps.length(), 2)
}

test "gd_find_overlapping_features - no overlap" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 200, "+", "GeneA", @src.gd_rectangle_shape())
  let d3 = @src.gd_add_feature(d2, 500, 600, "+", "GeneB", @src.gd_rectangle_shape())
  let overlaps = @src.gd_find_overlapping_features(d3, 300, 400)
  assert_eq(overlaps.length(), 0)
}

test "gd_find_overlapping_features - edge touching" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 200, "+", "GeneA", @src.gd_rectangle_shape())
  let overlaps = @src.gd_find_overlapping_features(d2, 200, 300)
  assert_eq(overlaps.length(), 0)
}

test "gd_set_style - basic style" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let style = @src.DiagramStyle::new()
  let d2 = @src.gd_set_style(d, style)
  assert_eq(d2.style().scale(), 1.0)
  assert_true(d2.style().linear())
  assert_true(d2.style().border())
}

test "gd_set_style - circular mode" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let style = @src.DiagramStyle::new(circular=true, linear=false)
  let d2 = @src.gd_set_style(d, style)
  assert_true(d2.style().circular())
  assert_false(d2.style().linear())
}

test "gd_diagram_length - different ranges" {
  let d1 = @src.gd_create_diagram("a", 0, 500)
  assert_eq(@src.gd_diagram_length(d1), 500)
  let d2 = @src.gd_create_diagram("b", 100, 200)
  assert_eq(@src.gd_diagram_length(d2), 100)
  let d3 = @src.gd_create_diagram("c", -50, 50)
  assert_eq(@src.gd_diagram_length(d3), 100)
}

test "gd_track_count - empty diagram" {
  let d = @src.gd_create_diagram("test", 0, 100)
  assert_eq(@src.gd_track_count(d), 0)
}

test "gd_feature_count - mixed features" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_feature(d2, 50, 150, "+", "G1", @src.gd_rectangle_shape())
  let d4 = @src.gd_add_track_feature(d3, "genes", 100, 200, "+", "G2")
  let d5 = @src.gd_add_feature(d4, 300, 400, "-", "G3", @src.gd_arrow_shape())
  assert_eq(@src.gd_feature_count(d5), 3)
}

test "gd_to_svg_string - generates valid SVG" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let d2 = @src.gd_add_track(d, "genes")
  let d3 = @src.gd_add_feature(d2, 50, 150, "+", "GeneA", @src.gd_rectangle_shape())
  let d4 = @src.gd_add_track_feature(d3, "genes", 200, 300, "+", "Exon1")
  let svg = @src.gd_to_svg_string(d4, 800, 400)
  assert_true(svg.length() > 0)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
}

test "gd_to_svg_string - empty diagram" {
  let d = @src.gd_create_diagram("empty", 0, 0)
  let svg = @src.gd_to_svg_string(d, 400, 200)
  assert_eq(svg, "")
}

test "gd_to_svg_string - contains feature labels" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let d2 = @src.gd_add_feature(d, 100, 200, "+", "MyGene", @src.gd_rectangle_shape())
  let svg = @src.gd_to_svg_string(d2, 800, 400)
  assert_true(svg.contains("MyGene"))
}

test "gd_set_feature_color - changes color" {
  let feat = @src.DiagramFeature::new(0, 100, color="#FF0000")
  let updated = @src.gd_set_feature_color(feat, "#00FF00")
  assert_eq(updated.color(), "#00FF00")
  assert_eq(feat.color(), "#FF0000")
}

test "gd_set_feature_color - shape preserved" {
  let feat = @src.DiagramFeature::new(0, 100, shape=@src.gd_arrow_shape(), color="#FF0000")
  let updated = @src.gd_set_feature_color(feat, "#00FF00")
  assert_eq(updated.shape().to_string(), "arrow")
}

test "gd_label_features - auto labels large features" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 300, "+", "", @src.gd_rectangle_shape())
  let d3 = @src.gd_add_feature(d2, 400, 450, "+", "", @src.gd_rectangle_shape())
  let d4 = @src.gd_label_features(d3, 50)
  let feats = d4.features()
  assert_eq(feats[0].label(), "feature_100")
  assert_eq(feats[1].label(), "")
}

test "gd_label_features - existing labels preserved" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 100, 300, "+", "Existing", @src.gd_rectangle_shape())
  let d3 = @src.gd_label_features(d2, 50)
  let feats = d3.features()
  assert_eq(feats[0].label(), "Existing")
}

test "FeatureShape - to_string conversions" {
  let r = @src.gd_rectangle_shape()
  assert_eq(r.to_string(), "rectangle")
  let a = @src.gd_arrow_shape()
  assert_eq(a.to_string(), "arrow")
  let d = @src.gd_diamond_shape()
  assert_eq(d.to_string(), "diamond")
  let c = @src.gd_crossed_arrow_shape()
  assert_eq(c.to_string(), "crossed_arrow")
  let t = @src.gd_terminators_shape()
  assert_eq(t.to_string(), "terminators")
}

test "FeatureShape - from_string conversions" {
  let r = @src.FeatureShape::from_string("rectangle")
  assert_eq(r.to_string(), "rectangle")
  let a = @src.FeatureShape::from_string("arrow")
  assert_eq(a.to_string(), "arrow")
  let d = @src.FeatureShape::from_string("diamond")
  assert_eq(d.to_string(), "diamond")
  let u = @src.FeatureShape::from_string("unknown")
  assert_eq(u.to_string(), "rectangle")
}

test "DiagramStyle - default values" {
  let style = @src.DiagramStyle::new()
  assert_eq(style.scale(), 1.0)
  assert_true(style.linear())
  assert_false(style.circular())
  assert_true(style.border())
  assert_eq(style.color_scheme(), "default")
}

test "DiagramStyle - setters" {
  let style = @src.DiagramStyle::new()
  let s2 = style.set_scale(2.5)
  assert_eq(s2.scale(), 2.5)
  let s3 = s2.set_border(false)
  assert_false(s3.border())
  let s4 = s3.set_color_scheme("grayscale")
  assert_eq(s4.color_scheme(), "grayscale")
}

test "DiagramStyle - circular and linear are mutually exclusive" {
  let style = @src.DiagramStyle::new()
  let circ = style.set_circular(true)
  assert_true(circ.circular())
  assert_false(circ.linear())
  let lin = circ.set_linear(true)
  assert_false(lin.circular())
  assert_true(lin.linear())
}

test "DiagramFeature - overlaps detection" {
  let feat = @src.DiagramFeature::new(100, 300)
  assert_true(feat.overlaps(50, 150))
  assert_true(feat.overlaps(250, 350))
  assert_true(feat.overlaps(150, 250))
  assert_false(feat.overlaps(300, 400))
  assert_false(feat.overlaps(0, 100))
}

test "DiagramFeature - length" {
  let feat = @src.DiagramFeature::new(100, 300)
  assert_eq(feat.length(), 200)
}

test "TrackFeature - length" {
  let feat = @src.TrackFeature::new(50, 150)
  assert_eq(feat.length(), 100)
}

test "Track - feature operations" {
  let t = @src.Track::new("test_track")
  assert_eq(t.name(), "test_track")
  assert_eq(t.feature_count(), 0)
  let f1 = @src.TrackFeature::new(100, 200, label="F1")
  let t2 = t.add_feature(f1)
  assert_eq(t2.feature_count(), 1)
  let f2 = @src.TrackFeature::new(300, 400, label="F2")
  let t3 = t2.add_feature(f2)
  assert_eq(t3.feature_count(), 2)
}

test "gd_feature_count - diagram features only" {
  let d = @src.gd_create_diagram("test", 0, 1000)
  let d2 = @src.gd_add_feature(d, 50, 100, "+", "F1", @src.gd_rectangle_shape())
  let d3 = @src.gd_add_feature(d2, 200, 300, "+", "F2", @src.gd_rectangle_shape())
  assert_eq(@src.gd_feature_count(d3), 2)
}

test "gd_to_svg_string - contains track names" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let d2 = @src.gd_add_track(d, "MyTrack")
  let d3 = @src.gd_add_track_feature(d2, "MyTrack", 100, 200, "+", "F1")
  let svg = @src.gd_to_svg_string(d3, 600, 300)
  assert_true(svg.contains("MyTrack"))
}

test "gd_to_svg_string - border toggle" {
  let d = @src.gd_create_diagram("test", 0, 500)
  let style = @src.DiagramStyle::new(border=false)
  let d2 = @src.gd_set_style(d, style)
  let svg = @src.gd_to_svg_string(d2, 600, 300)
  assert_true(svg.contains("<svg"))
}

test "DiagramFeature - default values" {
  let feat = @src.DiagramFeature::new(0, 100)
  assert_eq(feat.strand(), "+")
  assert_eq(feat.label(), "")
  assert_eq(feat.color(), "#45B7D1")
  assert_eq(feat.shape().to_string(), "rectangle")
  assert_eq(feat.depth(), 0)
}

test "TrackFeature - default values" {
  let feat = @src.TrackFeature::new(0, 100)
  assert_eq(feat.strand(), "+")
  assert_eq(feat.label(), "")
  assert_eq(feat.color(), "#FF6B6B")
  assert_eq(feat.shape().to_string(), "rectangle")
}