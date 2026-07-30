///|
/// Test file for Gviz module.

test "gviz_feature_creation" {
  let f = @src.gviz_feature("g1", "chr1", 100, 200, @src.track_strand_forward(), "exon", "GeneA")
  assert_eq(f.feature_id, "g1")
  assert_eq(f.chromosome, "chr1")
  assert_eq(f.start, 100)
  assert_eq(f.end_, 200)
  assert_eq(f.feature_type, "exon")
  assert_eq(f.label, "GeneA")
}

test "gviz_track_creation" {
  let t = @src.gviz_track("genes", @src.track_type_gene_region(), "chr1", 1000, 5000)
  assert_eq(t.track_name, "genes")
  assert_eq(t.chromosome, "chr1")
  assert_eq(t.start, 1000)
  assert_eq(t.end_, 5000)
  assert_eq(t.get_n_features(), 0)
}

test "gviz_track_add_feature" {
  let t = @src.gviz_track("genes", @src.track_type_gene_region(), "chr1", 1000, 5000)
  t.add_feature(@src.gviz_feature("g1", "chr1", 1200, 1800, @src.track_strand_forward(), "exon", "GeneA"))
  t.add_feature(@src.gviz_feature("g2", "chr1", 2000, 3500, @src.track_strand_reverse(), "exon", "GeneB"))
  assert_eq(t.get_n_features(), 2)
}

test "gviz_track_add_data_point" {
  let t = @src.gviz_track("coverage", @src.track_type_data(), "chr1", 1000, 5000)
  t.add_data_point(1000, 5.0)
  t.add_data_point(2000, 10.0)
  t.add_data_point(3000, 15.0)
  assert_eq(t.get_n_data_points(), 3)
}

test "gviz_track_set_color" {
  let t = @src.gviz_track("track1", @src.track_type_annotation(), "chr1", 1000, 5000)
  t.set_color("red")
  assert_eq(t.color, "red")
}

test "gviz_track_set_label" {
  let t = @src.gviz_track("track1", @src.track_type_annotation(), "chr1", 1000, 5000)
  t.set_label("My Track")
  assert_eq(t.display_label, "My Track")
}

test "gviz_track_get_features_in_region" {
  let t = @src.gviz_track("genes", @src.track_type_gene_region(), "chr1", 1000, 5000)
  t.add_feature(@src.gviz_feature("g1", "chr1", 1200, 1800, @src.track_strand_forward(), "exon", "A"))
  t.add_feature(@src.gviz_feature("g2", "chr1", 2000, 3500, @src.track_strand_reverse(), "exon", "B"))
  t.add_feature(@src.gviz_feature("g3", "chr1", 4000, 4500, @src.track_strand_forward(), "exon", "C"))

  let in_region = t.get_features_in_region(1500, 3000)
  assert_eq(in_region.length(), 2) // g1 (overlaps) and g2
}

test "gviz_track_get_data_in_region" {
  let t = @src.gviz_track("data", @src.track_type_data(), "chr1", 1000, 5000)
  t.add_data_point(1000, 5.0)
  t.add_data_point(2000, 10.0)
  t.add_data_point(3000, 15.0)
  t.add_data_point(4000, 20.0)

  let data = t.get_data_in_region(1500, 3500)
  assert_eq(data.length(), 2) // 2000 and 3000
}

test "gviz_region_creation" {
  let r = @src.gviz_region("chr1", 1000, 5000)
  assert_eq(r.chromosome, "chr1")
  assert_eq(r.start, 1000)
  assert_eq(r.end_, 5000)
}

test "gviz_plot_creation" {
  let r = @src.gviz_region("chr1", 1000, 5000)
  let p = @src.gviz_plot(r, title="Test Plot", width=80, height=25)
  assert_eq(p.get_n_tracks(), 0)
  assert_eq(p.get_region().chromosome, "chr1")
}

test "gviz_plot_add_track" {
  let r = @src.gviz_region("chr1", 1000, 5000)
  let p = @src.gviz_plot(r, title="Test")
  let t = @src.gviz_track("track1", @src.track_type_annotation(), "chr1", 1000, 5000)
  p.add_track(t)
  assert_eq(p.get_n_tracks(), 1)
}

test "gviz_plot_to_ascii" {
  let p = @src.gviz_sample_plot()
  let ascii = p.to_ascii()
  assert_true(ascii.contains("chr1"))
  assert_true(ascii.contains("GeneRegionTrack"))
}

test "gviz_plot_summary" {
  let p = @src.gviz_sample_plot()
  let s = p.summary()
  assert_true(s.contains("Gviz"))
  assert_true(s.contains("chr1"))
}

test "gviz_track_type_to_string" {
  assert_eq(@src.track_type_annotation().to_string(), "AnnotationTrack")
  assert_eq(@src.track_type_gene_region().to_string(), "GeneRegionTrack")
  assert_eq(@src.track_type_data().to_string(), "DataTrack")
  assert_eq(@src.gviz_track_type_ideogram().to_string(), "IdeogramTrack")
  assert_eq(@src.track_type_genome_axis().to_string(), "GenomeAxisTrack")
  assert_eq(@src.track_type_sequence().to_string(), "SequenceTrack")
}

test "gviz_strand_to_string" {
  assert_eq(@src.track_strand_forward().to_string(), "+")
  assert_eq(@src.track_strand_reverse().to_string(), "-")
  assert_eq(@src.track_strand_unstranded().to_string(), "*")
}

test "gviz_sample_plot" {
  let p = @src.gviz_sample_plot()
  assert_true(p.get_n_tracks() >= 3)
}
