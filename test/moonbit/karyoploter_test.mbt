///|
/// Test file for karyoploteR module.

test "karyotype_new" {
  let plot = @src.KaryotypePlot::new("hg38")
  assert_eq(plot.get_genome(), "hg38")
  assert_true(plot.get_chromosomes().length() > 0)
}

test "karyotype_chromosomes" {
  let plot = @src.KaryotypePlot::new("hg38")
  let chroms = plot.get_chromosomes()
  assert_true(chroms.contains("chr1"))
  assert_true(chroms.contains("chrX"))
  assert_true(chroms.contains("chrY"))
}

test "karyotype_chromosome_size" {
  let plot = @src.KaryotypePlot::new("hg38")
  let size = plot.get_chromosome_size("chr1")
  assert_eq(size, 248956422)
}

test "karyotype_unknown_chromosome" {
  let plot = @src.KaryotypePlot::new("hg38")
  let size = plot.get_chromosome_size("chrUnknown")
  assert_eq(size, 0)
}

test "karyotype_add_track" {
  let plot = @src.KaryotypePlot::new("hg38")
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 248956422, "Test Track")
  let plot2 = plot.add_track(track)
  assert_eq(plot2.get_n_tracks(), 1)
}

test "karyotype_add_region" {
  let plot = @src.KaryotypePlot::new("hg38")
  let region = @src.KaryotypeRegion::new("chr1", 1000000, 2000000, "Region 1", "#FF0000")
  let plot2 = plot.add_region(region)
  assert_eq(plot2.get_n_regions(), 1)
}

test "track_new" {
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 1000000, "My Track")
  assert_eq(track.track_id, "t1")
  assert_eq(track.chromosome, "chr1")
  assert_eq(track.start, 0)
  assert_eq(track.end, 1000000)
}

test "track_add_point" {
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 1000000, "Test")
  let point = @src.TrackPoint::new("chr1", 500000, 0.8, "peak")
  let track2 = track.add_point(point)
  assert_eq(track2.get_n_points(), 1)
  let data = track2.get_data()
  assert_eq(data[0].position, 500000)
}

test "track_set_color" {
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 1000000, "Test")
  let track2 = track.set_color("#FF0000")
  assert_eq(track2.color, "#FF0000")
}

test "track_set_y_range" {
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 1000000, "Test")
  let track2 = track.set_y_range(0.0, 2.0)
  assert_eq(track2.y_max, 2.0)
}

test "track_filter_chromosome" {
  let track = @src.KaryotypeTrack::new("t1", @src.track_type_points(), "chr1", 0, 248956422, "Test")
  let p1 = @src.TrackPoint::new("chr1", 1000000, 0.5, "")
  let p2 = @src.TrackPoint::new("chr2", 2000000, 0.8, "")
  let track2 = track.add_point(p1)
  let track3 = track2.add_point(p2)
  let filtered = track3.filter_chromosome("chr1")
  assert_eq(filtered.get_n_points(), 1)
}

test "track_point_new" {
  let pt = @src.TrackPoint::new("chr1", 500000, 0.75, "test_point")
  assert_eq(pt.chromosome, "chr1")
  assert_eq(pt.position, 500000)
  assert_eq(pt.value, 0.75)
  assert_eq(pt.label, "test_point")
}

test "karyotype_region_new" {
  let region = @src.KaryotypeRegion::new("chr1", 1000, 2000, "Region", "#FF0000")
  assert_eq(region.chromosome, "chr1")
  assert_eq(region.start, 1000)
  assert_eq(region.end, 2000)
}

test "ideogram_band_new" {
  let band = @src.IdeogramBand::new("chr1", 0, 5000000, "p36.3", "gneg")
  assert_eq(band.chromosome, "chr1")
  assert_eq(band.band_name, "p36.3")
}

test "track_type_to_string" {
  assert_eq(@src.track_type_points().to_string(), "points")
  assert_eq(@src.track_type_lines().to_string(), "lines")
  assert_eq(@src.track_type_bars().to_string(), "bars")
  assert_eq(@src.track_type_heatmap().to_string(), "heatmap")
  assert_eq(@src.track_type_ideogram().to_string(), "ideogram")
}

test "karyotype_to_ascii" {
  let plot = @src.karyotype_sample()
  let ascii = plot.to_ascii("chr1", 60)
  assert_true(ascii.length() > 0)
  assert_true(ascii.contains("hg38"))
  assert_true(ascii.contains("chr1"))
}

test "karyotype_summary" {
  let plot = @src.karyotype_sample()
  let summary = plot.summary()
  assert_true(summary.contains("Karyotype Plot Summary"))
  assert_true(summary.contains("hg38"))
  assert_true(summary.contains("chr1"))
}

test "karyotype_unknown_genome" {
  let plot = @src.KaryotypePlot::new("mm10")
  assert_eq(plot.get_genome(), "mm10")
  assert_eq(plot.get_chromosomes().length(), 0)
}

test "karyotype_sample" {
  let plot = @src.karyotype_sample()
  assert_true(plot.get_n_tracks() >= 2)
  assert_true(plot.get_n_regions() >= 2)
  assert_eq(plot.get_genome(), "hg38")
}

test "karyotype_ascii_unknown_chr" {
  let plot = @src.karyotype_sample()
  let ascii = plot.to_ascii("chr99", 60)
  assert_true(ascii.contains("not found"))
}