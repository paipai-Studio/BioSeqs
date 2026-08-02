///|
/// Tests for the methylSeekR module.
/// Mirrors the Bioconductor methylSeekR package API.

///|
test "msr_cytosine_site_creation" {
  let s = @src.CytosineSite::new(
    chr="chr1",
    position=1000,
    methylated=8,
    unmethylated=2,
  )
  assert_eq(@src.CytosineSite::chr(s), "chr1")
  assert_eq(@src.CytosineSite::position(s), 1000)
  assert_eq(@src.CytosineSite::methylated(s), 8)
  assert_eq(@src.CytosineSite::unmethylated(s), 2)
}

///|
test "msr_cytosine_site_coverage" {
  let s = @src.CytosineSite::new(
    chr="chr1",
    position=100,
    methylated=8,
    unmethylated=2,
  )
  assert_eq(@src.CytosineSite::coverage(s), 10)
  let zero = @src.CytosineSite::new(
    chr="chr1",
    position=200,
    methylated=0,
    unmethylated=0,
  )
  assert_eq(@src.CytosineSite::coverage(zero), 0)
}

///|
test "msr_cytosine_site_methylation_level" {
  // 8 methylated, 2 unmethylated -> 0.8
  let s = @src.CytosineSite::new(
    chr="chr1",
    position=100,
    methylated=8,
    unmethylated=2,
  )
  let level = @src.methylation_level(s)
  assert_true((level - 0.8).abs() < 0.000001)
  // 5/5 -> 0.5
  let half = @src.CytosineSite::new(
    chr="chr1",
    position=200,
    methylated=5,
    unmethylated=5,
  )
  assert_true((@src.methylation_level(half) - 0.5).abs() < 0.000001)
  // 0/0 -> 0.0 (zero coverage)
  let zero = @src.CytosineSite::new(
    chr="chr1",
    position=300,
    methylated=0,
    unmethylated=0,
  )
  assert_eq(@src.methylation_level(zero), 0.0)
}

///|
test "msr_cytosine_site_to_string" {
  let s = @src.CytosineSite::new(
    chr="chr1",
    position=100,
    methylated=8,
    unmethylated=2,
  )
  let str = @src.CytosineSite::to_string(s)
  assert_true(str.contains("CytosineSite"))
  assert_true(str.contains("chr1"))
  assert_true(str.contains("coverage=10"))
}

///|
test "msr_methyl_tile_creation" {
  let t = @src.MethylTile::new(
    chr="chr1",
    start=0,
    end=1000,
    methylation_level=0.25,
    coverage=40,
    region_type="LMR",
  )
  assert_eq(@src.MethylTile::chr(t), "chr1")
  assert_eq(@src.MethylTile::start(t), 0)
  assert_eq(@src.MethylTile::end(t), 1000)
  assert_true((@src.MethylTile::methylation_level(t) - 0.25).abs() < 0.000001)
  assert_eq(@src.MethylTile::coverage(t), 40)
  assert_eq(@src.MethylTile::region_type(t), "LMR")
}

///|
test "msr_methyl_tile_default_region_type" {
  // region_type omitted -> defaults to "FMR".
  let t = @src.MethylTile::new(
    chr="chr1",
    start=0,
    end=1000,
    methylation_level=0.9,
    coverage=100,
  )
  assert_eq(@src.MethylTile::region_type(t), "FMR")
}

///|
test "msr_methyl_region_creation" {
  let r = @src.MethylRegion::new(
    chr="chr2",
    start=5000,
    end=15000,
    region_type="PMD",
    n_tiles=10,
    mean_methylation=0.62,
    min_methylation=0.55,
    max_methylation=0.68,
  )
  assert_eq(@src.MethylRegion::chr(r), "chr2")
  assert_eq(@src.MethylRegion::start(r), 5000)
  assert_eq(@src.MethylRegion::end(r), 15000)
  assert_eq(@src.MethylRegion::region_type(r), "PMD")
  assert_eq(@src.MethylRegion::n_tiles(r), 10)
  assert_true((@src.MethylRegion::mean_methylation(r) - 0.62).abs() < 0.000001)
  assert_true((@src.MethylRegion::min_methylation(r) - 0.55).abs() < 0.000001)
  assert_true((@src.MethylRegion::max_methylation(r) - 0.68).abs() < 0.000001)
}

///|
test "msr_methyl_region_size" {
  let r = @src.MethylRegion::new(
    chr="chr1",
    start=1000,
    end=6000,
    region_type="UMR",
    n_tiles=5,
    mean_methylation=0.05,
    min_methylation=0.04,
    max_methylation=0.06,
  )
  assert_eq(@src.MethylRegion::size(r), 5000)
}

///|
test "msr_params_default" {
  let p = @src.MethylSeekRParams::default()
  assert_eq(@src.MethylSeekRParams::tile_size(p), 1000)
  assert_eq(@src.MethylSeekRParams::min_coverage(p), 5)
  assert_true((@src.MethylSeekRParams::alpha(p) - 0.1).abs() < 0.000001)
  assert_true((@src.MethylSeekRParams::beta(p) - 0.5).abs() < 0.000001)
  assert_true((@src.MethylSeekRParams::gamma(p) - 0.7).abs() < 0.000001)
  assert_eq(@src.MethylSeekRParams::min_region_size(p), 500)
  assert_eq(@src.MethylSeekRParams::pmd_window(p), 10000)
}

///|
test "msr_params_new_all_defaults" {
  let p = @src.MethylSeekRParams::new()
  assert_eq(@src.MethylSeekRParams::tile_size(p), 1000)
  assert_eq(@src.MethylSeekRParams::min_coverage(p), 5)
  assert_eq(@src.MethylSeekRParams::pmd_window(p), 10000)
}

///|
test "msr_params_new_custom" {
  let p = @src.MethylSeekRParams::new(
    tile_size=500,
    min_coverage=10,
    alpha=0.05,
    beta=0.4,
    gamma=0.6,
    min_region_size=1000,
    pmd_window=5000,
  )
  assert_eq(@src.MethylSeekRParams::tile_size(p), 500)
  assert_eq(@src.MethylSeekRParams::min_coverage(p), 10)
  assert_true((@src.MethylSeekRParams::alpha(p) - 0.05).abs() < 0.000001)
  assert_true((@src.MethylSeekRParams::beta(p) - 0.4).abs() < 0.000001)
  assert_true((@src.MethylSeekRParams::gamma(p) - 0.6).abs() < 0.000001)
  assert_eq(@src.MethylSeekRParams::min_region_size(p), 1000)
  assert_eq(@src.MethylSeekRParams::pmd_window(p), 5000)
}

///|
test "msr_params_new_partial_override" {
  let p = @src.MethylSeekRParams::new(tile_size=250, gamma=0.8)
  assert_eq(@src.MethylSeekRParams::tile_size(p), 250)
  assert_true((@src.MethylSeekRParams::gamma(p) - 0.8).abs() < 0.000001)
  // Untouched parameters keep their defaults.
  assert_eq(@src.MethylSeekRParams::min_coverage(p), 5)
  assert_eq(@src.MethylSeekRParams::pmd_window(p), 10000)
}

///|
test "msr_tiling_empty" {
  let params = @src.MethylSeekRParams::default()
  let tiles = @src.tile_methylation([], params)
  assert_eq(tiles.length(), 0)
}

///|
test "msr_tiling_single_site" {
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(
      chr="chr1",
      position=250,
      methylated=4,
      unmethylated=6,
    ),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 1)
  assert_eq(@src.MethylTile::chr(tiles[0]), "chr1")
  assert_eq(@src.MethylTile::start(tiles[0]), 0)
  assert_eq(@src.MethylTile::end(tiles[0]), 1000)
  assert_eq(@src.MethylTile::coverage(tiles[0]), 10)
  assert_true((@src.MethylTile::methylation_level(tiles[0]) - 0.4).abs() < 0.000001)
  // Newly created tiles default to "FMR".
  assert_eq(@src.MethylTile::region_type(tiles[0]), "FMR")
}

///|
test "msr_tiling_boundaries" {
  // A site at position 999 falls in tile [0,1000); at 1000 in [1000,2000).
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr1", position=999, methylated=1, unmethylated=0),
    @src.CytosineSite::new(chr="chr1", position=1000, methylated=1, unmethylated=0),
    @src.CytosineSite::new(chr="chr1", position=2500, methylated=0, unmethylated=1),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 3)
  assert_eq(@src.MethylTile::start(tiles[0]), 0)
  assert_eq(@src.MethylTile::end(tiles[0]), 1000)
  assert_eq(@src.MethylTile::start(tiles[1]), 1000)
  assert_eq(@src.MethylTile::end(tiles[1]), 2000)
  assert_eq(@src.MethylTile::start(tiles[2]), 2000)
  assert_eq(@src.MethylTile::end(tiles[2]), 3000)
}

///|
test "msr_tiling_methylation_computation" {
  // Two sites in the same tile: (3,7) and (2,8) -> meth=5, unmeth=15 -> 0.25.
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr1", position=100, methylated=3, unmethylated=7),
    @src.CytosineSite::new(chr="chr1", position=200, methylated=2, unmethylated=8),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 1)
  assert_eq(@src.MethylTile::coverage(tiles[0]), 20)
  assert_true(
    (@src.MethylTile::methylation_level(tiles[0]) - 0.25).abs() < 0.000001,
  )
}

///|
test "msr_tiling_coverage_aggregation" {
  // Coverage from multiple sites is summed per tile.
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr1", position=100, methylated=1, unmethylated=4),
    @src.CytosineSite::new(chr="chr1", position=200, methylated=2, unmethylated=3),
    @src.CytosineSite::new(chr="chr1", position=300, methylated=0, unmethylated=5),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 1)
  // Total coverage = 5 + 5 + 5 = 15, methylated = 3 -> level 0.2.
  assert_eq(@src.MethylTile::coverage(tiles[0]), 15)
  assert_true(
    (@src.MethylTile::methylation_level(tiles[0]) - 0.2).abs() < 0.000001,
  )
}

///|
test "msr_tiling_multiple_chromosomes" {
  // Tiles are sorted by chromosome name then start.
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr2", position=50, methylated=1, unmethylated=1),
    @src.CytosineSite::new(chr="chr1", position=100, methylated=1, unmethylated=1),
    @src.CytosineSite::new(chr="chr1", position=1500, methylated=1, unmethylated=1),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 3)
  assert_eq(@src.MethylTile::chr(tiles[0]), "chr1")
  assert_eq(@src.MethylTile::start(tiles[0]), 0)
  assert_eq(@src.MethylTile::chr(tiles[1]), "chr1")
  assert_eq(@src.MethylTile::start(tiles[1]), 1000)
  assert_eq(@src.MethylTile::chr(tiles[2]), "chr2")
  assert_eq(@src.MethylTile::start(tiles[2]), 0)
}

///|
test "msr_classify_umr" {
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
    ),
  ]
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "UMR")
}

///|
test "msr_classify_lmr" {
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.3,
      coverage=100,
    ),
  ]
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "LMR")
}

///|
test "msr_classify_fmr" {
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.9,
      coverage=100,
    ),
  ]
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "FMR")
}

///|
test "msr_classify_low_coverage" {
  // Same low methylation, but coverage below min_coverage -> "FMR".
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.02,
      coverage=2,
    ),
  ]
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "FMR")
}

///|
test "msr_classify_threshold_boundaries" {
  // level == alpha (0.1) -> not UMR (strict <), becomes LMR.
  // level == beta (0.5) -> not LMR (strict <), becomes FMR.
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.1,
      coverage=100,
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=1000,
      end=2000,
      methylation_level=0.5,
      coverage=100,
    ),
  ]
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "LMR")
  assert_eq(@src.MethylTile::region_type(out[1]), "FMR")
}

///|
test "msr_pmd_detection_broad_domain" {
  // 15 contiguous FMR tiles at 60% methylation -> broad domain -> PMD.
  let params = @src.MethylSeekRParams::default()
  let tiles : Array[@src.MethylTile] = Array::new()
  for k in 0..<15 {
    tiles.push(
      @src.MethylTile::new(
        chr="chr1",
        start=k * 1000,
        end=(k + 1) * 1000,
        methylation_level=0.6,
        coverage=100,
        region_type="FMR",
      ),
    )
  }
  let out = @src.detect_pmds(tiles, params)
  for k in 0..<15 {
    assert_eq(@src.MethylTile::region_type(out[k]), "PMD")
  }
}

///|
test "msr_pmd_detection_short_run" {
  // 3 contiguous FMR tiles at 60% methylation -> too short (< 10) -> stay FMR.
  let params = @src.MethylSeekRParams::default()
  let tiles : Array[@src.MethylTile] = Array::new()
  for k in 0..<3 {
    tiles.push(
      @src.MethylTile::new(
        chr="chr1",
        start=k * 1000,
        end=(k + 1) * 1000,
        methylation_level=0.6,
        coverage=100,
        region_type="FMR",
      ),
    )
  }
  let out = @src.detect_pmds(tiles, params)
  for k in 0..<3 {
    assert_eq(@src.MethylTile::region_type(out[k]), "FMR")
  }
}

///|
test "msr_pmd_does_not_override_umr_lmr" {
  // UMR and LMR tiles must not be reclassified as PMD.
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=1000,
      end=2000,
      methylation_level=0.3,
      coverage=100,
      region_type="LMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=2000,
      end=3000,
      methylation_level=0.6,
      coverage=100,
      region_type="FMR",
    ),
  ]
  let out = @src.detect_pmds(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "UMR")
  assert_eq(@src.MethylTile::region_type(out[1]), "LMR")
  // The single FMR tile is below gamma but the run is too short -> stays FMR.
  assert_eq(@src.MethylTile::region_type(out[2]), "FMR")
}

///|
test "msr_pmd_high_methylation_not_candidate" {
  // FMR tiles at 90% methylation (>= gamma) are not PMD candidates.
  let params = @src.MethylSeekRParams::default()
  let tiles : Array[@src.MethylTile] = Array::new()
  for k in 0..<15 {
    tiles.push(
      @src.MethylTile::new(
        chr="chr1",
        start=k * 1000,
        end=(k + 1) * 1000,
        methylation_level=0.9,
        coverage=100,
        region_type="FMR",
      ),
    )
  }
  let out = @src.detect_pmds(tiles, params)
  for k in 0..<15 {
    assert_eq(@src.MethylTile::region_type(out[k]), "FMR")
  }
}

///|
test "msr_merge_adjacent_same_class" {
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=1000,
      end=2000,
      methylation_level=0.06,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=2000,
      end=3000,
      methylation_level=0.04,
      coverage=100,
      region_type="UMR",
    ),
  ]
  let regions = @src.merge_regions(tiles, params)
  assert_eq(regions.length(), 1)
  assert_eq(@src.MethylRegion::start(regions[0]), 0)
  assert_eq(@src.MethylRegion::end(regions[0]), 3000)
  assert_eq(@src.MethylRegion::n_tiles(regions[0]), 3)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  // Mean of (0.05, 0.06, 0.04) = 0.05.
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[0]) - 0.05).abs() < 0.000001,
  )
  assert_true(
    (@src.MethylRegion::min_methylation(regions[0]) - 0.04).abs() < 0.000001,
  )
  assert_true(
    (@src.MethylRegion::max_methylation(regions[0]) - 0.06).abs() < 0.000001,
  )
}

///|
test "msr_merge_different_class_boundary" {
  // Adjacent tiles of different classes do not merge.
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=1000,
      end=2000,
      methylation_level=0.3,
      coverage=100,
      region_type="LMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=2000,
      end=3000,
      methylation_level=0.9,
      coverage=100,
      region_type="FMR",
    ),
  ]
  let regions = @src.merge_regions(tiles, params)
  assert_eq(regions.length(), 3)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  assert_eq(@src.MethylRegion::region_type(regions[1]), "LMR")
  assert_eq(@src.MethylRegion::region_type(regions[2]), "FMR")
  // Each region is a single tile.
  assert_eq(@src.MethylRegion::n_tiles(regions[0]), 1)
  assert_eq(@src.MethylRegion::n_tiles(regions[1]), 1)
  assert_eq(@src.MethylRegion::n_tiles(regions[2]), 1)
}

///|
test "msr_merge_size_filter" {
  // A single 1000 bp tile is smaller than min_region_size=2000 -> filtered.
  let params = @src.MethylSeekRParams::new(min_region_size=2000)
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
  ]
  let regions = @src.merge_regions(tiles, params)
  assert_eq(regions.length(), 0)
  // With a smaller threshold the region is retained.
  let params2 = @src.MethylSeekRParams::new(min_region_size=500)
  let regions2 = @src.merge_regions(tiles, params2)
  assert_eq(regions2.length(), 1)
}

///|
test "msr_merge_different_chromosomes" {
  // Touching coordinates but different chromosomes do not merge.
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr2",
      start=1000,
      end=2000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
  ]
  let regions = @src.merge_regions(tiles, params)
  assert_eq(regions.length(), 2)
  assert_eq(@src.MethylRegion::chr(regions[0]), "chr1")
  assert_eq(@src.MethylRegion::chr(regions[1]), "chr2")
}

///|
test "msr_merge_gap_breaks_region" {
  // A gap between two same-class tiles prevents merging.
  let params = @src.MethylSeekRParams::default()
  let tiles = [
    @src.MethylTile::new(
      chr="chr1",
      start=0,
      end=1000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
    @src.MethylTile::new(
      chr="chr1",
      start=3000,
      end=4000,
      methylation_level=0.05,
      coverage=100,
      region_type="UMR",
    ),
  ]
  let regions = @src.merge_regions(tiles, params)
  // Two separate regions (gap between 1000 and 3000).
  assert_eq(regions.length(), 2)
}

///|
test "msr_full_pipeline_sample_data" {
  let params = @src.MethylSeekRParams::default()
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  // Four regimes: one UMR, one LMR, one PMD, one FMR.
  assert_eq(regions.length(), 4)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  assert_eq(@src.MethylRegion::region_type(regions[1]), "LMR")
  assert_eq(@src.MethylRegion::region_type(regions[2]), "PMD")
  assert_eq(@src.MethylRegion::region_type(regions[3]), "FMR")
  // UMR: tiles 0-4, span 0-5000.
  assert_eq(@src.MethylRegion::start(regions[0]), 0)
  assert_eq(@src.MethylRegion::end(regions[0]), 5000)
  assert_eq(@src.MethylRegion::n_tiles(regions[0]), 5)
  // LMR: tiles 5-9, span 5000-10000.
  assert_eq(@src.MethylRegion::start(regions[1]), 5000)
  assert_eq(@src.MethylRegion::end(regions[1]), 10000)
  assert_eq(@src.MethylRegion::n_tiles(regions[1]), 5)
  // PMD: tiles 10-24, span 10000-25000.
  assert_eq(@src.MethylRegion::start(regions[2]), 10000)
  assert_eq(@src.MethylRegion::end(regions[2]), 25000)
  assert_eq(@src.MethylRegion::n_tiles(regions[2]), 15)
  // FMR: tiles 25-29, span 25000-30000.
  assert_eq(@src.MethylRegion::start(regions[3]), 25000)
  assert_eq(@src.MethylRegion::end(regions[3]), 30000)
  assert_eq(@src.MethylRegion::n_tiles(regions[3]), 5)
}

///|
test "msr_full_pipeline_region_methylation" {
  let params = @src.MethylSeekRParams::default()
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  // UMR mean ~ 0.05, PMD mean ~ 0.6, FMR mean ~ 0.9.
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[0]) - 0.05).abs() < 0.001,
  )
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[2]) - 0.6).abs() < 0.001,
  )
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[3]) - 0.9).abs() < 0.001,
  )
}

///|
test "msr_filter_regions" {
  let params = @src.MethylSeekRParams::default()
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  let umrs = @src.filter_regions(regions, "UMR")
  assert_eq(umrs.length(), 1)
  assert_eq(@src.MethylRegion::region_type(umrs[0]), "UMR")
  let pmds = @src.filter_regions(regions, "PMD")
  assert_eq(pmds.length(), 1)
  assert_eq(@src.MethylRegion::region_type(pmds[0]), "PMD")
  let fmrs = @src.filter_regions(regions, "FMR")
  assert_eq(fmrs.length(), 1)
  // Filtering a non-existent class returns an empty array.
  let empty = @src.filter_regions(regions, "UNKNOWN")
  assert_eq(empty.length(), 0)
}

///|
test "msr_region_summary" {
  let params = @src.MethylSeekRParams::default()
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  let s = @src.region_summary(regions)
  assert_true(s.contains("Methylation Regimes Summary"))
  assert_true(s.contains("Total regions: 4"))
  assert_true(s.contains("UMR: 1 regions"))
  assert_true(s.contains("LMR: 1 regions"))
  assert_true(s.contains("PMD: 1 regions"))
  assert_true(s.contains("FMR: 1 regions"))
  // UMR spans 5000 bp.
  assert_true(s.contains("UMR: 1 regions, 5000 bp"))
  // PMD spans 15000 bp.
  assert_true(s.contains("PMD: 1 regions, 15000 bp"))
}

///|
test "msr_region_summary_empty" {
  let s = @src.region_summary([])
  assert_true(s.contains("Total regions: 0"))
  assert_true(s.contains("UMR: 0 regions, 0 bp"))
}

///|
test "msr_sample_data_validation" {
  let sites = @src.methyl_seekr_sample_data()
  // 30 tiles * 10 sites = 300 sites.
  assert_eq(sites.length(), 300)
  // All sites are on chr1.
  let mut all_chr1 = true
  for s in sites {
    if @src.CytosineSite::chr(s) != "chr1" {
      all_chr1 = false
    }
  }
  assert_true(all_chr1)
  // First site is at position 50 in the UMR region (1 methylated, 19 unmethylated).
  assert_eq(@src.CytosineSite::position(sites[0]), 50)
  assert_eq(@src.CytosineSite::methylated(sites[0]), 1)
  assert_eq(@src.CytosineSite::unmethylated(sites[0]), 19)
  // Tiling produces 30 tiles.
  let params = @src.MethylSeekRParams::default()
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 30)
  // First tile [0,1000) and last tile [29000,30000).
  assert_eq(@src.MethylTile::start(tiles[0]), 0)
  assert_eq(@src.MethylTile::end(tiles[0]), 1000)
  assert_eq(@src.MethylTile::start(tiles[29]), 29000)
  assert_eq(@src.MethylTile::end(tiles[29]), 30000)
}

///|
test "msr_edge_empty_pipeline" {
  let params = @src.MethylSeekRParams::default()
  let regions = @src.call_methylation_regimes([], params)
  assert_eq(regions.length(), 0)
}

///|
test "msr_edge_single_site" {
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr1", position=100, methylated=0, unmethylated=20),
  ]
  let regions = @src.call_methylation_regimes(sites, params)
  // Single tile of 1000 bp >= min_region_size (500) -> one UMR region.
  assert_eq(regions.length(), 1)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  assert_eq(@src.MethylRegion::n_tiles(regions[0]), 1)
}

///|
test "msr_edge_all_methylated" {
  // All sites fully methylated -> all tiles FMR.
  let params = @src.MethylSeekRParams::default()
  let sites : Array[@src.CytosineSite] = Array::new()
  for k in 0..<12 {
    sites.push(
      @src.CytosineSite::new(
        chr="chr1",
        position=k * 1000 + 50,
        methylated=10,
        unmethylated=0,
      ),
    )
  }
  let regions = @src.call_methylation_regimes(sites, params)
  assert_eq(regions.length(), 1)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "FMR")
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[0]) - 1.0).abs() < 0.000001,
  )
}

///|
test "msr_edge_all_unmethylated" {
  // All sites fully unmethylated with sufficient coverage -> UMR (not PMD),
  // because UMR classification takes precedence over PMD detection.
  let params = @src.MethylSeekRParams::default()
  let sites : Array[@src.CytosineSite] = Array::new()
  for k in 0..<12 {
    sites.push(
      @src.CytosineSite::new(
        chr="chr1",
        position=k * 1000 + 50,
        methylated=0,
        unmethylated=10,
      ),
    )
  }
  let regions = @src.call_methylation_regimes(sites, params)
  assert_eq(regions.length(), 1)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  assert_true(
    (@src.MethylRegion::mean_methylation(regions[0]) - 0.0).abs() < 0.000001,
  )
}

///|
test "msr_edge_zero_coverage_sites" {
  // Sites with zero coverage produce a tile with coverage 0 and level 0.0,
  // which classifies as FMR (insufficient coverage).
  let params = @src.MethylSeekRParams::default()
  let sites = [
    @src.CytosineSite::new(chr="chr1", position=100, methylated=0, unmethylated=0),
  ]
  let tiles = @src.tile_methylation(sites, params)
  assert_eq(tiles.length(), 1)
  assert_eq(@src.MethylTile::coverage(tiles[0]), 0)
  assert_eq(@src.MethylTile::methylation_level(tiles[0]), 0.0)
  let out = @src.classify_tiles(tiles, params)
  assert_eq(@src.MethylTile::region_type(out[0]), "FMR")
}

///|
test "msr_pipeline_custom_params" {
  // With a smaller tile size and pmd_window the PMD domain still forms.
  let params = @src.MethylSeekRParams::new(tile_size=500, pmd_window=5000)
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  // Each original 1 kb tile splits into two 500 bp tiles; regimes are preserved.
  assert_eq(regions.length(), 4)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "UMR")
  assert_eq(@src.MethylRegion::region_type(regions[1]), "LMR")
  assert_eq(@src.MethylRegion::region_type(regions[2]), "PMD")
  assert_eq(@src.MethylRegion::region_type(regions[3]), "FMR")
}

///|
test "msr_pipeline_large_min_region_size" {
  // A large min_region_size keeps only the broad PMD domain and discards the
  // smaller UMR / LMR / FMR regions (each 5 kb).
  let params = @src.MethylSeekRParams::new(min_region_size=12000)
  let sites = @src.methyl_seekr_sample_data()
  let regions = @src.call_methylation_regimes(sites, params)
  assert_eq(regions.length(), 1)
  assert_eq(@src.MethylRegion::region_type(regions[0]), "PMD")
  // With an even larger threshold the PMD (15 kb) is also discarded.
  let params2 = @src.MethylSeekRParams::new(min_region_size=20000)
  let regions2 = @src.call_methylation_regimes(sites, params2)
  assert_eq(regions2.length(), 0)
}
