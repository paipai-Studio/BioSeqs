///|
/// methylSeekR Demo.
/// Demonstrates identifying methylation regimes (UMR / LMR / PMD / FMR)
/// from whole-genome bisulfite sequencing data.
fn main {
  println("=== methylSeekR Demo ===")

  // 1. Load sample cytosine sites.
  println("\n1. Load sample cytosine sites:")
  let sites = @bio.methyl_seekr_sample_data()
  println("  Loaded " + sites.length().to_string() + " cytosine sites")
  println("  First site: " + @bio.CytosineSite::to_string(sites[0]))
  println("  Last site:  " + @bio.CytosineSite::to_string(sites[sites.length() - 1]))

  // 2. Configure pipeline parameters.
  println("\n2. Pipeline parameters (defaults):")
  let params = @bio.MethylSeekRParams::default()
  println("  tile_size=" + @bio.MethylSeekRParams::tile_size(params).to_string())
  println("  min_coverage=" + @bio.MethylSeekRParams::min_coverage(params).to_string())
  println("  alpha=" + @bio.MethylSeekRParams::alpha(params).to_string())
  println("  beta=" + @bio.MethylSeekRParams::beta(params).to_string())
  println("  gamma=" + @bio.MethylSeekRParams::gamma(params).to_string())
  println(
    "  min_region_size=" +
    @bio.MethylSeekRParams::min_region_size(params).to_string(),
  )
  println("  pmd_window=" + @bio.MethylSeekRParams::pmd_window(params).to_string())

  // 3. Tile the genome and compute per-tile methylation.
  println("\n3. Tile the genome:")
  let tiles = @bio.tile_methylation(sites, params)
  println("  Produced " + tiles.length().to_string() + " tiles")
  let mut ti = 0
  while ti < tiles.length() && ti < 5 {
    println("    " + @bio.MethylTile::to_string(tiles[ti]))
    ti = ti + 1
  }
  if tiles.length() > 5 {
    println("    ... and " + (tiles.length() - 5).to_string() + " more")
  }

  // 4. Classify tiles into UMR / LMR / FMR.
  println("\n4. Classify tiles:")
  let classified = @bio.classify_tiles(tiles, params)
  let mut n_umr = 0
  let mut n_lmr = 0
  let mut n_fmr = 0
  for t in classified {
    match @bio.MethylTile::region_type(t) {
      "UMR" => n_umr = n_umr + 1
      "LMR" => n_lmr = n_lmr + 1
      _ => n_fmr = n_fmr + 1
    }
  }
  println("  UMR tiles: " + n_umr.to_string())
  println("  LMR tiles: " + n_lmr.to_string())
  println("  FMR tiles: " + n_fmr.to_string())

  // 5. Detect Partially Methylated Domains.
  println("\n5. Detect PMDs:")
  let with_pmds = @bio.detect_pmds(classified, params)
  let mut n_pmd = 0
  for t in with_pmds {
    if @bio.MethylTile::region_type(t) == "PMD" {
      n_pmd = n_pmd + 1
    }
  }
  println("  PMD tiles: " + n_pmd.to_string())

  // 6. Run the full pipeline in one call.
  println("\n6. Full pipeline (call_methylation_regimes):")
  let regions = @bio.call_methylation_regimes(sites, params)
  println("  Detected " + regions.length().to_string() + " regions")
  for r in regions {
    println("    " + @bio.MethylRegion::to_string(r))
  }

  // 7. Display summary statistics.
  println("\n7. Region summary:")
  println(@bio.region_summary(regions))

  // 8. Filter UMRs.
  println("8. Filter UMR regions:")
  let umrs = @bio.filter_regions(regions, "UMR")
  println("  Found " + umrs.length().to_string() + " UMR region(s)")
  for r in umrs {
    println(
      "    " +
      @bio.MethylRegion::chr(r) +
      ":" +
      @bio.MethylRegion::start(r).to_string() +
      "-" +
      @bio.MethylRegion::end(r).to_string() +
      " (n_tiles=" +
      @bio.MethylRegion::n_tiles(r).to_string() +
      ", mean_meth=" +
      @bio.MethylRegion::mean_methylation(r).to_string() +
      ")",
    )
  }

  // 9. Filter PMDs.
  println("\n9. Filter PMD regions:")
  let pmds = @bio.filter_regions(regions, "PMD")
  println("  Found " + pmds.length().to_string() + " PMD region(s)")
  for r in pmds {
    println(
      "    " +
      @bio.MethylRegion::chr(r) +
      ":" +
      @bio.MethylRegion::start(r).to_string() +
      "-" +
      @bio.MethylRegion::end(r).to_string() +
      " (size=" +
      @bio.MethylRegion::size(r).to_string() +
      " bp, mean_meth=" +
      @bio.MethylRegion::mean_methylation(r).to_string() +
      ")",
    )
  }

  // 10. Demonstrate a single cytosine methylation level.
  println("\n10. Single-site methylation level:")
  let site = @bio.CytosineSite::new(
    chr="chr1",
    position=1000,
    methylated=8,
    unmethylated=2,
  )
  println("  " + @bio.CytosineSite::to_string(site))
  println("  methylation_level = " + @bio.methylation_level(site).to_string())

  println("\nmethylSeekR demo completed!")
}
