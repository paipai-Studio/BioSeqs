///|
/// methylSeekR module - identify methylation regimes from whole-genome
/// bisulfite sequencing (WGBS) data.
/// Inspired by the Bioconductor methylSeekR package.
///
/// Classifies genomic regions into:
/// - UMR (Unmethylated Regions): regions with very low methylation
/// - LMR (Low Methylated Regions): regions with intermediate methylation
/// - PMD (Partially Methylated Domains): large domains with reduced methylation
/// - FMR (Fully Methylated Regions): highly methylated regions (default)
///
/// Pipeline:
/// 1. Tile the genome into fixed-size windows and compute per-tile methylation.
/// 2. Classify each tile as UMR / LMR / FMR using alpha / beta thresholds.
/// 3. Detect PMDs as broad contiguous domains with methylation below gamma.
/// 4. Merge adjacent same-class tiles into regions.
/// 5. Filter out regions smaller than a minimum size threshold.

///|
/// A single cytosine methylation call from whole-genome bisulfite sequencing.
pub struct CytosineSite {
  chr : String
  position : Int
  methylated : Int
  unmethylated : Int
}

///|
/// Construct a CytosineSite.
pub fn CytosineSite::new(
  chr~ : String,
  position~ : Int,
  methylated~ : Int,
  unmethylated~ : Int,
) -> CytosineSite {
  CytosineSite::{ chr, position, methylated, unmethylated }
}

///|
/// Chromosome name.
pub fn CytosineSite::chr(self : CytosineSite) -> String {
  self.chr
}

///|
/// Genomic position (1-based or 0-based, caller-defined).
pub fn CytosineSite::position(self : CytosineSite) -> Int {
  self.position
}

///|
/// Number of methylated reads covering this cytosine.
pub fn CytosineSite::methylated(self : CytosineSite) -> Int {
  self.methylated
}

///|
/// Number of unmethylated reads covering this cytosine.
pub fn CytosineSite::unmethylated(self : CytosineSite) -> Int {
  self.unmethylated
}

///|
/// Total read coverage (methylated + unmethylated).
pub fn CytosineSite::coverage(self : CytosineSite) -> Int {
  self.methylated + self.unmethylated
}

///|
/// Human-readable representation.
pub fn CytosineSite::to_string(self : CytosineSite) -> String {
  "CytosineSite(chr=" +
  self.chr +
  ", position=" +
  self.position.to_string() +
  ", methylated=" +
  self.methylated.to_string() +
  ", unmethylated=" +
  self.unmethylated.to_string() +
  ", coverage=" +
  self.coverage().to_string() +
  ")"
}

///|
/// Methylation level of a single cytosine site:
/// methylated / (methylated + unmethylated). Returns 0.0 when coverage is zero.
pub fn methylation_level(site : CytosineSite) -> Double {
  let total = site.methylated + site.unmethylated
  if total <= 0 {
    return 0.0
  }
  site.methylated.to_double() / total.to_double()
}

///|
/// A fixed-size genomic tile with aggregated methylation statistics.
pub struct MethylTile {
  chr : String
  start : Int
  end : Int
  methylation_level : Double
  coverage : Int
  mut region_type : String
}

///|
/// Construct a MethylTile. `region_type` defaults to "FMR".
pub fn MethylTile::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  methylation_level~ : Double,
  coverage~ : Int,
  region_type~ : String = "FMR",
) -> MethylTile {
  MethylTile::{
    chr,
    start,
    end,
    methylation_level,
    coverage,
    region_type,
  }
}

///|
/// Chromosome name.
pub fn MethylTile::chr(self : MethylTile) -> String {
  self.chr
}

///|
/// Tile start coordinate (inclusive).
pub fn MethylTile::start(self : MethylTile) -> Int {
  self.start
}

///|
/// Tile end coordinate (exclusive).
pub fn MethylTile::end(self : MethylTile) -> Int {
  self.end
}

///|
/// Methylation level in [0, 1].
pub fn MethylTile::methylation_level(self : MethylTile) -> Double {
  self.methylation_level
}

///|
/// Total read coverage within the tile.
pub fn MethylTile::coverage(self : MethylTile) -> Int {
  self.coverage
}

///|
/// Current region classification ("UMR", "LMR", "PMD", or "FMR").
pub fn MethylTile::region_type(self : MethylTile) -> String {
  self.region_type
}

///|
/// Human-readable representation.
pub fn MethylTile::to_string(self : MethylTile) -> String {
  "MethylTile(chr=" +
  self.chr +
  ", start=" +
  self.start.to_string() +
  ", end=" +
  self.end.to_string() +
  ", level=" +
  self.methylation_level.to_string() +
  ", coverage=" +
  self.coverage.to_string() +
  ", type=" +
  self.region_type +
  ")"
}

///|
/// A merged genomic region composed of one or more adjacent tiles of the
/// same methylation class.
pub struct MethylRegion {
  chr : String
  start : Int
  end : Int
  region_type : String
  n_tiles : Int
  mean_methylation : Double
  min_methylation : Double
  max_methylation : Double
}

///|
/// Construct a MethylRegion.
pub fn MethylRegion::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  region_type~ : String,
  n_tiles~ : Int,
  mean_methylation~ : Double,
  min_methylation~ : Double,
  max_methylation~ : Double,
) -> MethylRegion {
  MethylRegion::{
    chr,
    start,
    end,
    region_type,
    n_tiles,
    mean_methylation,
    min_methylation,
    max_methylation,
  }
}

///|
/// Chromosome name.
pub fn MethylRegion::chr(self : MethylRegion) -> String {
  self.chr
}

///|
/// Region start coordinate (inclusive).
pub fn MethylRegion::start(self : MethylRegion) -> Int {
  self.start
}

///|
/// Region end coordinate (exclusive).
pub fn MethylRegion::end(self : MethylRegion) -> Int {
  self.end
}

///|
/// Region class ("UMR", "LMR", "PMD", or "FMR").
pub fn MethylRegion::region_type(self : MethylRegion) -> String {
  self.region_type
}

///|
/// Number of tiles merged into this region.
pub fn MethylRegion::n_tiles(self : MethylRegion) -> Int {
  self.n_tiles
}

///|
/// Mean methylation level across the region's tiles.
pub fn MethylRegion::mean_methylation(self : MethylRegion) -> Double {
  self.mean_methylation
}

///|
/// Minimum methylation level among the region's tiles.
pub fn MethylRegion::min_methylation(self : MethylRegion) -> Double {
  self.min_methylation
}

///|
/// Maximum methylation level among the region's tiles.
pub fn MethylRegion::max_methylation(self : MethylRegion) -> Double {
  self.max_methylation
}

///|
/// Size of the region in base pairs (end - start).
pub fn MethylRegion::size(self : MethylRegion) -> Int {
  self.end - self.start
}

///|
/// Human-readable representation.
pub fn MethylRegion::to_string(self : MethylRegion) -> String {
  "MethylRegion(chr=" +
  self.chr +
  ", start=" +
  self.start.to_string() +
  ", end=" +
  self.end.to_string() +
  ", type=" +
  self.region_type +
  ", n_tiles=" +
  self.n_tiles.to_string() +
  ", mean_meth=" +
  self.mean_methylation.to_string() +
  ")"
}

///|
/// Configuration parameters for the methylSeekR pipeline.
pub struct MethylSeekRParams {
  tile_size : Int
  min_coverage : Int
  alpha : Double
  beta : Double
  gamma : Double
  min_region_size : Int
  pmd_window : Int
}

///|
/// Default parameters:
/// tile_size=1000, min_coverage=5, alpha=0.1, beta=0.5, gamma=0.7,
/// min_region_size=500, pmd_window=10000.
pub fn MethylSeekRParams::default() -> MethylSeekRParams {
  MethylSeekRParams::{
    tile_size: 1000,
    min_coverage: 5,
    alpha: 0.1,
    beta: 0.5,
    gamma: 0.7,
    min_region_size: 500,
    pmd_window: 10000,
  }
}

///|
/// Construct parameters with optional overrides. Any omitted argument takes
/// its default value (see `MethylSeekRParams::default`).
pub fn MethylSeekRParams::new(
  tile_size? : Int = 1000,
  min_coverage? : Int = 5,
  alpha? : Double = 0.1,
  beta? : Double = 0.5,
  gamma? : Double = 0.7,
  min_region_size? : Int = 500,
  pmd_window? : Int = 10000,
) -> MethylSeekRParams {
  MethylSeekRParams::{
    tile_size,
    min_coverage,
    alpha,
    beta,
    gamma,
    min_region_size,
    pmd_window,
  }
}

///|
/// Tile size in base pairs.
pub fn MethylSeekRParams::tile_size(self : MethylSeekRParams) -> Int {
  self.tile_size
}

///|
/// Minimum coverage required to confidently classify a tile.
pub fn MethylSeekRParams::min_coverage(self : MethylSeekRParams) -> Int {
  self.min_coverage
}

///|
/// UMR methylation threshold (tiles below this are UMR candidates).
pub fn MethylSeekRParams::alpha(self : MethylSeekRParams) -> Double {
  self.alpha
}

///|
/// LMR methylation threshold (tiles between alpha and beta are LMR).
pub fn MethylSeekRParams::beta(self : MethylSeekRParams) -> Double {
  self.beta
}

///|
/// PMD methylation threshold (broad domains below this are PMD candidates).
pub fn MethylSeekRParams::gamma(self : MethylSeekRParams) -> Double {
  self.gamma
}

///|
/// Minimum region size in base pairs; smaller merged regions are discarded.
pub fn MethylSeekRParams::min_region_size(self : MethylSeekRParams) -> Int {
  self.min_region_size
}

///|
/// Window size in base pairs used for PMD detection.
pub fn MethylSeekRParams::pmd_window(self : MethylSeekRParams) -> Int {
  self.pmd_window
}

///|
/// Tile the genome into fixed-size windows and compute per-tile methylation.
///
/// Each cytosine site is assigned to the tile containing its position, where
/// tile `k` spans `[k * tile_size, (k + 1) * tile_size)`. Methylated and
/// unmethylated counts are summed across all sites in a tile. Only tiles that
/// contain at least one site are emitted, sorted by chromosome then start.
/// Newly created tiles are assigned the default region type "FMR".
pub fn tile_methylation(
  sites : Array[CytosineSite],
  params : MethylSeekRParams,
) -> Array[MethylTile] {
  let tiles : Array[MethylTile] = Array::new()
  let n = sites.length()
  if n == 0 || params.tile_size <= 0 {
    return tiles
  }
  // Aggregate counts per (chr, tile_index).
  let agg : Map[String, (String, Int, Int, Int)] = Map([], capacity=256)
  let keys : Array[String] = Array::new()
  for site in sites {
    if site.position < 0 {
      continue
    }
    let tile_idx = site.position / params.tile_size
    let key = site.chr + ":" + tile_idx.to_string()
    match agg.get(key) {
      Some(entry) =>
        agg[key] = (
          entry.0,
          entry.1,
          entry.2 + site.methylated,
          entry.3 + site.unmethylated,
        )
      None => {
        agg[key] = (site.chr, tile_idx, site.methylated, site.unmethylated)
        keys.push(key)
      }
    }
  }
  // Build an entries list and sort by (chr, tile_index) using selection sort.
  let entries : Array[(String, Int, Int, Int)] = Array::new()
  for key in keys {
    entries.push(agg[key])
  }
  let m = entries.length()
  let mut i = 0
  while i < m - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < m {
      let lt = entries[j].0 < entries[min_idx].0 ||
        (entries[j].0 == entries[min_idx].0 && entries[j].1 < entries[min_idx].1)
      if lt {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tmp = entries[i]
      entries[i] = entries[min_idx]
      entries[min_idx] = tmp
    }
    i = i + 1
  }
  // Create tiles in sorted order.
  for entry in entries {
    let (chr, idx, meth, unmeth) = entry
    let total = meth + unmeth
    let level = if total > 0 {
      meth.to_double() / total.to_double()
    } else {
      0.0
    }
    let start = idx * params.tile_size
    tiles.push(
      MethylTile::new(
        chr=chr,
        start=start,
        end=start + params.tile_size,
        methylation_level=level,
        coverage=total,
        region_type="FMR",
      ),
    )
  }
  tiles
}

///|
/// Classify each tile as UMR, LMR, or FMR based on methylation level and
/// coverage. Tiles are modified in place and the input array is returned.
///
/// - Coverage below `min_coverage` -> "FMR" (insufficient evidence).
/// - Methylation < `alpha`           -> "UMR".
/// - Methylation in [`alpha`, `beta`) -> "LMR".
/// - Otherwise                       -> "FMR".
///
/// PMD detection is performed separately by `detect_pmds`.
pub fn classify_tiles(
  tiles : Array[MethylTile],
  params : MethylSeekRParams,
) -> Array[MethylTile] {
  for i in 0..<tiles.length() {
    if tiles[i].coverage < params.min_coverage {
      tiles[i].region_type = "FMR"
    } else if tiles[i].methylation_level < params.alpha {
      tiles[i].region_type = "UMR"
    } else if tiles[i].methylation_level < params.beta {
      tiles[i].region_type = "LMR"
    } else {
      tiles[i].region_type = "FMR"
    }
  }
  tiles
}

///|
/// Detect Partially Methylated Domains (PMDs) using a sliding window over
/// contiguous tiles. Tiles are modified in place and the input array returned.
///
/// A tile becomes a PMD candidate when it is currently classified as "FMR"
/// and its methylation level is below `gamma`. Maximal contiguous runs of
/// such candidates (on the same chromosome, with adjacent coordinates) that
/// span at least `pmd_window / tile_size` tiles are relabelled as "PMD".
/// Runs shorter than the minimum window are left unchanged.
pub fn detect_pmds(
  tiles : Array[MethylTile],
  params : MethylSeekRParams,
) -> Array[MethylTile] {
  let n = tiles.length()
  if n == 0 {
    return tiles
  }
  let min_pmd_tiles = if params.tile_size > 0 {
    params.pmd_window / params.tile_size
  } else {
    0
  }
  if min_pmd_tiles <= 0 {
    return tiles
  }
  let mut i = 0
  while i < n {
    // Skip tiles that are not PMD candidates.
    if tiles[i].region_type != "FMR" || tiles[i].methylation_level >= params.gamma {
      i = i + 1
      continue
    }
    let chr = tiles[i].chr
    // Collect a maximal contiguous run of candidate tiles.
    let run : Array[Int] = Array::new()
    run.push(i)
    let mut j = i + 1
    while j < n {
      if tiles[j].chr == chr &&
        tiles[j].start == tiles[j - 1].end &&
        tiles[j].region_type == "FMR" &&
        tiles[j].methylation_level < params.gamma {
        run.push(j)
        j = j + 1
      } else {
        break
      }
    }
    if run.length() >= min_pmd_tiles {
      for idx in run {
        tiles[idx].region_type = "PMD"
      }
    }
    // Advance past the processed run regardless of outcome.
    i = j
  }
  tiles
}

///|
/// Merge adjacent tiles of the same class into MethylRegions. Tiles must be
/// sorted by (chr, start); two tiles merge only when they share the same
/// chromosome, the same region type, and touch (next.start == current.end).
/// Regions smaller than `min_region_size` base pairs are discarded.
pub fn merge_regions(
  tiles : Array[MethylTile],
  params : MethylSeekRParams,
) -> Array[MethylRegion] {
  let result : Array[MethylRegion] = Array::new()
  let n = tiles.length()
  if n == 0 {
    return result
  }
  let mut i = 0
  while i < n {
    let chr = tiles[i].chr
    let rtype = tiles[i].region_type
    let start = tiles[i].start
    let mut end = tiles[i].end
    let mut n_tiles = 1
    let mut sum_level = tiles[i].methylation_level
    let mut min_level = tiles[i].methylation_level
    let mut max_level = tiles[i].methylation_level
    let mut j = i + 1
    while j < n {
      if tiles[j].chr == chr &&
        tiles[j].region_type == rtype &&
        tiles[j].start == end {
        end = tiles[j].end
        n_tiles = n_tiles + 1
        sum_level = sum_level + tiles[j].methylation_level
        if tiles[j].methylation_level < min_level {
          min_level = tiles[j].methylation_level
        }
        if tiles[j].methylation_level > max_level {
          max_level = tiles[j].methylation_level
        }
        j = j + 1
      } else {
        break
      }
    }
    if end - start >= params.min_region_size {
      result.push(
        MethylRegion::new(
          chr=chr,
          start=start,
          end=end,
          region_type=rtype,
          n_tiles=n_tiles,
          mean_methylation=sum_level / n_tiles.to_double(),
          min_methylation=min_level,
          max_methylation=max_level,
        ),
      )
    }
    i = j
  }
  result
}

///|
/// Run the full methylSeekR pipeline:
/// tile -> classify -> detect PMDs -> merge -> filter by size.
pub fn call_methylation_regimes(
  sites : Array[CytosineSite],
  params : MethylSeekRParams,
) -> Array[MethylRegion] {
  let tiles = tile_methylation(sites, params)
  let classified = classify_tiles(tiles, params)
  let with_pmds = detect_pmds(classified, params)
  merge_regions(with_pmds, params)
}

///|
/// Return only the regions whose type matches `region_type`
/// (e.g. "UMR", "LMR", "PMD", or "FMR").
pub fn filter_regions(
  regions : Array[MethylRegion],
  region_type : String,
) -> Array[MethylRegion] {
  let result : Array[MethylRegion] = Array::new()
  for r in regions {
    if r.region_type == region_type {
      result.push(r)
    }
  }
  result
}

///|
/// Build a human-readable summary of a set of methylation regions, reporting
/// the count and total base-pair span of each regime class.
pub fn region_summary(regions : Array[MethylRegion]) -> String {
  let mut n_umr = 0
  let mut n_lmr = 0
  let mut n_pmd = 0
  let mut n_fmr = 0
  let mut bp_umr = 0
  let mut bp_lmr = 0
  let mut bp_pmd = 0
  let mut bp_fmr = 0
  for r in regions {
    let span = r.end - r.start
    if r.region_type == "UMR" {
      n_umr = n_umr + 1
      bp_umr = bp_umr + span
    }
    if r.region_type == "LMR" {
      n_lmr = n_lmr + 1
      bp_lmr = bp_lmr + span
    }
    if r.region_type == "PMD" {
      n_pmd = n_pmd + 1
      bp_pmd = bp_pmd + span
    }
    if r.region_type == "FMR" {
      n_fmr = n_fmr + 1
      bp_fmr = bp_fmr + span
    }
  }
  let mut s = "Methylation Regimes Summary\n"
  s = s + "===========================\n"
  s = s + "Total regions: " + regions.length().to_string() + "\n"
  s = s + "UMR: " + n_umr.to_string() + " regions, " + bp_umr.to_string() + " bp\n"
  s = s + "LMR: " + n_lmr.to_string() + " regions, " + bp_lmr.to_string() + " bp\n"
  s = s + "PMD: " + n_pmd.to_string() + " regions, " + bp_pmd.to_string() + " bp\n"
  s = s + "FMR: " + n_fmr.to_string() + " regions, " + bp_fmr.to_string() + " bp\n"
  s
}

///|
/// Build a set of example cytosine sites for testing and demonstration.
///
/// The data spans 30 kb of chr1 with 10 sites per 1 kb tile (300 sites total)
/// arranged into four regimes:
/// - 0-5 kb:    5% methylation  -> UMR  (tiles 0-4)
/// - 5-10 kb:   30% methylation -> LMR  (tiles 5-9)
/// - 10-25 kb:  60% methylation -> PMD  (tiles 10-24, broad domain < gamma)
/// - 25-30 kb:  90% methylation -> FMR  (tiles 25-29)
pub fn methyl_seekr_sample_data() -> Array[CytosineSite] {
  let sites : Array[CytosineSite] = Array::new()
  // Each entry: (first_tile, last_tile, methylated, unmethylated) per site.
  let regions : Array[(Int, Int, Int, Int)] = [
    (0, 4, 1, 19), // UMR: 5% methylation
    (5, 9, 3, 7), // LMR: 30% methylation
    (10, 24, 6, 4), // PMD: 60% methylation
    (25, 29, 9, 1), // FMR: 90% methylation
  ]
  for region in regions {
    let (t_start, t_end, meth, unmeth) = region
    for t in t_start..=t_end {
      let base = t * 1000
      for k in 0..<10 {
        let pos = base + 50 + k * 100
        sites.push(
          CytosineSite::new(
            chr="chr1",
            position=pos,
            methylated=meth,
            unmethylated=unmeth,
          ),
        )
      }
    }
  }
  sites
}
