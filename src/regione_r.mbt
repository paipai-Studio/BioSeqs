// Region permutation testing - port of Bioconductor `regioneR`.
//
// Implements the permTest framework of Gel et al. (2016): a region set A is
// repeatedly randomized over a genome (optionally excluding masked regions),
// an association statistic is evaluated on every randomized copy against a
// second region set B, and the observed statistic is compared with the
// permuted distribution to obtain an empirical p-value and a z-score.
//
// Provided randomization strategies mirror regioneR's:
//   - randomizeRegions: places regions of the same widths at random valid
//     positions (per chromosome by default), never overlapping the mask.
//   - circularRandomizeRegions: applies a random circular shift to every
//     region, preserving local structure (mask is not used).
// Provided statistics mirror numOverlaps and meanDistance.

///|
/// Genome layout: chromosome names with their lengths. Coordinates are
/// 1-based, so valid positions on a chromosome run from 1 to its length.
pub struct GenomeLayout {
  chroms : Array[String]
  lengths : Array[Int]
}

///|
/// Build a genome layout from (chromosome, length) pairs.
pub fn GenomeLayout::new(pairs : Array[(String, Int)]) -> GenomeLayout raise SeqError {
  if pairs.length() == 0 {
    raise SeqError("genome layout must contain at least one chromosome")
  }
  let chroms : Array[String] = []
  let lengths : Array[Int] = []
  for pair in pairs {
    if pair.1 <= 0 {
      raise SeqError("chromosome length must be positive")
    }
    for name in chroms {
      if name == pair.0 {
        raise SeqError("duplicate chromosome \{pair.0} in genome layout")
      }
    }
    chroms.push(pair.0)
    lengths.push(pair.1)
  }
  { chroms, lengths }
}

///|
/// Look up a chromosome index by name; -1 when absent.
fn rr_chrom_index(genome : GenomeLayout, chrom : String) -> Int {
  let mut i = 0
  while i < genome.chroms.length() {
    if genome.chroms[i] == chrom {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Permutation test options mirroring regioneR::permTest defaults.
pub struct RegionPermOptions {
  ntimes : Int
  seed : Int
  alternative : String
  per_chromosome : Bool
  allow_overlaps : Bool
}

///|
/// Default options identical to `permTest()` defaults; `seed` makes the
/// randomizations reproducible.
pub fn RegionPermOptions::default() -> RegionPermOptions {
  {
    ntimes: 100,
    seed: 1,
    alternative: "greater",
    per_chromosome: true,
    allow_overlaps: true,
  }
}

///|
/// Construct options with explicit overrides; unset parameters keep the
/// regioneR defaults.
pub fn RegionPermOptions::of(
  ntimes? : Int = 100,
  seed? : Int = 1,
  alternative? : String = "greater",
  per_chromosome? : Bool = true,
  allow_overlaps? : Bool = true,
) -> RegionPermOptions {
  { ntimes, seed, alternative, per_chromosome, allow_overlaps }
}

///|
/// Result of a permutation test: the observed statistic, the permuted
/// statistics, the z-score relative to the permuted distribution and the
/// empirical p-value.
pub struct PermTestResult {
  statistic : Double
  permuted : Array[Double]
  z_score : Double
  p_value : Double
  alternative : String
  ntimes : Int
}

///|
/// Deterministic Lehmer (MINSTD) PRNG with Schrage decomposition for the
/// randomizations - the same generator scheme used across this library.
struct RrRng {
  mut state : Int
}

///|
fn RrRng::new(seed : Int) -> RrRng {
  // map any seed into [1, 2147483646], the valid multiplicative group
  let s = (seed - 1) % 2147483646 + 1
  let rng : RrRng = { state: if s > 0 { s } else { s + 2147483646 } }
  // warm up so that low-entropy seeds do not produce a biased first draw
  let mut i = 0
  while i < 8 {
    ignore(rng.next())
    i = i + 1
  }
  rng
}

///|
/// Advance the generator and return the raw value in [1, 2147483646].
fn RrRng::next(self : RrRng) -> Int {
  let normalized = (self.state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let nxt = if candidate > 0 { candidate } else { candidate + 2147483647 }
  self.state = nxt
  nxt
}

///|
/// Uniform double in [0, 1).
fn RrRng::next_uniform(self : RrRng) -> Double {
  self.next().to_double() / 2147483647.0
}

///|
/// Uniform integer in [0, bound).
fn RrRng::next_bound(self : RrRng, bound : Int) -> Int {
  let r = (self.next_uniform() * bound.to_double()).to_int()
  if r >= bound {
    bound - 1
  } else if r < 0 {
    0
  } else {
    r
  }
}

///|
/// Merge and sort a set of intervals (1-based, inclusive) per chromosome.
/// Returns a map from chromosome name to the sorted, merged interval list.
fn rr_mask_by_chrom(mask : GRanges) -> Map[String, Array[(Int, Int)]] {
  let raw : Map[String, Array[(Int, Int)]] = Map::new()
  let n = granges_length(mask)
  let mut i = 0
  while i < n {
    let chrom = mask.seqnames[i]
    match raw.get(chrom) {
      Some(ivs) => ivs.push((mask.starts[i], mask.ends[i]))
      None => raw.set(chrom, [(mask.starts[i], mask.ends[i])])
    }
    i = i + 1
  }
  let merged : Map[String, Array[(Int, Int)]] = Map::new()
  for chrom, ivs in raw {
    ivs.sort_by(fn(a : (Int, Int), b : (Int, Int)) -> Int {
      let c = a.0.compare(b.0)
      if c == 0 {
        a.1.compare(b.1)
      } else {
        c
      }
    })
    let out : Array[(Int, Int)] = []
    for iv in ivs {
      if out.length() > 0 && iv.0 <= out[out.length() - 1].1 + 1 {
        // contiguous or overlapping: extend the last interval
        let last = out[out.length() - 1]
        if iv.1 > last.1 {
          out[out.length() - 1] = (last.0, iv.1)
        }
      } else {
        out.push(iv)
      }
    }
    merged.set(chrom, out)
  }
  merged
}

///|
/// Validate that every region of `a` lies on a genome chromosome and within
/// its bounds. Returns an error message or "" when valid.
fn rr_validate(a : GRanges, genome : GenomeLayout) -> String {
  let n = granges_length(a)
  let mut i = 0
  while i < n {
    let ci = rr_chrom_index(genome, a.seqnames[i])
    if ci < 0 {
      return "region \{i} chromosome \{a.seqnames[i]} is not in the genome layout"
    }
    let len = genome.lengths[ci]
    if a.starts[i] < 1 || a.ends[i] > len || a.starts[i] > a.ends[i] {
      return "region \{i} (\{a.starts[i]}-\{a.ends[i]}) is out of bounds for \{a.seqnames[i]}"
    }
    i = i + 1
  }
  ""
}

///|
/// Number of valid start positions for a region of width `w` inside a
/// sorted, merged allowed-interval list, plus the interval list itself.
fn rr_total_starts(allowed : Array[(Int, Int)], w : Int) -> Int {
  let mut total = 0
  for iv in allowed {
    let avail = iv.1 - iv.0 + 1
    if avail >= w {
      total = total + avail - w + 1
    }
  }
  total
}

///|
/// Sample a start position for a region of width `w` from the allowed
/// interval list, using `r` uniform in [0, total_starts). Returns the
/// picked interval index together with the absolute start position.
fn rr_pick_start(
  allowed : Array[(Int, Int)],
  w : Int,
  r : Int,
) -> (Int, Int) raise SeqError {
  let mut remaining = r
  let mut idx = 0
  for iv in allowed {
    let avail = iv.1 - iv.0 + 1
    if avail >= w {
      let cnt = avail - w + 1
      if remaining < cnt {
        return (idx, iv.0 + remaining)
      }
      remaining = remaining - cnt
    }
    idx = idx + 1
  }
  raise SeqError("internal error: random position out of range")
}

///|
/// Allowed intervals for one cache entry plus, in global-pool mode, the
/// chromosome index each interval belongs to.
struct RrAllowed {
  ivs : Array[(Int, Int)]
  chrom_ids : Array[Int]
}

///|
/// Collect the allowed (unmasked) intervals of chromosome `ci` inside [1,
/// len], appending to `ivs`; `chrom_ids` receives `ci` per interval.
fn rr_collect_allowed(
  genome : GenomeLayout,
  ci : Int,
  mask_map : Map[String, Array[(Int, Int)]],
  ivs : Array[(Int, Int)],
  chrom_ids : Array[Int],
) -> Unit {
  let len = genome.lengths[ci]
  let masked = match mask_map.get(genome.chroms[ci]) {
    Some(m) => m
    None => []
  }
  let mut pos = 1
  for m in masked {
    if m.0 > pos {
      let hi = if m.0 - 1 < len { m.0 - 1 } else { len }
      if hi >= pos {
        ivs.push((pos, hi))
        chrom_ids.push(ci)
      }
    }
    if m.1 + 1 > pos {
      pos = m.1 + 1
    }
    if pos > len {
      // masked tail reaching the chromosome end
      return
    }
  }
  ivs.push((pos, len))
  chrom_ids.push(ci)
}

///|
/// Core randomization with an externally owned RNG so that repeated calls
/// (one per permutation) advance the same random stream.
fn rr_randomize_with_rng(
  a : GRanges,
  genome : GenomeLayout,
  mask_map : Map[String, Array[(Int, Int)]],
  per_chromosome : Bool,
  rng : RrRng,
) -> GRanges raise SeqError {
  let msg = rr_validate(a, genome)
  if msg != "" {
    raise SeqError(msg)
  }
  let n = granges_length(a)
  let new_seqnames : Array[String] = []
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  // cache of allowed-interval lists: per chromosome + width in
  // per-chromosome mode, per width across the whole genome otherwise
  let iv_cache : Map[String, RrAllowed] = Map::new()
  let mut i = 0
  while i < n {
    let w = a.ends[i] - a.starts[i] + 1
    let ci = rr_chrom_index(genome, a.seqnames[i])
    let cache_key = if per_chromosome {
      "\{ci}:\{w}"
    } else {
      "*:\{w}"
    }
    let allowed = match iv_cache.get(cache_key) {
      Some(entry) => entry
      None => {
        let entry : RrAllowed = { ivs: [], chrom_ids: [] }
        if per_chromosome {
          rr_collect_allowed(genome, ci, mask_map, entry.ivs, entry.chrom_ids)
        } else {
          // global pool: allowed intervals of every chromosome
          let mut c = 0
          while c < genome.chroms.length() {
            rr_collect_allowed(genome, c, mask_map, entry.ivs, entry.chrom_ids)
            c = c + 1
          }
        }
        iv_cache.set(cache_key, entry)
        entry
      }
    }
    let total = rr_total_starts(allowed.ivs, w)
    if total == 0 {
      raise SeqError(
        "not enough unmasked space for a region of width \{w}",
      )
    }
    let (iv_idx, start) = rr_pick_start(allowed.ivs, w, rng.next_bound(total))
    new_seqnames.push(genome.chroms[allowed.chrom_ids[iv_idx]])
    new_starts[i] = start
    new_ends[i] = start + w - 1
    i = i + 1
  }
  {
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: a.strands.copy(),
    widths: a.widths.copy(),
  }
}

///|
/// Randomize region set `a` over `genome`, keeping region widths and (by
/// default) their chromosomes, and never placing a region so that it
/// overlaps `mask`. With `per_chromosome = false` regions may move to any
/// chromosome of the genome (weighted by its unmasked space). Mirrors
/// `regioneR::randomizeRegions`.
pub fn randomize_regions(
  a : GRanges,
  genome : GenomeLayout,
  mask? : GRanges = granges([], [], []),
  options? : RegionPermOptions = RegionPermOptions::default(),
) -> GRanges raise SeqError {
  rr_randomize_with_rng(
    a,
    genome,
    rr_mask_by_chrom(mask),
    options.per_chromosome,
    RrRng::new(options.seed),
  )
}

///|
/// Randomize region set `a` by applying a random circular shift to every
/// region, which preserves the local structure of the region set. Regions
/// keep their chromosome and width; the shift wraps within the valid
/// placement space so all regions remain inside their chromosome. Mirrors
/// `regioneR::circularRandomizeRegions` (the mask is not used).
pub fn circular_randomize_regions(
  a : GRanges,
  genome : GenomeLayout,
  options? : RegionPermOptions = RegionPermOptions::default(),
) -> GRanges raise SeqError {
  let msg = rr_validate(a, genome)
  if msg != "" {
    raise SeqError(msg)
  }
  let rng = RrRng::new(options.seed)
  let n = granges_length(a)
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    let w = a.ends[i] - a.starts[i] + 1
    let span = genome.lengths[rr_chrom_index(genome, a.seqnames[i])] - w + 1
    // circular shift of the start within [1, span]
    let shift = rng.next_bound(span)
    let start = ((a.starts[i] - 1 + shift) % span) + 1
    new_starts[i] = start
    new_ends[i] = start + w - 1
    i = i + 1
  }
  {
    seqnames: a.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: a.strands.copy(),
    widths: a.widths.copy(),
  }
}

///|
/// Number of overlaps between two region sets, mirroring
/// `regioneR::numOverlaps`. With `count_once = false` (the default) every
/// pairwise overlap is counted; with `count_once = true` each region of `a`
/// that overlaps at least one region of `b` is counted once.
pub fn stat_num_overlaps(
  a : GRanges,
  b : GRanges,
  count_once? : Bool = false,
) -> Int {
  let counts = count_overlaps(a, b)
  let mut total = 0
  for c in counts {
    if count_once {
      if c > 0 {
        total = total + 1
      }
    } else {
      total = total + c
    }
  }
  total
}

///|
/// Mean distance from each region of `a` to its nearest region of `b`,
/// mirroring `regioneR::meanDistance`. Regions of `a` whose chromosome has
/// no region in `b` are skipped; returns NaN when no distance exists.
pub fn stat_mean_distance(a : GRanges, b : GRanges) -> Double {
  let ds = granges_distance(a, b)
  let mut sum = 0.0
  let mut cnt = 0
  for d in ds {
    if !d.is_nan() {
      sum = sum + d
      cnt = cnt + 1
    }
  }
  if cnt == 0 {
    @double.not_a_number
  } else {
    sum / cnt.to_double()
  }
}

///|
/// Empirical p-value of the observed statistic against the permuted
/// distribution (regioneR `compute.pvalue`).
fn rr_p_value(
  permuted : Array[Double],
  observed : Double,
  alternative : String,
) -> Double raise SeqError {
  let n = permuted.length()
  let mut ge = 0
  let mut le = 0
  for p in permuted {
    if p >= observed {
      ge = ge + 1
    }
    if p <= observed {
      le = le + 1
    }
  }
  let p_greater = (ge + 1).to_double() / (n + 1).to_double()
  let p_less = (le + 1).to_double() / (n + 1).to_double()
  if alternative == "greater" {
    p_greater
  } else if alternative == "less" {
    p_less
  } else if alternative == "two.sided" {
    let m = if p_greater < p_less { p_greater } else { p_less }
    let two = 2.0 * m
    if two > 1.0 {
      1.0
    } else {
      two
    }
  } else {
    raise SeqError("alternative must be \"greater\", \"less\" or \"two.sided\"")
  }
}

///|
/// z-score of the observed statistic relative to the permuted distribution
/// (regioneR `zscore`), using the sample standard deviation.
fn rr_z_score(permuted : Array[Double], observed : Double) -> Double {
  let n = permuted.length()
  if n == 0 {
    return @double.not_a_number
  }
  let mut sum = 0.0
  for p in permuted {
    sum = sum + p
  }
  let mean = sum / n.to_double()
  if n < 2 {
    return 0.0
  }
  let mut ss = 0.0
  for p in permuted {
    let d = p - mean
    ss = ss + d * d
  }
  let sd = (ss / (n - 1).to_double()).sqrt()
  if sd == 0.0 {
    if observed > mean {
      @double.infinity
    } else if observed < mean {
      -@double.infinity
    } else {
      0.0
    }
  } else {
    (observed - mean) / sd
  }
}

///|
/// Permutation test of the association between two region sets (port of
/// `regioneR::permTest`). Region set `a` is randomized `ntimes` times over
/// `genome` (excluding `mask`, when given) and `statistic` is evaluated on
/// every randomized copy against `b`. The observed statistic is compared
/// with the permuted distribution to obtain the empirical p-value and
/// z-score returned in a `PermTestResult`.
pub fn perm_test(
  a : GRanges,
  b : GRanges,
  genome : GenomeLayout,
  statistic : (GRanges, GRanges) -> Double,
  mask? : GRanges = granges([], [], []),
  options? : RegionPermOptions = RegionPermOptions::default(),
) -> PermTestResult raise SeqError {
  if options.ntimes <= 0 {
    raise SeqError("ntimes must be positive")
  }
  if options.alternative != "greater" &&
    options.alternative != "less" &&
    options.alternative != "two.sided" {
    raise SeqError("alternative must be \"greater\", \"less\" or \"two.sided\"")
  }
  let observed = statistic(a, b)
  let mask_map = rr_mask_by_chrom(mask)
  let rng = RrRng::new(options.seed)
  let permuted : Array[Double] = []
  let mut i = 0
  while i < options.ntimes {
    let ra = rr_randomize_with_rng(a, genome, mask_map, options.per_chromosome, rng)
    permuted.push(statistic(ra, b))
    i = i + 1
  }
  {
    statistic: observed,
    permuted,
    z_score: rr_z_score(permuted, observed),
    p_value: rr_p_value(permuted, observed, options.alternative),
    alternative: options.alternative,
    ntimes: options.ntimes,
  }
}

///|
/// Convenience wrapper of `perm_test` using the number of overlaps as the
/// statistic (port of `regioneR::overlapPermTest`).
pub fn overlap_perm_test(
  a : GRanges,
  b : GRanges,
  genome : GenomeLayout,
  mask? : GRanges = granges([], [], []),
  count_once? : Bool = false,
  options? : RegionPermOptions = RegionPermOptions::default(),
) -> PermTestResult raise SeqError {
  perm_test(
    a,
    b,
    genome,
    fn(x : GRanges, y : GRanges) { stat_num_overlaps(x, y, count_once=count_once).to_double() },
    mask=mask,
    options=options,
  )
}
