// Copy-number segmentation - port of Bioconductor `DNAcopy`.
//
// Implements Circular Binary Segmentation (CBS) of Olshen et al. (2004)
// with the Venkatraman & Olshen (2007) max-statistic formulation used by
// DNAcopy::segment(): the data are recursively split at the interval with
// maximal standardized difference from the rest, significance is assessed
// by permutation, and adjacent segments with nearly equal levels can be
// merged ("sdundo").

///|
/// Copy-number array (single sample), analogous to DNAcopy's `CNA` object.
pub struct Cna {
  chrom : Array[String]
  maploc : Array[Int]
  logratio : Array[Double]
}

///|
/// Build a copy-number array. `maploc` positions are used as reported
/// segment boundaries; probes are segmented in the given order.
pub fn Cna::new(
  chrom : Array[String],
  maploc : Array[Int],
  logratio : Array[Double],
) -> Cna raise SeqError {
  let n = chrom.length()
  if maploc.length() != n || logratio.length() != n {
    raise SeqError("chrom, maploc and logratio must have the same length")
  }
  { chrom, maploc, logratio }
}

///|
/// Segmentation options mirroring `DNAcopy::segment()` defaults.
pub struct DnaCopyOptions {
  alpha : Double
  nperm : Int
  min_width : Int
  trim : Double
  undo_splits : String
  undo_sd : Double
  seed : Int
}

///|
/// Default options identical to `segment()` defaults; `seed` makes the
/// permutation test reproducible.
pub fn DnaCopyOptions::default() -> DnaCopyOptions {
  {
    alpha: 0.01,
    nperm: 10000,
    min_width: 2,
    trim: 0.025,
    undo_splits: "none",
    undo_sd: 3.0,
    seed: 1,
  }
}

///|
/// Construct options with explicit overrides; unset parameters keep the
/// DNAcopy defaults.
pub fn DnaCopyOptions::of(
  alpha? : Double = 0.01,
  nperm? : Int = 10000,
  min_width? : Int = 2,
  trim? : Double = 0.025,
  undo_splits? : String = "none",
  undo_sd? : Double = 3.0,
  seed? : Int = 1,
) -> DnaCopyOptions {
  { alpha, nperm, min_width, trim, undo_splits, undo_sd, seed }
}

///|
/// One segmented copy-number region.
pub struct CnaSegment {
  chrom : String
  start : Int
  end : Int
  num_mark : Int
  seg_mean : Double
  start_row : Int
  end_row : Int
}

///|
/// Segmentation result, analogous to the `output` slot of a DNAcopy object.
pub struct DnaCopyResult {
  segments : Array[CnaSegment]
}

///|
/// Deterministic Lehmer (MINSTD) PRNG with Schrage decomposition for the
/// permutation tests - the same generator scheme used across this library.
struct DcRng {
  mut state : Int
}

///|
fn DcRng::new(seed : Int) -> DcRng {
  // map any seed into [1, 2147483646], the valid multiplicative group
  let s = (seed - 1) % 2147483646 + 1
  let rng : DcRng = { state: if s > 0 { s } else { s + 2147483646 } }
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
fn DcRng::next(self : DcRng) -> Int {
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
fn DcRng::next_uniform(self : DcRng) -> Double {
  self.next().to_double() / 2147483647.0
}

///|
/// In-place Fisher-Yates shuffle of a value slice [lo, hi).
fn dc_shuffle(values : Array[Double], lo : Int, hi : Int, rng : DcRng) -> Unit {
  let mut i = hi - 1
  while i > lo {
    let j = lo + (rng.next_uniform() * (i - lo + 1).to_double()).to_int()
    let tmp = values[i]
    values[i] = values[j]
    values[j] = tmp
    i = i - 1
  }
}

///|
/// Trimmed variance of first differences (DNAcopy `trimmed.variance`).
fn dc_trimmed_variance(x : Array[Double], lo : Int, hi : Int, trim : Double) -> Double {
  let n = hi - lo
  if n < 2 {
    return 0.0
  }
  let diffs : Array[Double] = []
  let mut i = lo + 1
  while i < hi {
    let d = x[i] - x[i - 1]
    diffs.push(if d < 0.0 { -d } else { d })
    i = i + 1
  }
  diffs.sort()
  let n_keep = ((1.0 - 2.0 * trim) * diffs.length().to_double() + 0.5).floor().to_int()
  if n_keep <= 0 {
    return 0.0
  }
  let a = stat_normal_quantile(1.0 - trim)
  let grid_n = 10000
  let mut sum1 = 0.0
  let mut k = 0
  while k < grid_n {
    let x1 = -a + (k.to_double() + 0.5) * (2.0 * a) / grid_n.to_double()
    let dens = 0.3989422804014327 * @math.exp(-x1 * x1 / 2.0)
    sum1 = sum1 + x1 * x1 * dens
    k = k + 1
  }
  let infl = 1.0 / (sum1 * (2.0 * a) / grid_n.to_double() / (1.0 - 2.0 * trim))
  let mut ss = 0.0
  let mut j = 0
  while j < n_keep {
    ss = ss + diffs[j] * diffs[j]
    j = j + 1
  }
  infl * ss / (2.0 * n_keep.to_double())
}

///|
/// Maximal CBS statistic within [lo, hi) and its interval [i1, i2) in
/// absolute probe indices. Port of DNAcopy `tmaxo`: over all intervals of
/// length in [min_width, n - min_width], maximize
/// |mean_int - mean_out| * sqrt(len * (n - len) / n).
fn dc_max_stat(
  x : Array[Double],
  lo : Int,
  hi : Int,
  min_width : Int,
) -> (Int, Int, Double) {
  let n = hi - lo
  // partial sums of centered data
  let s : Array[Double] = Array::make(n + 1, 0.0)
  let mut mean = 0.0
  let mut i = lo
  while i < hi {
    mean = mean + x[i]
    i = i + 1
  }
  mean = mean / n.to_double()
  let mut k = 1
  while k <= n {
    s[k] = s[k - 1] + (x[lo + k - 1] - mean)
    k = k + 1
  }
  let rn = n.to_double()
  let mut best = 0.0
  let mut best_i = 0
  let mut best_j = 0
  let mut i2 = min_width
  while i2 <= n - min_width {
    // i2 = interval end (exclusive), i1 = interval start (exclusive bound
    // of partial sums); interval = [i1, i2)
    let mut i1 = 0
    while i1 <= i2 - min_width {
      let len = (i2 - i1).to_double()
      let scale = rn / (len * (rn - len))
      let diff = s[i2] - s[i1]
      let stat = diff * diff * scale
      if stat > best {
        best = stat
        best_i = i1
        best_j = i2
      }
      i1 = i1 + 1
    }
    i2 = i2 + 1
  }
  (best_i, best_j, best.sqrt())
}

///|
/// Permutation p-value of the observed CBS statistic: the fraction of
/// permuted segments whose maximal statistic reaches the observed one.
/// Port of the `p.method = "perm"` branch of DNAcopy `fndcpt`.
fn dc_perm_test(
  x : Array[Double],
  lo : Int,
  hi : Int,
  min_width : Int,
  ostat : Double,
  alpha : Double,
  nperm : Int,
  rng : DcRng,
) -> Bool {
  let nrejc = (alpha * nperm.to_double()).to_int()
  let work = Array::make(hi - lo, 0.0)
  let mut i = lo
  while i < hi {
    work[i - lo] = x[i]
    i = i + 1
  }
  let mut np = 0
  let mut nrej = 0
  dc_shuffle(work, 0, work.length(), rng)
  while np < nperm {
    // shuffle then recompute max statistic on the permuted copy
    let (_, _, pstat) = dc_max_stat(work, 0, work.length(), min_width)
    if ostat <= pstat {
      nrej = nrej + 1
      if nrej > nrejc {
        // p-value clearly exceeds alpha: the observed statistic is not
        // extreme, so no changepoint (DNAcopy's permax early break).
        return false
      }
    }
    np = np + 1
    if np < nperm {
      dc_shuffle(work, 0, work.length(), rng)
    }
  }
  // few permutations reached the observed statistic: p <= alpha, significant
  nrej <= nrejc
}

///|
/// Binary-segmentation permutation p-value (DNAcopy `tpermp`): tests
/// whether the first n1 of [lo, lo+n1+n2) differ from the last n2.
fn dc_tpermp(
  x : Array[Double],
  lo : Int,
  n1 : Int,
  n2 : Int,
  nperm : Int,
  rng : DcRng,
) -> Double {
  let n = n1 + n2
  if n1 == 1 || n2 == 1 {
    return 1.0
  }
  let px : Array[Double] = Array::make(n, 0.0)
  let mut sum1 = 0.0
  let mut sum2 = 0.0
  let mut tss = 0.0
  let mut i = 0
  while i < n {
    px[i] = x[lo + i]
    if i < n1 {
      sum1 = sum1 + px[i]
    } else {
      sum2 = sum2 + px[i]
    }
    tss = tss + px[i] * px[i]
    i = i + 1
  }
  let xbar = (sum1 + sum2) / n.to_double()
  let m1 = if n1 <= n2 { n1 } else { n2 }
  let rm1 = m1.to_double()
  let mean1 = if n1 <= n2 { sum1 / n1.to_double() } else { sum2 / n2.to_double() }
  let ostat = 0.99999 * (mean1 - xbar).abs()
  let mut tstat = ostat * ostat * rm1 * n.to_double() / (n.to_double() - rm1)
  tstat = tstat / ((tss - n.to_double() * xbar * xbar - tstat) / (n - 2).to_double())
  if tstat > 25.0 && m1 >= 10 {
    return 0.0
  }
  let mut nrej = 0
  let mut np = 0
  while np < nperm {
    dc_shuffle(px, 0, n, rng)
    // resample mean of the smaller group from the tail of the shuffle
    let m1 = if n1 <= n2 { n1 } else { n2 }
    let mut psum = 0.0
    let mut i2 = n
    while i2 > n - m1 {
      psum = psum + px[i2 - 1]
      i2 = i2 - 1
    }
    let pstat = (psum / m1.to_double() - xbar).abs()
    if ostat <= pstat {
      nrej = nrej + 1
    }
    np = np + 1
  }
  nrej.to_double() / nperm.to_double()
}

///|
/// CBS changepoint search within one chromosome segment [lo, hi).
/// Returns probe indices (absolute) where segments end, excluding `hi`.
/// Port of DNAcopy `changepoints` + `fndcpt` (p.method = "perm").
fn dc_cbs_segment(
  x : Array[Double],
  lo : Int,
  hi : Int,
  options : DnaCopyOptions,
  rng : DcRng,
) -> Array[Int] {
  let min_width = options.min_width
  let seg_end : Array[Int] = [lo, hi]
  let change_loc : Array[Int] = []
  let mut k = seg_end.length()
  while k > 1 {
    let cur_lo = seg_end[k - 2]
    let cur_hi = seg_end[k - 1]
    let cur_n = cur_hi - cur_lo
    let mut ncpt = 0
    let mut icpt : Array[Int] = []
    if cur_n >= 2 * min_width {
      // skip constant segments
      let mut lo_v = x[cur_lo]
      let mut hi_v = x[cur_lo]
      let mut i = cur_lo + 1
      while i < cur_hi {
        if x[i] < lo_v {
          lo_v = x[i]
        }
        if x[i] > hi_v {
          hi_v = x[i]
        }
        i = i + 1
      }
      if hi_v > lo_v {
        let (i1, i2, ostat) = dc_max_stat(x, cur_lo, cur_hi, min_width)
        if ostat > 0.1 {
          let interval_len = i2 - i1
          let rest_len = cur_n - interval_len
          let flank = if interval_len < rest_len { interval_len } else { rest_len }
          let mut significant = false
          if ostat >= 7.0 && flank >= 10 {
            significant = true
          } else {
            significant = dc_perm_test(
              x,
              cur_lo,
              cur_hi,
              min_width,
              ostat,
              options.alpha,
              options.nperm,
              rng,
            )
          }
          if significant {
            // locate changepoints relative to cur_lo
            if i2 == cur_n {
              ncpt = 1
              icpt = [i1]
            } else if i1 == 0 {
              ncpt = 1
              icpt = [i2]
            } else {
              // validate the left and right edges by binary segmentation
              let pval_left = dc_tpermp(x, cur_lo, i1, i2 - i1, options.nperm, rng)
              if pval_left <= options.alpha {
                ncpt = ncpt + 1
                icpt.push(i1)
              }
              let pval_right = dc_tpermp(
                x,
                cur_lo + i1,
                cur_n - i1 - (cur_n - i2),
                cur_n - i2,
                options.nperm,
                rng,
              )
              if pval_right <= options.alpha {
                ncpt = ncpt + 1
                icpt.push(i2)
              }
            }
          }
        }
      }
    }
    if ncpt == 0 {
      change_loc.push(seg_end[k - 1])
      seg_end.pop() |> ignore
    } else if ncpt == 1 {
      seg_end[k - 1] = cur_lo + icpt[0]
      seg_end.push(cur_hi)
    } else {
      seg_end[k - 1] = cur_lo + icpt[0]
      seg_end.push(cur_lo + icpt[1])
      seg_end.push(cur_hi)
    }
    k = seg_end.length()
  }
  let reversed : Array[Int] = []
  let mut ri = change_loc.length() - 1
  while ri >= 0 {
    reversed.push(change_loc[ri])
    ri = ri - 1
  }
  reversed
}

///|
/// Merge adjacent segments whose medians differ by less than
/// `trimmed_sd * undo_sd` (DNAcopy `changepoints.sdundo`).
fn dc_sdundo(
  x : Array[Double],
  bounds : Array[Int],
  lo : Int,
  hi : Int,
  trimmed_sd : Double,
  undo_sd : Double,
) -> Array[Int] {
  let threshold = trimmed_sd * undo_sd
  let cpt_loc : Array[Int] = []
  for b in bounds {
    if b > lo && b < hi {
      cpt_loc.push(b)
    }
  }
  cpt_loc.push(hi)
  let mut merged = true
  while merged {
    merged = false
    let k = cpt_loc.length()
    if k > 1 {
      let seg_start : Array[Int] = Array::make(k, 0)
      let mut prev = lo
      let mut i = 0
      while i < k {
        seg_start[i] = prev
        prev = cpt_loc[i]
        i = i + 1
      }
      let med : Array[Double] = Array::make(k, 0.0)
      let mut i2 = 0
      while i2 < k {
        let slice : Array[Double] = []
        let mut j = seg_start[i2]
        while j < cpt_loc[i2] {
          slice.push(x[j])
          j = j + 1
        }
        slice.sort()
        med[i2] = lowess_median(slice)
        i2 = i2 + 1
      }
      let mut best_diff = 1.0e300
      let mut best_i = -1
      let mut i3 = 1
      while i3 < k {
        let d = (med[i3] - med[i3 - 1]).abs()
        if d < best_diff {
          best_diff = d
          best_i = i3 - 1
        }
        i3 = i3 + 1
      }
      if best_i >= 0 && best_diff < threshold {
        cpt_loc.remove(best_i) |> ignore
        merged = true
      }
    }
  }
  let out : Array[Int] = []
  for b in cpt_loc {
    if b < hi {
      out.push(b)
    }
  }
  out
}

///|
/// Segment a copy-number array with Circular Binary Segmentation.
///
/// Each chromosome is segmented independently; the result rows follow the
/// order of the chromosomes in the input. Set `undo_splits` to "none" or
/// "sdundo" (DNAcopy's "prune" mode is not supported).
pub fn segment(cna : Cna, options? : DnaCopyOptions = DnaCopyOptions::default()) -> DnaCopyResult raise SeqError {
  if options.min_width < 2 || options.min_width > 5 {
    raise SeqError("minimum segment width should be between 2 and 5")
  }
  if options.undo_splits != "none" && options.undo_splits != "sdundo" {
    raise SeqError("undo.splits must be \"none\" or \"sdundo\"")
  }
  if options.alpha <= 0.0 || options.alpha >= 1.0 {
    raise SeqError("alpha must be in (0, 1)")
  }
  let n = cna.chrom.length()
  let rng = DcRng::new(options.seed)
  // trimmed SD over the whole sample (as in segment())
  let trimmed_sd = dc_trimmed_variance(cna.logratio, 0, n, options.trim).sqrt()
  // chromosome order of first appearance
  let chrom_order : Array[String] = []
  for ch in cna.chrom {
    let mut found = false
    for c in chrom_order {
      if c == ch {
        found = true
        break
      }
    }
    if !found {
      chrom_order.push(ch)
    }
  }
  let segments : Array[CnaSegment] = []
  for ch in chrom_order {
    let idx : Array[Int] = []
    let mut i = 0
    while i < n {
      if cna.chrom[i] == ch {
        idx.push(i)
      }
      i = i + 1
    }
    let vals : Array[Double] = []
    for i in idx {
      vals.push(cna.logratio[i])
    }
    let bounds_raw = dc_cbs_segment(vals, 0, vals.length(), options, rng)
    // keep only interior boundaries
    let bounds : Array[Int] = []
    for b in bounds_raw {
      if b > 0 && b < vals.length() {
        bounds.push(b)
      }
    }
    let final_bounds = if options.undo_splits == "sdundo" && vals.length() > 0 {
      dc_sdundo(vals, bounds, 0, vals.length(), trimmed_sd, options.undo_sd)
    } else {
      bounds
    }
    // build segments
    let starts : Array[Int] = [0]
    for b in final_bounds {
      starts.push(b)
    }
    let m = starts.length()
    let mut si = 0
    while si < m {
      let s = starts[si]
      let e = if si + 1 < m { starts[si + 1] } else { vals.length() }
      let mut mean = 0.0
      let mut j = s
      while j < e {
        mean = mean + vals[j]
        j = j + 1
      }
      mean = mean / (e - s).to_double()
      let rounded = if mean >= 0.0 {
        ((mean * 10000.0) + 0.5).floor() / 10000.0
      } else {
        ((mean * 10000.0) - 0.5).ceil() / 10000.0
      }
      segments.push(
        {
          chrom: ch,
          start: cna.maploc[idx[s]],
          end: cna.maploc[idx[e - 1]],
          num_mark: e - s,
          seg_mean: rounded,
          start_row: idx[s],
          end_row: idx[e - 1],
        },
      )
      si = si + 1
    }
  }
  { segments }
}
