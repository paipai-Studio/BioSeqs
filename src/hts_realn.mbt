///|
/// htslib BAQ realignment and probabilistic alignment (realn.c + errmod.c + probaln.c port).
///
/// This module implements:
/// 1. probaln_glocal - probabilistic banded glocal (semi-global) alignment
/// 2. errmod - revised MAQ error model for genotype likelihoods
/// 3. sam_cap_mapq - mapping quality capping
/// 4. sam_prob_realn - BAQ calculation and realignment

///|
/// Alignment parameters for probaln_glocal.
pub struct ProbalnPar {
  d : Double
  e : Double
  bw : Int
}

///|
pub fn ProbalnPar::make(d : Double, e : Double, bw : Int) -> ProbalnPar {
  ProbalnPar::{ d, e, bw }
}

///|
/// Result of probabilistic glocal alignment.
pub struct ProbalnResult {
  score : Int
  state : Array[Int]
  q : Array[Int]
}

///|
let probaln_EI : Double = 0.25

///|
let probaln_EM : Double = 0.33333333333

///|
fn set_u(bw : Int, i : Int, k : Int) -> Int {
  let x = if i - bw > 0 { i - bw } else { 0 }
  (k - x + 1) * 3
}

///|
pub fn base_to_int(base : Char) -> Int {
  match base {
    'A' | 'a' => 0
    'C' | 'c' => 1
    'G' | 'g' => 2
    'T' | 't' => 3
    _ => 4
  }
}

///|
pub fn seq_to_ints(seq : String) -> Array[Int] {
  let result = Array::make(seq.length(), 4)
  for i = 0; i < seq.length(); i = i + 1 {
    result[i] = base_to_int(seq.unsafe_get(i).to_char().unwrap())
  }
  result
}

///|
/// Probabilistic banded glocal (semi-global) alignment.
///
/// Uses a profile HMM with forward-backward algorithm to compute
/// the most probable alignment of query to reference, allowing
/// for gaps at the ends of the reference (semi-global).
///
/// Parameters:
/// - ref: reference sequence (0-4 encoded, 4=ambiguous)
/// - query: query sequence (0-4 encoded)
/// - qual: base quality scores (phred scale, 0-255)
/// - params: alignment parameters (gap open, gap extend, band width)
///
/// Returns: ProbalnResult with score, state, and quality arrays
pub fn probaln_glocal(
  ref : Array[Int],
  query : Array[Int],
  qual : Array[Int],
  params : ProbalnPar,
) -> ProbalnResult {
  let l_ref = ref.length()
  let l_query = query.length()
  if l_ref == 0 || l_query == 0 {
    return ProbalnResult::{ score: 0, state: Array::make(l_query, 0), q: Array::make(l_query, 0) }
  }

  let mut bw = if l_ref > l_query { l_ref } else { l_query }
  if bw > params.bw {
    bw = params.bw
  }
  let diff = if l_ref - l_query > 0 { l_ref - l_query } else { l_query - l_ref }
  if bw < diff {
    bw = diff
  }
  let bw2 = bw * 2 + 1
  let i_dim = if bw2 < l_ref { bw2 * 3 + 6 } else { l_ref * 3 + 6 }

  let qual_prob = Array::make(l_query, 0.0)
  for i = 0; i < l_query; i = i + 1 {
    let q = if i < qual.length() { qual[i] } else { 30 }
    qual_prob[i] = @math.pow(10.0, -(q.to_double()) / 10.0)
  }

  let sM = 1.0 / (2.0 * l_query.to_double() + 2.0)
  let sI = sM
  let m = Array::make(9, 0.0)
  m[0] = (1.0 - params.d - params.d) * (1.0 - sM)
  m[1] = params.d * (1.0 - sM)
  m[2] = params.d * (1.0 - sM)
  m[3] = (1.0 - params.e) * (1.0 - sI)
  m[4] = params.e * (1.0 - sI)
  m[5] = 0.0
  m[6] = 1.0 - params.e
  m[7] = 0.0
  m[8] = params.e
  let bM = (1.0 - params.d) / l_ref.to_double()
  let bI = params.d / l_ref.to_double()

  let f = Array::make((l_query + 1) * i_dim, 0.0)
  let s = Array::make(l_query + 2, 0.0)

  let k0 = set_u(bw, 0, 0)
  f[k0] = 1.0
  s[0] = 1.0

  // f[1]
  {
    let fi_offset = i_dim
    let mut sum = 0.0
    let end1 = if l_ref < bw + 1 { l_ref } else { bw + 1 }
    for k = 1; k <= end1; k = k + 1 {
      let e = if ref[k - 1] > 3 || query[0] > 3 {
        1.0
      } else if ref[k - 1] == query[0] {
        1.0 - qual_prob[0]
      } else {
        qual_prob[0] * probaln_EM
      }
      let u = set_u(bw, 1, k)
      f[fi_offset + u] = e * bM
      f[fi_offset + u + 1] = probaln_EI * bI
      sum = sum + f[fi_offset + u] + f[fi_offset + u + 1]
    }
    s[1] = sum
  }

  // f[2..l_query]
  for i = 2; i <= l_query; i = i + 1 {
    let fi_offset = i * i_dim
    let fi1_offset = (i - 1) * i_dim
    let mut sum = 0.0
    let qli = qual_prob[i - 1]
    let qyi = query[i - 1]
    let mut beg = 1
    let mut end = l_ref
    let x_low = i - bw
    if beg < x_low {
      beg = x_low
    }
    let x_high = i + bw
    if end > x_high {
      end = x_high
    }

    let e_arr = [
      qli * probaln_EM,
      1.0 - qli,
      1.0,
      1.0,
    ]
    let m_scale = 1.0 / s[i - 1]

    let xm0 = m_scale * m[0]
    let xm1 = m_scale * m[3]
    let xm2 = m_scale * m[6]
    let xm3 = probaln_EI * m_scale * m[1]
    let xm4 = probaln_EI * m_scale * m[4]

    let u_start = set_u(bw, i, beg)
    let v11_start = set_u(bw, i - 1, beg - 1)
    let mut xi = u_start
    let mut yi = v11_start
    let mut l_x0 = m[2] * f[fi_offset + xi]
    let mut l_x2 = m[8] * f[fi_offset + xi + 2]

    for k = beg; k <= end; k = k + 1 {
      let cond = if ref[k - 1] > 3 || qyi > 3 {
        2 + (if ref[k - 1] == qyi { 1 } else { 0 })
      } else {
        if ref[k - 1] == qyi { 1 } else { 0 }
      }

      let z0 = xm0 * f[fi1_offset + yi]
      let z1 = xm1 * f[fi1_offset + yi + 1]
      let z2 = xm2 * f[fi1_offset + yi + 2]
      let z3 = xm3 * f[fi1_offset + yi + 3]
      let z4 = xm4 * f[fi1_offset + yi + 4]

      f[fi_offset + xi] = e_arr[cond] * (z0 + z1 + z2)
      f[fi_offset + xi + 1] = z3 + z4
      f[fi_offset + xi + 2] = l_x0 + l_x2
      sum = sum + f[fi_offset + xi] + f[fi_offset + xi + 1] + f[fi_offset + xi + 2]

      l_x0 = m[2] * f[fi_offset + xi]
      l_x2 = m[8] * f[fi_offset + xi + 2]
      xi = xi + 3
      yi = yi + 3
    }
    s[i] = sum
  }

  // Compute likelihood
  let mut pr_val = 0
  {
    let mut sum = 0.0
    let m_scale = 1.0 / s[l_query]
    for k = 1; k <= l_ref; k = k + 1 {
      let u = set_u(bw, l_query, k)
      if u < 3 || u >= i_dim {
        continue
      }
      let idx = l_query * i_dim + u
      sum = sum + m_scale * f[idx] * sM + m_scale * f[idx + 1] * sI
    }
    s[l_query + 1] = sum
  }
  {
    let mut p = 1.0
    let mut pr1 = 0.0
    for i = 0; i <= l_query + 1; i = i + 1 {
      p = p * s[i]
      if p < 0.0000000000000000000000000000001 {
        pr1 = pr1 + -4.343 * @math.ln(p)
        p = 1.0
      }
    }
    pr1 = pr1 + -4.343 * @math.ln(p * l_ref.to_double() * l_query.to_double())
    pr_val = (pr1 + 0.499).to_int()
  }

  // Backward
  let b = Array::make((l_query + 1) * i_dim, 0.0)

  for k = 1; k <= l_ref; k = k + 1 {
    let u = set_u(bw, l_query, k)
    if u < 3 || u >= i_dim {
      continue
    }
    let idx = l_query * i_dim + u
    b[idx] = sM / s[l_query] / s[l_query + 1]
    b[idx + 1] = sI / s[l_query] / s[l_query + 1]
  }

  for i = l_query - 1; i >= 1; i = i - 1 {
    let mut beg = 1
    let mut end = l_ref
    let x_low = i - bw
    if beg < x_low {
      beg = x_low
    }
    let x_high = i + bw
    if end > x_high {
      end = x_high
    }
    let bi_offset = i * i_dim
    let bi1_offset = (i + 1) * i_dim
    let y = if i > 1 { 1.0 } else { 0.0 }
    let qli1 = qual_prob[i]
    let qyi1 = query[i]

    let e_arr = [
      qli1 * probaln_EM,
      1.0 - qli1,
      1.0,
      1.0,
    ]

    let u_start = set_u(bw, i, end)
    let v10_start = set_u(bw, i + 1, end)
    let mut xi = u_start
    let mut yi = v10_start
    let mut xi_5 = b[bi_offset + xi + 5]
    let e1 = probaln_EI * m[1]
    let e4 = probaln_EI * m[4]
    let n = 1.0 / s[i]

    for k = end; k >= beg; k = k - 1 {
      let e = if k >= l_ref {
        0.0
      } else {
        let cond = if ref[k] > 3 || qyi1 > 3 {
          2 + (if ref[k] == qyi1 { 1 } else { 0 })
        } else {
          if ref[k] == qyi1 { 1 } else { 0 }
        }
        e_arr[cond] * b[bi1_offset + yi + 3]
      }

      b[bi_offset + xi + 1] = (e * m[3] + e4 * b[bi1_offset + yi + 1]) * n
      b[bi_offset + xi] = (e * m[0] + e1 * b[bi1_offset + yi + 1] + m[2] * xi_5) * n
      b[bi_offset + xi + 2] = (e * m[6] + m[8] * xi_5) * y * n
      xi_5 = b[bi_offset + xi + 2]

      xi = xi - 3
      yi = yi - 3
    }
  }

  // MAP
  let state = Array::make(l_query, 0)
  let q_arr = Array::make(l_query, 0)

  for i = 1; i <= l_query; i = i + 1 {
    let mut beg = 1
    let mut end = l_ref
    let x_low = i - bw
    if beg < x_low {
      beg = x_low
    }
    let x_high = i + bw
    if end > x_high {
      end = x_high
    }
    let fi_offset = i * i_dim
    let bi_offset = i * i_dim
    let m_scale = 1.0 / s[i]
    let mut sum = 0.0
    let mut max_val = 0.0
    let mut max_k = -1

    let mut u = set_u(bw, i, beg)
    for k = beg; k <= end; k = k + 1 {
      let z1 = m_scale * f[fi_offset + u] * b[bi_offset + u]
      let z2 = m_scale * f[fi_offset + u + 1] * b[bi_offset + u + 1]
      let which = if z2 >= z1 { 1 } else { 0 }
      let zm = if which == 1 { z2 } else { z1 }
      if zm > max_val {
        max_val = zm
        max_k = ((k - 1) << 2) | which
      }
      sum = sum + z1 + z2
      u = u + 3
    }
    if sum > 0.0 {
      max_val = max_val / sum
    }
    if max_val > 1.0 {
      max_val = 1.0
    }
    state[i - 1] = max_k
    let one_minus_max = 1.0 - max_val
    let qval = if one_minus_max <= 0.0 {
      99
    } else {
      (-4.343 * @math.ln(one_minus_max) + 0.499).to_int()
    }
    q_arr[i - 1] = if qval > 100 { 99 } else { qval }
  }

  ProbalnResult::{ score: pr_val, state, q: q_arr }
}

///|
/// Error model for genotype likelihood calculation (errmod.c port).
/// Uses revised MAQ error model with dependency correction.
pub struct ErrMod {
  depcorr : Double
  fk : Array[Double]
  lhet : Array[Double]
}

///|
fn lfact(n : Int) -> Double {
  if n <= 0 {
    return 0.0
  }
  kf_lgamma(n.to_double() + 1.0)
}

///|
/// Compute log binomial coefficient: log(n! / (k! * (n-k)!))
fn logbinom(n : Int, k : Int) -> Double {
  if n <= 0 || k <= 0 || k > n {
    return 0.0
  }
  lfact(n) - lfact(k) - lfact(n - k)
}

///|
/// Create an error model with given dependency correlation.
/// Precomputes fk and lhet tables.
pub fn errmod_init(depcorr : Double) -> ErrMod {
  let eta = 0.03
  let fk = Array::make(256, 0.0)
  fk[0] = 1.0
  for n = 1; n < 256; n = n + 1 {
    fk[n] = @math.pow(1.0 - depcorr, n.to_double()) * (1.0 - eta) + eta
  }

  // lhet table: 256*256 log heterozygote likelihoods
  let lhet = Array::make(256 * 256, 0.0)
  for n = 0; n < 256; n = n + 1 {
    for k = 0; k < 256; k = k + 1 {
      lhet[n * 256 + k] = logbinom(n, k) - 0.6931471805599453 * n.to_double()
    }
  }

  ErrMod::{ depcorr, fk, lhet }
}

///|
/// Compute beta value on demand (instead of precomputing 32MB table).
/// beta[q][n][k] = -10/ln(10) * log(sum of binomial terms)
fn errmod_beta(em : ErrMod, qual : Int, n : Int, k : Int) -> Double {
  if n <= 0 || k > n {
    return 0.0
  }
  if k == n {
    return 1.0e300 // HUGE_VAL equivalent
  }
  let e = @math.pow(10.0, -(qual.to_double()) / 10.0)
  let le = @math.ln(e)
  let le1 = @math.ln(1.0 - e)
  // sum1 starts at k=n: lC[n][n] + n*le
  let mut sum1 = logbinom(n, n) + n.to_double() * le
  // Iterate from k=n-1 down to k
  let mut beta_k = 0.0
  let mut j = n - 1
  while j >= k {
    let lC = logbinom(n, j)
    let val = lC + j.to_double() * le + (n - j).to_double() * le1
    let diff = val - sum1
    let sum = if diff > 700.0 {
      val
    } else if diff < -700.0 {
      sum1
    } else {
      sum1 + @math.ln(1.0 + @math.exp(diff))
    }
    if j == k {
      beta_k = -10.0 / 2.302585092994046 * (sum1 - sum)
    }
    sum1 = sum
    j = j - 1
  }
  beta_k
}

///|
/// Calculate genotype likelihoods given observed bases.
///
/// Parameters:
/// - em: error model
/// - bases: array of observed bases, each encoded as (quality<<5 | strand<<4 | base)
///   where base is 0-4 (A,C,G,T,N) and quality is phred score
/// - m: number of alleles (4 for ACGT)
///
/// Returns: m*m array of phred-scaled genotype likelihoods
pub fn errmod_cal(
  em : ErrMod,
  bases : Array[Int],
  m : Int,
) -> Array[Double] {
  let n_orig = bases.length()
  let n = if n_orig > 255 { 255 } else { n_orig }
  let q = Array::make(m * m, 0.0)
  if n == 0 {
    return q
  }

  // Sort bases by value (descending, as in C code)
  let sorted = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    sorted[i] = bases[i]
  }
  sorted.sort_by(fn(a, b) { b - a })

  // Per-base counts and sums
  let fsum = Array::make(16, 0.0)
  let bsum = Array::make(16, 0.0)
  let cnt = Array::make(16, 0)
  let w = Array::make(32, 0)

  for j = n - 1; j >= 0; j = j - 1 {
    let b = sorted[j]
    let qual_raw = b >> 5
    let qual = if qual_raw < 4 { 4 } else { qual_raw }
    let qual_capped = if qual > 63 { 63 } else { qual }
    let basestrand = b & 0x1f
    let base = b & 0xf
    fsum[base] = fsum[base] + em.fk[w[basestrand]]
    bsum[base] = bsum[base] + em.fk[w[basestrand]] * errmod_beta(em, qual_capped, n, cnt[base])
    cnt[base] = cnt[base] + 1
    w[basestrand] = w[basestrand] + 1
  }

  // Generate likelihood
  for j = 0; j < m; j = j + 1 {
    // Homozygous
    let mut tmp1 = 0.0
    let mut tmp2 = 0
    for k = 0; k < m; k = k + 1 {
      if k == j {
        continue
      }
      tmp1 = tmp1 + bsum[k]
      tmp2 = tmp2 + cnt[k]
    }
    if tmp2 > 0 {
      q[j * m + j] = tmp1
    }
    // Heterozygous
    for k = j + 1; k < m; k = k + 1 {
      let cjk = cnt[j] + cnt[k]
      let mut tmp1_h = 0.0
      let mut tmp2_h = 0
      for i = 0; i < m; i = i + 1 {
        if i == j || i == k {
          continue
        }
        tmp1_h = tmp1_h + bsum[i]
        tmp2_h = tmp2_h + cnt[i]
      }
      let val = if tmp2_h > 0 {
        -4.343 * em.lhet[cjk * 256 + cnt[k]] + tmp1_h
      } else {
        -4.343 * em.lhet[cjk * 256 + cnt[k]]
      }
      q[j * m + k] = val
      q[k * m + j] = val
    }
    // Clamp to >= 0
    for k = 0; k < m; k = k + 1 {
      if q[j * m + k] < 0.0 {
        q[j * m + k] = 0.0
      }
    }
  }
  q
}

///|
/// Mapping quality capping (sam_cap_mapq from realn.c).
///
/// Given a read's alignment to reference, compute a capped mapping quality
/// based on the number of mismatches and their quality scores.
///
/// Parameters:
/// - ref: reference sequence string
/// - query: query sequence string
/// - qual: quality scores (phred scale)
/// - cigar_ops: CIGAR operations as array of (op, length) pairs
///   where op: 0=M, 1=I, 2=D, 4=S, 5=H
/// - pos: 0-based reference position of alignment start
/// - thres: threshold (default 40 if negative)
///
/// Returns: capped mapping quality, or -1 if below threshold
pub fn sam_cap_mapq(
  ref : String,
  query : String,
  qual : Array[Int],
  cigar_ops : Array[(Int, Int)],
  pos : Int,
  thres : Int,
) -> Int {
  let thres = if thres < 0 { 40 } else { thres }
  let ref_len = ref.length()
  let mut mm = 0
  let mut q_sum = 0
  let mut len = 0
  let mut clip_l = 0
  let mut clip_q = 0
  let mut x = pos
  let mut y = 0

  for (op, l) in cigar_ops {
    if op == 0 || op == 7 || op == 8 {
      // Match/mismatch
      for j = 0; j < l; j = j + 1 {
        if x + j >= ref_len {
          break
        }
        let ref_base = base_to_int(ref.unsafe_get(x + j).to_char().unwrap())
        let query_base = base_to_int(query.unsafe_get(y + j).to_char().unwrap())
        let qval = if y + j < qual.length() { qual[y + j] } else { 30 }
        if ref_base != 4 && query_base != 4 && qval >= 13 {
          len = len + 1
          if query_base != ref_base && qval >= 13 {
            mm = mm + 1
            q_sum = q_sum + (if qval > 33 { 33 } else { qval })
          }
        }
      }
      x = x + l
      y = y + l
      len = len + l
    } else if op == 2 {
      // Deletion
      x = x + l
    } else if op == 4 {
      // Soft clip
      for j = 0; j < l; j = j + 1 {
        clip_q = clip_q + (if y + j < qual.length() { qual[y + j] } else { 13 })
      }
      clip_l = clip_l + l
      y = y + l
    } else if op == 5 {
      // Hard clip
      clip_q = clip_q + 13 * l
      clip_l = clip_l + l
    } else if op == 1 {
      // Insertion
      y = y + l
    }
  }

  // Compute t = q - 4.343 * log(prod(len/(i+1))) + clip_q/5
  let mut t = 1.0
  for i = 0; i < mm; i = i + 1 {
    t = t * len.to_double() / (i + 1).to_double()
  }
  let t_val = q_sum.to_double() - 4.343 * @math.ln(t) + clip_q.to_double() / 5.0
  if t_val > thres.to_double() {
    return -1
  }
  let t_val = if t_val < 0.0 { 0.0 } else { t_val }
  let result = ((thres.to_double() - t_val) / thres.to_double()).sqrt() * thres.to_double()
  (result + 0.499).to_int()
}
