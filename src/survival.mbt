///|
/// survival.mbt - Survival Analysis module for BioSeqs.
/// Inspired by Bioconductor's survival package.
/// Provides Kaplan-Meier estimation, the log-rank test, and Cox proportional
/// hazards model fitting.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A single survival observation record.
pub struct SurvRecord {
  time : Double
  event : Bool
  group : String
  covariates : Array[Double]
}

///|
/// A collection of survival records grouped for analysis.
pub struct SurvData {
  records : Array[SurvRecord]
  mut n_samples : Int
  mut n_groups : Int
  group_names : Array[String]
}

///|
/// Result of a Kaplan-Meier survival estimate.
pub struct KaplanMeierResult {
  times : Array[Double]
  survival : Array[Double]
  std_err : Array[Double]
  n_at_risk : Array[Int]
  n_events : Array[Int]
  group : String
}

///|
/// Result of a two-sample log-rank test.
pub struct LogRankResult {
  chi_square : Double
  df : Int
  p_value : Double
  group1 : String
  group2 : String
  n1 : Int
  n2 : Int
}

///|
/// Result of a Cox proportional hazards model fit.
pub struct CoxResult {
  coefficients : Array[Double]
  std_errors : Array[Double]
  hazard_ratios : Array[Double]
  z_values : Array[Double]
  p_values : Array[Double]
  covariate_names : Array[String]
  log_likelihood : Double
  n_samples : Int
  converged : Bool
}

// ===========================================================================
// SurvRecord methods
// ===========================================================================

///|
/// Create a new survival record.
pub fn SurvRecord::new(
  time : Double,
  event : Bool,
  group : String,
  covariates : Array[Double]
) -> SurvRecord {
  SurvRecord::{ time, event, group, covariates }
}

///|
/// Get the observed time.
pub fn SurvRecord::time(self : SurvRecord) -> Double {
  self.time
}

///|
/// Whether an event (e.g. death) was observed.
pub fn SurvRecord::event(self : SurvRecord) -> Bool {
  self.event
}

///|
/// Get the group label.
pub fn SurvRecord::group(self : SurvRecord) -> String {
  self.group
}

///|
/// Get the covariate vector.
pub fn SurvRecord::covariates(self : SurvRecord) -> Array[Double] {
  self.covariates
}

// ===========================================================================
// SurvData methods
// ===========================================================================

///|
/// Create an empty SurvData container.
pub fn SurvData::new() -> SurvData {
  SurvData::{
    records: Array::new(),
    n_samples: 0,
    n_groups: 0,
    group_names: Array::new()
  }
}

///|
/// Build a SurvData from a set of records, computing group metadata.
pub fn SurvData::from_records(records : Array[SurvRecord]) -> SurvData {
  let n = records.length()
  let group_names : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let g = records[i].group
    let mut found = false
    let mut j = 0
    while j < group_names.length() {
      if group_names[j] == g {
        found = true
        break
      }
      j = j + 1
    }
    if !found {
      group_names.push(g)
    }
    i = i + 1
  }
  SurvData::{
    records,
    n_samples: n,
    n_groups: group_names.length(),
    group_names
  }
}

///|
/// Number of samples.
pub fn SurvData::n_samples(self : SurvData) -> Int {
  self.n_samples
}

///|
/// Number of groups.
pub fn SurvData::n_groups(self : SurvData) -> Int {
  self.n_groups
}

///|
/// Group names in order of first appearance.
pub fn SurvData::group_names(self : SurvData) -> Array[String] {
  self.group_names
}

///|
/// Underlying records.
pub fn SurvData::records(self : SurvData) -> Array[SurvRecord] {
  self.records
}

///|
/// Append a record, updating derived counts.
pub fn SurvData::add_record(self : SurvData, record : SurvRecord) -> Unit {
  self.records.push(record)
  self.n_samples = self.n_samples + 1
  let g = record.group
  let mut found = false
  let mut j = 0
  while j < self.group_names.length() {
    if self.group_names[j] == g {
      found = true
      break
    }
    j = j + 1
  }
  if !found {
    self.group_names.push(g)
    self.n_groups = self.n_groups + 1
  }
}

///|
/// Return records belonging to a given group.
pub fn SurvData::filter_by_group(self : SurvData, group : String) -> Array[SurvRecord] {
  let result : Array[SurvRecord] = Array::new()
  let recs = self.records
  let n = recs.length()
  let mut i = 0
  while i < n {
    if recs[i].group == group {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

// ===========================================================================
// Math helpers
// ===========================================================================

///|
/// Wrapper for @math.exp.
pub fn survival_exp(x : Double) -> Double {
  @math.exp(x)
}

///|
/// Wrapper for @math.ln (safe against non-positive input).
pub fn survival_ln(x : Double) -> Double {
  let v = if x > 0.0 { x } else { 1.0e-300 }
  @math.ln(v)
}

///|
/// Wrapper for absolute value.
pub fn survival_abs(x : Double) -> Double {
  x.abs()
}

///|
/// Standard normal CDF (Abramowitz and Stegun approximation).
pub fn survival_normal_cdf(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  let t = 1.0 / (1.0 + 0.2316419 * abs_z)
  let d = 0.3989423 * @math.exp(-z * z / 2.0)
  let p = d *
    t *
    (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.3302744))))
  if z > 0.0 { 1.0 - p } else { p }
}

///|
/// Two-sided p-value from a standard normal variate.
pub fn survival_normal_p_value_two_sided(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  let p = 2.0 * (1.0 - survival_normal_cdf(abs_z))
  if p.is_nan() || p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Upper-tail p-value for a chi-square statistic with given degrees of freedom.
pub fn survival_chi_square_p_value(chi_sq : Double, df : Int) -> Double {
  if chi_sq <= 0.0 {
    return 1.0
  }
  if df <= 0 {
    return 1.0
  }
  if df == 1 {
    let p = 2.0 * (1.0 - survival_normal_cdf(chi_sq.sqrt()))
    return if p.is_nan() || p < 0.0 {
      0.0
    } else if p > 1.0 {
      1.0
    } else {
      p
    }
  }
  // Wilson-Hilferty approximation for general df
  let df_d = df.to_double()
  let z = (@math.pow(chi_sq / df_d, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_d))) /
    ((2.0 / (9.0 * df_d)).sqrt())
  let p = 1.0 - survival_normal_cdf(z)
  if p.is_nan() || p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

// ===========================================================================
// Linear algebra helpers (internal)
// ===========================================================================

///|
/// Solve the linear system a * x = b via Gauss-Jordan elimination.
fn survival_solve_linear(
  a : Array[Array[Double]],
  b : Array[Double],
  n : Int
) -> Array[Double] {
  if n == 0 {
    return Array::new()
  }
  let m : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(n + 1, 0.0)
    let mut j = 0
    while j < n {
      row[j] = a[i][j]
      j = j + 1
    }
    row[n] = b[i]
    m.push(row)
    i = i + 1
  }
  let mut col = 0
  while col < n {
    let mut pivot = col
    let mut r = col + 1
    while r < n {
      if m[r][col].abs() > m[pivot][col].abs() {
        pivot = r
      }
      r = r + 1
    }
    if pivot != col {
      let tmp = m[pivot]
      m[pivot] = m[col]
      m[col] = tmp
    }
    let piv_val = m[col][col]
    if piv_val.abs() >= 1.0e-30 {
      let mut rr = 0
      while rr < n {
        if rr != col {
          let factor = m[rr][col] / piv_val
          let mut jj = 0
          while jj <= n {
            m[rr][jj] = m[rr][jj] - factor * m[col][jj]
            jj = jj + 1
          }
        }
        rr = rr + 1
      }
    }
    col = col + 1
  }
  let x : Array[Double] = Array::make(n, 0.0)
  let mut k = 0
  while k < n {
    let piv = m[k][k]
    x[k] = if piv.abs() < 1.0e-30 { 0.0 } else { m[k][n] / piv }
    k = k + 1
  }
  x
}

///|
/// Invert a square matrix via Gauss-Jordan elimination with partial pivoting.
fn survival_invert_matrix(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = m.length()
  if n == 0 {
    return Array::new()
  }
  let aug : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(2 * n, 0.0)
    let mut j = 0
    while j < n {
      row[j] = m[i][j]
      j = j + 1
    }
    row[n + i] = 1.0
    aug.push(row)
    i = i + 1
  }
  let two_n = 2 * n
  let mut col = 0
  while col < n {
    let mut pivot = col
    let mut r = col + 1
    while r < n {
      if aug[r][col].abs() > aug[pivot][col].abs() {
        pivot = r
      }
      r = r + 1
    }
    if pivot != col {
      let tmp = aug[pivot]
      aug[pivot] = aug[col]
      aug[col] = tmp
    }
    let raw_piv = aug[col][col]
    let piv_val = if raw_piv.abs() < 1.0e-12 {
      if raw_piv >= 0.0 { 1.0e-12 } else { -1.0e-12 }
    } else {
      raw_piv
    }
    aug[col][col] = piv_val
    let mut j = 0
    while j < two_n {
      aug[col][j] = aug[col][j] / piv_val
      j = j + 1
    }
    let mut rr = 0
    while rr < n {
      if rr != col {
        let factor = aug[rr][col]
        let mut jj = 0
        while jj < two_n {
          aug[rr][jj] = aug[rr][jj] - factor * aug[col][jj]
          jj = jj + 1
        }
      }
      rr = rr + 1
    }
    col = col + 1
  }
  let inv : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(aug[i][n + j])
      j = j + 1
    }
    inv.push(row)
    i = i + 1
  }
  inv
}

///|
/// Format a double for summary output (4 decimal places).
fn survival_fmt(v : Double) -> String {
  if v.is_nan() {
    return "NA"
  }
  if v == @double.infinity || v == @double.neg_infinity {
    return "Inf"
  }
  let r = (v * 10000.0).round() / 10000.0
  "\{r}"
}

// ===========================================================================
// Kaplan-Meier estimator
// ===========================================================================

///|
/// Compute the Kaplan-Meier survival curve for a specific group.
/// Standard error uses Greenwood's formula.
pub fn survival_kaplan_meier(data : SurvData, group : String) -> KaplanMeierResult {
  let group_records = data.filter_by_group(group)
  let n = group_records.length()

  let times : Array[Double] = Array::new()
  let survival_arr : Array[Double] = Array::new()
  let std_err_arr : Array[Double] = Array::new()
  let n_at_risk_arr : Array[Int] = Array::new()
  let n_events_arr : Array[Int] = Array::new()

  // Start with S(0) = 1.0
  times.push(0.0)
  survival_arr.push(1.0)
  std_err_arr.push(0.0)
  n_at_risk_arr.push(n)
  n_events_arr.push(0)

  if n == 0 {
    return KaplanMeierResult::{
      times,
      survival: survival_arr,
      std_err: std_err_arr,
      n_at_risk: n_at_risk_arr,
      n_events: n_events_arr,
      group
    }
  }

  // Sort records by time (copy to avoid mutating input).
  let sorted = group_records.copy()
  sorted.sort_by(fn(a : SurvRecord, b : SurvRecord) -> Int {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      0
    }
  })

  let mut i = 0
  let mut surv = 1.0
  let mut cum_var = 0.0
  while i < n {
    let t = sorted[i].time
    let mut d = 0
    let mut j = i
    while j < n && sorted[j].time == t {
      if sorted[j].event {
        d = d + 1
      }
      j = j + 1
    }
    let n_at_risk = n - i
    if d > 0 && n_at_risk > 0 {
      surv = surv * (1.0 - d.to_double() / n_at_risk.to_double())
      if n_at_risk > d {
        cum_var = cum_var +
          d.to_double() /
          (n_at_risk.to_double() * (n_at_risk - d).to_double())
      }
      let se = surv * cum_var.sqrt()
      times.push(t)
      survival_arr.push(surv)
      std_err_arr.push(se)
      n_at_risk_arr.push(n_at_risk)
      n_events_arr.push(d)
    }
    i = j
  }

  KaplanMeierResult::{
    times,
    survival: survival_arr,
    std_err: std_err_arr,
    n_at_risk: n_at_risk_arr,
    n_events: n_events_arr,
    group
  }
}

// ===========================================================================
// Log-rank test
// ===========================================================================

///|
/// Two-sample log-rank test comparing survival curves of two groups.
pub fn survival_log_rank(
  data : SurvData,
  group1 : String,
  group2 : String
) -> LogRankResult {
  let r1 = data.filter_by_group(group1)
  let r2 = data.filter_by_group(group2)
  let n1_total = r1.length()
  let n2_total = r2.length()

  // (time, event_flag, group_idx)
  let combined : Array[(Double, Int, Int)] = Array::new()
  for r in r1 {
    combined.push((r.time, if r.event { 1 } else { 0 }, 0))
  }
  for r in r2 {
    combined.push((r.time, if r.event { 1 } else { 0 }, 1))
  }
  combined.sort_by(fn(a : (Double, Int, Int), b : (Double, Int, Int)) -> Int {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })

  let mut o1 = 0.0
  let mut o2 = 0.0
  let mut e1 = 0.0
  let mut e2 = 0.0
  let mut v = 0.0
  let mut n1 = n1_total.to_double()
  let mut n2 = n2_total.to_double()

  let n_total = combined.length()
  let mut i = 0
  while i < n_total {
    let t = combined[i].0
    let mut d1 = 0
    let mut d2 = 0
    let mut cnt1 = 0
    let mut cnt2 = 0
    let mut j = i
    while j < n_total && combined[j].0 == t {
      if combined[j].1 == 1 {
        if combined[j].2 == 0 {
          d1 = d1 + 1
        } else {
          d2 = d2 + 1
        }
      }
      if combined[j].2 == 0 {
        cnt1 = cnt1 + 1
      } else {
        cnt2 = cnt2 + 1
      }
      j = j + 1
    }
    let d = d1 + d2
    let at_risk = n1 + n2
    if d > 0 && at_risk > 0.0 {
      e1 = e1 + d.to_double() * n1 / at_risk
      e2 = e2 + d.to_double() * n2 / at_risk
      if at_risk > 1.0 {
        v = v +
          d.to_double() *
          n1 *
          n2 *
          (at_risk - d.to_double()) /
          (at_risk * at_risk * (at_risk - 1.0))
      }
      o1 = o1 + d1.to_double()
      o2 = o2 + d2.to_double()
    }
    n1 = n1 - cnt1.to_double()
    n2 = n2 - cnt2.to_double()
    i = j
  }

  let chi_sq = if v > 0.0 {
    (o1 - e1) * (o1 - e1) / v
  } else {
    0.0
  }
  let p = survival_chi_square_p_value(chi_sq, 1)
  LogRankResult::{
    chi_square: chi_sq,
    df: 1,
    p_value: p,
    group1,
    group2,
    n1: n1_total,
    n2: n2_total
  }
}

// ===========================================================================
// Cox proportional hazards model
// ===========================================================================

///|
/// Compute the partial-likelihood score, information matrix and log-likelihood
/// for a Cox model at the given coefficient vector. Ties use Breslow's method.
fn survival_cox_score_info(
  times : Array[Double],
  events : Array[Int],
  z : Array[Array[Double]],
  order : Array[Int],
  beta : Array[Double],
  n : Int,
  p : Int
) -> (Array[Double], Array[Array[Double]], Double) {
  let u : Array[Double] = Array::make(p, 0.0)
  let info : Array[Array[Double]] = Array::new()
  let mut kk = 0
  while kk < p {
    info.push(Array::make(p, 0.0))
    kk = kk + 1
  }
  let mut ll = 0.0

  let mut idx = 0
  while idx < n {
    let subj = order[idx]
    let t = times[subj]

    // Build the risk set (subjects with time >= t): order[idx .. n-1].
    let mut s_w = 0.0
    let sz : Array[Double] = Array::make(p, 0.0)
    let szz : Array[Array[Double]] = Array::new()
    let mut qi = 0
    while qi < p {
      szz.push(Array::make(p, 0.0))
      qi = qi + 1
    }
    let mut r = idx
    while r < n {
      let j = order[r]
      let mut dot = 0.0
      let mut c = 0
      while c < p {
        dot = dot + beta[c] * z[j][c]
        c = c + 1
      }
      let w = survival_exp(dot)
      s_w = s_w + w
      let mut c2 = 0
      while c2 < p {
        sz[c2] = sz[c2] + z[j][c2] * w
        let mut l = 0
        while l < p {
          szz[c2][l] = szz[c2][l] + z[j][c2] * z[j][l] * w
          l = l + 1
        }
        c2 = c2 + 1
      }
      r = r + 1
    }

    // Subjects with an event at this time.
    let event_subjs : Array[Int] = Array::new()
    let mut d = 0
    let mut m = idx
    while m < n && times[order[m]] == t {
      if events[order[m]] == 1 {
        event_subjs.push(order[m])
        d = d + 1
      }
      m = m + 1
    }

    if d > 0 && s_w > 0.0 {
      let ln_s_w = survival_ln(s_w)
      let mut ei = 0
      while ei < d {
        let ev = event_subjs[ei]
        let mut dot = 0.0
        let mut c = 0
        while c < p {
          dot = dot + beta[c] * z[ev][c]
          c = c + 1
        }
        ll = ll + dot - ln_s_w
        let mut c2 = 0
        while c2 < p {
          u[c2] = u[c2] + z[ev][c2] - sz[c2] / s_w
          c2 = c2 + 1
        }
        ei = ei + 1
      }
      let mut c = 0
      while c < p {
        let mut l = 0
        while l < p {
          let term = szz[c][l] / s_w - sz[c] * sz[l] / (s_w * s_w)
          info[c][l] = info[c][l] + d.to_double() * term
          l = l + 1
        }
        c = c + 1
      }
    }
    idx = m
  }

  (u, info, ll)
}

///|
/// Fit a Cox proportional hazards model using Newton-Raphson iteration on the
/// partial likelihood. Ties are handled with the Breslow approximation.
pub fn survival_cox_fit(
  data : SurvData,
  covariate_names : Array[String]
) -> CoxResult {
  let records = data.records()
  let n = records.length()
  let p = covariate_names.length()

  if n == 0 || p == 0 {
    return CoxResult::{
      coefficients: Array::new(),
      std_errors: Array::new(),
      hazard_ratios: Array::new(),
      z_values: Array::new(),
      p_values: Array::new(),
      covariate_names,
      log_likelihood: 0.0,
      n_samples: n,
      converged: true
    }
  }

  let times : Array[Double] = Array::make(n, 0.0)
  let events : Array[Int] = Array::make(n, 0)
  let z : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    times[i] = records[i].time
    events[i] = if records[i].event { 1 } else { 0 }
    z.push(records[i].covariates)
    i = i + 1
  }

  // Sort subject indices by time.
  let order : Array[Int] = Array::make(n, 0)
  i = 0
  while i < n {
    order[i] = i
    i = i + 1
  }
  order.sort_by(fn(a : Int, b : Int) -> Int {
    if times[a] < times[b] {
      -1
    } else if times[a] > times[b] {
      1
    } else {
      0
    }
  })

  let beta : Array[Double] = Array::make(p, 0.0)

  let mut iter = 0
  let mut converged = false
  let mut cont = true
  while cont && iter < 100 {
    let triple = survival_cox_score_info(times, events, z, order, beta, n, p)
    let u = triple.0
    let info = triple.1
    let delta = survival_solve_linear(info, u, p)
    let mut max_change = 0.0
    let mut k = 0
    while k < p {
      beta[k] = beta[k] + delta[k]
      let ch = delta[k].abs()
      if ch > max_change {
        max_change = ch
      }
      k = k + 1
    }
    iter = iter + 1
    if max_change < 1.0e-6 {
      converged = true
      cont = false
    }
  }

  // Final score/info/loglik at the converged beta.
  let final_triple = survival_cox_score_info(times, events, z, order, beta, n, p)
  let final_info = final_triple.1
  let final_ll = final_triple.2
  let inv_info = survival_invert_matrix(final_info)

  let std_errors : Array[Double] = Array::make(p, 0.0)
  let hazard_ratios : Array[Double] = Array::make(p, 0.0)
  let z_values : Array[Double] = Array::make(p, 0.0)
  let p_values : Array[Double] = Array::make(p, 0.0)
  let mut k = 0
  while k < p {
    let var_k = inv_info[k][k]
    let se = if var_k > 0.0 { var_k.sqrt() } else { 0.0 }
    std_errors[k] = se
    hazard_ratios[k] = survival_exp(beta[k])
    z_values[k] = if se > 0.0 { beta[k] / se } else { 0.0 }
    p_values[k] = survival_normal_p_value_two_sided(z_values[k])
    k = k + 1
  }

  CoxResult::{
    coefficients: beta,
    std_errors,
    hazard_ratios,
    z_values,
    p_values,
    covariate_names,
    log_likelihood: final_ll,
    n_samples: n,
    converged
  }
}

// ===========================================================================
// Summary and utility functions
// ===========================================================================

///|
/// Formatted summary of a Kaplan-Meier estimate.
pub fn survival_km_summary(result : KaplanMeierResult) -> String {
  let mut s = "Kaplan-Meier Survival Estimate\n"
  s = s + "Group: " + result.group + "\n"
  s = s + "Time points: " + result.times.length().to_string() + "\n"
  let med = survival_median_survival(result)
  let med_str = match med {
    Some(t) => survival_fmt(t)
    None => "NA"
  }
  s = s + "Median survival: " + med_str + "\n"
  s
}

///|
/// Formatted summary of a log-rank test.
pub fn survival_logrank_summary(result : LogRankResult) -> String {
  let mut s = "Log-Rank Test\n"
  s = s + "Group 1: " + result.group1 + " (n=" + result.n1.to_string() + ")\n"
  s = s + "Group 2: " + result.group2 + " (n=" + result.n2.to_string() + ")\n"
  s = s +
    "Chi-square = " +
    survival_fmt(result.chi_square) +
    ", df = " +
    result.df.to_string() +
    "\n"
  s = s + "p-value = " + survival_fmt(result.p_value) + "\n"
  s
}

///|
/// Formatted summary of a Cox model fit.
pub fn survival_cox_summary(result : CoxResult) -> String {
  let conv_str = if result.converged { "true" } else { "false" }
  let mut s = "Cox Proportional Hazards Model\n"
  s = s +
    "n = " +
    result.n_samples.to_string() +
    ", converged = " +
    conv_str +
    "\n"
  s = s + "Log-likelihood = " + survival_fmt(result.log_likelihood) + "\n"
  s = s + "Covariate\tcoef\texp(coef)\tse(coef)\tz\tPr(>|z|)\n"
  let p = result.covariate_names.length()
  let mut i = 0
  while i < p {
    s = s +
      result.covariate_names[i] +
      "\t" +
      survival_fmt(result.coefficients[i]) +
      "\t" +
      survival_fmt(result.hazard_ratios[i]) +
      "\t" +
      survival_fmt(result.std_errors[i]) +
      "\t" +
      survival_fmt(result.z_values[i]) +
      "\t" +
      survival_fmt(result.p_values[i]) +
      "\n"
    i = i + 1
  }
  s
}

///|
/// Time at which the survival curve first drops to or below 0.5.
pub fn survival_median_survival(result : KaplanMeierResult) -> Double? {
  let n = result.times.length()
  let mut i = 0
  while i < n {
    if result.survival[i] <= 0.5 {
      return Some(result.times[i])
    }
    i = i + 1
  }
  None
}

///|
/// Tabular string representation of a Kaplan-Meier estimate.
pub fn survival_km_to_string(result : KaplanMeierResult) -> String {
  let mut s = "Kaplan-Meier estimate for group: " + result.group + "\n"
  s = s + "time\tsurvival\tstd_err\tn_at_risk\tn_events\n"
  let n = result.times.length()
  let mut i = 0
  while i < n {
    s = s +
      survival_fmt(result.times[i]) +
      "\t" +
      survival_fmt(result.survival[i]) +
      "\t" +
      survival_fmt(result.std_err[i]) +
      "\t" +
      result.n_at_risk[i].to_string() +
      "\t" +
      result.n_events[i].to_string() +
      "\n"
    i = i + 1
  }
  s
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Sample survival data with two groups (control / treatment) and two
/// covariates per record (age, treatment indicator). 15 records total.
pub fn survival_sample_data() -> SurvData {
  let records : Array[SurvRecord] = Array::new()
  // Control group: covariates = [age, treatment_indicator]
  records.push(SurvRecord::new(6.0, true, "control", [55.0, 0.0]))
  records.push(SurvRecord::new(13.0, true, "control", [62.0, 0.0]))
  records.push(SurvRecord::new(21.0, true, "control", [48.0, 0.0]))
  records.push(SurvRecord::new(30.0, true, "control", [70.0, 0.0]))
  records.push(SurvRecord::new(37.0, true, "control", [59.0, 0.0]))
  records.push(SurvRecord::new(45.0, false, "control", [65.0, 0.0]))
  records.push(SurvRecord::new(46.0, true, "control", [53.0, 0.0]))
  // Treatment group: covariates = [age, treatment_indicator]
  records.push(SurvRecord::new(9.0, true, "treatment", [58.0, 1.0]))
  records.push(SurvRecord::new(18.0, true, "treatment", [64.0, 1.0]))
  records.push(SurvRecord::new(22.0, false, "treatment", [50.0, 1.0]))
  records.push(SurvRecord::new(26.0, true, "treatment", [60.0, 1.0]))
  records.push(SurvRecord::new(33.0, true, "treatment", [47.0, 1.0]))
  records.push(SurvRecord::new(50.0, false, "treatment", [68.0, 1.0]))
  records.push(SurvRecord::new(60.0, true, "treatment", [55.0, 1.0]))
  records.push(SurvRecord::new(80.0, false, "treatment", [72.0, 1.0]))
  SurvData::from_records(records)
}
