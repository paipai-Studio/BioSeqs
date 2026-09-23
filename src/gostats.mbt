///|
/// GOstats - Gene Ontology enrichment analysis using the hypergeometric test.
/// Inspired by Bioconductor GOstats package.
/// Provides hypergeometric test, conditional test, and Benjamini-Hochberg FDR correction.

///|
pub struct GOStatsTerm {
  go_id : String
  name : String
  ontology : String
  parent_ids : Array[String]
  gene_ids : Array[String]
}

///|
pub fn GOStatsTerm::new(
  go_id : String,
  name : String,
  ontology : String,
  parent_ids : Array[String],
  gene_ids : Array[String],
) -> GOStatsTerm {
  GOStatsTerm::{ go_id, name, ontology, parent_ids, gene_ids, }
}

///|
pub struct GOStatsGraph {
  terms : Map[String, GOStatsTerm]
}

///|
pub fn GOStatsGraph::new(terms : Map[String, GOStatsTerm]) -> GOStatsGraph {
  GOStatsGraph::{ terms, }
}

///|
pub fn GOStatsGraph::add_term(
  self : GOStatsGraph,
  term : GOStatsTerm,
) -> GOStatsGraph {
  let terms : Map[String, GOStatsTerm] = self.terms.copy()
  terms.set(term.go_id, term)
  GOStatsGraph::{ terms, }
}

///|
pub fn GOStatsGraph::get_term(
  self : GOStatsGraph,
  go_id : String,
) -> GOStatsTerm? {
  self.terms.get(go_id)
}

///|
pub fn GOStatsGraph::has_term(self : GOStatsGraph, go_id : String) -> Bool {
  self.terms.contains(go_id)
}

///|
pub struct GOHyperGeoResult {
  go_id : String
  name : String
  ontology : String
  p_value : Double
  fdr : Double
  odds_ratio : Double
  expected : Double
  observed : Int
  gene_count : Int
  universe_count : Int
}

///|
pub fn GOHyperGeoResult::new(
  go_id : String,
  name : String,
  ontology : String,
  p_value : Double,
  fdr : Double,
  odds_ratio : Double,
  expected : Double,
  observed : Int,
  gene_count : Int,
  universe_count : Int,
) -> GOHyperGeoResult {
  GOHyperGeoResult::{
    go_id,
    name,
    ontology,
    p_value,
    fdr,
    odds_ratio,
    expected,
    observed,
    gene_count,
    universe_count,
  }
}

///|
pub struct GOHyperGeoTest {
  universe : Array[String]
  selected : Array[String]
  go_graph : GOStatsGraph
}

///|
pub fn GOHyperGeoTest::new(
  universe : Array[String],
  selected : Array[String],
  go_graph : GOStatsGraph,
) -> GOHyperGeoTest {
  GOHyperGeoTest::{ universe, selected, go_graph, }
}

///|
/// Stirling's approximation for lgamma(x):
/// lgamma(x) ≈ 0.5*ln(2π) + (x-0.5)*ln(x) - x + 1/(12x) - 1/(360x^3)
fn gostats_lgamma(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x < 7.0 {
    let mut shift : Double = 0.0
    let mut xx : Double = x
    while xx < 7.0 {
      shift = shift - @math.ln(xx)
      xx = xx + 1.0
    }
    shift + gostats_lgamma_stirling(xx)
  } else {
    gostats_lgamma_stirling(x)
  }
}

///|
fn gostats_lgamma_stirling(x : Double) -> Double {
  let ln_2pi : Double = 0.9189385332046727
  let inv_12x : Double = 1.0 / (12.0 * x)
  let inv_360x3 : Double = 1.0 / (360.0 * x * x * x)
  ln_2pi + (x - 0.5) * @math.ln(x) - x + inv_12x - inv_360x3
}

///|
/// Log of binomial coefficient using lgamma:
/// log C(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
pub fn gostats_log_choose(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return 0.0
  }
  if k == 0 || k == n {
    return 0.0
  }
  let nk : Int = n - k
  gostats_lgamma(n.to_double() + 1.0) -
  gostats_lgamma(k.to_double() + 1.0) -
  gostats_lgamma(nk.to_double() + 1.0)
}

///|
/// Hypergeometric p-value: P(X >= k) = sum_{i=k}^{min(m,n)} C(m,i)*C(N-m,n-i)/C(N,n)
/// N = universe size, m = genes in GO term, n = selected genes, k = overlap
fn gostats_hypergeo_pvalue(n_univ : Int, m : Int, n : Int, k : Int) -> Double {
  let upper : Int = if m < n { m } else { n }
  if k > upper {
    return 1.0
  }
  if k <= 0 {
    return 1.0
  }

  let log_denom : Double = gostats_log_choose(n_univ, n)
  let mut log_max : Double = -1.0e300
  let mut ii : Int = k
  while ii <= upper {
    let log_term : Double = gostats_log_choose(m, ii) +
      gostats_log_choose(n_univ - m, n - ii) -
      log_denom
    if log_term > log_max {
      log_max = log_term
    }
    ii = ii + 1
  }

  let mut sum : Double = 0.0
  ii = k
  while ii <= upper {
    let log_term : Double = gostats_log_choose(m, ii) +
      gostats_log_choose(n_univ - m, n - ii) -
      log_denom
    sum = sum + @math.exp(log_term - log_max)
    ii = ii + 1
  }

  let pval : Double = @math.exp(log_max) * sum
  if pval > 1.0 {
    1.0
  } else {
    pval
  }
}

///|
/// Compute odds ratio = (k * (N-m-n+k)) / ((m-k) * (n-k))
fn gostats_odds_ratio(n_univ : Int, m : Int, n : Int, k : Int) -> Double {
  let aa : Int = k
  let bb : Int = m - k
  let cc : Int = n - k
  let dd : Int = n_univ - m - n + k
  if aa > 0 && bb > 0 && cc > 0 && dd > 0 {
    aa.to_double() * dd.to_double() / (bb.to_double() * cc.to_double())
  } else {
    0.0
  }
}

///|
/// Benjamini-Hochberg FDR correction for multiple testing.
pub fn gostats_bh_fdr(p_values : Array[Double]) -> Array[Double] {
  let nn : Int = p_values.length()
  if nn == 0 {
    return Array(capacity=0)
  }

  let sorted_with_idx : Array[(Double, Int)] = Array(capacity=0)
  let mut i : Int = 0
  while i < nn {
    sorted_with_idx.push((p_values[i], i))
    i = i + 1
  }

  let mut j : Int = 0
  while j < nn - 1 {
    let mut kk : Int = j + 1
    while kk < nn {
      if sorted_with_idx[j].0 > sorted_with_idx[kk].0 {
        let tmp : (Double, Int) = sorted_with_idx[j]
        sorted_with_idx[j] = sorted_with_idx[kk]
        sorted_with_idx[kk] = tmp
      }
      kk = kk + 1
    }
    j = j + 1
  }

  let corrected_sorted : Array[Double] = Array(capacity=0)
  i = 0
  while i < nn {
    let rank : Double = (i + 1).to_double()
    let pp : Double = sorted_with_idx[i].0
    let corrected : Double = pp * nn.to_double() / rank
    corrected_sorted.push(if corrected > 1.0 { 1.0 } else { corrected })
    i = i + 1
  }

  for i = nn - 2; i >= 0; i = i - 1 {
    if corrected_sorted[i] > corrected_sorted[i + 1] {
      corrected_sorted[i] = corrected_sorted[i + 1]
    }
  }

  let bh_result : Array[Double] = Array::make(nn, 0.0)
  i = 0
  while i < nn {
    let (_, original_idx) : (Double, Int) = sorted_with_idx[i]
    bh_result[original_idx] = corrected_sorted[i]
    i = i + 1
  }

  bh_result
}

///|
/// Build a set from an array of strings for fast membership testing.
fn gostats_string_set(arr : Array[String]) -> Map[String, Bool] {
  let s : Map[String, Bool] = Map([], capacity=arr.length())
  let mut i : Int = 0
  while i < arr.length() {
    s.set(arr[i], true)
    i = i + 1
  }
  s
}

///|
/// Compute intersection count between a string array and a string set.
fn gostats_intersection_count(
  arr : Array[String],
  set : Map[String, Bool],
) -> Int {
  let mut count : Int = 0
  let mut i : Int = 0
  while i < arr.length() {
    if set.contains(arr[i]) {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Compute intersection genes between a string array and a string set.
fn gostats_intersection_genes(
  arr : Array[String],
  set : Map[String, Bool],
) -> Array[String] {
  let inter_result : Array[String] = Array(capacity=0)
  let mut i : Int = 0
  while i < arr.length() {
    if set.contains(arr[i]) {
      inter_result.push(arr[i])
    }
    i = i + 1
  }
  inter_result
}

///|
/// GOstats hypergeometric test for all GO terms.
/// Optionally filter by ontology ("BP", "MF", "CC"). Empty string means all.
pub fn gostats_hypergeometric_test(
  hg_test : GOHyperGeoTest,
  ontology? : String = "",
) -> Array[GOHyperGeoResult] {
  let universe_set : Map[String, Bool] = gostats_string_set(hg_test.universe)
  let selected_set : Map[String, Bool] = gostats_string_set(hg_test.selected)
  let n_univ : Int = hg_test.universe.length()
  let n_sel : Int = hg_test.selected.length()

  let hyp_results : Array[GOHyperGeoResult] = Array(capacity=0)
  let hyp_p_values : Array[Double] = Array(capacity=0)

  let term_keys : Array[String] = hg_test.go_graph.terms.keys().collect()
  let mut i : Int = 0
  while i < term_keys.length() {
    let go_id : String = term_keys[i]
    let gterm : GOStatsTerm = hg_test.go_graph.terms.get(go_id).unwrap()

    if ontology != "" && gterm.ontology != ontology {
      i = i + 1
      continue
    }

    let term_genes_in_univ : Array[String] = gostats_intersection_genes(
      gterm.gene_ids,
      universe_set,
    )
    let m : Int = term_genes_in_univ.length()
    let k_overlap : Int = gostats_intersection_count(
      term_genes_in_univ, selected_set,
    )

    if k_overlap == 0 {
      i = i + 1
      continue
    }

    let p_value : Double = gostats_hypergeo_pvalue(n_univ, m, n_sel, k_overlap)
    hyp_p_values.push(p_value)

    let odds_val : Double = gostats_odds_ratio(n_univ, m, n_sel, k_overlap)
    let expected_val : Double = m.to_double() *
      n_sel.to_double() /
      n_univ.to_double()

    hyp_results.push(
      GOHyperGeoResult::new(
        gterm.go_id,
        gterm.name,
        gterm.ontology,
        p_value,
        0.0,
        odds_val,
        expected_val,
        k_overlap,
        m,
        n_univ,
      ),
    )

    i = i + 1
  }

  let fdr_vals : Array[Double] = gostats_bh_fdr(hyp_p_values)
  let mut j : Int = 0
  while j < hyp_results.length() {
    hyp_results[j] = GOHyperGeoResult::{
      go_id: hyp_results[j].go_id,
      name: hyp_results[j].name,
      ontology: hyp_results[j].ontology,
      p_value: hyp_results[j].p_value,
      fdr: fdr_vals[j],
      odds_ratio: hyp_results[j].odds_ratio,
      expected: hyp_results[j].expected,
      observed: hyp_results[j].observed,
      gene_count: hyp_results[j].gene_count,
      universe_count: hyp_results[j].universe_count,
    }
    j = j + 1
  }

  let mut a : Int = 0
  while a < hyp_results.length() - 1 {
    let mut b : Int = a + 1
    while b < hyp_results.length() {
      if hyp_results[a].p_value > hyp_results[b].p_value {
        let tmp : GOHyperGeoResult = hyp_results[a]
        hyp_results[a] = hyp_results[b]
        hyp_results[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }

  hyp_results
}

///|
/// Compute the depth (distance from root) of a GO term in the graph via parent traversal.
fn gostats_term_depth(go_graph : GOStatsGraph, go_id : String) -> Int {
  let visited : Map[String, Bool] = Map([], capacity=50)
  visited.set(go_id, true)
  let depth_queue : Array[(String, Int)] = Array(capacity=0)
  depth_queue.push((go_id, 0))
  let mut max_depth : Int = 0
  let mut idx : Int = 0
  while idx < depth_queue.length() {
    let (cid, depth) : (String, Int) = depth_queue[idx]
    if depth > max_depth {
      max_depth = depth
    }
    match go_graph.terms.get(cid) {
      Some(dterm) => {
        let mut pi : Int = 0
        while pi < dterm.parent_ids.length() {
          let pid : String = dterm.parent_ids[pi]
          if go_graph.terms.contains(pid) && !visited.contains(pid) {
            visited.set(pid, true)
            depth_queue.push((pid, depth + 1))
          }
          pi = pi + 1
        }
      }
      None => ()
    }
    idx = idx + 1
  }
  max_depth
}

///|
/// GOstats conditional test.
/// Uses the conditional method: for each GO term, remove genes mapped to
/// significant child terms before testing the parent term.
/// This reduces redundancy in the enrichment results.
pub fn gostats_conditional_test(
  hg_test : GOHyperGeoTest,
  ontology? : String = "",
  fdr_threshold? : Double = 0.05,
) -> Array[GOHyperGeoResult] {
  let cond_universe_set : Map[String, Bool] = gostats_string_set(
    hg_test.universe,
  )
  let cond_n_univ : Int = hg_test.universe.length()
  let cond_n : Int = hg_test.selected.length()

  let cond_term_keys : Array[String] = hg_test.go_graph.terms.keys().collect()

  let filtered_terms : Array[String] = Array(capacity=0)
  let mut fi : Int = 0
  while fi < cond_term_keys.length() {
    let fid : String = cond_term_keys[fi]
    let fterm : GOStatsTerm = hg_test.go_graph.terms.get(fid).unwrap()
    if ontology == "" || fterm.ontology == ontology {
      filtered_terms.push(fid)
    }
    fi = fi + 1
  }

  let term_level : Map[String, Int] = Map([], capacity=filtered_terms.length())
  let mut ti : Int = 0
  while ti < filtered_terms.length() {
    let tlid : String = filtered_terms[ti]
    let level : Int = gostats_term_depth(hg_test.go_graph, tlid)
    term_level.set(tlid, level)
    ti = ti + 1
  }

  let sorted_terms : Array[String] = Array(capacity=0)
  ti = 0
  while ti < filtered_terms.length() {
    sorted_terms.push(filtered_terms[ti])
    ti = ti + 1
  }

  let mut sa : Int = 0
  while sa < sorted_terms.length() - 1 {
    let mut sb : Int = sa + 1
    while sb < sorted_terms.length() {
      let la : Int = term_level.get(sorted_terms[sa]).unwrap_or(0)
      let lb : Int = term_level.get(sorted_terms[sb]).unwrap_or(0)
      if la > lb {
        let stmp : String = sorted_terms[sa]
        sorted_terms[sa] = sorted_terms[sb]
        sorted_terms[sb] = stmp
      }
      sb = sb + 1
    }
    sa = sa + 1
  }

  let removed_genes : Map[String, Bool] = Map([], capacity=cond_n)
  let cond_results : Array[GOHyperGeoResult] = Array(capacity=0)
  let cond_p_values : Array[Double] = Array(capacity=0)

  let mut si : Int = 0
  while si < sorted_terms.length() {
    let sgo_id : String = sorted_terms[si]
    let sterm : GOStatsTerm = hg_test.go_graph.terms.get(sgo_id).unwrap()

    let s_term_genes : Array[String] = gostats_intersection_genes(
      sterm.gene_ids,
      cond_universe_set,
    )
    let effective_genes : Array[String] = Array(capacity=0)
    let mut egi : Int = 0
    while egi < s_term_genes.length() {
      let egene : String = s_term_genes[egi]
      if !removed_genes.contains(egene) {
        effective_genes.push(egene)
      }
      egi = egi + 1
    }

    let em : Int = effective_genes.length()

    let mut eff_sel_count : Int = 0
    egi = 0
    while egi < effective_genes.length() {
      let mut ej : Int = 0
      while ej < hg_test.selected.length() {
        if effective_genes[egi] == hg_test.selected[ej] {
          eff_sel_count = eff_sel_count + 1
          break
        }
        ej = ej + 1
      }
      egi = egi + 1
    }

    let ek : Int = eff_sel_count

    if ek > 0 && em > 0 {
      let ep_val : Double = gostats_hypergeo_pvalue(cond_n_univ, em, cond_n, ek)
      cond_p_values.push(ep_val)

      let e_odds : Double = gostats_odds_ratio(cond_n_univ, em, cond_n, ek)
      let e_expected : Double = em.to_double() *
        cond_n.to_double() /
        cond_n_univ.to_double()

      cond_results.push(
        GOHyperGeoResult::new(
          sterm.go_id,
          sterm.name,
          sterm.ontology,
          ep_val,
          0.0,
          e_odds,
          e_expected,
          ek,
          em,
          cond_n_univ,
        ),
      )
    }

    si = si + 1
  }

  let cond_fdr_vals : Array[Double] = gostats_bh_fdr(cond_p_values)
  let mut ri : Int = 0
  while ri < cond_results.length() {
    cond_results[ri] = GOHyperGeoResult::{
      go_id: cond_results[ri].go_id,
      name: cond_results[ri].name,
      ontology: cond_results[ri].ontology,
      p_value: cond_results[ri].p_value,
      fdr: cond_fdr_vals[ri],
      odds_ratio: cond_results[ri].odds_ratio,
      expected: cond_results[ri].expected,
      observed: cond_results[ri].observed,
      gene_count: cond_results[ri].gene_count,
      universe_count: cond_results[ri].universe_count,
    }
    ri = ri + 1
  }

  ri = 0
  while ri < cond_results.length() {
    if cond_results[ri].fdr <= fdr_threshold {
      match hg_test.go_graph.terms.get(cond_results[ri].go_id) {
        Some(sig_term) => {
          let mut sgi : Int = 0
          while sgi < sig_term.gene_ids.length() {
            removed_genes.set(sig_term.gene_ids[sgi], true)
            sgi = sgi + 1
          }
        }
        None => ()
      }
    }
    ri = ri + 1
  }

  let mut ca : Int = 0
  while ca < cond_results.length() - 1 {
    let mut cb : Int = ca + 1
    while cb < cond_results.length() {
      if cond_results[ca].p_value > cond_results[cb].p_value {
        let ctmp : GOHyperGeoResult = cond_results[ca]
        cond_results[ca] = cond_results[cb]
        cond_results[cb] = ctmp
      }
      cb = cb + 1
    }
    ca = ca + 1
  }

  cond_results
}

///|
/// Filter GOHyperGeoResult array by FDR threshold and minimum observed count.
pub fn gostats_filter_results(
  results : Array[GOHyperGeoResult],
  fdr_cutoff : Double,
  min_observed : Int,
) -> Array[GOHyperGeoResult] {
  let filtered : Array[GOHyperGeoResult] = Array(capacity=0)
  let mut i4 : Int = 0
  while i4 < results.length() {
    let r : GOHyperGeoResult = results[i4]
    if r.fdr <= fdr_cutoff && r.observed >= min_observed {
      filtered.push(r)
    }
    i4 = i4 + 1
  }
  filtered
}

///|
/// Summarize GO enrichment results as a string.
pub fn gostats_summary(results : Array[GOHyperGeoResult]) -> String {
  let mut s : String = "GOstats Enrichment Results:\n"
  s = s + "Total terms tested: " + results.length().to_string() + "\n"

  let mut sig_count : Int = 0
  let mut i5 : Int = 0
  while i5 < results.length() {
    if results[i5].fdr <= 0.05 {
      sig_count = sig_count + 1
    }
    i5 = i5 + 1
  }
  s = s + "Significant terms (FDR <= 0.05): " + sig_count.to_string() + "\n"

  let mut bp_count : Int = 0
  let mut mf_count : Int = 0
  let mut cc_count : Int = 0
  i5 = 0
  while i5 < results.length() {
    let ont : String = results[i5].ontology
    if ont == "BP" {
      bp_count = bp_count + 1
    } else if ont == "MF" {
      mf_count = mf_count + 1
    } else if ont == "CC" {
      cc_count = cc_count + 1
    }
    i5 = i5 + 1
  }
  s = s + "BP terms: " + bp_count.to_string() + ", "
  s = s + "MF terms: " + mf_count.to_string() + ", "
  s = s + "CC terms: " + cc_count.to_string() + "\n"

  s
}
