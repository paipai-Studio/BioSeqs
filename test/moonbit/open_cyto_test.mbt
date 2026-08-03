///|
/// Tests for the openCyto module (Bioconductor openCyto).
///
/// Covers numeric helpers, KDE, gate construction, all built-in gating
/// methods (mindensity, tailgate, quantileGate, rangeGate, flowClust),
/// EM t-mixture fitting, t-pdf, lgamma, topological sort, channel
/// resolution, the full gating pipeline, population stats, and summaries.
///
/// Note: `oc_sort` is a private internal helper, so it is exercised
/// indirectly through the public `oc_median` / `oc_quantile` functions
/// that rely on it.

// ===========================================================================
// 1. Numeric helpers: oc_mean
// ===========================================================================

///|
test "oc_mean_basic" {
  let m = @src.oc_mean([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true((m - 3.0).abs() < 0.001)
}

///|
test "oc_mean_empty" {
  assert_eq(@src.oc_mean([]), 0.0)
}

///|
test "oc_mean_negative" {
  let m = @src.oc_mean([-1.0, -2.0, -3.0])
  assert_true((m - -2.0).abs() < 0.001)
}

// ===========================================================================
// 2. Numeric helpers: oc_sd
// ===========================================================================

///|
test "oc_sd_basic" {
  // Sample SD of [1,2,3,4,5]: sqrt(10/4) = sqrt(2.5) ~ 1.5811388
  let s = @src.oc_sd([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true((s - 1.5811388).abs() < 0.001)
}

///|
test "oc_sd_single_element" {
  assert_eq(@src.oc_sd([42.0]), 0.0)
}

///|
test "oc_sd_empty" {
  assert_eq(@src.oc_sd([]), 0.0)
}

// ===========================================================================
// 3. Numeric helpers: oc_median
// ===========================================================================

///|
test "oc_median_odd" {
  let m = @src.oc_median([3.0, 1.0, 2.0])
  assert_true((m - 2.0).abs() < 0.001)
}

///|
test "oc_median_even" {
  let m = @src.oc_median([1.0, 2.0, 3.0, 4.0])
  assert_true((m - 2.5).abs() < 0.001)
}

///|
test "oc_median_empty" {
  assert_eq(@src.oc_median([]), 0.0)
}

// ===========================================================================
// 4. Numeric helpers: oc_quantile
// ===========================================================================

///|
test "oc_quantile_median" {
  let q = @src.oc_quantile([1.0, 2.0, 3.0, 4.0], 0.5)
  assert_true((q - 2.5).abs() < 0.001)
}

///|
test "oc_quantile_min" {
  let q = @src.oc_quantile([1.0, 2.0, 3.0, 4.0], 0.0)
  assert_true((q - 1.0).abs() < 0.001)
}

///|
test "oc_quantile_max" {
  let q = @src.oc_quantile([1.0, 2.0, 3.0, 4.0], 1.0)
  assert_true((q - 4.0).abs() < 0.001)
}

///|
test "oc_quantile_quarter" {
  // pos = 0.25 * 3 = 0.75 -> lo=0, hi=1, frac=0.75 -> 1.0*0.25 + 2.0*0.75 = 1.75
  let q = @src.oc_quantile([1.0, 2.0, 3.0, 4.0], 0.25)
  assert_true((q - 1.75).abs() < 0.001)
}

// ===========================================================================
// 5. oc_sort (private) exercised via oc_quantile ordering
// ===========================================================================

///|
test "oc_sort_via_quantile_ordering" {
  // Shuffled input: oc_quantile relies on oc_sort internally.
  // q=0, 0.5, 1.0 must come out as min, median, max of the data.
  let data = [5.0, 1.0, 3.0, 2.0, 4.0]
  let lo = @src.oc_quantile(data, 0.0)
  let mid = @src.oc_quantile(data, 0.5)
  let hi = @src.oc_quantile(data, 1.0)
  assert_true((lo - 1.0).abs() < 0.001)
  assert_true((mid - 3.0).abs() < 0.001)
  assert_true((hi - 5.0).abs() < 0.001)
  assert_true(lo <= mid)
  assert_true(mid <= hi)
}

// ===========================================================================
// 6. KDE: oc_kde_evaluate
// ===========================================================================

///|
test "oc_kde_evaluate_positive" {
  let d = @src.oc_kde_evaluate([1.0, 2.0, 3.0], 2.0, 1.0)
  assert_true(d > 0.0)
}

///|
test "oc_kde_evaluate_empty" {
  assert_eq(@src.oc_kde_evaluate([], 2.0, 1.0), 0.0)
}

///|
test "oc_kde_evaluate_bad_bandwidth" {
  // Non-positive bandwidth should yield 0.
  assert_eq(@src.oc_kde_evaluate([1.0, 2.0], 1.5, 0.0), 0.0)
  assert_eq(@src.oc_kde_evaluate([1.0, 2.0], 1.5, -1.0), 0.0)
}

///|
test "oc_kde_evaluate_peak_near_data" {
  // Density at the data cluster center should exceed density far away.
  let data = [1.0, 1.1, 0.9, 1.2]
  let at_peak = @src.oc_kde_evaluate(data, 1.0, 0.5)
  let far = @src.oc_kde_evaluate(data, 10.0, 0.5)
  assert_true(at_peak > far)
}

// ===========================================================================
// 7. KDE: oc_kde_grid
// ===========================================================================

///|
test "oc_kde_grid_lengths" {
  let (gx, gy) = @src.oc_kde_grid([1.0, 2.0, 3.0, 4.0, 5.0], 20, 1.0)
  assert_eq(gx.length(), 20)
  assert_eq(gy.length(), 20)
  // Grid should span beyond the data range (10% padding each side).
  assert_true(gx[0] < 1.0)
  assert_true(gx[19] > 5.0)
}

///|
test "oc_kde_grid_empty" {
  let (gx, gy) = @src.oc_kde_grid([], 10, 1.0)
  assert_eq(gx.length(), 0)
  assert_eq(gy.length(), 0)
}

///|
test "oc_kde_grid_degenerate" {
  // All-equal data: range is 0, returns a single point.
  let (gx, gy) = @src.oc_kde_grid([3.0, 3.0, 3.0], 10, 1.0)
  assert_eq(gx.length(), 1)
  assert_true((gx[0] - 3.0).abs() < 0.001)
  assert_true(gy[0] > 0.0)
}

// ===========================================================================
// 8. KDE: oc_bandwidth_silverman
// ===========================================================================

///|
test "oc_bandwidth_silverman_positive" {
  let h = @src.oc_bandwidth_silverman([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true(h > 0.0)
}

///|
test "oc_bandwidth_silverman_single" {
  assert_eq(@src.oc_bandwidth_silverman([7.0]), 1.0)
}

///|
test "oc_bandwidth_silverman_constant" {
  // Constant data -> sd=0 and iqr=0 -> returns fallback 1.0.
  assert_eq(@src.oc_bandwidth_silverman([5.0, 5.0, 5.0, 5.0]), 1.0)
}

// ===========================================================================
// 9. Gate construction: OcGate::new1d / new2d
// ===========================================================================

///|
test "oc_gate_new1d_fields" {
  let g = @src.OcGate::new1d(dim="FSC", min=2.0, max=5.0)
  assert_eq(g.dim1, "FSC")
  assert_true((g.min1 - 2.0).abs() < 0.001)
  assert_true((g.max1 - 5.0).abs() < 0.001)
  let is_1d = match g.dim2 {
    None => true
    Some(_) => false
  }
  assert_true(is_1d)
}

///|
test "oc_gate_new2d_fields" {
  let g = @src.OcGate::new2d(
    dim1="FSC", dim2="SSC", min1=1.0, max1=2.0, min2=3.0, max2=4.0,
  )
  assert_eq(g.dim1, "FSC")
  assert_eq(g.dim2.unwrap(), "SSC")
  assert_true((g.min1 - 1.0).abs() < 0.001)
  assert_true((g.max1 - 2.0).abs() < 0.001)
  assert_true((g.min2 - 3.0).abs() < 0.001)
  assert_true((g.max2 - 4.0).abs() < 0.001)
}

// ===========================================================================
// 10. oc_event_in_gate
// ===========================================================================

///|
test "oc_event_in_gate_1d" {
  let g = @src.OcGate::new1d(dim="FSC", min=2.0, max=5.0)
  assert_true(@src.oc_event_in_gate(g, 3.0, 0.0))
  assert_false(@src.oc_event_in_gate(g, 1.0, 0.0))
  assert_false(@src.oc_event_in_gate(g, 6.0, 0.0))
}

///|
test "oc_event_in_gate_1d_boundary" {
  let g = @src.OcGate::new1d(dim="FSC", min=2.0, max=5.0)
  // Boundaries are inclusive.
  assert_true(@src.oc_event_in_gate(g, 2.0, 0.0))
  assert_true(@src.oc_event_in_gate(g, 5.0, 0.0))
}

///|
test "oc_event_in_gate_2d" {
  let g = @src.OcGate::new2d(
    dim1="FSC", dim2="SSC", min1=1.0, max1=3.0, min2=10.0, max2=20.0,
  )
  // Both dims in range.
  assert_true(@src.oc_event_in_gate(g, 2.0, 15.0))
  // First dim out of range.
  assert_false(@src.oc_event_in_gate(g, 5.0, 15.0))
  // Second dim out of range.
  assert_false(@src.oc_event_in_gate(g, 2.0, 99.0))
}

// ===========================================================================
// 11. mindensity
// ===========================================================================

///|
test "oc_mindensity_bimodal" {
  // Two clear clusters around 1.0 and 5.0; valley near 3.0.
  let data = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let args : Map[String, Double] = Map::new()
  args.set("bandwidth", 0.5)
  let g = @src.oc_mindensity(data, "FSC", args)
  assert_eq(g.dim1, "FSC")
  // Threshold should fall in the valley between the two clusters.
  assert_true(g.min1 > 2.0)
  assert_true(g.min1 < 4.0)
}

///|
test "oc_mindensity_empty" {
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_mindensity([], "FSC", args)
  assert_eq(g.dim1, "FSC")
  assert_true((g.min1 - 0.0).abs() < 0.001)
}

///|
test "oc_mindensity_single_peak" {
  // Unimodal data -> falls back to the median.
  let data = [1.0, 1.1, 1.0, 1.2, 0.9, 1.1]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_mindensity(data, "FSC", args)
  let med = @src.oc_median(data)
  assert_true((g.min1 - med).abs() < 0.001)
}

// ===========================================================================
// 12. tailgate
// ===========================================================================

///|
test "oc_tailgate_default_quantile" {
  // Default quantile = 0.99.
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_tailgate(data, "CD4", args)
  assert_eq(g.dim1, "CD4")
  // pos = 0.99 * 9 = 8.91 -> 9.0*0.09 + 10.0*0.91 = 9.91
  assert_true((g.min1 - 9.91).abs() < 0.01)
}

///|
test "oc_tailgate_custom_quantile" {
  let data = [10.0, 20.0, 30.0, 40.0, 50.0]
  let args : Map[String, Double] = Map::new()
  args.set("quantile", 0.5)
  let g = @src.oc_tailgate(data, "CD4", args)
  // pos = 0.5 * 4 = 2.0 -> threshold = 30.0
  assert_true((g.min1 - 30.0).abs() < 0.001)
}

///|
test "oc_tailgate_empty" {
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_tailgate([], "CD4", args)
  assert_true((g.min1 - 0.0).abs() < 0.001)
}

// ===========================================================================
// 13. quantileGate
// ===========================================================================

///|
test "oc_quantile_gate_default" {
  // Default probs = 0.95.
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_quantile_gate(data, "CD8", args)
  assert_eq(g.dim1, "CD8")
  // pos = 0.95 * 9 = 8.55 -> 9.0*0.45 + 10.0*0.55 = 9.55
  assert_true((g.min1 - 9.55).abs() < 0.01)
}

///|
test "oc_quantile_gate_custom" {
  let data = [10.0, 20.0, 30.0, 40.0, 50.0]
  let args : Map[String, Double] = Map::new()
  args.set("probs", 0.5)
  let g = @src.oc_quantile_gate(data, "CD8", args)
  assert_true((g.min1 - 30.0).abs() < 0.001)
}

///|
test "oc_quantile_gate_empty" {
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_quantile_gate([], "CD8", args)
  assert_true((g.min1 - 0.0).abs() < 0.001)
}

// ===========================================================================
// 14. rangeGate
// ===========================================================================

///|
test "oc_range_gate_bimodal" {
  let data = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_range_gate(data, "FSC", args)
  assert_eq(g.dim1, "FSC")
  // Threshold must lie within the data range.
  assert_true(g.min1 >= 0.9)
  assert_true(g.min1 <= 5.2)
}

///|
test "oc_range_gate_empty" {
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_range_gate([], "FSC", args)
  assert_true((g.min1 - 0.0).abs() < 0.001)
}

// ===========================================================================
// 15. flowClust (1D and 2D)
// ===========================================================================

///|
test "oc_flow_clust_1d_bimodal" {
  let data = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_flow_clust(data, [], "FSC", None, args, 2)
  assert_eq(g.dim1, "FSC")
  // Boundary is the midpoint of the two fitted means -> near 3.0.
  assert_true(g.min1 > 2.0)
  assert_true(g.min1 < 4.0)
}

///|
test "oc_flow_clust_2d_returns_2d_gate" {
  let v1 = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let v2 = [10.0, 10.1, 9.9, 10.2, 50.0, 50.1, 49.9, 50.2]
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_flow_clust(v1, v2, "FSC", Some("SSC"), args, 2)
  assert_eq(g.dim1, "FSC")
  assert_eq(g.dim2.unwrap(), "SSC")
  // Gate bounds should be ordered and finite.
  assert_true(g.max1 > g.min1)
  assert_true(g.max2 > g.min2)
}

///|
test "oc_flow_clust_1d_empty" {
  let args : Map[String, Double] = Map::new()
  let g = @src.oc_flow_clust([], [], "FSC", None, args, 2)
  assert_true((g.min1 - 0.0).abs() < 0.001)
}

// ===========================================================================
// 16. EM t-mixture (1D and 2D)
// ===========================================================================

///|
test "oc_em_t_mixture_1d_bimodal" {
  let data = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let (means, weights) = @src.oc_em_t_mixture_1d(data, 2, 50)
  assert_eq(means.length(), 2)
  assert_eq(weights.length(), 2)
  let mut lo = means[0]
  let mut hi = means[0]
  for m in means {
    if m < lo {
      lo = m
    }
    if m > hi {
      hi = m
    }
  }
  // One mean near the low cluster, one near the high cluster.
  assert_true(lo < 2.5)
  assert_true(hi > 2.5)
  assert_true(lo >= 0.5)
  assert_true(hi <= 5.5)
}

///|
test "oc_em_t_mixture_1d_single_cluster" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let (means, weights) = @src.oc_em_t_mixture_1d(data, 1, 10)
  assert_eq(means.length(), 1)
  assert_eq(weights.length(), 1)
  assert_true((means[0] - 3.0).abs() < 0.001)
  assert_true((weights[0] - 1.0).abs() < 0.001)
}

///|
test "oc_em_t_mixture_1d_empty" {
  let (means, weights) = @src.oc_em_t_mixture_1d([], 2, 10)
  assert_eq(means.length(), 0)
  assert_eq(weights.length(), 0)
}

///|
test "oc_em_t_mixture_2d_bimodal" {
  let v1 = [1.0, 1.1, 0.9, 1.2, 5.0, 5.1, 4.9, 5.2]
  let v2 = [10.0, 10.1, 9.9, 10.2, 50.0, 50.1, 49.9, 50.2]
  let (m1, m2) = @src.oc_em_t_mixture_2d(v1, v2, 2, 50)
  assert_eq(m1.length(), 2)
  assert_eq(m2.length(), 2)
  // dim1 means should straddle the valley.
  let mut lo1 = m1[0]
  let mut hi1 = m1[0]
  for m in m1 {
    if m < lo1 {
      lo1 = m
    }
    if m > hi1 {
      hi1 = m
    }
  }
  assert_true(lo1 < 2.5)
  assert_true(hi1 > 2.5)
}

// ===========================================================================
// 17. t-pdf
// ===========================================================================

///|
test "oc_t_pdf_at_mean_is_positive" {
  let p = @src.oc_t_pdf(0.0, 0.0, 1.0, 4.0)
  assert_true(p > 0.0)
}

///|
test "oc_t_pdf_symmetric" {
  // t-pdf is symmetric about its mean.
  let plus = @src.oc_t_pdf(3.0, 1.0, 1.0, 4.0)
  let minus = @src.oc_t_pdf(-1.0, 1.0, 1.0, 4.0)
  assert_true((plus - minus).abs() < 0.0001)
}

///|
test "oc_t_pdf_mode_at_mean" {
  // Density at the mean (t=0) should exceed density away from the mean.
  let at_mean = @src.oc_t_pdf(5.0, 5.0, 1.0, 4.0)
  let away = @src.oc_t_pdf(7.0, 5.0, 1.0, 4.0)
  assert_true(at_mean > away)
}

///|
test "oc_t_pdf_zero_variance" {
  assert_eq(@src.oc_t_pdf(1.0, 1.0, 0.0, 4.0), 0.0)
}

// ===========================================================================
// 18. lgamma
// ===========================================================================

///|
test "oc_lgamma_one" {
  // lgamma(1) = log(gamma(1)) = log(1) = 0
  assert_true((@src.oc_lgamma(1.0)).abs() < 0.001)
}

///|
test "oc_lgamma_two" {
  // lgamma(2) = log(gamma(2)) = log(1) = 0
  assert_true((@src.oc_lgamma(2.0)).abs() < 0.001)
}

///|
test "oc_lgamma_half" {
  // lgamma(0.5) = log(sqrt(pi)) = 0.5 * log(pi) ~ 0.572365
  let v = @src.oc_lgamma(0.5)
  assert_true((v - 0.572365).abs() < 0.01)
}

///|
test "oc_lgamma_five" {
  // lgamma(5) = log(4!) = log(24) ~ 3.178054
  let v = @src.oc_lgamma(5.0)
  assert_true((v - 3.178054).abs() < 0.01)
}

// ===========================================================================
// 19. Topological sort
// ===========================================================================

///|
test "oc_topo_sort_parent_before_child" {
  // Define rules out of dependency order; topo sort must put the parent
  // population before its children.
  let rule_cd4 = @src.OcGatingRule::new(
    parent="lymph",
    child="cd4",
    method="tailgate",
    dims=["CD4"],
    args=Map::new(),
  )
  let rule_cd8 = @src.OcGatingRule::new(
    parent="lymph",
    child="cd8",
    method="tailgate",
    dims=["CD8"],
    args=Map::new(),
  )
  let rule_lymph = @src.OcGatingRule::new(
    parent="root",
    child="lymph",
    method="mindensity",
    dims=["FSC"],
    args=Map::new(),
  )
  let sorted = @src.oc_topo_sort([rule_cd4, rule_cd8, rule_lymph])
  assert_eq(sorted.length(), 3)
  // "lymph" must come before "cd4" and "cd8".
  assert_eq(sorted[0].child, "lymph")
  let mut lymph_idx = 0
  let mut cd4_idx = 0
  for i in 0..<sorted.length() {
    if sorted[i].child == "lymph" {
      lymph_idx = i
    }
    if sorted[i].child == "cd4" {
      cd4_idx = i
    }
  }
  assert_true(lymph_idx < cd4_idx)
}

// ===========================================================================
// 20. Channel resolution
// ===========================================================================

///|
test "oc_resolve_channel_found" {
  let idx = @src.oc_resolve_channel(["FSC", "SSC", "FITC"], "SSC")
  assert_eq(idx, 1)
}

///|
test "oc_resolve_channel_not_found" {
  let idx = @src.oc_resolve_channel(["FSC", "SSC", "FITC"], "PE")
  assert_eq(idx, -1)
}

// ===========================================================================
// 21. Channel extraction helper
// ===========================================================================

///|
test "oc_extract_channel" {
  let events = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let col = @src.oc_extract_channel(events, 1)
  assert_eq(col.length(), 3)
  assert_true((col[0] - 2.0).abs() < 0.001)
  assert_true((col[1] - 4.0).abs() < 0.001)
  assert_true((col[2] - 6.0).abs() < 0.001)
}

// ===========================================================================
// 22. OcFlowSet and OcGatingRule construction
// ===========================================================================

///|
test "oc_flowset_construction" {
  let fs = @src.OcFlowSet::new(
    sample_names=["S1", "S2"],
    channel_names=["FSC", "SSC"],
    data=[[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]],
  )
  assert_eq(fs.sample_names.length(), 2)
  assert_eq(fs.channel_names.length(), 2)
  assert_eq(fs.data.length(), 2)
  assert_eq(fs.data[0].length(), 2)
  assert_true((fs.data[0][0][0] - 1.0).abs() < 0.001)
}

///|
test "oc_gating_rule_construction" {
  let args : Map[String, Double] = Map::new()
  args.set("probs", 0.95)
  let r = @src.OcGatingRule::new(
    parent="root",
    child="tcells",
    method="quantileGate",
    dims=["CD3"],
    args=args,
  )
  assert_eq(r.parent, "root")
  assert_eq(r.child, "tcells")
  assert_eq(r.method, "quantileGate")
  assert_eq(r.dims.length(), 1)
  assert_eq(r.dims[0], "CD3")
  assert_true((r.args.get("probs").unwrap() - 0.95).abs() < 0.001)
}

// ===========================================================================
// 23. Full pipeline: oc_gate_flow_set (single rule)
// ===========================================================================

///|
test "oc_gate_flow_set_single_rule" {
  // One sample, one channel (FSC), bimodal FSC values.
  let fs = @src.OcFlowSet::new(
    sample_names=["S1"],
    channel_names=["FSC"],
    data=[
      [
        [1.0],
        [1.1],
        [0.9],
        [1.2],
        [5.0],
        [5.1],
        [4.9],
        [5.2],
      ],
    ],
  )
  let args : Map[String, Double] = Map::new()
  args.set("bandwidth", 0.5)
  let rule = @src.OcGatingRule::new(
    parent="root",
    child="cells",
    method="mindensity",
    dims=["FSC"],
    args=args,
  )
  let results = @src.oc_gate_flow_set(fs, [rule])
  assert_eq(results.length(), 1)
  assert_eq(results[0].population, "cells")
  assert_eq(results[0].sample, "S1")
  assert_eq(results[0].parent_events, 8)
  // High-FSC cluster has 4 events (threshold ~3.0).
  assert_eq(results[0].child_events, 4)
}

// ===========================================================================
// 24. Full pipeline: oc_gate_flow_set (parent->child chain)
// ===========================================================================

///|
test "oc_gate_flow_set_chain" {
  // Two channels: FSC (bimodal) and CD4.
  let fs = @src.OcFlowSet::new(
    sample_names=["S1"],
    channel_names=["FSC", "CD4"],
    data=[
      [
        [1.0, 100.0],
        [1.1, 200.0],
        [0.9, 300.0],
        [1.2, 400.0],
        [5.0, 500.0],
        [5.1, 600.0],
        [4.9, 700.0],
        [5.2, 800.0],
      ],
    ],
  )
  let args1 : Map[String, Double] = Map::new()
  args1.set("bandwidth", 0.5)
  let rule_cells = @src.OcGatingRule::new(
    parent="root",
    child="cells",
    method="mindensity",
    dims=["FSC"],
    args=args1,
  )
  let args2 : Map[String, Double] = Map::new()
  args2.set("quantile", 0.5)
  let rule_cd4 = @src.OcGatingRule::new(
    parent="cells",
    child="cd4pos",
    method="tailgate",
    dims=["CD4"],
    args=args2,
  )
  let results = @src.oc_gate_flow_set(fs, [rule_cd4, rule_cells])
  assert_eq(results.length(), 2)
  // Topological order: "cells" before "cd4pos".
  assert_eq(results[0].population, "cells")
  assert_eq(results[1].population, "cd4pos")
  // The parent population of cd4pos is cells, so its parent_events
  // must equal the child_events of cells.
  assert_eq(results[1].parent_events, results[0].child_events)
}

// ===========================================================================
// 25. Population stats
// ===========================================================================

///|
test "oc_population_stats_basic" {
  // Build a gating result by hand: 8 parent events, 4 child events.
  let gate = @src.OcGate::new1d(dim="FSC", min=3.0, max=1.0e30)
  let indices = [
    false, false, false, false, true, true, true, true,
  ]
  let result = @src.OcGatingResult::new(
    population="cells",
    sample="S1",
    gate=gate,
    parent_events=8,
    child_events=4,
    indices=indices,
  )
  let stats = @src.oc_population_stats([result], 8)
  assert_eq(stats.length(), 1)
  assert_eq(stats[0].population, "cells")
  assert_eq(stats[0].sample, "S1")
  assert_eq(stats[0].count, 4)
  assert_eq(stats[0].parent_count, 8)
  assert_true((stats[0].percent_of_parent - 50.0).abs() < 0.001)
  assert_true((stats[0].percent_of_total - 50.0).abs() < 0.001)
}

///|
test "oc_population_stats_zero_parent" {
  // When parent_events == 0, percentages should be 0 (no division by zero).
  let gate = @src.OcGate::new1d(dim="FSC", min=0.0, max=0.0)
  let result = @src.OcGatingResult::new(
    population="dead",
    sample="S1",
    gate=gate,
    parent_events=0,
    child_events=0,
    indices=[],
  )
  let stats = @src.oc_population_stats([result], 10)
  assert_eq(stats[0].percent_of_parent, 0.0)
  assert_eq(stats[0].percent_of_total, 0.0)
}

// ===========================================================================
// 26. Gating summary
// ===========================================================================

///|
test "oc_gating_summary_string" {
  let gate = @src.OcGate::new1d(dim="FSC", min=3.0, max=1.0e30)
  let result = @src.OcGatingResult::new(
    population="cells",
    sample="S1",
    gate=gate,
    parent_events=8,
    child_events=4,
    indices=[false, false, false, false, true, true, true, true],
  )
  let s = @src.oc_gating_summary([result])
  // The summary header should be present.
  let has_header = s.contains("openCyto Gating Summary")
  assert_true(has_header)
  // The population and sample names should appear.
  assert_true(s.contains("cells"))
  assert_true(s.contains("S1"))
}

// ===========================================================================
// 27. oc_ln helper
// ===========================================================================

///|
test "oc_ln_one" {
  assert_true((@src.oc_ln(1.0)).abs() < 0.001)
}

///|
test "oc_ln_two" {
  // ln(2) ~ 0.693147
  assert_true((@src.oc_ln(2.0) - 0.693147).abs() < 0.001)
}

///|
test "oc_ln_non_positive" {
  // Non-positive input returns the sentinel -1e30.
  assert_true(@src.oc_ln(0.0) < -1.0e29)
  assert_true(@src.oc_ln(-1.0) < -1.0e29)
}
