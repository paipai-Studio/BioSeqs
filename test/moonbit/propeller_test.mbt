///|
/// Tests for the speckle propeller port (getTransformedProps, propeller
/// t-tests/ANOVA with limma empirical Bayes including robust variance
/// shrinkage, estimateBetaParam(FromCounts), normCounts). All expected values
/// were frozen from speckle 1.12.0 with limma 3.68.5 under R 4.3.3 using
/// /tmp/gen_speckle_expected.R and /tmp/gen_speckle_robust.R (R sources
/// sourced directly; %.17g full precision).

///|
fn spkt_approx(
  left : Double,
  right : Double,
  tolerance? : Double = 1.0e-10,
) -> Bool {
  let scale = 1.0 + right.abs()
  (left - right).abs() <= tolerance * scale
}

///|
fn spkt_vec(
  actual : Array[Double],
  want : Array[Double],
  tolerance? : Double = 1.0e-10,
) -> Unit {
  if actual.length() != want.length() {
    abort("length \{actual.length()} != \{want.length()}")
  }
  for i in 0..<want.length() {
    if !spkt_approx(actual[i], want[i], tolerance~) {
      abort("index \{i}: got \{actual[i]} want \{want[i]}")
    }
  }
}

///|
fn spkt_mat(
  actual : @src.Matrix,
  want : Array[Array[Double]],
  tolerance? : Double = 1.0e-10,
) -> Unit {
  if actual.rows != want.length() || actual.cols != want[0].length() {
    abort(
      "matrix dim \{actual.rows}x\{actual.cols} != \{want.length()}x\{want[0].length()}",
    )
  }
  for i in 0..<want.length() {
    for j in 0..<want[i].length() {
      if !spkt_approx(actual.get(i, j), want[i][j], tolerance~) {
        abort("cell (\{i},\{j}): got \{actual.get(i, j)} want \{want[i][j]}")
      }
    }
  }
}

///|
fn spkt_raises(operation : () -> Unit raise @src.PropellerError) -> Bool {
  try {
    operation()
    false
  } catch {
    _ => true
  }
}

///|
/// Expand a counts matrix (types x samples) into cell-level cluster/sample
/// vectors, as required by get_transformed_props and propeller.
fn spkt_expand(
  tab : @src.Matrix,
  cluster_prefix~ : String,
  sample_prefix~ : String,
) -> (Array[String], Array[String]) {
  let clusters = []
  let samples = []
  for i in 0..<tab.rows {
    for j in 0..<tab.cols {
      let count = tab.get(i, j).round().to_int()
      for _ in 0..<count {
        clusters.push("\{cluster_prefix}\{i}")
        samples.push("\{sample_prefix}\{j + 1}")
      }
    }
  }
  (clusters, samples)
}

///|
fn spkt_counts_a() -> @src.Matrix {
  @src.Matrix::from_2d(3, 4, [
    [500.0, 900.0, 270.0, 480.0],
    [300.0, 450.0, 360.0, 360.0],
    [200.0, 150.0, 270.0, 360.0],
  ])
}

///|
fn spkt_design_a() -> @src.Matrix {
  @src.Matrix::from_2d(4, 2, [[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
}

///|
fn spkt_design_asex() -> @src.Matrix {
  @src.Matrix::from_2d(4, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 1.0],
  ])
}

///|
fn spkt_counts_b() -> @src.Matrix {
  @src.Matrix::from_2d(3, 6, [
    [500.0, 900.0, 270.0, 480.0, 800.0, 600.0],
    [300.0, 450.0, 360.0, 360.0, 100.0, 160.0],
    [200.0, 150.0, 270.0, 360.0, 100.0, 40.0],
  ])
}

///|
fn spkt_design_b() -> @src.Matrix {
  @src.Matrix::from_2d(6, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
  ])
}

///|
fn spkt_counts_c() -> @src.Matrix {
  @src.Matrix::from_2d(6, 8, [
    [420.0, 510.0, 380.0, 460.0, 80.0, 120.0, 60.0, 100.0],
    [300.0, 260.0, 330.0, 290.0, 520.0, 480.0, 550.0, 500.0],
    [150.0, 180.0, 140.0, 170.0, 210.0, 230.0, 190.0, 220.0],
    [90.0, 70.0, 110.0, 80.0, 140.0, 120.0, 160.0, 130.0],
    [40.0, 55.0, 35.0, 45.0, 90.0, 75.0, 105.0, 85.0],
    [0.0, 5.0, 0.0, 3.0, 60.0, 45.0, 70.0, 55.0],
  ])
}

///|
fn spkt_design_c() -> @src.Matrix {
  @src.Matrix::from_2d(8, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [0.0, 1.0],
    [0.0, 1.0],
    [0.0, 1.0],
  ])
}

///|
fn spkt_counts_d() -> @src.Matrix {
  @src.Matrix::from_2d(5, 9, [
    [100.0, 120.0, 110.0, 300.0, 280.0, 320.0, 50.0, 60.0, 55.0],
    [200.0, 180.0, 210.0, 90.0, 110.0, 80.0, 150.0, 140.0, 160.0],
    [60.0, 70.0, 55.0, 40.0, 35.0, 45.0, 200.0, 190.0, 210.0],
    [30.0, 25.0, 35.0, 150.0, 160.0, 140.0, 20.0, 18.0, 22.0],
    [0.0, 5.0, 0.0, 20.0, 10.0, 30.0, 80.0, 90.0, 70.0],
  ])
}

///|
fn spkt_design_d() -> @src.Matrix {
  @src.Matrix::from_2d(9, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
  ])
}

// ---------------------------------------------------------------------------
// getTransformedProps
// ---------------------------------------------------------------------------

///|
test "spk_get_transformed_props_asin_counts_and_values" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  spkt_mat(
    pl.counts,
    [
      [500.0, 900.0, 270.0, 480.0],
      [300.0, 450.0, 360.0, 360.0],
      [200.0, 150.0, 270.0, 360.0],
    ],
    tolerance=1.0e-12,
  )
  spkt_mat(
    pl.transformed,
    [
      [
        0.78539816339744839, 0.88607712379261372, 0.57963974036370425, 0.68471920300228295,
      ],
      [
        0.57963974036370425, 0.57963974036370425, 0.68471920300228295, 0.57963974036370425,
      ],
      [
        0.46364760900080609, 0.32175055439664219, 0.57963974036370425, 0.57963974036370425,
      ],
    ],
    tolerance=1.0e-12,
  )
  assert_eq(pl.cluster_names, ["c0", "c1", "c2"])
  assert_eq(pl.sample_names, ["s1", "s2", "s3", "s4"])
}

///|
test "spk_get_transformed_props_logit_values" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples)
  spkt_mat(
    pl.proportions,
    [
      [0.5, 0.59999999999999998, 0.29999999999999999, 0.40000000000000002],
      [
        0.29999999999999999, 0.29999999999999999, 0.40000000000000002, 0.29999999999999999,
      ],
      [
        0.20000000000000001, 0.10000000000000001, 0.29999999999999999, 0.29999999999999999,
      ],
    ],
    tolerance=1.0e-12,
  )
  spkt_mat(
    pl.transformed,
    [
      [
        -0.0009985023295895493, 0.40435523008080321, -0.84703376425467702, -0.40581190887460183,
      ],
      [
        -0.8470601330593277, -0.84713929382855868, -0.40592732212147792, -0.84709970324794426,
      ],
      [
        -1.3850467003217353, -2.1946372537712682, -0.84703376425467702, -0.84709970324794426,
      ],
    ],
    tolerance=1.0e-12,
  )
}

// ---------------------------------------------------------------------------
// Fixture A: two-group moderated t-test
// ---------------------------------------------------------------------------

///|
test "spk_ttest_fixture_a_asin_matches_limma" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let design = spkt_design_a()
  let model = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    design,
    [1.0, -1.0],
  )
  for v in model.df_prior {
    assert_true(v > 1.0e200)
  }
  assert_true(spkt_approx(model.s2_prior, 0.0043628678423677722))
  spkt_vec(model.df_total, [6.0, 6.0, 6.0], tolerance=1.0e-12)
  spkt_vec(
    model.statistic,
    [3.0817862495555488, -0.7954297290756619, -2.8302039949815678],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.p_value,
    [0.021613588088540914, 0.4566687197787217, 0.029949742218953398],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.fdr,
    [0.044924613328430098, 0.4566687197787217, 0.044924613328430098],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.prop_ratio,
    [1.5714285714285712, 0.8571428571428571, 0.5],
    tolerance=1.0e-12,
  )
  spkt_mat(
    model.means,
    [
      [0.54999999999999993, 0.34999999999999998],
      [0.29999999999999999, 0.34999999999999998],
      [0.14999999999999997, 0.29999999999999993],
    ],
    tolerance=1.0e-12,
  )
}

///|
test "spk_ttest_fixture_a_robust_and_nonrobust_agree" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let design = spkt_design_a()
  let mr = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    design,
    [1.0, -1.0],
    robust=true,
  )
  let mn = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    design,
    [1.0, -1.0],
    robust=false,
  )
  spkt_vec(mr.statistic, mn.statistic, tolerance=1.0e-12)
  spkt_vec(mr.p_value, mn.p_value, tolerance=1.0e-12)
}

///|
test "spk_ttest_fixture_a_logit_matches_limma" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="logit")
  let design = spkt_design_a()
  let model = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    design,
    [1.0, -1.0],
    robust=false,
  )
  assert_true(spkt_approx(model.s2_prior, 0.10075485989683063))
  spkt_vec(
    model.statistic,
    [2.6088578279803816, -0.69493684623235819, -2.9701280135678827],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.p_value,
    [0.040181918397897719, 0.51309355987947658, 0.024953385783993101],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.fdr,
    [0.060272877596846575, 0.51309355987947658, 0.060272877596846575],
    tolerance=1.0e-10,
  )
}

///|
test "spk_ttest_rows_sorted_and_unsorted" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let sorted = @src.propeller_ttest(
    pl,
    spkt_design_a(),
    [1.0, -1.0],
    robust=false,
  )
  assert_eq(sorted.length(), 3)
  assert_eq(sorted[0].cluster, "c0")
  assert_eq(sorted[1].cluster, "c2")
  assert_eq(sorted[2].cluster, "c1")
  assert_true(spkt_approx(sorted[0].statistic, 3.0817862495555488))
  assert_true(sorted[0].p_value < sorted[1].p_value)
  assert_true(sorted[1].p_value < sorted[2].p_value)
  let raw = @src.propeller_ttest(
    pl,
    spkt_design_a(),
    [1.0, -1.0],
    robust=false,
    sort=false,
  )
  assert_eq(raw[0].cluster, "c0")
  assert_eq(raw[1].cluster, "c1")
  assert_eq(raw[2].cluster, "c2")
  spkt_vec(
    raw[1].means,
    [0.29999999999999999, 0.34999999999999998],
    tolerance=1.0e-12,
  )
}

///|
test "spk_ttest_fixture_asex_confounded_covariate" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let model = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    spkt_design_asex(),
    [1.0, -1.0, 0.0],
    robust=false,
  )
  assert_true(
    spkt_approx(model.df_prior[0], 0.67461890146092429, tolerance=1.0e-10),
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.00019497241148801591, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [1.6746189014609243, 1.6746189014609243, 1.6746189014609243],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.statistic,
    [22.557037110574001, -1.2642997246927716, -3.3660232816234048],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.p_value,
    [0.0044510296395541115, 0.35387577754376498, 0.099109422516639248],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.fdr,
    [0.013353088918662335, 0.35387577754376498, 0.14866413377495888],
    tolerance=1.0e-9,
  )
  // Zero-contrast sex column is dropped when estimating proportion means.
  spkt_mat(
    model.means,
    [
      [0.54999999999999993, 0.34999999999999998],
      [0.29999999999999999, 0.34999999999999998],
      [0.14999999999999997, 0.29999999999999993],
    ],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.prop_ratio,
    [1.5714285714285712, 0.8571428571428571, 0.5],
    tolerance=1.0e-12,
  )
}

// ---------------------------------------------------------------------------
// High-level wrapper
// ---------------------------------------------------------------------------

///|
test "spk_wrapper_fixture_a_asin_baseline_and_order" {
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let groups = []
  for cell in 0..<clusters.length() {
    groups.push(
      if samples[cell] == "s1" || samples[cell] == "s2" {
        "grp1"
      } else {
        "grp2"
      },
    )
  }
  let rows = @src.propeller(
    clusters,
    samples,
    groups,
    transform="asin",
    robust=false,
  )
  assert_eq(rows.length(), 3)
  assert_eq(rows[0].cluster, "c0")
  assert_eq(rows[1].cluster, "c2")
  assert_eq(rows[2].cluster, "c1")
  for r in rows {
    assert_eq(r.kind, "t")
  }
  spkt_vec(
    [rows[0].baseline_prop, rows[1].baseline_prop, rows[2].baseline_prop],
    [0.46739130434782611, 0.21304347826086956, 0.31956521739130433],
    tolerance=1.0e-12,
  )
  assert_true(spkt_approx(rows[0].statistic, 3.0817862495555488))
  assert_true(spkt_approx(rows[0].p_value, 0.021613588088540914))
  assert_true(spkt_approx(rows[0].fdr, 0.044924613328430098))
  spkt_vec(
    rows[0].means,
    [0.54999999999999993, 0.34999999999999998],
    tolerance=1.0e-12,
  )
  match rows[0].prop_ratio {
    Some(rr) => assert_true(spkt_approx(rr, 1.5714285714285712))
    None => assert_true(false)
  }
  assert_true(spkt_approx(rows[1].statistic, -2.8302039949815678))
}

// ---------------------------------------------------------------------------
// Fixture B: three-group moderated ANOVA
// ---------------------------------------------------------------------------

///|
test "spk_anova_fixture_b_asin_matches_limma" {
  let tab = spkt_counts_b()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let model = @src.propeller_anova_model(
    pl.transformed,
    pl.proportions,
    spkt_design_b(),
    robust=true,
  )
  assert_true(spkt_approx(model.s2_prior, 0.0047413844078486207))
  spkt_vec(model.df_total, [9.0, 9.0, 9.0], tolerance=1.0e-12)
  spkt_vec(
    model.statistic,
    [20.932450662069975, 6.6828646262656726, 10.036768558582171],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.p_value,
    [8.1124528327663366e-10, 0.00125218578559137, 4.3760955784341918e-05],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.fdr,
    [2.4337358498299009e-09, 0.00125218578559137, 6.5641433676512877e-05],
    tolerance=1.0e-10,
  )
  spkt_mat(
    model.means,
    [
      [0.54999999999999993, 0.34999999999999998, 0.77500000000000002],
      [0.29999999999999999, 0.34999999999999998, 0.14999999999999999],
      [0.14999999999999997, 0.29999999999999993, 0.074999999999999997],
    ],
    tolerance=1.0e-12,
  )
}

///|
test "spk_anova_fixture_b_logit_matches_limma" {
  let tab = spkt_counts_b()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="logit")
  let model = @src.propeller_anova_model(
    pl.transformed,
    pl.proportions,
    spkt_design_b(),
    robust=false,
  )
  assert_true(spkt_approx(model.s2_prior, 0.13850716957765846))
  spkt_vec(
    model.statistic,
    [12.60323722962165, 5.5048510558604917, 10.667836629194111],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.p_value,
    [3.361116896267234e-06, 0.004066994290735408, 2.3281846314724094e-05],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.fdr,
    [1.0083350688801703e-05, 0.004066994290735408, 3.4922769472086141e-05],
    tolerance=1.0e-10,
  )
}

///|
test "spk_anova_rows_sorted" {
  let tab = spkt_counts_b()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let rows = @src.propeller_anova(pl, spkt_design_b(), robust=true)
  assert_eq(rows[0].cluster, "c0")
  assert_eq(rows[1].cluster, "c2")
  assert_eq(rows[2].cluster, "c1")
  assert_true(rows[0].p_value < rows[1].p_value)
  assert_true(rows[1].p_value < rows[2].p_value)
  assert_eq(rows[0].means.length(), 3)
}

// ---------------------------------------------------------------------------
// Fixture C: six types, eight samples, zeros, finite (non-Inf) prior df
// ---------------------------------------------------------------------------

///|
test "spk_ttest_fixture_c_asin_matches_limma" {
  let tab = spkt_counts_c()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="K",
    sample_prefix="u",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  spkt_mat(
    pl.counts,
    [
      [420.0, 510.0, 380.0, 460.0, 80.0, 120.0, 60.0, 100.0],
      [300.0, 260.0, 330.0, 290.0, 520.0, 480.0, 550.0, 500.0],
      [150.0, 180.0, 140.0, 170.0, 210.0, 230.0, 190.0, 220.0],
      [90.0, 70.0, 110.0, 80.0, 140.0, 120.0, 160.0, 130.0],
      [40.0, 55.0, 35.0, 45.0, 90.0, 75.0, 105.0, 85.0],
      [0.0, 5.0, 0.0, 3.0, 60.0, 45.0, 70.0, 55.0],
    ],
    tolerance=1.0e-12,
  )
  let model = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    spkt_design_c(),
    [1.0, -1.0],
    robust=true,
  )
  assert_true(
    spkt_approx(model.df_prior[0], 73.283460728548178, tolerance=1.0e-9),
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.00086074673049024535, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [36.0, 36.0, 36.0, 36.0, 36.0, 36.0],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.statistic,
    [
      19.651473266353097, -8.9231413062557845, -2.5079962525562927, -3.177600742167487,
      -4.0049905423468308, -9.6417492045473487,
    ],
    tolerance=1.0e-11,
  )
  spkt_vec(
    model.p_value,
    [
      7.8261404668450984e-21, 1.1896866180887529e-10, 0.016788881594774225, 0.0030452553795657293,
      0.00029712444394838056, 1.6314700511748352e-11,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.fdr,
    [
      4.695684280107059e-20, 2.3793732361775057e-10, 0.016788881594774225, 0.0036543064554788751,
      0.00044568666592257082, 4.8944101535245055e-11,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.prop_ratio,
    [
      5.1992399961299887, 0.61627304769185842, 0.79936193936261701, 0.68387539769747763,
      0.52431290258133412, 0.035893605941521099,
    ],
    tolerance=1.0e-11,
  )
  spkt_mat(
    model.means,
    [
      [0.42826576691770979, 0.082370840206739043],
      [0.28727914736248383, 0.4661556244240066],
      [0.1548959811781554, 0.19377452634492953],
      [0.085425864124126449, 0.12491436950611849],
      [0.042260184155155193, 0.08060107608852822],
      [0.0018730562623692393, 0.052183563429678163],
    ],
    tolerance=1.0e-11,
  )
}

///|
test "spk_ttest_fixture_c_logit_matches_limma" {
  let tab = spkt_counts_c()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="K",
    sample_prefix="u",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="logit")
  let model = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    spkt_design_c(),
    [1.0, -1.0],
    robust=true,
  )
  assert_true(
    spkt_approx(model.df_prior[0], 1.844515520555013, tolerance=1.0e-10),
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.02973044875007123, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [
      7.8445155205550128, 7.8445155205550128, 7.8445155205550128, 7.8445155205550128,
      7.8445155205550128, 7.8445155205550128,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.statistic,
    [
      12.159034005159016, -7.389379342791103, -2.9394788920953312, -3.2917895391660403,
      -6.4415694039685594, -6.6650896920812439,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.p_value,
    [
      2.2903816571394258e-06, 8.513107380428836e-05, 0.019131977813676111, 0.011297598054278017,
      0.0002177843846437879, 0.00017297947223874011,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.fdr,
    [
      1.3742289942836555e-05, 0.00025539322141286507, 0.019131977813676111, 0.013557117665133621,
      0.00032667657696568186, 0.00032667657696568186,
    ],
    tolerance=1.0e-9,
  )
}

// ---------------------------------------------------------------------------
// Fixture D: five types, nine samples, three groups (two tested coefficients)
// ---------------------------------------------------------------------------

///|
test "spk_anova_fixture_d_asin_matches_limma" {
  let tab = spkt_counts_d()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="T",
    sample_prefix="v",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let model = @src.propeller_anova_model(
    pl.transformed,
    pl.proportions,
    spkt_design_d(),
    robust=true,
  )
  assert_true(
    spkt_approx(model.df_prior[0], 9.7475933924571603, tolerance=1.0e-10),
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.00076707780663314029, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [
      15.74759339245716, 15.74759339245716, 15.74759339245716, 15.74759339245716,
      15.74759339245716,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.statistic,
    [
      222.97332809691127, 122.23116218087978, 219.29323181544677, 137.03312741475563,
      75.202057739841152,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.p_value,
    [
      2.8056646948634276e-12, 2.56361167235081e-10, 3.1841121284556862e-12, 1.0974775749515196e-10,
      8.7663276588833175e-09,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.fdr,
    [
      7.9602803211392155e-12, 3.2045145904385123e-10, 7.9602803211392155e-12, 1.829129291585866e-10,
      8.7663276588833175e-09,
    ],
    tolerance=1.0e-9,
  )
  spkt_mat(
    model.means,
    [
      [0.27490097977902861, 0.49697114618205024, 0.10895496881141591],
      [0.49167187825724418, 0.15498508346428005, 0.29686741809274497],
      [0.15433083176985618, 0.06622030926191616, 0.39590521984779875],
      [0.07492964352720452, 0.24884994648265812, 0.039565923267538237],
      [0.0041666666666666683, 0.032973514609095679, 0.15870646998050231],
    ],
    tolerance=1.0e-11,
  )
}

///|
test "spk_anova_fixture_d_logit_matches_limma" {
  let tab = spkt_counts_d()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="T",
    sample_prefix="v",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="logit")
  let model = @src.propeller_anova_model(
    pl.transformed,
    pl.proportions,
    spkt_design_d(),
    robust=true,
  )
  assert_true(
    spkt_approx(model.df_prior[0], 1.4410231620951492, tolerance=1.0e-10),
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.017267190800887941, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [
      7.4410231620951492, 7.4410231620951492, 7.4410231620951492, 7.4410231620951492,
      7.4410231620951492,
    ],
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.statistic,
    [
      263.53661928404352, 94.794332184018543, 260.1633339790643, 206.59203588777513,
      22.00693394982488,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.p_value,
    [
      1.2402671297240438e-07, 5.0826212332502979e-06, 1.3002876155123735e-07, 3.0247326160902188e-07,
      0.00075082145192879694,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    model.fdr,
    [
      3.2507190387809339e-07, 6.3532765415628722e-06, 3.2507190387809339e-07, 5.0412210268170312e-07,
      0.00075082145192879694,
    ],
    tolerance=1.0e-8,
  )
}

// ---------------------------------------------------------------------------
// Fixture E: robust variance shrinkage actually distinguishes outlier rows
// ---------------------------------------------------------------------------

///|
fn spkt_etrans() -> @src.Matrix {
  @src.Matrix::from_2d(14, 8, [
    [
      -0.58958821854162957, -0.6215938152470456, -0.59721523769961449, -0.60169497569897146,
      -0.56333279230569189, -0.60032178064018904, -0.56470293594913346, -0.57040244526270178,
    ],
    [
      -0.48234429578943844, -0.49409120425218644, -0.49277967870214312, -0.49922405018725685,
      -0.43904411660350806, -0.43913171510396359, -0.4855530135379505, -0.40767535462944393,
    ],
    [
      -0.15348791818328003, -0.29308325082300851, -0.65088235987086129, -0.12781939706173248,
      -0.43426597416235307, -0.37686932115481331, -0.80644087715072565, 0.1938078971897178,
    ],
    [
      -0.25288201768164498, -0.21568447715331407, -0.2192272633014074, -0.25268323710096791,
      -0.19093959061065041, -0.16933357411602779, -0.214425492219527876, -0.20051831333390283,
    ],
    [
      -0.14355600995681733, -0.13180589127278022, -0.13481288220498572, -0.11516966347117838,
      -0.076591681446981058, -0.1062971330458858, -0.076916472254390356, -0.10495965008269927,
    ],
    [
      -0.046222000229391852, -0.0075169475967987263, -0.047161358796384341, -0.016942641033969395,
      0.010249906110235889, 0.055285900193151979, 0.041158799280091377, 0.018873122465829541,
    ],
    [
      0.10829493668198248, 0.10083738101883127, 0.068909718908979872, 0.088183689101831941,
      0.12369548612282007, 0.15129564362352704, 0.12941156321644603, 0.13767923851147315,
    ],
    [
      0.18296494746340033, 0.18847320172784743, 0.21037143332369895, 0.18771125323292207,
      0.25486289387604344, 0.26103890568528398, 0.25371653543887412, 0.25186816393046579,
    ],
    [
      0.28838898239193739, 0.3176223564422711, 0.28707833344298667, 0.34613174261298391,
      0.35301053336508303, 0.37604112396277389, 0.38346685669190667, 0.37318564896643835,
    ],
    [
      0.46540195515461735, 0.42152087011561867, 0.447329498807854, 0.41891854909507231,
      0.53095136117181063, 0.5022155176727694, 0.48159417405793065, 0.50416492791501744,
    ],
    [
      0.26495178865352459, 0.93107524044432854, 0.26497211378181507, 0.16464320129448295,
      0.59003417254026491, 0.81614444452183288, 0.6010719732202392, 0.4137676959009744,
    ],
    [
      0.69274675426429422, 0.71751976341146462, 0.63006051386624673, 0.6987025670993654,
      0.70973309000079488, 0.73885418890648691, 0.68274366484003823, 0.71401855288372407,
    ],
    [
      0.76543847836391143, 0.77480947041405412, 0.76306945034337148, 0.79200390282627353,
      0.85894071423817164, 0.82473851990310254, 0.83005969019584092, 0.81238240687443919,
    ],
    [
      0.91991666178819254, 0.91123657781892176, 0.86871442886968242, 0.9296215403821062,
      0.92661589699775537, 0.92895972925067294, 0.92056308950341958, 0.92613126158471637,
    ],
  ])
}

///|
fn spkt_eprops() -> @src.Matrix {
  @src.Matrix::from_2d(14, 8, [
    [
      0.35672934166995712, 0.34941904987452843, 0.35498106161804838, 0.35395600585304571,
      0.36277666713287071, 0.35427007888576456, 0.36245999094950343, 0.36114396789387315,
    ],
    [
      0.38169870693448105, 0.37893025918030687, 0.37923896535486501, 0.37772303720453643,
      0.39196876068744468, 0.39194788359999427, 0.38094172212844213, 0.39946965826939107,
    ],
    [
      0.46170317572549546, 0.42724920433760044, 0.34279072710215197, 0.46808858576274853,
      0.39310811775370375, 0.40688220009245757, 0.30864944132459632, 0.5483008810608001,
    ],
    [
      0.43711426391713937, 0.44628694612793651, 0.44541163833666469, 0.43716317357271534,
      0.45240960184053408, 0.45776747240685645, 0.44659808107888843, 0.45003771486402205,
    ],
    [
      0.46417250511823493, 0.4670961493767512, 0.46634773181623207, 0.47123936731578814,
      0.48086143474788251, 0.47345071578673697, 0.48078035651662454, 0.47378415036529503,
    ],
    [
      0.48844655683798638, 0.49812077194952897, 0.48821184514828853, 0.495764441106021134,
      0.50256245409318989, 0.51381795564369792, 0.51028824746333978, 0.50471814056921305,
    ],
    [
      0.52704730556843893, 0.5251880058614723, 0.51722061585446044, 0.52203164695179816,
      0.5308845023531501, 0.53775192539699967, 0.53230781418792472, 0.53436554189278818,
    ],
    [
      0.54561405913185701, 0.54697931551282619, 0.55239974983465734, 0.54679050330260504,
      0.56337306066234272, 0.56489166066076058, 0.56309105453818054, 0.5626362661630625,
    ],
    [
      0.57160168348198837, 0.57874469315922239, 0.571280771060554753, 0.5856792233432212,
      0.58734743454386507, 0.59291791644142033, 0.59470899373977892, 0.59222851867284931,
    ],
    [
      0.61429488507071994, 0.60384712338736735, 0.61000411108598585, 0.60322443922275659,
      0.62970497450972229, 0.62297984464556577, 0.61812424503024754, 0.62343760435594076,
    ],
    [
      0.56585315917165224, 0.71729337754290545, 0.56585815230430603, 0.5410680715448779,
      0.64337298642352669, 0.69341730033885862, 0.64590151849393607, 0.60199095662419266,
    ],
    [
      0.66657767710764004, 0.67206061842499054, 0.65250318340157121, 0.6679000518480952,
      0.67034217977384336, 0.67674524810326464, 0.66435077961117339, 0.67128850427939524,
    ],
    [
      0.68253331998427236, 0.68456035929277792, 0.68201977340191244, 0.68826144534053324,
      0.702439291690374, 0.69524126955181875, 0.6963675085393776, 0.6926169482769331,
    ],
    [
      0.71502512466159118, 0.71325313815845925, 0.70447812699494972, 0.71699849793383807,
      0.716388225029380245, 0.71686418946050767, 0.71515682517438384, 0.71628974361692388,
    ],
  ])
}

///|
fn spkt_design_e() -> @src.Matrix {
  @src.Matrix::from_2d(8, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [0.0, 1.0],
    [0.0, 1.0],
    [0.0, 1.0],
  ])
}

///|
test "spk_robust_fixture_e_nonrobust_matches_limma" {
  let model = @src.propeller_ttest_model(
    spkt_etrans(),
    spkt_eprops(),
    spkt_design_e(),
    [1.0, -1.0],
    robust=false,
  )
  spkt_vec(
    model.df_prior,
    Array::make(14, 1.1125003809951362),
    tolerance=1.0e-10,
  )
  assert_true(
    spkt_approx(model.s2_prior, 0.00028776224603485959, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    Array::make(14, 7.1125003809951366),
    tolerance=1.0e-10,
  )
  spkt_vec(
    model.s2_post,
    [
      0.00025083964324617458, 0.00050047116694901489, 0.096475982588783232, 0.0003728288491377978,
      0.00022183442342615771, 0.00039685614024864799, 0.00023077166294742233, 0.00011485710502368739,
      0.00044773754457483198, 0.00042500260635933179, 0.064035279951435098, 0.00087565445402206573,
      0.00028173725989884193, 0.00035428571305575938,
    ],
    tolerance=1.0e-12,
  )
  spkt_vec(
    model.statistic,
    [
      -2.485295727903023, -3.1139305176140328, 0.22594135926933251, -3.0259940200543536,
      -3.8118049128072338, -4.319943595426512, -4.0928088828622924, -8.3122277603627825,
      -4.1184118732108139, -4.5576485204209281, -1.1112654213393229, -1.2702908715721535,
      -4.8614824560758425, -1.3670828141970572,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.p_value,
    [
      0.041362021825932231, 0.016645535069548661, 0.82760270444733441, 0.018857877999684948,
      0.006417022207661955, 0.003353080094757023, 0.0044603174499085472, 6.5219179675231074e-05,
      0.004317466313456966, 0.0025080326473105236, 0.30260226228941123, 0.2439552217566367,
      0.001751335992502322, 0.21323618806986264,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.fdr,
    [
      0.057906830556305119, 0.029129686371710158, 0.82760270444733441, 0.029334476888398808,
      0.01283404441532391, 0.010407407383119944, 0.010407407383119944, 0.00091306851545323506,
      0.010407407383119944, 0.010407407383119944, 0.32587935938859669, 0.28461442538274284,
      0.010407407383119944, 0.2713915120889161,
    ],
    tolerance=1.0e-9,
  )
}

///|
test "spk_robust_fixture_e_shrinks_outlier_rows" {
  let model = @src.propeller_ttest_model(
    spkt_etrans(),
    spkt_eprops(),
    spkt_design_e(),
    [1.0, -1.0],
    robust=true,
  )
  let want_dfprior = [
    1.2088783448967213, 1.2088783448967213, 0.63585129118644823, 1.2088783448967213,
    1.2088783448967213, 1.2088783448967213, 1.2088783448967213, 1.2088783448967213,
    1.2088783448967213, 1.2088783448967213, 0.63585129118644823, 1.2088783448967213,
    1.2088783448967213, 1.2088783448967213,
  ]
  spkt_vec(model.df_prior, want_dfprior, tolerance=1.0e-9)
  assert_true(
    spkt_approx(model.s2_prior, 0.00030907557410523429, tolerance=1.0e-12),
  )
  spkt_vec(
    model.df_total,
    [
      7.2088783448967213, 7.2088783448967213, 6.635851291186448, 7.2088783448967213,
      7.2088783448967213, 7.2088783448967213, 7.2088783448967213, 7.2088783448967213,
      7.2088783448967213, 7.2088783448967213, 6.635851291186448, 7.2088783448967213,
      7.2088783448967213, 7.2088783448967213,
    ],
    tolerance=1.0e-9,
  )
  spkt_vec(
    model.s2_post,
    [
      0.00025490736982411503, 0.00050120148417850084, 0.10338716515925069, 0.00037526566023633799,
      0.00022628993072323737, 0.00039897172226208735, 0.00023510768522800568, 0.00012074282873419572,
      0.00044917287566506315, 0.00042674188862612707, 0.0686162669363891, 0.00087136881737245723,
      0.00028539190565617128, 0.00035697043369229336,
    ],
    tolerance=1.0e-11,
  )
  spkt_vec(
    model.statistic,
    [
      -2.4653862220446254, -3.1116609850789909, 0.21825893350561118, -3.0161532776979088,
      -3.7740923275222293, -4.3084749348300937, -4.0548920082000661, -8.1071030687471364,
      -4.1118264193026111, -4.5483511780836485, -1.0735293286954379, -1.2734108637214883,
      -4.8302547955150379, -1.3619323003307751,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    model.p_value,
    [
      0.042157530365192139, 0.016416154444466053, 0.83378662354809963, 0.018820222500862443,
      0.0065786730605406039, 0.0032961666835248002, 0.004551083099081866, 7.1260607926884321e-05,
      0.0042295462826059645, 0.0024508201653092453, 0.3205128309316696, 0.24239422305087643,
      0.0017488555885717791, 0.21424689411231851,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    model.fdr,
    [
      0.059020542511268988, 0.028728270277815594, 0.83378662354809963, 0.029275901668008246,
      0.013157346121081208, 0.010619193897857688, 0.010619193897857688, 0.00099764851097638042,
      0.010619193897857688, 0.010619193897857688, 0.34516766408025956, 0.28279326022602252,
      0.010619193897857688, 0.2726778652338599,
    ],
    tolerance=1.0e-8,
  )
}

// ---------------------------------------------------------------------------
// Beta moments, beta-binomial counts parameters and normCounts
// ---------------------------------------------------------------------------

///|
test "spk_estimate_beta_param_matches_speckle" {
  let xp = [
    0.094930328209210837, 0.14542877646897814, 0.17960197316287255, 0.11606734068730769,
    0.11840051102092848, 0.37613389739568731, 0.22039004141198454, 0.26264844139609522,
    0.29826180573142219, 0.2922928717115903, 0.31896000215301912, 0.12383869621468634,
  ]
  let mom = @src.estimate_beta_param(xp)
  assert_true(spkt_approx(mom.alpha, 3.6981234846105142, tolerance=1.0e-11))
  assert_true(spkt_approx(mom.beta, 13.725618707136791, tolerance=1.0e-11))
}

///|
test "spk_estimate_beta_params_from_counts_matches_speckle" {
  let p = @src.estimate_beta_params_from_counts(spkt_counts_c())
  assert_true(spkt_approx(p.n, 1075.0, tolerance=1.0e-9))
  spkt_vec(
    p.alpha,
    [
      1.3343277690175506, 10.112285721814995, 55.745904630659297, 19.493906614204978,
      9.6920858523974047, 1.090534240846722,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    p.beta,
    [
      3.8918066302828023, 16.730876023527522, 264.01675626470444, 165.86204053638281,
      148.08076839269322, 39.257317587822371,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    p.pi,
    [
      0.25531830356222435, 0.37671738589324522, 0.17433525376154249, 0.10517011681512244,
      0.061430630121841703, 0.027028309846023696,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    p.dispersion,
    [
      0.1913460166913952, 0.037253435697585338, 0.0031273194850202695, 0.005395025168453727,
      0.0063382259564535746, 0.02478446694625392,
    ],
    tolerance=1.0e-8,
  )
  spkt_vec(
    p.variance,
    [
      35461.514751964023, 9988.7322348751877, 672.84329522745998, 684.21289959376406,
      481.24543376580857, 762.57917506783804,
    ],
    tolerance=1.0e-6,
  )
}

///|
test "spk_norm_counts_matches_speckle" {
  let tab = spkt_counts_c()
  let normc = @src.propeller_norm_counts(tab)
  spkt_mat(
    normc,
    [
      [
        451.5, 507.63888888888886, 410.55276381909545, 471.85114503816794, 78.181818181818173,
        120.5607476635514, 56.828193832599119, 98.623853211009177,
      ],
      [
        322.5, 258.7962962962963, 356.5326633165829, 297.4713740458015, 508.18181818181819,
        482.24299065420558, 520.92511013215858, 493.1192660550459,
      ],
      [
        161.25, 179.16666666666666, 151.25628140703517, 174.37977099236642, 205.22727272727275,
        231.07476635514018, 179.95594713656388, 216.97247706422019,
      ],
      [
        96.75, 69.675925925925924, 118.84422110552765, 82.061068702290072, 136.81818181818181,
        120.5607476635514, 151.54185022026434, 128.21100917431193,
      ],
      [
        43.0, 54.745370370370367, 37.814070351758794, 46.159351145038173, 87.954545454545453,
        75.350467289719617, 99.449339207048453, 83.83027522935781,
      ],
      [
        0.0, 4.9768518518518512, 0.0, 3.0772900763358777, 58.636363636363633, 45.210280373831772,
        66.29955947136564, 54.243119266055047,
      ],
    ],
    tolerance=1.0e-10,
  )
  let lnormc = @src.propeller_norm_counts(tab, log=true)
  spkt_mat(
    lnormc,
    [
      [
        8.8188400001424583, 8.9877384109066956, 8.681842351842155, 8.8823763210595033,
        6.296692319049729, 6.9182903904944686, 5.8399397211793129, 6.6298760933373222,
      ],
      [
        8.334057528119045, 8.0171302443211534, 8.4785773799146398, 8.2176989074036477,
        8.9892790444753565, 8.9137715607634078, 9.0249752420198774, 8.945914636424229,
      ],
      [
        7.3363105074804649, 7.4878632864338721, 7.2443038552824612, 7.4489049889978967,
        7.6832686791044358, 7.8540099088087523, 7.4941870348516462, 7.7633666228541776,
      ],
      [
        6.6023434232209546, 6.1316489142011124, 6.8976887591311691, 6.3661195222826574,
        7.100075093763035, 6.9182903904944686, 7.2470160232291834, 7.0066912734345399,
      ],
      [
        5.4417487978485442, 5.7865530309461271, 5.2586293023400401, 5.5428888339478597,
        6.4655885541884714, 6.2438233859121697, 6.64184027305071709, 6.3967067751072424,
      ],
      [
        -0.98753239548711369, 2.4532598878453991, -0.98753239548711369, 1.83945266210009,
        5.8847367243957374, 5.5132443466997394, 6.0605163096236057, 5.7733782560066107,
      ],
    ],
    tolerance=1.0e-9,
  )
  let lnormc2 = @src.propeller_norm_counts(tab, log=true, prior_count=2.0)
  spkt_mat(
    lnormc2,
    [
      [
        8.8196102987467349, 8.9879763705257503, 8.6830924729334278, 8.8829390415414302,
        6.3201399167357577, 6.9321701214495075, 5.8734982368145303, 6.647694859436772,
      ],
      [
        8.3367505738290415, 8.0214763842751307, 8.480628529299004, 8.2209573540068561,
        8.9895124285656998, 8.9142350257977832, 9.0251039646758073, 8.9462786880818541,
      ],
      [
        7.345699629687226, 7.4959173473369027, 7.2545739376433538, 7.4572890810456949,
        7.6897944187120064, 7.8593584905817009, 7.5021883535014044, 7.7693228913099244,
      ],
      [
        6.6205804086414526, 6.1583945212703979, 6.9118248204348172, 6.388286599961722,
        7.1118431183280917, 6.9321701214495075, 7.2572593329200537, 7.0195112362988281,
      ],
      [
        5.48706053435538, 5.8215105125121873, 5.3104942681600535, 5.5849110108228608,
        6.486023341284989, 6.2682874693918391, 6.6594797405720909, 6.4183286167672264,
      ],
      [
        1.0084122398162001, 2.8011422925555336, 1.0084122398162001, 2.3441094239640612,
        5.9171596563416049, 5.5562076971231669, 6.0888013100986491, 5.8086887978716186,
      ],
    ],
    tolerance=1.0e-9,
  )
}

// ---------------------------------------------------------------------------
// Error paths and edge cases
// ---------------------------------------------------------------------------

///|
test "spk_error_empty_and_mismatched_inputs" {
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.get_transformed_props([], []) |> ignore
    }),
  )
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.get_transformed_props(["c0"], []) |> ignore
    }),
  )
  let tab = spkt_counts_a()
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let one_group = Array::make(clusters.length(), "g1")
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.get_transformed_props(clusters, samples, transform="sqrt") |> ignore
    }),
  )
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.propeller(clusters, samples, one_group, transform="asin") |> ignore
    }),
  )
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.propeller_ttest_model(tab, tab, spkt_design_b(), [1.0, -1.0, 0.0])
      |> ignore // sample / design mismatch
    }),
  )
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.propeller_ttest_model(tab, tab, spkt_design_a(), [1.0]) |> ignore // contrast length
    }),
  )
  assert_true(
    spkt_raises(fn() raise @src.PropellerError {
      @src.estimate_beta_param([0.3]) |> ignore
    }),
  )
}

///|
test "spk_two_types_forces_nonrobust_without_error" {
  let tab = @src.Matrix::from_2d(2, 4, [
    [500.0, 900.0, 270.0, 480.0],
    [500.0, 600.0, 630.0, 720.0],
  ])
  let (clusters, samples) = spkt_expand(
    tab,
    cluster_prefix="c",
    sample_prefix="s",
  )
  let pl = @src.get_transformed_props(clusters, samples, transform="asin")
  let mr = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    spkt_design_a(),
    [1.0, -1.0],
    robust=true,
  )
  let mn = @src.propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    spkt_design_a(),
    [1.0, -1.0],
    robust=false,
  )
  spkt_vec(mr.statistic, mn.statistic, tolerance=1.0e-12)
  spkt_vec(mr.p_value, mn.p_value, tolerance=1.0e-12)
}

///|
test "spk_levels_use_c_locale_lexicographic_order" {
  // R table()/sort() order is pure lexicographic (C locale), where "s10"
  // sorts before "s2"; MoonBit's native String order is length-first.
  let clusters = [
    "B cell", "B cell", "T cell", "Monocyte", "Monocyte", "B cell", "T cell", "T cell",
    "Monocyte",
  ]
  let samples = ["s1", "s1", "s1", "s10", "s10", "s10", "s2", "s2", "s2"]
  let pl = @src.get_transformed_props(clusters, samples)
  assert_eq(pl.cluster_names, ["B cell", "Monocyte", "T cell"])
  assert_eq(pl.sample_names, ["s1", "s10", "s2"])
  // rows = cluster_names, columns = sample_names
  assert_true(spkt_approx(pl.counts.get(0, 0), 2.0))
  assert_true(spkt_approx(pl.counts.get(0, 1), 1.0))
  assert_true(spkt_approx(pl.counts.get(0, 2), 0.0))
  assert_true(spkt_approx(pl.counts.get(1, 0), 0.0))
  assert_true(spkt_approx(pl.counts.get(1, 1), 2.0))
  assert_true(spkt_approx(pl.counts.get(1, 2), 1.0))
  assert_true(spkt_approx(pl.counts.get(2, 0), 1.0))
  assert_true(spkt_approx(pl.counts.get(2, 1), 0.0))
  assert_true(spkt_approx(pl.counts.get(2, 2), 2.0))
}
