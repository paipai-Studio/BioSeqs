///|
/// Tests for the Biopython Bio.motifs.thresholds ScoreDistribution port.
/// All expected values were frozen from Biopython 1.88
/// (generator: /tmp/gen_motif_threshold_expected.py).

///|
fn mscd_approx(left : Double, right : Double, tolerance~ : Double = 1.0e-9) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn mscd_infinity() -> Double {
  1.0e308 * 10.0
}

///|
fn mscd_nan() -> Double {
  mscd_infinity() - mscd_infinity()
}

///|
fn mscd_raises(operation : () -> Unit raise @src.MotifThresholdError) -> Bool {
  try {
    operation()
    false
  } catch {
    _ => true
  }
}

///|
fn mscd_uniform_pssm() -> Array[Array[Double]] {
  [
    [
      1.4854268271702415,
      -1.3219280948873622,
      1.4854268271702415,
      -0.3219280948873623,
      1.2630344058337937,
      -1.3219280948873622,
      1.0,
      0.2630344058337938,
    ],
    [
      -1.3219280948873622,
      1.2630344058337937,
      -1.3219280948873622,
      -1.3219280948873622,
      -0.3219280948873623,
      1.4854268271702415,
      -0.3219280948873623,
      0.6780719051126378,
    ],
    [
      -1.3219280948873622,
      -0.3219280948873623,
      -1.3219280948873622,
      1.2630344058337937,
      -1.3219280948873622,
      -1.3219280948873622,
      -0.3219280948873623,
      -0.3219280948873623,
    ],
    [
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
      -1.3219280948873622,
    ],
  ]
}

///|
fn mscd_biased_pssm() -> Array[Array[Double]] {
  [
    [
      1.222392421336448,
      -1.584962500721156,
      1.222392421336448,
      -0.5849625007211561,
      1.0,
      -1.584962500721156,
      0.7369655941662062,
      0.0,
    ],
    [
      -1.0,
      1.584962500721156,
      -1.0,
      -1.0,
      0.0,
      1.807354922057604,
      0.0,
      1.0,
    ],
    [
      -1.0,
      0.0,
      -1.0,
      1.584962500721156,
      -1.0,
      -1.0,
      0.0,
      0.0,
    ],
    [
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
      -1.584962500721156,
    ],
  ]
}

///|
fn mscd_counts8() -> Array[Array[Int]] {
  [
    [7, 1, 7, 2, 6, 1, 5, 3],
    [1, 6, 1, 1, 2, 7, 2, 4],
    [1, 2, 1, 6, 1, 1, 2, 2],
    [1, 1, 1, 1, 1, 1, 1, 1],
  ]
}

///|
fn mscd_argmax(values : Array[Double]) -> Int {
  let mut best = 0
  for i in 1..<values.length() {
    if values[i] > values[best] {
      best = i
    }
  }
  best
}

test "mt uniform background grid and information content" {
  let dist = @src.MotifScoreDistribution::from_pssm(mscd_uniform_pssm())
  inspect!(
    mscd_approx(dist.min_score(), -10.575424759098896),
    content="true"
  )
  inspect!(
    mscd_approx(dist.step(), 0.0025626803804505113),
    content="true"
  )
  inspect!(dist.n_points(), content="8000")
  inspect!(
    mscd_approx(dist.information_content(), 3.6094058761801775),
    content="true"
  )
  let mo = dist.motif_density()
  let bg = dist.background_density()
  let mo_sum = mo.fold(init=0.0, fn(acc : Double, x : Double) -> Double { acc + x })
  let bg_sum = bg.fold(init=0.0, fn(acc : Double, x : Double) -> Double { acc + x })
  inspect!(mscd_approx(mo_sum, 1.0, tolerance=1.0e-12), content="true")
  inspect!(mscd_approx(bg_sum, 1.0, tolerance=1.0e-12), content="true")
  // Frozen Biopython argmax bins and their probability mass.
  inspect!(mscd_argmax(mo), content="5380")
  inspect!(
    mscd_approx(dist.grid_score(5380), 3.2117956877248552),
    content="true"
  )
  inspect!(
    mscd_approx(mo[5380], 0.019051199999999994, tolerance=1.0e-14),
    content="true"
  )
  inspect!(mscd_argmax(bg), content="1398")
  inspect!(
    mscd_approx(dist.grid_score(1398), -6.992797587229081),
    content="true"
  )
  inspect!(
    mscd_approx(bg[1398], 0.032135009765625, tolerance=1.0e-14),
    content="true"
  )
}

test "mt uniform background thresholds" {
  let dist = @src.MotifScoreDistribution::from_pssm(mscd_uniform_pssm())
  inspect!(
    mscd_approx(dist.threshold_fpr(0.001), 6.604784511441332),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.01), 4.019040007566765),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.05), 1.6254965322259878),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.1), 0.4030979907510943),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.001), -6.5776433655960975),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.01), -3.994461542101983),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.05), -1.6393582724679625),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.1), -0.37595684490586123),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_balanced(), -0.012056230881889007),
    content="true"
  )
  let (balanced, rate) = dist.threshold_balanced_with_rate()
  inspect!(
    mscd_approx(balanced, -0.012056230881889007),
    content="true"
  )
  inspect!(mscd_approx(rate, 0.1300048828125, tolerance=1.0e-12), content="true")
  inspect!(
    mscd_approx(dist.threshold_balanced(rate_proportion=0.5), -0.639912924092263),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_balanced(rate_proportion=2.0), 0.6260511838502882),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_patser(), 0.7234330383074088),
    content="true"
  )
}

test "mt biased background thresholds" {
  let bg = [0.30, 0.20, 0.20, 0.30]
  let dist = @src.MotifScoreDistribution::from_pssm(
    mscd_biased_pssm(),
    background=bg
  )
  inspect!(
    mscd_approx(dist.min_score(), -12.679700005769249),
    content="true"
  )
  inspect!(
    mscd_approx(dist.step(), 0.0028551981955379753),
    content="true"
  )
  inspect!(
    mscd_approx(dist.information_content(), 3.8449806323944546),
    content="true"
  )
  let mo = dist.motif_density()
  let background = dist.background_density()
  inspect!(mscd_argmax(mo), content="7648")
  inspect!(
    mscd_approx(dist.grid_score(7648), 9.156855793705187),
    content="true"
  )
  inspect!(
    mscd_approx(mo[7648], 0.023461200000000012, tolerance=1.0e-14),
    content="true"
  )
  inspect!(mscd_argmax(background), content="2854")
  inspect!(
    mscd_approx(background[2854], 0.021179879999999998, tolerance=1.0e-14),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.001), 6.572901426743318),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.01), 4.028919834518984),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.1), 0.37426614423037385),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.001), -6.600983047468899),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.01), -3.9456487256185824),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.1), -0.3623749902184237),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_balanced(), 0.011655973397051156),
    content="true"
  )
  let (_, rate) = dist.threshold_balanced_with_rate()
  inspect!(mscd_approx(rate, 0.11796088, tolerance=1.0e-12), content="true")
  inspect!(
    mscd_approx(dist.threshold_balanced(rate_proportion=0.5), -0.7249851610517464),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_balanced(rate_proportion=2.0), 0.7882698825833803),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_patser(), 1.0109753418353424),
    content="true"
  )
}

test "mt pseudocount small motif" {
  let pssm = [
    [
      1.222392421336448,
      0.0,
      0.0,
      0.7369655941662062,
    ],
    [
      0.0,
      -1.5849625007211563,
      0.7369655941662062,
      0.0,
    ],
    [
      -1.5849625007211563,
      0.7369655941662062,
      -1.5849625007211563,
      0.0,
    ],
    [
      -1.5849625007211563,
      0.0,
      0.0,
      -1.5849625007211563,
    ],
  ]
  let dist = @src.MotifScoreDistribution::from_pssm(
    pssm,
    precision=500
  )
  inspect!(
    mscd_approx(dist.min_score(), -6.339850002884625),
    content="true"
  )
  inspect!(
    mscd_approx(dist.step(), 0.004889014110415053),
    content="true"
  )
  inspect!(dist.n_points(), content="2000")
  inspect!(
    mscd_approx(dist.information_content(), 0.973868196520204),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fpr(0.01), 2.699937087272807),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_fnr(0.01), -3.1668798452252562),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_balanced(), 0.0012012983236981967),
    content="true"
  )
  inspect!(
    mscd_approx(dist.threshold_patser(), -0.8445981427781062),
    content="true"
  )
  let mo = dist.motif_density()
  let bg = dist.background_density()
  inspect!(
    mscd_approx(mo.fold(init=0.0, fn(acc : Double, x : Double) -> Double { acc + x }), 1.0, tolerance=1.0e-12),
    content="true"
  )
  inspect!(
    mscd_approx(bg.fold(init=0.0, fn(acc : Double, x : Double) -> Double { acc + x }), 1.0, tolerance=1.0e-12),
    content="true"
  )
}

test "mt from_counts reproduces pssm-based distribution" {
  let from_counts = @src.MotifScoreDistribution::from_counts(mscd_counts8())
  let from_pssm = @src.MotifScoreDistribution::from_pssm(mscd_uniform_pssm())
  inspect!(
    mscd_approx(from_counts.min_score(), from_pssm.min_score()),
    content="true"
  )
  inspect!(
    mscd_approx(from_counts.step(), from_pssm.step(), tolerance=1.0e-15),
    content="true"
  )
  inspect!(
    mscd_approx(from_counts.information_content(), from_pssm.information_content()),
    content="true"
  )
  inspect!(
    mscd_approx(from_counts.threshold_fpr(0.01), 4.019040007566765),
    content="true"
  )
  inspect!(
    mscd_approx(from_counts.threshold_balanced(), -0.012056230881889007),
    content="true"
  )
  // Zero counts without pseudocount are rejected (would yield -infinite scores).
  let sparse_counts = [
    [3, 1, 1, 2],
    [1, 0, 2, 1],
    [0, 2, 0, 1],
    [0, 1, 1, 0],
  ]
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_counts(sparse_counts)
    }),
    content="true"
  )
  // With a 0.5 scalar pseudocount the same counts match the frozen Biopython
  // pseudocount distribution up to one discretization bin.
  let with_pc = @src.MotifScoreDistribution::from_counts(
    sparse_counts,
    pseudocount=0.5,
    precision=500
  )
  inspect!(
    mscd_approx(with_pc.min_score(), -6.339850002884625, tolerance=1.0e-12),
    content="true"
  )
  inspect!(
    mscd_approx(with_pc.threshold_fpr(0.01), 2.699937087272807),
    content="true"
  )
  inspect!(
    mscd_approx(with_pc.threshold_patser(), -0.8445981427781062),
    content="true"
  )
}

test "mt default background and density copies" {
  let explicit = @src.MotifScoreDistribution::from_pssm(
    mscd_uniform_pssm(),
    background=[0.25, 0.25, 0.25, 0.25]
  )
  let implicit = @src.MotifScoreDistribution::from_pssm(mscd_uniform_pssm())
  inspect!(
    mscd_approx(implicit.threshold_fpr(0.01), explicit.threshold_fpr(0.01)),
    content="true"
  )
  // Unnormalized-but-uniform background is renormalized to the same result.
  let scaled = @src.MotifScoreDistribution::from_pssm(
    mscd_uniform_pssm(),
    background=[1.0, 1.0, 1.0, 1.0]
  )
  inspect!(
    mscd_approx(scaled.threshold_fnr(0.05), explicit.threshold_fnr(0.05)),
    content="true"
  )
  // Accessors return copies.
  let mo = implicit.motif_density()
  let before = mo[0]
  mo[0] += 1.0
  inspect!(implicit.motif_density()[0] == before, content="true")
}

test "mt invalid input validation" {
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm([[1.0], [1.0], [1.0]])
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[1.0], [1.0], [1.0], [1.0], [1.0]]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[1.0, 2.0], [1.0], [1.0, 2.0], [1.0, 2.0]]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm([[], [], [], []])
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[mscd_infinity()], [1.0], [1.0], [1.0]]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[mscd_nan()], [1.0], [1.0], [1.0]]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[1.0]],
        background=[0.25, 0.25, 0.5]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[1.0]],
        background=[0.25, 0.25, -0.25, 0.75]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm(
        [[1.0]],
        background=[0.0, 0.0, 0.0, 0.0]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_pssm([[1.0]], precision=0)
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_counts(
        [[1], [1], [1], [-1]]
      )
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = @src.MotifScoreDistribution::from_counts(
        [[1], [1], [1], [1]],
        pseudocount=-0.1
      )
    }),
    content="true"
  )
  let dist = @src.MotifScoreDistribution::from_pssm(mscd_uniform_pssm())
  inspect!(
    mscd_raises(fn() {
      let _ = dist.threshold_fpr(1.5)
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = dist.threshold_fnr(-0.1)
    }),
    content="true"
  )
  inspect!(
    mscd_raises(fn() {
      let _ = dist.threshold_balanced(rate_proportion=-1.0)
    }),
    content="true"
  )
}
