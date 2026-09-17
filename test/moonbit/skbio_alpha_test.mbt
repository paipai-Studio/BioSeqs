///|
/// Tests for the scikit-bio alpha diversity metrics port. Expected values come
/// from running scikit-bio 0.7.3's `skbio.diversity.alpha` functions on the
/// same vectors.

///|
fn skb_close(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-6
}

///|
test "skbio_alpha_counts_and_richness" {
  // singletons / doubletons
  let c = [1.0, 1.0, 2.0, 2.0, 5.0, 8.0]
  assert_eq(@src.skbio_alpha_sobs(c), 6)
  assert_eq(@src.skbio_alpha_singles(c), 2)
  assert_eq(@src.skbio_alpha_doubles(c), 2)
  let osd = @src.skbio_alpha_osd(c)
  assert_eq(osd.0, 6)
  assert_eq(osd.1, 2)
  assert_eq(osd.2, 2)

  // zero counts are dropped
  let z = [1.0, 2.0, 0.0, 3.0, 0.0]
  assert_eq(@src.skbio_alpha_sobs(z), 3)
  assert_eq(@src.skbio_alpha_singles(z), 1)
  assert_eq(@src.skbio_alpha_doubles(z), 1)

  // single taxon
  let s = [7.0]
  assert_eq(@src.skbio_alpha_sobs(s), 1)
  assert_eq(@src.skbio_alpha_singles(s), 0)
  assert_eq(@src.skbio_alpha_doubles(s), 0)
}

///|
test "skbio_alpha_dominance_evenness" {
  let c = [1.0, 2.0, 3.0, 4.0]
  assert_true(skb_close(@src.skbio_alpha_dominance(c), 0.3))
  assert_true(skb_close(@src.skbio_alpha_simpson_d(c), 0.3))
  assert_true(skb_close(@src.skbio_alpha_inv_simpson(c), 3.3333333333))
  assert_true(skb_close(@src.skbio_alpha_enspie(c), 3.3333333333))
  assert_true(skb_close(@src.skbio_alpha_berger_parker_d(c), 0.4))
  assert_true(skb_close(@src.skbio_alpha_simpson_e(c), 0.8333333333))
  assert_true(skb_close(@src.skbio_alpha_heip_e(c), 0.8653718222))
  assert_true(skb_close(@src.skbio_alpha_goods_coverage(c), 0.9))
}

///|
test "skbio_alpha_richness_indices" {
  let c = [1.0, 2.0, 3.0, 4.0]
  assert_true(skb_close(@src.skbio_alpha_margalef(c), 1.3028834457))
  assert_true(skb_close(@src.skbio_alpha_menhinick(c), 1.26491106407))
  assert_true(skb_close(@src.skbio_alpha_mcintosh_d(c), 0.6614445864))
  assert_true(skb_close(@src.skbio_alpha_mcintosh_e(c), 0.7595545253))
  assert_true(skb_close(@src.skbio_alpha_brillouin_d(c), 0.9441452093))
  assert_true(skb_close(@src.skbio_alpha_strong(c), 0.2))
}

///|
test "skbio_alpha_generalized_entropy" {
  let c = [1.0, 2.0, 3.0, 4.0]
  assert_true(skb_close(@src.skbio_alpha_renyi(c, 2.0), 1.2039728043))
  assert_true(skb_close(@src.skbio_alpha_renyi(c, 0.0), 1.3862943611))
  assert_true(skb_close(@src.skbio_alpha_renyi(c, 1.0), 1.2798542258))
  assert_true(skb_close(@src.skbio_alpha_tsallis(c, 2.0), 0.7))
  assert_true(skb_close(@src.skbio_alpha_tsallis(c, 1.0), 1.2798542258))
  assert_true(skb_close(@src.skbio_alpha_hill(c, 2.0), 3.3333333333))
  assert_true(skb_close(@src.skbio_alpha_hill(c, 0.0), 4.0))
  assert_true(skb_close(@src.skbio_alpha_hill(c, 1.0), 3.5961154666))
}

///|
test "skbio_alpha_single_taxon_edge" {
  let s = [7.0]
  assert_true(skb_close(@src.skbio_alpha_dominance(s), 1.0))
  assert_true(skb_close(@src.skbio_alpha_simpson_e(s), 1.0))
  assert_true(skb_close(@src.skbio_alpha_heip_e(s), 1.0))
  assert_true(skb_close(@src.skbio_alpha_mcintosh_e(s), 1.0))
  assert_true(skb_close(@src.skbio_alpha_mcintosh_d(s), 0.0))
  assert_true(skb_close(@src.skbio_alpha_brillouin_d(s), 0.0))
  assert_true(skb_close(@src.skbio_alpha_strong(s), 0.0))
  assert_true(skb_close(@src.skbio_alpha_renyi(s, 2.0), 0.0))
  assert_true(skb_close(@src.skbio_alpha_tsallis(s, 2.0), 0.0))
  assert_true(skb_close(@src.skbio_alpha_hill(s, 2.0), 1.0))
}