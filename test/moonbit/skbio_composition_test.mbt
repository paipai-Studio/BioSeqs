///|
/// Tests for the scikit-bio composition port. Expected values are the exact
/// doctest outputs from scikit-bio 0.7.3 `stats/composition/_base.py`.

///|
test "skbio_closure" {
  let r = @src.skbio_closure([2.0, 2.0, 6.0])
  assert_true(assert_close(r[0], 0.2, 1.0e-9))
  assert_true(assert_close(r[1], 0.2, 1.0e-9))
  assert_true(assert_close(r[2], 0.6, 1.0e-9))
}

///|
test "skbio_multi_replace" {
  // doctest row: [.2, .4, .4, 0] -> [.1875, .375, .375, .0625]
  let r = @src.skbio_multi_replace([0.2, 0.4, 0.4, 0.0])
  assert_true(assert_close(r[0], 0.1875, 1.0e-9))
  assert_true(assert_close(r[1], 0.375, 1.0e-9))
  assert_true(assert_close(r[2], 0.375, 1.0e-9))
  assert_true(assert_close(r[3], 0.0625, 1.0e-9))
  // explicit delta keeps the same default (delta = (1/4)^2 = 0.0625)
  let r2 = @src.skbio_multi_replace_delta([0.2, 0.4, 0.4, 0.0], 0.0625)
  assert_true(assert_close(r2[3], 0.0625, 1.0e-9))
  // no zeros -> unchanged (and still closes to 1)
  let r3 = @src.skbio_multi_replace([0.2, 0.4, 0.4])
  assert_true(assert_close(r3[0], 0.2, 1.0e-9))
}

///|
test "skbio_perturb" {
  // doctest: perturb([1/3,1/3,1/3], [1/2,1/2,1]) = [.25,.25,.5]
  let r = @src.skbio_perturb(
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [0.5, 0.5, 1.0],
  )
  assert_true(assert_close(r[0], 0.25, 1.0e-9))
  assert_true(assert_close(r[1], 0.25, 1.0e-9))
  assert_true(assert_close(r[2], 0.5, 1.0e-9))
}

///|
test "skbio_perturb_inv" {
  // doctest: perturb_inv([.1,.3,.4,.2], [1/6,1/6,1/3,1/3]) = [.14285714,.42857143,.28571429,.14285714]
  let r = @src.skbio_perturb_inv(
    [0.1, 0.3, 0.4, 0.2],
    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0],
  )
  assert_true(assert_close(r[0], 0.14285714, 1.0e-8))
  assert_true(assert_close(r[1], 0.42857143, 1.0e-8))
  assert_true(assert_close(r[2], 0.28571429, 1.0e-8))
  assert_true(assert_close(r[3], 0.14285714, 1.0e-8))
}

///|
test "skbio_power" {
  // doctest: power([.1,.3,.4,.2], .1) = [.23059566,.25737316,.26488486,.24714631]
  let r = @src.skbio_power([0.1, 0.3, 0.4, 0.2], 0.1)
  assert_true(assert_close(r[0], 0.23059566, 1.0e-8))
  assert_true(assert_close(r[1], 0.25737316, 1.0e-8))
  assert_true(assert_close(r[2], 0.26488486, 1.0e-8))
  assert_true(assert_close(r[3], 0.24714631, 1.0e-8))
}

///|
test "skbio_inner" {
  // doctest: inner([.1,.3,.4,.2], [.2,.4,.2,.2]) = 0.2107852473...
  let r = @src.skbio_inner([0.1, 0.3, 0.4, 0.2], [0.2, 0.4, 0.2, 0.2])
  assert_true(assert_close(r, 0.2107852473, 1.0e-9))
}

///|
test "skbio_clr" {
  // clr([.1,.3,.4,.2]) = [-.7945134, .3040988, .5917809, -.1013663]
  let r = @src.skbio_clr([0.1, 0.3, 0.4, 0.2])
  assert_true(assert_close(r[0], -0.7945134, 1.0e-6))
  assert_true(assert_close(r[1], 0.3040988, 1.0e-6))
  assert_true(assert_close(r[2], 0.5917809, 1.0e-6))
  assert_true(assert_close(r[3], -0.1013663, 1.0e-6))
  // clr_inv(clr(x)) == closure(x) == x (already closed)
  let back = @src.skbio_clr_inv(r)
  assert_true(assert_close(back[0], 0.1, 1.0e-9))
  assert_true(assert_close(back[1], 0.3, 1.0e-9))
  assert_true(assert_close(back[2], 0.4, 1.0e-9))
  assert_true(assert_close(back[3], 0.2, 1.0e-9))
}

///|
test "skbio_clr_inv" {
  // doctest: clr_inv([.1,.3,.4,.2]) = [.21383822,.26118259,.28865141,.23632778]
  let r = @src.skbio_clr_inv([0.1, 0.3, 0.4, 0.2])
  assert_true(assert_close(r[0], 0.21383822, 1.0e-8))
  assert_true(assert_close(r[1], 0.26118259, 1.0e-8))
  assert_true(assert_close(r[2], 0.28865141, 1.0e-8))
  assert_true(assert_close(r[3], 0.23632778, 1.0e-8))
}

///|
test "skbio_alr" {
  // alr([.1,.3,.6], ref=0) = [ln3, ln6] = [1.0986123, 1.7917595]
  let r = @src.skbio_alr([0.1, 0.3, 0.6], 0)
  assert_true(assert_close(r[0], 1.0986123, 1.0e-6))
  assert_true(assert_close(r[1], 1.7917595, 1.0e-6))
  // round-trip
  let back = @src.skbio_alr_inv(r, 0)
  assert_true(assert_close(back[0], 0.1, 1.0e-9))
  assert_true(assert_close(back[1], 0.3, 1.0e-9))
  assert_true(assert_close(back[2], 0.6, 1.0e-9))
}

///|
test "skbio_ilr" {
  // doctest: ilr_inv([.1,.3,.6], default basis(4)) = [.34180297,.29672718,.22054469,.14092516]
  let basis = @src.skbio_gram_schmidt_basis(4)
  let r = @src.skbio_ilr_inv([0.1, 0.3, 0.6], basis)
  assert_true(assert_close(r[0], 0.34180297, 1.0e-8))
  assert_true(assert_close(r[1], 0.29672718, 1.0e-8))
  assert_true(assert_close(r[2], 0.22054469, 1.0e-8))
  assert_true(assert_close(r[3], 0.14092516, 1.0e-8))
  // round-trip: ilr_inv(ilr(x)) == closure(x)
  let x = [0.1, 0.3, 0.4, 0.2]
  let back = @src.skbio_ilr_inv(@src.skbio_ilr(x, basis), basis)
  assert_true(assert_close(back[0], 0.1, 1.0e-9))
  assert_true(assert_close(back[1], 0.3, 1.0e-9))
  assert_true(assert_close(back[2], 0.4, 1.0e-9))
  assert_true(assert_close(back[3], 0.2, 1.0e-9))
}

///|
test "skbio_centralize" {
  // doctest: centralize([[.1,.3,.4,.2],[.2,.2,.2,.4]])
  //   = [[.17445763,.30216948,.34891526,.17445763],[.32495488,.18761279,.16247744,.32495488]]
  let r = @src.skbio_centralize([[0.1, 0.3, 0.4, 0.2], [0.2, 0.2, 0.2, 0.4]])
  assert_true(assert_close(r[0][0], 0.17445763, 1.0e-8))
  assert_true(assert_close(r[0][1], 0.30216948, 1.0e-8))
  assert_true(assert_close(r[0][2], 0.34891526, 1.0e-8))
  assert_true(assert_close(r[0][3], 0.17445763, 1.0e-8))
  assert_true(assert_close(r[1][0], 0.32495488, 1.0e-8))
  assert_true(assert_close(r[1][1], 0.18761279, 1.0e-8))
  assert_true(assert_close(r[1][2], 0.16247744, 1.0e-8))
  assert_true(assert_close(r[1][3], 0.32495488, 1.0e-8))
}