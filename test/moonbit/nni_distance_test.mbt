///|
/// Tests for NNI distance (minimum NNI moves between two trees).

///|
fn nnid_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn nnid_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children, branch_length=Some(1.0))
}

///|
/// ((A,B),(C,D))
fn nnid_balanced4() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("B")]),
    nnid_internal([nnid_leaf("C"), nnid_leaf("D")]),
  ])
}

///|
/// ((A,C),(B,D))
fn nnid_cross4() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("C")]),
    nnid_internal([nnid_leaf("B"), nnid_leaf("D")]),
  ])
}

///|
/// ((A,D),(B,C))
fn nnid_cross4b() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("D")]),
    nnid_internal([nnid_leaf("B"), nnid_leaf("C")]),
  ])
}

///|
/// ((A,B),(C,(D,E)))
fn nnid_5a() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("B")]),
    nnid_internal([
      nnid_leaf("C"),
      nnid_internal([nnid_leaf("D"), nnid_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,C),(B,(D,E)))
fn nnid_5b() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("C")]),
    nnid_internal([
      nnid_leaf("B"),
      nnid_internal([nnid_leaf("D"), nnid_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,D),(B,(C,E)))
fn nnid_5c() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("D")]),
    nnid_internal([
      nnid_leaf("B"),
      nnid_internal([nnid_leaf("C"), nnid_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,B),((C,D),(E,F)))
fn nnid_6a() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("B")]),
    nnid_internal([
      nnid_internal([nnid_leaf("C"), nnid_leaf("D")]),
      nnid_internal([nnid_leaf("E"), nnid_leaf("F")]),
    ]),
  ])
}

///|
/// ((A,C),((B,D),(E,F)))
fn nnid_6b() -> @src.Clade {
  nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("C")]),
    nnid_internal([
      nnid_internal([nnid_leaf("B"), nnid_leaf("D")]),
      nnid_internal([nnid_leaf("E"), nnid_leaf("F")]),
    ]),
  ])
}

///|
/// ((A,B),(C,D)) — same topology, different branch lengths
fn nnid_balanced4_diffbl() -> @src.Clade {
  @src.Clade::new(
    name=Some("internal"),
    clades=[
      @src.Clade::new(
        name=Some("internal"),
        clades=[
          @src.Clade::new(name=Some("A"), branch_length=Some(3.0)),
          @src.Clade::new(name=Some("B"), branch_length=Some(2.0)),
        ],
        branch_length=Some(5.0),
      ),
      @src.Clade::new(
        name=Some("internal"),
        clades=[
          @src.Clade::new(name=Some("C"), branch_length=Some(0.5)),
          @src.Clade::new(name=Some("D"), branch_length=Some(1.5)),
        ],
        branch_length=Some(4.0),
      ),
    ],
    branch_length=Some(0.0),
  )
}

///|
test "nnid: identical trees distance 0" {
  let t = nnid_balanced4()
  assert_eq(@src.dendropy_nni_distance(t, t), 0)
}

///|
test "nnid: same topology different branch lengths distance 0" {
  let t1 = nnid_balanced4()
  let t2 = nnid_balanced4_diffbl()
  assert_eq(@src.dendropy_nni_distance(t1, t2), 0)
}

///|
test "nnid: different taxon count returns -1" {
  let t1 = nnid_balanced4()
  let t2 = nnid_5a()
  assert_eq(@src.dendropy_nni_distance(t1, t2), -1)
}

///|
test "nnid: 4-leaf different topology distance 1" {
  // ((A,B),(C,D)) vs ((A,C),(B,D)) — one NNI move
  let t1 = nnid_balanced4()
  let t2 = nnid_cross4()
  let d = @src.dendropy_nni_distance(t1, t2)
  assert_eq(d, 1)
}

///|
test "nnid: 4-leaf all three topologies distance 1" {
  let t1 = nnid_balanced4()
  let t2 = nnid_cross4()
  let t3 = nnid_cross4b()
  // All three 4-leaf topologies are one NNI move apart
  assert_eq(@src.dendropy_nni_distance(t1, t2), 1)
  assert_eq(@src.dendropy_nni_distance(t1, t3), 1)
  assert_eq(@src.dendropy_nni_distance(t2, t3), 1)
}

///|
test "nnid: 5-leaf same topology distance 0" {
  let t = nnid_5a()
  assert_eq(@src.dendropy_nni_distance(t, t), 0)
}

///|
test "nnid: 5-leaf one NNI move" {
  // ((A,B),(C,(D,E))) vs ((A,C),(B,(D,E))) — one NNI move
  let t1 = nnid_5a()
  let t2 = nnid_5b()
  let d = @src.dendropy_nni_distance(t1, t2)
  assert_eq(d, 1)
}

///|
test "nnid: 5-leaf distance symmetric" {
  let t1 = nnid_5a()
  let t2 = nnid_5b()
  let t3 = nnid_5c()
  let d12 = @src.dendropy_nni_distance(t1, t2)
  let d21 = @src.dendropy_nni_distance(t2, t1)
  assert_eq(d12, d21)
  let d13 = @src.dendropy_nni_distance(t1, t3)
  let d31 = @src.dendropy_nni_distance(t3, t1)
  assert_eq(d13, d31)
}

///|
test "nnid: 3-leaf trees always distance 0" {
  let t1 = nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("B")]),
    nnid_leaf("C"),
  ])
  let t2 = nnid_internal([
    nnid_internal([nnid_leaf("A"), nnid_leaf("C")]),
    nnid_leaf("B"),
  ])
  assert_eq(@src.dendropy_nni_distance(t1, t2), 0)
}

///|
test "nnid: 6-leaf same topology distance 0" {
  let t = nnid_6a()
  assert_eq(@src.dendropy_nni_distance(t, t), 0)
}

///|
test "nnid: 6-leaf NNI distance" {
  // ((A,B),((C,D),(E,F))) vs ((A,C),((B,D),(E,F)))
  // Not 1 NNI move: must swap A<->C across multiple edges
  let t1 = nnid_6a()
  let t2 = nnid_6b()
  let d = @src.dendropy_nni_distance(t1, t2)
  assert_true(d >= 2)
  assert_true(d <= 4)
}

///|
test "nnid: triangle inequality" {
  let t1 = nnid_5a()
  let t2 = nnid_5b()
  let t3 = nnid_5c()
  let d12 = @src.dendropy_nni_distance(t1, t2)
  let d23 = @src.dendropy_nni_distance(t2, t3)
  let d13 = @src.dendropy_nni_distance(t1, t3)
  assert_true(d13 <= d12 + d23)
}

///|
test "nnid: non-negativity" {
  let t1 = nnid_balanced4()
  let t2 = nnid_cross4()
  let d = @src.dendropy_nni_distance(t1, t2)
  assert_true(d >= 0)
}

///|
test "nnid: NNI distance >= SPR distance" {
  // NNI is a subset of SPR, so NNI distance >= SPR distance
  let t1 = nnid_5a()
  let t2 = nnid_5c()
  let d_nni = @src.dendropy_nni_distance(t1, t2)
  let d_spr = @src.dendropy_spr_distance(t1, t2)
  assert_true(d_nni >= d_spr)
}
