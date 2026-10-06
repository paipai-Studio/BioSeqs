///|
/// Tests for SPR distance (minimum SPR moves between two trees).

///|
fn sprd_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn sprd_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children, branch_length=Some(1.0))
}

///|
/// ((A,B),(C,D))
fn sprd_balanced4() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("B")]),
    sprd_internal([sprd_leaf("C"), sprd_leaf("D")]),
  ])
}

///|
/// ((A,C),(B,D))
fn sprd_cross4() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("C")]),
    sprd_internal([sprd_leaf("B"), sprd_leaf("D")]),
  ])
}

///|
/// ((A,D),(B,C))
fn sprd_cross4b() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("D")]),
    sprd_internal([sprd_leaf("B"), sprd_leaf("C")]),
  ])
}

///|
/// ((A,B),(C,(D,E)))
fn sprd_5a() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("B")]),
    sprd_internal([
      sprd_leaf("C"),
      sprd_internal([sprd_leaf("D"), sprd_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,C),(B,(D,E)))
fn sprd_5b() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("C")]),
    sprd_internal([
      sprd_leaf("B"),
      sprd_internal([sprd_leaf("D"), sprd_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,D),(B,(C,E)))
fn sprd_5c() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("D")]),
    sprd_internal([
      sprd_leaf("B"),
      sprd_internal([sprd_leaf("C"), sprd_leaf("E")]),
    ]),
  ])
}

///|
/// ((A,B),((C,D),(E,F)))
fn sprd_6a() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("B")]),
    sprd_internal([
      sprd_internal([sprd_leaf("C"), sprd_leaf("D")]),
      sprd_internal([sprd_leaf("E"), sprd_leaf("F")]),
    ]),
  ])
}

///|
/// ((A,C),((B,D),(E,F)))
fn sprd_6b() -> @src.Clade {
  sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("C")]),
    sprd_internal([
      sprd_internal([sprd_leaf("B"), sprd_leaf("D")]),
      sprd_internal([sprd_leaf("E"), sprd_leaf("F")]),
    ]),
  ])
}

///|
/// ((A,B),(C,D)) — same topology, different branch lengths
fn sprd_balanced4_diffbl() -> @src.Clade {
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

test "sprd: identical trees distance 0" {
  let t = sprd_balanced4()
  assert_eq(@src.dendropy_spr_distance(t, t), 0)
}

///|
test "sprd: same topology different branch lengths distance 0" {
  let t1 = sprd_balanced4()
  let t2 = sprd_balanced4_diffbl()
  assert_eq(@src.dendropy_spr_distance(t1, t2), 0)
}

///|
test "sprd: different taxon count returns -1" {
  let t1 = sprd_balanced4()
  let t2 = sprd_5a()
  assert_eq(@src.dendropy_spr_distance(t1, t2), -1)
}

///|
test "sprd: 4-leaf different topology distance 1" {
  // ((A,B),(C,D)) vs ((A,C),(B,D)) — one SPR move
  let t1 = sprd_balanced4()
  let t2 = sprd_cross4()
  let d = @src.dendropy_spr_distance(t1, t2)
  assert_true(d >= 1)
  assert_true(d <= 2)
}

///|
test "sprd: 4-leaf all three topologies" {
  let t1 = sprd_balanced4()
  let t2 = sprd_cross4()
  let t3 = sprd_cross4b()
  // All three 4-leaf topologies are within SPR distance 1 of each other
  let d12 = @src.dendropy_spr_distance(t1, t2)
  let d13 = @src.dendropy_spr_distance(t1, t3)
  let d23 = @src.dendropy_spr_distance(t2, t3)
  assert_true(d12 >= 1 && d12 <= 2)
  assert_true(d13 >= 1 && d13 <= 2)
  assert_true(d23 >= 1 && d23 <= 2)
}

///|
test "sprd: 5-leaf same topology distance 0" {
  let t = sprd_5a()
  assert_eq(@src.dendropy_spr_distance(t, t), 0)
}

///|
test "sprd: 5-leaf one SPR move" {
  // ((A,B),(C,(D,E))) vs ((A,C),(B,(D,E))) — should be 1 SPR move
  let t1 = sprd_5a()
  let t2 = sprd_5b()
  let d = @src.dendropy_spr_distance(t1, t2)
  assert_true(d >= 1)
  assert_true(d <= 2)
}

///|
test "sprd: 5-leaf distance symmetric" {
  let t1 = sprd_5a()
  let t2 = sprd_5b()
  let t3 = sprd_5c()
  let d12 = @src.dendropy_spr_distance(t1, t2)
  let d21 = @src.dendropy_spr_distance(t2, t1)
  assert_eq(d12, d21)
  let d13 = @src.dendropy_spr_distance(t1, t3)
  let d31 = @src.dendropy_spr_distance(t3, t1)
  assert_eq(d13, d31)
}

///|
test "sprd: 3-leaf trees always distance 0" {
  // Only one unrooted topology for 3 taxa
  let t1 = sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("B")]),
    sprd_leaf("C"),
  ])
  let t2 = sprd_internal([
    sprd_internal([sprd_leaf("A"), sprd_leaf("C")]),
    sprd_leaf("B"),
  ])
  assert_eq(@src.dendropy_spr_distance(t1, t2), 0)
}

///|
test "sprd: 6-leaf same topology distance 0" {
  let t = sprd_6a()
  assert_eq(@src.dendropy_spr_distance(t, t), 0)
}

///|
test "sprd: 6-leaf one SPR move" {
  // ((A,B),((C,D),(E,F))) vs ((A,C),((B,D),(E,F)))
  let t1 = sprd_6a()
  let t2 = sprd_6b()
  let d = @src.dendropy_spr_distance(t1, t2)
  assert_true(d >= 1)
  assert_true(d <= 2)
}

///|
test "sprd: triangle inequality" {
  let t1 = sprd_5a()
  let t2 = sprd_5b()
  let t3 = sprd_5c()
  let d12 = @src.dendropy_spr_distance(t1, t2)
  let d23 = @src.dendropy_spr_distance(t2, t3)
  let d13 = @src.dendropy_spr_distance(t1, t3)
  assert_true(d13 <= d12 + d23)
}

///|
test "sprd: non-negativity" {
  let t1 = sprd_balanced4()
  let t2 = sprd_cross4()
  let d = @src.dendropy_spr_distance(t1, t2)
  assert_true(d >= 0)
}

///|
test "sprd: 4-leaf distance to self is 0" {
  let t1 = sprd_cross4()
  let t2 = sprd_cross4b()
  // These are different topologies
  let d = @src.dendropy_spr_distance(t1, t2)
  assert_true(d >= 1)
}
