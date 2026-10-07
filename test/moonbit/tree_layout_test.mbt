///|
/// Tests for tree layout (rectangular and circular).

///|
fn tl_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn tl_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("i"), clades=children, branch_length=Some(1.0))
}

test "layout: rectangular basic 2 leaves" {
  let tree = tl_internal([tl_leaf("A"), tl_leaf("B")])
  let layout = @src.phylo_layout_rectangular(tree)
  // 3 nodes: root + 2 leaves
  assert_eq(layout.nodes.length(), 3)
  assert_eq(layout.n_leaves, 2)
  // Leaves should have y = 0 and y = 1
  let mut y0 = -1.0
  let mut y1 = -1.0
  for n in layout.nodes {
    if n.is_leaf && n.name == "A" {
      y0 = n.y
    }
    if n.is_leaf && n.name == "B" {
      y1 = n.y
    }
  }
  assert_eq(y0, 0.0)
  assert_eq(y1, 1.0)
}

///|
test "layout: rectangular x positions" {
  // ((A:1, B:1):1, C:1):1
  // Root x=1, left child x=2, A x=3, B x=3, C x=2
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  // Find root (first node)
  let root = layout.nodes[0]
  assert_true(root.x > 0.0)
  // Width should be max x = 3.0
  assert_true(layout.width > 0.0)
}

///|
test "layout: rectangular internal y is average" {
  // ((A, B), C) — internal node (A,B) should have y = (0+1)/2 = 0.5
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  // Find the internal node that is not root
  for n in layout.nodes {
    if not(n.is_leaf) && n.children.length() == 2 {
      // This should be the (A,B) node with y = 0.5
      let mut has_a = false
      let mut has_b = false
      for ci in n.children {
        if layout.nodes[ci].name == "A" {
          has_a = true
        }
        if layout.nodes[ci].name == "B" {
          has_b = true
        }
      }
      if has_a && has_b {
        assert_true((n.y - 0.5).abs() < 0.001)
      }
    }
  }
}

///|
test "layout: rectangular 5 leaves" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_internal([
      tl_internal([tl_leaf("C"), tl_leaf("D")]),
      tl_leaf("E"),
    ]),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  assert_eq(layout.n_leaves, 5)
  // Leaves should have y = 0, 1, 2, 3, 4
  let mut leaf_count = 0
  for n in layout.nodes {
    if n.is_leaf {
      leaf_count = leaf_count + 1
    }
  }
  assert_eq(leaf_count, 5)
}

///|
test "layout: rectangular single leaf" {
  let tree = tl_leaf("A")
  let layout = @src.phylo_layout_rectangular(tree)
  assert_eq(layout.nodes.length(), 1)
  assert_eq(layout.n_leaves, 1)
}

///|
test "layout: circular basic" {
  let tree = tl_internal([tl_leaf("A"), tl_leaf("B")])
  let layout = @src.phylo_layout_circular(tree)
  assert_eq(layout.nodes.length(), 3)
  assert_eq(layout.n_leaves, 2)
  assert_true(layout.max_radius > 0.0)
}

///|
test "layout: circular angles distributed" {
  // 4 leaves → angles at 0, π/2, π, 3π/2
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_internal([tl_leaf("C"), tl_leaf("D")]),
  ])
  let layout = @src.phylo_layout_circular(tree)
  // Check that leaf angles are distributed
  let angles : Array[Double] = Array(capacity=0)
  for n in layout.nodes {
    if n.is_leaf {
      angles.push(n.angle)
    }
  }
  assert_eq(angles.length(), 4)
  // First angle should be 0
  assert_true(angles[0].abs() < 0.01)
}

///|
test "layout: circular cartesian coords" {
  let tree = tl_internal([tl_leaf("A"), tl_leaf("B")])
  let layout = @src.phylo_layout_circular(tree)
  // Each node should have x, y computed
  for n in layout.nodes {
    // x = radius * cos(angle), y = radius * sin(angle)
    let expected_x = n.radius * tl_cos(n.angle)
    let expected_y = n.radius * tl_sin(n.angle)
    assert_true((n.x - expected_x).abs() < 0.01)
    assert_true((n.y - expected_y).abs() < 0.01)
  }
}

///|
/// Approximate cos for test verification.
fn tl_cos(x : Double) -> Double {
  let pi = 3.141592653589793
  let two_pi = 2.0 * pi
  let mut y = x
  while y > pi {
    y = y - two_pi
  }
  while y < -pi {
    y = y + two_pi
  }
  let y2 = y * y
  1.0 - y2 / 2.0 + y2 * y2 / 24.0 - y2 * y2 * y2 / 720.0
}

///|
fn tl_sin(x : Double) -> Double {
  let pi = 3.141592653589793
  let two_pi = 2.0 * pi
  let mut y = x
  while y > pi {
    y = y - two_pi
  }
  while y < -pi {
    y = y + two_pi
  }
  let y2 = y * y
  y * (1.0 - y2 / 6.0 + y2 * y2 / 120.0 - y2 * y2 * y2 / 5040.0)
}

///|
test "layout: circular custom range" {
  let tree = tl_internal([
    tl_leaf("A"),
    tl_leaf("B"),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_circular_range(
    tree, 0.0, 3.141592653589793
  )
  // Half circle: angles in [0, π]
  for n in layout.nodes {
    if n.is_leaf {
      assert_true(n.angle >= 0.0 - 0.01)
      assert_true(n.angle <= 3.141592653589793 + 0.01)
    }
  }
}

///|
test "layout: render rectangular non-empty" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  let s = @src.phylo_render_rectangular(layout, 30, 10)
  assert_true(s.length() > 0)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
  assert_true(s.contains("C"))
}

///|
test "layout: render circular non_empty" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_internal([tl_leaf("C"), tl_leaf("D")]),
  ])
  let layout = @src.phylo_layout_circular(tree)
  let s = @src.phylo_render_circular(layout, 30, 30)
  assert_true(s.length() > 0)
}

///|
test "layout: render rectangular empty tree" {
  let layout : @src.TreeLayoutRect = @src.TreeLayoutRect::{
    nodes: [],
    width: 0.0,
    height: 0.0,
    n_leaves: 0,
  }
  let s = @src.phylo_render_rectangular(layout, 10, 10)
  assert_eq(s, "")
}

///|
test "layout: render circular empty tree" {
  let layout : @src.TreeLayoutCirc = @src.TreeLayoutCirc::{
    nodes: [],
    max_radius: 0.0,
    n_leaves: 0,
    start_angle: 0.0,
    end_angle: 6.283185307179586,
  }
  let s = @src.phylo_render_circular(layout, 10, 10)
  assert_eq(s, "")
}

///|
test "layout: rectangular depth tracking" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  // Root depth = 0, children depth = 1, grandchildren depth = 2
  assert_eq(layout.nodes[0].depth, 0)
  for n in layout.nodes {
    if n.is_leaf && (n.name == "A" || n.name == "B") {
      assert_eq(n.depth, 2)
    }
    if n.is_leaf && n.name == "C" {
      assert_eq(n.depth, 1)
    }
  }
}

///|
test "layout: circular depth tracking" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_leaf("C"),
  ])
  let layout = @src.phylo_layout_circular(tree)
  assert_eq(layout.nodes[0].depth, 0)
  for n in layout.nodes {
    if n.is_leaf && (n.name == "A" || n.name == "B") {
      assert_eq(n.depth, 2)
    }
  }
}

///|
test "layout: rectangular branch lengths preserved" {
  let tree = @src.Clade::new(
    name=Some("root"),
    clades=[
      @src.Clade::new(name=Some("A"), branch_length=Some(2.5)),
      @src.Clade::new(name=Some("B"), branch_length=Some(1.5)),
    ],
  )
  let layout = @src.phylo_layout_rectangular(tree)
  for n in layout.nodes {
    if n.name == "A" {
      assert_true((n.branch_length - 2.5).abs() < 0.001)
    }
    if n.name == "B" {
      assert_true((n.branch_length - 1.5).abs() < 0.001)
    }
  }
}

///|
test "layout: circular radius equals cumulative branch" {
  // Root bl=1, children bl=1 → leaves at radius 2.0
  let tree = tl_internal([tl_leaf("A"), tl_leaf("B")])
  let layout = @src.phylo_layout_circular(tree)
  for n in layout.nodes {
    if n.is_leaf {
      assert_true((n.radius - 2.0).abs() < 0.001)
    }
  }
}

///|
test "layout: render rectangular large tree" {
  let tree = tl_internal([
    tl_internal([
      tl_internal([tl_leaf("A"), tl_leaf("B")]),
      tl_internal([tl_leaf("C"), tl_leaf("D")]),
    ]),
    tl_internal([
      tl_leaf("E"),
      tl_internal([tl_leaf("F"), tl_leaf("G")]),
    ]),
  ])
  let layout = @src.phylo_layout_rectangular(tree)
  let s = @src.phylo_render_rectangular(layout, 40, 15)
  assert_true(s.contains("A"))
  assert_true(s.contains("G"))
}

///|
test "layout: circular 6 leaves" {
  let tree = tl_internal([
    tl_internal([tl_leaf("A"), tl_leaf("B")]),
    tl_internal([
      tl_internal([tl_leaf("C"), tl_leaf("D")]),
      tl_internal([tl_leaf("E"), tl_leaf("F")]),
    ]),
  ])
  let layout = @src.phylo_layout_circular(tree)
  assert_eq(layout.n_leaves, 6)
  // All 6 leaf angles should be distinct
  let angles : Array[Double] = Array(capacity=0)
  for n in layout.nodes {
    if n.is_leaf {
      angles.push(n.angle)
    }
  }
  for i in 0..<angles.length() {
    for j in (i + 1)..<angles.length() {
      assert_true((angles[i] - angles[j]).abs() > 0.01)
    }
  }
}
