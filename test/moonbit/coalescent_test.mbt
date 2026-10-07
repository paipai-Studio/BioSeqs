///|
/// Tests for Kingman coalescent simulation and tree visualization.

///|
fn coal_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn coal_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("i"), clades=children, branch_length=Some(1.0))
}

test "coalescent: basic simulation n=2" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(2, rng)
  assert_eq(result.n_lineages, 2)
  assert_eq(result.coalescent_times.length(), 1)
  assert_true(result.total_depth > 0.0)
  assert_eq(result.tree.count_terminals(), 2)
}

///|
test "coalescent: basic simulation n=5" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(5, rng)
  assert_eq(result.n_lineages, 5)
  assert_eq(result.coalescent_times.length(), 4)
  assert_true(result.total_depth > 0.0)
  assert_eq(result.tree.count_terminals(), 5)
}

///|
test "coalescent: n=1 returns single leaf" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(1, rng)
  assert_eq(result.n_lineages, 1)
  assert_eq(result.coalescent_times.length(), 0)
  assert_eq(result.total_depth, 0.0)
}

///|
test "coalescent: n=0 returns empty" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(0, rng)
  assert_eq(result.n_lineages, 0)
  assert_eq(result.total_depth, 0.0)
}

///|
test "coalescent: times are monotonically increasing" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(10, rng)
  for i in 1..<result.coalescent_times.length() {
    assert_true(result.coalescent_times[i] > result.coalescent_times[i - 1])
  }
}

///|
test "coalescent: branch lengths are positive" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(5, rng)
  assert_true(coal_check_positive_branches(result.tree))
}

///|
fn coal_check_positive_branches(node : @src.Clade) -> Bool {
  match node.branch_length {
    Some(bl) => if bl <= 0.0 { return false }
    None => ()
  }
  for ch in node.clades {
    if not(coal_check_positive_branches(ch)) {
      return false
    }
  }
  true
}

///|
test "coalescent: with Ne scaling" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent_ne(5, 100.0, rng)
  assert_true(result.total_depth > 0.0)
  assert_eq(result.n_lineages, 5)
}

///|
test "coalescent: larger Ne gives larger depth" {
  let rng1 = @src.Rng::new(42UL)
  let rng2 = @src.Rng::new(42UL)
  let r1 = @src.phylo_simulate_coalescent_ne(5, 1.0, rng1)
  let r2 = @src.phylo_simulate_coalescent_ne(5, 10.0, rng2)
  assert_true(r2.total_depth > r1.total_depth)
}

///|
test "coalescent: reproducible with same seed" {
  let rng1 = @src.Rng::new(123UL)
  let rng2 = @src.Rng::new(123UL)
  let r1 = @src.phylo_simulate_coalescent(5, rng1)
  let r2 = @src.phylo_simulate_coalescent(5, rng2)
  assert_eq(r1.total_depth, r2.total_depth)
  assert_eq(r1.coalescent_times.length(), r2.coalescent_times.length())
}

///|
test "coalescent: expected depth formula" {
  let exp2 = @src.phylo_coalescent_expected_depth(2)
  assert_true((exp2 - 1.0).abs() < 0.001)
  let exp5 = @src.phylo_coalescent_expected_depth(5)
  assert_true((exp5 - 1.6).abs() < 0.001)
  let exp10 = @src.phylo_coalescent_expected_depth(10)
  assert_true((exp10 - 1.8).abs() < 0.001)
}

///|
test "coalescent: mean depth approximates expected" {
  let rng = @src.Rng::new(42UL)
  let mean_depth = @src.phylo_coalescent_mean_depth(10, 1.0, 1000, rng)
  let expected = @src.phylo_coalescent_expected_depth(10)
  let ratio = mean_depth / expected
  assert_true(ratio > 0.7 && ratio < 1.3)
}

///|
test "coalescent: watterson expected" {
  let s2 = @src.phylo_watterson_expected(2, 1.0)
  assert_true((s2 - 1.0).abs() < 0.001)
  let s3 = @src.phylo_watterson_expected(3, 1.0)
  assert_true((s3 - 1.5).abs() < 0.001)
}

///|
test "coalescent: summary output" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(5, rng)
  let summary = @src.phylo_coalescent_summary(result)
  assert_true(summary.contains("Kingman"))
  assert_true(summary.contains("Lineages: 5"))
  assert_true(summary.contains("Expected depth"))
}

///|
test "coalescent: tree has correct number of internal nodes" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(7, rng)
  assert_eq(result.coalescent_times.length(), 6)
}

// =========================================================================
// Visualization tests
// =========================================================================

test "viz: draw ascii enhanced basic" {
  let tree = coal_internal([coal_leaf("A"), coal_leaf("B")])
  let s = @src.phylo_draw_ascii_enhanced(tree)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
  assert_true(s.contains("i"))
}

///|
test "viz: draw ascii enhanced with branch lengths" {
  let tree = @src.Clade::new(
    name=Some("root"),
    clades=[
      @src.Clade::new(name=Some("A"), branch_length=Some(1.5)),
      @src.Clade::new(name=Some("B"), branch_length=Some(2.3)),
    ],
  )
  let s = @src.phylo_draw_ascii_enhanced(tree)
  assert_true(s.contains("1.5"))
  assert_true(s.contains("2.3"))
}

///|
test "viz: draw ascii enhanced with confidence" {
  let tree = @src.Clade::new(
    name=Some("node"),
    confidence=Some(95.0),
    clades=[
      @src.Clade::new(name=Some("A"), branch_length=Some(1.0)),
      @src.Clade::new(name=Some("B"), branch_length=Some(1.0)),
    ],
  )
  let s = @src.phylo_draw_ascii_enhanced(tree)
  assert_true(s.contains("95"))
}

///|
test "viz: draw cladogram basic" {
  let tree = coal_internal([
    coal_internal([coal_leaf("A"), coal_leaf("B")]),
    coal_leaf("C"),
  ])
  let s = @src.phylo_draw_cladogram(tree)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
  assert_true(s.contains("C"))
}

///|
test "viz: draw rectangular basic" {
  let tree = coal_internal([coal_leaf("A"), coal_leaf("B")])
  let s = @src.phylo_draw_rectangular(tree)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
}

///|
test "viz: draw barplot basic" {
  let tree = coal_internal([
    @src.Clade::new(name=Some("A"), branch_length=Some(1.0)),
    @src.Clade::new(name=Some("B"), branch_length=Some(2.0)),
  ])
  let s = @src.phylo_draw_barplot(tree, 20)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
}

///|
test "viz: draw barplot proportional" {
  let tree = coal_internal([
    @src.Clade::new(name=Some("short"), branch_length=Some(0.5)),
    @src.Clade::new(name=Some("long"), branch_length=Some(5.0)),
  ])
  let s = @src.phylo_draw_barplot(tree, 20)
  let short_count = coal_count_blocks(s, "short")
  let long_count = coal_count_blocks(s, "long")
  assert_true(long_count > short_count)
}

///|
fn coal_count_blocks(s : String, prefix : String) -> Int {
  // Find the line starting with prefix, count block chars on that line
  let mut i = 0
  let n = s.length()
  while i < n {
    // Check if this line starts with prefix
    let plen = prefix.length()
    if i + plen <= n && s.substring(start=i, end=i + plen) == prefix {
      // Count block chars until newline
      let mut count = 0
      let mut j = i
      while j < n {
        let c = s.substring(start=j, end=j + 1)
        if c == "\n" {
          break
        }
        if c == "█" {
          count = count + 1
        }
        j = j + 1
      }
      return count
    }
    // Skip to next line
    while i < n && s.substring(start=i, end=i + 1) != "\n" {
      i = i + 1
    }
    i = i + 1
  }
  0
}

///|
test "viz: draw ascii enhanced single leaf" {
  let tree = coal_leaf("A")
  let s = @src.phylo_draw_ascii_enhanced(tree)
  assert_true(s.contains("A"))
}

///|
test "viz: draw cladogram single leaf" {
  let tree = coal_leaf("A")
  let s = @src.phylo_draw_cladogram(tree)
  assert_true(s.contains("A"))
}

///|
test "viz: coalescent tree can be drawn" {
  let rng = @src.Rng::new(42UL)
  let result = @src.phylo_simulate_coalescent(5, rng)
  let s = @src.phylo_draw_ascii_enhanced(result.tree)
  assert_true(s.length() > 0)
  let s2 = @src.phylo_draw_cladogram(result.tree)
  assert_true(s2.length() > 0)
}

///|
test "viz: deep tree visualization" {
  let tree = coal_internal([
    coal_internal([
      coal_internal([coal_leaf("A"), coal_leaf("B")]),
      coal_leaf("C"),
    ]),
    coal_internal([coal_leaf("D"), coal_leaf("E")]),
  ])
  let s = @src.phylo_draw_ascii_enhanced(tree)
  assert_true(s.contains("A"))
  assert_true(s.contains("B"))
  assert_true(s.contains("C"))
  assert_true(s.contains("D"))
  assert_true(s.contains("E"))
}
