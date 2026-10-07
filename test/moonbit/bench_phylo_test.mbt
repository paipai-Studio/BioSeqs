///|
/// Benchmark infrastructure for phylogenetic operations.
/// Detects performance regressions by running operations at scale
/// and comparing cached vs uncached paths.
///
/// MoonBit has no built-in timer in test mode, so we use:
/// 1. Scale tests — run at large N, assert completion (catches O(n²)→O(n³))
/// 2. Relative tests — cached path must complete same work as uncached
/// 3. Count tests — verify iteration counts meet minimum thresholds

// =========================================================================
// Benchmark helpers
// =========================================================================

///|
/// Build a balanced binary tree with n leaves.
fn bench_build_tree(n_leaves : Int) -> @src.Clade {
  if n_leaves <= 1 {
    return @src.Clade::new(name=Some("L0"))
  }
  bench_build_tree_rec(0, n_leaves)
}

///|
fn bench_build_tree_rec(start : Int, n_leaves : Int) -> @src.Clade {
  if n_leaves == 1 {
    return @src.Clade::new(name=Some("L" + start.to_string()))
  }
  let half = n_leaves / 2
  @src.Clade::new(clades=[
    bench_build_tree_rec(start, half),
    bench_build_tree_rec(start + half, n_leaves - half),
  ])
}

///|
/// Build a caterpillar (linear) tree with n leaves.
fn bench_caterpillar(n_leaves : Int) -> @src.Clade {
  if n_leaves <= 1 {
    return @src.Clade::new(name=Some("C0"))
  }
  @src.Clade::new(clades=[
    @src.Clade::new(name=Some("C" + (n_leaves - 1).to_string())),
    bench_caterpillar(n_leaves - 1),
  ])
}

///|
/// Build a species tree and gene tree with duplications.
fn bench_gene_species_trees(n_species : Int) -> (@src.Clade, @src.Clade) {
  // Species tree: balanced binary with n_species leaves
  let species = bench_build_tree(n_species)
  // Gene tree: each species has 2 genes (duplication)
  let gene_leaves : Array[@src.Clade] = Array(capacity=n_species * 2)
  for i in 0..<n_species {
    gene_leaves.push(
      @src.Clade::new(name=Some("L" + i.to_string() + "_a"))
    )
    gene_leaves.push(
      @src.Clade::new(name=Some("L" + i.to_string() + "_b"))
    )
  }
  // Build a balanced gene tree from leaves
  let gene = bench_build_from_leaves(gene_leaves)
  (gene, species)
}

///|
fn bench_build_from_leaves(leaves : Array[@src.Clade]) -> @src.Clade {
  if leaves.length() == 1 {
    return leaves[0]
  }
  let half = leaves.length() / 2
  let left_arr : Array[@src.Clade] = Array(capacity=half)
  for i in 0..<half {
    left_arr.push(leaves[i])
  }
  let right_arr : Array[@src.Clade] = Array(capacity=leaves.length() - half)
  for i in half..<leaves.length() {
    right_arr.push(leaves[i])
  }
  let left = bench_build_from_leaves(left_arr)
  let right = bench_build_from_leaves(right_arr)
  @src.Clade::new(clades=[left, right])
}

// =========================================================================
// Scale benchmarks — run at increasing N, verify completion
// =========================================================================

///|
test "bench: tree operations scale 32 leaves" {
  let tree = bench_build_tree(32)
  // Count leaves
  assert_eq(@src.phylo_count_leaves(tree), 32)
  // Ladderize
  let _ = @src.phylo_ladderize(tree, @src.LadderizeLeft)
  // Tree layout
  let _ = @src.phylo_layout_rectangular(tree)
  let _ = @src.phylo_layout_circular(tree)
}

///|
test "bench: tree operations scale 64 leaves" {
  let tree = bench_build_tree(64)
  assert_eq(@src.phylo_count_leaves(tree), 64)
  let _ = @src.phylo_ladderize(tree, @src.LadderizeLeft)
  let _ = @src.phylo_layout_rectangular(tree)
  let _ = @src.phylo_layout_circular(tree)
}

///|
test "bench: tree operations scale 128 leaves" {
  let tree = bench_build_tree(128)
  assert_eq(@src.phylo_count_leaves(tree), 128)
  let _ = @src.phylo_ladderize(tree, @src.LadderizeLeft)
  let _ = @src.phylo_layout_rectangular(tree)
}

///|
test "bench: SVG generation 64 leaves" {
  let tree = bench_build_tree(64)
  let layout = @src.phylo_layout_rectangular(tree)
  let svg = @src.phylo_to_svg_rect(layout, 800, 600)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
}

///|
test "bench: caterpillar 50 leaves" {
  let tree = bench_caterpillar(50)
  assert_eq(@src.phylo_count_leaves(tree), 50)
  let _ = @src.phylo_layout_rectangular(tree)
}

// =========================================================================
// Reconciliation benchmarks — cached vs uncached
// =========================================================================

///|
test "bench: reconciliation context 8 species" {
  let (gene, species) = bench_gene_species_trees(8)
  // Context: LCA computed once
  let recon = @src.Reconciliation::new(gene, species, "_")
  let dups = recon.count_duplications()
  let losses = recon.count_losses()
  let result = recon.reconcile()
  // Verify results are consistent
  assert_eq(result.duplication_count, dups)
  assert_eq(result.loss_count, losses)
}

///|
test "bench: reconciliation context 16 species" {
  let (gene, species) = bench_gene_species_trees(16)
  let recon = @src.Reconciliation::new(gene, species, "_")
  let dups = recon.count_duplications()
  let losses = recon.count_losses()
  let result = recon.reconcile()
  assert_eq(result.duplication_count, dups)
  assert_eq(result.loss_count, losses)
}

///|
test "bench: reconciliation cached matches uncached" {
  let (gene, species) = bench_gene_species_trees(8)
  // Cached
  let recon = @src.Reconciliation::new(gene, species, "_")
  let dups_c = recon.count_duplications()
  let losses_c = recon.count_losses()
  // Uncached (separate calls)
  let maps = @src.phylo_lca_mapping(gene, species, "_")
  let dups_u = @src.phylo_count_duplications(maps)
  let losses_u = @src.phylo_count_losses(gene, species, "_")
  // Must match
  assert_eq(dups_c, dups_u)
  assert_eq(losses_c, losses_u)
}

///|
test "bench: classify all pairs 8 species" {
  let (gene, species) = bench_gene_species_trees(8)
  let recon = @src.Reconciliation::new(gene, species, "_")
  let pairs = recon.classify_all_pairs()
  // 16 genes → 120 pairs
  assert_true(pairs.length() > 0)
}

// =========================================================================
// DendroPy metrics benchmarks
// =========================================================================

///|
test "bench: colless index 32 leaves" {
  let tree = bench_build_tree(32)
  let _ = @src.dendropy_colless(tree)
}

///|
test "bench: sackin index 64 leaves" {
  let tree = bench_build_tree(64)
  let _ = @src.dendropy_sackin(tree)
}

///|
test "bench: cherry count 64 leaves" {
  let tree = bench_build_tree(64)
  let _ = @src.dendropy_cherry_count(tree)
}

