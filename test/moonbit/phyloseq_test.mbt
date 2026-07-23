///|
/// Tests for phyloseq module.

test "phyloseq creation" {
  let ps = @src.create_example_phyloseq()
  assert_eq(@src.ps_num_otus(ps), 3)
  assert_eq(@src.ps_num_samples(ps), 3)
}

test "phyloseq total abundance" {
  let ps = @src.create_example_phyloseq()
  let total = @src.ps_total_abundance(ps)
  assert_eq(total, 975.0)
}

test "phyloseq filter by abundance" {
  let ps = @src.create_example_phyloseq()
  let filtered = @src.ps_filter_by_abundance(ps, 300.0)
  assert_eq(@src.ps_num_otus(filtered), 2)
}

test "phyloseq filter by taxonomy" {
  let ps = @src.create_example_phyloseq()
  let filtered = @src.ps_filter_by_taxonomy(ps, "phylum", "Proteobacteria")
  assert_eq(@src.ps_num_otus(filtered), 1)
}

test "phyloseq taxa summary" {
  let ps = @src.create_example_phyloseq()
  let summary = @src.ps_get_taxa_summary(ps, "phylum")
  assert_true(summary.contains("Proteobacteria"))
}