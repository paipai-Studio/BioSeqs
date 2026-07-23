///|
/// Tests for topGO module.

test "TopGOTerm creation" {
  let term = @src.TopGOTerm::new("GO:0008150", "biological_process", "BP", 0)
  assert_eq(term.go_id, "GO:0008150")
  assert_eq(term.name, "biological_process")
  assert_eq(term.namespace, "BP")
  assert_eq(term.level, 0)
  assert_eq(term.gene_count, 0)
}

test "TopGOGraph construction" {
  let mut graph = @src.TopGOGraph::new("GO:0008150")
  let root = @src.TopGOTerm::new("GO:0008150", "biological_process", "BP", 0)
  graph = graph.add_term(root)
  assert_true(graph.terms.contains("GO:0008150"))
}

test "TopGOGraph add edge" {
  let mut graph = @src.TopGOGraph::new("GO:0008150")
  let root = @src.TopGOTerm::new("GO:0008150", "biological_process", "BP", 0)
  let child = @src.TopGOTerm::new("GO:0009987", "cellular_process", "BP", 1)
  graph = graph.add_term(root)
  graph = graph.add_term(child)
  graph = graph.add_edge("GO:0008150", "GO:0009987")
  assert_true(graph.terms.contains("GO:0008150"))
  assert_true(graph.terms.contains("GO:0009987"))
}

test "topgo_elim_algorithm" {
  let graph = @src.create_example_topgo_graph()
  let genes_of_interest = ["gene1", "gene2", "gene3"]
  let all_genes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"]
  let results = @src.topgo_elim_algorithm(graph, genes_of_interest, all_genes, "BP")
  assert_true(results.length() > 0)
}

test "topgo_weight01_algorithm" {
  let graph = @src.create_example_topgo_graph()
  let genes_of_interest = ["gene1", "gene2", "gene3"]
  let all_genes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"]
  let results = @src.topgo_weight01_algorithm(graph, genes_of_interest, all_genes, "BP")
  assert_true(results.length() > 0)
}

test "topgo_fisher_exact" {
  let p = @src.topgo_fisher_exact(3, 10, 5, 100)
  assert_true(p >= 0.0 && p <= 1.0)
}
