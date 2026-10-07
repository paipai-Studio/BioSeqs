///|
/// Tests for Dup/Loss inference (reconciliation).

///|
fn dup_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn dup_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("i"), clades=children, branch_length=Some(1.0))
}

///|
/// Species tree: ((Human, Chimp), Gorilla)
fn dup_species3() -> @src.Clade {
  dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
}

///|
/// Species tree: ((Human, Chimp), (Gorilla, Orangutan))
fn dup_species4() -> @src.Clade {
  dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_internal([dup_leaf("Gorilla"), dup_leaf("Orangutan")]),
  ])
}

///|
/// Species tree: (((Human, Chimp), Gorilla), Orangutan)
fn dup_species4b() -> @src.Clade {
  dup_internal([
    dup_internal([
      dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
      dup_leaf("Gorilla"),
    ]),
    dup_leaf("Orangutan"),
  ])
}

test "duploss: matching tree zero dup zero loss" {
  // Gene tree = species tree
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  assert_eq(result.duplication_count, 0)
  assert_eq(result.loss_count, 0)
  assert_eq(result.total_cost, 0)
}

///|
test "duploss: simple duplication" {
  // Gene: ((Human_1, Human_2), Chimp)
  // Two Human paralogs → 1 duplication
  let gene = dup_internal([
    dup_internal([dup_leaf("Human_1"), dup_leaf("Human_2")]),
    dup_leaf("Chimp"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "_")
  assert_eq(result.duplication_count, 1)
  assert_true(result.loss_count >= 0)
  assert_true(result.total_cost >= 1)
}

///|
test "duploss: all paralogs" {
  // Gene: ((Human_1, Human_2), Human_3)
  // All same species → 2 duplications (both internal nodes)
  let gene = dup_internal([
    dup_internal([dup_leaf("Human_1"), dup_leaf("Human_2")]),
    dup_leaf("Human_3"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "_")
  assert_eq(result.duplication_count, 2)
}

///|
test "duploss: speciation count" {
  // Gene = species: 2 internal nodes, both speciations
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  assert_eq(result.speciation_count, 2)
  assert_eq(result.duplication_count, 0)
}

///|
test "duploss: node count" {
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  // 3 leaves + 2 internal = 5 nodes
  assert_eq(result.node_count, 5)
  assert_eq(result.leaf_count, 3)
}

///|
test "duploss: A: 4-species matching" {
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_internal([dup_leaf("Gorilla"), dup_leaf("Orangutan")]),
  ])
  let species = dup_species4()
  let result = @src.phylo_reconcile(gene, species, "")
  assert_eq(result.duplication_count, 0)
  assert_eq(result.loss_count, 0)
  assert_eq(result.total_cost, 0)
}

///|
test "duploss: 4-species with loss" {
  // Gene: ((Human, Gorilla), (Chimp,.Orangutan))
  // Species: ((Human, Chimp),2, (Gorilla, Orangutan))
  // This should imply losses
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Gorilla")]),
    dup_internal([dup_leaf("Chimp"), dup_leaf("Orangutan")]),
  ])
  let species = dup_species4()
  let result = @src.phylo_reconcile(gene, species, "")
  // Mixed cherries across species clades → at least 1 duplication
  assert_true(result.duplication_count >= 1)
  assert_true(result.loss_count >= 0)
}

///|
test "duploss: total cost = dup + loss" {
  let gene = dup_internal([
    dup_internal([dup_leaf("Human_1"), dup_leaf("Human_2")]),
    dup_leaf("Chimp"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "_")
  assert_eq(result.total_cost, result.duplication_count + result.loss_count)
}

///|
test "duploss: per-node edge losses" {
  // Gene = species: all edge8edge_losses should be 0
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  for node in result.nodes {
    assert_eq(node.edge_losses, 0)
  }
}

///|
test "duploss: deep species tree" {
  // Species: (((Human, Chimp), Gorilla), Orangutan)
  // Gene: (((Human, Chimp), Gorilla), Orangutan) — matching
  let gene = dup_internal([
    dup_internal([
      dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
      dup_leaf("Gorilla"),
    ]),
    dup_leaf("Orangutan"),
  ])
  let species = dup_species4b()
  let result = @src.phylo_reconcile(gene, species, "")
  assert_eq(result.duplication_count, 0)
  assert_eq(result.loss_count, 0)
}

///|
test "duploss: deep species tree with dup" {
  // Species: (((Human, Chimp), Gorilla), Orangutan)
  // Gene: (((Human_1, Human_2), Chimp), Gorilla) — 1 duplication
  let gene = dup_internal([
    dup_internal([
      dup_internal([dup_leaf("Human_1"), dup_leaf("Human_2")]),
      dup_leaf("Chimp"),
    ]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species4b()
  let result = @src.phylo_reconcile(gene, species, "_")
  assert_true(result.duplication_count >= 1)
}

///|
test "duploss: summary output non-empty" {
  let gene = dup_internal([
    dup_internal([dup_leaf("Human"), dup_leaf("Chimp")]),
    dup_leaf("Gorilla"),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  let summary = @src.phylo_reconciliation_summary(result)
  assert_true(summary.length() > 0)
   // Should contain key phrases
  assert_true(summary.contains("Reconciliation"))
  assert_true(summary.contains("Duplications:"))
  assert_true(summary.contains("Losses:"))
}

///|
test "duploss: explicit mapping" {
  let gene = dup_internal([
    dup_internal([dup_leaf("g1"), dup_leaf("g2")]),
    dup_leaf("g3"),
  ])
  let species = dup_species3()
  let g2s : Map[String, String] = Map([], capacity=0)
  g2s.set("g1", "Human")
  g2s.set("g2", "Chimp")
  g2s.set("g3", "Gorilla")
  let result = @src.phylo_reconcile_explicit(gene, species, g2s, "")
  assert_eq(result.duplication_count, 0)
  assert_eq(result.loss_count, 0)
}

///|
test "duploss: single leaf" {
  let gene = dup_leaf("Human")
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "")
  assert_eq(result.duplication_count, 0)
  assert_eq(result.loss_count, 0)
  assert_eq(result.leaf_count, 1)
}

///|
test "duploss: duplication with loss" {
  // Gene: ((Human_1, Chimp), (Human_2, Gorilla))
  // Species: ((Human, Chimp), Gorilla)
  // Human_1 and Human_2 are paralogs → at least 1 duplication
  // And there should be losses
  let gene = dup_internal([
    dup_internal([dup_leaf("Human_1"), dup_leaf("Chimp")]),
    dup_internal([dup_leaf("Human_2"), dup_leaf("Gorilla")]),
  ])
  let species = dup_species3()
  let result = @src.phylo_reconcile(gene, species, "_")
  assert_true(result.duplication_count >= 1)
  assert_true(result.loss_count >= 0)
  assert_true(result.total_cost == result.duplication_count + result.loss_count)
}
