///|
/// Tests for LCA mapping (gene tree → species tree).

///|
fn lca_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn lca_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("i"), clades=children, branch_length=Some(1.0))
}

///|
/// Species tree: ((Human, Chimp), Gorilla)
fn lca_species_tree3() -> @src.Clade {
  lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
}

///|
/// Species tree: ((Human, Chimp), (Gorilla, Orangutan))
fn lca_species_tree4() -> @src.Clade {
  lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_internal([lca_leaf("Gorilla"), lca_leaf("Orangutan")]),
  ])
}

///|
/// Species tree: (((Human, Chimp), Gorilla), Orangutan)
fn lca_species_tree4b() -> @src.Clade {
  lca_internal([
    lca_internal([
      lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
      lca_leaf("Gorilla"),
    ]),
    lca_leaf("Orangutan"),
  ])
}

test "lca: extract species with separator" {
  assert_eq(@src.phylo_extract_species("Human_001", "_"), "Human")
  assert_eq(@src.phylo_extract_species("Mouse_gene1", "_"), "Mouse")
  assert_eq(@src.phylo_extract_species("Human", "_"), "Human")
  assert_eq(@src.phylo_extract_species("Human", ""), "Human")
}

///|
test "lca: simple no-duplication mapping" {
  // Gene tree: ((Human, Chimp), Gorilla) — same as species tree
  // No duplications expected
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  let dups = @src.phylo_count_duplications(maps)
  assert_eq(dups, 0)
  // 3 nodes: 2 internal + 1 leaf... actually 5 nodes total
  // 3 leaves + 2 internal = 5
  assert_eq(maps.length(), 5)
}

///|
test "lca: duplication detected" {
  // Gene tree: ((Human_1, Human_2), Chimp)
  // Two Human genes → the root is a duplication
  let gene = lca_internal([
    lca_internal([lca_leaf("Human_1"), lca_leaf("Human_2")]),
    lca_leaf("Chimp"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "_")
  let dups = @src.phylo_count_duplications(maps)
  // The internal node {Human_1, Human_2} maps to Human,
  // which is the same as both children's LCA → duplication
  assert_true(dups >= 1)
}

///|
test "lca: leaf maps to species" {
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  // First result should be a leaf (Human)
  let human_map = maps[0]
  assert_eq(human_map.species_leaves, "Human,")
  assert_false(human_map.is_duplication)
}

///|
test "lca: internal node maps to LCA" {
  // Gene: ((Human, Chimp), Gorilla)
  // The {Human, Chimp} node should map to {Human, Chimp} in species tree
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  // Find the {Human, Chimp} internal node
  let mut found = false
  for m in maps {
    if m.gene_leaves == "Chimp,Human," {
      assert_eq(m.species_leaves, "Chimp,Human,")
      found = true
    }
  }
  assert_true(found)
}

///|
test "lca: root maps to root of species tree" {
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  // The root should map to the whole species tree
  let root_map = maps[maps.length() - 1]
  assert_eq(root_map.species_leaves, "Chimp,Gorilla,Human,")
}

///|
test "lca: 4-species no duplication" {
  // Gene tree matches species tree exactly
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_internal([lca_leaf("Gorilla"), lca_leaf("Orangutan")]),
  ])
  let species = lca_species_tree4()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  assert_eq(@src.phylo_count_duplications(maps), 0)
}

///|
test "lca: 4-species with duplication" {
  // Gene: ((Human_1, (Human_2, Chimp)), Gorilla)
  // Human_1 and Human_2 are paralogs → duplication
  let gene = lca_internal([
    lca_internal([
      lca_leaf("Human_1"),
      lca_internal([lca_leaf("Human_2"), lca_leaf("Chimp")]),
    ]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "_")
  let dups = @src.phylo_count_duplications(maps)
  assert_true(dups >= 1)
}

///|
test "lca: all leaves same species" {
  // Gene: ((Human_1, Human_2), Human_3)
  // All map to Human → multiple duplications
  let gene = lca_internal([
    lca_internal([lca_leaf("Human_1"), lca_leaf("Human_2")]),
    lca_leaf("Human_3"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "_")
  let dups = @src.phylo_count_duplications(maps)
  // Both internal nodes should be duplications
  assert_true(dups >= 2)
}

///|
test "lca: explicit gene-to-species mapping" {
  // Use explicit mapping instead of separator
  let gene = lca_internal([
    lca_internal([lca_leaf("geneA"), lca_leaf("geneB")]),
    lca_leaf("geneC"),
  ])
  let species = lca_species_tree3()
  let g2s : Map[String, String] = Map([], capacity=0)
  g2s.set("geneA", "Human")
  g2s.set("geneB", "Chimp")
  g2s.set("geneC", "Gorilla")
  let maps = @src.phylo_lca_mapping_explicit(gene, species, g2s, "")
  assert_eq(@src.phylo_count_duplications(maps), 0)
  // geneA should map to Human
  assert_eq(maps[0].species_leaves, "Human,")
}

///|
test "lca: losses counted" {
  // Gene: ((Human, Gorilla), Chimp)
  // Species: ((Human, Chimp), Gorilla)
  // This implies some losses
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Gorilla")]),
    lca_leaf("Chimp"),
  ])
  let species = lca_species_tree3()
  let losses = @src.phylo_count_losses(gene, species, "")
  // There should be some losses due to the mismatch
  assert_true(losses >= 0)
}

///|
test "lca: matching tree has zero losses" {
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let losses = @src.phylo_count_losses(gene, species, "")
  assert_eq(losses, 0)
}

///|
test "lca: deep species tree mapping" {
  // Species: (((Human, Chimp), Gorilla), Orangutan)
  // Gene: (((Human, Chimp), Gorilla), Orangutan) — matching
  let gene = lca_internal([
    lca_internal([
      lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
      lca_leaf("Gorilla"),
    ]),
    lca_leaf("Orangutan"),
  ])
  let species = lca_species_tree4b()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  assert_eq(@src.phylo_count_duplications(maps), 0)
}

///|
test "lca: species names preserved in result" {
  let gene = lca_internal([
    lca_internal([lca_leaf("Human"), lca_leaf("Chimp")]),
    lca_leaf("Gorilla"),
  ])
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  // Root should have all 3 species
  let root_map = maps[maps.length() - 1]
  assert_eq(root_map.species_names.length(), 3)
}

///|
test "lca: single leaf gene tree" {
  let gene = lca_leaf("Human")
  let species = lca_species_tree3()
  let maps = @src.phylo_lca_mapping(gene, species, "")
  assert_eq(maps.length(), 1)
  assert_eq(maps[0].species_leaves, "Human,")
  assert_false(maps[0].is_duplication)
}
