///|
/// Tests for Ortholog/Paralog classification.

///|
fn orthopara_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(1.0))
}

///|
fn orthopara_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("i"), clades=children, branch_length=Some(1.0))
}

///|
/// Species tree: ((Human, Chimp), Gorilla)
fn orthopara_species3() -> @src.Clade {
  orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
}

///|
/// Species tree: ((Human, Chimp), (Gorilla, Orangutan))
fn orthopara_species4() -> @src.Clade {
  orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_internal([orthopara_leaf("Gorilla"), orthopara_leaf("Orangutan")]),
  ])
}

test "orthoparalog: simple orthologs" {
  // Gene tree = species tree: all pairs are orthologs
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "")
  // 3 leaves → 3 pairs, all orthologs
  assert_eq(pairs.length(), 3)
  for p in pairs {
    assert_true(match p.relationship {
      Ortholog => true
      _ => false
    })
  }
}

///|
test "orthoparalog: human-chimp ortholog" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human", "Chimp", "")
  match cls {
    Some(c) => assert_true(match c.relationship {
      Ortholog => true
      _ => false
    })
    None => assert_true(false)
  }
}

///|
test "orthoparalog: human-gorilla ortholog" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human", "Gorilla", "")
  match cls {
    Some(c) => assert_true(match c.relationship {
      Ortholog => true
      _ => false
    })
    None => assert_true(false)
  }
}

///|
test "orthoparalog: simple paralog" {
  // Gene: ((Human_1, Human_2), Chimp)
  // Human_1 and Human_2 are paralogs (duplication)
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human_1", "Human_2", "_")
  match cls {
    Some(c) => assert_true(match c.relationship {
      SuperParalog => true
      UltraParalog => true
      _ => false
    })
    None => assert_true(false)
  }
}

///|
test "orthoparalog: paralog vs ortholog mix" {
  // Gene: ((Human_1, Human_2), Chimp)
  // Human_1-Human_2: paralog (same species)
  // Human_1-Chimp: ortholog
  // Human_2-Chimp: ortholog
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "_")
  assert_eq(pairs.length(), 3)
  let mut ortho_count = 0
  let mut para_count = 0
  for p in pairs {
    match p.relationship {
      Ortholog => ortho_count = ortho_count + 1
      SuperParalog => para_count = para_count + 1
      UltraParalog => para_count = para_count + 1
      Paralog => para_count = para_count + 1
    }
  }
  assert_eq(ortho_count, 2)
  assert_eq(para_count, 1)
}

///|
test "orthoparalog: out-paralog" {
  // Gene: ((Human_1, Chimp_1), (Human_2, Gorilla_1))
  // Species: ((Human, Chimp), Gorilla)
  // Root is duplication (right child LCA = root = parent LCA)
  // Human_1 (Human) and Gorilla_1 (Gorilla): different species → out-paralog
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Chimp_1")]),
    orthopara_internal([orthopara_leaf("Human_2"), orthopara_leaf("Gorilla_1")]),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human_1", "Gorilla_1", "_")
  match cls {
    Some(c) => assert_true(match c.relationship {
      Paralog => true
      _ => false
    })
    None => assert_true(false)
  }
}

///|
test "orthoparalog: all pairs count" {
  // 4 leaves → 6 pairs
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_internal([orthopara_leaf("Gorilla"), orthopara_leaf("Orangutan")]),
  ])
  let species = orthopara_species4()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "")
  assert_eq(pairs.length(), 6)
}

///|
test "orthoparalog: 4-species all orthologs" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_internal([orthopara_leaf("Gorilla"), orthopara_leaf("Orangutan")]),
  ])
  let species = orthopara_species4()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "")
  for p in pairs {
    assert_true(match p.relationship {
      Ortholog => true
      _ => false
    })
  }
}

///|
test "orthoparalog: extract orthologs" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "_")
  let orthos = @src.phylo_orthologs(pairs)
  assert_eq(orthos.length(), 2)
}

///|
test "orthoparalog: extract paralogs" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "_")
  let paras = @src.phylo_paralogs(pairs)
  assert_eq(paras.length(), 1)
}

///|
test "orthoparalog: pair not found" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human", "Mouse", "")
  match cls {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "orthoparalog: species names in result" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp_1"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human_1", "Chimp_1", "_")
  match cls {
    Some(c) => {
      assert_eq(c.species1, "Human")
      assert_eq(c.species2, "Chimp")
    }
    None => assert_true(false)
  }
}

///|
test "orthoparalog: summary output" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "_")
  let summary = @src.phylo_pairs_summary(pairs)
  assert_true(summary.contains("Ortholog"))
  assert_true(summary.contains("Paralog"))
  assert_true(summary.contains("Total pairs: 3"))
}

///|
test "orthoparalog: pair to string" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("Human"), orthopara_leaf("Chimp")]),
    orthopara_leaf("Gorilla"),
  ])
  let species = orthopara_species3()
  let cls = @src.phylo_classify_pair(gene, species, "Human", "Chimp", "")
  match cls {
    Some(c) => {
      let s = @src.phylo_pair_to_string(c)
      assert_true(s.contains("Human"))
      assert_true(s.contains("Chimp"))
      assert_true(s.contains("Ortholog"))
    }
    None => assert_true(false)
  }
}

///|
test "orthoparalog: explicit mapping" {
  let gene = orthopara_internal([
    orthopara_internal([orthopara_leaf("g1"), orthopara_leaf("g2")]),
    orthopara_leaf("g3"),
  ])
  let species = orthopara_species3()
  let g2s : Map[String, String] = Map([], capacity=0)
  g2s.set("g1", "Human")
  g2s.set("g2", "Chimp")
  g2s.set("g3", "Gorilla")
  let pairs = @src.phylo_classify_all_pairs_explicit(gene, species, g2s, "")
  for p in pairs {
    assert_true(match p.relationship {
      Ortholog => true
      _ => false
    })
  }
}

///|
test "orthoparalog: nested duplication" {
  // Gene: (((Human_1, Human_2), Human_3), Chimp)
  // Two nested duplications, all Human pairs are super/ultra paralogs
  let gene = orthopara_internal([
    orthopara_internal([
      orthopara_internal([orthopara_leaf("Human_1"), orthopara_leaf("Human_2")]),
      orthopara_leaf("Human_3"),
    ]),
    orthopara_leaf("Chimp"),
  ])
  let species = orthopara_species3()
  let pairs = @src.phylo_classify_all_pairs(gene, species, "_")
  // 4 leaves → 6 pairs
  assert_eq(pairs.length(), 6)
  // Human_1-Human_2, Human_1-Human_3, Human_2-Human_3: paralogs
  // Human_*-Chimp: orthologs
  let mut ortho_count = 0
  let mut para_count = 0
  for p in pairs {
    match p.relationship {
      Ortholog => ortho_count = ortho_count + 1
      _ => para_count = para_count + 1
    }
  }
  assert_eq(ortho_count, 3)
  assert_eq(para_count, 3)
}
