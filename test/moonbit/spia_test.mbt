///|
/// Tests for SPIA module.

///|
test "spia_pathway_node_creation" {
  let n = @src.PathwayNode::new("Gene1", "EGFR")
  assert_eq(n.gene_id, "Gene1")
  assert_eq(n.gene_name, "EGFR")
  assert_eq(n.downstream.length(), 0)
  assert_eq(n.upstream.length(), 0)
  assert_eq(n.activation.length(), 0)
}

///|
test "spia_pathway_node_add_downstream_activation" {
  let n = @src.PathwayNode::new("Gene1", "EGFR")
  n.add_downstream("Gene2", true)
  assert_eq(n.downstream.length(), 1)
  assert_eq(n.downstream[0], "Gene2")
  assert_eq(n.activation[0], 1)
}

///|
test "spia_pathway_node_add_downstream_inhibition" {
  let n = @src.PathwayNode::new("Gene1", "BCL2")
  n.add_downstream("Gene2", false)
  assert_eq(n.downstream.length(), 1)
  assert_eq(n.downstream[0], "Gene2")
  assert_eq(n.activation[0], -1)
}

///|
test "spia_signaling_pathway_creation" {
  let p = @src.SignalingPathway::new("PW1", "MAPK Cascade")
  assert_eq(p.id, "PW1")
  assert_eq(p.name, "MAPK Cascade")
  assert_eq(p.nodes.length(), 0)
}

///|
test "spia_signaling_pathway_add_node" {
  let p = @src.SignalingPathway::new("PW1", "Test Pathway")
  let n1 = @src.PathwayNode::new("Gene1", "EGFR")
  let n2 = @src.PathwayNode::new("Gene2", "RAS")
  p.add_node(n1)
  p.add_node(n2)
  assert_eq(p.nodes.length(), 2)
  assert_eq(p.nodes[0].gene_id, "Gene1")
  assert_eq(p.nodes[1].gene_id, "Gene2")
}

///|
test "spia_signaling_pathway_get_node" {
  let p = @src.SignalingPathway::new("PW1", "Test Pathway")
  let n1 = @src.PathwayNode::new("Gene1", "EGFR")
  let n2 = @src.PathwayNode::new("Gene2", "RAS")
  p.add_node(n1)
  p.add_node(n2)
  match p.get_node("Gene2") {
    Some(n) => assert_eq(n.gene_name, "RAS")
    None => assert_true(false)
  }
  match p.get_node("Gene3") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "spia_sample_pathways" {
  let pathways = @src.spia_sample_pathways()
  assert_eq(pathways.length(), 2)
  assert_eq(pathways[0].id, "PW1")
  assert_eq(pathways[0].name, "MAPK Cascade")
  assert_eq(pathways[1].id, "PW2")
  assert_eq(pathways[1].name, "Apoptosis")
  // PW1 should have 5 nodes
  assert_eq(pathways[0].nodes.length(), 5)
  // PW2 should have 3 nodes
  assert_eq(pathways[1].nodes.length(), 3)
}

///|
test "spia_sample_de_genes" {
  let de = @src.spia_sample_de_genes()
  assert_eq(de.size(), 4)
  assert_true(de.contains("Gene1"))
  assert_true(de.contains("Gene2"))
  assert_true(de.contains("Gene5"))
  assert_true(de.contains("Gene6"))
}

///|
test "spia_run_basic" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  assert_eq(results.n_pathways, 2)
  assert_true(results.results.length() > 0)
  // Each result should have valid structure
  let mut i = 0
  while i < results.results.length() {
    let r = results.results[i]
    assert_true(r.p_or >= 0.0)
    assert_true(r.p_or <= 1.0)
    assert_true(r.p_pert >= 0.0)
    assert_true(r.p_pert <= 1.0)
    assert_true(r.p_combined >= 0.0)
    assert_true(r.p_combined <= 1.0)
    assert_true(r.p_adj >= 0.0)
    assert_true(r.p_adj <= 1.0)
    i = i + 1
  }
}

///|
test "spia_run_pw1_de_genes_count" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  // PW1 (MAPK) contains Gene1, Gene2, Gene3, Gene4, Gene5
  // DE genes: Gene1, Gene2, Gene5 (3 out of 5)
  let mut pw1_result : @src.SpiaResult? = None
  let mut i = 0
  while i < results.results.length() {
    if results.results[i].pathway_id == "PW1" {
      pw1_result = Some(results.results[i])
    }
    i = i + 1
  }
  match pw1_result {
    Some(r) => {
      assert_eq(r.n_genes, 5)
      assert_eq(r.n_de, 3)
    }
    None => assert_true(false)
  }
}

///|
test "spia_run_pw2_de_genes_count" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  // PW2 (Apoptosis) contains Gene6, Gene7, Gene8
  // DE genes: Gene6 (1 out of 3)
  let mut pw2_result : @src.SpiaResult? = None
  let mut i = 0
  while i < results.results.length() {
    if results.results[i].pathway_id == "PW2" {
      pw2_result = Some(results.results[i])
    }
    i = i + 1
  }
  match pw2_result {
    Some(r) => {
      assert_eq(r.n_genes, 3)
      assert_eq(r.n_de, 1)
    }
    None => assert_true(false)
  }
}

///|
test "spia_activation_status" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let mut i = 0
  while i < results.results.length() {
    let r = results.results[i]
    // Activation status must be one of -1, 0, 1
    assert_true(r.activation_status == -1 || r.activation_status == 0 || r.activation_status == 1)
    i = i + 1
  }
}

///|
test "spia_significance_threshold" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  // Use very permissive threshold
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=1.0)
  let sig = results.get_significant()
  assert_eq(sig.length(), results.n_significant)
  // With threshold=1.0, all pathways should be significant
  assert_eq(results.n_significant, results.n_pathways)
}

///|
test "spia_strict_threshold" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  // Use very strict threshold
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.0)
  // No pathway should be significant at threshold 0
  assert_eq(results.n_significant, 0)
}

///|
test "spia_get_activated" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let activated = results.get_activated()
  // Each activated pathway should have activation_status > 0
  let mut i = 0
  while i < activated.length() {
    assert_true(activated[i].activation_status > 0)
    i = i + 1
  }
}

///|
test "spia_get_inhibited" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let inhibited = results.get_inhibited()
  // Each inhibited pathway should have activation_status < 0
  let mut i = 0
  while i < inhibited.length() {
    assert_true(inhibited[i].activation_status < 0)
    i = i + 1
  }
}

///|
test "spia_get_top_pathways" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let top = results.get_top_pathways(1)
  assert_eq(top.length(), 1)
  // Top pathway should have p_combined <= the other's p_combined
  let mut i = 0
  while i < results.results.length() {
    assert_true(top[0].p_combined <= results.results[i].p_combined)
    i = i + 1
  }
}

///|
test "spia_get_top_pathways_more_than_available" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let top = results.get_top_pathways(10)
  // Should return at most n_pathways
  assert_eq(top.length(), results.n_pathways)
}

///|
test "spia_summary" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let s = results.summary()
  assert_true(s.length() > 0)
  // Summary should contain header text
  let has_header = s.contains("SPIA")
  assert_true(has_header)
}

///|
test "spia_empty_pathways" {
  let pathways : Array[@src.SignalingPathway] = Array::new()
  let de : Map[String, Double] = Map::new()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  assert_eq(results.n_pathways, 0)
  assert_eq(results.n_significant, 0)
}

///|
test "spia_no_de_genes" {
  let pathways = @src.spia_sample_pathways()
  let de : Map[String, Double] = Map::new()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  // With no DE genes, all p-values should be 1.0 (not significant)
  assert_eq(results.n_significant, 0)
  let mut i = 0
  while i < results.results.length() {
    assert_true(results.results[i].n_de == 0)
    i = i + 1
  }
}

///|
test "spia_pert_factor_finite" {
  let pathways = @src.spia_sample_pathways()
  let de = @src.spia_sample_de_genes()
  let results = @src.run_spia(pathways, de, 100, fdr_threshold=0.05)
  let mut i = 0
  while i < results.results.length() {
    let r = results.results[i]
    // Perturbation factor should be a finite number (not NaN)
    assert_true(r.pert_factor == r.pert_factor) // NaN != NaN check
    i = i + 1
  }
}
