///|
/// Tests for the velociraptor module (Bioconductor velociraptor / scVelo).
///
/// Covers gene data construction, steady-state model, kinetic EM model,
/// velocity computation, embedding projection, transition matrix,
/// root cell identification, and sample data validation.

// ===========================================================================
// 1. VelocityGeneData construction
// ===========================================================================

///|
test "vr_gene_data_creation" {
  let gd = @src.VelocityGeneData::new("GeneA", [1.0, 2.0, 3.0], [0.5, 1.0, 1.5])
  assert_eq(gd.gene_name(), "GeneA")
  assert_eq(gd.spliced().length(), 3)
  assert_eq(gd.unspliced().length(), 3)
  assert_eq(gd.spliced()[0], 1.0)
  assert_eq(gd.unspliced()[0], 0.5)
}

// ===========================================================================
// 2. Sample data
// ===========================================================================

///|
test "vr_sample_data_dimensions" {
  let gene_data = @src.velocity_sample_data()
  assert_eq(gene_data.length(), 4) // 4 genes
  for gd in gene_data {
    assert_eq(gd.spliced().length(), 20) // 20 cells
    assert_eq(gd.unspliced().length(), 20)
  }
}

///|
test "vr_sample_data_nonneg" {
  let gene_data = @src.velocity_sample_data()
  for gd in gene_data {
    for s in gd.spliced() {
      assert_true(s >= 0.0)
    }
    for u in gd.unspliced() {
      assert_true(u >= 0.0)
    }
  }
}

///|
test "vr_sample_embeddings" {
  let emb = @src.velocity_sample_embeddings(10)
  assert_eq(emb.length(), 10)
  for e in emb {
    assert_eq(e.length(), 2) // 2D embedding
  }
}

// ===========================================================================
// 3. Steady-state model
// ===========================================================================

///|
test "vr_steady_state_returns_ratio" {
  let spliced = [1.0, 2.0, 3.0, 4.0, 5.0]
  let unspliced = [2.0, 4.0, 5.0, 6.0, 7.0]
  let (ratio, ss_s, ss_u) = @src.velocity_steady_state(spliced, unspliced)
  assert_true(ratio > 0.0)
  assert_true(ss_s > 0.0)
  assert_true(ss_u > 0.0)
}

///|
test "vr_steady_state_perfect_linear" {
  // Perfect linear relationship: u = 2 * s.
  let spliced = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let unspliced = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let (ratio, _, _) = @src.velocity_steady_state(spliced, unspliced)
  // Ratio should be close to 2.0.
  assert_true(ratio > 1.5)
  assert_true(ratio < 2.5)
}

///|
test "vr_steady_state_too_few_cells" {
  let (ratio, ss_s, ss_u) = @src.velocity_steady_state([1.0, 2.0], [1.0, 2.0])
  // With < 3 cells, returns defaults.
  assert_eq(ratio, 1.0)
  assert_eq(ss_s, 0.0)
  assert_eq(ss_u, 0.0)
}

// ===========================================================================
// 4. Kinetic model (EM)
// ===========================================================================

///|
test "vr_kinetic_model_returns_params" {
  let spliced = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let unspliced = [5.0, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0]
  let k = @src.velocity_kinetic_model("GeneA", spliced, unspliced, 50)
  assert_eq(k.gene_name(), "GeneA")
  assert_true(k.alpha() >= 0.0)
  assert_true(k.beta() > 0.0)
  assert_true(k.gamma() > 0.0)
}

///|
test "vr_kinetic_model_sample_data" {
  let gene_data = @src.velocity_sample_data()
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    assert_true(k.beta() > 0.0)
    assert_true(k.gamma() > 0.0)
  }
}

///|
test "vr_kinetic_model_too_few_cells" {
  let k = @src.velocity_kinetic_model("Gene", [1.0, 2.0], [1.0, 2.0], 10)
  assert_eq(k.gene_name(), "Gene")
  assert_eq(k.alpha(), 0.0)
  assert_eq(k.beta(), 1.0)
  assert_eq(k.gamma(), 1.0)
}

// ===========================================================================
// 5. Velocity computation
// ===========================================================================

///|
test "vr_compute_returns_velocities" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  assert_eq(velocities.length(), 20) // 20 cells
  for v in velocities {
    assert_eq(v.velocity().length(), 4) // 4 genes
    assert_true(v.speed() >= 0.0)
  }
}

///|
test "vr_gene_velocity" {
  let spliced = [1.0, 2.0, 3.0]
  let unspliced = [2.0, 4.0, 6.0]
  let vel = @src.velocity_gene(spliced, unspliced, 1.0, 0.5)
  assert_eq(vel.length(), 3)
  // v = beta*u - gamma*s = 1*2 - 0.5*1 = 1.5
  assert_eq(vel[0], 1.5)
}

// ===========================================================================
// 6. Velocity embedding
// ===========================================================================

///|
test "vr_velocity_embedding" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  let embeddings = @src.velocity_sample_embeddings(20)
  let emb_vel = @src.velocity_embedding(velocities, embeddings, 5)
  assert_eq(emb_vel.length(), 20)
  for ev in emb_vel {
    assert_eq(ev.embedding().length(), 2)
    assert_eq(ev.velocity_embedding().length(), 2)
  }
}

///|
test "vr_embedding_empty_input" {
  let velocities : Array[@src.CellVelocity] = []
  let emb_vel = @src.velocity_embedding(velocities, [], 5)
  assert_eq(emb_vel.length(), 0)
}

// ===========================================================================
// 7. Root cell identification
// ===========================================================================

///|
test "vr_find_root_cells" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  let roots = @src.velocity_find_root_cells(velocities, 3)
  assert_eq(roots.length(), 3)
  // Root indices should be valid.
  for r in roots {
    assert_true(r >= 0)
    assert_true(r < 20)
  }
}

///|
test "vr_root_cells_sorted_by_speed" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  let roots = @src.velocity_find_root_cells(velocities, 5)
  // Speeds should be sorted descending.
  for i in 1..<roots.length() {
    assert_true(velocities[roots[i - 1]].speed() >= velocities[roots[i]].speed())
  }
}

// ===========================================================================
// 8. Transition matrix
// ===========================================================================

///|
test "vr_transition_matrix" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  let embeddings = @src.velocity_sample_embeddings(20)
  let matrix = @src.velocity_transition_matrix(velocities, embeddings, 5, 1.0)
  assert_eq(matrix.length(), 20)
  for row in matrix {
    assert_eq(row.length(), 20)
  }
  // Each row should sum to approximately 1.0 (probabilities).
  for row in matrix {
    let mut sum = 0.0
    for v in row {
      sum = sum + v
    }
    // Sum should be close to 1.0 (within floating point tolerance).
    assert_true(sum > 0.99 || sum == 0.0) // 0.0 if no neighbors
  }
}

///|
test "vr_transition_matrix_nonneg" {
  let gene_data = @src.velocity_sample_data()
  let kinetics : Array[@src.GeneKinetics] = []
  for gd in gene_data {
    let k = @src.velocity_kinetic_model(gd.gene_name(), gd.spliced(), gd.unspliced(), 30)
    kinetics.push(k)
  }
  let velocities = @src.velocity_compute(gene_data, kinetics)
  let embeddings = @src.velocity_sample_embeddings(20)
  let matrix = @src.velocity_transition_matrix(velocities, embeddings, 5, 1.0)
  for row in matrix {
    for v in row {
      assert_true(v >= 0.0)
    }
  }
}

// ===========================================================================
// 9. Full pipeline
// ===========================================================================

///|
test "vr_run_pipeline" {
  let gene_data = @src.velocity_sample_data()
  let embeddings = @src.velocity_sample_embeddings(20)
  let (kinetics, velocities, emb_vel) = @src.velocity_run(gene_data, embeddings, 5)
  assert_eq(kinetics.length(), 4) // 4 genes
  assert_eq(velocities.length(), 20) // 20 cells
  assert_eq(emb_vel.length(), 20)
  // Check kinetics are valid.
  for k in kinetics {
    assert_true(k.beta() > 0.0)
    assert_true(k.gamma() > 0.0)
  }
}

///|
test "vr_run_without_embedding" {
  let gene_data = @src.velocity_sample_data()
  let (kinetics, velocities, emb_vel) = @src.velocity_run(gene_data, [], 5)
  assert_eq(kinetics.length(), 4)
  assert_eq(velocities.length(), 20)
  assert_eq(emb_vel.length(), 0) // no embedding => no projected velocity
}

// ===========================================================================
// 10. GeneKinetics accessors
// ===========================================================================

///|
test "vr_kinetics_accessors" {
  let k = @src.velocity_kinetic_model(
    "TestGene",
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 3.0, 2.0, 1.0],
    20,
  )
  assert_eq(k.gene_name(), "TestGene")
  assert_true(k.alpha() >= 0.0)
  assert_true(k.beta() > 0.0)
  assert_true(k.gamma() > 0.0)
  assert_true(k.likelihood() <= 0.0) // negative log-likelihood
}

// ===========================================================================
// 11. Edge cases
// ===========================================================================

///|
test "vr_empty_gene_data" {
  let gene_data : Array[@src.VelocityGeneData] = []
  let (kinetics, velocities, _) = @src.velocity_run(gene_data, [], 5)
  assert_eq(kinetics.length(), 0)
  assert_eq(velocities.length(), 0)
}

///|
test "vr_single_gene_single_cell" {
  let gd = @src.VelocityGeneData::new("G1", [5.0], [3.0])
  let k = @src.velocity_kinetic_model("G1", [5.0], [3.0], 10)
  // With < 3 cells, returns defaults.
  assert_eq(k.alpha(), 0.0)
}

///|
test "vr_velocity_gene_zero_params" {
  let vel = @src.velocity_gene([1.0, 2.0], [3.0, 4.0], 0.0, 0.0)
  assert_eq(vel[0], 0.0)
  assert_eq(vel[1], 0.0)
}
