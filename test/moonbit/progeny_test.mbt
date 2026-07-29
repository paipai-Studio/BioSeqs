///|
test "progeny_matrix_transpose" {
  let matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let result = @src.progeny_transpose(matrix)
  assert_eq(result.length(), 3)
  assert_eq(result[0].length(), 2)
  assert_eq(result[0][0], 1.0)
  assert_eq(result[0][1], 4.0)
  assert_eq(result[1][0], 2.0)
  assert_eq(result[1][1], 5.0)
}

///|
test "progeny_matrix_transpose_empty" {
  let matrix : Array[Array[Double]] = []
  let result = @src.progeny_transpose(matrix)
  assert_eq(result.length(), 0)
}

///|
test "progeny_mat_vec_multiply" {
  let matrix = [[1.0, 2.0], [3.0, 4.0]]
  let vector = [1.0, 2.0]
  let result = @src.progeny_mat_vec(matrix, vector)
  assert_eq(result.length(), 2)
  assert_eq(result[0], 5.0)   // 1*1 + 2*2
  assert_eq(result[1], 11.0)  // 3*1 + 4*2
}

///|
test "progeny_mat_mul" {
  let a = [[1.0, 2.0], [3.0, 4.0]]
  let b = [[5.0, 6.0], [7.0, 8.0]]
  let result = @src.progeny_mat_mul(a, b)
  assert_eq(result.length(), 2)
  assert_eq(result[0][0], 19.0)  // 1*5 + 2*7
  assert_eq(result[0][1], 22.0)  // 1*6 + 2*8
  assert_eq(result[1][0], 43.0)  // 3*5 + 4*7
  assert_eq(result[1][1], 50.0)  // 3*6 + 4*8
}

///|
test "progeny_identity_matrix" {
  let id = @src.progeny_identity(3, 1.0)
  assert_eq(id.length(), 3)
  assert_eq(id[0][0], 1.0)
  assert_eq(id[1][1], 1.0)
  assert_eq(id[2][2], 1.0)
  assert_eq(id[0][1], 0.0)
  assert_eq(id[0][2], 0.0)
  assert_eq(id[1][0], 0.0)
}

///|
test "progeny_mat_add" {
  let a = [[1.0, 2.0], [3.0, 4.0]]
  let b = [[5.0, 6.0], [7.0, 8.0]]
  let result = @src.progeny_mat_add(a, b)
  assert_eq(result[0][0], 6.0)
  assert_eq(result[0][1], 8.0)
  assert_eq(result[1][0], 10.0)
  assert_eq(result[1][1], 12.0)
}

///|
test "progeny_invert_2x2" {
  let matrix = [[2.0, 1.0], [5.0, 3.0]]
  let inv = @src.progeny_invert(matrix)
  assert_eq(inv.length(), 2)
  // Check: A * A^(-1) ≈ I
  let product = @src.progeny_mat_mul(matrix, inv)
  // Allow small tolerance
  assert_eq(product[0][0] > 0.99 && product[0][0] < 1.01, true)
  assert_eq(product[1][1] > 0.99 && product[1][1] < 1.01, true)
  assert_eq(product[0][1] > -0.01 && product[0][1] < 0.01, true)
  assert_eq(product[1][0] > -0.01 && product[1][0] < 0.01, true)
}

///|
test "progeny_invert_3x3" {
  let matrix = [[4.0, 7.0, 2.0], [2.0, 6.0, 3.0], [9.0, 1.0, 5.0]]
  let inv = @src.progeny_invert(matrix)
  assert_eq(inv.length(), 3)
  // Check: A * A^(-1) ≈ I
  let product = @src.progeny_mat_mul(matrix, inv)
  let mut i = 0
  while i < 3 {
    let mut j = 0
    while j < 3 {
      let val = product[i][j]
      if i == j {
        assert_eq(val > 0.99 && val < 1.01, true)
      } else {
        assert_eq(val > -0.01 && val < 0.01, true)
      }
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "progeny_pathway_gene_set" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let weights = [1.0, -1.0, 0.5]
  let pgs = @src.PathwayGeneSet::new("Pathway1", genes, weights)
  assert_eq(pgs.pathway_name, "Pathway1")
  assert_eq(pgs.genes.length(), 3)
  assert_eq(pgs.weights.length(), 3)
}

///|
test "progeny_build_pathway_matrix" {
  let gene_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let pathway1_genes = ["GeneA", "GeneB"]
  let pathway1_weights = [1.0, -1.0]
  let pgs1 = @src.PathwayGeneSet::new("Pathway1", pathway1_genes, pathway1_weights)

  let pathway2_genes = ["GeneB", "GeneC", "GeneD"]
  let pathway2_weights = [0.5, 1.0, -0.5]
  let pgs2 = @src.PathwayGeneSet::new("Pathway2", pathway2_genes, pathway2_weights)

  let pathways = [pgs1, pgs2]
  let matrix = @src.build_pathway_matrix(gene_names, pathways)

  // Matrix should be (genes × pathways) = 4 × 2
  assert_eq(matrix.length(), 4)
  assert_eq(matrix[0].length(), 2)

  // GeneA is in Pathway1 with weight 1.0
  assert_eq(matrix[0][0], 1.0)
  // GeneA is not in Pathway2
  assert_eq(matrix[0][1], 0.0)

  // GeneB is in Pathway1 with weight -1.0
  assert_eq(matrix[1][0], -1.0)
  // GeneB is in Pathway2 with weight 0.5
  assert_eq(matrix[1][1], 0.5)

  // GeneC is not in Pathway1
  assert_eq(matrix[2][0], 0.0)
  // GeneC is in Pathway2 with weight 1.0
  assert_eq(matrix[2][1], 1.0)
}

///|
test "progeny_run_basic" {
  // Expression: 2 samples × 4 genes
  let expression = [
    [5.0, 10.0, 3.0, 8.0],
    [2.0, 6.0, 7.0, 4.0],
  ]
  let sample_names = ["Sample1", "Sample2"]
  let gene_names = ["GeneA", "GeneB", "GeneC", "GeneD"]

  // Pathway1: GeneA + GeneB (activation)
  let pathway1 = @src.PathwayGeneSet::new(
    "Pathway1",
    ["GeneA", "GeneB"],
    [1.0, 1.0],
  )
  // Pathway2: GeneC + GeneD (activation)
  let pathway2 = @src.PathwayGeneSet::new(
    "Pathway2",
    ["GeneC", "GeneD"],
    [1.0, 1.0],
  )

  let data = @src.ProgenyData::new(
    expression,
    sample_names,
    gene_names,
    [pathway1, pathway2],
    0.01,
  )

  let result = @src.run_progeny(data)
  assert_eq(result.sample_names.length(), 2)
  assert_eq(result.pathway_names.length(), 2)
  assert_eq(result.activity_matrix.length(), 2)
  assert_eq(result.activity_matrix[0].length(), 2)

  // Pathway1 activity should be higher for Sample1 (higher GeneA + GeneB)
  let pathway1_activity = @src.get_pathway_activity(result, "Pathway1")
  assert_eq(pathway1_activity.length(), 2)
  assert_eq(pathway1_activity[0] > pathway1_activity[1], true)
}

///|
test "progeny_get_pathway_activity" {
  let expression = [[10.0, 5.0], [3.0, 8.0]]
  let sample_names = ["S1", "S2"]
  let gene_names = ["G1", "G2"]
  let pathway1 = @src.PathwayGeneSet::new("P1", ["G1"], [1.0])
  let data = @src.ProgenyData::new(expression, sample_names, gene_names, [pathway1], 0.01)
  let result = @src.run_progeny(data)

  let activity = @src.get_pathway_activity(result, "P1")
  assert_eq(activity.length(), 2)
  assert_eq(activity[0] > activity[1], true)
}

///|
test "progeny_get_sample_profile" {
  let expression = [[5.0, 3.0], [2.0, 8.0]]
  let sample_names = ["S1", "S2"]
  let gene_names = ["G1", "G2"]
  let p1 = @src.PathwayGeneSet::new("P1", ["G1"], [1.0])
  let p2 = @src.PathwayGeneSet::new("P2", ["G2"], [1.0])
  let data = @src.ProgenyData::new(expression, sample_names, gene_names, [p1, p2], 0.01)
  let result = @src.run_progeny(data)

  let profile = @src.get_sample_profile(result, "S1")
  assert_eq(profile.length(), 2)
}
