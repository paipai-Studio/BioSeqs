///|
fn main {
  println("=== MAST Demo ===")

  let n_genes = 10
  let n_group1 = 6
  let n_group2 = 6
  let n_cells = n_group1 + n_group2

  let expression : Array[Array[Double]] = Array::make(
    n_cells,
    Array::make(n_genes, 0.0),
  )
  let gene_names : Array[String] = Array::make(n_genes, "")
  let cell_names : Array[String] = Array::make(n_cells, "")
  let groups : Array[Int] = Array::make(n_cells, 0)

  let mut g = 0
  while g < n_genes {
    gene_names[g] = "Gene" + g.to_string()
    g = g + 1
  }

  let mut i = 0
  while i < n_cells {
    cell_names[i] = "Cell" + i.to_string()
    groups[i] = if i < n_group1 { 0 } else { 1 }

    let row : Array[Double] = Array::make(n_genes, 0.0)
    let mut j = 0
    while j < n_genes {
      let base_expr = if j < 5 {
        if i < n_group1 {
          10.0
        } else {
          50.0
        }
      } else if i < n_group1 {
        30.0
      } else {
        25.0
      }
      let noise = (i * 7 + j * 13).to_double() % 5.0
      let zero_prob = if j < 3 {
        if i < n_group1 {
          0.4
        } else {
          0.1
        }
      } else {
        0.2
      }
      let is_zero = ((i * 31 + j * 17) % 100).to_double() / 100.0 < zero_prob
      row[j] = if is_zero { 0.0 } else { base_expr + noise }
      j = j + 1
    }
    expression[i] = row
    i = i + 1
  }

  println("\n1. Create MAST Data:")
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  println("  Number of cells: " + data.n_cells.to_string())
  println("  Number of genes: " + data.n_genes.to_string())
  println("  Number of groups: " + data.n_groups.to_string())

  println("\n2. Detection rate (cngeneson) per cell:")
  let mut c = 0
  while c < 5 && c < data.n_cells {
    println(
      "  " +
      cell_names[c] +
      ": " +
      (data.cngeneson[c] * 100.0).to_string() +
      "%",
    )
    c = c + 1
  }

  println("\n3. Run Hurdle Model (group 0 vs group 1):")
  let result = @src.zlm_hurdle(data, 0, 1)
  println(result.to_string())

  println("\n4. Detailed results for all genes:")
  let mut k = 0
  while k < result.gene.length() {
    println("  " + result.gene[k] + ":")
    println(
      "    Detection rate: " +
      (result.detection_rate[k] * 100.0).to_string() +
      "%",
    )
    println("    Discrete p-value: " + result.discrete_pval[k].to_string())
    println("    Continuous p-value: " + result.continuous_pval[k].to_string())
    println("    Combined p-value: " + result.combined_pval[k].to_string())
    println("    Combined FDR: " + result.combined_fdr[k].to_string())
    println("    Log2 FC: " + result.log2_fc[k].to_string())
    k = k + 1
  }

  println("\n5. MAST Summary (FDR < 0.05):")
  let summary = @src.MastSummary::new(result, 0.05)
  println(summary.to_string())

  println("\n6. Top genes by combined p-value:")
  let indices : Array[Int] = Array::make(result.gene.length(), 0)
  let mut m = 0
  while m < result.gene.length() {
    indices[m] = m
    m = m + 1
  }

  let mut p = 0
  while p < 5 && p < result.gene.length() {
    let mut best_idx = p
    let mut best_pval = result.combined_pval[p]
    let mut q = p + 1
    while q < result.gene.length() {
      if result.combined_pval[q] < best_pval {
        best_pval = result.combined_pval[q]
        best_idx = q
      }
      q = q + 1
    }
    let temp = indices[p]
    indices[p] = indices[best_idx]
    indices[best_idx] = temp

    let idx = indices[p]
    println(
      "  " +
      (p + 1).to_string() +
      ". " +
      result.gene[idx] +
      " (p=" +
      result.combined_pval[idx].to_string() +
      ", FDR=" +
      result.combined_fdr[idx].to_string() +
      ")",
    )
    p = p + 1
  }

  println("\nMAST demo completed!")
}
