///|
/// dorothea Demo - TF activity prediction using VIPER.
/// Inspired by Bioconductor dorothea package.

///|
fn main {
  println("=== dorothea TF Activity Prediction Demo ===")
  println("Inspired by Bioconductor dorothea package")
  println("")
  
  println("1. Loading regulons (TF → target gene sets)...")
  let regulons = @bio.dorothea_get_regulons()
  println("   Total regulons: \{regulons.length().to_string()}")
  for regulon in regulons {
    let mut n_activation = 0
    let mut n_repression = 0
    for target in regulon.targets {
      if target.direction == "activation" {
        n_activation = n_activation + 1
      } else {
        n_repression = n_repression + 1
      }
    }
    println("     - \{regulon.tf}: \{regulon.targets.length().to_string()} targets (\{n_activation.to_string()} activating, \{n_repression.to_string()} repressing)")
  }
  println("")
  
  println("2. Creating synthetic single-cell data...")
  let test_data = @bio.dorothea_create_test_data()
  let (expr_data, cell_ids, gene_names) = test_data
  println("   Cells: \{cell_ids.length().to_string()}")
  println("   Genes: \{gene_names.length().to_string()}")
  println("")
  
  println("3. Running dorothea TF activity analysis...")
  let params = @bio.DorotheaParams::with_params(50, 3, 0.05)
  let results = @bio.dorothea_analyze_cells(expr_data, cell_ids, gene_names, regulons, params)
  println("   Analyzed TF activity for \{results.length().to_string()} TFs")
  println("")
  
  println("4. TF activity results:")
  let sorted = @bio.dorothea_sort_by_activity(results)
  for result in sorted {
    println("   \{result.tf}: activity=\{result.activity.to_string()}, z-score=\{result.z_score.to_string()}, targets=\{result.n_targets.to_string()}")
  }
  println("")
  
  println("5. Top 5 TFs by absolute activity:")
  let top5 = @bio.dorothea_get_top_tfs(results, 5)
  for result in top5 {
    println("   \{result.tf}: z-score=\{result.z_score.to_string()}")
  }
  println("")
  
  println("6. Analyzing single cell in detail...")
  let cell_expr : Array[Double] = Array::new()
  let n_genes = expr_data.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(expr_data[g][0])
    g = g + 1
  }
  let cell_results = @bio.dorothea_analyze_cell(cell_expr, gene_names, regulons, params)
  println("   Single cell TF activities:")
  for result in cell_results {
    println("     \{result.tf}: activity=\{result.activity.to_string()}")
  }
  println("")
  
  println("=== Demo completed successfully! ===")
  println("")
  println("Key takeaways:")
  println("- dorothea predicts TF activity from gene expression data")
  println("- Regulons define target genes and their direction (activation/repression)")
  println("- VIPER scoring uses the collective expression of target genes")
  println("- Permutation testing provides significance estimates")
  println("- The z-score indicates how far the activity is from the null expectation")
}