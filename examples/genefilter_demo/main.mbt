///|
/// genefilter module demo.

///|
pub fn main() {
  println("=== genefilter Demo ===")
  
  println("\n1. Creating GeneExpression data")
  let gene_ids = ["gene1", "gene2", "gene3", "gene4", "gene5"]
  let expression = [
    [10.5, 12.3, 11.1, 25.0, 27.5, 26.2],
    [5.8, 6.2, 5.5, 7.0, 6.8, 7.2],
    [100.2, 98.5, 102.1, 50.0, 48.5, 52.0],
    [15.7, 16.2, 15.9, 14.8, 15.2, 14.5],
    [3.2, 3.5, 3.1, 3.0, 3.3, 2.9],
  ]
  let groups = ["control", "control", "control", "treatment", "treatment", "treatment"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  println("Genes: {}", ge.gene_ids.length())
  println("Samples: {}", ge.groups.length())
  
  println("\n2. t-test filtering")
  let t_result = @src.row_ttest(ge, "control", "treatment")
  println("Genes passed t-test: {}", t_result.passing_genes.length())
  
  println("\n3. Wilcoxon rank sum test")
  let w_result = @src.row_wilcoxon(ge, "control", "treatment")
  println("Genes passed Wilcoxon: {}", w_result.passing_genes.length())
  
  println("\n4. Variance filter")
  let var_passing = @src.variance_filter(ge, 10.0)
  println("Genes with variance >= 10: {}", var_passing.length())
  
  println("\n5. CV filter")
  let cv_passing = @src.cv_filter(ge, 0.3)
  println("Genes with CV >= 0.3: {}", cv_passing.length())
  
  println("\n=== Demo complete ===")
}