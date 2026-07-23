///|
test "genefilter_row_ttest" {
  let gene_ids = ["gene1", "gene2", "gene3"]
  let expression = [
    [10.0, 12.0, 8.0, 20.0, 22.0, 18.0],
    [5.0, 6.0, 5.5, 7.0, 6.5, 7.5],
    [100.0, 98.0, 102.0, 50.0, 48.0, 52.0],
  ]
  let groups = ["control", "control", "control", "treatment", "treatment", "treatment"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  let result = @src.row_ttest(ge, "control", "treatment")
  
  assert_eq(result.passing_genes.length(), 3)
}

///|
test "genefilter_row_wilcoxon" {
  let gene_ids = ["gene1", "gene2"]
  let expression = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    [10.0, 11.0, 12.0, 1.0, 2.0, 3.0],
  ]
  let groups = ["groupA", "groupA", "groupA", "groupB", "groupB", "groupB"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  let result = @src.row_wilcoxon(ge, "groupA", "groupB")
  
  assert_eq(result.passing_genes.length(), 2)
}

///|
test "genefilter_variance_filter" {
  let gene_ids = ["gene1", "gene2", "gene3"]
  let expression = [
    [1.0, 1.0, 1.0],
    [1.0, 2.0, 3.0],
    [10.0, 20.0, 30.0],
  ]
  let groups = ["g1", "g1", "g1"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  let passing = @src.variance_filter(ge, 1.0)
  
  assert_eq(passing.length() >= 2, true)
}

///|
test "genefilter_cv_filter" {
  let gene_ids = ["gene1", "gene2"]
  let expression = [
    [1.0, 1.0, 1.0],
    [1.0, 2.0, 4.0],
  ]
  let groups = ["g1", "g1", "g1"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  let passing = @src.cv_filter(ge, 0.5)
  
  assert_eq(passing.length() >= 1, true)
}

///|
test "genefilter_row_quantile_filter" {
  let gene_ids = ["gene1", "gene2", "gene3"]
  let expression = [
    [5.0, 6.0, 7.0],
    [1.0, 2.0, 3.0],
    [9.0, 10.0, 11.0],
  ]
  let groups = ["g1", "g1", "g1"]
  
  let ge = @src.GeneExpression::new(gene_ids, expression, groups)
  let passing = @src.row_quantile_filter(ge, 0.25, 0.75)
  
  assert_eq(passing.length() >= 1, true)
}