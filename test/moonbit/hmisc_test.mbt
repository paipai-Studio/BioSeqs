///|
/// Test file for hmisc module.

// hmisc_describe tests

test "hmisc_describe basic valid data" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_describe(data, name="test_var")
  assert_eq(result.name, "test_var")
  assert_eq(result.n, 5)
  assert_eq(result.n_missing, 0)
  assert_true(result.mean > 2.9 && result.mean < 3.1)
  assert_true(result.sd > 1.4 && result.sd < 1.6)
  assert_eq(result.min, 1.0)
  assert_eq(result.max, 5.0)
  assert_true(result.median > 2.9 && result.median < 3.1)
}

test "hmisc_describe with NaN values" {
  let data = [1.0, @double.not_a_number, 3.0, @double.not_a_number, 5.0]
  let result = @src.hmisc_describe(data)
  assert_eq(result.n, 3)
  assert_eq(result.n_missing, 2)
  assert_true(result.mean > 2.9 && result.mean < 3.1)
  assert_eq(result.min, 1.0)
  assert_eq(result.max, 5.0)
}

test "hmisc_describe empty data" {
  let data : Array[Double] = Array::new()
  let result = @src.hmisc_describe(data)
  assert_eq(result.n, 0)
  assert_eq(result.n_missing, 0)
  assert_eq(result.mean, 0.0)
  assert_eq(result.sd, 0.0)
  assert_eq(result.min, 0.0)
  assert_eq(result.max, 0.0)
}

test "hmisc_describe all NaN data" {
  let data = [@double.not_a_number, @double.not_a_number, @double.not_a_number]
  let result = @src.hmisc_describe(data)
  assert_eq(result.n, 0)
  assert_eq(result.n_missing, 3)
  assert_eq(result.mean, 0.0)
  assert_eq(result.sd, 0.0)
}

test "hmisc_describe single value" {
  let data = [42.0]
  let result = @src.hmisc_describe(data)
  assert_eq(result.n, 1)
  assert_eq(result.mean, 42.0)
  assert_eq(result.sd, 0.0)
  assert_eq(result.min, 42.0)
  assert_eq(result.max, 42.0)
}

test "hmisc_describe default name" {
  let data = [1.0, 2.0, 3.0]
  let result = @src.hmisc_describe(data)
  assert_eq(result.name, "")
}

test "hmisc_describe quartiles" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.hmisc_describe(data)
  assert_true(result.q1 > 3.0 && result.q1 < 4.0)
  assert_true(result.median > 5.0 && result.median < 6.0)
  assert_true(result.q3 > 7.0 && result.q3 < 8.0)
}

// hmisc_pearson tests

test "hmisc_pearson perfect positive correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

test "hmisc_pearson perfect negative correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [10.0, 8.0, 6.0, 4.0, 2.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result < -0.99 && result > -1.01)
}

test "hmisc_pearson no correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result > -0.5 && result < 0.5)
}

test "hmisc_pearson mismatched lengths" {
  let x = [1.0, 2.0, 3.0, 4.0]
  let y = [1.0, 2.0, 3.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result.is_nan())
}

test "hmisc_pearson small arrays less than 3" {
  let x = [1.0, 2.0]
  let y = [3.0, 4.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result.is_nan())
}

test "hmisc_pearson constant values" {
  let x = [5.0, 5.0, 5.0, 5.0, 5.0]
  let y = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result.is_nan())
}

// hmisc_spearman tests

test "hmisc_spearman perfect positive correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

test "hmisc_spearman perfect negative correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [10.0, 8.0, 6.0, 4.0, 2.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result < -0.99 && result > -1.01)
}

test "hmisc_spearman with ties" {
  let x = [1.0, 2.0, 2.0, 3.0, 4.0]
  let y = [2.0, 3.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result > 0.9 && result < 1.01)
}

test "hmisc_spearman mismatched lengths" {
  let x = [1.0, 2.0, 3.0, 4.0]
  let y = [1.0, 2.0, 3.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result.is_nan())
}

test "hmisc_spearman small arrays less than 3" {
  let x = [1.0, 2.0]
  let y = [3.0, 4.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result.is_nan())
}

test "hmisc_spearman monotonic but non-linear" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

// hmisc_rcorr tests

test "hmisc_rcorr basic pearson correlation" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
    [5.0, 4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.hmisc_rcorr(data, cor_type=@src.hmisc_cor_type_pearson())
  assert_eq(result.names.length(), 3)
  assert_eq(result.matrix.length(), 3)
  assert_eq(result.p_values.length(), 3)
  assert_eq(result.matrix[0][0], 1.0)
  assert_eq(result.matrix[1][1], 1.0)
  assert_eq(result.matrix[2][2], 1.0)
  assert_true(result.matrix[0][1] > 0.99)
  assert_true(result.matrix[0][2] < -0.99)
  assert_true(result.matrix[1][2] < -0.99)
}

test "hmisc_rcorr basic spearman correlation" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_rcorr(data, cor_type=@src.hmisc_cor_type_spearman())
  assert_eq(@src.hmisc_cor_type_to_string(result.cor_type), "Spearman")
  assert_true(result.matrix[0][1] > 0.99)
}

test "hmisc_rcorr with custom names" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_rcorr(data, names=["Height", "Weight"])
  assert_eq(result.names[0], "Height")
  assert_eq(result.names[1], "Weight")
}

test "hmisc_rcorr default names" {
  let data = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 4.0, 6.0, 8.0],
    [3.0, 6.0, 9.0, 12.0],
  ]
  let result = @src.hmisc_rcorr(data)
  assert_eq(result.names[0], "V1")
  assert_eq(result.names[1], "V2")
  assert_eq(result.names[2], "V3")
}

test "hmisc_rcorr with NaN pairwise deletion" {
  let data = [
    [1.0, 2.0, @double.not_a_number, 4.0, 5.0],
    [2.0, @double.not_a_number, 6.0, 8.0, 10.0],
    [5.0, 4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.hmisc_rcorr(data)
  assert_true(!result.matrix[0][1].is_nan())
  assert_true(!result.matrix[0][2].is_nan())
  assert_true(!result.matrix[1][2].is_nan())
}

test "hmisc_rcorr p-values computed" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
  ]
  let result = @src.hmisc_rcorr(data)
  assert_true(!result.p_values[0][1].is_nan())
  assert_true(result.p_values[0][1] >= 0.0 && result.p_values[0][1] <= 1.0)
}

test "hmisc_rcorr symmetry" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
    [5.0, 4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.hmisc_rcorr(data)
  assert_true(result.matrix[0][1] == result.matrix[1][0])
  assert_true(result.matrix[0][2] == result.matrix[2][0])
  assert_true(result.matrix[1][2] == result.matrix[2][1])
}

test "hmisc_rcorr n matrix" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_rcorr(data)
  assert_eq(result.n[0][0], 5)
  assert_eq(result.n[1][1], 5)
  assert_eq(result.n[0][1], 5)
  assert_eq(result.n[1][0], 5)
}

// hmisc_varclus tests

test "hmisc_varclus basic clustering" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
    [5.0, 4.0, 6.0, 3.0, 7.0, 5.0, 8.0, 4.0, 9.0, 6.0],
  ]
  let result = @src.hmisc_varclus(data, names=["X", "Y", "Z"], min_cluster_size=2)
  assert_eq(result.names.length(), 3)
  assert_true(result.n_clusters >= 1)
  assert_eq(result.cluster_assignments.length(), 3)
}

test "hmisc_varclus all variables in one cluster" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
    [3.0, 6.0, 9.0, 12.0, 15.0],
  ]
  let result = @src.hmisc_varclus(data, min_cluster_size=3)
  assert_eq(result.n_clusters, 1)
  assert_eq(result.cluster_assignments[0], 0)
  assert_eq(result.cluster_assignments[1], 0)
  assert_eq(result.cluster_assignments[2], 0)
}

test "hmisc_varclus with two variables" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_varclus(data, names=["A", "B"])
  assert_eq(result.n_clusters, 1)
  assert_eq(result.cluster_assignments.length(), 2)
  assert_eq(result.merges.length(), 0)
  assert_eq(result.heights.length(), 0)
}

test "hmisc_varclus merge steps" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
    [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0],
    [5.0, 4.0, 6.0, 3.0, 7.0, 5.0, 8.0, 4.0],
    [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0],
  ]
  let result = @src.hmisc_varclus(data, min_cluster_size=2)
  assert_true(result.merges.length() > 0)
  assert_true(result.heights.length() > 0)
  assert_eq(result.merges.length(), result.heights.length())
}

// hmisc_somers_d tests

test "hmisc_somers_d basic positive" {
  let g1 = [3.0, 5.0, 7.0, 9.0, 11.0]
  let g2 = [1.0, 2.0, 4.0, 6.0, 8.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_eq(result.n, 10)
  assert_true(result.d > 0.0)
  assert_true(result.lower <= result.d)
  assert_true(result.upper >= result.d)
}

test "hmisc_somers_d basic negative" {
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [4.0, 5.0, 6.0, 7.0, 8.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.d < 0.0)
}

test "hmisc_somers_d completely separated" {
  let g1 = [10.0, 20.0, 30.0, 40.0, 50.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.d > 0.9 && result.d < 1.01)
}

test "hmisc_somers_d equal distributions" {
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.d > -0.1 && result.d < 0.1)
}

test "hmisc_somers_d empty group1" {
  let g1 : Array[Double] = Array::new()
  let g2 = [1.0, 2.0, 3.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.d.is_nan())
  assert_true(result.lower.is_nan())
  assert_true(result.upper.is_nan())
  assert_eq(result.n, 0)
}

test "hmisc_somers_d empty group2" {
  let g1 = [1.0, 2.0, 3.0]
  let g2 : Array[Double] = Array::new()
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.d.is_nan())
  assert_eq(result.n, 0)
}

test "hmisc_somers_d confidence interval" {
  let g1 = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_true(result.upper > result.lower)
}

// hmisc_impute tests

test "hmisc_impute basic NaN imputation" {
  let data = [
    [1.0, @double.not_a_number, 3.0],
    [4.0, 5.0, @double.not_a_number],
    [@double.not_a_number, 8.0, 9.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_eq(result.length(), 3)
  assert_eq(result[0].length(), 3)
  assert_true(!result[0][1].is_nan())
  assert_true(!result[1][2].is_nan())
  assert_true(!result[2][0].is_nan())
  assert_eq(result[0][0], 1.0)
  assert_eq(result[1][1], 5.0)
}

test "hmisc_impute no NaN values" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_eq(result[0][0], 1.0)
  assert_eq(result[0][1], 2.0)
  assert_eq(result[1][2], 6.0)
}

test "hmisc_impute all NaN column" {
  let data = [
    [1.0, @double.not_a_number, 3.0],
    [4.0, @double.not_a_number, 6.0],
    [7.0, @double.not_a_number, 9.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_eq(result[0][1], 0.0)
  assert_eq(result[1][1], 0.0)
  assert_eq(result[2][1], 0.0)
}

test "hmisc_impute empty data" {
  let data : Array[Array[Double]] = Array::new()
  let result = @src.hmisc_impute(data)
  assert_eq(result.length(), 0)
}

test "hmisc_impute single row" {
  let data = [
    [1.0, @double.not_a_number, 3.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_eq(result.length(), 1)
  assert_eq(result[0][0], 1.0)
  assert_eq(result[0][1], 0.0)
  assert_eq(result[0][2], 3.0)
}

// hmisc_rcorr_summary tests

test "hmisc_rcorr_summary basic formatting" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_rcorr(data, names=["X", "Y"])
  let summary = @src.hmisc_rcorr_summary(result)
  assert_true(summary.contains("Correlation Matrix"))
  assert_true(summary.contains("Pearson"))
  assert_true(summary.contains("Variables: 2"))
  assert_true(summary.contains("X"))
  assert_true(summary.contains("Y"))
}

test "hmisc_rcorr_summary spearman formatting" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
  ]
  let result = @src.hmisc_rcorr(data, names=["A", "B"], cor_type=@src.hmisc_cor_type_spearman())
  let summary = @src.hmisc_rcorr_summary(result)
  assert_true(summary.contains("Spearman"))
}

// hmisc_describe_summary tests

test "hmisc_describe_summary formatting" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let stats = @src.hmisc_describe(data, name="Height")
  let summary = @src.hmisc_describe_summary(stats)
  assert_true(summary.contains("Variable: Height"))
  assert_true(summary.contains("n: 5"))
  assert_true(summary.contains("missing: 0"))
  assert_true(summary.contains("Mean:"))
  assert_true(summary.contains("SD:"))
  assert_true(summary.contains("Min:"))
  assert_true(summary.contains("Median:"))
  assert_true(summary.contains("Max:"))
}

test "hmisc_describe_summary with NaN" {
  let data = [1.0, @double.not_a_number, 3.0]
  let stats = @src.hmisc_describe(data, name="Test")
  let summary = @src.hmisc_describe_summary(stats)
  assert_true(summary.contains("Variable: Test"))
  assert_true(summary.contains("n: 2"))
  assert_true(summary.contains("missing: 1"))
}

// hmisc_sample_data and hmisc_sample_names tests

test "hmisc_sample_data structure" {
  let data = @src.hmisc_sample_data()
  assert_eq(data.length(), 3)
  assert_eq(data[0].length(), 10)
  assert_eq(data[1].length(), 10)
  assert_eq(data[2].length(), 10)
}

test "hmisc_sample_data values" {
  let data = @src.hmisc_sample_data()
  assert_eq(data[0][0], 1.0)
  assert_eq(data[0][9], 10.0)
  assert_eq(data[1][0], 2.1)
  assert_eq(data[1][9], 20.0)
  assert_eq(data[2][0], 5.0)
  assert_eq(data[2][9], 6.0)
}

test "hmisc_sample_names structure" {
  let names = @src.hmisc_sample_names()
  assert_eq(names.length(), 3)
  assert_eq(names[0], "X")
  assert_eq(names[1], "Y")
  assert_eq(names[2], "Z")
}

test "hmisc_sample_data_with_describe" {
  let data = @src.hmisc_sample_data()
  let stats = @src.hmisc_describe(data[0], name="X")
  assert_eq(stats.n, 10)
  assert_eq(stats.min, 1.0)
  assert_eq(stats.max, 10.0)
}

test "hmisc_sample_data_with_rcorr" {
  let data = @src.hmisc_sample_data()
  let names = @src.hmisc_sample_names()
  let result = @src.hmisc_rcorr(data, names=names)
  assert_eq(result.names.length(), 3)
  assert_eq(result.matrix.length(), 3)
  assert_true(!result.matrix[0][1].is_nan())
  assert_true(!result.matrix[0][2].is_nan())
  assert_true(!result.matrix[1][2].is_nan())
}

test "hmisc_sample_data_with_varclus" {
  let data = @src.hmisc_sample_data()
  let names = @src.hmisc_sample_names()
  let result = @src.hmisc_varclus(data, names=names)
  assert_eq(result.names.length(), 3)
  assert_true(result.n_clusters >= 1)
}

// Edge case tests

test "hmisc_pearson with NaN values" {
  let x = [1.0, 2.0, @double.not_a_number, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result.is_nan())
}

test "hmisc_spearman with NaN values" {
  let x = [1.0, 2.0, @double.not_a_number, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result.is_nan())
}

test "hmisc_describe large dataset" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
              11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
  let result = @src.hmisc_describe(data)
  assert_eq(result.n, 20)
  assert_true(result.mean > 10.0 && result.mean < 11.0)
  assert_eq(result.min, 1.0)
  assert_eq(result.max, 20.0)
}

test "hmisc_somers_d large groups" {
  let g1 = [6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0, 54.0, 60.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_eq(result.n, 20)
  assert_true(result.d > 0.9)
}

test "hmisc_impute preserves non-NaN values" {
  let data = [
    [10.0, 20.0, 30.0],
    [40.0, @double.not_a_number, 60.0],
    [70.0, 80.0, 90.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_eq(result[0][0], 10.0)
  assert_eq(result[0][1], 20.0)
  assert_eq(result[0][2], 30.0)
  assert_eq(result[1][0], 40.0)
  assert_eq(result[1][2], 60.0)
  assert_eq(result[2][0], 70.0)
  assert_eq(result[2][1], 80.0)
  assert_eq(result[2][2], 90.0)
}

test "hmisc_impute column mean computation" {
  let data = [
    [@double.not_a_number, 10.0],
    [20.0, @double.not_a_number],
    [30.0, 30.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_true(result[0][0] > 24.0 && result[0][0] < 26.0)
  assert_true(result[1][1] > 19.0 && result[1][1] < 21.0)
}

test "hmisc_rcorr empty data" {
  let data : Array[Array[Double]] = Array::new()
  let result = @src.hmisc_rcorr(data)
  assert_eq(result.names.length(), 0)
  assert_eq(result.matrix.length(), 0)
}

test "hmisc_varclus single variable" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
  ]
  let result = @src.hmisc_varclus(data, names=["X"])
  assert_eq(result.n_clusters, 1)
  assert_eq(result.cluster_assignments[0], 0)
}

test "hmisc_varclus two correlated variables" {
  let data = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
    [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
    [20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0],
  ]
  let result = @src.hmisc_varclus(data, min_cluster_size=2)
  assert_true(result.n_clusters >= 1)
}

test "hmisc_pearson_identical arrays" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

test "hmisc_pearson zero correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
  let result = @src.hmisc_pearson(x, y)
  assert_true(result > -0.5 && result < 0.5)
}

test "hmisc_spearman_identical arrays" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.hmisc_spearman(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

test "hmisc_somers_d single element groups" {
  let g1 = [5.0]
  let g2 = [3.0]
  let result = @src.hmisc_somers_d(g1, g2)
  assert_eq(result.n, 2)
  assert_true(result.d > 0.0)
}

test "hmisc_describe negative values" {
  let data = [-5.0, -3.0, -1.0, 0.0, 2.0]
  let result = @src.hmisc_describe(data)
  assert_eq(result.min, -5.0)
  assert_eq(result.max, 2.0)
  assert_true(result.mean > -1.5 && result.mean < -0.5)
}

test "hmisc_impute all NaN row" {
  let data = [
    [1.0, 2.0, 3.0],
    [@double.not_a_number, @double.not_a_number, @double.not_a_number],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.hmisc_impute(data)
  assert_true(!result[1][0].is_nan())
  assert_true(!result[1][1].is_nan())
  assert_true(!result[1][2].is_nan())
}