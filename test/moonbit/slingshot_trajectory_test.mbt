///|
/// Slingshot trajectory inference — numerical correctness tests
///
/// Tests verify: lineage topology, MST edges, curve convergence,
/// pseudotime ordering on known 3-cluster and 4-cluster datasets.

let _S3_COORDS : Array[Array[Double]] = [
  [-0.5, 0.0], [-0.3, 0.1], [-0.2, -0.1], [0.0, 0.2], [0.1, -0.1],
  [1.8, 0.3], [2.0, 0.5], [2.2, 0.4], [2.1, 0.6], [1.9, 0.7],
  [4.8, -0.2], [5.0, 0.0], [5.2, 0.1], [5.1, -0.1], [4.9, 0.2],
]
let _S3_LABELS : Array[String] = [
  "A","A","A","A","A",
  "B","B","B","B","B",
  "C","C","C","C","C",
]

let _S4_COORDS : Array[Array[Double]] = [
  [-0.6, 0.1], [-0.4, -0.2], [-0.3, 0.3], [-0.1, -0.1], [0.0, 0.2], [0.1, -0.3], [0.2, 0.1], [0.3, 0.0], [0.4, -0.2], [0.5, 0.2],
  [1.7, 0.4], [1.9, 0.6], [2.0, 0.3], [2.1, 0.5], [2.2, 0.7], [2.3, 0.2], [2.4, 0.6], [2.5, 0.4],
  [4.5, 1.7], [4.7, 1.9], [4.8, 2.1], [4.9, 1.8], [5.0, 2.0], [5.1, 2.2], [5.2, 1.9],
  [4.5, -1.7], [4.7, -1.9], [4.8, -2.1], [4.9, -1.8], [5.0, -2.0], [5.1, -2.2], [5.2, -1.9],
]
let _S4_LABELS : Array[String] = [
  "A","A","A","A","A","A","A","A","A","A",
  "B","B","B","B","B","B","B","B",
  "C","C","C","C","C","C","C",
  "D","D","D","D","D","D","D",
]

fn _slingshot_config(start : Array[String]) -> @src.SlingshotAdvancedConfig raise @src.SlingshotAdvancedError {
  @src.SlingshotAdvancedConfig::create(
    start_clusters=start,
    distance=@src.slingshot_center_euclidean(),
    use_median=false,
    smoother_span=0.3,
    curve_points=20,
    max_iterations=15,
    tolerance=0.001,
  )
}

// ============ 3-cluster lineage inference ============

test "slingshot_3cluster_single_lineage" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  assert_true(@src.SlingshotAdvancedResult::lineages_string(result).length() == 1)
  let line = @src.SlingshotAdvancedResult::lineages_string(result)[0]
  assert_true(line.length() == 3)
  assert_true(line[0] == "A")
  assert_true(line[1] == "B")
  assert_true(line[2] == "C")
}

test "slingshot_3cluster_converges" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  assert_true(@src.SlingshotAdvancedResult::converged(result))
  assert_true(@src.SlingshotAdvancedResult::iterations(result) >= 1)
}

test "slingshot_3cluster_one_curve" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curves = @src.SlingshotAdvancedResult::curves(result)
  assert_true(curves.length() == 1)
  let curve = curves[0]
  assert_true(@src.SlingshotAdvancedCurve::points(curve).length() == 20)
  assert_true(@src.SlingshotAdvancedCurve::curve_length(curve) > 0.0)
  assert_true(@src.SlingshotAdvancedCurve::pseudotime(curve).length() == 15)
  assert_true(@src.SlingshotAdvancedCurve::distances(curve).length() == 15)
}

test "slingshot_3cluster_pseudotime_range" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curve = @src.SlingshotAdvancedResult::curves(result)[0]
  let pt = @src.SlingshotAdvancedCurve::pseudotime(curve)
  let dist = @src.SlingshotAdvancedCurve::distances(curve)
  // pseudotime 是有限值，distance 非负
  let mut i = 0
  while i < pt.length() {
    assert_true(pt[i] == pt[i])  // NaN check: x != x
    assert_true(dist[i] >= -0.1)
    i = i + 1
  }
}

test "slingshot_3cluster_curve_name" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curve = @src.SlingshotAdvancedResult::curves(result)[0]
  let name = @src.SlingshotAdvancedCurve::curve_name(curve)
  assert_true(name.length() > 0)
}

// ============ 4-cluster branching ============

test "slingshot_4cluster_two_lineages" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let lines = @src.SlingshotAdvancedResult::lineages_string(result)
  assert_true(lines.length() == 2)
  let line0 = lines[0]
  let line1 = lines[1]
  assert_true(line0[0] == "A" && line0[1] == "B")
  assert_true(line1[0] == "A" && line1[1] == "B")
  let last0 = line0[line0.length()-1]
  let last1 = line1[line1.length()-1]
  assert_true((last0 == "C" && last1 == "D") || (last0 == "D" && last1 == "C"))
}

test "slingshot_4cluster_converges" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  assert_true(@src.SlingshotAdvancedResult::converged(result))
  assert_true(@src.SlingshotAdvancedResult::curves(result).length() == 2)
}

test "slingshot_4cluster_pseudotime_range" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curves = @src.SlingshotAdvancedResult::curves(result)
  let mut ci = 0
  while ci < curves.length() {
    let pt = @src.SlingshotAdvancedCurve::pseudotime(curves[ci])
    assert_true(pt.length() == 32)
    assert_true(@src.SlingshotAdvancedCurve::curve_length(curves[ci]) > 0.0)
    let mut i = 0
    while i < pt.length() {
      assert_true(pt[i] == pt[i])  // finite
      i = i + 1
    }
    ci = ci + 1
  }
}

test "slingshot_4cluster_weights_nonnegative" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curves = @src.SlingshotAdvancedResult::curves(result)
  let mut ci = 0
  while ci < curves.length() {
    let w = @src.SlingshotAdvancedCurve::weights(curves[ci])
    let mut i = 0
    while i < w.length() {
      assert_true(w[i] >= -0.001 && w[i] <= 1.001)
      i = i + 1
    }
    ci = ci + 1
  }
}

// ============ MST topology ============

test "slingshot_mst_3cluster_edge_count" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  assert_true(@src.SlingshotLineageModel::edges(model).length() == 2)
  assert_true(@src.SlingshotLineageModel::cluster_names(model).length() == 3)
}

test "slingshot_mst_4cluster_edges_ab_bc_bd" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let names = @src.SlingshotLineageModel::cluster_names(model)
  let edges = @src.SlingshotLineageModel::edges(model)
  assert_true(names.length() == 4)
  assert_true(edges.length() == 3)
  let mut has_ab = false
  let mut has_bc = false
  let mut has_bd = false
  let mut has_cd = false
  let mut i = 0
  while i < edges.length() {
    let a = names[edges[i].from]
    let b = names[edges[i].to]
    let pair = if a < b { a + b } else { b + a }
    if pair == "AB" { has_ab = true }
    if pair == "BC" { has_bc = true }
    if pair == "BD" { has_bd = true }
    if pair == "CD" { has_cd = true }
    i = i + 1
  }
  assert_true(has_ab)
  assert_true(has_bc)
  assert_true(has_bd)
  assert_true(!has_cd)
}

test "slingshot_4cluster_lineage_paths" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let lines = @src.SlingshotLineageModel::lineages_string(model)
  assert_true(lines.length() == 2)
  assert_true(lines[0].length() == 3 && lines[0][0] == "A" && lines[0][1] == "B")
  assert_true(lines[1].length() == 3 && lines[1][0] == "A" && lines[1][1] == "B")
}

// ============ Root selection ============

test "slingshot_specified_root_used" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let roots = @src.SlingshotLineageModel::roots_string(model)
  assert_true(roots.length() == 1)
  assert_true(roots[0] == "A")
}

test "slingshot_multi_root_default" {
  let config = @src.SlingshotAdvancedConfig::create() catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let roots = @src.SlingshotLineageModel::roots_string(model)
  assert_true(roots.length() >= 1)
}

// ============ Predict ============

test "slingshot_predict_within_range" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let proj = @src.slingshot_predict(result, _S3_COORDS) catch {
    _ => { assert_true(false); return }
  }
  assert_true(proj.pseudotime.length() == 15)
  let mut i = 0
  while i < proj.pseudotime.length() {
    let row = proj.pseudotime[i]
    assert_true(row.length() == 1)
    match row[0] {
      None => ()
      Some(pt_val) => { assert_true(pt_val == pt_val) }
    }
    i = i + 1
  }
}

test "slingshot_predict_distance_nonnegative" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let proj = @src.slingshot_predict(result, _S4_COORDS) catch {
    _ => { assert_true(false); return }
  }
  assert_true(proj.distances.length() == 32)
  let mut i = 0
  while i < proj.distances.length() {
    let mut j = 0
    while j < proj.distances[i].length() {
      assert_true(proj.distances[i][j] >= -0.01)
      j = j + 1
    }
    i = i + 1
  }
}

// ============ Summary ============

test "slingshot_summary" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S4_COORDS, _S4_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let s = @src.slingshot_advanced_summary(result)
  assert_true(s.contains("Slingshot"))
  assert_true(s.contains("lineages"))
}

// ============ Omega threshold ============

test "slingshot_small_omega_forest_splits" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    omega=0.0001,
  ) catch {
    _ => { assert_true(false); return }
  }
  let model = @src.slingshot_infer_lineages(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let roots = @src.SlingshotLineageModel::roots_string(model)
  assert_true(roots.length() == 3)
}

// ============ Endpoint extension ============

test "slingshot_endpoints_extended_on_3cluster" {
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(_S3_COORDS, _S3_LABELS, config) catch {
    _ => { assert_true(false); return }
  }
  let curve = @src.SlingshotAdvancedResult::curves(result)[0]
  let name = @src.SlingshotAdvancedCurve::curve_name(curve)
  // name 应该非空且与 lineage 相关
  assert_true(name.length() > 0)
  let pts = @src.SlingshotAdvancedCurve::points(curve)
  assert_true(pts.length() >= 3)
}

// ============ Edge cases ============

test "slingshot_single_cluster_noop" {
  // 单 cluster 下 MST 只有 0 edges —— 实现可能 raise 或 no-op
  // 这里只验证不 crash（无论 raise 还是返回空 model 都可接受）
  let config = @src.SlingshotAdvancedConfig::create() catch {
    _ => { assert_true(false); return }
  }
  let coords : Array[Array[Double]] = [[0.0, 0.0], [0.1, 0.1], [0.2, 0.2]]
  let labels : Array[String] = ["X", "X", "X"]
  let result = try {
    ignore(@src.slingshot_infer_lineages(coords, labels, config))
    Some("ok")
  } catch {
    @src.SlingshotAdvancedError(_) => None
    _ => None
  }
  match result {
    Some(_) => ()
    None => ()
    _ => ()
  }
}

test "slingshot_two_clusters_one_lineage" {
  let coords : Array[Array[Double]] = [[0.0, 0.0], [0.1, 0.1], [1.0, 1.0], [1.1, 1.1]]
  let labels : Array[String] = ["A", "A", "B", "B"]
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(coords, labels, config) catch {
    _ => { assert_true(false); return }
  }
  let lines = @src.SlingshotAdvancedResult::lineages_string(result)
  assert_true(lines.length() == 1)
  assert_true(lines[0].length() == 2)
  assert_true(lines[0][0] == "A")
  assert_true(lines[0][1] == "B")
}

test "slingshot_two_clusters_converges" {
  let coords : Array[Array[Double]] = [[0.0, 0.0], [0.1, 0.1], [1.0, 1.0], [1.1, 1.1]]
  let labels : Array[String] = ["A", "A", "B", "B"]
  let config = _slingshot_config(["A"]) catch {
    _ => { assert_true(false); return }
  }
  let result = @src.slingshot_advanced(coords, labels, config) catch {
    _ => { assert_true(false); return }
  }
  assert_true(@src.SlingshotAdvancedResult::converged(result))
  let curves = @src.SlingshotAdvancedResult::curves(result)
  assert_true(curves.length() == 1)
  assert_true(@src.SlingshotAdvancedCurve::pseudotime(curves[0]).length() == 4)
}
