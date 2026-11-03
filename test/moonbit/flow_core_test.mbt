///|
/// Tests for the flow_core module.
/// Mirrors the Bioconductor flowCore package API.

///|
test "flow_core_sqrt_basic" {
  let result = @src.flow_core_sqrt(4.0)
  assert_true(result > 1.99 && result < 2.01)
}

///|
test "flow_core_sqrt_zero" {
  assert_eq(@src.flow_core_sqrt(0.0), 0.0)
}

///|
test "flow_core_sqrt_one" {
  assert_eq(@src.flow_core_sqrt(1.0), 1.0)
}

///|
test "flow_core_create_example_fcs" {
  let fcs = @src.flow_core_create_example_fcs()
  assert_eq(fcs.n_events, 5)
  assert_eq(fcs.n_parameters, 4)
  assert_eq(fcs.parameter_names.length(), 4)
  assert_eq(fcs.data_segment.length(), 5)
  assert_eq(fcs.parameter_names[0], "FSC-A")
  assert_eq(fcs.parameter_names[3], "PE-A")
}

///|
test "flow_core_create_example_frame" {
  let frame = @src.flow_core_create_example_frame()
  assert_eq(frame.n_events, 5)
  assert_eq(frame.n_params, 4)
  assert_eq(frame.exprs.length(), 5)
  assert_eq(frame.params.length(), 4)
  assert_eq(frame.exprs[0][0], 100.0)
  assert_eq(frame.exprs[4][3], 45.0)
}

///|
test "flow_core_create_example_set" {
  let set = @src.flow_core_create_example_set()
  assert_eq(set.frames.length(), 2)
  assert_eq(set.sample_names.length(), 2)
  assert_eq(set.sample_names[0], "sample1")
  assert_eq(set.sample_names[1], "sample2")
}

///|
test "flow_core_create_example_gate" {
  let gate = @src.flow_core_create_example_gate()
  assert_eq(gate.parameters.length(), 2)
  assert_eq(gate.coordinates.length(), 2)
  assert_eq(gate.filter_id, "example_rect_gate")
  assert_eq(gate.parameters[0], "FITC-A")
  assert_eq(gate.parameters[1], "PE-A")
}

///|
test "flow_core_create_example_compensation" {
  let comp = @src.flow_core_create_example_compensation()
  assert_eq(comp.fluorochromes.length(), 3)
  assert_eq(comp.matrix.length(), 3)
  assert_eq(comp.spillover.length(), 9)
  assert_eq(comp.fluorochromes[0], "FITC-A")
}

///|
test "flow_core_parameter_description_new" {
  let pd = @src.ParameterDescription::new("FITC-A", 1000.0, 0.0, 1000.0)
  assert_eq(pd.name, "FITC-A")
  assert_eq(pd.range, 1000.0)
  assert_eq(pd.min_range, 0.0)
  assert_eq(pd.max_range, 1000.0)
}

///|
test "flow_core_flow_frame_new" {
  let exprs = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let params = [
    @src.ParameterDescription::new("FSC-A", 1000.0, 0.0, 1000.0),
    @src.ParameterDescription::new("SSC-A", 1000.0, 0.0, 1000.0),
  ]
  let frame = @src.FlowFrame::new(exprs, params)
  assert_eq(frame.n_events, 3)
  assert_eq(frame.n_params, 2)
  assert_eq(frame.exprs.length(), 3)
  assert_eq(frame.params.length(), 2)
}

///|
test "flow_core_compensation_matrix_new" {
  let matrix = [[1.0, 0.1], [0.1, 1.0]]
  let fluorochromes = ["FITC-A", "PE-A"]
  let comp = @src.CompensationMatrix::new(matrix, fluorochromes)
  assert_eq(comp.matrix.length(), 2)
  assert_eq(comp.fluorochromes.length(), 2)
  assert_eq(comp.spillover.length(), 4)
}

///|
test "flow_core_get_markers" {
  let frame = @src.flow_core_create_example_frame()
  let markers = @src.flow_core_get_markers(frame)
  assert_eq(markers.length(), 4)
  assert_eq(markers[0], "FSC-A")
  assert_eq(markers[1], "SSC-A")
  assert_eq(markers[2], "FITC-A")
  assert_eq(markers[3], "PE-A")
}

///|
test "flow_core_log_transform_basic" {
  // log10(100) = 2 with a=1, b=10, c=0
  let result = @src.flow_core_log_transform(100.0, 1.0, 10.0, 0.0)
  assert_true(result > 1.99 && result < 2.01)
}

///|
test "flow_core_log_transform_zero_x" {
  // x <= 0 returns c
  let result = @src.flow_core_log_transform(0.0, 1.0, 10.0, 5.0)
  assert_eq(result, 5.0)
}

///|
test "flow_core_asinh_transform_zero" {
  let result = @src.flow_core_asinh_transform(0.0, 1.0, 1.0, 0.0)
  assert_eq(result, 0.0)
}

///|
test "flow_core_asinh_transform_basic" {
  // asinh(0) = 0, asinh(1) = ln(1 + sqrt(2)) ≈ 0.881
  let result = @src.flow_core_asinh_transform(1.0, 1.0, 1.0, 0.0)
  assert_true(result > 0.87 && result < 0.89)
}

///|
test "flow_core_logicle_transform_at_top" {
  // x >= t -> returns m
  let result = @src.flow_core_logicle_transform(10000.0, 1000.0, 0.5, 4.5, 0.0)
  assert_eq(result, 4.5)
}

///|
test "flow_core_logicle_transform_zero" {
  let result = @src.flow_core_logicle_transform(0.0, 1000.0, 0.5, 4.5, 0.0)
  assert_eq(result, 0.0)
}

///|
test "flow_core_compensate_basic" {
  let frame = @src.flow_core_create_example_frame()
  let comp = @src.flow_core_create_example_compensation()
  let result = @src.flow_core_compensate(frame, comp)
  assert_eq(result.n_events, frame.n_events)
  assert_eq(result.n_params, frame.n_params)
}

///|
test "flow_core_rectangle_gate_basic" {
  let frame = @src.flow_core_create_example_frame()
  let gate = @src.flow_core_create_example_gate()
  let result = @src.flow_core_rectangle_gate(frame, gate)
  assert_eq(result.n_total, 5)
  assert_true(result.n_filtered >= 0)
  assert_eq(result.percentages.length(), 2)
}

///|
test "flow_core_polygon_gate_insufficient_coords" {
  let frame = @src.flow_core_create_example_frame()
  // Example gate has only 2 coordinates; polygon requires 3+
  let gate = @src.flow_core_create_example_gate()
  let result = @src.flow_core_polygon_gate(frame, gate)
  assert_eq(result.n_total, 5)
  assert_eq(result.n_filtered, 0)
}

///|
test "flow_core_ellipsoid_gate_basic" {
  let frame = @src.flow_core_create_example_frame()
  let gate = @src.flow_core_create_example_gate()
  let result = @src.flow_core_ellipsoid_gate(frame, gate)
  assert_eq(result.n_total, 5)
  assert_eq(result.percentages.length(), 2)
}

///|
test "flow_core_quadrant_gate_basic" {
  let frame = @src.flow_core_create_example_frame()
  let gate = @src.flow_core_create_example_gate()
  let result = @src.flow_core_quadrant_gate(frame, gate)
  assert_eq(result.n_total, 5)
  assert_eq(result.percentages.length(), 4)
}

///|
test "flow_core_filter_events_dispatch" {
  let frame = @src.flow_core_create_example_frame()
  let gate = @src.flow_core_create_example_gate()
  let result = @src.flow_core_filter_events(frame, gate)
  assert_eq(result.n_total, 5)
}

///|
test "flow_core_point_in_polygon_inside" {
  let polygon = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
  assert_true(@src.flow_core_point_in_polygon(5.0, 5.0, polygon))
}

///|
test "flow_core_point_in_polygon_outside" {
  let polygon = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
  assert_true(!@src.flow_core_point_in_polygon(20.0, 20.0, polygon))
}

///|
test "flow_core_point_in_polygon_too_few_vertices" {
  let polygon = [(0.0, 0.0), (10.0, 10.0)]
  assert_true(!@src.flow_core_point_in_polygon(5.0, 5.0, polygon))
}

///|
test "flow_core_point_in_ellipse_inside" {
  assert_true(
    @src.flow_core_point_in_ellipse(0.0, 0.0, (0.0, 0.0), (10.0, 10.0)),
  )
}

///|
test "flow_core_point_in_ellipse_outside" {
  assert_true(
    !@src.flow_core_point_in_ellipse(20.0, 0.0, (0.0, 0.0), (10.0, 10.0)),
  )
}

///|
test "flow_core_compute_spillover_basic" {
  let matrix = [[1.0, 0.1], [0.1, 1.0]]
  let spillover = @src.flow_core_compute_spillover(matrix)
  assert_eq(spillover.length(), 4)
  assert_eq(spillover[0], 1.0)
  assert_eq(spillover[3], 1.0)
}

///|
test "flow_core_normalize_scale" {
  let frame = @src.flow_core_create_example_frame()
  let result = @src.flow_core_normalize(frame, "scale")
  assert_eq(result.exprs.length(), frame.exprs.length())
  // After scale normalization, values should be in [0, 1]
  assert_true(result.exprs[0][0] >= 0.0 && result.exprs[0][0] <= 1.0)
}

///|
test "flow_core_normalize_log" {
  let frame = @src.flow_core_create_example_frame()
  let result = @src.flow_core_normalize(frame, "log")
  assert_eq(result.exprs.length(), frame.exprs.length())
  // log normalization: ln(v) if v > 0, else 0
  assert_true(result.exprs[0][0] > 0.0)
}

///|
test "flow_core_merge_frames_basic" {
  let frame1 = @src.flow_core_create_example_frame()
  let frame2 = @src.flow_core_create_example_frame()
  let merged = @src.flow_core_merge_frames([frame1, frame2])
  assert_eq(merged.n_events, 10)
  assert_eq(merged.exprs.length(), 10)
}

///|
test "flow_core_summary_basic" {
  let frame = @src.flow_core_create_example_frame()
  let summary = @src.flow_core_summary(frame)
  assert_true(summary.length() > 0)
}

///|
test "flow_core_gate_type_to_string" {
  // Get a GateType from the example gate (RectangleGate).
  let gate = @src.flow_core_create_example_gate()
  let s = @src.GateType::to_string(gate.gate_type)
  assert_eq(s, "RectangleGate")
}
