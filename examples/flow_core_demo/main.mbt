///|
/// flowCore Demo.
fn main {
  println("=== flowCore Demo ===")

  // 1. Create example FlowFrame
  println("\n1. Create example FlowFrame:")
  let frame = @src.flow_core_create_example_frame()
  println(frame.to_string())

  // 2. Demonstrate data transformations (log, asinh, logicle)
  println("\n2. Data Transformations:")
  let log_val = @src.flow_core_log_transform(100.0, 1.0, 10.0, 0.0)
  let asinh_val = @src.flow_core_asinh_transform(100.0, 1.0, 1.0, 0.0)
  let logicle_val = @src.flow_core_logicle_transform(
    100.0, 262144.0, 0.5, 4.5, 0.0,
  )
  println("  log(100) = " + log_val.to_string())
  println("  asinh(100) = " + asinh_val.to_string())
  println("  logicle(100) = " + logicle_val.to_string())

  // Apply log transform to each value in the frame manually
  println("  Apply log transform to first row of frame:")
  if frame.exprs.length() > 0 {
    let mut j = 0
    while j < frame.exprs[0].length() {
      let v = frame.exprs[0][j]
      let tv = @src.flow_core_log_transform(v, 1.0, 10.0, 0.0)
      println(
        "    [" + j.to_string() + "] " + v.to_string() + " -> " + tv.to_string(),
      )
      j = j + 1
    }
  }

  // 3. Demonstrate fluorescence compensation
  println("\n3. Fluorescence Compensation:")
  let comp_matrix = @src.flow_core_create_example_compensation()
  println(comp_matrix.to_string())
  let comped_frame = @src.flow_core_compensate(frame, comp_matrix)
  println("  Compensated frame events: " + comped_frame.n_events.to_string())

  // 4. Demonstrate rectangle gate
  println("\n4. Rectangle Gate:")
  let rect_gate = @src.flow_core_create_example_gate()
  println(rect_gate.to_string())
  let rect_result = @src.flow_core_rectangle_gate(frame, rect_gate)
  println(rect_result.to_string())

  // 5. Demonstrate polygon gate
  println("\n5. Polygon Gate:")
  // Use Gate::new with gate_type extracted from the example gate (enum variants are read-only externally).
  let poly_gate = @src.Gate::new(
    rect_gate.gate_type,
    ["FITC-A", "PE-A"],
    [(0.0, 0.0), (50.0, 0.0), (50.0, 50.0), (0.0, 50.0)],
    "poly_gate_1",
  )
  println(poly_gate.to_string())
  let poly_result = @src.flow_core_polygon_gate(frame, poly_gate)
  println(poly_result.to_string())

  // 6. Demonstrate ellipsoid gate
  println("\n6. Ellipsoid Gate:")
  let ellipse_gate = @src.Gate::new(
    rect_gate.gate_type,
    ["FITC-A", "PE-A"],
    [(25.0, 25.0), (20.0, 20.0)],
    "ellipse_gate_1",
  )
  println(ellipse_gate.to_string())
  let ellipse_result = @src.flow_core_ellipsoid_gate(frame, ellipse_gate)
  println(ellipse_result.to_string())

  // 7. Demonstrate quadrant gate
  println("\n7. Quadrant Gate:")
  let quad_gate = @src.Gate::new(
    rect_gate.gate_type,
    ["FITC-A", "PE-A"],
    [(25.0, 25.0)],
    "quad_gate_1",
  )
  println(quad_gate.to_string())
  let quad_result = @src.flow_core_quadrant_gate(frame, quad_gate)
  println(quad_result.to_string())

  // 8. Demonstrate main gating entry point
  println("\n8. Main Gating Entry Point (filter_events):")
  let dispatch_result = @src.flow_core_filter_events(frame, rect_gate)
  println(dispatch_result.to_string())

  // 9. Demonstrate normalization
  println("\n9. Normalization:")
  let scaled_frame = @src.flow_core_normalize(frame, "scale")
  println("  Scaled frame events: " + scaled_frame.n_events.to_string())
  let log_norm_frame = @src.flow_core_normalize(frame, "log")
  println(
    "  Log-normalized frame events: " + log_norm_frame.n_events.to_string(),
  )

  // 10. Demonstrate merging multiple FlowFrames
  println("\n10. Merge Multiple FlowFrames:")
  let flow_set = @src.flow_core_create_example_set()
  println(flow_set.to_string())
  let merged_frame = @src.flow_core_merge_frames(flow_set.frames)
  println(
    "  Merged " +
    flow_set.frames.length().to_string() +
    " frames into one with " +
    merged_frame.n_events.to_string() +
    " events",
  )

  // 11. Demonstrate data summary
  println("\n11. FlowFrame Summary:")
  println(@src.flow_core_summary(frame))

  // 12. Bonus: example FCSFile, ParameterDescription, CompensationMatrix, and markers
  println(
    "\n12. Bonus: FCSFile, ParameterDescription, CompensationMatrix, markers:",
  )
  let fcs = @src.flow_core_create_example_fcs()
  println(fcs.to_string())
  let markers = @src.flow_core_get_markers(frame)
  println("  Markers: " + markers.length().to_string())
  let mut i = 0
  while i < markers.length() {
    println("    " + (i + 1).to_string() + ". " + markers[i])
    i = i + 1
  }
  let param_desc = @src.ParameterDescription::new("Custom", 1000.0, 0.0, 1000.0)
  println("  " + param_desc.to_string())
  let custom_comp = @src.CompensationMatrix::new([[1.0, 0.1], [0.2, 1.0]], [
    "FITC-A", "PE-A",
  ])
  println("  " + custom_comp.to_string())

  println("\nflowCore demo completed!")
}
