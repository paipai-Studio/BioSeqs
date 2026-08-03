///|
/// openCyto automated flow-cytometry gating demo.
///
/// Demonstrates template-driven gating on synthetic flow cytometry data
/// using mindensity and flowClust gating methods.

fn main() -> Unit {
  println("=== openCyto Gating Pipeline Demo ===\n")

  // Create synthetic flow cytometry data
  // Channels: FSC-A, SSC-A, FITC-A, PE-A
  let channels = ["FSC-A", "SSC-A", "FITC-A", "PE-A"]
  let sample_names = ["sample_01", "sample_02"]

  // Generate events for 2 samples
  let all_data : Array[Array[Array[Double]]] = Array::new()

  for _s in 0..<2 {
    let events : Array[Array[Double]] = Array::new()
    // 500 events: mix of negative and positive populations
    for i in 0..<500 {
      let fsc = 100.0 + (i.to_double() * 0.3)
      let ssc = 80.0 + (i.to_double() * 0.2)
      // First 250 events: FITC-negative (low FITC), PE-negative
      // Last 250 events: FITC-positive (high FITC), PE-positive
      let fitc = if i < 250 {
        10.0 + (i.to_double() % 50.0) * 0.3
      } else {
        200.0 + (i.to_double() % 50.0) * 0.5
      }
      let pe = if i < 250 {
        5.0 + (i.to_double() % 30.0) * 0.2
      } else {
        150.0 + (i.to_double() % 30.0) * 0.4
      }
      events.push([fsc, ssc, fitc, pe])
    }
    all_data.push(events)
  }

  let flowset = @bio.OcFlowSet::new(
    sample_names=sample_names,
    channel_names=channels,
    data=all_data,
  )

  println("FlowSet: " + flowset.data.length().to_string() + " samples")
  println("Channels: " + channels.length().to_string())
  println("Events per sample: " + flowset.data[0].length().to_string())
  println("")

  // Define gating template
  let args1 : Map[String, Double] = Map([], capacity=4)
  let args2 : Map[String, Double] = Map([], capacity=4)
  args2["K"] = 2.0

  let template : Array[@bio.OcGatingRule] = Array::new()
  template.push(
    @bio.OcGatingRule::new(
      parent="root",
      child="cells",
      method="mindensity",
      dims=["FSC-A"],
      args=args1,
    ),
  )
  template.push(
    @bio.OcGatingRule::new(
      parent="cells",
      child="fitc_pos",
      method="mindensity",
      dims=["FITC-A"],
      args=args1,
    ),
  )
  template.push(
    @bio.OcGatingRule::new(
      parent="fitc_pos",
      child="double_pos",
      method="flowClust",
      dims=["FITC-A", "PE-A"],
      args=args2,
    ),
  )

  println("Gating Template:")
  let sorted_template = @bio.oc_topo_sort(template)
  for rule in sorted_template {
    let mut dims_str = rule.dims[0]
    if rule.dims.length() > 1 {
      dims_str = rule.dims[0] + ", " + rule.dims[1]
    }
    println("  " + rule.parent + " -> " + rule.child + " (" + rule.method + " on " + dims_str + ")")
  }
  println("")

  // Run gating pipeline
  let results = @bio.oc_gate_flow_set(flowset, template)

  // Display results
  println("Gating Results:")
  for r in results {
    let pct = if r.parent_events > 0 {
      r.child_events.to_double() / r.parent_events.to_double() * 100.0
    } else {
      0.0
    }
    println(
      "  [" +
      r.sample +
      "] " +
      r.population +
      ": " +
      r.child_events.to_string() +
      "/" +
      r.parent_events.to_string() +
      " events (" +
      pct.to_string() +
      "%)",
    )
  }
  println("")

  // Population statistics
  println("Population Statistics:")
  let stats = @bio.oc_population_stats(results, 500)
  for s in stats {
    println(
      "  " +
      s.population +
      " (" +
      s.sample +
      "): count=" +
      s.count.to_string() +
      ", %parent=" +
      s.percent_of_parent.to_string() +
      ", %total=" +
      s.percent_of_total.to_string(),
    )
  }
  println("")

  // Demonstrate individual gating methods
  println("--- Individual Gating Method Demos ---\n")

  // mindensity on FITC channel
  let fitc_vals = @bio.oc_extract_channel(flowset.data[0], 2)
  let mindensity_gate = @bio.oc_mindensity(fitc_vals, "FITC-A", Map([], capacity=4))
  println("mindensity gate on FITC-A: [" + mindensity_gate.min1.to_string() + ", " + mindensity_gate.max1.to_string() + "]")

  // tailgate on PE channel
  let pe_vals = @bio.oc_extract_channel(flowset.data[0], 3)
  let tail_args : Map[String, Double] = Map([], capacity=4)
  tail_args["quantile"] = 0.95
  let tail_gate = @bio.oc_tailgate(pe_vals, "PE-A", tail_args)
  println("tailgate gate on PE-A (q=0.95): [" + tail_gate.min1.to_string() + ", " + tail_gate.max1.to_string() + "]")

  // flowClust 1D on FITC
  let fc_gate = @bio.oc_flow_clust(fitc_vals, [], "FITC-A", None, Map([], capacity=4), 2)
  println("flowClust 1D gate on FITC-A (K=2): [" + fc_gate.min1.to_string() + ", " + fc_gate.max1.to_string() + "]")

  // KDE demonstration
  let h = @bio.oc_bandwidth_silverman(fitc_vals)
  let (gx, gy) = @bio.oc_kde_grid(fitc_vals, 10, h)
  println("\nKDE density (10 grid points, bandwidth=" + h.to_string() + "):")
  for i in 0..<gx.length() {
    println("  mz=" + gx[i].to_string() + " density=" + gy[i].to_string())
  }

  println("\n=== Demo Complete ===")
}
