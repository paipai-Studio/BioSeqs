///|

test "isoform_create_expression" {
  let isoform = @src.IsoformExpression::new(
    "iso1", "gene1", [100.0, 150.0, 120.0],
    [50.0, 75.0, 60.0], [40.0, 60.0, 48.0]
  )
  
  assert_eq(isoform.isoform_id, "iso1")
  assert_eq(isoform.gene_id, "gene1")
  assert_eq(isoform.counts.length(), 3)
  assert_eq(isoform.tpm.length(), 3)
  assert_eq(isoform.fpkm.length(), 3)
}

test "isoform_calculate_usage" {
  let isoform = @src.IsoformExpression::new(
    "iso1", "gene1", [100.0, 150.0], [50.0, 75.0], [40.0, 60.0]
  )
  
  let usage = @src.bio_isoform_calculate_usage(isoform)
  assert_eq(usage.length(), 2)
  assert_eq(usage[0], 50.0)
  assert_eq(usage[1], 75.0)
}

test "isoform_calculate_dpsi" {
  let iso1 = @src.IsoformExpression::new("iso1", "gene1", [100.0, 200.0], [50.0, 150.0], [40.0, 120.0])
  let iso2 = @src.IsoformExpression::new("iso2", "gene1", [100.0, 50.0], [50.0, 50.0], [40.0, 40.0])
  
  let dpsi = @src.bio_isoform_calculate_dpsi(iso1, iso2, [0], [1])
  
  assert_true((dpsi - 0.25).abs() < 0.01)
}

test "isoform_find_switches" {
  let data = @src.SwitchAnalyzeRlist::new(["s1", "s2", "s3", "s4"], ["control", "control", "treatment", "treatment"])
  
  let iso1 = @src.IsoformExpression::new("iso1", "gene1", [100.0, 100.0, 200.0, 200.0], [50.0, 50.0, 150.0, 150.0], [40.0, 40.0, 120.0, 120.0])
  let iso2 = @src.IsoformExpression::new("iso2", "gene1", [100.0, 100.0, 50.0, 50.0], [50.0, 50.0, 50.0, 50.0], [40.0, 40.0, 40.0, 40.0])
  
  let data2 = data.add_isoform(iso1).add_isoform(iso2).add_gene_annotation("gene1", "Gene1")
  
  let switches = @src.bio_isoform_find_switches(data2, 0.1, 0.05)
  
  assert_true(switches.length() >= 1)
}

test "isoform_switch_summary" {
  let switches = [
    @src.IsoformSwitch::new("gene1", "Gene1", "iso1", "iso2", 0.3, 0.25, 0.001, 0.005, "isoform1_up", ["significant_isoform_switch"]),
    @src.IsoformSwitch::new("gene2", "Gene2", "iso3", "iso4", -0.25, -0.2, 0.002, 0.01, "isoform2_up", ["significant_isoform_switch"]),
  ]
  
  let summary = @src.bio_isoform_switch_summary(switches)
  
  assert_true(summary.contains("Total switches"))
  assert_true(summary.contains("Gene1"))
  assert_true(summary.contains("Gene2"))
}

test "isoform_plot_psi" {
  let data = @src.SwitchAnalyzeRlist::new(["s1", "s2"], ["control", "treatment"])
  
  let iso1 = @src.IsoformExpression::new("iso1", "gene1", [100.0, 200.0], [50.0, 100.0], [40.0, 80.0])
  let iso2 = @src.IsoformExpression::new("iso2", "gene1", [100.0, 100.0], [50.0, 50.0], [40.0, 40.0])
  
  let data2 = data.add_isoform(iso1).add_isoform(iso2).add_gene_annotation("gene1", "Gene1")
  
  let plot = @src.bio_isoform_plot_psi(data2, "gene1")
  
  assert_true(plot.contains("Gene1"))
  assert_true(plot.contains("iso1"))
  assert_true(plot.contains("iso2"))
}