///|

fn main {
  println("=== IsoformSwitchAnalyzeR Demo - Isoform Switch Analysis ===\n")
  
  let samples = ["control_1", "control_2", "control_3", "treatment_1", "treatment_2", "treatment_3"]
  let conditions = ["control", "control", "control", "treatment", "treatment", "treatment"]
  
  let mut data = @src.SwitchAnalyzeRlist::new(samples, conditions)
  
  data = data.add_gene_annotation("TP53", "Tumor Protein P53")
  data = data.add_gene_annotation("EGFR", "Epidermal Growth Factor Receptor")
  data = data.add_gene_annotation("BRCA1", "Breast Cancer 1")
  
  let tp53_iso1 = @src.IsoformExpression::new("TP53-001", "TP53", [100.0, 120.0, 110.0, 200.0, 220.0, 210.0], [50.0, 60.0, 55.0, 100.0, 110.0, 105.0], [40.0, 48.0, 44.0, 80.0, 88.0, 84.0])
  let tp53_iso2 = @src.IsoformExpression::new("TP53-002", "TP53", [100.0, 120.0, 110.0, 50.0, 55.0, 52.0], [50.0, 60.0, 55.0, 25.0, 27.5, 26.0], [40.0, 48.0, 44.0, 20.0, 22.0, 21.0])
  
  let egfr_iso1 = @src.IsoformExpression::new("EGFR-001", "EGFR", [80.0, 90.0, 85.0, 150.0, 160.0, 155.0], [40.0, 45.0, 42.5, 75.0, 80.0, 77.5], [32.0, 36.0, 34.0, 60.0, 64.0, 62.0])
  let egfr_iso2 = @src.IsoformExpression::new("EGFR-002", "EGFR", [80.0, 90.0, 85.0, 80.0, 85.0, 82.0], [40.0, 45.0, 42.5, 40.0, 42.5, 41.0], [32.0, 36.0, 34.0, 32.0, 34.0, 33.0])
  
  let brca1_iso1 = @src.IsoformExpression::new("BRCA1-001", "BRCA1", [60.0, 70.0, 65.0, 65.0, 72.0, 68.0], [30.0, 35.0, 32.5, 32.5, 36.0, 34.0], [24.0, 28.0, 26.0, 26.0, 29.0, 27.0])
  let brca1_iso2 = @src.IsoformExpression::new("BRCA1-002", "BRCA1", [60.0, 70.0, 65.0, 60.0, 65.0, 62.0], [30.0, 35.0, 32.5, 30.0, 32.5, 31.0], [24.0, 28.0, 26.0, 24.0, 26.0, 25.0])
  
  data = data.add_isoform(tp53_iso1).add_isoform(tp53_iso2)
  data = data.add_isoform(egfr_iso1).add_isoform(egfr_iso2)
  data = data.add_isoform(brca1_iso1).add_isoform(brca1_iso2)
  
  println("1. PSI Plot for TP53:")
  println("----------------------")
  println(@src.bio_isoform_plot_psi(data, "TP53"))
  
  println("\n2. PSI Plot for EGFR:")
  println("----------------------")
  println(@src.bio_isoform_plot_psi(data, "EGFR"))
  
  println("\n3. Find Isoform Switches (dPSI > 0.1, padj < 0.05):")
  println("------------------------------------------------------")
  
  let switches = @src.bio_isoform_find_switches(data, 0.1, 0.05)
  
  if switches.length() > 0 {
    println(@src.bio_isoform_switch_summary(switches))
  } else {
    println("No significant isoform switches found.")
  }
  
  println("\n4. Detailed Switch Analysis:")
  println("------------------------------")
  
  let mut i = 0
  while i < switches.length() {
    let s = switches[i]
    println("Gene: " + s.gene_name)
    println("  Switch: " + s.isoform_id_1 + " vs " + s.isoform_id_2)
    println("  dPSI: " + s.dpsi.to_string())
    println("  DIF: " + s.dif.to_string())
    println("  Adjusted p-value: " + s.adjusted_pvalue.to_string())
    println("  Direction: " + s.switch_direction)
    println("  Consequences: " + s.consequence_type.join(", "))
    println("")
    i = i + 1
  }
  
  println("\n=== IsoformSwitchAnalyzeR Demo Complete ===")
}