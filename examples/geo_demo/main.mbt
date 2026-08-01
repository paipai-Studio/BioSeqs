///|
/// GEO SOFT module demo - Gene Expression Omnibus format parser.
///
/// Run: moon run examples/geo_demo
fn main {
  println("=== GEO SOFT Module Demo ===")
  println("Gene Expression Omnibus Simple Omnibus Format parser")
  println("Reference: Barrett T et al. (2013) Nucleic Acids Res 41:D991-5")
  println("")

  // ------------------------------------------------------------------
  // 1. Generate and parse sample SOFT file
  // ------------------------------------------------------------------
  println("=== 1. Generate and parse sample SOFT file ===")
  let text = @bio.geo_sample_text()
  let f = @bio.geo_parse(text)
  println("Parsed " + f.n_records().to_string() + " record(s)")
  println("")

  // ------------------------------------------------------------------
  // 2. Display platform record
  // ------------------------------------------------------------------
  println("=== 2. Platform record ===")
  let platforms = f.get_by_type(@bio.geo_platform())
  println("Found " + platforms.length().to_string() + " platform(s)")
  let p = platforms[0]
  println("Accession: " + p.accession())
  println("Title: " + p.title())
  println("Data table: " + p.data_table().n_columns().to_string() + " columns, " + p.data_table().n_rows().to_string() + " rows")
  println("")

  // ------------------------------------------------------------------
  // 3. Display platform data table
  // ------------------------------------------------------------------
  println("=== 3. Platform data table ===")
  let table = p.data_table()
  let ids = table.get_column("ID")
  let vals = table.get_column("VALUE")
  for i in 0..<ids.length() {
    println("  " + ids[i] + ": " + vals[i])
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Display sample records
  // ------------------------------------------------------------------
  println("=== 4. Sample records ===")
  let samples = f.get_by_type(@bio.geo_sample())
  println("Found " + samples.length().to_string() + " sample(s)")
  for s in samples {
    println("  " + s.accession() + ": " + s.title())
    let source = s.get_attribute("Sample_source")
    match source {
      Some(v) => println("    Source: " + v)
      None => ()
    }
    let chars = s.get_attribute("Sample_characteristics")
    match chars {
      Some(v) => println("    Characteristics: " + v)
      None => ()
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Compare expression values between samples
  // ------------------------------------------------------------------
  println("=== 5. Compare expression values ===")
  let gsm1 = f.get_by_accession("GSM1").unwrap()
  let gsm2 = f.get_by_accession("GSM2").unwrap()
  let ctrl_vals = gsm1.data_table().get_column("VALUE")
  let treat_vals = gsm2.data_table().get_column("VALUE")
  let probe_ids = gsm1.data_table().get_column("ID")
  println("  Probe ID\tControl\tTreatment\tFold Change")
  for i in 0..<probe_ids.length() {
    let ctrl = ctrl_vals[i]
    let treat = treat_vals[i]
    match (@bio.parse_double(ctrl), @bio.parse_double(treat)) {
      (Some(c), Some(t)) =>
        println("  " + probe_ids[i] + "\t" + ctrl + "\t" + treat + "\t" + (t / c).to_string())
      _ => println("  " + probe_ids[i] + "\t" + ctrl + "\t" + treat + "\tN/A")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Full file summary
  // ------------------------------------------------------------------
  println("=== 6. Full file summary ===")
  println(f.to_string())

  println("=== Demo complete ===")
}
