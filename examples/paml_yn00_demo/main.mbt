///|
fn optional_value(value : Double?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "undefined"
  }
}

///|
fn yn00_demo_control(text : String) -> @bio.Yn00Control {
  @bio.yn00_parse_control(text) catch {
    Yn00Error(message) => abort("failed to parse YN00 control: " + message)
  }
}

///|
fn yn00_demo_results() -> @bio.Yn00Results {
  @bio.yn00_parse_results(@bio.yn00_example_results_text()) catch {
    Yn00Error(message) => abort("failed to parse YN00 results: " + message)
  }
}

///|
fn yn00_demo_matrix(results : @bio.Yn00Results) -> @bio.Yn00DistanceMatrix {
  results.matrix("YN00", "dS") catch {
    Yn00Error(message) => abort("failed to build YN00 matrix: " + message)
  }
}

///|
fn yn00_demo_mean(results : @bio.Yn00Results) -> Double {
  results.mean("YN00", "dS") catch {
    Yn00Error(message) => abort("failed to calculate YN00 mean: " + message)
  }
}

///|
fn main {
  println("=== Biopython Bio.Phylo.PAML.yn00 Offline Demo ===")
  let control = yn00_demo_control(@bio.yn00_example_control_text())
  let canonical = @bio.yn00_write_control(control)
  println("\n1. Strict YN00 control parsing")
  println("  sequence file: " + control.sequence_file())
  println("  output file: " + control.output_file())
  println("  genetic code: " + control.option("icode").unwrap_or("default"))
  println(
    "  canonical round-trip: " +
    (@bio.yn00_write_control(yn00_demo_control(canonical)) == canonical).to_string(),
  )

  let results = yn00_demo_results()
  println("\n2. Complete pairwise result")
  println(
    "  alignment=" +
    results.alignment_file() +
    ", sequences=" +
    results.sequence_count().to_string() +
    ", codons=" +
    results.codon_count().to_string() +
    ", pairs=" +
    results.pairs().length().to_string(),
  )
  println("  names: " + results.sequence_names().join(", "))

  let pair = results.pair("Homo_sapie", "Pan_troglo").unwrap()
  let yn = pair.yang_nielsen()
  println("\n3. Five methods for Homo_sapie vs Pan_troglo")
  println(
    "  NG86: omega=" +
    pair.ng86().omega().to_string() +
    ", dN=" +
    pair.ng86().dn().to_string() +
    ", dS=" +
    pair.ng86().ds().to_string(),
  )
  println(
    "  YN00: kappa=" +
    yn.kappa().to_string() +
    ", omega=" +
    yn.omega().to_string() +
    ", dN=" +
    yn.dn().to_string() +
    " +- " +
    yn.dn_standard_error().to_string() +
    ", dS=" +
    yn.ds().to_string() +
    " +- " +
    yn.ds_standard_error().to_string(),
  )
  println(
    "  LWL85: omega=" +
    optional_value(pair.lwl85().omega()) +
    ", dS=" +
    optional_value(pair.lwl85().ds()),
  )
  println(
    "  LWL85m: omega=" +
    optional_value(pair.lwl85_modified().omega()) +
    ", rho=" +
    optional_value(pair.lwl85_modified().rho()),
  )
  println(
    "  LPB93: omega=" +
    optional_value(pair.lpb93().omega()) +
    ", dS=" +
    optional_value(pair.lpb93().ds()),
  )

  let matrix = yn00_demo_matrix(results)
  println("\n4. Symmetric Yang-Nielsen dS matrix")
  let names = matrix.names()
  for first in names {
    let row : Array[String] = []
    for second in names {
      row.push(optional_value(matrix.value(first, second)))
    }
    println("  " + first + ": " + row.join(", "))
  }
  println("  mean pairwise dS: " + yn00_demo_mean(results).to_string())

  println("\nThe demo parses existing YN00 artifacts; it does not launch PAML.")
}
