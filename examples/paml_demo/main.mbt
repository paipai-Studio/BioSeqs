///|
/// PAML Demo - Demonstrates molecular evolution analysis.

fn main {
  println("=== PAML Demo ===")
  
  let alignment = @bio.paml_create_example_alignment()
  println("Alignment: \{alignment.sequences.length()} sequences")
  
  println("")
  println("Parameter estimation:")
  let result = @bio.paml_estimate_parameters(alignment)
  println("  ln(Likelihood): \{result.ln_likelihood}")
  println("  pi_A: \{result.parameters.get("pi_A").unwrap_or(0.0)}")
  println("  pi_T: \{result.parameters.get("pi_T").unwrap_or(0.0)}")
  println("  pi_C: \{result.parameters.get("pi_C").unwrap_or(0.0)}")
  println("  pi_G: \{result.parameters.get("pi_G").unwrap_or(0.0)}")
  
  println("")
  println("dN/dS calculation:")
  let dnds = @bio.paml_calculate_dnds(
    "ATGGCCGGTAAATGA",
    "ATGGCCGGTAATAGA",
    "Nei-Gojobori"
  )
  println("  dN: \{dnds.dN}")
  println("  dS: \{dnds.dS}")
  println("  omega (dN/dS): \{dnds.omega}")
  
  println("")
  println("Likelihood analysis:")
  let likelihood_result = @bio.paml_run_likelihood(alignment)
  println("  dN/dS ratios: \{likelihood_result.dnds_ratios.length()} pairs")
}