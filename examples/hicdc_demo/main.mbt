///|
/// hicdc module demo - Hi-C/HiChIP significant interaction detection.
///
/// Demonstrates NB background model fitting, significance testing,
/// TAD calling, and A/B compartment calling on simulated Hi-C data.
///
/// Reference: Sahin M et al. (2021) HiC-DC+. Nat Commun 12:3366.
/// Run: moon run examples/hicdc_demo
fn main {
  println("=== hicdc Module Demo ===")
  println("Hi-C/HiChIP significant interaction detection")
  println("Reference: Bioconductor HiC-DC+")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let contacts = @bio.hicdc_sample_data()
  println("Total contacts: " + contacts.length().to_string())
  let mut max_dist = 0
  for c in contacts {
    if c.distance() > max_dist {
      max_dist = c.distance()
    }
  }
  println("Max genomic distance: " + max_dist.to_string() + " bins")
  println("Chromosome: " + contacts[0].chr())
  println("")

  // ------------------------------------------------------------------
  // 2. Fit background model
  // ------------------------------------------------------------------
  println("=== 2. Fit NB Background Model ===")
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @bio.hicdc_fit_background(contacts)
  println("Intercept: " + intercept.to_string())
  println("Beta (distance): " + beta_dist.to_string())
  println("Beta (GC): " + beta_gc.to_string())
  println("Beta (mappability): " + beta_map.to_string())
  println("Dispersion: " + dispersion.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Predict expected counts
  // ------------------------------------------------------------------
  println("=== 3. Expected Counts by Distance ===")
  let test_distances = [1, 5, 10, 15, 19]
  for d in test_distances {
    let c = @bio.HiCContact::new("chr1", 0, d, 0, 1)
    let expected = @bio.hicdc_predict_expected(
      [c], intercept, beta_dist, beta_gc, beta_map
    )
    println("  Distance " + d.to_string() + ": expected=" + expected[0].to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Test significance
  // ------------------------------------------------------------------
  println("=== 4. Significance Testing ===")
  let results = @bio.hicdc_run(contacts, 0.1)
  let mut n_sig = 0
  for r in results {
    if r.significant() {
      n_sig = n_sig + 1
    }
  }
  println("Total contacts tested: " + results.length().to_string())
  println("Significant interactions (FDR < 0.1): " + n_sig.to_string())
  // Show top significant contacts.
  let mut shown = 0
  for r in results {
    if r.significant() && shown < 5 {
      println(
        "  chr1:" +
          r.bin1().to_string() +
          "-" +
          r.bin2().to_string() +
          " obs=" +
          r.observed().to_string() +
          " exp=" +
          r.expected().to_string() +
          " FDR=" +
          r.fdr().to_string(),
      )
      shown = shown + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. TAD calling
  // ------------------------------------------------------------------
  println("=== 5. TAD Boundary Calling ===")
  let boundaries = @bio.hicdc_call_tads(contacts, 20, 5, "chr1")
  println("TAD boundaries detected: " + boundaries.length().to_string())
  for b in boundaries {
    println(
      "  Boundary at bin " +
        b.bin().to_string() +
        " (DI=" +
        b.directionality_index().to_string() +
        ")",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. A/B compartment calling
  // ------------------------------------------------------------------
  println("=== 6. A/B Compartment Calling ===")
  let compartments = @bio.hicdc_call_compartments(contacts, 10, "chr1")
  let mut n_a = 0
  let mut n_b = 0
  for c in compartments {
    if c.compartment() == "A" {
      n_a = n_a + 1
    } else {
      n_b = n_b + 1
    }
  }
  println("A compartment bins: " + n_a.to_string())
  println("B compartment bins: " + n_b.to_string())
  for c in compartments {
    println(
      "  Bin " +
        c.bin().to_string() +
        ": " +
        c.compartment() +
        " (PC1=" +
        c.pc1().to_string() +
        ")",
    )
  }
  println("")

  println("=== Demo complete ===")
}
