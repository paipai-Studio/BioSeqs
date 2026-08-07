// Advanced protein-sequence prediction demo.
//
// Showcases six empirical protein-analysis algorithms:
//   1. Chou–Fasman secondary-structure prediction.
//   2. IUPred intrinsically-disordered-region prediction.
//   3. COILS coiled-coil prediction.
//   4. Kolaskar–Tongaonkar antigenicity prediction.
//   5. Emini surface-accessibility prediction.
//   6. Karplus–Schulz flexibility prediction.

///|
fn paa_demo_round(value : Double, digits : Int) -> Double {
  let scale = @math.pow(10.0, digits.to_double())
  (value * scale).round() / scale
}

///|
fn paa_demo_format_ss(prediction : Array[String]) -> String {
  let mut result = ""
  for p in prediction {
    result = result + p
  }
  result
}

///|
fn main {
  println("=== Protein Analysis Advanced Demo ===")

  // -----------------------------------------------------------------------
  // 1. Chou-Fasman secondary-structure prediction
  // -----------------------------------------------------------------------
  println("\n1. Chou-Fasman secondary-structure prediction")
  let helical = @src.demo_helical_sequence()
  let cf = @src.chou_fasman_predict(helical) catch {
    ProteinAdvancedError(msg) => abort("Chou-Fasman failed: " + msg)
  }
  println("  Sequence: " + helical)
  println("  SS:       " + paa_demo_format_ss(cf.prediction()))
  println(
    "  Helices: " +
    cf.helices().length().to_string() +
    ", Sheets: " +
    cf.sheets().length().to_string() +
    ", Turns: " +
    cf.turns().length().to_string(),
  )
  for region in cf.helices() {
    let (start, end) = region
    println("  helix " + start.to_string() + "-" + end.to_string())
  }

  // -----------------------------------------------------------------------
  // 2. IUPred disorder prediction
  // -----------------------------------------------------------------------
  println("\n2. IUPred disorder prediction")
  let disorder_seq = @src.demo_disorder_sequence()
  let iu = @src.iupred(disorder_seq, window_size=25) catch {
    ProteinAdvancedError(msg) => abort("IUPred failed: " + msg)
  }
  let iu_scores = iu.scores()
  let mut max_score = 0.0
  for s in iu_scores {
    if s > max_score {
      max_score = s
    }
  }
  println("  Sequence: " + disorder_seq)
  println("  Max IUPred score: " + paa_demo_round(max_score, 4).to_string())
  let regions = iu.disordered_regions()
  if regions.length() > 0 {
    for region in regions {
      let (start, end) = region
      println(
        "  Disordered region: " +
        start.to_string() +
        "-" +
        end.to_string(),
      )
    }
  } else {
    println("  No disordered regions detected")
  }

  // -----------------------------------------------------------------------
  // 3. COILS coiled-coil prediction
  // -----------------------------------------------------------------------
  println("\n3. COILS coiled-coil prediction")
  let cc_seq = @src.demo_coiled_coil_sequence()
  let cc = @src.predict_coiled_coils(cc_seq, window=14, threshold=0.5) catch {
    ProteinAdvancedError(msg) => abort("COILS failed: " + msg)
  }
  println("  Sequence: " + cc_seq)
  let cc_scores = cc.scores()
  let mut cc_max = 0.0
  for s in cc_scores {
    if s > cc_max {
      cc_max = s
    }
  }
  println("  Max COILS score: " + paa_demo_round(cc_max, 4).to_string())
  let cc_regions = cc.coiled_coil_regions()
  if cc_regions.length() > 0 {
    for region in cc_regions {
      let (start, end) = region
      println(
        "  Coiled-coil region: " + start.to_string() + "-" + end.to_string()
      )
    }
  } else {
    println("  No coiled-coil regions detected")
  }

  // -----------------------------------------------------------------------
  // 4. Kolaskar-Tongaonkar antigenicity prediction
  // -----------------------------------------------------------------------
  println("\n4. Kolaskar-Tongaonkar antigenicity prediction")
  let ag_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
  let kt = @src.kolaskar_tongaonkar_antigenicity(ag_seq, window=7) catch {
    ProteinAdvancedError(msg) => abort("Kolaskar failed: " + msg)
  }
  println("  Sequence: " + ag_seq)
  let kt_scores = kt.scores()
  let mut kt_max = 0.0
  for s in kt_scores {
    if s > kt_max {
      kt_max = s
    }
  }
  println("  Max antigenicity score: " + paa_demo_round(kt_max, 4).to_string())
  let sites = kt.antigenic_sites()
  println("  Antigenic sites: " + sites.length().to_string())
  for site in sites {
    let (start, end) = site
    println("  site " + start.to_string() + "-" + end.to_string())
  }

  // -----------------------------------------------------------------------
  // 5. Emini surface accessibility prediction
  // -----------------------------------------------------------------------
  println("\n5. Emini surface accessibility prediction")
  let em_scores = @src.emini_surface_accessibility(ag_seq, window=6) catch {
    ProteinAdvancedError(msg) => abort("Emini failed: " + msg)
  }
  let mut em_max = 0.0
  let mut em_sum = 0.0
  for s in em_scores {
    em_sum = em_sum + s
    if s > em_max {
      em_max = s
    }
  }
  println(
    "  Mean: " +
    paa_demo_round(em_sum / em_scores.length().to_double(), 4).to_string() +
    ", Max: " +
    paa_demo_round(em_max, 4).to_string(),
  )

  // -----------------------------------------------------------------------
  // 6. Karplus-Schulz flexibility prediction
  // -----------------------------------------------------------------------
  println("\n6. Karplus-Schulz flexibility prediction")
  let ks_scores = @src.karplus_schulz_flexibility(ag_seq, window=5) catch {
    ProteinAdvancedError(msg) => abort("Karplus-Schulz failed: " + msg)
  }
  let mut ks_sum = 0.0
  let mut ks_max = 0.0
  for s in ks_scores {
    ks_sum = ks_sum + s
    if s > ks_max {
      ks_max = s
    }
  }
  println(
    "  Mean flexibility: " +
    paa_demo_round(ks_sum / ks_scores.length().to_double(), 4).to_string() +
    ", Max: " +
    paa_demo_round(ks_max, 4).to_string(),
  )

  println("\n=== Demo complete ===")
}
