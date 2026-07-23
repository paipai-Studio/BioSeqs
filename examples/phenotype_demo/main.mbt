///|
/// Bio.phenotype Demo - Phenotype microarray analysis.
/// Inspired by Biopython's Bio.phenotype package.

fn main {
  println("=== Bio.phenotype Demo ===")

  // 1. Create example plate data
  println("")
  println("1. Load example phenotype plate data:")
  let plate = @bio.phen_create_example_plate()
  println("   Plate ID: " + plate.plate_id)
  println("   Wells: " + plate.phen_n_wells().to_string())
  println("   Time points: " + plate.phen_n_times().to_string())

  // 2. Inspect individual wells
  println("")
  println("2. Inspect wells:")
  let mut i = 0
  while i < plate.phen_n_wells() {
    let well = plate.phen_get_well_at(i)
    let max_sig = well.phen_max_signal()
    let area = well.phen_area()
    let slope = well.phen_max_slope()
    println("   " + well.well_id +
      ": points=" + well.phen_n_points().to_string() +
      ", max_signal=" + max_sig.to_string() +
      ", area=" + area.to_string() +
      ", max_slope=" + slope.to_string())
    i = i + 1
  }

  // 3. Fit growth curve (logistic model)
  println("")
  println("3. Fit logistic growth curve to well A02:")
  let a02 = plate.phen_get_well("A02").unwrap()
  let fit = @bio.phen_fit_logistic(a02)
  println("   Model: " + fit.model)
  println("   Asymptote: " + fit.asymptote.to_string())
  println("   Inflection point: " + fit.inflection.to_string())
  println("   Slope: " + fit.slope.to_string())
  println("   Lag time: " + fit.lag_time.to_string())
  println("   Initial: " + fit.initial.to_string())
  println("   R-squared: " + fit.r_squared.to_string())

  // 4. Fit Gompertz model
  println("")
  println("4. Fit Gompertz growth curve to well A02:")
  let fit_g = @bio.phen_fit_gompertz(a02)
  println("   Model: " + fit_g.model)
  println("   Asymptote: " + fit_g.asymptote.to_string())
  println("   Slope: " + fit_g.slope.to_string())
  println("   R-squared: " + fit_g.r_squared.to_string())

  // 5. Predict values
  println("")
  println("5. Predict signal values at sample time points:")
  let times = [0.0, 5.0, 10.0, 15.0, 20.0, 24.0]
  let preds = @bio.phen_predict(fit, times)
  let mut j = 0
  while j < times.length() {
    println("   t=" + times[j].to_string() + " -> predicted=" + preds[j].to_string())
    j = j + 1
  }

  // 6. Control subtraction
  println("")
  println("6. Subtract control well (A01) from sample wells:")
  let sub = @bio.phen_subtract_control(plate, "A01")
  println("   Control-subtracted wells: " + sub.well_ids.length().to_string())
  let mut k = 0
  while k < sub.well_ids.length() {
    let mut max_sub = 0.0
    let mut m = 0
    while m < sub.signals[k].length() {
      if sub.signals[k][m] > max_sub {
        max_sub = sub.signals[k][m]
      }
      m = m + 1
    }
    println("   " + sub.well_ids[k] + ": max_subtracted_signal=" + max_sub.to_string())
    k = k + 1
  }

  // 7. Parse CSV data
  println("")
  println("7. Parse CSV phenotype data:")
  let csv = @bio.phen_example_csv()
  let csv_plate = @bio.phen_parse_csv(csv)
  println("   Parsed plate: " + csv_plate.plate_id)
  println("   Wells: " + csv_plate.phen_n_wells().to_string())
  let csv_well = csv_plate.phen_get_well("A02").unwrap()
  println("   A02 signals: " + csv_well.signals[0].to_string() + ", " + csv_well.signals[1].to_string() + ", " + csv_well.signals[2].to_string())

  // 8. Parse JSON data
  println("")
  println("8. Parse JSON phenotype data:")
  let json = "{\"plate\":\"P2\",\"times\":[0,2,4],\"wells\":[{\"id\":\"A01\",\"signals\":[0.05,0.06,0.08]},{\"id\":\"A02\",\"signals\":[0.05,0.1,0.2]}]}"
  let json_plate = @bio.phen_parse_json(json)
  println("   Parsed plate: " + json_plate.plate_id)
  println("   Wells: " + json_plate.phen_n_wells().to_string())
  println("   Time points: " + json_plate.phen_n_times().to_string())

  println("")
  println("=== Bio.phenotype Demo Completed ===")
}
