///|
/// Bio.phenotype module tests - Phenotype microarray analysis

///|
test "phen_create_example_plate" {
  let plate = @src.phen_create_example_plate()
  assert_eq(plate.phen_n_wells(), 4)
  assert_eq(plate.phen_n_times(), 13)
  assert_eq(plate.plate_id, "ExamplePlate")
}

///|
test "phen_create_example_well" {
  let well = @src.phen_create_example_well()
  assert_eq(well.well_id, "A01")
  assert_eq(well.phen_n_points(), 11)
  assert_eq(well.plate_id, "ExamplePlate")
}

///|
test "phen_well_methods" {
  let well = @src.phen_create_example_well()
  assert_eq(well.phen_get_signal(0), 0.05)
  assert_eq(well.phen_get_time(0), 0.0)
  // Max signal is 1.0 (at the end)
  assert_true((well.phen_max_signal() - 1.0).abs() < 1.0e-9)
  // Min signal is 0.05
  assert_true((well.phen_min_signal() - 0.05).abs() < 1.0e-9)
  // Area should be positive
  assert_true(well.phen_area() > 0.0)
  // Max slope should be positive
  assert_true(well.phen_max_slope() > 0.0)
}

///|
test "phen_plate_methods" {
  let plate = @src.phen_create_example_plate()
  // Get well by id
  let a01_opt = plate.phen_get_well("A01")
  assert_true(a01_opt is Some(_))
  let a01 = a01_opt.unwrap()
  assert_eq(a01.well_id, "A01")
  // Get well at index
  let w0 = plate.phen_get_well_at(0)
  assert_eq(w0.well_id, "A01")
  // Mean signal at time 0 should be the average of all wells' initial signals
  let mean0 = plate.phen_mean_signal_at(0)
  // A01=0.05, A02=0.05, A03=0.04, B01=0.05 -> mean = 0.0475
  assert_true((mean0 - 0.0475).abs() < 1.0e-9)
  // Missing well
  let missing = plate.phen_get_well("Z99")
  assert_true(missing is None)
}

///|
test "phen_fit_logistic" {
  let well = @src.phen_create_example_well()
  let fit = @src.phen_fit_logistic(well)
  assert_eq(fit.model, "logistic")
  // Asymptote should be near 1.0 (max signal)
  assert_true(fit.asymptote > 0.9)
  assert_true(fit.asymptote < 1.1)
  // Slope should be positive
  assert_true(fit.slope > 0.0)
  // Lag time should be non-negative
  assert_true(fit.lag_time >= 0.0)
  // Initial should be near 0.05
  assert_true(fit.initial > 0.0)
  assert_true(fit.initial < 0.2)
  // R-squared should be in [0, 1]
  assert_true(fit.r_squared >= 0.0)
  assert_true(fit.r_squared <= 1.0)
}

///|
test "phen_fit_gompertz" {
  let well = @src.phen_create_example_well()
  let fit = @src.phen_fit_gompertz(well)
  assert_eq(fit.model, "gompertz")
  // Asymptote should be near 1.0
  assert_true(fit.asymptote > 0.9)
  assert_true(fit.asymptote < 1.1)
  // Slope should be positive
  assert_true(fit.slope > 0.0)
  // R-squared should be in [0, 1]
  assert_true(fit.r_squared >= 0.0)
  assert_true(fit.r_squared <= 1.0)
}

///|
test "phen_predict" {
  let well = @src.phen_create_example_well()
  let fit = @src.phen_fit_logistic(well)
  let times = [0.0, 5.0, 10.0, 20.0]
  let preds = @src.phen_predict(fit, times)
  assert_eq(preds.length(), 4)
  // Predicted values should be increasing for a growth curve
  assert_true(preds[3] >= preds[0])
  // Predicted value at time 0 should be near initial
  assert_true(preds[0] < 0.3)
  // Predicted value at time 20 should be near asymptote
  assert_true(preds[3] > 0.5)
}

///|
test "phen_control_subtraction" {
  let plate = @src.phen_create_example_plate()
  let result = @src.phen_subtract_control(plate, "A01")
  // A01 is the control, so it should not be in the result
  // Result should have 3 wells (A02, A03, B01)
  assert_eq(result.well_ids.length(), 3)
  // Times should match plate
  assert_eq(result.times.length(), plate.phen_n_times())
  // Signals should have 3 arrays
  assert_eq(result.signals.length(), 3)
  // Each signal array should have 13 points
  let mut i = 0
  while i < result.signals.length() {
    assert_eq(result.signals[i].length(), 13)
    i = i + 1
  }
  // A02 at time 0: 0.05 - 0.05 = 0.0
  let a02_idx = find_well_idx(result.well_ids, "A02")
  assert_true(a02_idx >= 0)
  if a02_idx >= 0 {
    assert_true((result.signals[a02_idx][0] - 0.0).abs() < 1.0e-9)
  }
}

///|
fn find_well_idx(well_ids : Array[String], target : String) -> Int {
  let mut idx = -1
  let mut i = 0
  while i < well_ids.length() {
    if well_ids[i] == target {
      idx = i
      i = well_ids.length()
    }
    i = i + 1
  }
  idx
}

///|
test "phen_parse_csv" {
  let csv = @src.phen_example_csv()
  let plate = @src.phen_parse_csv(csv)
  assert_eq(plate.plate_id, "P1")
  // Should have 2 wells (A01 and A02)
  assert_eq(plate.phen_n_wells(), 2)
  // A01 should have 3 time points
  let a01 = plate.phen_get_well("A01").unwrap()
  assert_eq(a01.phen_n_points(), 3)
  assert_eq(a01.well_id, "A01")
  // A02 should have 3 time points
  let a02 = plate.phen_get_well("A02").unwrap()
  assert_eq(a02.phen_n_points(), 3)
  // Check a parsed signal value
  assert_true((a02.phen_get_signal(2) - 0.2).abs() < 1.0e-9)
}

///|
test "phen_parse_json" {
  let json = "{\"plate\":\"P2\",\"times\":[0,2,4],\"wells\":[{\"id\":\"A01\",\"signals\":[0.05,0.06,0.08]},{\"id\":\"A02\",\"signals\":[0.05,0.1,0.2]}]}"
  let plate = @src.phen_parse_json(json)
  assert_eq(plate.plate_id, "P2")
  assert_eq(plate.phen_n_wells(), 2)
  assert_eq(plate.phen_n_times(), 3)
  let a01 = plate.phen_get_well("A01").unwrap()
  assert_eq(a01.phen_n_points(), 3)
  assert_true((a01.phen_get_signal(2) - 0.08).abs() < 1.0e-9)
  let a02 = plate.phen_get_well("A02").unwrap()
  assert_eq(a02.phen_n_points(), 3)
  assert_true((a02.phen_get_signal(1) - 0.1).abs() < 1.0e-9)
}

///|
test "phen_well_record_new" {
  let times = [0.0, 1.0, 2.0, 3.0]
  let signals = [0.1, 0.2, 0.3, 0.4]
  let well = @src.WellRecord::new("B05", "TestPlate", times, signals)
  assert_eq(well.well_id, "B05")
  assert_eq(well.plate_id, "TestPlate")
  assert_eq(well.phen_n_points(), 4)
  assert_eq(well.phen_get_time(2), 2.0)
  assert_eq(well.phen_get_signal(2), 0.3)
}

///|
test "phen_plate_add_well" {
  let times = [0.0, 1.0, 2.0]
  let plate = @src.PlateRecord::new("P3", times)
  assert_eq(plate.phen_n_wells(), 0)
  let well = @src.WellRecord::new("A01", "P3", times, [0.1, 0.2, 0.3])
  let plate2 = plate.phen_add_well(well)
  assert_eq(plate2.phen_n_wells(), 1)
  // Original should be unchanged (immutable)
  assert_eq(plate.phen_n_wells(), 0)
}
