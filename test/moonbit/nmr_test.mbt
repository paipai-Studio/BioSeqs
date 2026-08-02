///|
/// Tests for the Bio.NMR module (nmr.mbt).
///
/// Covers NOE distance restraints, dihedral-angle restraints, chemical
/// shifts, and NMR peak-list parsing.

// ============================================================================
// NoeRestraint construction and violation
// ============================================================================

///|
test "noe_restraint_construction" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="15ALA_HA",
    atom2="16GLY_HN",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.2),
    energy=Some(2.3),
  )
  assert_eq(r.restraint_id, 1)
  assert_eq(r.atom1, "15ALA_HA")
  assert_eq(r.atom2, "16GLY_HN")
  assert_true((r.lower_bound - 1.8).abs() < 0.001)
  assert_true((r.upper_bound - 5.0).abs() < 0.001)
  assert_true((r.observed.unwrap_or(0.0) - 3.2).abs() < 0.001)
  assert_true((r.energy.unwrap_or(0.0) - 2.3).abs() < 0.001)
}

///|
test "noe_restraint_construction_no_observed" {
  let r = @src.NoeRestraint::new(
    restraint_id=2,
    atom1="A",
    atom2="B",
    lower_bound=2.0,
    upper_bound=4.0,
    observed=None,
    energy=None,
  )
  assert_eq(r.observed.is_empty(), true)
  assert_eq(r.energy.is_empty(), true)
}

///|
test "noe_violation_within_range" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=None,
    energy=None,
  )
  assert_true((@src.noe_violation(r, 3.0) - 0.0).abs() < 0.001)
  assert_true((@src.noe_violation(r, 1.8) - 0.0).abs() < 0.001)
  assert_true((@src.noe_violation(r, 5.0) - 0.0).abs() < 0.001)
}

///|
test "noe_violation_below_range" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=None,
    energy=None,
  )
  // observed = 1.5, lower = 1.8, violation = 0.3
  assert_true((@src.noe_violation(r, 1.5) - 0.3).abs() < 0.001)
}

///|
test "noe_violation_above_range" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=None,
    energy=None,
  )
  // observed = 5.5, upper = 5.0, violation = 0.5
  assert_true((@src.noe_violation(r, 5.5) - 0.5).abs() < 0.001)
}

///|
test "noe_self_violation_with_observed" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(6.0),
    energy=None,
  )
  let v = @src.noe_self_violation(r)
  assert_true(v.is_some())
  assert_true((v.unwrap() - 1.0).abs() < 0.001)
}

///|
test "noe_self_violation_without_observed" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=None,
    energy=None,
  )
  assert_eq(@src.noe_self_violation(r).is_empty(), true)
}

///|
test "filter_violated_restraints_basic" {
  let r1 = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.0), // within
    energy=None,
  )
  let r2 = @src.NoeRestraint::new(
    restraint_id=2,
    atom1="C",
    atom2="D",
    lower_bound=1.8,
    upper_bound=3.5,
    observed=Some(4.5), // violation = 1.0
    energy=None,
  )
  let r3 = @src.NoeRestraint::new(
    restraint_id=3,
    atom1="E",
    atom2="F",
    lower_bound=1.8,
    upper_bound=6.0,
    observed=Some(2.1), // within
    energy=None,
  )
  let violated = @src.filter_violated_restraints([r1, r2, r3], 0.5)
  assert_eq(violated.length(), 1)
  assert_eq(violated[0].restraint_id, 2)
}

///|
test "filter_violated_restraints_threshold_zero" {
  // Threshold of 0.0 keeps all restraints with any positive violation.
  let r1 = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.0), // within
    energy=None,
  )
  let r2 = @src.NoeRestraint::new(
    restraint_id=2,
    atom1="C",
    atom2="D",
    lower_bound=1.8,
    upper_bound=3.5,
    observed=Some(3.6), // violation = 0.1
    energy=None,
  )
  let violated = @src.filter_violated_restraints([r1, r2], 0.0)
  assert_eq(violated.length(), 1)
  assert_eq(violated[0].restraint_id, 2)
}

// ============================================================================
// NOE restraint parsing
// ============================================================================

///|
test "parse_noe_restraints_basic" {
  let text =
    #|1 15ALA_HA 16GLY_HN 1.8 5.0 3.2 2.3
    #|2 20VAL_HN 21LYS_HN 1.8 3.5 4.5 5.1
    #|3 30LEU_HB 32ASP_HD 1.8 6.0 2.1 0.5
  let restraints = @src.parse_noe_restraints(text)
  assert_eq(restraints.length(), 3)
  assert_eq(restraints[0].restraint_id, 1)
  assert_eq(restraints[0].atom1, "15ALA_HA")
  assert_eq(restraints[0].atom2, "16GLY_HN")
  assert_true((restraints[0].lower_bound - 1.8).abs() < 0.001)
  assert_true((restraints[0].upper_bound - 5.0).abs() < 0.001)
  assert_true((restraints[0].observed.unwrap_or(0.0) - 3.2).abs() < 0.001)
}

///|
test "parse_noe_restraints_with_comments" {
  let text =
    #|# This is a comment
    #|! Another comment style
    #|1 A B 1.8 5.0
    #|# mid-file comment
    #|2 C D 1.8 3.5
  let restraints = @src.parse_noe_restraints(text)
  assert_eq(restraints.length(), 2)
  assert_eq(restraints[0].restraint_id, 1)
  assert_eq(restraints[1].restraint_id, 2)
}

///|
test "parse_noe_restraints_short_lines_skipped" {
  let text =
    #|1 A B 1.8 5.0
    #|too short
    #|2 C D 1.8 3.5
  let restraints = @src.parse_noe_restraints(text)
  assert_eq(restraints.length(), 2)
}

///|
test "parse_noe_restraints_no_observed" {
  let text = "1 A B 1.8 5.0"
  let restraints = @src.parse_noe_restraints(text)
  assert_eq(restraints.length(), 1)
  assert_eq(restraints[0].observed.is_empty(), true)
  assert_eq(restraints[0].energy.is_empty(), true)
}

///|
test "format_noe_restraint_roundtrip" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A_HA",
    atom2="B_HN",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.2),
    energy=Some(2.3),
  )
  let s = @src.format_noe_restraint(r)
  assert_true(s.contains("A_HA"))
  assert_true(s.contains("B_HN"))
}

///|
test "format_noe_restraints_table" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=None,
    energy=None,
  )
  let s = @src.format_noe_restraints([r])
  assert_true(s.contains("id"))
  assert_true(s.contains("atom1"))
  assert_true(s.contains("A"))
}

///|
test "noe_restraint_summary_basic" {
  let r1 = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.0), // within
    energy=None,
  )
  let r2 = @src.NoeRestraint::new(
    restraint_id=2,
    atom1="C",
    atom2="D",
    lower_bound=1.8,
    upper_bound=3.5,
    observed=Some(4.5), // violation = 1.0
    energy=None,
  )
  let summary = @src.noe_restraint_summary([r1, r2])
  assert_eq(summary.total, 2)
  assert_eq(summary.violated, 1)
  assert_true((summary.max_violation - 1.0).abs() < 0.001)
}

///|
test "noe_restraint_summary_no_violations" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(3.0),
    energy=None,
  )
  let summary = @src.noe_restraint_summary([r])
  assert_eq(summary.total, 1)
  assert_eq(summary.violated, 0)
  assert_true((summary.max_violation - 0.0).abs() < 0.001)
}

// ============================================================================
// Dihedral restraints
// ============================================================================

///|
test "dihedral_restraint_construction" {
  let d = @src.DihedralRestraint::new(
    restraint_id=1,
    angle_name="PHI",
    residue=15,
    lower_bound=-120.0,
    upper_bound=-60.0,
    observed=Some(-85.0),
    energy=Some(0.3),
  )
  assert_eq(d.restraint_id, 1)
  assert_eq(d.angle_name, "PHI")
  assert_eq(d.residue, 15)
  assert_true((d.lower_bound - (-120.0)).abs() < 0.001)
  assert_true((d.upper_bound - (-60.0)).abs() < 0.001)
}

///|
test "dihedral_violation_within_range" {
  let d = @src.DihedralRestraint::new(
    restraint_id=1,
    angle_name="PHI",
    residue=15,
    lower_bound=-120.0,
    upper_bound=-60.0,
    observed=None,
    energy=None,
  )
  assert_true((@src.dihedral_violation(d, -85.0) - 0.0).abs() < 0.001)
  assert_true((@src.dihedral_violation(d, -120.0) - 0.0).abs() < 0.001)
  assert_true((@src.dihedral_violation(d, -60.0) - 0.0).abs() < 0.001)
}

///|
test "dihedral_violation_above_range" {
  let d = @src.DihedralRestraint::new(
    restraint_id=1,
    angle_name="PHI",
    residue=15,
    lower_bound=-120.0,
    upper_bound=-60.0,
    observed=None,
    energy=None,
  )
  // observed = -40, upper = -60, violation = 20
  assert_true((@src.dihedral_violation(d, -40.0) - 20.0).abs() < 0.001)
}

///|
test "dihedral_violation_below_range" {
  let d = @src.DihedralRestraint::new(
    restraint_id=1,
    angle_name="PSI",
    residue=15,
    lower_bound=-50.0,
    upper_bound=10.0,
    observed=None,
    energy=None,
  )
  // observed = -70, lower = -50, violation = 20
  assert_true((@src.dihedral_violation(d, -70.0) - 20.0).abs() < 0.001)
}

///|
test "parse_dihedral_restraints_basic" {
  let text =
    #|1 PHI 15 -120.0 -60.0 -85.0 0.3
    #|2 PSI 15 -50.0 10.0 -30.0 0.2
  let restraints = @src.parse_dihedral_restraints(text)
  assert_eq(restraints.length(), 2)
  assert_eq(restraints[0].angle_name, "PHI")
  assert_eq(restraints[0].residue, 15)
  assert_true((restraints[0].lower_bound - (-120.0)).abs() < 0.001)
  assert_true((restraints[0].observed.unwrap_or(0.0) - (-85.0)).abs() < 0.001)
}

// ============================================================================
// Chemical shifts
// ============================================================================

///|
test "chemical_shift_construction" {
  let cs = @src.ChemicalShift::new(
    id=1,
    residue=15,
    atom_name="HA",
    element="H",
    shift_value=4.32,
    error=Some(0.02),
  )
  assert_eq(cs.id, 1)
  assert_eq(cs.residue, 15)
  assert_eq(cs.atom_name, "HA")
  assert_eq(cs.element, "H")
  assert_true((cs.shift_value - 4.32).abs() < 0.001)
  assert_true((cs.error.unwrap_or(0.0) - 0.02).abs() < 0.001)
}

///|
test "chemical_shift_table_lookup" {
  let table = @src.nmr_sample_chemical_shifts()
  let v = @src.lookup_shift(table, 15, "HA")
  assert_true(v.is_some())
  assert_true((v.unwrap() - 4.32).abs() < 0.001)
  let missing = @src.lookup_shift(table, 999, "X")
  assert_eq(missing.is_empty(), true)
}

///|
test "shifts_for_residue" {
  let table = @src.nmr_sample_chemical_shifts()
  let shifts = @src.shifts_for_residue(table, 15)
  // Residue 15 has HA, HN, N
  assert_eq(shifts.length(), 3)
}

///|
test "shifts_for_element" {
  let table = @src.nmr_sample_chemical_shifts()
  let h_shifts = @src.shifts_for_element(table, "H")
  // H: residue 15 HA, residue 15 HN, residue 16 HN = 3
  assert_eq(h_shifts.length(), 3)
  let c_shifts = @src.shifts_for_element(table, "C")
  assert_eq(c_shifts.length(), 1)
  let n_shifts = @src.shifts_for_element(table, "N")
  assert_eq(n_shifts.length(), 1)
}

///|
test "mean_shift_by_element" {
  let table = @src.nmr_sample_chemical_shifts()
  let mean_h = @src.mean_shift_by_element(table, "H")
  // H shifts: 4.32, 8.34, 8.21 -> mean = 6.9567
  assert_true((mean_h - 6.9567).abs() < 0.01)
}

///|
test "mean_shift_missing_element_returns_zero" {
  let table = @src.nmr_sample_chemical_shifts()
  let mean = @src.mean_shift_by_element(table, "P")
  assert_true((mean - 0.0).abs() < 0.001)
}

///|
test "parse_chemical_shifts_basic" {
  let text =
    #|1 15 HA H 4.32 0.02
    #|2 15 HN H 8.34 0.01
    #|3 15 N N 118.5 0.3
  let table = @src.parse_chemical_shifts(text)
  let v = @src.lookup_shift(table, 15, "HA")
  assert_true((v.unwrap() - 4.32).abs() < 0.001)
  let n = @src.lookup_shift(table, 15, "N")
  assert_true((n.unwrap() - 118.5).abs() < 0.001)
}

///|
test "parse_chemical_shifts_with_comments" {
  let text =
    #|# comment
    #|1 15 HA H 4.32
  let table = @src.parse_chemical_shifts(text)
  // The table has no public shifts field; use lookup.
  let v = @src.lookup_shift(table, 15, "HA")
  assert_true(v.is_some())
}

///|
test "format_chemical_shifts_table" {
  let table = @src.nmr_sample_chemical_shifts()
  let s = @src.format_chemical_shifts(table)
  assert_true(s.contains("id"))
  assert_true(s.contains("residue"))
  assert_true(s.contains("15"))
  assert_true(s.contains("HA"))
}

// ============================================================================
// Peak lists
// ============================================================================

///|
test "nmr_peak_construction" {
  let p = @src.NmrPeak::new(
    peak_id=1,
    dimensions=["15.HN", "15.N"],
    chemical_shifts=[8.34, 118.5],
    intensities=[1.0e6, 1.0e6],
    volume=Some(2.3e6),
    merity=Some(0.95),
  )
  assert_eq(p.peak_id, 1)
  assert_eq(p.dimensions.length(), 2)
  assert_eq(p.dimensions[0], "15.HN")
  assert_true((p.chemical_shifts[0] - 8.34).abs() < 0.001)
  assert_true((p.volume.unwrap() - 2.3e6).abs() < 1.0)
  assert_true((p.merity.unwrap() - 0.95).abs() < 0.001)
}

///|
test "peak_list_construction" {
  let list = @src.nmr_sample_peak_list()
  assert_eq(list.label, "15N-HSQC")
  assert_eq(list.n_dim, 2)
  assert_eq(list.peaks.length(), 3)
}

///|
test "peaks_by_residue" {
  let list = @src.nmr_sample_peak_list()
  let peaks = @src.peaks_by_residue(list, 15, 0)
  assert_eq(peaks.length(), 1)
  assert_eq(peaks[0].peak_id, 1)
}

///|
test "peaks_by_atom" {
  let list = @src.nmr_sample_peak_list()
  let peaks = @src.peaks_by_atom(list, "15.HN", 0)
  assert_eq(peaks.length(), 1)
  let none_peaks = @src.peaks_by_atom(list, "NONEXISTENT", 0)
  assert_eq(none_peaks.length(), 0)
}

///|
test "max_peak_intensity" {
  let list = @src.nmr_sample_peak_list()
  let m = @src.max_peak_intensity(list)
  // Max intensity across all peaks and dimensions: 1.2e6
  assert_true((m - 1.2e6).abs() < 1.0)
}

///|
test "max_peak_volume" {
  let list = @src.nmr_sample_peak_list()
  let m = @src.max_peak_volume(list)
  // Max volume: 2.8e6
  assert_true((m - 2.8e6).abs() < 1.0)
}

///|
test "parse_peak_list_basic" {
  let text =
    #|label 15N-HSQC
    #|1 15.HN 8.34 1.0e6 15.N 118.5 1.0e6 2.3e6 0.95
    #|2 16.HN 8.21 8.0e5 16.N 121.0 8.0e5 1.7e6 0.91
  let list = @src.parse_peak_list(text)
  assert_eq(list.label, "15N-HSQC")
  assert_eq(list.peaks.length(), 2)
  assert_eq(list.peaks[0].peak_id, 1)
  assert_eq(list.peaks[0].dimensions[0], "15.HN")
  assert_true((list.peaks[0].chemical_shifts[0] - 8.34).abs() < 0.001)
  assert_true((list.peaks[0].volume.unwrap() - 2.3e6).abs() < 1.0)
}

///|
test "parse_peak_list_with_comments" {
  let text =
    #|# comment line
    #|label test
    #|1 A 1.0 100 B 2.0 100
  let list = @src.parse_peak_list(text)
  assert_eq(list.label, "test")
  assert_eq(list.peaks.length(), 1)
}

///|
test "format_peak_list_output" {
  let list = @src.nmr_sample_peak_list()
  let s = @src.format_peak_list(list)
  assert_true(s.contains("label"))
  assert_true(s.contains("15N-HSQC"))
  assert_true(s.contains("id"))
}

// ============================================================================
// Sample data and summary
// ============================================================================

///|
test "nmr_sample_noe_restraints" {
  let restraints = @src.nmr_sample_noe_restraints()
  assert_eq(restraints.length(), 3)
  assert_eq(restraints[0].atom1, "15ALA_HA")
}

///|
test "nmr_sample_dihedral_restraints" {
  let restraints = @src.nmr_sample_dihedral_restraints()
  assert_eq(restraints.length(), 3)
  assert_eq(restraints[0].angle_name, "PHI")
}

///|
test "format_restraint_summary" {
  let r = @src.NoeRestraint::new(
    restraint_id=1,
    atom1="A",
    atom2="B",
    lower_bound=1.8,
    upper_bound=5.0,
    observed=Some(6.0),
    energy=None,
  )
  let summary = @src.noe_restraint_summary([r])
  let s = @src.format_restraint_summary(summary)
  assert_true(s.contains("NMR Restraint Summary"))
  assert_true(s.contains("Total restraints: 1"))
  assert_true(s.contains("Violated: 1"))
}
