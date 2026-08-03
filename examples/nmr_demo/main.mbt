///|
/// Bio.NMR demo: NMR restraint and chemical-shift analysis.
///
/// Demonstrates NOE distance restraints, dihedral-angle restraints,
/// chemical-shift table queries, and NMR peak-list parsing.
fn main {
  println("=== Bio.NMR Demo ===\n")

  // 1. NOE distance restraints.
  println("1. NOE Distance Restraints:")
  let restraints = @bio.nmr_sample_noe_restraints()
  for i in 0..<restraints.length() {
    let r = restraints[i]
    println(
      "   #\{r.restraint_id}: \{r.atom1} - \{r.atom2}, range [\{r.lower_bound}, \{r.upper_bound}] A, observed=\{r.observed.unwrap_or(-1.0)}",
    )
  }
  // Compute violations.
  println("\n2. Restraint Violations:")
  for i in 0..<restraints.length() {
    let r = restraints[i]
    match r.observed {
      Some(d) => {
        let v = @bio.noe_violation(r, d)
        let status = if v > 0.0 { "VIOLATED" } else { "ok" }
        println(
          "   #\{r.restraint_id}: observed=\{d} A, violation=\{v} A [\{status}]",
        )
      }
      None => continue
    }
  }
  // Summary.
  let summary = @bio.noe_restraint_summary(restraints)
  println(
    "\n   Summary: total=\{summary.total}, violated=\{summary.violated}, max_violation=\{summary.max_violation} A",
  )

  // 2. Dihedral-angle restraints.
  println("\n3. Dihedral Angle Restraints:")
  let dihedrals = @bio.nmr_sample_dihedral_restraints()
  for i in 0..<dihedrals.length() {
    let d = dihedrals[i]
    println(
      "   #\{d.restraint_id}: \{d.angle_name}(\{d.residue}), range [\{d.lower_bound}, \{d.upper_bound}] deg, observed=\{d.observed.unwrap_or(0.0)}",
    )
  }
  println("\n4. Dihedral Violations:")
  for i in 0..<dihedrals.length() {
    let d = dihedrals[i]
    match d.observed {
      Some(v) => {
        let viol = @bio.dihedral_violation(d, v)
        let status = if viol > 0.0 { "VIOLATED" } else { "ok" }
        println(
          "   #\{d.restraint_id}: observed=\{v} deg, violation=\{viol} deg [\{status}]",
        )
      }
      None => continue
    }
  }

  // 3. Chemical shifts.
  println("\n5. Chemical Shifts:")
  let table = @bio.nmr_sample_chemical_shifts()
  println("   Total shifts: \{table.shifts.length()}")
  let h_shifts = @bio.shifts_for_element(table, "H")
  println("   H shifts: \{h_shifts.length()}")
  let n_shifts = @bio.shifts_for_element(table, "N")
  println("   N shifts: \{n_shifts.length()}")
  let c_shifts = @bio.shifts_for_element(table, "C")
  println("   C shifts: \{c_shifts.length()}")
  let mean_h = @bio.mean_shift_by_element(table, "H")
  println("   Mean H shift: \{mean_h} ppm")
  // Lookup example.
  let ha_15 = @bio.lookup_shift(table, 15, "HA")
  println("   15.HA shift: \{ha_15.unwrap_or(-1.0)} ppm")

  // 4. Peak list.
  println("\n6. Peak List:")
  let list = @bio.nmr_sample_peak_list()
  println("   Label: \{list.label}")
  println("   Dimensions: \{list.n_dim}")
  println("   Peaks: \{list.peaks.length()}")
  for i in 0..<list.peaks.length() {
    let p = list.peaks[i]
    println(
      "   - #\{p.peak_id}: \{p.dimensions[0]} / \{p.dimensions[1]}, shifts=\{p.chemical_shifts[0]}, \{p.chemical_shifts[1]}, volume=\{p.volume.unwrap_or(0.0)}",
    )
  }
  let max_int = @bio.max_peak_intensity(list)
  println("   Max intensity: \{max_int}")
  let max_vol = @bio.max_peak_volume(list)
  println("   Max volume: \{max_vol}")

  // 5. Parse NOE restraints from text.
  println("\n7. Parse NOE Restraints from Text:")
  let text =
    #|1 A_HA B_HN 1.8 5.0 3.2 2.3
    #|2 C_HN D_HN 1.8 3.5 4.5 5.1
  let parsed = @bio.parse_noe_restraints(text)
  println("   Parsed \{parsed.length()} restraints")
  for i in 0..<parsed.length() {
    let r = parsed[i]
    println("   - #\{r.restraint_id}: \{r.atom1} - \{r.atom2}")
  }

  // 6. Parse a peak list from text.
  println("\n8. Parse Peak List from Text:")
  let peak_text =
    #|label demo
    #|1 15.HN 8.34 1.0e6 15.N 118.5 1.0e6 2.3e6 0.95
    #|2 16.HN 8.21 8.0e5 16.N 121.0 8.0e5 1.7e6 0.91
  let parsed_list = @bio.parse_peak_list(peak_text)
  println("   Label: \{parsed_list.label}")
  println("   Parsed \{parsed_list.peaks.length()} peaks")

  // 7. Format restraints.
  println("\n9. Format Restraints:")
  let fmt = @bio.format_noe_restraints(restraints)
  println(fmt)

  // 8. Format summary.
  println("10. Restraint Summary:")
  let summary_str = @bio.format_restraint_summary(summary)
  println(summary_str)

  println("=== Demo Complete ===")
}
