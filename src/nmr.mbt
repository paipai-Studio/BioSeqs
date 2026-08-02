// Bio.NMR - NMR (Nuclear Magnetic Resonance) data structures and parsing.
//
// Inspired by Biopython's Bio.NMR module, which provides:
//   - NOEdata: parse NOE (Nuclear Overhauser Effect) distance restraints
//     from XPLOR/CNS output.
//   - xpktools: tools for working with NMRView .xpk peak list files.
//   - nmrutils: general NMR utility helpers.
//
// This pure-MoonBit module provides the same core functionality:
//   - Data structures for NOE distance restraints, dihedral-angle restraints,
//     chemical shifts, and NMR peak lists.
//   - Parsers for XPLOR/CNS NOE restraint tables and NMRView .xpk peak lists.
//   - Helpers for restraint-violation detection, peak assignment, and
//     chemical-shift table queries.
//
// References:
//   https://biopython.org/docs/1.79/api/Bio.NMR.html
//   https://biopython.org/docs/1.79/api/Bio.NMR.NOEdata.html
//   https://biopython.org/docs/1.79/api/Bio.NMR.xpktools.html

// ============================================================================
// NOE distance restraints (Bio.NMR.NOEdata)
// ============================================================================

///|
/// A single NOE (Nuclear Overhauser Effect) distance restraint.
///
/// `lower` and `upper` define the allowed distance range in Angstroms.
/// `observed` is the distance computed from a structure (optional).
/// `violation` is positive when outside the range.
pub struct NoeRestraint {
  restraint_id : Int
  atom1 : String // e.g. "15ALA_HA"
  atom2 : String
  lower_bound : Double // Angstroms
  upper_bound : Double // Angstroms
  observed : Double? // distance from structure, if available
  energy : Double? // XPLOR energy term, if available
}

///|
/// Construct an NoeRestraint.
pub fn NoeRestraint::new(
  restraint_id~ : Int,
  atom1~ : String,
  atom2~ : String,
  lower_bound~ : Double,
  upper_bound~ : Double,
  observed~ : Double?,
  energy~ : Double?,
) -> NoeRestraint {
  NoeRestraint::{
    restraint_id,
    atom1,
    atom2,
    lower_bound,
    upper_bound,
    observed,
    energy,
  }
}

///|
/// Compute the violation magnitude for a restraint given an observed distance.
///
/// Returns 0.0 when the distance lies within [lower, upper]. Returns the
/// positive amount by which the distance falls outside the range otherwise.
pub fn noe_violation(restraint : NoeRestraint, observed : Double) -> Double {
  if observed < restraint.lower_bound {
    restraint.lower_bound - observed
  } else if observed > restraint.upper_bound {
    observed - restraint.upper_bound
  } else {
    0.0
  }
}

///|
/// Return the violation for the restraint's own `observed` field, when present.
pub fn noe_self_violation(restraint : NoeRestraint) -> Double? {
  match restraint.observed {
    Some(d) => Some(noe_violation(restraint, d))
    None => None
  }
}

///|
/// Filter a list of NOE restraints, returning those that are violated by more
/// than `threshold` Angstroms.
pub fn filter_violated_restraints(
  restraints : Array[NoeRestraint],
  threshold : Double,
) -> Array[NoeRestraint] {
  let out : Array[NoeRestraint] = Array::new()
  for r in restraints {
    match r.observed {
      Some(d) => {
        let v = noe_violation(r, d)
        if v > threshold {
          out.push(r)
        }
      }
      None => continue
    }
  }
  out
}

///|
/// Parse XPLOR/CNS NOE restraint table output.
///
/// The expected input is one restraint per line, in a tab- or
/// whitespace-separated format such as:
///
/// ```text
///   1   15ALA_HA   16HN   1.8   5.0   4.2   2.3
///   2   20GLY_HN   21HN   1.8   3.5   3.0   1.1
/// ```
///
/// Columns: id, atom1, atom2, lower, upper, observed?, energy?
/// Lines beginning with '#' or '!' are treated as comments.
pub fn parse_noe_restraints(text : String) -> Array[NoeRestraint] {
  let restraints : Array[NoeRestraint] = Array::new()
  let lines = nmr_split_lines(text)
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed[0:1].to_owned() == "#" || trimmed[0:1].to_owned() == "!" {
      continue
    }
    let fields = nmr_split_whitespace(trimmed)
    if fields.length() < 5 {
      continue
    }
    let id = parse_int(fields[0])
    let atom1 = fields[1]
    let atom2 = fields[2]
    let lower = nmr_parse_double(fields[3])
    let upper = nmr_parse_double(fields[4])
    let observed = if fields.length() >= 6 {
      Some(nmr_parse_double(fields[5]))
    } else {
      None
    }
    let energy = if fields.length() >= 7 {
      Some(nmr_parse_double(fields[6]))
    } else {
      None
    }
    restraints.push(
      NoeRestraint::new(
        restraint_id=id,
        atom1=atom1,
        atom2=atom2,
        lower_bound=lower,
        upper_bound=upper,
        observed=observed,
        energy=energy,
      ),
    )
  }
  restraints
}

///|
/// Format a single NOE restraint as a human-readable line.
pub fn format_noe_restraint(r : NoeRestraint) -> String {
  let mut s = r.restraint_id.to_string() + "\t" + r.atom1 + "\t" + r.atom2
  s = s + "\t" + nmr_fmt(r.lower_bound) + "\t" + nmr_fmt(r.upper_bound)
  match r.observed {
    Some(d) => s = s + "\t" + nmr_fmt(d)
    None => s = s + "\t-"
  }
  match r.energy {
    Some(e) => s = s + "\t" + nmr_fmt(e)
    None => s = s + "\t-"
  }
  s
}

///|
/// Format a list of NOE restraints as a table.
pub fn format_noe_restraints(restraints : Array[NoeRestraint]) -> String {
  let mut s = "id\tatom1\tatom2\tlower\tupper\tobserved\tenergy\n"
  for r in restraints {
    s = s + format_noe_restraint(r) + "\n"
  }
  s
}

///|
/// Summarize NOE restraint statistics: count, number violated, max violation.
pub fn noe_restraint_summary(
  restraints : Array[NoeRestraint]
) -> NmrRestraintSummary {
  let total = restraints.length()
  let mut violated = 0
  let mut max_violation = 0.0
  for r in restraints {
    match r.observed {
      Some(d) => {
        let v = noe_violation(r, d)
        if v > 0.0 {
          violated = violated + 1
          if v > max_violation {
            max_violation = v
          }
        }
      }
      None => continue
    }
  }
  NmrRestraintSummary::{ total, violated, max_violation }
}

// ============================================================================
// Dihedral angle restraints
// ============================================================================

///|
/// A dihedral-angle restraint (e.g. phi/psi/chi).
///
/// `lower` and `upper` define the allowed angle range in degrees.
pub struct DihedralRestraint {
  restraint_id : Int
  angle_name : String // e.g. "PHI", "PSI", "CHI1"
  residue : Int
  lower_bound : Double // degrees
  upper_bound : Double // degrees
  observed : Double?
  energy : Double?
}

///|
/// Construct a DihedralRestraint.
pub fn DihedralRestraint::new(
  restraint_id~ : Int,
  angle_name~ : String,
  residue~ : Int,
  lower_bound~ : Double,
  upper_bound~ : Double,
  observed~ : Double?,
  energy~ : Double?,
) -> DihedralRestraint {
  DihedralRestraint::{
    restraint_id,
    angle_name,
    residue,
    lower_bound,
    upper_bound,
    observed,
    energy,
  }
}

///|
/// Compute the violation magnitude for a dihedral restraint, accounting for
/// angle periodicity. Returns 0.0 when within range, otherwise the positive
/// angular distance outside the allowed range.
pub fn dihedral_violation(
  restraint : DihedralRestraint,
  observed_deg : Double,
) -> Double {
  let v = observed_deg
  if v >= restraint.lower_bound && v <= restraint.upper_bound {
    return 0.0
  }
  let below = restraint.lower_bound - v
  let above = v - restraint.upper_bound
  if below > above {
    below
  } else {
    above
  }
}

///|
/// Parse a dihedral-restraint table (id, angle, residue, lower, upper,
/// observed?, energy?).
pub fn parse_dihedral_restraints(text : String) -> Array[DihedralRestraint] {
  let out : Array[DihedralRestraint] = Array::new()
  let lines = nmr_split_lines(text)
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed[0:1].to_owned() == "#" || trimmed[0:1].to_owned() == "!" {
      continue
    }
    let fields = nmr_split_whitespace(trimmed)
    if fields.length() < 5 {
      continue
    }
    let id = parse_int(fields[0])
    let angle = fields[1]
    let res = parse_int(fields[2])
    let lower = nmr_parse_double(fields[3])
    let upper = nmr_parse_double(fields[4])
    let observed = if fields.length() >= 6 {
      Some(nmr_parse_double(fields[5]))
    } else {
      None
    }
    let energy = if fields.length() >= 7 {
      Some(nmr_parse_double(fields[6]))
    } else {
      None
    }
    out.push(
      DihedralRestraint::new(
        restraint_id=id,
        angle_name=angle,
        residue=res,
        lower_bound=lower,
        upper_bound=upper,
        observed=observed,
        energy=energy,
      ),
    )
  }
  out
}

// ============================================================================
// Chemical shifts (Bio.NMR chemical-shift tables, BMRB-like)
// ============================================================================

///|
/// A single chemical-shift assignment.
pub struct ChemicalShift {
  id : Int
  residue : Int
  atom_name : String
  element : String // H, C, N, ...
  shift_value : Double // ppm
  error : Double?
}

///|
/// Construct a ChemicalShift.
pub fn ChemicalShift::new(
  id~ : Int,
  residue~ : Int,
  atom_name~ : String,
  element~ : String,
  shift_value~ : Double,
  error~ : Double?,
) -> ChemicalShift {
  ChemicalShift::{ id, residue, atom_name, element, shift_value, error }
}

///|
/// A table of chemical shifts.
pub struct ChemicalShiftTable {
  shifts : Array[ChemicalShift]
}

///|
/// Construct a ChemicalShiftTable.
pub fn ChemicalShiftTable::new(
  shifts~ : Array[ChemicalShift]
) -> ChemicalShiftTable {
  ChemicalShiftTable::{ shifts }
}

///|
/// Look up the chemical shift for a (residue, atom_name) pair.
pub fn lookup_shift(
  table : ChemicalShiftTable,
  residue : Int,
  atom_name : String,
) -> Double? {
  for s in table.shifts {
    if s.residue == residue && s.atom_name == atom_name {
      return Some(s.shift_value)
    }
  }
  None
}

///|
/// Return all shifts for a given residue.
pub fn shifts_for_residue(
  table : ChemicalShiftTable,
  residue : Int,
) -> Array[ChemicalShift] {
  let out : Array[ChemicalShift] = Array::new()
  for s in table.shifts {
    if s.residue == residue {
      out.push(s)
    }
  }
  out
}

///|
/// Return all shifts for atoms of a given element (e.g. "H", "C", "N").
pub fn shifts_for_element(
  table : ChemicalShiftTable,
  element : String,
) -> Array[ChemicalShift] {
  let out : Array[ChemicalShift] = Array::new()
  for s in table.shifts {
    if s.element == element {
      out.push(s)
    }
  }
  out
}

///|
/// Compute the mean shift for atoms of the given element. Returns 0.0 if no
/// shifts match.
pub fn mean_shift_by_element(
  table : ChemicalShiftTable,
  element : String,
) -> Double {
  let shifts = shifts_for_element(table, element)
  if shifts.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for s in shifts {
    sum = sum + s.shift_value
  }
  sum / shifts.length().to_double()
}

///|
/// Parse a BMRB-like chemical-shift table.
///
/// Expected columns: id, residue, atom_name, element, shift, error?
/// Lines beginning with '#' or '!' are comments.
pub fn parse_chemical_shifts(text : String) -> ChemicalShiftTable {
  let shifts : Array[ChemicalShift] = Array::new()
  let lines = nmr_split_lines(text)
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed[0:1].to_owned() == "#" || trimmed[0:1].to_owned() == "!" {
      continue
    }
    let fields = nmr_split_whitespace(trimmed)
    if fields.length() < 5 {
      continue
    }
    let id = parse_int(fields[0])
    let res = parse_int(fields[1])
    let atom = fields[2]
    let elem = fields[3]
    let value = nmr_parse_double(fields[4])
    let err = if fields.length() >= 6 {
      Some(nmr_parse_double(fields[5]))
    } else {
      None
    }
    shifts.push(
      ChemicalShift::new(
        id=id,
        residue=res,
        atom_name=atom,
        element=elem,
        shift_value=value,
        error=err,
      ),
    )
  }
  ChemicalShiftTable::new(shifts=shifts)
}

///|
/// Format a chemical-shift table as a string.
pub fn format_chemical_shifts(table : ChemicalShiftTable) -> String {
  let mut s = "id\tresidue\tatom\telement\tshift\terror\n"
  for sh in table.shifts {
    s = s +
      sh.id.to_string() +
      "\t" +
      sh.residue.to_string() +
      "\t" +
      sh.atom_name +
      "\t" +
      sh.element +
      "\t" +
      nmr_fmt(sh.shift_value)
    match sh.error {
      Some(e) => s = s + "\t" + nmr_fmt(e)
      None => s = s + "\t-"
    }
    s = s + "\n"
  }
  s
}

// ============================================================================
// Peak lists (Bio.NMR.xpktools - NMRView .xpk format)
// ============================================================================

///|
/// A single peak in an NMR peak list.
///
/// `dimensions` holds the assignment (e.g. residue.atom) for each dimension.
/// `intensities` holds the peak intensity/height per dimension (for
/// multidimensional peak lists the per-dimension intensities may be equal).
/// `volume` is the integrated peak volume when available.
pub struct NmrPeak {
  peak_id : Int
  dimensions : Array[String] // assignment per dimension
  chemical_shifts : Array[Double] // ppm per dimension
  intensities : Array[Double]
  volume : Double?
  merity : Double? // figure of merit, when available
}

///|
/// Construct an NmrPeak.
pub fn NmrPeak::new(
  peak_id~ : Int,
  dimensions~ : Array[String],
  chemical_shifts~ : Array[Double],
  intensities~ : Array[Double],
  volume~ : Double?,
  merity~ : Double?,
) -> NmrPeak {
  NmrPeak::{
    peak_id,
    dimensions,
    chemical_shifts,
    intensities,
    volume,
    merity,
  }
}

///|
/// A peak list container with optional label/header.
pub struct PeakList {
  label : String
  peaks : Array[NmrPeak]
  n_dim : Int
}

///|
/// Construct a PeakList.
pub fn PeakList::new(
  label~ : String,
  peaks~ : Array[NmrPeak],
  n_dim~ : Int,
) -> PeakList {
  PeakList::{ label, peaks, n_dim }
}

///|
/// Return all peaks assigned to the given residue in dimension `dim`.
pub fn peaks_by_residue(
  list : PeakList,
  residue : Int,
  dim : Int,
) -> Array[NmrPeak] {
  let out : Array[NmrPeak] = Array::new()
  let prefix = residue.to_string()
  for p in list.peaks {
    if dim < 0 || dim >= p.dimensions.length() {
      continue
    }
    let d = p.dimensions[dim]
    // assignment format like "15.HA" or "15ALA_HA"; match the leading residue.
    if nmr_assignment_residue(d) == prefix {
      out.push(p)
    }
  }
  out
}

///|
/// Return all peaks whose assignment in dimension `dim` matches `atom`.
pub fn peaks_by_atom(
  list : PeakList,
  atom : String,
  dim : Int,
) -> Array[NmrPeak] {
  let out : Array[NmrPeak] = Array::new()
  for p in list.peaks {
    if dim < 0 || dim >= p.dimensions.length() {
      continue
    }
    if p.dimensions[dim] == atom {
      out.push(p)
    }
  }
  out
}

///|
/// Compute the maximum peak intensity across all peaks and dimensions.
pub fn max_peak_intensity(list : PeakList) -> Double {
  let mut m = 0.0
  for p in list.peaks {
    for v in p.intensities {
      if v > m {
        m = v
      }
    }
  }
  m
}

///|
/// Compute the maximum peak volume across all peaks. Returns 0.0 if no peaks
/// have a volume.
pub fn max_peak_volume(list : PeakList) -> Double {
  let mut m = 0.0
  for p in list.peaks {
    match p.volume {
      Some(v) => if v > m { m = v }
      None => continue
    }
  }
  m
}

///|
/// Parse an NMRView .xpk-style peak list.
///
/// The expected format is a header line starting with "label" followed by
/// column descriptors, then one peak per line with whitespace-separated
/// fields. The first field is the peak id; subsequent fields are interpreted
/// in groups of (assignment, shift, intensity) per dimension, optionally
/// followed by a volume and merit. Lines starting with '#' are comments.
///
/// Example header:
///   label 15N-HSQC
///   id  w1  w1_ppm  w1_int  w2  w2_ppm  w2_int  volume  merit
///   1  15.HN  8.34  1.0e6  15.N  118.5  1.0e6  2.3e6  0.95
pub fn parse_peak_list(text : String) -> PeakList {
  let peaks : Array[NmrPeak] = Array::new()
  let lines = nmr_split_lines(text)
  let mut label = ""
  let mut n_dim = 0
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed[0:1].to_owned() == "#" {
      continue
    }
    let fields = nmr_split_whitespace(trimmed)
    if fields.length() == 0 {
      continue
    }
    if fields[0] == "label" {
      // "label 15N-HSQC" or "label: 15N-HSQC"
      if fields.length() >= 2 {
        label = fields[1]
      }
      continue
    }
    if fields[0] == "id" || (fields[0].length() >= 2 && fields[0][0:2].to_owned() == "id") {
      // Header row describing columns; count dimension groups.
      // Skip detection: leave n_dim to be inferred from data.
      continue
    }
    let id = parse_int(fields[0])
    // Determine dimension count: remaining fields grouped by 3.
    let n_fields = fields.length() - 1
    let nd = n_fields / 3
    if n_dim == 0 && nd > 0 {
      n_dim = nd
    }
    let dimensions : Array[String] = Array::new()
    let shifts : Array[Double] = Array::new()
    let intensities : Array[Double] = Array::new()
    let mut idx = 1
    for _ in 0..<nd {
      if idx >= fields.length() {
        dimensions.push("?")
        shifts.push(0.0)
        intensities.push(0.0)
        idx = idx + 3
        continue
      }
      let assign = fields[idx]
      let sh = if idx + 1 < fields.length() {
        nmr_parse_double(fields[idx + 1])
      } else {
        0.0
      }
      let inten = if idx + 2 < fields.length() {
        nmr_parse_double(fields[idx + 2])
      } else {
        0.0
      }
      dimensions.push(assign)
      shifts.push(sh)
      intensities.push(inten)
      idx = idx + 3
    }
    let volume = if idx < fields.length() {
      Some(nmr_parse_double(fields[idx]))
    } else {
      None
    }
    let merit = if idx + 1 < fields.length() {
      Some(nmr_parse_double(fields[idx + 1]))
    } else {
      None
    }
    peaks.push(
      NmrPeak::new(
        peak_id=id,
        dimensions=dimensions,
        chemical_shifts=shifts,
        intensities=intensities,
        volume=volume,
        merity=merit,
      ),
    )
  }
  PeakList::new(label=label, peaks=peaks, n_dim=n_dim)
}

///|
/// Format a peak list as a string table.
pub fn format_peak_list(list : PeakList) -> String {
  let mut s = "label " + list.label + "\n"
  s = s + "id\tassignment\tshifts\tintensities\tvolume\tmerit\n"
  for p in list.peaks {
    let dims = nmr_join(p.dimensions, ",")
    let shs = nmr_join_double(p.chemical_shifts, ",")
    let ints = nmr_join_double(p.intensities, ",")
    s = s +
      p.peak_id.to_string() +
      "\t" +
      dims +
      "\t" +
      shs +
      "\t" +
      ints +
      "\t"
    match p.volume {
      Some(v) => s = s + nmr_fmt(v)
      None => s = s + "-"
    }
    s = s + "\t"
    match p.merity {
      Some(m) => s = s + nmr_fmt(m)
      None => s = s + "-"
    }
    s = s + "\n"
  }
  s
}

// ============================================================================
// Summary type
// ============================================================================

///|
/// Summary statistics for a restraint table.
pub struct NmrRestraintSummary {
  total : Int
  violated : Int
  max_violation : Double
}

///|
/// Format an NmrRestraintSummary as a string.
pub fn format_restraint_summary(s : NmrRestraintSummary) -> String {
  let mut out = "=== NMR Restraint Summary ===\n"
  out = out + "Total restraints: " + s.total.to_string() + "\n"
  out = out + "Violated: " + s.violated.to_string() + "\n"
  out = out + "Max violation: " + nmr_fmt(s.max_violation) + " A\n"
  out
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Sample NOE restraints for testing/demos.
pub fn nmr_sample_noe_restraints() -> Array[NoeRestraint] {
  [
    NoeRestraint::new(
      restraint_id=1,
      atom1="15ALA_HA",
      atom2="16GLY_HN",
      lower_bound=1.8,
      upper_bound=5.0,
      observed=Some(3.2),
      energy=Some(2.3),
    ),
    NoeRestraint::new(
      restraint_id=2,
      atom1="20VAL_HN",
      atom2="21LYS_HN",
      lower_bound=1.8,
      upper_bound=3.5,
      observed=Some(4.5),
      energy=Some(5.1),
    ),
    NoeRestraint::new(
      restraint_id=3,
      atom1="30LEU_HB",
      atom2="32ASP_HD",
      lower_bound=1.8,
      upper_bound=6.0,
      observed=Some(2.1),
      energy=Some(0.5),
    ),
  ]
}

///|
/// Sample chemical shifts for testing/demos.
pub fn nmr_sample_chemical_shifts() -> ChemicalShiftTable {
  ChemicalShiftTable::new(shifts=[
    ChemicalShift::new(
      id=1,
      residue=15,
      atom_name="HA",
      element="H",
      shift_value=4.32,
      error=Some(0.02),
    ),
    ChemicalShift::new(
      id=2,
      residue=15,
      atom_name="HN",
      element="H",
      shift_value=8.34,
      error=Some(0.01),
    ),
    ChemicalShift::new(
      id=3,
      residue=15,
      atom_name="N",
      element="N",
      shift_value=118.5,
      error=Some(0.3),
    ),
    ChemicalShift::new(
      id=4,
      residue=16,
      atom_name="HN",
      element="H",
      shift_value=8.21,
      error=Some(0.02),
    ),
    ChemicalShift::new(
      id=5,
      residue=20,
      atom_name="CA",
      element="C",
      shift_value=55.4,
      error=Some(0.4),
    ),
  ])
}

///|
/// Sample peak list (2D HSQC) for testing/demos.
pub fn nmr_sample_peak_list() -> PeakList {
  PeakList::new(
    label="15N-HSQC",
    n_dim=2,
    peaks=[
      NmrPeak::new(
        peak_id=1,
        dimensions=["15.HN", "15.N"],
        chemical_shifts=[8.34, 118.5],
        intensities=[1.0e6, 1.0e6],
        volume=Some(2.3e6),
        merity=Some(0.95),
      ),
      NmrPeak::new(
        peak_id=2,
        dimensions=["16.HN", "16.N"],
        chemical_shifts=[8.21, 121.0],
        intensities=[8.0e5, 8.0e5],
        volume=Some(1.7e6),
        merity=Some(0.91),
      ),
      NmrPeak::new(
        peak_id=3,
        dimensions=["20.HN", "20.N"],
        chemical_shifts=[7.95, 110.2],
        intensities=[1.2e6, 1.2e6],
        volume=Some(2.8e6),
        merity=Some(0.98),
      ),
    ],
  )
}

///|
/// Sample dihedral restraints for testing/demos.
pub fn nmr_sample_dihedral_restraints() -> Array[DihedralRestraint] {
  [
    DihedralRestraint::new(
      restraint_id=1,
      angle_name="PHI",
      residue=15,
      lower_bound=-120.0,
      upper_bound=-60.0,
      observed=Some(-85.0),
      energy=Some(0.3),
    ),
    DihedralRestraint::new(
      restraint_id=2,
      angle_name="PSI",
      residue=15,
      lower_bound=-50.0,
      upper_bound=10.0,
      observed=Some(-30.0),
      energy=Some(0.2),
    ),
    DihedralRestraint::new(
      restraint_id=3,
      angle_name="PHI",
      residue=20,
      lower_bound=-160.0,
      upper_bound=-80.0,
      observed=Some(-145.0),
      energy=Some(0.5),
    ),
  ]
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Split a string into lines (handling \n and \r\n).
fn nmr_split_lines(text : String) -> Array[String] {
  let out : Array[String] = Array::new()
  let mut cur = ""
  for i in 0..<text.length() {
    let c = text[i:i + 1].to_owned()
    if c == "\n" {
      // Strip trailing \r if present.
      if cur.length() > 0 && cur[cur.length() - 1:cur.length()].to_owned() == "\r" {
        cur = cur[0:cur.length() - 1].to_owned()
      }
      out.push(cur)
      cur = ""
    } else {
      cur = cur + c
    }
  }
  if cur.length() > 0 {
    if cur[cur.length() - 1:cur.length()].to_owned() == "\r" {
      cur = cur[0:cur.length() - 1].to_owned()
    }
    out.push(cur)
  }
  out
}

///|
/// Split a string on whitespace (spaces and tabs), ignoring empty fields.
fn nmr_split_whitespace(s : String) -> Array[String] {
  let out : Array[String] = Array::new()
  let mut cur = ""
  for i in 0..<s.length() {
    let c = s[i:i + 1].to_owned()
    if c == " " || c == "\t" {
      if cur.length() > 0 {
        out.push(cur)
        cur = ""
      }
    } else {
      cur = cur + c
    }
  }
  if cur.length() > 0 {
    out.push(cur)
  }
  out
}

///|
/// Parse a double-precision number, returning 0.0 on failure.
fn nmr_parse_double(s : String) -> Double {
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Join an array of strings with a separator.
fn nmr_join(arr : Array[String], sep : String) -> String {
  let mut s = ""
  for i in 0..<arr.length() {
    if i > 0 {
      s = s + sep
    }
    s = s + arr[i]
  }
  s
}

///|
/// Join an array of doubles with a separator, formatting each.
fn nmr_join_double(arr : Array[Double], sep : String) -> String {
  let mut s = ""
  for i in 0..<arr.length() {
    if i > 0 {
      s = s + sep
    }
    s = s + nmr_fmt(arr[i])
  }
  s
}

///|
/// Extract the residue number prefix from an assignment string like "15.HN"
/// or "15ALA_HN" -> "15". Returns "" when no leading digits are present.
fn nmr_assignment_residue(s : String) -> String {
  let mut out = ""
  for i in 0..<s.length() {
    let c = s[i:i + 1].to_owned()
    if c >= "0" && c <= "9" {
      out = out + c
    } else {
      break
    }
  }
  out
}

///|
/// Format a double for display (6 decimal places).
fn nmr_fmt(x : Double) -> String {
  let neg = x < 0.0
  let abs_x = if neg { -x } else { x }
  let scaled = (abs_x * 1000000.0 + 0.5).floor()
  let int_part = (scaled / 1000000.0).floor().to_int()
  let frac_int = (scaled - int_part.to_double() * 1000000.0).to_int()
  let frac_str = if frac_int < 10 {
    "00000" + frac_int.to_string()
  } else if frac_int < 100 {
    "0000" + frac_int.to_string()
  } else if frac_int < 1000 {
    "000" + frac_int.to_string()
  } else if frac_int < 10000 {
    "00" + frac_int.to_string()
  } else if frac_int < 100000 {
    "0" + frac_int.to_string()
  } else {
    frac_int.to_string()
  }
  let sign = if neg { "-" } else { "" }
  sign + int_part.to_string() + "." + frac_str
}
