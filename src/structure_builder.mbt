///|
/// StructureBuilder — build a full peptide AtomArray from scratch
/// using internal-coordinate (φ/ψ/ω/χ) specification.
///
/// Wraps the existing BackboneAngleChain::build_backbone /
/// ChainAngleChain::build_full inverse-geometry engine with a
/// method-chaining API, default secondary-structure constants,
/// and per-residue sidechain χ overrides.
///
/// Inspired by Biopython's Bio.PDB.InternalCoordinateBuilder.

// ============================================================================
// Default secondary-structure backbone angles (degrees)
// ============================================================================

///|
/// Classical α-helix (right-handed 3.6₁₃): φ = -60°, ψ = -45°.
pub const SS_ALPHA_HELIX_PHI : Double = -60.0

///|
pub const SS_ALPHA_HELIX_PSI : Double = -45.0

///|
/// β-strand (extended): φ = -135°, ψ = 135°.
pub const SS_BETA_STRAND_PHI : Double = -135.0

///|
pub const SS_BETA_STRAND_PSI : Double = 135.0

///|
/// π-helix (5.4₁₀): φ ≈ -57°, ψ ≈ -70°.
pub const SS_PI_HELIX_PHI : Double = -57.0

///|
pub const SS_PI_HELIX_PSI : Double = -70.0

///|
/// 3₁₀ helix: φ ≈ -60°, ψ ≈ -30°.
pub const SS_3TO10_HELIX_PHI : Double = -60.0

///|
pub const SS_3TO10_HELIX_PSI : Double = -30.0

///|
/// Polyproline II helix: φ ≈ -75°, ψ ≈ 150°.
pub const SS_PPII_PHI : Double = -75.0

///|
pub const SS_PPII_PSI : Double = 150.0

///|
/// Generic "coil" default: φ = -80°, ψ = 120°.
pub const SS_COIL_PHI : Double = -80.0

///|
pub const SS_COIL_PSI : Double = 120.0

///|
/// Default peptide bond ω (trans, ~180°).
pub const DEFAULT_OMEGA : Double = 180.0

// ============================================================================
// Per-residue specification
// ============================================================================

///|
/// One residue's full internal-coordinate specification.
/// All angles in degrees; None means "use standard default".
pub struct BuilderResidue {
  name : String
  mut phi : Double?
  mut psi : Double?
  omega : Double?
  mut chi1 : Double?
  mut chi2 : Double?
  mut chi3 : Double?
  mut chi4 : Double?
}

// ============================================================================
// StructureBuilder
// ============================================================================

///|
/// Builder that accumulates residues and produces an AtomArray via
/// inverse internal-coordinate reconstruction.
pub struct StructureBuilder {
  residues : Array[BuilderResidue]
  chain_id : Char
  start_res_id : Int32
  mut default_phi : Double
  mut default_psi : Double
  default_omega : Double
  default_chi : Double // default χ angle when sidechain present
}

///|
/// Create a new empty StructureBuilder.
pub fn StructureBuilder::new(
  chain_id~ : Char = 'A',
  start_res_id~ : Int32 = 1,
) -> StructureBuilder {
  StructureBuilder::{
    residues: [],
    chain_id,
    start_res_id,
    default_phi: SS_COIL_PHI,
    default_psi: SS_COIL_PSI,
    default_omega: DEFAULT_OMEGA,
    default_chi: 180.0,
  }
}

// --- Configuration --------------------------------------------------------

///|
/// Override the default (φ, ψ) used when adding residues via
/// simple `add_residue(name)` without explicit angles.
pub fn StructureBuilder::set_default_backbone(
  self : StructureBuilder,
  phi : Double,
  psi : Double,
) -> StructureBuilder {
  self.default_phi = phi
  self.default_psi = psi
  self
}

///|
/// Set defaults to α-helix φ/ψ.
pub fn StructureBuilder::use_alpha_helix(
  self : StructureBuilder,
) -> StructureBuilder {
  self.default_phi = SS_ALPHA_HELIX_PHI
  self.default_psi = SS_ALPHA_HELIX_PSI
  self
}

///|
/// Set defaults to β-strand φ/ψ.
pub fn StructureBuilder::use_beta_strand(
  self : StructureBuilder,
) -> StructureBuilder {
  self.default_phi = SS_BETA_STRAND_PHI
  self.default_psi = SS_BETA_STRAND_PSI
  self
}

///|
/// Set defaults to π-helix φ/ψ.
pub fn StructureBuilder::use_pi_helix(
  self : StructureBuilder,
) -> StructureBuilder {
  self.default_phi = SS_PI_HELIX_PHI
  self.default_psi = SS_PI_HELIX_PSI
  self
}

///|
/// Set defaults to coil φ/ψ.
pub fn StructureBuilder::use_coil(self : StructureBuilder) -> StructureBuilder {
  self.default_phi = SS_COIL_PHI
  self.default_psi = SS_COIL_PSI
  self
}

// --- Add residues ---------------------------------------------------------

///|
/// Append one residue. Uses current default (φ, ψ) unless overridden.
pub fn StructureBuilder::add_residue(
  self : StructureBuilder,
  name : String,
  phi~ : Double = -999.0,
  psi~ : Double = -999.0,
  omega~ : Double = -999.0,
  chi1~ : Double = -999.0,
  chi2~ : Double = -999.0,
  chi3~ : Double = -999.0,
  chi4~ : Double = -999.0,
) -> StructureBuilder {
  let r : BuilderResidue = BuilderResidue::{
    name,
    phi: opt_or_default(phi, self.default_phi),
    psi: opt_or_default(psi, self.default_psi),
    omega: opt_or_default(omega, self.default_omega),
    chi1: opt_or_none(chi1),
    chi2: opt_or_none(chi2),
    chi3: opt_or_none(chi3),
    chi4: opt_or_none(chi4),
  }
  self.residues.push(r)
  self
}

///|
/// Append `n` residues of `name` with α-helix geometry.
pub fn StructureBuilder::add_helix(
  self : StructureBuilder,
  name : String,
  n~ : Int32 = 1,
) -> StructureBuilder {
  let mut i : Int32 = 0
  while i < n {
    self.residues.push(BuilderResidue::{
      name,
      phi: Some(SS_ALPHA_HELIX_PHI),
      psi: Some(SS_ALPHA_HELIX_PSI),
      omega: Some(DEFAULT_OMEGA),
      chi1: self.default_chi_or_none(name),
      chi2: None,
      chi3: None,
      chi4: None,
    })
    i = i + 1
  }
  self
}

///|
/// Append `n` residues of `name` with β-strand (extended) geometry.
pub fn StructureBuilder::add_strand(
  self : StructureBuilder,
  name : String,
  n~ : Int32 = 1,
) -> StructureBuilder {
  let mut i : Int32 = 0
  while i < n {
    self.residues.push(BuilderResidue::{
      name,
      phi: Some(SS_BETA_STRAND_PHI),
      psi: Some(SS_BETA_STRAND_PSI),
      omega: Some(DEFAULT_OMEGA),
      chi1: self.default_chi_or_none(name),
      chi2: None,
      chi3: None,
      chi4: None,
    })
    i = i + 1
  }
  self
}

///|
/// Append `n` residues of `name` with 3₁₀ helix geometry.
pub fn StructureBuilder::add_3to10_helix(
  self : StructureBuilder,
  name : String,
  n~ : Int32 = 1,
) -> StructureBuilder {
  let mut i : Int32 = 0
  while i < n {
    self.residues.push(BuilderResidue::{
      name,
      phi: Some(SS_3TO10_HELIX_PHI),
      psi: Some(SS_3TO10_HELIX_PSI),
      omega: Some(DEFAULT_OMEGA),
      chi1: self.default_chi_or_none(name),
      chi2: None,
      chi3: None,
      chi4: None,
    })
    i = i + 1
  }
  self
}

///|
/// Append `n` residues of `name` with polyproline II helix geometry.
pub fn StructureBuilder::add_ppii(
  self : StructureBuilder,
  name : String,
  n~ : Int32 = 1,
) -> StructureBuilder {
  let mut i : Int32 = 0
  while i < n {
    self.residues.push(BuilderResidue::{
      name,
      phi: Some(SS_PPII_PHI),
      psi: Some(SS_PPII_PSI),
      omega: Some(DEFAULT_OMEGA),
      chi1: self.default_chi_or_none(name),
      chi2: None,
      chi3: None,
      chi4: None,
    })
    i = i + 1
  }
  self
}

///|
/// Append `n` residues of `name` with generic coil geometry.
pub fn StructureBuilder::add_coil(
  self : StructureBuilder,
  name : String,
  n~ : Int32 = 1,
) -> StructureBuilder {
  let mut i : Int32 = 0
  while i < n {
    self.residues.push(BuilderResidue::{
      name,
      phi: Some(SS_COIL_PHI),
      psi: Some(SS_COIL_PSI),
      omega: Some(DEFAULT_OMEGA),
      chi1: self.default_chi_or_none(name),
      chi2: None,
      chi3: None,
      chi4: None,
    })
    i = i + 1
  }
  self
}

// --- Per-residue edits ----------------------------------------------------

///|
/// Override φ for residue index `i`.
pub fn StructureBuilder::set_phi(
  self : StructureBuilder,
  i : Int32,
  phi : Double,
) -> StructureBuilder {
  if i >= 0 && i < self.residues.length() {
    self.residues[i].phi = Some(phi)
  }
  self
}

///|
/// Override ψ for residue index `i`.
pub fn StructureBuilder::set_psi(
  self : StructureBuilder,
  i : Int32,
  psi : Double,
) -> StructureBuilder {
  if i >= 0 && i < self.residues.length() {
    self.residues[i].psi = Some(psi)
  }
  self
}

///|
/// Override sidechain χ for residue index `i`.
/// chi_idx is 1-based (1=χ¹, 2=χ², 3=χ³, 4=χ⁴).
pub fn StructureBuilder::set_chi(
  self : StructureBuilder,
  i : Int32,
  chi_idx : Int32,
  angle : Double,
) -> StructureBuilder {
  if i >= 0 && i < self.residues.length() {
    let r = self.residues[i]
    match chi_idx {
      1 => r.chi1 = Some(angle)
      2 => r.chi2 = Some(angle)
      3 => r.chi3 = Some(angle)
      4 => r.chi4 = Some(angle)
      _ => ()
    }
  }
  self
}

///|
/// Number of residues currently accumulated.
pub fn StructureBuilder::len(self : StructureBuilder) -> Int32 {
  self.residues.length()
}

// --- Build ---------------------------------------------------------------

///|
/// Assemble the accumulated residues into a full AtomArray
/// (backbone + sidechain) via inverse internal-coordinate reconstruction.
pub fn StructureBuilder::build(self : StructureBuilder) -> AtomArray {
  let n = self.residues.length()
  if n == 0 {
    return AtomArray::new()
  }

  let res_names : Array[String] = Array(capacity=n)
  let bb : Array[BackboneAngles] = Array(capacity=n)
  let sc_arr : Array[SidechainAngles] = Array(capacity=n)

  let mut i = 0
  while i < n {
    let r = self.residues[i]
    res_names.push(r.name)

    let phi = r.phi
    let psi = r.psi
    let omega = if i < n - 1 { r.omega } else { None }
    bb.push(BackboneAngles::{ phi, psi, omega, })

    sc_arr.push(SidechainAngles::{
      chi1: r.chi1,
      chi2: r.chi2,
      chi3: r.chi3,
      chi4: r.chi4,
    })

    i = i + 1
  }

  let chain = ChainAngleChain::from_arrays(bb, sc_arr)
  chain.build_full(
    res_names,
    chain_id=self.chain_id,
    start_res_id=self.start_res_id,
  )
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Interpret a sentinel (-999.0) as "use default", otherwise Some(value).
fn opt_or_default(val : Double, def : Double) -> Double? {
  if val < -900.0 {
    Some(def)
  } else {
    Some(val)
  }
}

///|
/// Interpret a sentinel (-999.0) as None, otherwise Some(value).
fn opt_or_none(val : Double) -> Double? {
  if val < -900.0 {
    None
  } else {
    Some(val)
  }
}

///|
/// Decide default χ¹ for a residue name.
/// GLY/ALA have no χ¹; others get 180° (trans) as standard start.
fn StructureBuilder::default_chi_or_none(
  self : StructureBuilder,
  name : String,
) -> Double? {
  match name {
    "GLY" | "ALA" => None
    _ => Some(self.default_chi)
  }
}
