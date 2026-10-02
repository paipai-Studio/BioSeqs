///|
/// Tier 2 ProDy tests: normal mode analysis, perturbation response,
/// domain decomposition, rotamer classification, stereochemistry.

///|
fn make_t2_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (1.0, 0.0, 0.0),
    (0.31, 0.95, 1.5),
    (-0.81, 0.59, 3.0),
    (-0.81, -0.59, 4.5),
    (0.31, -0.95, 6.0),
    (1.0, 0.0, 7.5),
    (0.31, 0.95, 9.0),
    (-0.81, 0.59, 10.5),
  ]
  for i = 0; i < 8; i = i + 1 {
    let (x, y, z) = coords[i]
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, z), resname="ALA", chainid='A', resseq=i + 1, element="C")
    let n_atom = @src.Atom::new(name="N", coord=@src.Vector3::new(x - 0.5, y, z - 0.3), resname="ALA", chainid='A', resseq=i + 1, element="N")
    let c_atom = @src.Atom::new(name="C", coord=@src.Vector3::new(x + 0.5, y, z + 0.3), resname="ALA", chainid='A', resseq=i + 1, element="C")
    let o_atom = @src.Atom::new(name="O", coord=@src.Vector3::new(x + 0.7, y + 0.2, z + 0.5), resname="ALA", chainid='A', resseq=i + 1, element="O")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=i + 1, atoms=[n_atom, ca, c_atom, o_atom])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

// === Normal Mode Analysis ===

test "mode_analysis_from_gnm" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      assert_true(ma.num_modes() > 0)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_from_anm" {
  let chain = make_t2_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let ma = @src.ModeAnalysis::from_anm(anm)
      assert_true(ma.num_modes() > 0)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_eigenvalue" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let ev = ma.eigenvalue(0)
      assert_true(ev >= 0.0)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_eigenvector_length" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let vec = ma.eigenvector(0)
      assert_eq(vec.length(), 8)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_collectivity" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let c = ma.collectivity(0)
      assert_true(c >= 0.0)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_participation_ratio" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let pr = ma.participation_ratio(0)
      assert_true(pr >= 0.0 && pr <= 1.01)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_project" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let disp = Array::make(8, 1.0)
      let proj = ma.project(disp, 0)
      assert_true(proj == proj)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_deformation" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let defm = ma.deformation(0)
      assert_eq(defm.length(), 8)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_bfactors_correlation" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let bf = Array::make(8, 1.0)
      let corr = ma.bfactors_correlation(bf)
      assert_true(corr >= -1.01 && corr <= 1.01)
    }
    None => assert_true(false)
  }
}

test "mode_analysis_mode_energy" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let ma = @src.ModeAnalysis::from_gnm(gnm)
      let e = ma.mode_energy(Array::make(8, 1.0), 0)
      assert_true(e >= 0.0)
    }
    None => assert_true(false)
  }
}

// === Perturbation Response ===

test "perturbation_from_gnm" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      assert_eq(pr.num_residues(), 8)
    }
    None => assert_true(false)
  }
}

test "perturbation_from_anm" {
  let chain = make_t2_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let pr = @src.PerturbationResponse::from_anm(anm)
      assert_eq(pr.num_residues(), 8)
    }
    None => assert_true(false)
  }
}

test "perturbation_response_matrix_shape" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let mat = pr.response_matrix()
      assert_eq(mat.length(), 8)
      assert_eq(mat[0].length(), 8)
    }
    None => assert_true(false)
  }
}

test "perturbation_sensitivity" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let s = pr.sensitivity(0)
      assert_true(s >= 0.0)
    }
    None => assert_true(false)
  }
}

test "perturbation_effectiveness" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let e = pr.effectiveness(0)
      assert_true(e >= 0.0)
    }
    None => assert_true(false)
  }
}

test "perturbation_profile_length" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let prof = pr.perturbation_profile()
      assert_eq(prof.length(), 8)
    }
    None => assert_true(false)
  }
}

test "perturbation_hotspot_residues" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let hotspots = pr.hotspot_residues(3)
      assert_true(hotspots.length() <= 3)
    }
    None => assert_true(false)
  }
}

test "perturbation_dfi" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let d = pr.dfi(0)
      assert_true(d >= 0.0)
    }
    None => assert_true(false)
  }
}

test "perturbation_dfi_profile" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let pr = @src.PerturbationResponse::from_gnm(gnm)
      let prof = pr.dfi_profile()
      assert_eq(prof.length(), 8)
    }
    None => assert_true(false)
  }
}

// === Domain Decomposition ===

test "domain_decomp_from_gnm" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      assert_true(dd.num_domains() >= 1)
    }
    None => assert_true(false)
  }
}

test "domain_decomp_from_anm" {
  let chain = make_t2_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let dd = @src.DomainDecomposition::from_anm(anm, 2)
      assert_true(dd.num_domains() >= 1)
    }
    None => assert_true(false)
  }
}

test "domain_decomp_domains" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      let domains = dd.domains()
      assert_eq(domains.length(), dd.num_domains())
    }
    None => assert_true(false)
  }
}

test "domain_decomp_domain_sizes" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      let sizes = dd.domain_sizes()
      let total = {
        let mut s = 0
        for sz in sizes {
          s = s + sz
        }
        s
      }
      assert_eq(total, 8)
    }
    None => assert_true(false)
  }
}

test "domain_decomp_rigidity_profile" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      let rp = dd.rigidity_profile()
      assert_eq(rp.length(), 8)
    }
    None => assert_true(false)
  }
}

test "domain_decomp_largest_domain" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      let largest = dd.largest_domain_size()
      assert_true(largest >= 1)
    }
    None => assert_true(false)
  }
}

test "domain_decomp_domain_of" {
  let chain = make_t2_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let dd = @src.DomainDecomposition::from_gnm(gnm, 2)
      let d = dd.domain_of(0)
      assert_true(d >= 0)
    }
    None => assert_true(false)
  }
}

// === Rotamer Analysis ===

test "classify_rotamer_trans" {
  let bin = @src.classify_rotamer_angle(180.0)
  assert_true(bin is @src.Trans)
}

test "classify_rotamer_gauche_plus" {
  let bin = @src.classify_rotamer_angle(60.0)
  assert_true(bin is @src.GauchePlus)
}

test "classify_rotamer_gauche_minus" {
  let bin = @src.classify_rotamer_angle(-60.0)
  assert_true(bin is @src.GaucheMinus)
}

test "classify_rotamer_near_trans" {
  let bin = @src.classify_rotamer_angle(170.0)
  assert_true(bin is @src.Trans)
}

test "rotamer_deviation_zero_at_center" {
  let d = @src.rotamer_deviation(180.0)
  assert_true(d < 0.01)
}

test "rotamer_deviation_nonzero" {
  let d = @src.rotamer_deviation(90.0)
  assert_true(d > 20.0)
}

test "rotamer_bin_name" {
  assert_eq(@src.rotamer_bin_name(@src.Trans), "trans")
  assert_eq(@src.rotamer_bin_name(@src.GauchePlus), "gauche+")
  assert_eq(@src.rotamer_bin_name(@src.GaucheMinus), "gauche-")
}

test "has_chi1_serine" {
  assert_true(@src.has_chi1("SER"))
  assert_true(not(@src.has_chi1("GLY")))
  assert_true(not(@src.has_chi1("ALA")))
}

test "has_chi2_valine" {
  assert_true(not(@src.has_chi2("VAL")))
  assert_true(@src.has_chi2("PHE"))
}

test "has_chi3_lysine" {
  assert_true(@src.has_chi3("LYS"))
  assert_true(not(@src.has_chi3("VAL")))
}

test "has_chi4_arginine" {
  assert_true(@src.has_chi4("ARG"))
  assert_true(not(@src.has_chi4("MET")))
}

test "compute_chi_angles_basic" {
  let n = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let ca = Some(@src.Vector3::new(1.5, 0.0, 0.0))
  let cb = Some(@src.Vector3::new(2.0, 1.0, 0.0))
  let cg = Some(@src.Vector3::new(2.0, 2.0, 1.0))
  let (chi1, _, _, _) = @src.compute_chi_angles(n, ca, cb, cg, None, None, None)
  assert_true(chi1 == chi1)
}

test "compute_chi_angles_missing" {
  let nan = 0.0 / 0.0
  let (chi1, chi2, chi3, chi4) = @src.compute_chi_angles(None, None, None, None, None, None, None)
  assert_true(chi1 != chi1)
  assert_true(chi2 != chi2)
  assert_true(chi3 != chi3)
  assert_true(chi4 != chi4)
  ignore(nan)
}

test "classify_rotamer_outlier" {
  let nan = 0.0 / 0.0
  let r = @src.classify_rotamer(1, "SER", 0.0, nan, nan, nan)
  assert_true(r.is_outlier)
}

test "classify_rotamer_normal" {
  let nan = 0.0 / 0.0
  let r = @src.classify_rotamer(1, "SER", 60.0, nan, nan, nan)
  assert_true(not(r.is_outlier))
}

test "rotamer_outliers_filter" {
  let nan = 0.0 / 0.0
  let r1 = @src.classify_rotamer(1, "SER", 60.0, nan, nan, nan)
  let r2 = @src.classify_rotamer(2, "SER", 0.0, nan, nan, nan)
  let outliers = @src.rotamer_outliers([r1, r2])
  assert_eq(outliers.length(), 1)
}

test "rotamer_summary_counts" {
  let nan = 0.0 / 0.0
  let r1 = @src.classify_rotamer(1, "SER", 60.0, nan, nan, nan)
  let r2 = @src.classify_rotamer(2, "SER", 180.0, nan, nan, nan)
  let r3 = @src.classify_rotamer(3, "SER", -60.0, nan, nan, nan)
  let (total, outliers, trans, gp, gm) = @src.rotamer_summary([r1, r2, r3])
  assert_eq(total, 3)
  assert_eq(outliers, 0)
  assert_eq(trans, 1)
  assert_eq(gp, 1)
  assert_eq(gm, 1)
}

test "rotamer_probability_proline" {
  let p = @src.rotamer_probability("PRO", @src.GaucheMinus, @src.Trans)
  assert_true(p > 0.0)
}

// === Stereochemistry ===

test "detect_chirality_l_amino" {
  let n = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let ca = Some(@src.Vector3::new(1.0, 0.0, 0.0))
  let c = Some(@src.Vector3::new(1.0, 1.0, 0.0))
  let cb = Some(@src.Vector3::new(1.0, 0.0, 1.0))
  let r = @src.detect_chirality(1, "ALA", n, ca, c, cb)
  assert_true(r.has_chirality)
}

test "detect_chirality_missing" {
  let r = @src.detect_chirality(1, "GLY", None, None, None, None)
  assert_true(not(r.has_chirality))
}

test "detect_cis_peptide_trans" {
  let cp = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let n = Some(@src.Vector3::new(1.0, 0.0, 0.0))
  let ca = Some(@src.Vector3::new(1.5, 0.5, 0.0))
  let c = Some(@src.Vector3::new(2.0, 0.0, 0.0))
  let r = @src.detect_cis_peptide(1, 2, cp, n, ca, c)
  assert_true(not(r.is_cis))
}

test "detect_cis_peptide_missing" {
  let r = @src.detect_cis_peptide(1, 2, None, None, None, None)
  assert_true(not(r.is_cis))
}

test "count_cis_peptides" {
  let nan = 0.0 / 0.0
  ignore(nan)
  let r1 = @src.CisPeptideResult::{ resseq1: 1, resseq2: 2, omega: 180.0, is_cis: false }
  let r2 = @src.CisPeptideResult::{ resseq1: 2, resseq2: 3, omega: 0.0, is_cis: true }
  let r3 = @src.CisPeptideResult::{ resseq1: 3, resseq2: 4, omega: 10.0, is_cis: true }
  let count = @src.count_cis_peptides([r1, r2, r3])
  assert_eq(count, 2)
}

test "validate_bond_length_normal" {
  let p1 = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let p2 = Some(@src.Vector3::new(1.5, 0.0, 0.0))
  let r = @src.validate_bond_length("CA", "CB", 1, p1, p2, 1.54, 0.1)
  assert_true(not(r.is_outlier))
}

test "validate_bond_length_outlier" {
  let p1 = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let p2 = Some(@src.Vector3::new(3.0, 0.0, 0.0))
  let r = @src.validate_bond_length("CA", "CB", 1, p1, p2, 1.54, 0.1)
  assert_true(r.is_outlier)
}

test "validate_bond_length_missing" {
  let r = @src.validate_bond_length("CA", "CB", 1, None, None, 1.54, 0.1)
  assert_true(not(r.is_outlier))
}

test "validate_bond_angle_normal" {
  let p1 = Some(@src.Vector3::new(1.0, 0.0, 0.0))
  let p2 = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let p3 = Some(@src.Vector3::new(0.0, 1.0, 0.0))
  let r = @src.validate_bond_angle("N", "CA", "C", 1, p1, p2, p3, 90.0, 15.0)
  assert_true(not(r.is_outlier))
}

test "validate_bond_angle_outlier" {
  let p1 = Some(@src.Vector3::new(1.0, 0.0, 0.0))
  let p2 = Some(@src.Vector3::new(0.0, 0.0, 0.0))
  let p3 = Some(@src.Vector3::new(-0.9, 0.1, 0.0))
  let r = @src.validate_bond_angle("N", "CA", "C", 1, p1, p2, p3, 109.5, 5.0)
  assert_true(r.is_outlier)
}

test "peptide_planarity_trans" {
  let r = @src.peptide_planarity(1, 180.0, 10.0)
  assert_true(r.is_planar)
}

test "peptide_planarity_nonplanar" {
  let r = @src.peptide_planarity(1, 90.0, 10.0)
  assert_true(not(r.is_planar))
}

test "standard_bond_length_peptide" {
  let l = @src.standard_bond_length("C", "N")
  assert_true((l - 1.33).abs() < 0.001)
}

test "standard_bond_length_ca_cb" {
  let l = @src.standard_bond_length("CA", "CB")
  assert_true((l - 1.54).abs() < 0.001)
}

test "standard_bond_angle_ca" {
  let a = @src.standard_bond_angle("N", "CA", "C")
  assert_true((a - 110.0).abs() < 0.001)
}

test "validate_bond_lengths_filter" {
  let r1 = @src.BondValidationResult::{ atom1: "CA", atom2: "CB", resseq: 1, observed: 1.54, expected: 1.54, deviation: 0.0, is_outlier: false }
  let r2 = @src.BondValidationResult::{ atom1: "CA", atom2: "CB", resseq: 2, observed: 3.0, expected: 1.54, deviation: 1.46, is_outlier: true }
  let outliers = @src.validate_bond_lengths([r1, r2])
  assert_eq(outliers.length(), 1)
}

test "validate_bond_angles_filter" {
  let r1 = @src.AngleValidationResult::{ atom1: "N", atom2: "CA", atom3: "C", resseq: 1, observed: 110.0, expected: 110.0, deviation: 0.0, is_outlier: false }
  let r2 = @src.AngleValidationResult::{ atom1: "N", atom2: "CA", atom3: "C", resseq: 2, observed: 90.0, expected: 110.0, deviation: -20.0, is_outlier: true }
  let outliers = @src.validate_bond_angles([r1, r2])
  assert_eq(outliers.length(), 1)
}

test "stereo_summary_basic" {
  let chir = [@src.ChiralityResult::{ resseq: 1, resname: "ALA", has_chirality: true, is_correct: true, volume: 1.0 }]
  let cis = [@src.CisPeptideResult::{ resseq1: 1, resseq2: 2, omega: 180.0, is_cis: false }]
  let bonds = [@src.BondValidationResult::{ atom1: "CA", atom2: "CB", resseq: 1, observed: 1.54, expected: 1.54, deviation: 0.0, is_outlier: false }]
  let angles = [@src.AngleValidationResult::{ atom1: "N", atom2: "CA", atom3: "C", resseq: 1, observed: 110.0, expected: 110.0, deviation: 0.0, is_outlier: false }]
  let planarity = [@src.PeptidePlanarityResult::{ resseq: 1, omega: 180.0, deviation_from_180: 0.0, is_planar: true }]
  let s = @src.stereo_summary(chir, cis, bonds, angles, planarity)
  assert_eq(s.num_residues, 1)
  assert_eq(s.num_chiral, 1)
  assert_eq(s.num_wrong_chirality, 0)
  assert_eq(s.num_cis_peptides, 0)
  assert_eq(s.num_bond_outliers, 0)
  assert_eq(s.num_angle_outliers, 0)
  assert_eq(s.num_nonplanar_peptides, 0)
}

test "stereo_summary_with_issues" {
  let chir = [
    @src.ChiralityResult::{ resseq: 1, resname: "ALA", has_chirality: true, is_correct: true, volume: 1.0 },
    @src.ChiralityResult::{ resseq: 2, resname: "ALA", has_chirality: true, is_correct: false, volume: -1.0 },
  ]
  let cis = [
    @src.CisPeptideResult::{ resseq1: 1, resseq2: 2, omega: 0.0, is_cis: true },
  ]
  let bonds = [
    @src.BondValidationResult::{ atom1: "CA", atom2: "CB", resseq: 1, observed: 3.0, expected: 1.54, deviation: 1.46, is_outlier: true },
  ]
  let angles = [
    @src.AngleValidationResult::{ atom1: "N", atom2: "CA", atom3: "C", resseq: 1, observed: 90.0, expected: 110.0, deviation: -20.0, is_outlier: true },
  ]
  let planarity = [
    @src.PeptidePlanarityResult::{ resseq: 1, omega: 90.0, deviation_from_180: 90.0, is_planar: false },
  ]
  let s = @src.stereo_summary(chir, cis, bonds, angles, planarity)
  assert_eq(s.num_residues, 2)
  assert_eq(s.num_chiral, 2)
  assert_eq(s.num_wrong_chirality, 1)
  assert_eq(s.num_cis_peptides, 1)
  assert_eq(s.num_bond_outliers, 1)
  assert_eq(s.num_angle_outliers, 1)
  assert_eq(s.num_nonplanar_peptides, 1)
}
