///|
/// Tests for ProtParam protein analysis module.

///|
/// Standard ProtParam docstring example sequence (152 residues).
let pp_long_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"

///|
fn pp_close(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-9
}

///|
fn pp_vec_close(a : Array[Double], b : Array[Double]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let mut i = 0
  while i < a.length() {
    if (a[i] - b[i]).abs() >= 1.0e-9 {
      return false
    }
    i = i + 1
  }
  true
}

///|
test "protparam_molecular_weight" {
  let protein = @src.ProteinAnalysis::new("ACDEFGHIKLMNPQRSTVWY")
  let mw = protein.molecular_weight()
  assert_true(mw > 2000.0 && mw < 2500.0)
}

///|
test "protparam_count_amino_acids" {
  let protein = @src.ProteinAnalysis::new("AABBCC")
  let counts = protein.count_amino_acids()
  assert_eq(counts.get('A'), Some(2))
  assert_eq(counts.get('B'), Some(2))
  assert_eq(counts.get('C'), Some(2))
}

///|
test "protparam_instability_index" {
  let stable = @src.ProteinAnalysis::new("AVGLEALVLGLLAVGLLPGL")
  let index = stable.instability_index()
  assert_true(index < 40.0)
}

///|
test "protparam_gravy" {
  let hydrophobic = @src.ProteinAnalysis::new("LLLLLL")
  let gravy = hydrophobic.gravy()
  assert_true(gravy > 3.0)

  let hydrophilic = @src.ProteinAnalysis::new("RRRRRR")
  let gravy2 = hydrophilic.gravy()
  assert_true(gravy2 < -3.0)
}

///|
test "protparam_aromaticity" {
  let protein = @src.ProteinAnalysis::new("FFWWYY")
  let aroma = protein.aromaticity()
  assert_eq(aroma, 1.0)
}

///|
test "protparam_isoelectric_point" {
  let acidic = @src.ProteinAnalysis::new("EEEEEE")
  let pi = acidic.isoelectric_point()
  assert_true(pi < 5.0)

  let basic = @src.ProteinAnalysis::new("KKKKKK")
  let pi2 = basic.isoelectric_point()
  assert_true(pi2 > 9.0)
}

///|
test "protparam_charge_at_ph" {
  let neutral = @src.ProteinAnalysis::new("G")
  let charge = neutral.charge_at_ph(7.0)
  assert_true(charge > -0.5 && charge < 0.5)
}

///|
test "protparam_secondary_structure" {
  let protein = @src.ProteinAnalysis::new("AEAEAE")
  let (h, s, c) = protein.secondary_structure()
  assert_true(h > 1.0)
}

///|
test "protparam_amino_acid_composition" {
  let protein = @src.ProteinAnalysis::new("AABB")
  let comp = protein.amino_acid_composition()
  assert_eq(comp.get('A'), Some(0.5))
  assert_eq(comp.get('B'), Some(0.5))
}

///|
test "protparam_signal_peptide" {
  let sp_seq = "MRVLLFLLLIAAPVLAAPVTLSRGA"
  let protein = @src.ProteinAnalysis::new(sp_seq)
  let result = protein.predict_signal_peptide()
  assert_true(result.has_signal_peptide)
}

///|
test "protparam_length" {
  let protein = @src.ProteinAnalysis::new("ABCDE")
  assert_eq(protein.length(), 5)
}

///|
test "protparam_get_amino_acids_percent" {
  let protein = @src.ProteinAnalysis::new("AABBCC")
  let percent = protein.get_amino_acids_percent()
  assert_eq(percent.get('A'), Some(33.33333333333333))
}

///|
test "protparam_empty_sequence" {
  let protein = @src.ProteinAnalysis::new("")
  assert_eq(protein.length(), 0)
  assert_eq(protein.molecular_weight(), 0.0)
  assert_eq(protein.instability_index(), 0.0)
}

///|
test "protparam_small_sequence" {
  let protein = @src.ProteinAnalysis::new("A")
  assert_eq(protein.length(), 1)
  assert_true(protein.molecular_weight() > 70.0)
}

///|
test "protparam_instability_index_small" {
  let protein = @src.ProteinAnalysis::new("AB")
  let index = protein.instability_index()
  assert_true(index >= 0.0)
}

///|
test "protparam_secondary_structure_fraction_long" {
  // Authoritative values from Biopython 1.88 ProteinAnalysis(X)
  let protein = @src.ProteinAnalysis::new(pp_long_seq)
  let (helix, turn, sheet) = protein.secondary_structure_fraction()
  assert_true(pp_close(helix, 0.32894736842105265))
  assert_true(pp_close(turn, 0.2894736842105263))
  assert_true(pp_close(sheet, 0.368421052631579))
}

///|
test "protparam_secondary_structure_fraction_groups" {
  // L belongs to both helix and sheet.
  let (h1, t1, s1) = @src.ProteinAnalysis::new("MMMMM").secondary_structure_fraction()
  assert_true(pp_close(h1, 1.0) && pp_close(t1, 0.0) && pp_close(s1, 0.0))
  let (h2, t2, s2) = @src.ProteinAnalysis::new("NPGSD").secondary_structure_fraction()
  assert_true(pp_close(h2, 0.0) && pp_close(t2, 1.0) && pp_close(s2, 0.0))
  let (h3, t3, s3) = @src.ProteinAnalysis::new("VIYFWLT").secondary_structure_fraction()
  assert_true(pp_close(h3, 1.0 / 7.0) && pp_close(t3, 0.0) && pp_close(s3, 1.0))
  let (h4, t4, s4) = @src.ProteinAnalysis::new("").secondary_structure_fraction()
  assert_true(h4 == 0.0 && t4 == 0.0 && s4 == 0.0)
}

///|
test "protparam_molar_extinction_coefficient" {
  let (red_x, ox_x) = @src.ProteinAnalysis::new(pp_long_seq).molar_extinction_coefficient()
  assert_eq(red_x, 17420)
  assert_eq(ox_x, 17545)

  // Single residues and cystine pairing (odd Cys -> one unpaired).
  assert_eq(
    @src.ProteinAnalysis::new("W").molar_extinction_coefficient(),
    (5500, 5500),
  )
  assert_eq(
    @src.ProteinAnalysis::new("Y").molar_extinction_coefficient(),
    (1490, 1490),
  )
  assert_eq(
    @src.ProteinAnalysis::new("C").molar_extinction_coefficient(),
    (0, 0),
  )
  assert_eq(
    @src.ProteinAnalysis::new("CC").molar_extinction_coefficient(),
    (0, 125),
  )
  assert_eq(
    @src.ProteinAnalysis::new("CCC").molar_extinction_coefficient(),
    (0, 125),
  )
  assert_eq(
    @src.ProteinAnalysis::new("WYYCCC").molar_extinction_coefficient(),
    (8480, 8605),
  )
  assert_eq(
    @src.ProteinAnalysis::new("MMMMM").molar_extinction_coefficient(),
    (0, 0),
  )

  // Lower-case input is upper-cased, mirroring Biopython's constructor.
  assert_eq(
    @src.ProteinAnalysis::new("wwwyyy").molar_extinction_coefficient(),
    (20970, 20970),
  )
}

///|
test "protparam_flexibility_long" {
  // flexibility() returns length - 9 scores (152 - 9 = 143).
  let flex = @src.ProteinAnalysis::new(pp_long_seq).flexibility()
  assert_eq(flex.length(), 143)
  assert_true(pp_close(flex[0], 0.9825119047619049))
  assert_true(pp_close(flex[1], 1.0166904761904763))
  assert_true(pp_close(flex[2], 0.9947857142857144))
  assert_true(pp_close(flex[3], 0.9660238095238095))
  assert_true(pp_close(flex[4], 0.9890714285714285))
  assert_true(pp_close(flex[71], 1.0289761904761903))
  assert_true(pp_close(flex[140], 0.9785833333333332))
  assert_true(pp_close(flex[141], 0.9609642857142857))
  assert_true(pp_close(flex[142], 0.9650833333333334))
  let mut sum = 0.0
  for v in flex {
    sum = sum + v
  }
  assert_true(pp_close(sum, 143.82985714285715))
}

///|
test "protparam_flexibility_small_full_vector" {
  let small = "ACDEFGHIKLMNPQRSTVWYK"
  let expected = [
    1.0175714285714286, 0.9786785714285715, 0.9849880952380953, 1.019107142857143,
    0.9594642857142857, 0.9993928571428573, 1.0194880952380951, 1.0140357142857144,
    1.0185714285714287, 1.010904761904762, 1.0184761904761905, 0.9913809523809524,
  ]
  let flex = @src.ProteinAnalysis::new(small).flexibility()
  assert_eq(flex.length(), 12)
  assert_true(pp_vec_close(flex, expected))
}

///|
test "protparam_flexibility_short_and_unknown" {
  // n <= 9 yields an empty profile (Python range(n - 9) is empty).
  assert_eq(@src.ProteinAnalysis::new("ACDEFGHIK").flexibility().length(), 0)
  assert_eq(@src.ProteinAnalysis::new("ACDE").flexibility().length(), 0)
  // n = 10 yields exactly one window.
  assert_eq(@src.ProteinAnalysis::new("ACDEFGHIKL").flexibility().length(), 1)

  // Non-standard residues at a covered window position raise.
  let raised = try {
    ignore(@src.ProteinAnalysis::new("XAAAAAAAAA").flexibility())
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "protparam_protein_scale_kd9_long" {
  let kd = @src.gravy_scale_table("KyteDoolitle")
  let profile = @src.ProteinAnalysis::new(pp_long_seq).protein_scale(kd, 9)
  assert_eq(profile.length(), 144)
  assert_true(pp_close(profile[0], 0.24444444444444427))
  assert_true(pp_close(profile[1], -0.04444444444444448))
  assert_true(pp_close(profile[2], -0.04444444444444446))
  assert_true(pp_close(profile[3], 0.7666666666666666))
  assert_true(pp_close(profile[70], -1.1444444444444446))
  assert_true(pp_close(profile[71], -1.2222222222222223))
  assert_true(pp_close(profile[72], -0.588888888888889))
  assert_true(pp_close(profile[73], -0.288888888888889))
  assert_true(pp_close(profile[142], 1.4888888888888887))
  assert_true(pp_close(profile[143], 2.388888888888889))
  let mut sum = 0.0
  for v in profile {
    sum = sum + v
  }
  assert_true(pp_close(sum, -99.14444444444445))
}

///|
test "protparam_protein_scale_kd5_edge04" {
  let kd = @src.gravy_scale_table("KyteDoolitle")
  // weights for window 5, edge 0.4: [0.4, 0.7, 1, 0.7, 0.4]
  let profile = @src.ProteinAnalysis::new(pp_long_seq).protein_scale(
    kd,
    5,
    edge=0.4,
  )
  assert_eq(profile.length(), 148)
  assert_true(pp_close(profile[0], -0.9875))
  assert_true(pp_close(profile[1], -0.8687499999999998))
  assert_true(pp_close(profile[2], -0.7218749999999999))
  assert_true(pp_close(profile[3], 0.35000000000000014))
  assert_true(pp_close(profile[146], 1.3124999999999998))
  assert_true(pp_close(profile[147], 1.4874999999999998))
  let mut sum = 0.0
  for v in profile {
    sum = sum + v
  }
  assert_true(pp_close(sum, -96.121875))
}

///|
test "protparam_protein_scale_hw19" {
  let hw = @src.hydrophilicity_table()
  let profile = @src.ProteinAnalysis::new(pp_long_seq).protein_scale(hw, 19)
  assert_eq(profile.length(), 134)
  assert_true(pp_close(profile[0], -0.11052631578947364))
  assert_true(pp_close(profile[1], -0.0421052631578947))
  assert_true(pp_close(profile[2], -0.01578947368421053))
  assert_true(pp_close(profile[66], 0.005263157894736888))
  assert_true(pp_close(profile[133], -0.30000000000000004))
  let mut sum = 0.0
  for v in profile {
    sum = sum + v
  }
  assert_true(pp_close(sum, 16.18947368421053))
}

///|
test "protparam_protein_scale_small_full_vectors" {
  let small = "ACDEFGHIKLMNPQRSTVWYK"
  let kd = @src.gravy_scale_table("KyteDoolitle")
  let p5 = @src.ProteinAnalysis::new(small).protein_scale(kd, 5)
  let expected5 = [
    0.019999999999999928, -0.42000000000000004, -1.56, 0.039999999999999925, -0.040000000000000126,
    0.16000000000000006, 0.62, 0.5599999999999999, -0.66, -0.5800000000000001, -2.2399999999999998,
    -2.78, -2.2199999999999998, -1.06, -0.54, 0.1, -0.5199999999999998,
  ]
  assert_true(pp_vec_close(p5, expected5))

  let p7 = @src.ProteinAnalysis::new(small).protein_scale(kd, 7, edge=0.5)
  let expected7 = [
    -0.6766666666666667, -0.2833333333333334, -0.7533333333333333, -0.07333333333333343,
    0.64, 0.039999999999999945, 0.07999999999999999, -0.2033333333333333, -1.45,
    -1.5333333333333332, -2.17, -1.853333333333333, -1.1833333333333331, -0.7733333333333332,
    -0.5466666666666665,
  ]
  assert_true(pp_vec_close(p7, expected7))
}

///|
test "protparam_protein_scale_unknown_and_errors" {
  let kd = @src.gravy_scale_table("KyteDoolitle")
  // Unknown residues: their weighted pairs and centre are skipped, like
  // Biopython's stderr-warning behaviour: AXBXCX window 5 -> [0.86, 0.0].
  let weird = @src.ProteinAnalysis::new("AXBXCX").protein_scale(kd, 5)
  assert_eq(weird.length(), 2)
  assert_true(pp_close(weird[0], 0.86))
  assert_true(pp_close(weird[1], 0.0))

  // Window larger than the sequence yields no windows.
  assert_eq(@src.ProteinAnalysis::new("ACDE").protein_scale(kd, 9).length(), 0)

  // Window < 1 raises SeqError.
  let raised = try {
    ignore(@src.ProteinAnalysis::new("ACDEFG").protein_scale(kd, 0))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "protparam_gravy_named_scales" {
  let protein = @src.ProteinAnalysis::new(pp_long_seq)
  // Authoritative GRAVY values for every registered scale (Biopython 1.88).
  let expected : Array[(String, Double)] = [
    ("KyteDoolitle", -0.5973684210526317),
    ("Aboderin", 4.567105263157895),
    ("AbrahamLeo", 0.23776315789473684),
    ("Argos", 0.8606578947368421),
    ("BlackMould", 0.507375),
    ("BullBreese", -0.044473684210526304),
    ("Casari", -0.24144736842105266),
    ("Cid", -0.06776315789473684),
    ("Cowan3.4", 0.02335526315789476),
    ("Cowan7.5", -0.07328947368421052),
    ("Eisenberg", -0.04348684210526315),
    ("Engelman", 1.5999999999999999),
    ("Fasman", -0.36144736842105263),
    ("Fauchere", 0.32703947368421055),
    ("GoldSack", 1.1563815789473686),
    ("Guy", 0.0675),
    ("Jones", 1.2230263157894736),
    ("Juretic", -0.6672368421052631),
    ("Kidera", 0.13828947368421052),
    ("Miyazawa", 5.310921052631579),
    ("Parker", 1.7486842105263158),
    ("Ponnuswamy", 0.3491447368421053),
    ("Rose", 0.7146710526315789),
    ("Roseman", -0.47289473684210526),
    ("Sweet", -0.07907894736842107),
    ("Tanford", 0.0625),
    ("Wilson", 1.549342105263158),
    ("Zimmerman", 1.284078947368421),
  ]
  assert_eq(expected.length(), 28)
  for pair in expected {
    assert_true(pp_close(protein.gravy_scale(pair.0), pair.1))
  }
  // The named Kyte-Doolittle scale agrees with the classic gravy() helper.
  assert_true(pp_close(protein.gravy_scale("KyteDoolitle"), protein.gravy()))
  assert_eq(@src.gravy_scale_names().length(), 28)

  // Unknown scale names raise SeqError (Biopython raises ValueError).
  let raised = try {
    ignore(protein.gravy_scale("NoSuchScale"))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
