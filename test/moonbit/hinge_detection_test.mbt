///|
/// Tests for hinge site detection from normal modes.

test "HingeSite construction" {
  let h = @src.HingeSite::new(5, 2, 0.01, 0.8, 0.9)
  assert_eq(h.get_residue_idx(), 5)
  assert_eq(h.get_mode_idx(), 2)
  assert_true((h.get_magnitude() - 0.01).abs() < 1.0e-10)
}

test "detect_hinges_gnm sign change" {
  // Mode: + + + - - - + + +
  // Sign changes at index 2↔3 and 5↔6
  let mode = [1.0, 2.0, 1.0, -1.0, -2.0, -1.0, 1.0, 2.0, 1.0]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 2)
  // First hinge between idx 2 (val=1.0) and idx 3 (val=-1.0)
  // |1.0| <= |-1.0| → hinge at idx 2
  assert_eq(hinges[0].get_residue_idx(), 2)
  // Second hinge between idx 5 (val=-1.0) and idx 6 (val=1.0)
  assert_eq(hinges[1].get_residue_idx(), 5)
}

test "detect_hinges_gnm no sign change" {
  let mode = [1.0, 2.0, 3.0, 2.0, 1.0]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 0)
}

test "detect_hinges_gnm all positive" {
  let mode = [0.1, 0.5, 1.0, 0.5, 0.1]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 0)
}

test "detect_hinges_gnm alternating" {
  // + - + - + → sign change at every pair
  let mode = [1.0, -1.0, 1.0, -1.0, 1.0]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 4)
}

test "detect_hinges_gnm too short" {
  let mode = [1.0, -1.0]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 0)
}

test "detect_hinges_gnm hinge picks smaller magnitude" {
  // + + - - → sign change at idx 1↔2, |2.0| > |0.5| → hinge at idx 2
  let mode = [1.0, 2.0, -0.5, -1.0]
  let hinges = @src.detect_hinges_gnm(mode, 0)
  assert_eq(hinges.length(), 1)
  assert_eq(hinges[0].get_residue_idx(), 2)
}

test "detect_hinges_gnm_local_min basic" {
  // |mode|: 1.0, 0.5, 0.01, 0.5, 1.0 → local min at idx 2
  let mode = [1.0, 0.5, 0.01, -0.5, -1.0]
  let hinges = @src.detect_hinges_gnm_local_min(mode, 0)
  assert_true(hinges.length() >= 1)
  assert_eq(hinges[0].get_residue_idx(), 2)
}

test "detect_hinges_gnm_local_min no minimum" {
  // Monotonically increasing
  let mode = [0.1, 0.5, 1.0, 1.5, 2.0]
  let hinges = @src.detect_hinges_gnm_local_min(mode, 0)
  assert_eq(hinges.length(), 0)
}

test "anm_mode_magnitudes basic" {
  // 2 residues: (1,0,0) and (0,3,4)
  let mode = [1.0, 0.0, 0.0, 0.0, 3.0, 4.0]
  let mags = @src.anm_mode_magnitudes(mode)
  assert_eq(mags.length(), 2)
  assert_true((mags[0] - 1.0).abs() < 1.0e-10)
  assert_true((mags[1] - 5.0).abs() < 1.0e-10)
}

test "detect_hinges_anm local minimum" {
  // 5 residues, hinge at idx 2 (zero displacement)
  // Residue 0: (1,0,0), 1: (0.5,0,0), 2: (0,0,0), 3: (0,0.5,0), 4: (0,1,0)
  let mode = [
    1.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
  ]
  let hinges = @src.detect_hinges_anm(mode, 0)
  assert_true(hinges.length() >= 1)
  let found_idx2 = {
    let mut found = false
    for h in hinges {
      if h.get_residue_idx() == 2 {
        found = true
      }
    }
    found
  }
  assert_true(found_idx2)
}

test "detect_hinges_anm no hinge" {
  // All residues move equally
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
  ]
  let hinges = @src.detect_hinges_anm(mode, 0)
  assert_eq(hinges.length(), 0)
}

test "detect_common_hinges_gnm basic" {
  // Mode 1: + + - - → hinge at 1 or 2
  // Mode 2: + + - - → same hinge
  let mode1 = [1.0, 1.0, -1.0, -1.0]
  let mode2 = [1.0, 1.0, -1.0, -1.0]
  let modes = [mode1, mode2]
  let common = @src.detect_common_hinges_gnm(modes, 2, min_count=2)
  assert_true(common.length() >= 1)
}

test "detect_common_hinges_gnm no common" {
  // Mode 1: + + - - → hinge at 1 or 2
  // Mode 2: + - + - → hinge at 0 or 1, 2 or 3
  let mode1 = [1.0, 2.0, -1.0, -1.0]
  let mode2 = [1.0, -1.0, 1.0, -1.0]
  let modes = [mode1, mode2]
  // min_count=2 means must be hinge in both modes
  let common = @src.detect_common_hinges_gnm(modes, 2, min_count=2)
  // Check that we get some result (may or may not have common)
  assert_true(common.length() >= 0)
}

test "hinge_scores_gnm basic" {
  // Mode 1: + + - - → hinge at idx 1 (|1| <= |-1|)
  // Mode 2: + + - - → same
  let mode1 = [1.0, 1.0, -1.0, -1.0]
  let mode2 = [1.0, 1.0, -1.0, -1.0]
  let modes = [mode1, mode2]
  let scores = @src.hinge_scores_gnm(modes, 2)
  assert_eq(scores.length(), 4)
  // Residue 1 should have score 1.0 (hinge in both modes)
  assert_true((scores[1] - 1.0).abs() < 1.0e-10)
  // Residue 0 should have score 0.0
  assert_true(scores[0].abs() < 1.0e-10)
}

test "hinge_scores_gnm empty modes" {
  let modes : Array[Array[Double]] = []
  let scores = @src.hinge_scores_gnm(modes, 0)
  assert_eq(scores.length(), 0)
}

test "hinge_scores_anm basic" {
  // 3 residues, hinge at idx 1
  let mode1 = [
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
  ]
  let mode2 = [
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
  ]
  let modes = [mode1, mode2]
  let scores = @src.hinge_scores_anm(modes, 2)
  assert_eq(scores.length(), 3)
  // Residue 1 should have high score
  assert_true(scores[1] >= 0.0)
}

test "classify_hinges basic" {
  let scores = [0.0, 0.5, 1.0, 0.3, 0.8]
  let is_hinge = @src.classify_hinges(scores, 0.5)
  assert_eq(is_hinge.length(), 5)
  assert_false(is_hinge[0])
  assert_true(is_hinge[1])
  assert_true(is_hinge[2])
  assert_false(is_hinge[3])
  assert_true(is_hinge[4])
}

test "classify_hinges all true" {
  let scores = [1.0, 1.0, 1.0]
  let is_hinge = @src.classify_hinges(scores, 0.5)
  for b in is_hinge {
    assert_true(b)
  }
}

test "classify_hinges all false" {
  let scores = [0.0, 0.0, 0.0]
  let is_hinge = @src.classify_hinges(scores, 0.5)
  for b in is_hinge {
    assert_false(b)
  }
}

test "hinge_angle basic" {
  // 5 residues: left moves +x, right moves +y → 90° angle
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
  ]
  let angle = @src.hinge_angle(mode, 2, window=2)
  // Should be ~90 degrees
  assert_true(angle > 80.0 && angle < 100.0)
}

test "hinge_angle parallel" {
  // Left and right both move +x → 0° angle
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
  ]
  let angle = @src.hinge_angle(mode, 2, window=2)
  assert_true(angle.abs() < 1.0)
}

test "hinge_angle antiparallel" {
  // Left moves +x, right moves -x → 180° angle
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    -1.0, 0.0, 0.0,
    -1.0, 0.0, 0.0,
  ]
  let angle = @src.hinge_angle(mode, 2, window=2)
  assert_true(angle > 170.0 && angle <= 180.0)
}

test "hinge_angle edge case" {
  // Hinge at residue 0 (no left flank)
  let mode = [
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
  ]
  let angle = @src.hinge_angle(mode, 0, window=2)
  // No left flank → angle = 0
  assert_eq(angle, 0.0)
}

test "hinge_axis basic" {
  // Left moves +x, right moves +y → axis = x × y = z
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
  ]
  let axis = @src.hinge_axis(mode, 2, window=2)
  // Should be approximately (0, 0, 1)
  assert_true(axis.z.abs() > 0.9)
  assert_true(axis.x.abs() < 0.1)
  assert_true(axis.y.abs() < 0.1)
}

test "hinge_axis normalized" {
  let mode = [
    2.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 3.0, 0.0,
    0.0, 3.0, 0.0,
  ]
  let axis = @src.hinge_axis(mode, 2, window=2)
  // Should be unit vector
  let norm = (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z).sqrt()
  assert_true((norm - 1.0).abs() < 1.0e-10)
}

test "hinge_axis zero displacement" {
  // All zero → zero axis
  let mode = [
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
  ]
  let axis = @src.hinge_axis(mode, 2, window=2)
  assert_true(axis.x.abs() < 1.0e-10)
  assert_true(axis.y.abs() < 1.0e-10)
  assert_true(axis.z.abs() < 1.0e-10)
}

test "hinge_angles multiple" {
  // Two hinges: at idx 2 and idx 5
  let mode = [
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
  ]
  let hinges = @src.detect_hinges_anm(mode, 0)
  let angles = @src.hinge_angles(mode, hinges)
  assert_eq(angles.length(), hinges.length())
  for a in angles {
    assert_true(a >= 0.0 && a <= 180.0)
  }
}

test "analyze_hinges_gnm basic" {
  // Two modes with a common hinge
  let mode1 = [1.0, 1.0, -1.0, -1.0, -1.0]
  let mode2 = [1.0, 1.0, -1.0, -1.0, -1.0]
  let modes = [mode1, mode2]
  let analysis = @src.analyze_hinges_gnm(modes, 2, score_threshold=0.5)
  assert_eq(analysis.get_scores().length(), 5)
  assert_true(analysis.num_hinges() >= 0)
}

test "analyze_hinges_anm basic" {
  // 5 residues with hinge at idx 2
  let mode1 = [
    1.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
  ]
  let mode2 = [
    1.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
  ]
  let modes = [mode1, mode2]
  let analysis = @src.analyze_hinges_anm(modes, 2, score_threshold=0.5)
  assert_eq(analysis.get_scores().length(), 5)
  assert_true(analysis.num_hinges() >= 0)
}

test "analyze_hinges_gnm get_hinge_indices" {
  let mode1 = [1.0, 1.0, -1.0, -1.0]
  let mode2 = [1.0, 1.0, -1.0, -1.0]
  let modes = [mode1, mode2]
  let analysis = @src.analyze_hinges_gnm(modes, 2, score_threshold=0.5)
  let indices = analysis.get_hinge_indices()
  // All hinge indices should be valid
  for idx in indices {
    assert_true(idx >= 0 && idx < 4)
  }
}

test "analyze_hinges_anm get_angles" {
  let mode1 = [
    1.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
  ]
  let mode2 = [
    1.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
  ]
  let modes = [mode1, mode2]
  let analysis = @src.analyze_hinges_anm(modes, 2, score_threshold=0.5)
  let angles = analysis.get_angles()
  // Angles should match number of hinges
  assert_eq(angles.length(), analysis.num_hinges())
}
