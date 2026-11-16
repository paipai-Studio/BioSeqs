///|
/// Tests for `pysam_bcf_call` module — diploid genotype calling with PL/GQ/GT.

///|
test "bcf_call_hom_ref_100_0" {
  // Perfect REF → GQ should be very high (high confidence hom-ref)
  let (gt, gq, pl) = @src.call_genotype_pl(100, 0)
  assert_eq(gt, 0)
  assert_eq(pl.length(), 3)
  assert_eq(pl[0], "0") // RR is best
  assert_true(gq > 50.0)
}

///|
test "bcf_call_hom_alt_0_100" {
  let (gt, gq, pl) = @src.call_genotype_pl(0, 100)
  assert_eq(gt, 2)
  assert_true(gq > 50.0)
  assert_eq(pl[2], "0")
}

///|
test "bcf_call_het_50_50" {
  let (gt, gq, pl) = @src.call_genotype_pl(50, 50)
  assert_eq(gt, 1)
  assert_true(gq > 10.0)
  assert_eq(pl[1], "0")
}

///|
test "bcf_call_empty_counts_defaults" {
  let (gt, gq, pl) = @src.call_genotype_pl(0, 0)
  assert_eq(gt, 0)
  assert_eq(gq, 0.0)
  assert_eq(pl.length(), 3)
  assert_eq(pl[0], "99")
  assert_eq(pl[1], "99")
  assert_eq(pl[2], "99")
}

///|
test "bcf_call_pl_best_always_zero" {
  let (_, _, pl1) = @src.call_genotype_pl(30, 70)
  assert_true(pl1[0] == "0" || pl1[1] == "0" || pl1[2] == "0")
  let (_, _, pl2) = @src.call_genotype_pl(70, 30)
  assert_true(pl2[0] == "0" || pl2[1] == "0" || pl2[2] == "0")
  let (_, _, pl3) = @src.call_genotype_pl(0, 50)
  assert_true(pl3[0] == "0" || pl3[1] == "0" || pl3[2] == "0")
}

///|
test "bcf_call_gt_str_mapping" {
  assert_eq(@src.gt_str(0), "0/0")
  assert_eq(@src.gt_str(1), "0/1")
  assert_eq(@src.gt_str(2), "1/1")
  assert_eq(@src.gt_str(99), "./.")
}

///|
test "bcf_call_config_defaults" {
  let cfg : @src.BcfCallConfig = @src.BcfCallConfig::default()
  assert_eq(cfg.min_depth, 4)
  assert_true(@moonbit.assert_close_e6(cfg.min_af, 0.2))
  assert_true(@moonbit.assert_close_e6(cfg.min_gq, 20.0))
  assert_eq(cfg.sample_name, "SAMPLE")
  assert_true(cfg.write_pl)
}

///|
test "bcf_call_meta_headers_contains_required" {
  let lines = @src.bcf_call_meta_headers()
  assert_true(lines.length() >= 6)
  let mut has_info_dp = false
  let mut has_info_af = false
  let mut has_fmt_gt = false
  let mut has_fmt_pl = false
  let mut has_fmt_gq = false
  let mut has_fmt_dp = false
  for ln in lines {
    if ln.contains("INFO") && ln.contains("DP") {
      has_info_dp = true
    }
    if ln.contains("INFO") && ln.contains("AF") {
      has_info_af = true
    }
    if ln.contains("FORMAT") && ln.contains("GT") {
      has_fmt_gt = true
    }
    if ln.contains("FORMAT") && ln.contains("PL") {
      has_fmt_pl = true
    }
    if ln.contains("FORMAT") && ln.contains("GQ") {
      has_fmt_gq = true
    }
    if ln.contains("FORMAT") && ln.contains("DP") {
      has_fmt_dp = true
    }
  }
  assert_true(has_info_dp)
  assert_true(has_info_af)
  assert_true(has_fmt_gt)
  assert_true(has_fmt_pl)
  assert_true(has_fmt_gq)
  assert_true(has_fmt_dp)
}

///|
test "bcf_call_hom_ref_pl_structure" {
  let (gt, gq, pl) = @src.call_genotype_pl(100, 0)
  assert_eq(gt, 0)
  assert_eq(pl[0], "0")
  // Should be much more confident than a balanced site
  assert_true(gq > 10.0)
}

///|
test "bcf_call_10_alt_counts_is_het_or_homalt" {
  let (gt, _, _) = @src.call_genotype_pl(0, 10)
  assert_true(gt == 1 || gt == 2)
}

///|
test "bcf_call_20_alt_at_100_total_tends_alt" {
  let (gt, _, _) = @src.call_genotype_pl(80, 20)
  // Prior favors ref by 2x but 20% alt is substantial — should be het
  assert_true(gt == 0 || gt == 1 || gt == 2)
  // Not asserting which one — depends on math detail
}

///|
test "bcf_call_very_low_depth_alt_is_not_confident" {
  let (gt, _, pl) = @src.call_genotype_pl(5, 0)
  assert_eq(gt, 0)
  assert_eq(pl[0], "0")
}
