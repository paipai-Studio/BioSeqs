///|
/// Tests for motifbreakR — motif disruption variant prediction.

///|
test "motifbreakr_variant construction" {
  let v = @src.MotifbreakrVariant::new(
    name="rs1",
    chrom="chr1",
    position=5,
    ref_allele="A",
    alt_allele="G",
  )
  assert_eq(v.name, "rs1")
  assert_eq(v.position, 5)
  assert_eq(v.ref_allele, "A")
  assert_eq(v.alt_allele, "G")
}

///|
test "motifbreakr_effect enum equality" {
  assert_eq(@src.MotifbreakrEffect::Gain, @src.MotifbreakrEffect::Gain)
  assert_true(@src.MotifbreakrEffect::Gain != @src.MotifbreakrEffect::Loss)
  assert_true(@src.MotifbreakrEffect::Gain != @src.MotifbreakrEffect::Neutral)
}

///|
test "motifbreakr_motif_break basic" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let seq = @src.motifbreakr_sample_sequence()
  let variants = @src.motifbreakr_sample_variants()
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", seq)
    m
  }
  let results = @src.motif_break(
    variants,
    [("GATA1", pwm)],
    sequences,
  )
  assert_true(results.length() > 0)
  for r in results {
    assert_eq(r.chrom, "chr1")
    assert_true(r.ref_score > -1.0e8)
    assert_true(r.alt_score > -1.0e8)
  }
}

///|
test "motifbreakr_motif_break gain variant" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let seq = "GGGGAGATAAGGGG"
  let variants = [
    @src.MotifbreakrVariant::new(
      name="rs_gain",
      chrom="chr1",
      position=4,
      ref_allele="G",
      alt_allele="A",
    ),
  ]
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", seq)
    m
  }
  let results = @src.motif_break(variants, [("GATA1", pwm)], sequences)
  assert_true(results.length() >= 1)
}

///|
test "motifbreakr_significant filters neutral" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let seq = @src.motifbreakr_sample_sequence()
  let variants = @src.motifbreakr_sample_variants()
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", seq)
    m
  }
  let results = @src.motif_break(variants, [("GATA1", pwm)], sequences)
  let sig = @src.motifbreakr_significant(results)
  for r in sig {
    assert_true(r.effect != @src.MotifbreakrEffect::Neutral)
  }
}

///|
test "motifbreakr_sort_by_delta descending" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let seq = @src.motifbreakr_sample_sequence()
  let variants = @src.motifbreakr_sample_variants()
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", seq)
    m
  }
  let results = @src.motif_break(variants, [("GATA1", pwm)], sequences)
  let sorted = @src.motifbreakr_sort_by_delta(results)
  for i in 0..<(sorted.length() - 1) {
    let abs_i = if sorted[i].delta < 0.0 {
      0.0 - sorted[i].delta
    } else {
      sorted[i].delta
    }
    let abs_next = if sorted[i + 1].delta < 0.0 {
      0.0 - sorted[i + 1].delta
    } else {
      sorted[i + 1].delta
    }
    assert_true(abs_i >= abs_next)
  }
}

///|
test "motifbreakr_motif_break missing chrom" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let variants = [
    @src.MotifbreakrVariant::new(
      name="rs1",
      chrom="chr2",
      position=5,
      ref_allele="A",
      alt_allele="G",
    ),
  ]
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", "ACGTACGT")
    m
  }
  let results = @src.motif_break(variants, [("M", pwm)], sequences)
  assert_eq(results.length(), 0)
}

///|
test "motifbreakr_motif_break variant outside motif width" {
  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let variants = [
    @src.MotifbreakrVariant::new(
      name="rs_far",
      chrom="chr1",
      position=0,
      ref_allele="G",
      alt_allele="T",
    ),
  ]
  let sequences = {
    let m : Map[String, String] = Map([])
    m.set("chr1", "GGGGAGATAAGGGG")
    m
  }
  let results = @src.motif_break(variants, [("GATA1", pwm)], sequences)
  for r in results {
    assert_eq(r.variant_name, "rs_far")
  }
}