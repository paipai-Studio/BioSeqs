///|
fn main {
  println("=== motifbreakR Demo ===")

  let pwm = @src.motifbreakr_sample_gata1_pwm()
  let seq = @src.motifbreakr_sample_sequence()
  let variants = @src.motifbreakr_sample_variants()
  let sequences = {
    let m : Map[String, String] = Map::new()
    m.set("chr1", seq)
    m
  }

  println("\n1. Sequence: " + seq)
  println("   Variants: " + variants.length().to_string())
  for v in variants {
    println(
      "   " + v.name + ": " + v.chrom + ":" + v.position.to_string() +
      " " + v.ref_allele + ">" + v.alt_allele,
    )
  }

  println("\n2. Motif break analysis")
  let results = @src.motif_break(variants, [("GATA1", pwm)], sequences)
  for r in results {
    let effect_str = match r.effect {
      @src.MotifbreakrEffect::Gain => "GAIN"
      @src.MotifbreakrEffect::Loss => "LOSS"
      @src.MotifbreakrEffect::Neutral => "NEUTRAL"
    }
    println(
      "   " + r.variant_name + " vs " + r.motif_name +
      ": ref=" + r.ref_score.to_string() +
      ", alt=" + r.alt_score.to_string() +
      ", delta=" + r.delta.to_string() +
      " -> " + effect_str,
    )
  }

  println("\n3. Significant effects (non-neutral)")
  let sig = @src.motifbreakr_significant(results)
  println("   " + sig.length().to_string() + " significant variant-motif pairs")

  println("\n4. Sorted by |delta| (most disruptive first)")
  let sorted = @src.motifbreakr_sort_by_delta(results)
  for r in sorted {
    let abs_d = if r.delta < 0.0 { 0.0 - r.delta } else { r.delta }
    println("   " + r.variant_name + ": |delta|=" + abs_d.to_string())
  }
}