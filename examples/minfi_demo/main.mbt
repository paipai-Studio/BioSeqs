///|
fn main {
  println("=== minfi Demo ===")

  // 1. Create example RGChannelSet
  println("\n1. Create example RGChannelSet:")
  let rgset = @src.minfi_create_example_rgset()
  println(rgset.to_string())

  // 2. NOOB preprocessing
  println("\n2. preprocess_noob (background correction):")
  let mset_noob = @src.minfi_preprocess_noob(rgset)
  println(mset_noob.to_string())

  // 3. Illumina preprocessing
  println("\n3. preprocess_illumina:")
  let mset_illumina = @src.minfi_preprocess_illumina(rgset)
  println(mset_illumina.to_string())

  // 4. quantile normalization
  println("\n4. preprocess_quantile:")
  let mset_quantile = @src.minfi_preprocess_quantile(rgset)
  println(mset_quantile.to_string())

  // 5. funnorm functional normalization
  println("\n5. preprocess_funnorm:")
  let grset_funnorm = @src.minfi_preprocess_funnorm(rgset)
  println(grset_funnorm.to_string())

  // 6. ratio_convert to GenomicRatioSet
  println("\n6. ratio_convert (MethylSet -> GenomicRatioSet):")
  let grset = @src.minfi_ratio_convert(mset_noob)
  println(grset.to_string())

  // 7. get_beta and get_mvalues
  println("\n7. get_beta / get_mvalues:")
  let beta = @src.minfi_get_beta(grset)
  let mvalues = @src.minfi_get_mvalues(grset)
  println(
    "  Beta matrix: " +
    beta.length().to_string() +
    " probes x " +
    (if beta.length() > 0 {
      beta[0].length().to_string()
    } else {
      "0".to_string()
    }) +
    " samples",
  )
  println(
    "  M-value matrix: " +
    mvalues.length().to_string() +
    " probes x " +
    (if mvalues.length() > 0 {
      mvalues[0].length().to_string()
    } else {
      "0".to_string()
    }) +
    " samples",
  )

  // 8. detection_pvalue
  println("\n8. detection_pvalue:")
  let det_p = @src.minfi_detection_pvalue(rgset)
  println("  Computed " + det_p.length().to_string() + " detection p-values")
  let mut i = 0
  while i < det_p.length() && i < 5 {
    println("  " + rgset.probe_names[i] + ": " + det_p[i].to_string())
    i = i + 1
  }

  // 9. get_qc - QC report
  println("\n9. get_qc (quality control report):")
  let qc = @src.minfi_get_qc(rgset)
  println(qc.to_string())

  // 10. estimate_sex
  println("\n10. estimate_sex:")
  let sexes = @src.minfi_estimate_sex(grset)
  let mut j = 0
  while j < sexes.length() {
    println("  " + grset.sample_names[j] + ": " + sexes[j])
    j = j + 1
  }

  // 11. find_dmps - differential methylation analysis
  println("\n11. find_dmps (groups: 0,0,1,1):")
  let groups : Array[Int] = [0, 0, 1, 1]
  let dmp_result = @src.minfi_find_dmps(grset, groups)
  println(dmp_result.to_string())

  // 12. find_dmrs - differential methylation regions
  println("\n12. find_dmrs (max_gap=500):")
  let dmr_result = @src.minfi_find_dmrs(dmp_result, 500)
  println(dmr_result.to_string())

  // 13. filter_probes
  println("\n13. filter_probes (sd >= 0.01):")
  let filtered = @src.minfi_filter_probes(grset, 0.01)
  println(
    "  Before: " +
    grset.n_probes().to_string() +
    " probes, after: " +
    filtered.n_probes().to_string() +
    " probes",
  )

  // 14. drop_snp_probes
  println("\n14. drop_snp_probes:")
  let no_snp = @src.minfi_drop_snp_probes(grset)
  println(
    "  Before: " +
    grset.n_probes().to_string() +
    " probes, after: " +
    no_snp.n_probes().to_string() +
    " probes",
  )

  // 15. map_to_genome
  println("\n15. map_to_genome:")
  let loc_map = @src.minfi_map_to_genome(grset)
  let keys = loc_map.keys().collect()
  println("  Mapped " + keys.length().to_string() + " probes to genome")
  let mut k = 0
  while k < keys.length() && k < 5 {
    let probe = keys[k]
    let loc = loc_map.get(probe)
    if loc.is_some() {
      let (chr_s, pos_i) = loc.unwrap()
      println("  " + probe + " -> " + chr_s + ":" + pos_i.to_string())
    }
    k = k + 1
  }

  // Also demo with example GenomicRatioSet
  println("\n16. Bonus: example GenomicRatioSet summary:")
  let example_grset = @src.minfi_create_example_grset()
  println(example_grset.to_string())

  println("\nminfi demo completed!")
}
