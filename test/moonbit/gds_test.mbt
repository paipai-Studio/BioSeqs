///|
/// Tests for the SeqArray GDS container and SeqVarData filtering.

///|
/// Sample VCF text used to build GDS containers.
fn gds_test_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n" +
  "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n" +
  "1\t200\t.\tC\tA,T\t.\tPASS\t.\tGT\t0/1\t./.\t2/2\n" +
  "2\t150\trs3\tG\tC\t.\tPASS\t.\tGT\t1\t0/0\t0|1\n" +
  "2\t300\trs4\tT\tA\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
}

///|
/// Build a SeqVarData from the test VCF.
fn gds_test_data() -> @src.SeqVarData {
  let vcf = @src.parse_vcf(gds_test_vcf())
  @src.seq_vcf2gds(vcf)
}

///|
test "gds_file_node_tree" {
  let file = @src.gds_file_new()
  let _ = file.add_node(
    "annotation/info",
    @src.GdsNode::leaf("DP", "depth", @src.GdsValue::GdsIntArray([10, 20])),
  )
  let _ = file.add_node(
    "",
    @src.GdsNode::leaf("sample.id", "samples", @src.GdsValue::GdsStringArray(["A", "B"])),
  )
  assert_true(file.has_node("annotation/info/DP"))
  assert_true(file.has_node("sample.id"))
  assert_true(!file.has_node("annotation/missing"))
  let node = file.get_node("annotation/info/DP")
  assert_eq(node.name, "DP")
  match node.value {
    Some(@src.GdsValue::GdsIntArray(items)) => assert_eq(items, [10, 20])
    _ => assert_true(false)
  }
  // missing node raises
  let absent = try {
    let _ = file.get_node("nope")
    false
  } catch {
    _ => true
  }
  assert_true(absent)
}

///|
test "gds_parse_gt_tokens" {
  // haploid, diploid phased, unphased and missing calls all map to
  // (first, second) allele indices
  let vcf = @src.parse_vcf(gds_test_vcf())
  let data = @src.seq_vcf2gds(vcf)
  let dosages = data.genotype_dosages()
  // variant 2 (haploid 1 for S1, 0/0 for S2, 0|1 for S3)
  assert_eq(dosages[2][0], 1)
  assert_eq(dosages[2][1], 0)
  assert_eq(dosages[2][2], 1)
  // variant 1 has a missing call for S2
  assert_eq(dosages[1][1], -1)
  assert_eq(dosages[1][2], 2)
}

///|
test "gds_vcf_import_counts" {
  let data = gds_test_data()
  assert_eq(data.n_samples(), 3)
  assert_eq(data.n_variants_total(), 4)
  assert_eq(data.n_variants(), 4)
  assert_eq(data.sample_ids(), ["S1", "S2", "S3"])
  assert_eq(data.variant_ids(), ["rs1", "1:200:C", "rs3", "rs4"])
  assert_eq(data.chromosomes(), ["1", "1", "2", "2"])
  assert_eq(data.positions(), [100, 200, 150, 300])
  assert_eq(data.alleles(), ["A,G", "C,A,T", "G,C", "T,A"])
}

///|
test "gds_variant_index_filter" {
  let data = gds_test_data()
  data.set_variant_filter([1, 3])
  assert_eq(data.n_variants(), 2)
  assert_eq(data.variant_ids(), ["1:200:C", "rs4"])
  assert_eq(data.positions(), [200, 300])
  data.reset_filter()
  assert_eq(data.n_variants(), 4)
  // out of range raises
  let failed = try {
    data.set_variant_filter([9])
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "gds_variant_chrom_and_range_filters" {
  let data = gds_test_data()
  data.set_variant_chrom_filter(["2"])
  assert_eq(data.n_variants(), 2)
  assert_eq(data.positions(), [150, 300])
  data.reset_filter()
  data.set_variant_range_filter(100, 200)
  assert_eq(data.n_variants(), 2)
  assert_eq(data.positions(), [100, 150])
}

///|
test "gds_sample_filter" {
  let data = gds_test_data()
  data.set_sample_filter([0, 2])
  assert_eq(data.n_samples(), 2)
  assert_eq(data.sample_ids(), ["S1", "S3"])
  let dosages = data.genotype_dosages()
  assert_eq(dosages[0].length(), 2)
  assert_eq(dosages[0][0], 0)
  // S3 carries 1/1 at variant 0 -> dosage 2
  assert_eq(dosages[0][1], 2)
}

///|
test "gds_genotype_dosages" {
  let data = gds_test_data()
  let dosages = data.genotype_dosages()
  assert_eq(dosages.length(), 4)
  // variant 0: 0/0, 0/1, 1/1
  assert_eq(dosages[0], [0, 1, 2])
  // variant 1: 0/1, ., 2/2 (multi-allelic alt index 2 counts both copies)
  assert_eq(dosages[1], [1, -1, 2])
  // variant 3: 0/1, 1/1, 0/0
  assert_eq(dosages[3], [1, 2, 0])
}

///|
test "gds_allele_frequencies" {
  let data = gds_test_data()
  let freqs = data.allele_frequencies()
  assert_eq(freqs.length(), 4)
  assert_true((freqs[0] - 0.5).abs() < 1.0e-12)
  // variant 1: dosages 1 + 2 over 4 called alleles (missing excluded)
  assert_true((freqs[1] - 0.75).abs() < 1.0e-12)
  // variant 2: haploid 1 + 0/0 + 0|1 -> 2 alt over 5 called alleles
  assert_true((freqs[2] - 0.4).abs() < 1.0e-12)
}

///|
test "gds_missing_rates" {
  let data = gds_test_data()
  let rates = data.missing_rates()
  assert_eq(rates.length(), 4)
  assert_eq(rates[0], 0.0)
  assert_true((rates[1] - (1.0 / 3.0)).abs() < 1.0e-12)
}

///|
test "gds_summary" {
  let data = gds_test_data()
  let summary = data.summary()
  assert_eq(summary.n_samples, 3)
  assert_eq(summary.n_variants, 4)
  assert_eq(summary.n_chromosomes, 2)
  assert_true(summary.mean_allele_frequency > 0.0)
  assert_true(summary.mean_missing_rate > 0.0)
}

///|
test "gds_apply_per_variant" {
  let data = gds_test_data()
  let means = data.apply_per_variant(fn (row, _idx) {
    let mut total = 0.0
    let mut n = 0
    for dosage in row {
      if dosage >= 0 {
        total = total + dosage.to_double()
        n = n + 1
      }
    }
    if n > 0 {
      total / n.to_double()
    } else {
      0.0
    }
  })
  assert_eq(means.length(), 4)
  assert_true((means[0] - 1.0).abs() < 1.0e-12)
  assert_true((means[1] - 1.5).abs() < 1.0e-12)
}
