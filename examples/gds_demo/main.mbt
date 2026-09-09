///|
/// Bioconductor SeqArray Demo - GDS-style variant data container.
///
/// Demonstrates the port of the Bioconductor `SeqArray` package
/// (Zheng et al.): VCF import via seqVCF2GDS, the hierarchical GDS
/// node tree, sample/variant filters (seqSetFilter semantics) and
/// dosage / allele-frequency / missing-rate summaries.

///|
/// A small multi-sample VCF with a haploid call, a multi-allelic site
/// and a missing genotype.
fn demo_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n" +
  "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n" +
  "1\t200\t.\tC\tA,T\t.\tPASS\t.\tGT\t0/1\t./.\t2/2\n" +
  "2\t150\trs3\tG\tC\t.\tPASS\t.\tGT\t1\t0/0\t0|1\n" +
  "2\t300\trs4\tT\tA\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
}

///|
fn mean(values : Array[Double]) -> Double {
  let mut total = 0.0
  for v in values {
    total = total + v
  }
  total / values.length().to_double()
}

///|
fn mean_int(values : Array[Int]) -> Double {
  let mut total = 0
  for v in values {
    total = total + v
  }
  total.to_double() / values.length().to_double()
}

///|
fn main raise {
  println("=== SeqArray: VCF -> GDS import (seqVCF2GDS) ===")
  let vcf = @src.parse_vcf(demo_vcf())
  let data = @src.seq_vcf2gds(vcf)
  println("samples  : \{data.sample_ids()}")
  println("variants : \{data.variant_ids()}")
  println("chroms   : \{data.chromosomes()}")
  println("positions: \{data.positions()}")
  println("alleles  : \{data.alleles()}")

  println("\n=== GDS node tree ===")
  let file = @src.gds_file_new()
  let _ = file.add_node(
    "",
    @src.GdsNode::leaf(
      "sample.id",
      "sample identifiers",
      @src.GdsValue::GdsStringArray(["S1", "S2", "S3"]),
    ),
  )
  let _ = file.add_node(
    "annotation/info/DP",
    @src.GdsNode::leaf("DP", "read depth", @src.GdsValue::GdsIntArray([10, 20, 30])),
  )
  println("has annotation/info/DP : \{file.has_node("annotation/info/DP")}")
  println("has sample.id          : \{file.has_node("sample.id")}")
  println("has annotation/missing : \{file.has_node("annotation/missing")}")

  println("\n=== Genotype dosages (filtered view) ===")
  let dosages = data.genotype_dosages()
  for v in 0..<dosages.length() {
    println("variant \{v}: \{dosages[v]}")
  }

  println("\n=== Per-variant summaries ===")
  let freqs = data.allele_frequencies()
  let rates = data.missing_rates()
  for v in 0..<freqs.length() {
    println(
      "variant \{v}: alt freq=\{freqs[v]}, missing rate=\{rates[v]}",
    )
  }

  println("\n=== seqSetFilter: keep only chromosome 2 ===")
  data.set_variant_chrom_filter(["2"])
  println("variants : \{data.variant_ids()}")
  println("positions: \{data.positions()}")
  let dos2 = data.genotype_dosages()
  for v in 0..<dos2.length() {
    println("variant \{v} dosages: \{dos2[v]}")
  }

  println("\n=== Sample filter: keep S1 and S3 ===")
  data.reset_filter()
  data.set_sample_filter([0, 2])
  println("samples  : \{data.sample_ids()}")
  println("dosages[0]: \{data.genotype_dosages()[0]}")

  println("\n=== Summary (seqSummary) ===")
  data.reset_filter()
  let summary = data.summary()
  println("n_samples = \{summary.n_samples}, n_variants = \{summary.n_variants}")
  println("n_chromosomes = \{summary.n_chromosomes}")
  println("mean allele frequency = \{summary.mean_allele_frequency}")
  println("mean missing rate     = \{summary.mean_missing_rate}")
  println("mean dosage per variant = \{mean(data.apply_per_variant(fn(row, _i) { mean_int(row) }))}")

  println("\n=== Range filter [100, 200) ===")
  data.set_variant_range_filter(100, 200)
  println("positions: \{data.positions()}")
}
