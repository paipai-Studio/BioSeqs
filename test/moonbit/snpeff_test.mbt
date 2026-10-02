///|
/// Helper: build a simple gene with one transcript for testing.
fn build_test_annotation() -> @src.GenomeAnnotation {
  let ann = @src.GenomeAnnotation::new("test", "1.0", "test")
  let gene = @src.GFGene::new("g1", "GENE1", "protein_coding", "chr1", 100, 500, "+", "")
  let tx = @src.GFTranscript::new("t1", "TX1", "g1", "chr1", 100, 500, "+")
  tx.exons.push(@src.GFExon::new("e1", 1, "chr1", 100, 500, "+"))
  tx.cds.push(@src.GFExon::new("cds1", 1, "chr1", 103, 498, "+"))
  tx.utr5.push(@src.GFExon::new("utr5_1", 1, "chr1", 100, 102, "+"))
  tx.utr3.push(@src.GFExon::new("utr3_1", 1, "chr1", 499, 500, "+"))
  gene.transcripts.push(tx)
  ann.genes.push(gene)
  ann
}

///|
/// Helper: build reference sequence for testing.
fn build_test_reference() -> String {
  let mut seq = "AAA".repeat(34) // 102 bp 5'UTR
  seq += "ATG" // 103-105 start codon
  seq += "GCT".repeat(10) // 106-135 Ala x10
  seq += "GAA".repeat(10) // 136-165 Glu x10
  seq += "TCT".repeat(10) // 166-195 Ser x10
  seq += "CGT".repeat(10) // 196-225 Arg x10
  seq += "GGT".repeat(10) // 226-255 Gly x10
  seq += "CCT".repeat(10) // 256-285 Pro x10
  seq += "ACT".repeat(10) // 286-315 Thr x10
  seq += "TGT".repeat(10) // 316-345 Cys x10
  seq += "GTT".repeat(10) // 346-375 Val x10
  seq += "ATT".repeat(10) // 376-405 Ile x10
  seq += "TTT".repeat(10) // 406-435 Phe x10
  seq += "CAT".repeat(10) // 436-465 His x10
  seq += "AAC".repeat(10) // 466-495 Asn x10
  seq += "TAA" // 496-498 stop codon
  seq += "AA" // 499-500 3'UTR
  seq
}

///|
/// Helper: build a VCF from text.
fn make_vcf(text : String) -> @src.VcfFile {
  @src.parse_vcf(text)
}

///|
/// Helper: make a single-variant VCF.
fn make_single_variant_vcf(chrom : String, pos : Int, ref_ : String, alt : String) -> @src.VcfFile {
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    chrom + "\t" + pos.to_string() + "\t.\t" + ref_ + "\t" + alt + "\t.\t.\t.\n"
  make_vcf(text)
}

///|
test "snpeff_so_string" {
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::MissenseVariant), "missense_variant")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::StopGained), "stop_gained")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::FrameshiftVariant), "frameshift_variant")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::SpliceDonorVariant), "splice_donor_variant")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::SynonymousVariant), "synonymous_variant")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::UpstreamVariant), "upstream_gene_variant")
  assert_eq(@src.snpeff_so_string(@src.SnpEffSO::FivePrimeUTRVariant), "5_prime_UTR_variant")
}

///|
test "snpeff_impact_string" {
  assert_eq(@src.snpeff_impact_string(@src.SnpEffImpact::High), "HIGH")
  assert_eq(@src.snpeff_impact_string(@src.SnpEffImpact::Moderate), "MODERATE")
  assert_eq(@src.snpeff_impact_string(@src.SnpEffImpact::Low), "LOW")
  assert_eq(@src.snpeff_impact_string(@src.SnpEffImpact::Modifier), "MODIFIER")
}

///|
test "snpeff_so_impact mapping" {
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::StopGained) is @src.SnpEffImpact::High)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::StopLost) is @src.SnpEffImpact::High)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::StartLost) is @src.SnpEffImpact::High)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::FrameshiftVariant) is @src.SnpEffImpact::High)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::SpliceDonorVariant) is @src.SnpEffImpact::High)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::MissenseVariant) is @src.SnpEffImpact::Moderate)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::SynonymousVariant) is @src.SnpEffImpact::Low)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::UpstreamVariant) is @src.SnpEffImpact::Modifier)
  assert_true(@src.snpeff_so_impact(@src.SnpEffSO::IntergenicVariant) is @src.SnpEffImpact::Modifier)
}

///|
test "snpeff predict missense" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 107, "C", "T", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_missense = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::MissenseVariant {
      has_missense = true
      assert_true(eff.hgvs_c.contains("c."))
      assert_true(eff.hgvs_p.contains("p."))
    }
  }
  assert_true(has_missense)
}

///|
test "snpeff predict synonymous" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 108, "T", "C", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_synonymous = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::SynonymousVariant {
      has_synonymous = true
    }
  }
  assert_true(has_synonymous)
}

///|
test "snpeff predict stop_gained" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 136, "G", "T", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_stop = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::StopGained {
      has_stop = true
    }
  }
  assert_true(has_stop)
}

///|
test "snpeff predict start_lost" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 103, "A", "C", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_start_lost = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::StartLost {
      has_start_lost = true
    }
  }
  assert_true(has_start_lost)
}

///|
test "snpeff predict frameshift" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 110, "A", "AT", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_frameshift = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::FrameshiftVariant {
      has_frameshift = true
    }
  }
  assert_true(has_frameshift)
}

///|
test "snpeff predict 5prime_UTR" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 101, "A", "G", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_utr5 = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::FivePrimeUTRVariant {
      has_utr5 = true
    }
  }
  assert_true(has_utr5)
}

///|
test "snpeff predict 3prime_UTR" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 499, "A", "G", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_utr3 = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::ThreePrimeUTRVariant {
      has_utr3 = true
    }
  }
  assert_true(has_utr3)
}

///|
test "snpeff predict intergenic downstream" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 1000, "A", "G", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_downstream = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::DownstreamVariant {
      has_downstream = true
      assert_true(eff.distance > 0)
    }
  }
  assert_true(has_downstream)
}

///|
test "snpeff predict upstream" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 50, "A", "G", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_upstream = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::UpstreamVariant {
      has_upstream = true
      assert_true(eff.distance > 0)
    }
  }
  assert_true(has_upstream)
}

///|
test "snpeff HGVS c. notation" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 107, "C", "T", ann, ref_seq)
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::MissenseVariant {
      assert_true(eff.hgvs_c.contains("c.5"))
      assert_true(eff.hgvs_c.contains("C>T"))
    }
  }
}

///|
test "snpeff HGVS p. notation" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 107, "C", "T", ann, ref_seq)
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::MissenseVariant {
      assert_true(eff.hgvs_p.contains("p."))
      assert_true(eff.hgvs_p.contains("Ala"))
      assert_true(eff.hgvs_p.contains("Val"))
    }
  }
}

///|
test "snpeff predict stop_lost" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 496, "T", "C", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let mut has_stop_lost = false
  for eff in effects {
    if eff.so_term is @src.SnpEffSO::StopLost {
      has_stop_lost = true
    }
  }
  assert_true(has_stop_lost)
}

///|
test "snpeff ANN field format" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 107, "C", "T", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let ann_str = @src.snpeff_to_ann(effects[0])
  // ANN should have 16 fields separated by |
  let fields = ann_str.split("|").to_array()
  assert_eq(fields.length(), 16)
  assert_eq(fields[0], "T")
  assert_true(fields[1].length() > 0)
  assert_true(fields[2].length() > 0)
}

///|
test "snpeff ANN field multiple effects" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let effects = @src.snpeff_predict_variant("chr1", 107, "C", "T", ann, ref_seq)
  assert_true(effects.length() >= 1)
  let ann_str = @src.snpeff_to_ann_field(effects)
  assert_true(ann_str.contains("missense_variant"))
}

///|
test "snpeff_split_multiallelic" {
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tT,G,C\t.\t.\t.\n"
  let vcf = make_vcf(text)
  let split = @src.snpeff_split_multiallelic(vcf)
  assert_eq(split.records.length(), 3)
  assert_eq(split.records[0].alt.length(), 1)
  assert_eq(split.records[1].alt.length(), 1)
  assert_eq(split.records[2].alt.length(), 1)
}

///|
test "snpeff end-to-end annotation" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t107\t.\tC\tT\t.\t.\t.\n" +
    "chr1\t1000\t.\tA\tG\t.\t.\t.\n"
  let vcf = make_vcf(text)
  let ref_map : Map[String, String] = Map([], capacity=0)
  ref_map["chr1"] = ref_seq
  let result = @src.snpeff_annotate(vcf, ann, ref_map)
  assert_eq(result.num_variants, 2)
  assert_true(result.num_annotated >= 1)
  assert_true(result.summary.size() > 0)
}

///|
test "snpeff summary text" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t107\t.\tC\tT\t.\t.\t.\n"
  let vcf = make_vcf(text)
  let ref_map : Map[String, String] = Map([], capacity=0)
  ref_map["chr1"] = ref_seq
  let result = @src.snpeff_annotate(vcf, ann, ref_map)
  let summary = @src.snpeff_summary_text(result)
  assert_true(summary.contains("SnpEff annotation summary"))
  assert_true(summary.contains("Total variants: 1"))
}

///|
test "snpeff count by impact" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t107\t.\tC\tT\t.\t.\t.\n" +
    "chr1\t1000\t.\tA\tG\t.\t.\t.\n"
  let vcf = make_vcf(text)
  let ref_map : Map[String, String] = Map([], capacity=0)
  ref_map["chr1"] = ref_seq
  let result = @src.snpeff_annotate(vcf, ann, ref_map)
  let counts = @src.snpeff_count_by_impact(result)
  assert_true(counts.size() > 0)
}

///|
test "snpeff annotate VCF output" {
  let ann = build_test_annotation()
  let ref_seq = build_test_reference()
  let text =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t107\t.\tC\tT\t.\t.\t.\n"
  let vcf = make_vcf(text)
  let ref_map : Map[String, String] = Map([], capacity=0)
  ref_map["chr1"] = ref_seq
  let output = @src.snpeff_annotate_vcf(vcf, ann, ref_map)
  assert_true(output.contains("##fileformat=VCFv4.2"))
  assert_true(output.contains("ANN="))
  assert_true(output.contains("missense_variant"))
}
