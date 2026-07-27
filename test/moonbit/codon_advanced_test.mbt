///|
test "codon_advanced_codon_to_aa" {
  assert_eq(@src.codon_to_aa_map.get("ATG").unwrap_or(""), "M")
  assert_eq(@src.codon_to_aa_map.get("TGG").unwrap_or(""), "W")
  assert_eq(@src.codon_to_aa_map.get("TAA").unwrap_or(""), "*")
  assert_eq(@src.codon_to_aa_map.get("TAG").unwrap_or(""), "*")
  assert_eq(@src.codon_to_aa_map.get("TGA").unwrap_or(""), "*")
  assert_eq(@src.codon_to_aa_map.get("GCA").unwrap_or(""), "A")
  assert_eq(@src.codon_to_aa_map.get("GCG").unwrap_or(""), "A")
  assert_eq(@src.codon_to_aa_map.get("TTC").unwrap_or(""), "F")
  assert_eq(@src.codon_to_aa_map.get("TTT").unwrap_or(""), "F")
  assert_eq(@src.codon_to_aa_map.get("CTG").unwrap_or(""), "L")
}

///|
test "codon_advanced_aa_to_codons" {
  let codons = @src.aa_to_codons.get("A").unwrap_or([])
  assert_eq(codons.length(), 4)
  assert_eq(@src.aa_to_codons.get("M").unwrap_or([]).length(), 1)
  assert_eq(@src.aa_to_codons.get("W").unwrap_or([]).length(), 1)
  assert_eq(@src.aa_to_codons.get("L").unwrap_or([]).length(), 6)
  assert_eq(@src.aa_to_codons.get("S").unwrap_or([]).length(), 6)
  assert_eq(@src.aa_to_codons.get("*").unwrap_or([]).length(), 3)
}

///|
test "codon_advanced_usage_from_sequence" {
  let seq = "ATGGCCATGAAATGGTGA"
  let usage = @src.codon_usage_from_sequence(seq)

  assert_eq(@src.CodonUsageTable::total(usage), 6)
  assert_eq(@src.CodonUsageTable::get_count(usage, "ATG"), 2)
  assert_eq(@src.CodonUsageTable::get_count(usage, "GCC"), 1)
  assert_eq(@src.CodonUsageTable::get_count(usage, "AAA"), 1)
  assert_eq(@src.CodonUsageTable::get_count(usage, "TGG"), 1)
  assert_eq(@src.CodonUsageTable::get_count(usage, "TGA"), 1)
}

///|
test "codon_advanced_usage_from_sequence_empty" {
  let usage = @src.codon_usage_from_sequence("")
  assert_eq(@src.CodonUsageTable::total(usage), 0)
  assert_eq(@src.CodonUsageTable::get_count(usage, "ATG"), 0)
}

///|
test "codon_advanced_usage_from_sequence_short" {
  let usage = @src.codon_usage_from_sequence("AT")
  assert_eq(@src.CodonUsageTable::total(usage), 0)
}

///|
test "codon_advanced_usage_from_sequences" {
  let seqs = ["ATGGCCATG", "GCCGCCAAA", "ATGGCC"]
  let usage = @src.codon_usage_from_sequences(seqs)

  assert_eq(@src.CodonUsageTable::total(usage), 8)
  assert_eq(@src.CodonUsageTable::get_count(usage, "ATG"), 3)
  assert_eq(@src.CodonUsageTable::get_count(usage, "GCC"), 4)
  assert_eq(@src.CodonUsageTable::get_count(usage, "AAA"), 1)
}

///|
test "codon_advanced_usage_frequencies" {
  let seq = "ATGGCCATG"
  let usage = @src.codon_usage_from_sequence(seq)

  let freq_atg = @src.CodonUsageTable::get_frequency(usage, "ATG")
  let freq_gcc = @src.CodonUsageTable::get_frequency(usage, "GCC")
  assert_true(freq_atg > 0.0)
  assert_true(freq_gcc > 0.0)
  assert_true(@src.CodonUsageTable::get_frequency(usage, "XXX") == 0.0)
}

///|
test "codon_advanced_get_amino_acid" {
  assert_eq(@src.get_amino_acid("ATG"), "M")
  assert_eq(@src.get_amino_acid("TGG"), "W")
  assert_eq(@src.get_amino_acid("GCA"), "A")
  assert_eq(@src.get_amino_acid("TAA"), "*")
  assert_eq(@src.get_amino_acid("XXX"), "X")
  assert_eq(@src.get_amino_acid("atg"), "M")
}

///|
test "codon_advanced_get_codons_for_amino_acid" {
  let codons_a = @src.get_codons_for_amino_acid("A")
  assert_eq(codons_a.length(), 4)
  let codons_m = @src.get_codons_for_amino_acid("M")
  assert_eq(codons_m.length(), 1)
  assert_eq(codons_m[0], "ATG")
}

///|
test "codon_advanced_is_stop_codon" {
  assert_true(@src.is_stop_codon("TAA"))
  assert_true(@src.is_stop_codon("TAG"))
  assert_true(@src.is_stop_codon("TGA"))
  assert_false(@src.is_stop_codon("ATG"))
  assert_false(@src.is_stop_codon("GCA"))
}

///|
test "codon_advanced_is_start_codon" {
  assert_true(@src.is_start_codon("ATG"))
  assert_false(@src.is_start_codon("TTG"))
  assert_false(@src.is_start_codon("CTG"))
}

///|
test "codon_advanced_calculate_rscu" {
  let seq = "ATGGCCATGAAATGGTGA"
  let usage = @src.codon_usage_from_sequence(seq)
  let rscu = @src.calculate_rscu(usage)

  let val_atg = @src.RSCUResult::get_value(rscu, "ATG")
  assert_true(val_atg >= 0.0)

  let val_gcc = @src.RSCUResult::get_value(rscu, "GCC")
  assert_true(val_gcc >= 0.0)

  let values = @src.RSCUResult::values(rscu)
  assert_true(values.length() > 0)
}

///|
test "codon_advanced_calculate_rscu_equal_counts" {
  let counts = Map([], capacity=4)
  counts.set("GCA", 5)
  counts.set("GCC", 5)
  counts.set("GCG", 5)
  counts.set("GCT", 5)
  let usage = @src.CodonUsageTable::new(counts)
  let rscu = @src.calculate_rscu(usage)

  let val_gca = @src.RSCUResult::get_value(rscu, "GCA")
  let val_gcc = @src.RSCUResult::get_value(rscu, "GCC")
  assert_true(val_gca == 1.0)
  assert_true(val_gcc == 1.0)
}

///|
test "codon_advanced_calculate_cai" {
  let seq = "ATGGCCATGAAATGGTGA"
  let ref_seq = "ATGATGATGGCCGCCGCCAAAAAA"
  let ref_usage = @src.codon_usage_from_sequence(ref_seq)
  let cai_val = @src.calculate_cai(seq, ref_usage)

  assert_true(cai_val >= 0.0 && cai_val <= 1.0)
}

///|
test "codon_advanced_calculate_cai_empty" {
  let usage = @src.codon_usage_from_sequence("ATG")
  let cai_val = @src.calculate_cai("", usage)
  assert_eq(cai_val, 0.0)
}

///|
test "codon_advanced_calculate_cai_detailed" {
  let seq = "ATGGCCATGAAATGGTGA"
  let ref_seq = "ATGATGATGGCCGCCGCCAAAAAA"
  let ref_usage = @src.codon_usage_from_sequence(ref_seq)
  let result = @src.calculate_cai_detailed(seq, ref_usage)

  assert_true(@src.CAIResult::score(result) >= 0.0)
  assert_true(@src.CAIResult::codon_scores(result).length() > 0)
}

///|
test "codon_advanced_calculate_enc" {
  let seq = "ATGGCCATGAAATGGTGA"
  let usage = @src.codon_usage_from_sequence(seq)
  let enc_val = @src.calculate_enc(usage)

  assert_true(enc_val >= 20.0 && enc_val <= 61.0)
}

///|
test "codon_advanced_calculate_enc_low_bias" {
  let counts = Map([], capacity=24)
  counts.set("GCA", 1); counts.set("GCC", 1)
  counts.set("GCG", 1); counts.set("GCT", 1)
  counts.set("AAA", 1); counts.set("AAG", 1)
  counts.set("GAA", 1); counts.set("GAG", 1)
  counts.set("CTT", 1); counts.set("CTC", 1)
  counts.set("CTA", 1); counts.set("CTG", 1)
  counts.set("TTA", 1); counts.set("TTG", 1)
  counts.set("ATT", 1); counts.set("ATC", 1)
  counts.set("ATA", 1)
  counts.set("CGT", 1); counts.set("CGC", 1)
  counts.set("CGA", 1); counts.set("CGG", 1)
  counts.set("AGA", 1); counts.set("AGG", 1)
  let usage = @src.CodonUsageTable::new(counts)
  let enc_val = @src.calculate_enc(usage)

  assert_true(enc_val > 50.0)
}

///|
test "codon_advanced_calculate_gc3" {
  let seq = "ATGCCCAAA"
  let gc3 = @src.calculate_gc3(seq)
  assert_eq(gc3, 2.0 / 3.0)
}

///|
test "codon_advanced_calculate_gc3_empty" {
  let gc3 = @src.calculate_gc3("")
  assert_eq(gc3, 0.0)
}

///|
test "codon_advanced_optimal_codons" {
  let seq = "ATGGCCATGAAAGCCGCC"
  let usage = @src.codon_usage_from_sequence(seq)
  let optimal = @src.optimal_codons(usage)

  assert_true(optimal.length() > 0)
  let mut found_atg = false
  let mut found_gcc = false
  let mut i = 0
  while i < optimal.length() {
    if optimal[i].1 == "ATG" {
      found_atg = true
    }
    if optimal[i].1 == "GCC" {
      found_gcc = true
    }
    i = i + 1
  }
  assert_true(found_atg)
  assert_true(found_gcc)
}

///|
test "codon_advanced_rare_codons" {
  let counts = Map([], capacity=4)
  counts.set("GCA", 100)
  counts.set("GCC", 100)
  counts.set("GCG", 100)
  counts.set("GCT", 1)
  let usage = @src.CodonUsageTable::new(counts)
  let rare = @src.rare_codons(usage)

  assert_true(rare.length() > 0)
  let mut found_gct = false
  let mut i = 0
  while i < rare.length() {
    if rare[i] == "GCT" {
      found_gct = true
    }
    i = i + 1
  }
  assert_true(found_gct)
}

///|
test "codon_advanced_rare_codons_no_rare" {
  let counts = Map([], capacity=4)
  counts.set("GCA", 5)
  counts.set("GCC", 5)
  counts.set("GCG", 5)
  counts.set("GCT", 5)
  let usage = @src.CodonUsageTable::new(counts)
  let rare = @src.rare_codons(usage)

  assert_eq(rare.length(), 0)
}

///|
test "codon_advanced_bias_analysis" {
  let seq = "ATGGCCATGAAATGGTGAGCCGCCGCCAAAAAA"
  let result = @src.codon_bias_analysis(seq)

  assert_true(result.enc >= 20.0 && result.enc <= 61.0)
  assert_true(result.cai >= 0.0 && result.cai <= 1.0)
  assert_true(result.gc3 >= 0.0 && result.gc3 <= 1.0)
  assert_true(result.bias_score >= 0.0 && result.bias_score <= 1.0)
  assert_true(result.optimal_fraction >= 0.0 && result.optimal_fraction <= 1.0)
  assert_true(result.rare_fraction >= 0.0 && result.rare_fraction <= 1.0)
  assert_true(result.total_codons_analyzed > 0)
}

///|
test "codon_advanced_get_reference_table" {
  let high_gc = @src.get_reference_table("high_gc")
  assert_true(high_gc.length() > 0)
  assert_true(high_gc.contains("GCC"))

  let low_gc = @src.get_reference_table("low_gc")
  assert_true(low_gc.length() > 0)
  assert_true(low_gc.contains("GCT"))

  let general = @src.get_reference_table("general")
  assert_true(general.length() > 0)

  let unknown = @src.get_reference_table("unknown")
  assert_true(unknown.length() > 0)
}

///|
test "codon_advanced_calculate_cai_with_org" {
  let seq = "ATGGCCATGAAATGGTGA"
  let cai_high = @src.calculate_cai_with_org(seq, "high_gc")
  let cai_low = @src.calculate_cai_with_org(seq, "low_gc")

  assert_true(cai_high >= 0.0 && cai_high <= 1.0)
  assert_true(cai_low >= 0.0 && cai_low <= 1.0)
}

///|
test "codon_advanced_amino_acid_composition" {
  let seq = "ATGGCCATGAAATGGTGA"
  let usage = @src.codon_usage_from_sequence(seq)
  let comp = @src.get_amino_acid_composition(usage)

  assert_true(comp.length() > 0)
  assert_true(comp.contains("M"))
  assert_true(comp.contains("A"))
  assert_true(comp.contains("K"))
}

///|
test "codon_advanced_count_amino_acid_codons" {
  let seq = "ATGGCCATGAAATGGTGA"
  let count_m = @src.count_amino_acid_codons(seq, "M")
  let count_a = @src.count_amino_acid_codons(seq, "A")

  assert_eq(count_m, 2)
  assert_eq(count_a, 1)
}

///|
test "codon_advanced_at_third_fraction" {
  let seq = "ATGAAA"
  let frac = @src.at_third_fraction(seq)
  assert_eq(frac, 0.5)
}

///|
test "codon_advanced_gc_first_fraction" {
  let seq = "GCCGCA"
  let frac = @src.gc_first_fraction(seq)
  assert_eq(frac, 1.0)
}

///|
test "codon_advanced_gc_second_fraction" {
  let seq = "GCTGCA"
  let frac = @src.gc_second_fraction(seq)
  assert_eq(frac, 1.0)
}

///|
test "codon_advanced_codon_usage_to_string" {
  let seq = "ATGGCC"
  let usage = @src.codon_usage_from_sequence(seq)
  let str = @src.codon_usage_to_string(usage)

  assert_true(str.length() > 0)
  assert_true(str.contains("CodonUsageTable"))
  assert_true(@src.CodonUsageTable::total(usage) > 0)
}

///|
test "codon_advanced_compute_codon_adaptiveness" {
  let seq = "ATGGCCATGAAAGCCGCC"
  let usage = @src.codon_usage_from_sequence(seq)
  let w = @src.compute_codon_adaptiveness(usage)

  assert_true(w.length() > 0)
  let w_atg = w.get("ATG").unwrap_or(0.0)
  assert_true(w_atg >= 0.0 && w_atg <= 1.0)
}

///|
test "codon_advanced_multiple_synonymous" {
  let seq = "CTTCTCCTACTGCTTTTA"
  let usage = @src.codon_usage_from_sequence(seq)
  let rscu = @src.calculate_rscu(usage)

  assert_true(@src.RSCUResult::get_value(rscu, "CTT") > 0.0)
  assert_true(@src.RSCUResult::get_value(rscu, "CTC") > 0.0)
  assert_true(@src.RSCUResult::get_value(rscu, "CTG") > 0.0)
  assert_true(@src.RSCUResult::get_value(rscu, "TTA") > 0.0)
  assert_true(@src.RSCUResult::get_value(rscu, "CTA") > 0.0)
}

///|
test "codon_advanced_all_stop_codons" {
  assert_true(@src.is_stop_codon("TAA"))
  assert_true(@src.is_stop_codon("TAG"))
  assert_true(@src.is_stop_codon("TGA"))
  assert_eq(@src.get_amino_acid("TAA"), "*")
  assert_eq(@src.get_amino_acid("TAG"), "*")
  assert_eq(@src.get_amino_acid("TGA"), "*")
}