///|
/// Tests for Bio.PopGen.GenePop - GenePop file format parsing and analysis.

///|
test "gene_pop_parse_basic" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(record.title, "Sample GenePop File")
  assert_eq(record.populations.length(), 2)
  assert_eq(record.populations[0].individuals.length(), 3)
  assert_eq(record.populations[1].individuals.length(), 3)
}

///|
test "gene_pop_parse_empty" {
  let record = @src.parse_genepop("")
  assert_eq(record.title, "")
  assert_eq(record.populations.length(), 0)
}

///|
test "gene_pop_parse_loci_auto_generated" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  // 3 loci per individual -> auto-generated names Locus_1, Locus_2, Locus_3.
  assert_eq(record.loci_names.length(), 3)
  assert_eq(record.loci_names[0], "Locus_1")
  assert_eq(record.loci_names[1], "Locus_2")
  assert_eq(record.loci_names[2], "Locus_3")
}

///|
test "gene_pop_parse_loci_from_comment" {
  let content = "Test Title\n#Loci: LocusA LocusB\nPop\nind1, 0102 0304\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.loci_names.length(), 2)
  assert_eq(record.loci_names[0], "LocusA")
  assert_eq(record.loci_names[1], "LocusB")
}

///|
test "gene_pop_parse_digit_size_2" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(record.locus_digit_size, 2)
}

///|
test "gene_pop_parse_digit_size_3" {
  let content = "Title\nPop\nind1, 001002 003004\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.locus_digit_size, 3)
}

///|
test "gene_pop_parse_individual_names" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(record.populations[0].individuals[0].name, "ind1")
  assert_eq(record.populations[0].individuals[1].name, "ind2")
  assert_eq(record.populations[1].individuals[0].name, "ind4")
}

///|
test "gene_pop_parse_genotypes_diploid_2digit" {
  let content = "Title\nPop\nind1, 0102 0304\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes
  assert_eq(g.length(), 2)
  assert_eq(g[0].allele1, 1)
  assert_eq(g[0].allele2, Some(2))
  assert_eq(g[1].allele1, 3)
  assert_eq(g[1].allele2, Some(4))
}

///|
test "gene_pop_parse_genotypes_diploid_3digit" {
  let content = "Title\nPop\nind1, 001002 003004\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes
  assert_eq(g.length(), 2)
  assert_eq(g[0].allele1, 1)
  assert_eq(g[0].allele2, Some(2))
  assert_eq(g[1].allele1, 3)
  assert_eq(g[1].allele2, Some(4))
}

///|
test "gene_pop_parse_missing_data_2digit" {
  let content = "Title\nPop\nind1, 0000 0102\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes
  assert_eq(g[0].is_missing, true)
  assert_eq(g[1].is_missing, false)
}

///|
test "gene_pop_parse_missing_data_3digit" {
  let content = "Title\nPop\nind1, 000000 001002\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes
  assert_eq(g[0].is_missing, true)
  assert_eq(g[1].is_missing, false)
}

///|
test "gene_pop_parse_homozygous" {
  let content = "Title\nPop\nind1, 0101\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes[0]
  assert_eq(g.allele1, 1)
  assert_eq(g.allele2, Some(1))
}

///|
test "gene_pop_parse_pop_case_insensitive" {
  let content = "Title\nPOP\nind1, 0102\npop\nind2, 0304\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.populations.length(), 2)
}

///|
test "gene_pop_parse_multiple_populations" {
  let content = "Title\n" +
    "Pop\nind1, 0102\nind2, 0304\n" +
    "Pop\nind3, 0506\n" +
    "Pop\nind4, 0708\nind5, 0101\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.populations.length(), 3)
  assert_eq(record.populations[0].individuals.length(), 2)
  assert_eq(record.populations[1].individuals.length(), 1)
  assert_eq(record.populations[2].individuals.length(), 2)
}

///|
test "gene_pop_to_string_roundtrip" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  let text = @src.genepop_to_string(record)
  // Re-parse and verify structure.
  let record2 = @src.parse_genepop(text)
  assert_eq(record2.populations.length(), 2)
  assert_eq(record2.populations[0].individuals.length(), 3)
  assert_eq(record2.populations[1].individuals.length(), 3)
}

///|
test "gene_pop_total_individuals" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(@src.genepop_total_individuals(record), 6)
}

///|
test "gene_pop_num_populations" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(@src.genepop_num_populations(record), 2)
}

///|
test "gene_pop_num_loci" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  assert_eq(@src.genepop_num_loci(record), 3)
}

///|
test "gene_pop_alleles_at_locus" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  // Locus 0: ind1=01/02, ind2=02/02, ind3=01/01, ind4=01/02, ind5=02/02, ind6=01/01
  // Alleles: 1, 2
  let alleles = @src.genepop_alleles_at_locus(record, 0)
  assert_eq(alleles.length(), 2)
  assert_eq(alleles[0], 1)
  assert_eq(alleles[1], 2)
}

///|
test "gene_pop_allele_frequencies" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  // Population 0, locus 0: ind1=01/02, ind2=02/02, ind3=01/01
  // Allele counts: 01 appears 3 times, 02 appears 3 times. Total = 6.
  // Frequencies: 01=0.5, 02=0.5
  let pop = record.populations[0]
  let freqs = @src.genepop_allele_frequencies(pop, 0)
  let f1 = freqs.get_or_default(1, 0.0)
  let f2 = freqs.get_or_default(2, 0.0)
  assert_eq(f1, 0.5)
  assert_eq(f2, 0.5)
}

///|
test "gene_pop_observed_heterozygosity" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  // Population 0, locus 0: ind1=01/02 (het), ind2=02/02 (hom), ind3=01/01 (hom)
  // H_obs = 1/3
  let pop = record.populations[0]
  let h_obs = @src.genepop_observed_heterozygosity(pop, 0)
  // 1/3 ≈ 0.3333
  assert_eq(h_obs > 0.3, true)
  assert_eq(h_obs < 0.4, true)
}

///|
test "gene_pop_expected_heterozygosity" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  // Population 0, locus 0: p(01)=0.5, p(02)=0.5
  // H_exp = 1 - (0.5^2 + 0.5^2) = 1 - 0.5 = 0.5
  let pop = record.populations[0]
  let h_exp = @src.genepop_expected_heterozygosity(pop, 0)
  assert_eq(h_exp > 0.45, true)
  assert_eq(h_exp < 0.55, true)
}

///|
test "gene_pop_observed_heterozygosity_all_homozygous" {
  let content = "Title\nPop\nind1, 0101\nind2, 0202\nind3, 0303\n"
  let record = @src.parse_genepop(content)
  let pop = record.populations[0]
  let h_obs = @src.genepop_observed_heterozygosity(pop, 0)
  assert_eq(h_obs, 0.0)
}

///|
test "gene_pop_observed_heterozygosity_all_heterozygous" {
  let content = "Title\nPop\nind1, 0102\nind2, 0103\nind3, 0203\n"
  let record = @src.parse_genepop(content)
  let pop = record.populations[0]
  let h_obs = @src.genepop_observed_heterozygosity(pop, 0)
  assert_eq(h_obs, 1.0)
}

///|
test "gene_pop_summary_not_empty" {
  let sample = @src.genepop_sample()
  let record = @src.parse_genepop(sample)
  let s = @src.genepop_summary(record)
  assert_eq(s.length() > 0, true)
}

///|
test "gene_pop_sample_not_empty" {
  let sample = @src.genepop_sample()
  assert_eq(sample.length() > 0, true)
}

///|
test "gene_pop_parse_skip_blank_lines" {
  let content = "Title\n\nPop\n\nind1, 0102\n\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.populations.length(), 1)
  assert_eq(record.populations[0].individuals.length(), 1)
}

///|
test "gene_pop_parse_no_title" {
  let content = "Pop\nind1, 0102\n"
  let record = @src.parse_genepop(content)
  assert_eq(record.title, "")
  assert_eq(record.populations.length(), 1)
}

///|
test "gene_pop_genotype_diploid_constructor" {
  let g = @src.GenePopGenotype::diploid(1, 2)
  assert_eq(g.allele1, 1)
  assert_eq(g.allele2, Some(2))
  assert_eq(g.is_missing, false)
}

///|
test "gene_pop_genotype_haploid_constructor" {
  let g = @src.GenePopGenotype::haploid(5)
  assert_eq(g.allele1, 5)
  assert_eq(g.allele2, None)
  assert_eq(g.is_missing, false)
}

///|
test "gene_pop_genotype_missing_diploid" {
  let g = @src.GenePopGenotype::diploid(0, 0)
  assert_eq(g.is_missing, true)
}

///|
test "gene_pop_parse_3digit_high_allele" {
  let content = "Title\nPop\nind1, 100200\n"
  let record = @src.parse_genepop(content)
  let g = record.populations[0].individuals[0].genotypes[0]
  assert_eq(g.allele1, 100)
  assert_eq(g.allele2, Some(200))
}

///|
test "gene_pop_allele_frequencies_missing_excluded" {
  let content = "Title\nPop\nind1, 0102\nind2, 0000\nind3, 0101\n"
  let record = @src.parse_genepop(content)
  let pop = record.populations[0]
  // ind1=01/02, ind2=missing, ind3=01/01
  // Allele counts: 01 appears 3 times, 02 appears 1 time. Total = 4.
  // Frequencies: 01=0.75, 02=0.25
  let freqs = @src.genepop_allele_frequencies(pop, 0)
  let f1 = freqs.get_or_default(1, 0.0)
  let f2 = freqs.get_or_default(2, 0.0)
  assert_eq(f1, 0.75)
  assert_eq(f2, 0.25)
}

///|
test "gene_pop_parse_extra_whitespace" {
  let content = "Title\nPop\n   ind1  ,   0102   0304   \n"
  let record = @src.parse_genepop(content)
  assert_eq(record.populations[0].individuals[0].name, "ind1")
  assert_eq(record.populations[0].individuals[0].genotypes.length(), 2)
}
