///|
test "variant_annotation_detect_variant_type" {
  assert_true(@src.detect_variant_type("A", "T") is @src.VariantType::SNP)
  assert_true(
    @src.detect_variant_type("A", "AT") is @src.VariantType::Insertion,
  )
  assert_true(@src.detect_variant_type("AT", "A") is @src.VariantType::Deletion)
  assert_true(@src.detect_variant_type("AT", "CG") is @src.VariantType::Indel)
}

///|
test "variant_annotation_locate_variants" {
  let gene = @src.Gene::new("gene1", "chr1", 100, 500, @src.strand_plus())
  let tx = @src.Transcript::new(
    "tx1",
    "gene1",
    "chr1",
    100,
    500,
    @src.strand_plus(),
  )
  let exon = @src.Exon::new(
    "exon1",
    "tx1",
    "gene1",
    "chr1",
    100,
    200,
    @src.strand_plus(),
    1,
  )
  let cds = @src.CDS::new(
    "cds1",
    "tx1",
    "gene1",
    "chr1",
    120,
    180,
    @src.strand_plus(),
    1,
  )
  let txdb = @src.TxDb::new([gene], [tx], [exon], [cds], "hg38", "Homo sapiens")

  let vcf_content = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t150\t.\tA\tT\t.\t.\t.\nchr1\t600\t.\tG\tC\t.\t.\t.\n"
  let vcf = @src.parse_vcf(vcf_content)
  let rec1 = vcf.records[0]
  let rec2 = vcf.records[1]

  let annotations = @src.locate_variants([rec1, rec2], txdb)
  assert_eq(annotations.length(), 2)
}

///|
test "variant_annotation_predict_coding" {
  let gene = @src.Gene::new("gene1", "chr1", 100, 500, @src.strand_plus())
  let tx = @src.Transcript::new(
    "tx1",
    "gene1",
    "chr1",
    100,
    500,
    @src.strand_plus(),
  )
  let cds = @src.CDS::new(
    "cds1",
    "tx1",
    "gene1",
    "chr1",
    100,
    106,
    @src.strand_plus(),
    1,
  )
  let txdb = @src.TxDb::new([gene], [tx], [], [cds], "hg38", "Homo sapiens")

  let vcf_content = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\t.\tT\tC\t.\t.\t.\n"
  let vcf = @src.parse_vcf(vcf_content)
  let rec = vcf.records[0]

  let predictions = @src.predict_coding([rec], txdb, @hashmap.HashMap([]))
  assert_eq(predictions.length(), 0)
}

///|
test "variant_annotation_variant_summary" {
  let vcf_content = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t150\t.\tG\tGA\t.\t.\t.\nchr1\t200\t.\tCT\tC\t.\t.\t.\n"
  let vcf = @src.parse_vcf(vcf_content)
  let rec1 = vcf.records[0]
  let rec2 = vcf.records[1]
  let rec3 = vcf.records[2]

  let va1 = @src.VariantAnnotation::new(
    "chr1",
    100,
    "A",
    "T",
    @src.detect_variant_type("A", "T"),
    "gene1",
    "gene1",
    "tx1",
    "exon1",
    0,
    "",
    "",
    "exon",
  )
  let va2 = @src.VariantAnnotation::new(
    "chr1",
    150,
    "G",
    "GA",
    @src.detect_variant_type("G", "GA"),
    "gene1",
    "gene1",
    "tx1",
    "exon1",
    0,
    "",
    "",
    "coding",
  )

  let summary = @src.variant_summary([rec1, rec2, rec3], [va1, va2])
  assert_eq(summary.total_variants, 3)
  assert_eq(summary.snp_count, 1)
  assert_eq(summary.insertion_count, 1)
  assert_eq(summary.deletion_count, 1)
  assert_eq(summary.coding_variants, 1)
  assert_eq(summary.non_coding_variants, 1)
}

///|
test "variant_annotation_variant_summary_with_coding" {
  let vcf_content = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t101\t.\tT\tA\t.\t.\t.\nchr1\t102\t.\tG\tA\t.\t.\t.\n"
  let vcf = @src.parse_vcf(vcf_content)
  let rec1 = vcf.records[0]
  let rec2 = vcf.records[1]
  let rec3 = vcf.records[2]

  let cp1 = @src.CodingPrediction::new(
    "chr1", 100, "A", "T", "gene1", "tx1", 1, 0, "ATG", "TTG", "M", "L", "missense",
    "missense_variant",
  )
  let cp2 = @src.CodingPrediction::new(
    "chr1", 101, "T", "A", "gene1", "tx1", 2, 1, "ATG", "AAG", "M", "K", "missense",
    "missense_variant",
  )
  let cp3 = @src.CodingPrediction::new(
    "chr1", 102, "G", "A", "gene1", "tx1", 3, 2, "ATG", "ATA", "M", "I", "missense",
    "missense_variant",
  )

  let summary = @src.variant_summary_with_coding([rec1, rec2, rec3], [
    cp1, cp2, cp3,
  ])
  assert_eq(summary.total_variants, 3)
  assert_eq(summary.snp_count, 3)
  assert_eq(summary.missense_count, 3)
}

///|
test "variant_annotation_translate_codon" {
  assert_eq(@src.translate_codon("ATG"), "M")
  assert_eq(@src.translate_codon("TAA"), "*")
  assert_eq(@src.translate_codon("GCC"), "A")
  assert_eq(@src.translate_codon("TTT"), "F")
}
