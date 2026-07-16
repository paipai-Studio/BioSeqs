///|
fn main {
  println("=== VariantAnnotation Demo ===")
  println("Inspired by Bioconductor VariantAnnotation package")
  println("")

  println("1. Detect variant types")
  println("  SNP: \{@bio.detect_variant_type("A", "T").show()}")
  println("  Insertion: \{@bio.detect_variant_type("A", "AT").show()}")
  println("  Deletion: \{@bio.detect_variant_type("AT", "A").show()}")
  println("  Indel: \{@bio.detect_variant_type("AT", "CG").show()}")
  println("")

  println("2. Create TxDb for gene annotations")
  let gene = @bio.Gene::new("gene1", "chr1", 100, 500, @bio.strand_plus())
  let tx = @bio.Transcript::new(
    "tx1",
    "gene1",
    "chr1",
    100,
    500,
    @bio.strand_plus(),
  )
  let exon = @bio.Exon::new(
    "exon1",
    "tx1",
    "gene1",
    "chr1",
    100,
    200,
    @bio.strand_plus(),
    1,
  )
  let cds = @bio.CDS::new(
    "cds1",
    "tx1",
    "gene1",
    "chr1",
    120,
    180,
    @bio.strand_plus(),
    1,
  )
  let txdb = @bio.TxDb::new([gene], [tx], [exon], [cds], "hg38", "Homo sapiens")
  println("TxDb created with 1 gene, 1 transcript, 1 exon, 1 CDS")
  println("")

  println("3. Locate variants in genes")
  let vcf_content1 = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t150\t.\tA\tT\t.\t.\t.\nchr1\t600\t.\tG\tC\t.\t.\t.\n"
  let vcf1 = @bio.parse_vcf(vcf_content1)
  let rec1 = vcf1.records[0]
  let rec2 = vcf1.records[1]
  let annotations = @bio.locate_variants([rec1, rec2], txdb)
  println("Found \{annotations.length()} variant annotations")
  let mut i = 0
  while i < annotations.length() {
    println(
      "  Variant at \{annotations[i].chrom}:\{annotations[i].pos}: \{annotations[i].ref_allele} -> \{annotations[i].alt_allele} in \{annotations[i].gene_id}",
    )
    i = i + 1
  }
  println("")

  println("4. Predict coding effects")
  let cds_short = @bio.CDS::new(
    "cds1",
    "tx1",
    "gene1",
    "chr1",
    100,
    106,
    @bio.strand_plus(),
    1,
  )
  let txdb_short = @bio.TxDb::new(
    [gene],
    [tx],
    [],
    [cds_short],
    "hg38",
    "Homo sapiens",
  )
  let vcf_content2 = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\t.\tT\tC\t.\t.\t.\n"
  let vcf2 = @bio.parse_vcf(vcf_content2)
  let rec_coding = vcf2.records[0]
  let predictions = @bio.predict_coding(
    [rec_coding],
    txdb_short,
    @hashmap.HashMap([]),
  )
  println("Found \{predictions.length()} coding predictions")
  i = 0
  while i < predictions.length() {
    let pred = predictions[i]
    println(
      "  \{pred.chrom}:\{pred.pos}: \{pred.ref_allele}->\{pred.alt_allele}",
    )
    println("    Gene: \{pred.gene_id}, Tx: \{pred.transcript_id}")
    println("    Codon: \{pred.original_codon} -> \{pred.mutated_codon}")
    println(
      "    AA: \{pred.original_aa} -> \{pred.mutated_aa} (\{pred.aa_change})",
    )
    println("    Consequence: \{pred.consequence}")
    i = i + 1
  }
  println("")

  println("5. Generate variant summary")
  let vcf_content3 = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t150\t.\tG\tGA\t.\t.\t.\nchr1\t200\t.\tCT\tC\t.\t.\t.\n"
  let vcf3 = @bio.parse_vcf(vcf_content3)
  let rec_snp = vcf3.records[0]
  let rec_ins = vcf3.records[1]
  let rec_del = vcf3.records[2]
  let va1 = @bio.VariantAnnotation::new(
    "chr1",
    100,
    "A",
    "T",
    @bio.detect_variant_type("A", "T"),
    "gene1",
    "gene1",
    "tx1",
    "exon1",
    0,
    "",
    "",
    "exon",
  )
  let va2 = @bio.VariantAnnotation::new(
    "chr1",
    150,
    "G",
    "GA",
    @bio.detect_variant_type("G", "GA"),
    "gene1",
    "gene1",
    "tx1",
    "exon1",
    0,
    "",
    "",
    "coding",
  )
  let summary = @bio.variant_summary([rec_snp, rec_ins, rec_del], [va1, va2])
  println("Total variants: \{summary.total_variants}")
  println("SNPs: \{summary.snp_count}")
  println("Insertions: \{summary.insertion_count}")
  println("Deletions: \{summary.deletion_count}")
  println("Coding variants: \{summary.coding_variants}")
  println("Non-coding variants: \{summary.non_coding_variants}")
}
