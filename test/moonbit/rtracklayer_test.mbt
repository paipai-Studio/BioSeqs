///|
test "rtracklayer_parse_bed_minimal" {
  let bed_content = "chr1\t100\t200\nchr2\t300\t400"
  let records = @src.parse_bed(bed_content)
  
  assert_eq(records.length(), 2)
  assert_eq(records[0].chrom, "chr1")
  assert_eq(records[0].chrom_start, 100)
  assert_eq(records[0].chrom_end, 200)
  assert_eq(records[0].name, "")
}

///|
test "rtracklayer_parse_bed_full" {
  let bed_content = "chr1\t100\t200\tgene1\t500\t+\t120\t180\t255,0,0\t2\t10,20\t0,50"
  let records = @src.parse_bed(bed_content)
  
  assert_eq(records.length(), 1)
  assert_eq(records[0].chrom, "chr1")
  assert_eq(records[0].chrom_start, 100)
  assert_eq(records[0].chrom_end, 200)
  assert_eq(records[0].name, "gene1")
  assert_true(records[0].score is Some(_))
}

///|
test "rtracklayer_parse_bed_with_strand" {
  let bed_content = "chr1\t100\t200\t\t\t+\nchr2\t300\t400\t\t\t-\nchr3\t500\t600"
  let records = @src.parse_bed(bed_content)
  
  assert_eq(records.length(), 3)
}

///|
test "rtracklayer_parse_bed_skip_comments" {
  let bed_content = "# This is a comment\nchr1\t100\t200\n# Another comment\nchr2\t300\t400"
  let records = @src.parse_bed(bed_content)
  
  assert_eq(records.length(), 2)
}

///|
test "rtracklayer_write_bed" {
  let bed_content = "chr1\t100\t200\tgene1\t500\t+\t120\t180\t255,0,0"
  let records = @src.parse_bed(bed_content)
  let output = @src.write_bed(records)
  
  assert_true(output.length() > 0)
}

///|
test "rtracklayer_bed_to_granges" {
  let bed_content = "chr1\t100\t200\t\t\t+\nchr2\t300\t400\t\t\t-"
  let records = @src.parse_bed(bed_content)
  let gr = @src.bed_to_granges(records)
  
  assert_eq(gr.seqnames.length(), 2)
  assert_eq(gr.starts[0], 100)
  assert_eq(gr.ends[0], 200)
}

///|
test "rtracklayer_parse_wig_variable_step" {
  let wig_content = "variableStep chrom=chr1 span=5\n100\t1.5\n105\t2.0\n110\t1.8"
  let records = @src.parse_wig(wig_content)
  
  assert_eq(records.length(), 1)
  assert_eq(records[0].positions.length(), 3)
  assert_eq(records[0].values.length(), 3)
}

///|
test "rtracklayer_parse_wig_fixed_step" {
  let wig_content = "fixedStep chrom=chr1 start=100 step=10\n1.5\n2.0\n1.8"
  let records = @src.parse_wig(wig_content)
  
  assert_eq(records.length(), 1)
  assert_eq(records[0].positions.length(), 3)
}

///|
test "rtracklayer_write_wig" {
  let wig_content = "variableStep chrom=chr1\n100\t1.5\n105\t2.0"
  let records = @src.parse_wig(wig_content)
  let output = @src.write_wig(records)
  
  assert_true(output.length() > 0)
}

///|
test "rtracklayer_parse_bedgraph" {
  let bg_content = "track type=bedGraph\nchr1\t100\t200\t1.5\nchr1\t200\t300\t2.0\nchr2\t100\t150\t0.5"
  let records = @src.parse_bedgraph(bg_content)
  
  assert_eq(records.length(), 3)
  assert_eq(records[0].chrom, "chr1")
  assert_eq(records[0].chrom_start, 100)
  assert_eq(records[0].chrom_end, 200)
}

///|
test "rtracklayer_write_bedgraph" {
  let bg_content = "track type=bedGraph\nchr1\t100\t200\t1.5\nchr2\t300\t400\t2.0"
  let records = @src.parse_bedgraph(bg_content)
  let output = @src.write_bedgraph(records)
  
  assert_true(output.length() > 0)
}

///|
test "rtracklayer_parse_gff" {
  let gff_content = "##gff-version 3\nchr1\tsource\tgene\t100\t200\t.\t+\t0\tID=gene1;Name=GeneA"
  let records = @src.parse_gff(gff_content)
  
  assert_eq(records.length(), 1)
  assert_eq(records[0].seqid, "chr1")
  assert_eq(records[0].source, "source")
  assert_eq(records[0].feature_type, "gene")
  assert_eq(records[0].start, 100)
  assert_eq(records[0].end, 200)
}

///|
test "rtracklayer_parse_gff_attributes" {
  let gff_content = "chr1\tsource\tgene\t100\t200\t5.5\t-\t1\tID=gene1;Name=GeneA"
  let records = @src.parse_gff(gff_content)
  
  assert_eq(records.length(), 1)
}

///|
test "rtracklayer_write_gff" {
  let gff_content = "chr1\tsource\tgene\t100\t200\t.\t+\t0\tID=gene1"
  let records = @src.parse_gff(gff_content)
  let output = @src.write_gff(records)
  
  assert_true(output.length() > 0)
}

///|
test "rtracklayer_gff_to_granges" {
  let gff_content = "chr1\tsource\tgene\t100\t200\t.\t+\t0\nchr2\tsource\tgene\t300\t400\t.\t-\t0"
  let records = @src.parse_gff(gff_content)
  let gr = @src.gff_to_granges(records)
  
  assert_eq(gr.seqnames.length(), 2)
  assert_eq(gr.starts[0], 100)
  assert_eq(gr.ends[0], 200)
}

///|
test "rtracklayer_parse_bed_empty" {
  let records = @src.parse_bed("")
  assert_eq(records.length(), 0)
}

///|
test "rtracklayer_parse_wig_empty" {
  let records = @src.parse_wig("")
  assert_eq(records.length(), 0)
}

///|
test "rtracklayer_parse_bedgraph_empty" {
  let records = @src.parse_bedgraph("")
  assert_eq(records.length(), 0)
}

///|
test "rtracklayer_parse_gff_empty" {
  let records = @src.parse_gff("")
  assert_eq(records.length(), 0)
}