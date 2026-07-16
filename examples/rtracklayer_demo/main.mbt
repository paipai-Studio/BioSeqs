///|
fn main {
  println("=== rtracklayer Demo ===")
  println("")

  println("1. Parse BED format")
  let bed_content = "chr1\t100\t200\tgene1\t500\t+\nchr2\t300\t400\tgene2\t700\t-\nchr1\t500\t600\tgene3\t300\t."
  let bed_records = @bio.parse_bed(bed_content)
  
  let mut i = 0
  while i < bed_records.length() {
    let r = bed_records[i]
    println("  \{r.chrom}:\{r.chrom_start}-\{r.chrom_end} \{r.name}")
    i = i + 1
  }
  println("")

  println("2. BED to GRanges")
  let gr_from_bed = @bio.bed_to_granges(bed_records)
  i = 0
  while i < gr_from_bed.seqnames.length() {
    println("  \{gr_from_bed.seqnames[i]}:\{gr_from_bed.starts[i]}-\{gr_from_bed.ends[i]}")
    i = i + 1
  }
  println("")

  println("3. Parse WIG format (variableStep)")
  let wig_content = "variableStep chrom=chr1 span=5\n100\t1.5\n105\t2.0\n110\t1.8"
  let wig_records = @bio.parse_wig(wig_content)
  
  i = 0
  while i < wig_records.length() {
    let r = wig_records[i]
    let mut j = 0
    while j < r.positions.length() {
      println("  pos=\{r.positions[j]}, value=\{r.values[j]}")
      j = j + 1
    }
    i = i + 1
  }
  println("")

  println("4. Parse WIG format (fixedStep)")
  let wig_fixed = "fixedStep chrom=chr1 start=100 step=10\n1.5\n2.0\n1.8"
  let wig_fixed_records = @bio.parse_wig(wig_fixed)
  
  i = 0
  while i < wig_fixed_records.length() {
    let r = wig_fixed_records[i]
    let mut j = 0
    while j < r.positions.length() {
      println("  pos=\{r.positions[j]}, value=\{r.values[j]}")
      j = j + 1
    }
    i = i + 1
  }
  println("")

  println("5. Parse BEDGraph format")
  let bg_content = "track type=bedGraph\nchr1\t100\t200\t1.5\nchr1\t200\t300\t2.0\nchr2\t100\t150\t0.5"
  let bg_records = @bio.parse_bedgraph(bg_content)
  
  i = 0
  while i < bg_records.length() {
    let r = bg_records[i]
    println("  \{r.chrom}:\{r.chrom_start}-\{r.chrom_end} value=\{r.value}")
    i = i + 1
  }
  println("")

  println("6. Parse GFF format")
  let gff_content = "##gff-version 3\nchr1\tensembl\tgene\t100\t200\t5.5\t+\t0\tID=gene1;Name=GeneA\nchr2\tensembl\tgene\t300\t400\t.\t-\t1\tID=gene2;Name=GeneB"
  let gff_records = @bio.parse_gff(gff_content)
  
  i = 0
  while i < gff_records.length() {
    let r = gff_records[i]
    println("  \{r.seqid}:\{r.start}-\{r.end} \{r.feature_type} \{r.source}")
    i = i + 1
  }
  println("")

  println("7. GFF to GRanges")
  let gr_from_gff = @bio.gff_to_granges(gff_records)
  i = 0
  while i < gr_from_gff.seqnames.length() {
    println("  \{gr_from_gff.seqnames[i]}:\{gr_from_gff.starts[i]}-\{gr_from_gff.ends[i]}")
    i = i + 1
  }
  println("")

  println("8. Write BED format")
  let bed_output = @bio.write_bed(bed_records)
  println("BED output:")
  println(bed_output)
  println("")

  println("9. Write WIG format")
  let wig_output = @bio.write_wig(wig_records)
  println("WIG output:")
  println(wig_output)
  println("")

  println("10. Write BEDGraph format")
  let bg_output = @bio.write_bedgraph(bg_records)
  println("BEDGraph output:")
  println(bg_output)
  println("")

  println("11. Write GFF format")
  let gff_output = @bio.write_gff(gff_records)
  println("GFF output:")
  println(gff_output)
}