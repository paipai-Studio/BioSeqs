///|
/// SAM & VCF parsing example.
///
/// Demonstrates parsing SAM alignment files and VCF variant files,
/// accessing records, and querying variant properties.
///
/// Run: moon run examples/sam_vcf_demo/main.mbt
fn main {
  // === Parse SAM ===
  let sam_content = "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:249250621\n" +
    "@SQ\tSN:chr2\tLN:243199373\n" +
    "r001\t99\tchr1\t7\t60\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\n" +
    "r002\t0\tchr1\t9\t30\t1S2I6M\t*\t0\t0\tAGCTAA\t*\n" +
    "r003\t0\tchr1\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\n" +
    "r004\t0\tchr1\t16\t30\t6M14N5M\t*\t0\t0\tATAGATTTATCT\t*\n"

  let sam = @bio.parse_sam(sam_content)
  println("=== SAM File ===")
  println("  Headers: \{sam.num_headers()}")
  println("  Records: \{sam.len()}")

  // Sequence names from header
  let seq_names = sam.get_sequence_names()
  println("  Sequences: \{seq_names}")
  for i = 0; i < seq_names.length(); i = i + 1 {
    let len = sam.get_sequence_length(seq_names[i])
    match len {
      Some(l) => println("    \{seq_names[i]}: \{l} bp")
      None => println("    \{seq_names[i]}: unknown length")
    }
  }

  // Record details
  println("\n=== Alignment Records ===")
  for i = 0; i < sam.len(); i = i + 1 {
    let rec = sam.records[i]
    println("  [\{i}] \{rec.qname}")
    println("       flag=\{rec.flag} rname=\{rec.rname} pos=\{rec.pos}")
    println("       mapq=\{rec.mapq} cigar_len=\{rec.cigar.length()}")
    println("       seq=\{rec.seq}")
    println("       align_len=\{rec.get_alignment_length()}")
    println("       end=\{rec.get_end()}")
    println("       is_paired=\{rec.is_paired()}")
    println("       is_reverse=\{rec.is_reverse()}")
    println("       is_unmapped=\{rec.is_unmapped()}")
  }

  // === Parse VCF ===
  let vcf_content = "##fileformat=VCFv4.2\n" +
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples\">\n" +
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
    "##contig=<ID=chr1,length=249250621>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\n" +
    "chr1\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3\tGT\t0|0\t1|0\n" +
    "chr1\t17330\t.\tT\tA\t3\tq10\tNS=3\tGT\t0|0\t1/1\n" +
    "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2\tGT\t1|2\t2/1\n" +
    "chr1\t1230237\t.\tT\t.\t47\tPASS\tNS=3\tGT\t0|0\t0/1\n"

  let vcf = @bio.parse_vcf(vcf_content)
  println("\n=== VCF File ===")
  let version = vcf.get_version()
  match version {
    Some(v) => println("  Version: \{v}")
    None => println("  Version: unknown")
  }
  println("  Meta headers: \{vcf.num_meta_headers()}")
  println("  Records: \{vcf.len()}")

  // Sample names
  let samples = vcf.get_sample_names()
  println("  Samples: \{samples}")

  // Variant details
  println("\n=== Variant Records ===")
  for i = 0; i < vcf.len(); i = i + 1 {
    let rec = vcf.records[i]
    println("  [\{i}] \{rec.chrom}:\{rec.pos}")
    println("       ID: \{rec.id}")
    println("       REF: \{rec.ref}  ALT: \{rec.alt}")
    match rec.qual {
      Some(q) => println("       QUAL: \{q}")
      None => println("       QUAL: .")
    }
    println("       FILTER: \{rec.filter}")
    println("       is_snp: \{rec.is_snp()}")
    println("       is_insertion: \{rec.is_insertion()}")
    println("       is_deletion: \{rec.is_deletion()}")

    // Genotypes
    for j = 0; j < samples.length(); j = j + 1 {
      let gt = rec.get_genotype(samples[j])
      match gt {
        Some(g) => println("       \{samples[j]}: \{g}")
        None => println("       \{samples[j]}: .")
      }
    }
  }
}
