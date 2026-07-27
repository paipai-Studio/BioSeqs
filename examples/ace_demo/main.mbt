///|
/// ACE contig assembly format parsing example.
///
/// Run: moon run examples/ace_demo/main.mbt
fn main {
  println("=== ACE Contig Assembly Format Parser Demo ===")
  println("")

  let ace_content = "AF contig1 60\n" +
    "RD read1 60 chemistry1 dye1\n" +
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" +
    "q 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30\n" +
    "RD read2 60 chemistry2 dye2\n" +
    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\n" +
    "q 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig1 60 0 60\n" +
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" +
    "q 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 60 1 60 1 60 +\n" +
    "AF read2 1 60 1 60 1 60 -\n"

  println("1. Parsing ACE content...")
  let data = @bio.ace_parse(ace_content)
  println("   Total contigs: \{data.total_contigs}")
  println("   Total reads: \{data.total_reads}")

  println("")
  println("2. Contig details:")
  let mut ci = 0
  while ci < data.contigs.length() {
    let contig = data.contigs[ci]
    println("   Contig name: \{contig.contig_name}")
    println("   Length: \{contig.length}")
    println("   Reads in contig: \{contig.reads.length()}")
    println("   GC content: \{@bio.ace_contig_gc_content(contig)}")
    ci = ci + 1
  }

  println("")
  println("3. Getting contig by name:")
  let contig = @bio.ace_get_contig(data, "contig1")
  match contig {
    Some(c) => {
      println("   Found contig: \{c.contig_name}")
      println("   Length: \{@bio.ace_contig_length(data, "contig1")}")
    }
    None => println("   Contig not found")
  }

  println("")
  println("4. Consensus sequence:")
  let cs = @bio.ace_consensus_sequence(@bio.ace_get_contig(data, "contig1").unwrap())
  println("   First 20 chars: \{cs[0:20]}...")

  println("")
  println("5. Read coverage (first 10 positions):")
  let cov = @bio.ace_read_coverage(@bio.ace_get_contig(data, "contig1").unwrap())
  let mut pos = 1
  while pos <= 10 {
    println("   Position \{pos}: depth \{cov.get(pos).unwrap_or(0)}")
    pos = pos + 1
  }

  println("")
  println("6. Reads in contig1:")
  let reads = @bio.ace_contig_reads(data, "contig1")
  let mut ri = 0
  while ri < reads.length() {
    let read = reads[ri]
    println("   Read: \{read.read_id}, length: \{read.sequence_length}, strand: \{read.strand}")
    println("     Chemistry: \{read.chemistry}, Dye: \{read.dye}")
    println("     Quality: \{read.quality.length()} values")
    ri = ri + 1
  }

  println("")
  println("7. Serializing back to ACE format...")
  let output = @bio.ace_to_string(data)
  let data2 = @bio.ace_parse(output)
  println("   Roundtrip contigs: \{data2.contigs.length()}")
  println("   Roundtrip reads: \{data2.total_reads}")

  println("")
  println("=== Done ===")
}