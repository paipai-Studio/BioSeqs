///|
/// MAF (Multiple Alignment Format) Demo.
/// Demonstrates parsing, writing, and analyzing MAF alignment data.

fn main {
  println("=== MAF (Multiple Alignment Format) Demo ===\n")

  println("1. Parsing a MAF string with sequences...")
  let maf_content = "##maf version=1 score=100\n\n" +
    "a score=80\n" +
    "s human.chr1 0 10 + 10000 ACGTACGTAC\n" +
    "s mouse.chr1 0 10 + 8000 ACGTACGTAC\n" +
    "s rat.chr1 0 10 + 6000 AC--ACGTAC\n\n" +
    "a score=50\n" +
    "s human.chr1 100 14 + 10000 TTTTAAAACCCCCC\n" +
    "s mouse.chr1 100 14 + 8000 TTTTAAAACCCCCC\n" +
    "s rat.chr1 100 14 - 6000 GGGAATTCCGGGGG\n\n" +
    "a score=60\n" +
    "s human.chr1 5000 400 + 250000000\n" +
    "s mouse.chr1 5200 380 + 240000000\n" +
    "s rat.chr1 5100 360 + 230000000\n\n"
  let ali = @src.parse_maf(maf_content)
  println("   Version: \{ali.version}")
  println("   Number of blocks: \{@src.maf_num_blocks(ali)}")
  println("   Total alignment length: \{@src.maf_total_length(ali)} columns")
  let seq_names = @src.maf_all_seq_names(ali)
  println("   Sequences: \{seq_names}")

  println("\n2. Calculating percent identity per block...")
  let mut i = 0
  while i < ali.blocks.length() {
    let block = ali.blocks[i]
    let identity = @src.maf_block_percent_identity(block)
    let n_seqs = @src.maf_block_num_seqs(block)
    let aln_len = @src.maf_block_alignment_length(block)
    println("   Block \{i}: %identity = \{identity}%, seqs = \{n_seqs}, columns = \{aln_len}")
    i = i + 1
  }

  println("\n3. Computing alignment statistics...")
  let stats = @src.maf_compute_stats(ali)
  let stats_str = @src.maf_stats_to_string(stats)
  println("\{stats_str}")

  println("4. Selecting sequences (human.chr1 and mouse.chr1 only)...")
  let selected = @src.maf_select_seqs(ali, ["human.chr1", "mouse.chr1"])
  println("   Selected blocks: \{@src.maf_num_blocks(selected)}")
  let sel_names = @src.maf_all_seq_names(selected)
  println("   Selected sequences: \{sel_names}")
  let avg_id_sel = @src.maf_average_percent_identity(selected)
  println("   Average percent identity (selected): \{avg_id_sel}%")

  println("\n5. Filtering by minimum alignment length...")
  let filtered = @src.maf_filter_by_length(ali, 10)
  println("   Original blocks: \{@src.maf_num_blocks(ali)}")
  println("   Filtered blocks (min length 10): \{@src.maf_num_blocks(filtered)}")
  let mut j = 0
  while j < filtered.blocks.length() {
    let block = filtered.blocks[j]
    let len = @src.maf_block_alignment_length(block)
    println("   Block \{j}: length = \{len}")
    j = j + 1
  }

  println("\n6. Writing MAF alignment back to string...")
  let output = @src.write_maf(ali)
  let lines = output.split("\n").to_array()
  println("   Output has \{lines.length()} lines")
  println("   First line: \{lines[0]}")
  let has_human = output.contains("human.chr1")
  let has_mouse = output.contains("mouse.chr1")
  println("   Contains human.chr1: \{has_human}")
  println("   Contains mouse.chr1: \{has_mouse}")

  println("\n7. Block-level operations...")
  let block0 = ali.blocks[0]
  let coverage0 = @src.maf_block_coverage(block0)
  println("   Block 0 coverage: \{coverage0}")
  let seq_cov = @src.maf_sequence_coverage(ali)
  let keys = seq_cov.keys().collect()
  let mut k = 0
  while k < keys.length() {
    let key = keys[k]
    println("   Sequence '\{key}' coverage: \{seq_cov[key]}%")
    k = k + 1
  }

  println("\n=== MAF Demo Complete ===")
}