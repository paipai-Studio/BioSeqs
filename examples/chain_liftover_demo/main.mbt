///|
/// Chain file parsing and liftOver demo.
///
/// Demonstrates parsing a UCSC chain file, performing liftOver position
/// conversion between genome assemblies, finding chains for positions,
/// computing chain summaries, and serializing chains back to text.
///
/// Run: moon run examples/chain_liftover_demo/main.mbt
fn main {
  println("=== Chain File Parsing & liftOver Demo ===")
  println("")

  // A sample UCSC chain file mapping chr1 (hg38) to chr2 (mm10)
  // Chain format: alignment between genome assemblies
  let chain_content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  400 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  300 0 0"

  println("1. Parsing a UCSC chain file")
  let chains = @bio.cl_parse_chain(chain_content)
  println("   Chain alignments parsed: \{chains.alignments.length()}")
  println("   Target assembly: \{chains.target_name}")
  println("   Query assembly: \{chains.query_name}")

  println("")
  println("2. Inspecting chain alignments")
  let mut i = 0
  while i < chains.alignments.length() {
    let ch = chains.alignments[i]
    println(
      "   Chain \{i}: id=\{ch.id}, score=\{ch.score}, \{ch.t_name}:\{ch.t_start}-\{ch.t_end} -> \{ch.q_name}:\{ch.q_start}-\{ch.q_end} (strand=\{ch.q_strand})",
    )
    println("     Blocks: \{ch.blocks.length()}")
    let mut b = 0
    while b < ch.blocks.length() {
      let blk = ch.blocks[b]
      println(
        "       Block \{b}: size=\{blk.size}, dt=\{blk.dt}, dq=\{blk.dq}",
      )
      b = b + 1
    }
    i = i + 1
  }

  println("")
  println("3. LiftOver position conversion")
  // Lift over chr1:100 -> should map to chr2:100 (first chain, 0 mismatches)
  let pos_result = @bio.cl_liftover_pos(chains, "chr1", 100)
  match pos_result {
    Some(qpos) =>
      println("   chr1:100 -> chr2/3:\{qpos} (successful liftOver)")
    None => println("   chr1:100 -> no chain covers this position")
  }

  // Position in the first gap (dt=10 after block 0, pos 501-510 are gap)
  let gap_result = @bio.cl_liftover_pos(chains, "chr1", 501)
  match gap_result {
    Some(_) => println("   chr1:501 -> mapped (within alignment block)")
    None => println("   chr1:501 -> gap in alignment, unmappable")
  }

  // Position in the second block (pos 600 = 500+10+90, offset into block 1)
  let block2_result = @bio.cl_liftover_pos(chains, "chr1", 600)
  match block2_result {
    Some(qpos) =>
      println("   chr1:600 -> query pos:\{qpos} (second block, offset by gap)")
    None => println("   chr1:600 -> unmappable")
  }

  println("")
  println("4. LiftOver range conversion")
  let range_result = @bio.cl_liftover_range(chains, "chr1", 100, 200)
  match range_result {
    Some((qs, qe)) =>
      println("   chr1:100-200 -> query range \{qs}-\{qe}")
    None => println("   chr1:100-200 -> range unmappable")
  }

  println("")
  println("5. Finding chain for a specific position")
  let found_chain = @bio.cl_find_chain_for_pos(chains, "chr1", 6500)
  match found_chain {
    Some(ch) =>
      println("   chr1:6500 found in chain id=\{ch.id} (\{ch.t_name}:\{ch.t_start}-\{ch.t_end})")
    None => println("   chr1:6500 -> no chain found")
  }

  println("")
  println("6. Chain summary statistics")
  i = 0
  while i < chains.alignments.length() {
    let summary = @bio.cl_get_chain_summary(chains.alignments[i])
    println("   Chain \{chains.alignments[i].id}:")
    println("     Total block size: \{summary["total_block_size"]}")
    println("     Number of blocks: \{summary["num_blocks"]}")
    println("     Total gap (target): \{summary["total_gap_target"]}")
    println("     Total gap (query): \{summary["total_gap_query"]}")
    println("     Aligned target: \{summary["aligned_target"]} bp")
    println("     Aligned query: \{summary["aligned_query"]} bp")
    i = i + 1
  }

  println("")
  println("7. Serializing chain back to text")
  let serialized = @bio.cl_chain_to_string(chains.alignments[0])
  println("   First chain serialized:")
  let lines = serialized.split("\n").to_array()
  let mut l = 0
  while l < lines.length() && l < 4 {
    println("     \{lines[l]}")
    l = l + 1
  }
  if lines.length() > 4 {
    println("     ... (\{lines.length()} total lines)")
  }

  println("")
  println("=== Demo Complete ===")
}