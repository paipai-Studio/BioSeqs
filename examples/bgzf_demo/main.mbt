///|
/// BGZF demo (port of Biopython Bio.bgzf):
///  1. compress text with BgzfWriter (stored DEFLATE blocks + EOF marker),
///  2. scan the stream with bgzf_blocks to inspect block layout,
///  3. read it back with BgzfReader: read/readline and virtual-offset seek.
fn main {
  let text =
    ">chr1 demo chromosome\n" +
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    "TTTTAAAACCCCGGGGACGTACGTACGTACGTACGTACGTACGTACGT\n" +
    ">chr2\n" +
    "NNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"

  // ---- Compress -----------------------------------------------------------
  let compressed = try {
    @src.bgzf_compress(text)
  } catch {
    @src.BgzfError(msg) => {
      println("compression failed: " + msg)
      return
    }
  }
  println("original   : " + text.length().to_string() + " bytes")
  println("compressed : " + compressed.length().to_string() + " bytes")

  // ---- Inspect blocks -----------------------------------------------------
  let blocks = try {
    @src.bgzf_blocks(compressed)
  } catch {
    @src.BgzfError(msg) => {
      println("block scan failed: " + msg)
      return
    }
  }
  println("blocks     : " + blocks.length().to_string())
  let mut i = 0
  while i < blocks.length() {
    let b = blocks[i]
    println(
      "  block " + i.to_string() +
      ": file offset " + b.start_offset.to_string() +
      ", on-disk " + b.block_length.to_string() +
      " bytes, payload " + b.data_length.to_string() + " bytes",
    )
    i = i + 1
  }

  // ---- Random-access read -------------------------------------------------
  try {
    let reader = @src.BgzfReader::new(compressed)

    // Remember a virtual offset after the first line, then rewind to it.
    let line1 = reader.readline()
    let checkpoint = reader.tell()
    println("\nfirst line : " + line1)
    println("checkpoint virtual offset: " + checkpoint.to_string())

    let skipped = reader.read(size=20)
    println("skipped 20 bytes: " + skipped)

    reader.seek(checkpoint)
    let (coff, uoff) = @src.bgzf_split_virtual_offset(checkpoint)
    println(
      "seeked back to block " + coff.to_string() +
      " + " + uoff.to_string() + " bytes",
    )
    let line2 = reader.readline()
    println("next line  : " + line2)

    // Read everything from the start.
    reader.seek(@src.bgzf_make_virtual_offset(0, 0))
    let all = reader.read()
    println("full round-trip matches: " + (all == text).to_string())
  } catch {
    @src.BgzfError(msg) => println("reading failed: " + msg)
  }
}
