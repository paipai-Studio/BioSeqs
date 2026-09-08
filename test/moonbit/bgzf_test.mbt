///|
/// Test BGZF decompression functionality.

///|
/// Test bgzf_decompress with empty data.
test "bgzf_decompress_empty" {
  let result = @src.bgzf_decompress("")
  assert_eq(result, "")
}

///|
/// Test is_bgzf function.
test "bgzf_is_bgzf" {
  let bgzf_magic = "\x1f\x8b\x08\x04"
  assert_true(@src.is_bgzf(bgzf_magic))

  let not_bgzf = "BAM\x01"
  assert_true(!@src.is_bgzf(not_bgzf))
}

///|
/// The 51-byte stored block + 28-byte EOF block produced by the Biopython
/// reference implementation for the 20-byte payload "Hello BGZF\nLine two\n".
fn bgzf_reference_block() -> String {
  // 51-byte data block (stored DEFLATE)
  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x32\x00\x01\x14\x00\xeb\xff\x48\x65\x6c\x6c\x6f\x20\x42\x47\x5a\x46\x0a\x4c\x69\x6e\x65\x20\x74\x77\x6f\x0a\x90\xef\xcb\x4d\x14\x00\x00\x00" +
  // 28-byte canonical EOF block
  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
}

///|
/// Our stored-block writer must produce byte-identical output to the
/// Biopython reference for the same payload.
test "bgzf_compress_matches_reference_bytes" {
  let compressed = try {
    @src.bgzf_compress("Hello BGZF\nLine two\n")
  } catch {
    _ => ""
  }
  assert_eq(compressed, bgzf_reference_block())
  assert_eq(compressed.length(), 79)
}

///|
/// Decompress the reference block back to the original payload.
test "bgzf_decompress_reference_block" {
  let result = try {
    @src.bgzf_decompress(bgzf_reference_block())
  } catch {
    _ => ""
  }
  assert_eq(result, "Hello BGZF\nLine two\n")
}

///|
/// Compress then decompress round-trips arbitrary text.
test "bgzf_compress_decompress_roundtrip" {
  let text =
    ">seq1 desc\nACGTACGTACGT\n>seq2\nNNNNACGT\n"
  let compressed = try {
    @src.bgzf_compress(text)
  } catch {
    _ => ""
  }
  assert_true(@src.is_bgzf(compressed))
  let restored = try {
    @src.bgzf_decompress(compressed)
  } catch {
    _ => ""
  }
  assert_eq(restored, text)
}

///|
/// Data larger than one block (max payload 65505 bytes) is split into
/// multiple blocks and round-trips losslessly.
test "bgzf_compress_large_multiblock" {
  let size = 140000
  let big = String::make(size, 'A')
  let compressed = try {
    @src.bgzf_compress(big)
  } catch {
    _ => ""
  }
  let blocks = try {
    @src.bgzf_blocks(compressed)
  } catch {
    _ => []
  }
  // 65505 + 65505 + 8990 = 140000, then the 28-byte EOF block.
  assert_eq(blocks.length(), 4)
  assert_eq(blocks[0].data_length, 65505)
  assert_eq(blocks[0].block_length, 65536)
  assert_eq(blocks[1].data_length, 65505)
  assert_eq(blocks[2].data_length, 8990)
  assert_eq(blocks[3].data_length, 0)
  let restored = try {
    @src.bgzf_decompress(compressed)
  } catch {
    _ => ""
  }
  assert_eq(restored.length(), size)
  assert_eq(restored, big)
}

///|
/// bgzf_blocks reports the metadata of the reference stream correctly.
test "bgzf_blocks_reference" {
  let blocks = try {
    @src.bgzf_blocks(bgzf_reference_block())
  } catch {
    _ => []
  }
  assert_eq(blocks.length(), 2)
  assert_eq(blocks[0].start_offset, 0)
  assert_eq(blocks[0].block_length, 51)
  assert_eq(blocks[0].data_offset, 0)
  assert_eq(blocks[0].data_length, 20)
  assert_eq(blocks[1].start_offset, 51)
  assert_eq(blocks[1].block_length, 28)
  assert_eq(blocks[1].data_offset, 20)
  assert_eq(blocks[1].data_length, 0)
}

///|
/// Virtual offsets pack/unpack as (block offset << 16 | within-block offset).
test "bgzf_virtual_offset" {
  let voff = @src.bgzf_make_virtual_offset(51, 12)
  assert_eq(voff, (51L << 16) | 12L)
  let (coff, uoff) = @src.bgzf_split_virtual_offset(voff)
  assert_eq(coff, 51)
  assert_eq(uoff, 12)
  // Offset 0 within block 0 is just the within-block offset.
  let (c0, u0) = @src.bgzf_split_virtual_offset(
    @src.bgzf_make_virtual_offset(0, 5),
  )
  assert_eq(c0, 0)
  assert_eq(u0, 5)
}

///|
/// BgzfReader read(size), readline, tell and seek over a single block.
test "bgzf_reader_read_readline_seek" {
  try {
    let reader = @src.BgzfReader::new(bgzf_reference_block())
    // read(size)
    let head = reader.read(size=6)
    assert_eq(head, "Hello ")
    let (c1, u1) = @src.bgzf_split_virtual_offset(reader.tell())
    assert_eq(c1, 0)
    assert_eq(u1, 6)
    // Seek back to the position and read the rest.
    reader.seek(@src.bgzf_make_virtual_offset(0, 6))
    let rest = reader.read()
    assert_eq(rest, "BGZF\nLine two\n")

    // Fresh reader for line reading.
    let lines = @src.BgzfReader::new(bgzf_reference_block())
    assert_eq(lines.readline(), "Hello BGZF\n")
    assert_eq(lines.readline(), "Line two\n")
    assert_eq(lines.readline(), "")
  } catch {
    @src.BgzfError(_) => assert_true(false)
    _ => assert_true(false)
  }
}

///|
/// BgzfReader reads across block boundaries and seeks into later blocks.
test "bgzf_reader_multiblock" {
  let size = 140000
  let big = String::make(size, 'A')
  try {
    let compressed = @src.bgzf_compress(big)
    let blocks = @src.bgzf_blocks(compressed)

    // Reading across the first block boundary.
    let reader = @src.BgzfReader::new(compressed)
    let part = reader.read(size=70000)
    assert_eq(part.length(), 70000)
    let (coff, uoff) = @src.bgzf_split_virtual_offset(reader.tell())
    assert_eq(coff, blocks[1].start_offset)
    assert_eq(uoff, 70000 - 65505)
    assert_eq(reader.read().length(), size - 70000)

    // Seek into the second block and read 10 bytes.
    let reader2 = @src.BgzfReader::new(compressed)
    reader2.seek(
      @src.bgzf_make_virtual_offset(blocks[1].start_offset, 100),
    )
    let chunk = reader2.read(size=10)
    assert_eq(chunk.length(), 10)
    assert_eq(chunk, "AAAAAAAAAA")

    // Seek into the final (short) block.
    let reader3 = @src.BgzfReader::new(compressed)
    reader3.seek(
      @src.bgzf_make_virtual_offset(blocks[2].start_offset, 5),
    )
    assert_eq(reader3.read().length(), 8990 - 5)
  } catch {
    @src.BgzfError(_) => assert_true(false)
    _ => assert_true(false)
  }
}

///|
/// Seeking to an unknown compressed offset raises.
test "bgzf_reader_seek_bad_offset" {
  let raised = try {
    let reader = @src.BgzfReader::new(bgzf_reference_block())
    reader.seek(@src.bgzf_make_virtual_offset(999999, 0))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
/// Writing to a closed writer raises.
test "bgzf_writer_write_after_close" {
  let raised = try {
    let writer = @src.BgzfWriter::new()
    writer.write("x")
    writer.close()
    writer.write("y")
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
/// tell() on a writer reflects buffered bytes before flushing.
test "bgzf_writer_tell" {
  try {
    let writer = @src.BgzfWriter::new()
    writer.write("abc")
    let (coff, uoff) = @src.bgzf_split_virtual_offset(writer.tell())
    assert_eq(coff, 0)
    assert_eq(uoff, 3)
    writer.close()
    // After close everything is flushed to blocks: the virtual offset is
    // the end of the compressed file (31 overhead + 3 payload = 34 byte
    // block, then the 28 byte EOF marker = 62 bytes total).
    let (c2, u2) = @src.bgzf_split_virtual_offset(writer.tell())
    assert_eq(c2, 62)
    assert_eq(u2, 0)
  } catch {
    @src.BgzfError(_) => assert_true(false)
    _ => assert_true(false)
  }
}
