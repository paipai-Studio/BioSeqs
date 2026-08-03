///|
/// Bio.SeqIO.TwoBitIO module tests - UCSC 2bit binary genome format.
// Tests for the TwoBitIO parser/writer mirroring BioPython's Bio.SeqIO.TwoBitIO.
// Run `moon test` to execute.

// ============================================================================
// Constants (verified indirectly via the serialized header, since the
// `twobit_magic` / `twobit_version` lets are package-private in `src`).
// ============================================================================

///|
test "twobit_magic_constant" {
  // The first 4 bytes of a serialized 2bit file are the little-endian magic
  // number 0x1A412743 (0x43, 0x27, 0x41, 0x1A).
  let bytes = @src.twobit_to_bytes(@src.TwoBitFile::new())
  assert_eq(bytes[0], 0x43)
  assert_eq(bytes[1], 0x27)
  assert_eq(bytes[2], 0x41)
  assert_eq(bytes[3], 0x1A)
  let magic = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24)
  assert_eq(magic, 0x1A412743)
}

///|
test "twobit_version_constant" {
  // Bytes 4-7 of the header encode the version (0).
  let bytes = @src.twobit_to_bytes(@src.TwoBitFile::new())
  assert_eq(bytes[4], 0)
  assert_eq(bytes[5], 0)
  assert_eq(bytes[6], 0)
  assert_eq(bytes[7], 0)
}

// ============================================================================
// Base encoding / decoding
// ============================================================================

///|
test "twobit_base_to_code_basic" {
  assert_eq(@src.twobit_base_to_code('A'), 0)
  assert_eq(@src.twobit_base_to_code('C'), 1)
  assert_eq(@src.twobit_base_to_code('T'), 2)
  assert_eq(@src.twobit_base_to_code('G'), 3)
  assert_eq(@src.twobit_base_to_code('N'), -1)
}

///|
test "twobit_base_to_code_invalid" {
  // Any character that is not A/C/T/G maps to -1.
  assert_eq(@src.twobit_base_to_code('X'), -1)
  assert_eq(@src.twobit_base_to_code('Z'), -1)
  assert_eq(@src.twobit_base_to_code('5'), -1)
  // twobit_base_to_code only handles uppercase A/C/T/G; lowercase is not mapped.
  assert_eq(@src.twobit_base_to_code('a'), -1)
}

///|
test "twobit_code_to_base_basic" {
  assert_eq(@src.twobit_code_to_base(0), 'A')
  assert_eq(@src.twobit_code_to_base(1), 'C')
  assert_eq(@src.twobit_code_to_base(2), 'T')
  assert_eq(@src.twobit_code_to_base(3), 'G')
}

// ============================================================================
// TwoBitNBlock
// ============================================================================

///|
test "twobit_nblock_creation" {
  let b = @src.TwoBitNBlock::new(5, 10)
  assert_eq(b.start(), 5)
  assert_eq(b.size(), 10)
  assert_eq(b.end(), 15)
}

///|
test "twobit_nblock_zero_size" {
  let b = @src.TwoBitNBlock::new(0, 0)
  assert_eq(b.start(), 0)
  assert_eq(b.size(), 0)
  assert_eq(b.end(), 0)
}

///|
test "twobit_nblock_to_string" {
  let b = @src.TwoBitNBlock::new(2, 4)
  let s = b.to_string()
  assert_true(s.contains("NBlock"))
  assert_true(s.contains("start=2"))
  assert_true(s.contains("size=4"))
}

// ============================================================================
// TwoBitMaskBlock
// ============================================================================

///|
test "twobit_maskblock_creation" {
  let b = @src.TwoBitMaskBlock::new(3, 7)
  assert_eq(b.start(), 3)
  assert_eq(b.size(), 7)
  assert_eq(b.end(), 10)
}

///|
test "twobit_maskblock_to_string" {
  let b = @src.TwoBitMaskBlock::new(1, 5)
  let s = b.to_string()
  assert_true(s.contains("MaskBlock"))
  assert_true(s.contains("start=1"))
  assert_true(s.contains("size=5"))
}

// ============================================================================
// TwoBitSequence
// ============================================================================

///|
test "twobit_sequence_creation" {
  let seq = @src.TwoBitSequence::new("myseq", 100)
  assert_eq(seq.name(), "myseq")
  assert_eq(seq.dna_size(), 100)
  assert_eq(seq.packed_dna().length(), 0)
  assert_eq(seq.n_block_count(), 0)
  assert_eq(seq.mask_block_count(), 0)
  assert_eq(seq.n_blocks().length(), 0)
  assert_eq(seq.mask_blocks().length(), 0)
}

///|
test "twobit_sequence_add_n_blocks" {
  let seq = @src.TwoBitSequence::new("s", 50)
  assert_eq(seq.n_block_count(), 0)
  seq.add_n_block(5, 10)
  assert_eq(seq.n_block_count(), 1)
  seq.add_n_block(20, 5)
  assert_eq(seq.n_block_count(), 2)
  let nbs = seq.n_blocks()
  assert_eq(nbs.length(), 2)
  assert_eq(nbs[0].start(), 5)
  assert_eq(nbs[0].size(), 10)
  assert_eq(nbs[0].end(), 15)
  assert_eq(nbs[1].start(), 20)
  assert_eq(nbs[1].size(), 5)
}

///|
test "twobit_sequence_add_mask_blocks" {
  let seq = @src.TwoBitSequence::new("s", 50)
  assert_eq(seq.mask_block_count(), 0)
  seq.add_mask_block(0, 3)
  assert_eq(seq.mask_block_count(), 1)
  seq.add_mask_block(10, 4)
  assert_eq(seq.mask_block_count(), 2)
  let mbs = seq.mask_blocks()
  assert_eq(mbs.length(), 2)
  assert_eq(mbs[0].start(), 0)
  assert_eq(mbs[0].size(), 3)
  assert_eq(mbs[1].start(), 10)
  assert_eq(mbs[1].end(), 14)
}

// ============================================================================
// TwoBitFile
// ============================================================================

///|
test "twobit_file_creation" {
  let tbf = @src.TwoBitFile::new()
  assert_eq(tbf.n_sequences(), 0)
  assert_eq(tbf.version(), 0)
  assert_eq(tbf.sequences().length(), 0)
  assert_eq(tbf.sequence_names().length(), 0)
}

///|
test "twobit_file_add_and_lookup" {
  let tbf = @src.TwoBitFile::new()
  let s1 = @src.TwoBitSequence::new("seq1", 4)
  let s2 = @src.TwoBitSequence::new("seq2", 8)
  tbf.add_sequence(s1)
  tbf.add_sequence(s2)
  assert_eq(tbf.n_sequences(), 2)
  assert_eq(tbf.version(), 0)
  let names = tbf.sequence_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "seq1")
  assert_eq(names[1], "seq2")
  let seqs = tbf.sequences()
  assert_eq(seqs.length(), 2)
  assert_eq(seqs[0].name(), "seq1")
  assert_eq(seqs[1].name(), "seq2")
  let found = tbf.get_sequence("seq2")
  assert_true(found is Some(_))
  assert_eq(found.unwrap().name(), "seq2")
  assert_eq(found.unwrap().dna_size(), 8)
  let missing = tbf.get_sequence("nope")
  assert_true(missing is None)
}

// ============================================================================
// Pack / unpack
// ============================================================================

///|
test "twobit_pack_sequence_simple" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("ATGCATGC")
  assert_eq(packed.length(), 2)
  // A=0, T=2, G=3, C=1 packed 4-per-byte: 00 10 11 01 = 0x2D = 45.
  assert_eq(packed[0], 45)
  assert_eq(packed[1], 45)
  assert_eq(n_blocks.length(), 0)
  assert_eq(mask_blocks.length(), 0)
}

///|
test "twobit_unpack_sequence_simple" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("ATGCATGC")
  let result = @src.twobit_unpack_sequence(packed, 8, n_blocks, mask_blocks)
  assert_eq(result, "ATGCATGC")
}

///|
test "twobit_pack_sequence_with_n_blocks" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("ATGCNNNNATGC")
  assert_eq(packed.length(), 3) // ceil(12/4) = 3 bytes
  assert_eq(n_blocks.length(), 1)
  assert_eq(n_blocks[0].start(), 4)
  assert_eq(n_blocks[0].size(), 4)
  assert_eq(n_blocks[0].end(), 8)
  assert_eq(mask_blocks.length(), 0)
}

///|
test "twobit_unpack_sequence_with_n_blocks" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("ATGCNNNNATGC")
  let result = @src.twobit_unpack_sequence(packed, 12, n_blocks, mask_blocks)
  assert_eq(result, "ATGCNNNNATGC")
}

///|
test "twobit_pack_sequence_with_mask" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("atgcATGC")
  assert_eq(n_blocks.length(), 0)
  assert_eq(mask_blocks.length(), 1)
  assert_eq(mask_blocks[0].start(), 0)
  assert_eq(mask_blocks[0].size(), 4)
  assert_eq(mask_blocks[0].end(), 4)
  // Lowercase bases pack as their uppercase equivalents.
  let (upper_packed, _, _) = @src.twobit_pack_sequence("ATGCATGC")
  assert_eq(packed.length(), upper_packed.length())
  for i in 0..<packed.length() {
    assert_eq(packed[i], upper_packed[i])
  }
}

///|
test "twobit_unpack_sequence_with_mask" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("atgcATGC")
  let result = @src.twobit_unpack_sequence(packed, 8, n_blocks, mask_blocks)
  assert_eq(result, "atgcATGC")
}

///|
test "twobit_pack_sequence_with_n_and_mask" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence(
    "atgcNNNNATGC",
  )
  assert_eq(packed.length(), 3)
  assert_eq(mask_blocks.length(), 1)
  assert_eq(mask_blocks[0].start(), 0)
  assert_eq(mask_blocks[0].size(), 4)
  assert_eq(n_blocks.length(), 1)
  assert_eq(n_blocks[0].start(), 4)
  assert_eq(n_blocks[0].size(), 4)
}

///|
test "twobit_unpack_sequence_with_n_and_mask" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence(
    "atgcNNNNATGC",
  )
  let result = @src.twobit_unpack_sequence(packed, 12, n_blocks, mask_blocks)
  assert_eq(result, "atgcNNNNATGC")
}

///|
test "twobit_pack_unpack_odd_length" {
  // 5 bases (not divisible by 4).
  let (p1, n1, m1) = @src.twobit_pack_sequence("ATGCA")
  assert_eq(p1.length(), 2) // ceil(5/4) = 2
  assert_eq(n1.length(), 0)
  assert_eq(m1.length(), 0)
  assert_eq(@src.twobit_unpack_sequence(p1, 5, n1, m1), "ATGCA")
  // 3 bases.
  let (p2, n2, m2) = @src.twobit_pack_sequence("ATG")
  assert_eq(p2.length(), 1)
  assert_eq(@src.twobit_unpack_sequence(p2, 3, n2, m2), "ATG")
  // 6 bases.
  let (p3, n3, m3) = @src.twobit_pack_sequence("ATGCAT")
  assert_eq(p3.length(), 2)
  assert_eq(@src.twobit_unpack_sequence(p3, 6, n3, m3), "ATGCAT")
}

// ============================================================================
// Serialization (byte array)
// ============================================================================

///|
test "twobit_serialization_roundtrip" {
  let tbf = @src.twobit_sample_file()
  let bytes = @src.twobit_to_bytes(tbf)
  assert_true(bytes.length() > 16)
  let restored = @src.twobit_from_bytes(bytes)
  assert_eq(restored.n_sequences(), 2)
  assert_eq(restored.version(), 0)
  let names = restored.sequence_names()
  assert_eq(names[0], "seq1")
  assert_eq(names[1], "seq2")
  let s1 = @src.twobit_get_sequence_string(restored, "seq1")
  assert_true(s1 is Some(_))
  assert_eq(s1.unwrap(), "ATGCATGCATGC")
  let s2 = @src.twobit_get_sequence_string(restored, "seq2")
  assert_true(s2 is Some(_))
  assert_eq(s2.unwrap(), "ATGCNNNNATGC")
  // N-block metadata survives the round trip.
  let seq2 = restored.get_sequence("seq2").unwrap()
  assert_eq(seq2.n_block_count(), 1)
  assert_eq(seq2.mask_block_count(), 0)
  assert_eq(seq2.dna_size(), 12)
}

///|
test "twobit_serialization_custom_seq" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ATGCNNNNATGC"), id="x1")
  let tbf = @src.twobit_from_seq_records([rec])
  let bytes = @src.twobit_to_bytes(tbf)
  let restored = @src.twobit_from_bytes(bytes)
  assert_eq(restored.n_sequences(), 1)
  let s = @src.twobit_get_sequence_string(restored, "x1")
  assert_true(s is Some(_))
  assert_eq(s.unwrap(), "ATGCNNNNATGC")
  let seq = restored.get_sequence("x1").unwrap()
  assert_eq(seq.n_block_count(), 1)
  assert_eq(seq.dna_size(), 12)
}

///|
test "twobit_from_bytes_invalid_magic" {
  // A header whose magic does not match yields an empty file.
  let bad = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  let tbf = @src.twobit_from_bytes(bad)
  assert_eq(tbf.n_sequences(), 0)
}

// ============================================================================
// Hex string I/O
// ============================================================================

///|
test "twobit_write_produces_hex" {
  let tbf = @src.twobit_sample_file()
  let hex = @src.twobit_write(tbf)
  // Each byte produces two hex characters, so the length is even.
  assert_eq(hex.length() % 2, 0)
  assert_true(hex.length() > 0)
  // The first serialized byte is the low byte of the magic (0x43 -> "43").
  assert_eq(hex[0].to_int(), '4'.to_int())
  assert_eq(hex[1].to_int(), '3'.to_int())
  // Every character must be a lowercase hex digit.
  let code_0 = '0'.to_int()
  let code_9 = '9'.to_int()
  let code_a = 'a'.to_int()
  let code_f = 'f'.to_int()
  for i in 0..<hex.length() {
    let c = hex[i].to_int()
    let valid = (c >= code_0 && c <= code_9) || (c >= code_a && c <= code_f)
    assert_true(valid)
  }
}

///|
test "twobit_hex_io_roundtrip" {
  let tbf = @src.twobit_sample_file()
  let hex = @src.twobit_write(tbf)
  let restored = @src.twobit_read(hex)
  assert_eq(restored.n_sequences(), 2)
  assert_eq(restored.version(), 0)
  let s1 = @src.twobit_get_sequence_string(restored, "seq1")
  assert_true(s1 is Some(_))
  assert_eq(s1.unwrap(), "ATGCATGCATGC")
  let s2 = @src.twobit_get_sequence_string(restored, "seq2")
  assert_true(s2 is Some(_))
  assert_eq(s2.unwrap(), "ATGCNNNNATGC")
}

// ============================================================================
// Get sequence string
// ============================================================================

///|
test "twobit_get_sequence_string" {
  let tbf = @src.twobit_sample_file()
  let s1 = @src.twobit_get_sequence_string(tbf, "seq1")
  assert_true(s1 is Some(_))
  assert_eq(s1.unwrap(), "ATGCATGCATGC")
  let s2 = @src.twobit_get_sequence_string(tbf, "seq2")
  assert_true(s2 is Some(_))
  assert_eq(s2.unwrap(), "ATGCNNNNATGC")
  let missing = @src.twobit_get_sequence_string(tbf, "missing")
  assert_true(missing is None)
}

// ============================================================================
// SeqRecord conversion
// ============================================================================

///|
test "twobit_to_seq_records" {
  let tbf = @src.twobit_sample_file()
  let records = @src.twobit_to_seq_records(tbf)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].name, "seq1")
  assert_eq(records[0].description, "")
  assert_eq(records[0].seq.to_string(), "ATGCATGCATGC")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "ATGCNNNNATGC")
}

///|
test "twobit_from_seq_records" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("ATGCATGC"), id="r1"),
    @src.SeqRecord::new(@src.Seq::new("ATGCNNNNATGC"), id="r2"),
  ]
  let tbf = @src.twobit_from_seq_records(records)
  assert_eq(tbf.n_sequences(), 2)
  let s1 = @src.twobit_get_sequence_string(tbf, "r1")
  assert_true(s1 is Some(_))
  assert_eq(s1.unwrap(), "ATGCATGC")
  let s2 = @src.twobit_get_sequence_string(tbf, "r2")
  assert_true(s2 is Some(_))
  assert_eq(s2.unwrap(), "ATGCNNNNATGC")
  // N-block detection from the SeqRecord sequence.
  let seq2 = tbf.get_sequence("r2").unwrap()
  assert_eq(seq2.n_block_count(), 1)
  assert_eq(seq2.mask_block_count(), 0)
  assert_eq(seq2.dna_size(), 12)
}

///|
test "twobit_seq_records_roundtrip" {
  let tbf1 = @src.twobit_sample_file()
  let records = @src.twobit_to_seq_records(tbf1)
  let tbf2 = @src.twobit_from_seq_records(records)
  assert_eq(tbf2.n_sequences(), 2)
  let s1 = @src.twobit_get_sequence_string(tbf2, "seq1")
  assert_true(s1 is Some(_))
  assert_eq(s1.unwrap(), "ATGCATGCATGC")
  let s2 = @src.twobit_get_sequence_string(tbf2, "seq2")
  assert_true(s2 is Some(_))
  assert_eq(s2.unwrap(), "ATGCNNNNATGC")
}

// ============================================================================
// Sample file validation
// ============================================================================

///|
test "twobit_sample_file" {
  let tbf = @src.twobit_sample_file()
  assert_eq(tbf.n_sequences(), 2)
  assert_eq(tbf.version(), 0)
  let names = tbf.sequence_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "seq1")
  assert_eq(names[1], "seq2")
  // seq1: ATGCATGCATGC, no N or mask blocks.
  let seq1 = tbf.get_sequence("seq1").unwrap()
  assert_eq(seq1.name(), "seq1")
  assert_eq(seq1.dna_size(), 12)
  assert_eq(seq1.n_block_count(), 0)
  assert_eq(seq1.mask_block_count(), 0)
  // seq2: ATGCNNNNATGC, one N-block at position 4 size 4.
  let seq2 = tbf.get_sequence("seq2").unwrap()
  assert_eq(seq2.name(), "seq2")
  assert_eq(seq2.dna_size(), 12)
  assert_eq(seq2.n_block_count(), 1)
  assert_eq(seq2.mask_block_count(), 0)
  let nb = seq2.n_blocks()
  assert_eq(nb.length(), 1)
  assert_eq(nb[0].start(), 4)
  assert_eq(nb[0].size(), 4)
  assert_eq(nb[0].end(), 8)
}

// ============================================================================
// Edge cases
// ============================================================================

///|
test "twobit_empty_file" {
  let tbf = @src.TwoBitFile::new()
  assert_eq(tbf.n_sequences(), 0)
  assert_eq(tbf.sequences().length(), 0)
  assert_eq(tbf.sequence_names().length(), 0)
  let bytes = @src.twobit_to_bytes(tbf)
  assert_eq(bytes.length(), 16) // header only, no index or sequence data
  let restored = @src.twobit_from_bytes(bytes)
  assert_eq(restored.n_sequences(), 0)
  // Hex round trip of an empty file (16 bytes -> 32 hex chars).
  let hex = @src.twobit_write(tbf)
  assert_eq(hex.length(), 32)
  let restored2 = @src.twobit_read(hex)
  assert_eq(restored2.n_sequences(), 0)
}

///|
test "twobit_single_sequence_file" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ACGT"), id="only")
  let tbf = @src.twobit_from_seq_records([rec])
  assert_eq(tbf.n_sequences(), 1)
  let names = tbf.sequence_names()
  assert_eq(names.length(), 1)
  assert_eq(names[0], "only")
  let s = @src.twobit_get_sequence_string(tbf, "only")
  assert_true(s is Some(_))
  assert_eq(s.unwrap(), "ACGT")
}

///|
test "twobit_single_base_pack_unpack" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("A")
  assert_eq(packed.length(), 1)
  assert_eq(n_blocks.length(), 0)
  assert_eq(mask_blocks.length(), 0)
  let result = @src.twobit_unpack_sequence(packed, 1, n_blocks, mask_blocks)
  assert_eq(result, "A")
}

///|
test "twobit_all_n_pack_unpack" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("NNNN")
  assert_eq(packed.length(), 1)
  assert_eq(n_blocks.length(), 1)
  assert_eq(n_blocks[0].start(), 0)
  assert_eq(n_blocks[0].size(), 4)
  assert_eq(n_blocks[0].end(), 4)
  assert_eq(mask_blocks.length(), 0)
  let result = @src.twobit_unpack_sequence(packed, 4, n_blocks, mask_blocks)
  assert_eq(result, "NNNN")
}

///|
test "twobit_pack_empty_sequence" {
  let (packed, n_blocks, mask_blocks) = @src.twobit_pack_sequence("")
  assert_eq(packed.length(), 0)
  assert_eq(n_blocks.length(), 0)
  assert_eq(mask_blocks.length(), 0)
  let result = @src.twobit_unpack_sequence(packed, 0, n_blocks, mask_blocks)
  assert_eq(result, "")
}
