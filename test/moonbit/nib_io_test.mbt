///|
/// Tests for the NibIO module (Bio.SeqIO.NibIO equivalent).
///
/// These tests exercise the public API exposed in `src/nib_io.mbt`:
/// the NibSequence struct, the 2-bit base encoding helpers
/// (`nib_base_to_code` / `nib_base_to_code_ci` / `nib_code_to_base`),
/// the pack/unpack functions (`nib_encode` / `nib_decode`),
/// the high-level operations (`nib_to_string`, `nib_get_base`,
/// `nib_subsequence`, `nib_reverse_complement`, `nib_gc_content`,
/// `nib_compression_ratio`), the hex I/O (`nib_to_hex` / `nib_from_hex`),
/// and the sample data (`nib_sample_data`).

// ============================================================================
// 1. Base encoding - nib_base_to_code / nib_code_to_base
// ============================================================================

///|
test "nib_base_to_code_valid" {
  // T=0, C=1, A=2, G=3 per the nib specification.
  assert_eq(@src.nib_base_to_code('T'), 0)
  assert_eq(@src.nib_base_to_code('C'), 1)
  assert_eq(@src.nib_base_to_code('A'), 2)
  assert_eq(@src.nib_base_to_code('G'), 3)
}

///|
test "nib_base_to_code_invalid" {
  assert_eq(@src.nib_base_to_code('N'), -1)
  assert_eq(@src.nib_base_to_code('X'), -1)
  assert_eq(@src.nib_base_to_code('Z'), -1)
  assert_eq(@src.nib_base_to_code('5'), -1)
  // The strict variant does not accept lowercase.
  assert_eq(@src.nib_base_to_code('a'), -1)
  assert_eq(@src.nib_base_to_code('t'), -1)
}

///|
test "nib_base_to_code_ci_accepts_lowercase" {
  assert_eq(@src.nib_base_to_code_ci('T'), 0)
  assert_eq(@src.nib_base_to_code_ci('t'), 0)
  assert_eq(@src.nib_base_to_code_ci('C'), 1)
  assert_eq(@src.nib_base_to_code_ci('c'), 1)
  assert_eq(@src.nib_base_to_code_ci('A'), 2)
  assert_eq(@src.nib_base_to_code_ci('a'), 2)
  assert_eq(@src.nib_base_to_code_ci('G'), 3)
  assert_eq(@src.nib_base_to_code_ci('g'), 3)
}

///|
test "nib_base_to_code_ci_rejects_non_acgt" {
  assert_eq(@src.nib_base_to_code_ci('N'), -1)
  assert_eq(@src.nib_base_to_code_ci('n'), -1)
  assert_eq(@src.nib_base_to_code_ci('X'), -1)
  assert_eq(@src.nib_base_to_code_ci('5'), -1)
}

///|
test "nib_code_to_base_valid" {
  assert_eq(@src.nib_code_to_base(0), 'T')
  assert_eq(@src.nib_code_to_base(1), 'C')
  assert_eq(@src.nib_code_to_base(2), 'A')
  assert_eq(@src.nib_code_to_base(3), 'G')
}

///|
test "nib_code_to_base_invalid_returns_N" {
  assert_eq(@src.nib_code_to_base(4), 'N')
  assert_eq(@src.nib_code_to_base(-1), 'N')
  assert_eq(@src.nib_code_to_base(99), 'N')
}

///|
test "nib_base_code_round_trip" {
  // code -> base -> code should be identity for all four valid codes.
  for c in 0..<4 {
    let base = @src.nib_code_to_base(c)
    assert_eq(@src.nib_base_to_code(base), c)
  }
}

// ============================================================================
// 2. nib_encode - packing DNA strings
// ============================================================================

///|
test "nib_encode_known_values" {
  // ATGC: A=2, T=0, G=3, C=1 -> (2<<6)|(0<<4)|(3<<2)|(1<<0) = 141 = 0x8D
  assert_eq(@src.nib_encode("ATGC"), [141])
  // TTTT: all 0 -> 0x00
  assert_eq(@src.nib_encode("TTTT"), [0])
  // GGGG: all 3 -> 0xFF = 255
  assert_eq(@src.nib_encode("GGGG"), [255])
  // AAAA: all 2 -> 0xAA = 170
  assert_eq(@src.nib_encode("AAAA"), [170])
  // CCCC: all 1 -> 0x55 = 85
  assert_eq(@src.nib_encode("CCCC"), [85])
}

///|
test "nib_encode_empty" {
  assert_eq(@src.nib_encode("").length(), 0)
}

///|
test "nib_encode_single_base" {
  // One base occupies the high 2 bits of the single byte.
  // T -> 0 << 6 = 0; A -> 2 << 6 = 128; C -> 1 << 6 = 64; G -> 3 << 6 = 192.
  assert_eq(@src.nib_encode("T"), [0])
  assert_eq(@src.nib_encode("C"), [64])
  assert_eq(@src.nib_encode("A"), [128])
  assert_eq(@src.nib_encode("G"), [192])
}

///|
test "nib_encode_odd_length" {
  // 5 bases -> 2 bytes; the last byte holds a single base plus zero padding.
  let packed = @src.nib_encode("ATGCA")
  assert_eq(packed.length(), 2)
  // byte 0 = ATGC = 141; byte 1 = A<<6 = 128.
  assert_eq(packed, [141, 128])
}

///|
test "nib_encode_lowercase_same_as_uppercase" {
  // Lowercase input packs identically to uppercase input.
  assert_eq(@src.nib_encode("atgc"), @src.nib_encode("ATGC"))
  assert_eq(@src.nib_encode("aaaa"), @src.nib_encode("AAAA"))
}

// ============================================================================
// 3. nib_decode - unpacking to DNA strings
// ============================================================================

///|
test "nib_decode_known_values" {
  assert_eq(@src.nib_decode([141], 4), "ATGC")
  assert_eq(@src.nib_decode([0], 4), "TTTT")
  assert_eq(@src.nib_decode([255], 4), "GGGG")
  assert_eq(@src.nib_decode([170], 4), "AAAA")
  assert_eq(@src.nib_decode([85], 4), "CCCC")
}

///|
test "nib_decode_empty" {
  assert_eq(@src.nib_decode([], 0), "")
}

///|
test "nib_decode_odd_length" {
  // 5 bases from [141, 128] -> "ATGCA".
  assert_eq(@src.nib_decode([141, 128], 5), "ATGCA")
}

// ============================================================================
// 4. Encode / decode round-trip
// ============================================================================

///|
test "nib_encode_decode_round_trip_simple" {
  let seq = "ATGC"
  let packed = @src.nib_encode(seq)
  let decoded = @src.nib_decode(packed, seq.length())
  assert_eq(decoded, seq)
}

///|
test "nib_encode_decode_round_trip_all_bases" {
  for base in ["A", "C", "G", "T"] {
    let packed = @src.nib_encode(base)
    let decoded = @src.nib_decode(packed, 1)
    assert_eq(decoded, base)
  }
}

///|
test "nib_encode_decode_round_trip_single_base" {
  let seq = "G"
  let packed = @src.nib_encode(seq)
  assert_eq(packed.length(), 1)
  assert_eq(@src.nib_decode(packed, seq.length()), seq)
}

///|
test "nib_encode_decode_round_trip_odd" {
  let seq = "ATGCA"
  let packed = @src.nib_encode(seq)
  assert_eq(packed.length(), 2)
  assert_eq(@src.nib_decode(packed, seq.length()), seq)
}

///|
test "nib_encode_decode_round_trip_long" {
  let seq = "ATGCATGCATGCATGCATGC"
  let packed = @src.nib_encode(seq)
  assert_eq(packed.length(), 5)
  assert_eq(@src.nib_decode(packed, seq.length()), seq)
}

///|
test "nib_encode_decode_round_trip_mixed" {
  let seq = "GGGCCCAAATTTGGGCCC"
  let packed = @src.nib_encode(seq)
  assert_eq(@src.nib_decode(packed, seq.length()), seq)
}

// ============================================================================
// 5. NibSequence - construction and accessors
// ============================================================================

///|
test "nib_sequence_new_uppercase" {
  let nib = @src.NibSequence::new("seq1", "ATGCATGC")
  assert_eq(nib.seq_id(), "seq1")
  assert_eq(nib.length(), 8)
  assert_true(nib.is_uppercase())
  // 8 bases -> 2 packed bytes.
  assert_eq(nib.packed_data().length(), 2)
}

///|
test "nib_sequence_new_lowercase" {
  let nib = @src.NibSequence::new("seq2", "atgcatgc")
  assert_eq(nib.length(), 8)
  // Lowercase input sets is_uppercase to false.
  assert_false(nib.is_uppercase())
  // Lowercase packs identically to uppercase.
  let upper = @src.NibSequence::new("u", "ATGCATGC")
  assert_eq(nib.packed_data(), upper.packed_data())
}

///|
test "nib_sequence_new_mixed_case_treated_as_lowercase" {
  // Any lowercase character flips is_uppercase to false.
  let nib = @src.NibSequence::new("m", "ATGCa")
  assert_false(nib.is_uppercase())
}

///|
test "nib_sequence_new_empty" {
  let nib = @src.NibSequence::new("empty", "")
  assert_eq(nib.length(), 0)
  assert_eq(nib.packed_data().length(), 0)
  assert_true(nib.is_uppercase())
}

///|
test "nib_sequence_packed_data_returns_copy" {
  let nib = @src.NibSequence::new("s", "ATGC")
  let arr = nib.packed_data()
  arr.push(99)
  // Mutating the returned copy must not affect the original.
  assert_eq(nib.packed_data().length(), 1)
}

///|
test "nib_sequence_to_string" {
  let nib = @src.NibSequence::new("seq1", "ATGC")
  let s = nib.to_string()
  assert_true(s.contains("NibSequence(seq_id=seq1"))
  assert_true(s.contains("length=4"))
  assert_true(s.contains("is_uppercase=true"))
}

// ============================================================================
// 6. nib_to_string - decode with case handling
// ============================================================================

///|
test "nib_to_string_uppercase_round_trip" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  assert_eq(@src.nib_to_string(nib), "ATGCATGC")
}

///|
test "nib_to_string_lowercase_round_trip" {
  let nib = @src.NibSequence::new("s", "atgcatgc")
  // Lowercase input is decoded back to lowercase.
  assert_eq(@src.nib_to_string(nib), "atgcatgc")
}

///|
test "nib_to_string_empty" {
  let nib = @src.NibSequence::new("e", "")
  assert_eq(@src.nib_to_string(nib), "")
}

// ============================================================================
// 7. nib_size / nib_compressed_size / nib_compression_ratio
// ============================================================================

///|
test "nib_size" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  assert_eq(@src.nib_size(nib), 8)
}

///|
test "nib_compressed_size" {
  // 8 bases -> 2 bytes; 9 bases -> 3 bytes; 1 base -> 1 byte; 0 bases -> 0.
  assert_eq(@src.nib_compressed_size(@src.NibSequence::new("a", "ATGCATGC")), 2)
  assert_eq(@src.nib_compressed_size(@src.NibSequence::new("b", "ATGCATGCA")), 3)
  assert_eq(@src.nib_compressed_size(@src.NibSequence::new("c", "A")), 1)
  assert_eq(@src.nib_compressed_size(@src.NibSequence::new("d", "")), 0)
}

///|
test "nib_compression_ratio" {
  // 8 bases / 2 bytes = 4.0; 4 bases / 1 byte = 4.0; 1 base / 1 byte = 1.0.
  let r1 = @src.nib_compression_ratio(@src.NibSequence::new("a", "ATGCATGC"))
  assert_true(r1 > 3.99 && r1 < 4.01)
  let r2 = @src.nib_compression_ratio(@src.NibSequence::new("b", "ATGC"))
  assert_true(r2 > 3.99 && r2 < 4.01)
  let r3 = @src.nib_compression_ratio(@src.NibSequence::new("c", "A"))
  assert_true(r3 > 0.99 && r3 < 1.01)
  // Empty sequence returns 0.0 to avoid division by zero.
  let r4 = @src.nib_compression_ratio(@src.NibSequence::new("d", ""))
  assert_true(r4 == 0.0)
}

///|
test "nib_compression_ratio_odd_length" {
  // 5 bases / 2 bytes = 2.5.
  let r = @src.nib_compression_ratio(@src.NibSequence::new("a", "ATGCA"))
  assert_true(r > 2.49 && r < 2.51)
}

// ============================================================================
// 8. nib_get_base
// ============================================================================

///|
test "nib_get_base_valid" {
  let nib = @src.NibSequence::new("s", "ATGC")
  assert_eq(@src.nib_get_base(nib, 0), "A")
  assert_eq(@src.nib_get_base(nib, 1), "T")
  assert_eq(@src.nib_get_base(nib, 2), "G")
  assert_eq(@src.nib_get_base(nib, 3), "C")
}

///|
test "nib_get_base_out_of_bounds" {
  let nib = @src.NibSequence::new("s", "ATGC")
  assert_eq(@src.nib_get_base(nib, -1), "")
  assert_eq(@src.nib_get_base(nib, 4), "")
  assert_eq(@src.nib_get_base(nib, 100), "")
}

///|
test "nib_get_base_lowercase" {
  let nib = @src.NibSequence::new("s", "atgc")
  assert_eq(@src.nib_get_base(nib, 0), "a")
  assert_eq(@src.nib_get_base(nib, 3), "c")
}

// ============================================================================
// 9. nib_subsequence
// ============================================================================

///|
test "nib_subsequence_basic" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  assert_eq(@src.nib_subsequence(nib, 0, 4), "ATGC")
  assert_eq(@src.nib_subsequence(nib, 4, 4), "ATGC")
  assert_eq(@src.nib_subsequence(nib, 2, 4), "GCAT")
}

///|
test "nib_subsequence_clamped" {
  let nib = @src.NibSequence::new("s", "ATGC")
  // Request more than available; should clamp to the remaining bases.
  assert_eq(@src.nib_subsequence(nib, 2, 10), "GC")
}

///|
test "nib_subsequence_empty_cases" {
  let nib = @src.NibSequence::new("s", "ATGC")
  assert_eq(@src.nib_subsequence(nib, -1, 2), "")
  assert_eq(@src.nib_subsequence(nib, 4, 2), "")
  assert_eq(@src.nib_subsequence(nib, 0, 0), "")
  assert_eq(@src.nib_subsequence(nib, 0, -1), "")
}

///|
test "nib_subsequence_lowercase" {
  let nib = @src.NibSequence::new("s", "atgcatgc")
  assert_eq(@src.nib_subsequence(nib, 0, 4), "atgc")
}

///|
test "nib_subsequence_single_base" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  assert_eq(@src.nib_subsequence(nib, 3, 1), "C")
}

// ============================================================================
// 10. nib_reverse_complement
// ============================================================================

///|
test "nib_reverse_complement_basic" {
  let nib = @src.NibSequence::new("s", "ATGC")
  // Reverse of "ATGC" is "CGTA"; complement of each is "GCAT".
  assert_eq(@src.nib_to_string(@src.nib_reverse_complement(nib)), "GCAT")
}

///|
test "nib_reverse_complement_palindrome" {
  // "ATAT" reverse-complements to itself.
  let nib = @src.NibSequence::new("s", "ATAT")
  assert_eq(@src.nib_to_string(@src.nib_reverse_complement(nib)), "ATAT")
}

///|
test "nib_reverse_complement_double_round_trip" {
  // Reversing the complement twice yields the original sequence.
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  let rc1 = @src.nib_reverse_complement(nib)
  let rc2 = @src.nib_reverse_complement(rc1)
  assert_eq(@src.nib_to_string(rc2), "ATGCATGC")
}

///|
test "nib_reverse_complement_preserves_case" {
  let nib = @src.NibSequence::new("s", "atgc")
  let rc = @src.nib_reverse_complement(nib)
  assert_false(rc.is_uppercase())
  assert_eq(@src.nib_to_string(rc), "gcat")
}

///|
test "nib_reverse_complement_uses_rc_suffix" {
  let nib = @src.NibSequence::new("myseq", "ATGC")
  let rc = @src.nib_reverse_complement(nib)
  assert_eq(rc.seq_id(), "myseq_rc")
}

///|
test "nib_reverse_complement_empty" {
  let nib = @src.NibSequence::new("s", "")
  let rc = @src.nib_reverse_complement(nib)
  assert_eq(@src.nib_to_string(rc), "")
  assert_eq(@src.nib_size(rc), 0)
}

// ============================================================================
// 11. nib_gc_content
// ============================================================================

///|
test "nib_gc_content_all_gc" {
  let nib = @src.NibSequence::new("s", "GCGCGCGC")
  let gc = @src.nib_gc_content(nib)
  assert_true(gc > 0.99 && gc < 1.01)
}

///|
test "nib_gc_content_no_gc" {
  let nib = @src.NibSequence::new("s", "ATATATAT")
  let gc = @src.nib_gc_content(nib)
  assert_true(gc == 0.0)
}

///|
test "nib_gc_content_half" {
  let nib = @src.NibSequence::new("s", "ATGC")
  // 2 of 4 bases are G/C -> 0.5.
  let gc = @src.nib_gc_content(nib)
  assert_true(gc > 0.49 && gc < 0.51)
}

///|
test "nib_gc_content_empty" {
  let nib = @src.NibSequence::new("s", "")
  assert_true(@src.nib_gc_content(nib) == 0.0)
}

///|
test "nib_gc_content_quarter" {
  let nib = @src.NibSequence::new("s", "ATTTG")
  // 1 G out of 5 bases -> 0.2.
  let gc = @src.nib_gc_content(nib)
  assert_true(gc > 0.19 && gc < 0.21)
}

// ============================================================================
// 12. Hex string I/O - nib_to_hex / nib_from_hex
// ============================================================================

///|
test "nib_to_hex_known_value" {
  // "ATGC" packs to [141] = [0x8D] -> hex "8d".
  let nib = @src.NibSequence::new("s", "ATGC")
  assert_eq(@src.nib_to_hex(nib), "8d")
}

///|
test "nib_to_hex_empty" {
  let nib = @src.NibSequence::new("s", "")
  assert_eq(@src.nib_to_hex(nib), "")
}

///|
test "nib_to_hex_lowercase_output" {
  let nib = @src.NibSequence::new("s", "GGGG")
  // GGGG packs to [255] = [0xFF] -> hex "ff".
  let hex = @src.nib_to_hex(nib)
  assert_eq(hex, "ff")
  // Every character must be a lowercase hex digit.
  for c in hex {
    let code = c.to_int()
    let is_digit = code >= '0'.to_int() && code <= '9'.to_int()
    let is_lower = code >= 'a'.to_int() && code <= 'f'.to_int()
    assert_true(is_digit || is_lower)
  }
}

///|
test "nib_to_hex_long_sequence" {
  let nib = @src.NibSequence::new("s", "ATGCATGCATGC")
  // 3 packed bytes of 0x8D each -> "8d8d8d".
  assert_eq(@src.nib_to_hex(nib), "8d8d8d")
}

///|
test "nib_from_hex_round_trip" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  let hex = @src.nib_to_hex(nib)
  let nib2 = @src.nib_from_hex("s2", hex, 8)
  assert_eq(nib2.length(), 8)
  assert_eq(@src.nib_to_string(nib2), "ATGCATGC")
}

///|
test "nib_from_hex_preserves_packed_data" {
  let nib = @src.NibSequence::new("s", "GGGGCCCC")
  let hex = @src.nib_to_hex(nib)
  let nib2 = @src.nib_from_hex("s2", hex, 8)
  assert_eq(nib2.packed_data(), nib.packed_data())
}

///|
test "nib_from_hex_defaults_uppercase" {
  // Hex storage does not preserve case; the result is always uppercase.
  let nib = @src.nib_from_hex("s", "8d", 4)
  assert_true(nib.is_uppercase())
}

///|
test "nib_from_hex_odd_length" {
  // 5 bases from "8d80" (byte 1 = 0x80 = A<<6).
  let nib = @src.nib_from_hex("s", "8d80", 5)
  assert_eq(nib.length(), 5)
  assert_eq(@src.nib_to_string(nib), "ATGCA")
}

// ============================================================================
// 13. N handling
// ============================================================================

///|
test "nib_encode_n_packed_as_t" {
  // N maps to code 0 (T) during packing; this is a lossy transformation.
  let packed = @src.nib_encode("NNNN")
  assert_eq(packed, [0])
  assert_eq(@src.nib_decode(packed, 4), "TTTT")
}

///|
test "nib_encode_mixed_with_n" {
  // "AN" -> A=2, N->0; packed = (2<<6)|(0<<4) = 128.
  let packed = @src.nib_encode("AN")
  assert_eq(packed, [128])
  // Decoding gives "AT" (N is lossy -> T).
  assert_eq(@src.nib_decode(packed, 2), "AT")
}

///|
test "nib_sequence_with_n_round_trip_lossy" {
  let nib = @src.NibSequence::new("s", "ATNGC")
  // N becomes T after a round-trip.
  assert_eq(@src.nib_to_string(nib), "ATTGC")
}

// ============================================================================
// 14. Sample data
// ============================================================================

///|
test "nib_sample_data_structure" {
  let nib = @src.nib_sample_data()
  assert_eq(nib.seq_id(), "sample1")
  assert_eq(nib.length(), 12)
  assert_true(nib.is_uppercase())
}

///|
test "nib_sample_data_decodes_correctly" {
  let nib = @src.nib_sample_data()
  assert_eq(@src.nib_to_string(nib), "ATGCATGCATGC")
}

///|
test "nib_sample_data_packed" {
  let nib = @src.nib_sample_data()
  // 12 bases -> 3 packed bytes, each 0x8D = 141.
  assert_eq(nib.packed_data(), [141, 141, 141])
  assert_eq(@src.nib_compressed_size(nib), 3)
}

///|
test "nib_sample_data_hex" {
  let nib = @src.nib_sample_data()
  assert_eq(@src.nib_to_hex(nib), "8d8d8d")
}

///|
test "nib_sample_data_gc_content" {
  let nib = @src.nib_sample_data()
  // "ATGCATGCATGC": 6 G/C out of 12 -> 0.5.
  let gc = @src.nib_gc_content(nib)
  assert_true(gc > 0.49 && gc < 0.51)
}

///|
test "nib_sample_data_reverse_complement" {
  let nib = @src.nib_sample_data()
  let rc = @src.nib_reverse_complement(nib)
  // Reverse complement of ATGCATGCATGC is GCATGCATGCAT.
  assert_eq(@src.nib_to_string(rc), "GCATGCATGCAT")
}

// ============================================================================
// 15. Edge cases
// ============================================================================

///|
test "edge_case_single_base_each" {
  for base in ["A", "C", "G", "T"] {
    let nib = @src.NibSequence::new("s", base)
    assert_eq(@src.nib_size(nib), 1)
    assert_eq(@src.nib_compressed_size(nib), 1)
    assert_eq(@src.nib_to_string(nib), base)
    assert_eq(@src.nib_get_base(nib, 0), base)
  }
}

///|
test "edge_case_empty_sequence_full_round_trip" {
  let nib = @src.NibSequence::new("empty", "")
  assert_eq(@src.nib_size(nib), 0)
  assert_eq(@src.nib_compressed_size(nib), 0)
  assert_eq(@src.nib_compression_ratio(nib), 0.0)
  assert_eq(@src.nib_to_string(nib), "")
  assert_eq(@src.nib_to_hex(nib), "")
  assert_eq(@src.nib_gc_content(nib), 0.0)
  assert_eq(@src.nib_subsequence(nib, 0, 0), "")
  assert_eq(@src.nib_get_base(nib, 0), "")
}

///|
test "edge_case_two_bases" {
  let nib = @src.NibSequence::new("s", "AG")
  assert_eq(@src.nib_compressed_size(nib), 1)
  assert_eq(@src.nib_to_string(nib), "AG")
  assert_eq(@src.nib_subsequence(nib, 0, 1), "A")
  assert_eq(@src.nib_subsequence(nib, 1, 1), "G")
}

///|
test "edge_case_three_bases" {
  let nib = @src.NibSequence::new("s", "ATG")
  assert_eq(@src.nib_compressed_size(nib), 1)
  assert_eq(@src.nib_to_string(nib), "ATG")
}

///|
test "edge_case_hex_round_trip_various_lengths" {
  for seq in ["A", "AT", "ATG", "ATGC", "ATGCA", "ATGCAT", "ATGCATG", "ATGCATGC"] {
    let nib = @src.NibSequence::new("s", seq)
    let hex = @src.nib_to_hex(nib)
    let nib2 = @src.nib_from_hex("s2", hex, seq.length())
    assert_eq(@src.nib_to_string(nib2), seq)
  }
}

///|
test "edge_case_long_sequence_round_trip" {
  let sb = StringBuilder::new()
  for i in 0..<100 {
    let bases = ["A", "C", "G", "T"]
    sb.write_string(bases[i % 4])
  }
  let seq = sb.to_string()
  let nib = @src.NibSequence::new("long", seq)
  assert_eq(@src.nib_size(nib), 100)
  assert_eq(@src.nib_compressed_size(nib), 25)
  assert_eq(@src.nib_to_string(nib), seq)
  assert_eq(@src.nib_compression_ratio(nib), 4.0)
}

///|
test "edge_case_subsequence_full" {
  let nib = @src.NibSequence::new("s", "ATGCATGC")
  // Taking the whole sequence via subsequence equals the decoded string.
  assert_eq(@src.nib_subsequence(nib, 0, 8), "ATGCATGC")
}
