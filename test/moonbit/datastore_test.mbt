///|
test "md5_hex_rfc1321_vectors" {
  // Test vectors from RFC 1321.
  assert_eq(@src.md5_hex(""), "d41d8cd98f00b204e9800998ecf8427e")
  assert_eq(@src.md5_hex("a"), "0cc175b9c0f1b6a831c399e269772661")
  assert_eq(@src.md5_hex("abc"), "900150983cd24fb0d6963f7d28e17f72")
  assert_eq(@src.md5_hex("message digest"), "f96b697d7cb7938d525a2f31aaf161d0")
  assert_eq(
    @src.md5_hex("abcdefghijklmnopqrstuvwxyz"),
    "c3fcd3d76192e4007dfb496cca67e13b",
  )
  assert_eq(
    @src.md5_hex(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
    ),
    "d174ab98d277d9f5a5611c2c9f419d9f",
  )
  assert_eq(
    @src.md5_hex("12345678901234567890123456789012345678901234567890123456789012345678901234567890"),
    "57edf4a22be3c955ac49da2e2107b67a",
  )
}

///|
test "md5_hex_block_boundaries" {
  // Message whose padded length spans several 64-byte blocks.
  assert_eq(
    @src.md5_hex("The quick brown fox jumps over the lazy dog"),
    "9e107d9d372bb6826bd81d3542a419d6",
  )
  let mut seq = ""
  for _ in 0..<250 {
    seq = seq + "ACGT"
  }
  assert_eq(seq.length(), 1000)
  assert_eq(@src.md5_hex(seq), "dee764dc7924160b161de01b8c77167c")
  // Multi-byte UTF-8 characters are hashed over their UTF-8 bytes.
  assert_eq(
    @src.md5_hex("café résumé"),
    "10871a98898bc4c0980b3224dddbc33f",
  )
}

///|
test "datastore_add_get_roundtrip" {
  let bag = @src.Datastore::new(root="demo-bag")
  assert_eq(bag.root(), "demo-bag")
  assert_eq(bag.length(), 0)
  let sum = bag.add("reads.fasta", ">r1\nACGT\n")
  assert_eq(sum, @src.md5_hex(">r1\nACGT\n"))
  assert_eq(bag.length(), 1)
  assert_true(bag.contains("reads.fasta"))
  assert_eq(bag.get("reads.fasta"), Some(">r1\nACGT\n"))
  assert_eq(bag.checksum("reads.fasta"), Some(sum))
  assert_eq(bag.get("missing.fasta"), None)
  // Re-adding replaces the payload.
  let sum2 = bag.add("reads.fasta", ">r1\nACGTA\n")
  assert_eq(bag.length(), 1)
  assert_eq(bag.checksum("reads.fasta"), Some(sum2))
  // Removal.
  assert_true(bag.remove("reads.fasta"))
  assert_false(bag.remove("reads.fasta"))
  assert_eq(bag.length(), 0)
}

///|
test "datastore_names_sorted" {
  let bag = @src.Datastore::new(root="b")
  let _ = bag.add("gamma", "3")
  let _ = bag.add("alpha", "1")
  let _ = bag.add("beta", "2")
  assert_eq(bag.names(), ["alpha", "beta", "gamma"])
}

///|
test "datastore_manifest_and_verify" {
  let bag = @src.Datastore::new(root="b")
  let s1 = bag.add("tree.newick", "((a,b),c);")
  let s2 = bag.add("align.fasta", ">a\nAC\n>b\nAG\n")
  let lines = bag.manifest_lines()
  assert_eq(lines.length(), 2)
  // Manifest lines are sorted by payload name.
  assert_eq(lines[0], s2 + "  data/align.fasta")
  assert_eq(lines[1], s1 + "  data/tree.newick")
  assert_true(
    bag.manifest_md5_txt().contains(s2 + "  data/align.fasta\n"),
  )
  assert_true(bag.verify())
  // Entry accessor.
  let entry = bag.entry("tree.newick")
  match entry {
    Some(e) => {
      assert_eq(e.name(), "tree.newick")
      assert_eq(e.data(), "((a,b),c);")
      assert_eq(e.checksum(), s1)
    }
    None => fail("expected entry")
  }
  assert_true(bag.entry("nope") is None)
}

///|
test "datastore_load_bag_verify_detects_corruption" {
  let good = "((a,b),c);"
  let bad_checksum_manifest = "00000000000000000000000000000000  data/tree.newick\n"
  // Declared checksum disagrees with payload contents -> verify fails.
  let bag = @src.Datastore::load_bag(
    root="corrupt",
    manifest_text=bad_checksum_manifest,
    loader=fn(_) { Some(good) },
  )
  assert_eq(bag.length(), 1)
  assert_eq(bag.get("tree.newick"), Some(good))
  assert_false(bag.verify())
  // With the correct checksum the bag verifies.
  let sum = @src.md5_hex(good)
  let ok_manifest = sum + "  data/tree.newick\n"
  let ok_bag = @src.Datastore::load_bag(
    root="ok",
    manifest_text=ok_manifest,
    loader=fn(_) { Some(good) },
  )
  assert_true(ok_bag.verify())
}

///|
test "datastore_load_bag_errors" {
  // Missing payload.
  let result = try {
    let _ = @src.Datastore::load_bag(
      root="b",
      manifest_text="d41d8cd98f00b204e9800998ecf8427e  data/gone.fasta\n",
      loader=fn(_) { None },
    )
    false
  } catch {
    DatastoreError(msg) => msg.contains("Missing payload")
  }
  assert_true(result)
  // Malformed line.
  let result2 = try {
    let _ = @src.Datastore::load_bag(
      root="b",
      manifest_text="no-hash-here\n",
      loader=fn(_) { Some("x") },
    )
    false
  } catch {
    DatastoreError(msg) => msg.contains("Malformed")
  }
  assert_true(result2)
  // Non-MD5 checksum.
  let result3 = try {
    let _ = @src.Datastore::load_bag(
      root="b",
      manifest_text="abc123  data/x.txt\n",
      loader=fn(_) { Some("x") },
    )
    false
  } catch {
    DatastoreError(msg) => msg.contains("not MD5")
  }
  assert_true(result3)
  // Payload outside data/.
  let result4 = try {
    let _ = @src.Datastore::load_bag(
      root="b",
      manifest_text="d41d8cd98f00b204e9800998ecf8427e  etc/passwd\n",
      loader=fn(_) { Some("x") },
    )
    false
  } catch {
    DatastoreError(msg) => msg.contains("data/")
  }
  assert_true(result4)
  // Empty manifest.
  let result5 = try {
    let _ = @src.Datastore::load_bag(
      root="b",
      manifest_text="\n",
      loader=fn(_) { Some("x") },
    )
    false
  } catch {
    DatastoreError(msg) => msg.contains("Empty manifest")
  }
  assert_true(result5)
}

///|
test "datastore_fetch_cache" {
  let bag = @src.Datastore::new(root="cache")
  let mut loads = 0
  let fetch = fn(name : String) -> String {
    loads = loads + 1
    "content-of-" + name
  }
  let first = bag.fetch("entrez-1", fetch)
  assert_eq(first, "content-of-entrez-1")
  let second = bag.fetch("entrez-1", fetch)
  assert_eq(second, first)
  assert_eq(loads, 1)
  let _ = bag.fetch("entrez-2", fetch)
  assert_eq(loads, 2)
  assert_eq(bag.hits(), 1)
  assert_eq(bag.misses(), 2)
  // Fetched payloads are content-addressed like added ones.
  assert_eq(bag.checksum("entrez-1"), Some(@src.md5_hex("content-of-entrez-1")))
}

///|
test "datastore_bagit_and_bag_info" {
  let bag = @src.Datastore::new(root="b")
  assert_eq(
    bag.bagit_txt(),
    "BagIt-Version: 1.0\nTag-File-Character-Encoding: UTF-8\n",
  )
  let _ = bag.add("a.txt", "hello") // 5 octets
  let _ = bag.add("b.txt", "hi") // 2 octets
  bag.set_tag("External-Description", "synthetic test bag")
  bag.set_tag("Contact-Name", "BioSeqs")
  let info = bag.bag_info_txt()
  // Tags are sorted, followed by the agent and the Payload-Oxum.
  assert_true(info.contains("Contact-Name: BioSeqs\n"))
  assert_true(info.contains("External-Description: synthetic test bag\n"))
  assert_true(info.contains("Bag-Software-Agent: BioSeqs Datastore\n"))
  assert_true(info.contains("Payload-Oxum: 7.2\n"))
  assert_eq(bag.payload_oxum(), "7.2")
  assert_eq(bag.total_octets(), 7)
}

///|
test "datastore_total_octets_utf8" {
  let bag = @src.Datastore::new(root="b")
  // "café" is 5 UTF-8 bytes.
  let _ = bag.add("cafe.txt", "café")
  assert_eq(bag.total_octets(), 5)
  assert_eq(bag.payload_oxum(), "5.1")
}
