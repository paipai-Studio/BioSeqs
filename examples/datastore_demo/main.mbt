///|
/// Bio.Datastore Demo.
///
/// Demonstrates the MD5-keyed BagIt-style local data cache: payload
/// registration with MD5 manifests, fetch-through caching with hit/miss
/// accounting, bag-info tags (Payload-Oxum), manifest serialization,
/// verification and bag loading.

///|
fn main raise {
  println("=== 1. MD5 (RFC 1321 test vectors) ===")
  println("md5(\"\")               = \{@src.md5_hex("")}")
  println("md5(\"a\")              = \{@src.md5_hex("a")}")
  println("md5(\"abc\")            = \{@src.md5_hex("abc")}")
  println("md5(\"message digest\") = \{@src.md5_hex("message digest")}")

  println("\n=== 2. Create a bag and register payloads ===")
  let bag = @src.Datastore::new(root="demo-bag")
  let _ = bag.add("sequences/fasta/hsa.fas", ">hsa\nMVLSPADKTNVKAAWGK\n")
  let _ = bag.add("sequences/fasta/mmu.fas", ">mmu\nMVLSGEDKSNIKAAWGK\n")
  let _ = bag.add("notes/readme.txt", "synthetic demo payload\n")
  bag.set_tag("Contact-Name", "Demo User")
  bag.set_tag("External-Identifier", "demo:sequences-001")
  println("root:      \{bag.root()}")
  println("payloads:  \{bag.names()}")
  println("tags:      \{bag.tags()}")

  println("\n=== 3. Fetch-through cache ===")
  // first fetch: cache miss -> the loader callback provides the payload
  let first = bag.fetch("structures/pdb/1sct.pdb", fn(name) {
    println("  (loading \{name} from the remote source)")
    "ATOM      1  N   ALA A   1\nEND\n"
  })
  // second fetch: cache hit -> no loader call
  let second = bag.fetch("structures/pdb/1sct.pdb", fn(name) {
    println("  (this loader call should never happen: \{name})")
    ""
  })
  println("fetched payload stable: \{first == second}")
  println("cache hits: \{bag.hits()}, cache misses: \{bag.misses()}")

  println("\n=== 4. Manifest and BagIt files ===")
  println("payload-oxum: \{bag.payload_oxum()}")
  for line in bag.manifest_lines() {
    println("manifest: \{line}")
  }
  println("--- bagit.txt ---")
  println(bag.bagit_txt())
  println("--- bag-info.txt ---")
  println(bag.bag_info_txt())

  println("\n=== 5. Verify ===")
  println("bag verifies: \{bag.verify()}")
  let entry = bag.entry("notes/readme.txt")
  let checksum = match entry {
    Some(e) => e.checksum()
    None => "(missing)"
  }
  println("stored checksum of notes/readme.txt: \{checksum}")

  println("\n=== 6. Load a bag from an external manifest ===")
  // a consistent manifest: declare the real md5 of the provided payload
  let data_one = ">one\nAAAA\n"
  let good_manifest = @src.md5_hex(data_one) + "  data/seqs/one.fas\n"
  let loaded = @src.Datastore::load_bag(
    root="external-bag",
    manifest_text=good_manifest,
    loader=fn(_name) { Some(data_one) },
  )
  println("loaded payload names: \{loaded.names()}")
  println("consistent bag verifies: \{loaded.verify()}")
  // a manifest with a stale checksum makes verification fail
  let stale_manifest = "0cc175b9c0f1b6a831c399e269772661  data/seqs/one.fas\n"
  let mismatched = @src.Datastore::load_bag(
    root="stale-bag",
    manifest_text=stale_manifest,
    loader=fn(_name) { Some(data_one) },
  )
  println("stale manifest verifies: \{mismatched.verify()}")

  println("\nDone.")
}
