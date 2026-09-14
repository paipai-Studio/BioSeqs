///|
/// Tests for BiocFileCache — local file cache management.

///|
test "bfc add and get roundtrip" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t1")
  let entry = cache.add(rname="reads.fastq", data="@r1\nACGT\n")
  assert_true(entry.rid != "")
  assert_eq(entry.rname, "reads.fastq")
  assert_eq(entry.rtype, @src.BfcResourceType::Local)
  match cache.get(entry.rid) {
    Some(e) => assert_eq(e.rid, entry.rid)
    None => assert_true(false)
  }
}

///|
test "bfc hit and miss statistics" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t2")
  let e1 = cache.add(rname="a.txt", data="hello")
  let _ = cache.add(rname="b.txt", data="world")
  let _ = cache.get(e1.rid)
  let _ = cache.get("BFC999")
  assert_eq(cache.hits, 1)
  assert_eq(cache.misses, 3)
}

///|
test "bfc bfcrpath resolves path" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t3")
  let e = cache.add(rname="vcf.vcf.gz", data="##fileformat=VCFv4.2\n")
  match cache.bfcrpath(e.rid) {
    Some(p) => assert_true(p.contains("/tmp/t3/"))
    None => assert_true(false)
  }
  assert_true(cache.bfcrpath("NOPE") is None)
}

///|
test "bfc query by keyword" {
  let cache = @src.biocfilecache_sample_cache()
  let found = cache.bfcquery("genome")
  assert_eq(found.length(), 1)
  assert_eq(found[0].rname, "genome.fa")
  let none = cache.bfcquery("missing_file")
  assert_eq(none.length(), 0)
}

///|
test "bfc needs update compares etag" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t4")
  let web = cache.add(rname="remote.h5", data="HDF5", rtype=@src.BfcResourceType::Web, etag="r1")
  match cache.bfcneedsupdate(web.rid, "r1") {
    Some(false) => ()
    _ => assert_true(false)
  }
  match cache.bfcneedsupdate(web.rid, "r2") {
    Some(true) => ()
    _ => assert_true(false)
  }
  assert_true(cache.bfcneedsupdate("BFC404", "x") is None)
}

///|
test "bfc update changes checksum" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t5")
  let e = cache.add(rname="matrix.tsv", data="1\t2\n3\t4\n")
  let old = e.checksum
  let _ = cache.update(e.rid, data="1\t2\n3\t5\n")
  match cache.get(e.rid) {
    Some(updated) => assert_true(updated.checksum != old)
    None => assert_true(false)
  }
}

///|
test "bfc remove" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t6")
  let e = cache.add(rname="temp.txt", data="x")
  assert_true(cache.remove(e.rid))
  assert_false(cache.remove(e.rid))
  assert_eq(cache.length(), 0)
}

///|
test "bfc clean_cache sweeps expired" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t7")
  let _ = cache.add(rname="fresh.fa", data="A", expires="never")
  let _ = cache.add(rname="stale.fa", data="B", expires="2020-01-01")
  let removed = cache.clean_cache("2020-01-01")
  assert_eq(removed, 1)
  assert_eq(cache.length(), 1)
}

///|
test "bfc checksum deterministic" {
  let c1 = @src.bfc_checksum("ACGTACGT")
  let c2 = @src.bfc_checksum("ACGTACGT")
  let c3 = @src.bfc_checksum("ACGTACGA")
  assert_eq(c1, c2)
  assert_true(c1 != c3)
}

///|
test "bfc tags" {
  let cache = @src.BfcCache::new(cache_root="/tmp/t8")
  let e = cache.add(rname="anno.gff", data="##gff-version 3\n")
  assert_true(cache.set_tag(e.rid, "project", "demo"))
  match cache.get(e.rid) {
    Some(got) => match got.tags.get("project") {
      Some(v) => assert_eq(v, "demo")
      None => assert_true(false)
    }
    None => assert_true(false)
  }
}

///|
test "bfc info summary" {
  let cache = @src.biocfilecache_sample_cache()
  let s = cache.info()
  assert_true(s.contains("2 resources"))
}

///|
test "bfc resource type enum" {
  assert_eq(@src.BfcResourceType::Web, @src.BfcResourceType::Web)
  assert_true(@src.BfcResourceType::Local != @src.BfcResourceType::Remote)
}