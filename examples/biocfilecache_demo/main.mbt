///|
fn main {
  println("=== BiocFileCache Demo ===")

  let cache = @src.BfcCache::new(cache_root="/tmp/demo_bfc")
  println("\n1. Add a local FASTA resource")
  let e1 = cache.add(
    rname="genome.fa",
    data=">chr1\nACGTACGTACGT\n",
    rtype=@src.BfcResourceType::Local,
  )
  println("   rid=" + e1.rid + ", rname=" + e1.rname)

  println("\n2. Add a remote resource with etag")
  let e2 = cache.add(
    rname="annotation.gtf.gz",
    data="chr1\tsource\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n",
    rtype=@src.BfcResourceType::Web,
    etag="v2024-01",
  )
  println("   rid=" + e2.rid + ", etag=" + e2.etag)

  println("\n3. Resolve cached file path (bfcrpath)")
  match cache.bfcrpath(e1.rid) {
    Some(p) => println("   " + e1.rid + " -> " + p)
    None => println("   not found")
  }

  println("\n4. Hit / miss statistics")
  let _ = cache.get(e1.rid)
  let _ = cache.get("BFC999")
  println("   " + cache.info())

  println("\n5. Query by keyword")
  let found = cache.bfcquery("genome")
  println("   found " + found.length().to_string() + " resource(s)")

  println("\n6. Check staleness (bfcneedsupdate)")
  match cache.bfcneedsupdate(e2.rid, "v2024-02") {
    Some(true) => println("   annotation.gtf.gz needs update")
    Some(false) => println("   annotation.gtf.gz is fresh")
    None => println("   not found")
  }

  println("\n7. Update resource content")
  let _ = cache.update(e2.rid, data="chr1\tsource\tgene\t1\t100\t.\t+\t.\tgene_id \"g2\";")
  match cache.get(e2.rid) {
    Some(updated) => println("   new checksum=" + updated.checksum)
    None => ()
  }

  println("\n8. Tag a resource")
  let _ = cache.set_tag(e1.rid, "project", "demo")
  match cache.get(e1.rid) {
    Some(t) => match t.tags.get("project") {
      Some(v) => println("   tag project=" + v)
      None => ()
    }
    None => ()
  }

  println("\n9. Remove a resource")
  println("   removed: " + cache.remove(e1.rid).to_string())
  println("   " + cache.info())

  println("\n10. Content checksum")
  println("   checksum('ACGT') = " + @src.bfc_checksum("ACGT"))
}