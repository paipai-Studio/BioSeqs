///|
/// Local file cache management, ported from Bioconductor's
/// `BiocFileCache` package.
///
/// The reference package manages a per-user cache directory: each
/// resource is registered with a unique `rid`, an `rname`, the local
/// file path (`fpath`), a resource type (local / relative / remote /
/// web) and web metadata (`etag`, `expires`).  Lookups go through
/// `bfcrpath`; staleness is detected with `bfcneedsupdate`.  This
/// port keeps the same vocabulary over an in-memory store with
/// content checksums.

///|
/// Where a cached resource came from.
pub(all) enum BfcResourceType {
  Local
  Relative
  Remote
  Web
} derive(Eq, Show, Debug)
pub extend BfcResourceType with Eq::{equal, not_equal}
pub extend BfcResourceType with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A single cached resource record.
pub struct BfcEntry {
  rid : String
  rname : String
  fpath : String
  rtype : BfcResourceType
  etag : String
  expires : String
  checksum : String
  size : Int
  create_time : Int
  access_time : Int
  tags : Map[String, String]
} derive(Eq, Show, Debug)
pub extend BfcEntry with Eq::{equal, not_equal}
pub extend BfcEntry with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// The cache itself: entries keyed by rid plus hit/miss statistics.
pub struct BfcCache {
  cache_root : String
  entries : Map[String, BfcEntry]
  mut hits : Int
  mut misses : Int
  mut next_rid : Int
} derive(Eq, Show, Debug)
pub extend BfcCache with Eq::{equal, not_equal}
pub extend BfcCache with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn BfcCache::new(cache_root~ : String = "bfc_cache") -> BfcCache {
  BfcCache::{
    cache_root,
    entries: Map([], capacity=8),
    hits: 0,
    misses: 0,
    next_rid: 1,
  }
}

///|
/// FNV-1a 64-bit content checksum, rendered as 16 hex characters.
pub fn bfc_checksum(data : String) -> String {
  let mut hash : UInt64 = 0xcbf29ce484222325UL
  let mut i = 0
  let n = data.length()
  while i < n {
    let b = data[i]
    hash = hash ^ b.to_uint64()
    hash = hash * 0x100000001b3UL
    i += 1
  }
  hash.to_string()
}

///|
/// Add a resource to the cache.  Returns the created entry.
pub fn BfcCache::add(
  self : BfcCache,
  rname~ : String,
  data~ : String,
  rtype~ : BfcResourceType = BfcResourceType::Local,
  etag~ : String = "",
  expires~ : String = "",
) -> BfcEntry {
  let rid = "BFC" + self.next_rid.to_string()
  self.next_rid += 1
  let fpath = self.cache_root + "/" + rid
  let entry = BfcEntry::{
    rid,
    rname,
    fpath,
    rtype,
    etag,
    expires,
    checksum: bfc_checksum(data),
    size: data.length(),
    create_time: self.next_rid,
    access_time: self.next_rid,
    tags: Map([], capacity=4),
  }
  self.entries.set(rid, entry)
  self.misses += 1
  entry
}

///|
/// Look up an entry by rid; records a hit or miss.
pub fn BfcCache::get(self : BfcCache, rid : String) -> BfcEntry? {
  match self.entries.get(rid) {
    Some(e) => {
      self.hits += 1
      Some(e)
    }
    None => {
      self.misses += 1
      None
    }
  }
}

///|
/// Resolve a rid to its cached file path (`bfcrpath`).
pub fn BfcCache::bfcrpath(self : BfcCache, rid : String) -> String? {
  match self.entries.get(rid) {
    Some(e) => Some(e.fpath)
    None => None
  }
}

///|
/// Remove an entry from the cache.
pub fn BfcCache::remove(self : BfcCache, rid : String) -> Bool {
  match self.entries.get(rid) {
    Some(_) => {
      self.entries.remove(rid)
      true
    }
    None => false
  }
}

///|
/// Fuzzy search over rname and fpath (`bfcquery`).
pub fn BfcCache::bfcquery(
  self : BfcCache,
  keyword : String,
) -> Array[BfcEntry] {
  let results : Array[BfcEntry] = Array::new()
  for rid in self.entries.keys() {
    match self.entries.get(rid) {
      Some(e) => {
        if e.rname.contains(keyword) || e.fpath.contains(keyword) {
          results.push(e)
        }
      }
      None => ()
    }
  }
  results
}

///|
/// Check whether a resource needs refresh: true when the remote etag
/// differs from the cached one.
pub fn BfcCache::bfcneedsupdate(
  self : BfcCache,
  rid : String,
  remote_etag : String,
) -> Bool? {
  match self.entries.get(rid) {
    Some(e) => {
      if e.rtype is Web {
        Some(e.etag != remote_etag)
      } else {
        Some(false)
      }
    }
    None => None
  }
}

///|
/// Update the stored content (and checksum) for a resource.
pub fn BfcCache::update(
  self : BfcCache,
  rid : String,
  data~ : String,
  etag~ : String = "",
) -> Bool {
  match self.entries.get(rid) {
    Some(_) => {
      let e = self.entries[rid]
      let updated = BfcEntry::{
        rid: e.rid,
        rname: e.rname,
        fpath: e.fpath,
        rtype: e.rtype,
        etag: if etag != "" { etag } else { e.etag },
        expires: e.expires,
        checksum: bfc_checksum(data),
        size: data.length(),
        create_time: e.create_time,
        access_time: e.access_time,
        tags: e.tags,
      }
      self.entries.set(rid, updated)
      true
    }
    None => false
  }
}

///|
/// Attach a tag to a resource.
pub fn BfcCache::set_tag(
  self : BfcCache,
  rid : String,
  key : String,
  value : String,
) -> Bool {
  match self.entries.get(rid) {
    Some(_) => {
      let e = self.entries[rid]
      e.tags.set(key, value)
      true
    }
    None => false
  }
}

///|
/// Remove all entries whose `expires` field equals `marker` (used to
/// simulate `cleanbfc` expiry sweeping).
pub fn BfcCache::clean_cache(self : BfcCache, marker : String) -> Int {
  let doomed : Array[String] = Array::new()
  for rid in self.entries.keys() {
    match self.entries.get(rid) {
      Some(e) => {
        if e.expires == marker {
          doomed.push(rid)
        }
      }
      None => ()
    }
  }
  for rid in doomed {
    self.entries.remove(rid)
  }
  doomed.length()
}

///|
/// Number of cached resources.
pub fn BfcCache::length(self : BfcCache) -> Int {
  self.entries.length()
}

///|
/// One-line cache info summary.
pub fn BfcCache::info(self : BfcCache) -> String {
  "BiocFileCache at '" +
  self.cache_root +
  "': " +
  self.length().to_string() +
  " resources, " +
  self.hits.to_string() +
  " hits, " +
  self.misses.to_string() +
  " misses"
}

///|
/// Sample cache populated with two FASTA resources for demos/tests.
pub fn biocfilecache_sample_cache() -> BfcCache {
  let cache = BfcCache::new(cache_root="/tmp/demo_bfc")
  let _ = cache.add(
    rname="genome.fa",
    data=">chr1\nACGTACGT\n",
    rtype=BfcResourceType::Local,
  )
  let _ = cache.add(
    rname="annotation.gtf.gz",
    data="chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n",
    rtype=BfcResourceType::Web,
    etag="v1",
  )
  cache
}