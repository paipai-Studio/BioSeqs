///|
/// Bio.Datastore - MD5-keyed, BagIt-style local data cache.
///
/// Inspired by the bag ("BagIt") packaging conventions used to store and
/// transfer arbitrary digital content: a bag has just enough structure to
/// enclose descriptive metadata "tags" and a file "payload", and every
/// payload file is recorded in a checksum manifest (MD5 by default).
///
/// This module provides an in-memory bag so that code can:
///   - add payloads and obtain their MD5 checksum (content addressing),
///   - build/parse a `manifest-md5.txt` and `bag-info.txt`/`bagit.txt`,
///   - verify that payload contents still match the manifest,
///   - use the bag as a fetch-through cache (like Bio.Entrez style caching:
///     look up the local store first, download on miss, keep for next time).

///|
/// Error raised by Datastore operations for invalid bags or manifests.
pub suberror DatastoreError {
  DatastoreError(String)
}

///|
fn[T] datastore_fail(msg : String) -> T raise DatastoreError {
  raise DatastoreError(msg)
}

///|
/// One payload entry inside a bag.
pub struct DatastoreEntry {
  name : String
  data : String
  checksum : String
}

///|
/// Construct a DatastoreEntry.
pub fn DatastoreEntry::new(
  name~ : String,
  data~ : String,
  checksum~ : String,
) -> DatastoreEntry {
  DatastoreEntry::{ name, data, checksum, }
}

///|
/// Get the payload name.
pub fn DatastoreEntry::name(self : DatastoreEntry) -> String {
  self.name
}

///|
/// Get the payload text.
pub fn DatastoreEntry::data(self : DatastoreEntry) -> String {
  self.data
}

///|
/// Get the MD5 checksum recorded for the payload.
pub fn DatastoreEntry::checksum(self : DatastoreEntry) -> String {
  self.checksum
}

///|
/// MD5-keyed bag cache with tag metadata and fetch statistics.
pub struct Datastore {
  priv root : String
  priv entries : Map[String, String]
  priv manifest : Map[String, String]
  priv tags : Map[String, String]
  priv mut cache_hits : Int
  priv mut cache_misses : Int
}

///|
/// Create an empty bag rooted at `root`.
pub fn Datastore::new(root~ : String) -> Datastore {
  {
    root,
    entries: Map([], capacity=0),
    manifest: Map([], capacity=0),
    tags: Map([], capacity=0),
    cache_hits: 0,
    cache_misses: 0,
  }
}

///|
/// Get the bag root name.
pub fn Datastore::root(self : Datastore) -> String {
  self.root
}

///|
/// Number of payloads held in the bag.
pub fn Datastore::length(self : Datastore) -> Int {
  self.entries.size()
}

///|
/// Sorted payload names.
pub fn Datastore::names(self : Datastore) -> Array[String] {
  let keys : Array[String] = []
  for key in self.entries.keys() {
    keys.push(key)
  }
  datastore_sort_strings(keys)
}

///|
/// Set a bag-info tag (e.g. "Bag-Software-Agent", "External-Description").
pub fn Datastore::set_tag(
  self : Datastore,
  key : String,
  value : String,
) -> Unit {
  self.tags[key] = value
}

///|
/// Look up a bag-info tag.
pub fn Datastore::get_tag(self : Datastore, key : String) -> String? {
  self.tags.get(key)
}

///|
/// Sorted (key, value) tag pairs.
pub fn Datastore::tags(self : Datastore) -> Array[(String, String)] {
  let pairs : Array[(String, String)] = []
  for key in self.tags.keys() {
    match self.tags.get(key) {
      Some(v) => pairs.push((key, v))
      None => ()
    }
  }
  datastore_sort_pairs(pairs)
}

///|
/// Store a payload under `name`; returns its MD5 checksum.
/// Re-adding the same name replaces the payload and manifest entry.
pub fn Datastore::add(self : Datastore, name : String, data : String) -> String {
  let checksum = md5_hex(data)
  self.entries[name] = data
  self.manifest[name] = checksum
  checksum
}

///|
/// Look up a payload by name.
pub fn Datastore::get(self : Datastore, name : String) -> String? {
  self.entries.get(name)
}

///|
/// Whether a payload with this name is held.
pub fn Datastore::contains(self : Datastore, name : String) -> Bool {
  self.entries.contains(name)
}

///|
/// Remove a payload; returns true when it existed.
pub fn Datastore::remove(self : Datastore, name : String) -> Bool {
  let existed = self.entries.contains(name)
  let _ = self.entries.remove(name)
  let _ = self.manifest.remove(name)
  existed
}

///|
/// Fetch-through cache: return the stored payload on hit; otherwise call
/// `load(name)` to obtain the payload, store it, and return it.
pub fn Datastore::fetch(
  self : Datastore,
  name : String,
  load : (String) -> String,
) -> String {
  match self.entries.get(name) {
    Some(data) => {
      self.cache_hits = self.cache_hits + 1
      data
    }
    None => {
      self.cache_misses = self.cache_misses + 1
      let data = load(name)
      let _ = self.add(name, data)
      data
    }
  }
}

///|
/// Number of cache hits recorded so far.
pub fn Datastore::hits(self : Datastore) -> Int {
  self.cache_hits
}

///|
/// Number of cache misses recorded so far.
pub fn Datastore::misses(self : Datastore) -> Int {
  self.cache_misses
}

///|
/// The MD5 checksum recorded for a payload.
pub fn Datastore::checksum(self : Datastore, name : String) -> String? {
  self.manifest.get(name)
}

///|
/// A single manifest entry with its payload (when present).
pub fn Datastore::entry(self : Datastore, name : String) -> DatastoreEntry? {
  match (self.entries.get(name), self.manifest.get(name)) {
    (Some(data), Some(checksum)) =>
      Some(DatastoreEntry::{ name, data, checksum, })
    _ => None
  }
}

///|
/// Total payload size in bytes (UTF-8), for BagIt's Payload-Oxum.
pub fn Datastore::total_octets(self : Datastore) -> Int {
  let mut octets = 0
  for name in self.entries.keys() {
    match self.entries.get(name) {
      Some(data) => octets = octets + datastore_utf8_length(data)
      None => ()
    }
  }
  octets
}

///|
/// Payload-Oxum tag value "octets.streams" (BagIt bag-info.txt).
pub fn Datastore::payload_oxum(self : Datastore) -> String {
  datastore_int_to_string(self.total_octets()) +
  "." +
  datastore_int_to_string(self.entries.size())
}

///|
/// The `bagit.txt` tag file contents.
pub fn Datastore::bagit_txt(self : Datastore) -> String {
  "BagIt-Version: 1.0\nTag-File-Character-Encoding: UTF-8\n"
}

///|
/// The `bag-info.txt` tag file contents: user tags (sorted), the
/// software agent, and Payload-Oxum.
pub fn Datastore::bag_info_txt(self : Datastore) -> String {
  let mut text = ""
  for tag in self.tags() {
    let (key, value) = tag
    text = text + key + ": " + value + "\n"
  }
  text = text + "Bag-Software-Agent: BioSeqs Datastore\n"
  text = text + "Payload-Oxum: " + self.payload_oxum() + "\n"
  text
}

///|
/// Manifest lines "checksum  data/name" in BagIt manifest-md5.txt format,
/// sorted by payload name.
pub fn Datastore::manifest_lines(self : Datastore) -> Array[String] {
  let names : Array[String] = []
  for name in self.manifest.keys() {
    names.push(name)
  }
  let sorted = datastore_sort_strings(names)
  let lines : Array[String] = []
  for name in sorted {
    match self.manifest.get(name) {
      Some(checksum) => lines.push(checksum + "  data/" + name)
      None => ()
    }
  }
  lines
}

///|
/// The `manifest-md5.txt` file contents.
pub fn Datastore::manifest_md5_txt(self : Datastore) -> String {
  let mut text = ""
  for line in self.manifest_lines() {
    text = text + line + "\n"
  }
  text
}

///|
/// Verify that every payload still matches the checksum recorded in the
/// manifest. Payloads loaded from an external manifest whose contents do
/// not match the declared checksums make verification fail.
pub fn Datastore::verify(self : Datastore) -> Bool {
  for name in self.manifest.keys() {
    let declared = match self.manifest.get(name) {
      Some(v) => v
      None => return false
    }
    let data = match self.entries.get(name) {
      Some(v) => v
      None => return false
    }
    if md5_hex(data) != declared {
      return false
    }
  }
  self.manifest.size() == self.entries.size()
}

///|
/// Load a bag from a manifest-md5.txt text, fetching payloads through
/// `loader` (payload path -> contents). The declared checksums are kept
/// verbatim in the manifest, so a later `verify()` detects payloads whose
/// contents do not match the manifest.
pub fn Datastore::load_bag(
  root~ : String,
  manifest_text~ : String,
  loader~ : (String) -> String?,
) -> Datastore raise DatastoreError {
  let bag = Datastore::new(root~)
  let lines = datastore_split_lines(manifest_text)
  let mut loaded = 0
  for line in lines {
    let trimmed = datastore_rstrip(line)
    if trimmed.length() == 0 {
      continue
    }
    let parts = split_by_char(trimmed, ' '.to_int())
    let mut split = -1
    for i in 0..<parts.length() {
      if parts[i].length() == 0 {
        split = i
        break
      }
    }
    if split <= 0 {
      datastore_fail("Malformed manifest line: " + trimmed)
    }
    let checksum = parts[0]
    let path = parts[split + 1]
    if checksum.length() != 32 {
      datastore_fail("Manifest checksum is not MD5: " + checksum)
    }
    let prefix_len = "data/".length()
    if path.length() <= prefix_len || path[0:prefix_len].to_string() != "data/" {
      datastore_fail("Manifest payload must live under data/: " + path)
    }
    let name = path[prefix_len:path.length()].to_string()
    match loader(name) {
      Some(data) => {
        bag.entries[name] = data
        bag.manifest[name] = checksum
        loaded = loaded + 1
      }
      None => datastore_fail("Missing payload for manifest entry: " + name)
    }
  }
  if loaded == 0 {
    datastore_fail("Empty manifest")
  }
  bag
}

///|
/// Compute the MD5 checksum of a UTF-8 text (32 lowercase hex digits).
///
/// Reference values from RFC 1321:
///   md5_hex("")                      = "d41d8cd98f00b204e9800998ecf8427e"
///   md5_hex("a")                     = "0cc175b9c0f1b6a831c399e269772661"
///   md5_hex("abc")                   = "900150983246c27d70659a7e88e33d69"
///   md5_hex("message digest")        = "f96b697d7cb7938d525a2f31aaf161d0"
pub fn md5_hex(text : String) -> String {
  let bytes = datastore_utf8_bytes(text)
  let msg_len = bytes.length()
  // Pad: message + 0x80 + zeros, length congruent to 56 mod 64,
  // followed by the 64-bit little-endian bit length.
  let mut padded_len = msg_len + 9
  while padded_len % 64 != 0 {
    padded_len = padded_len + 1
  }
  let padded : Array[Int] = Array::make(padded_len, 0)
  for i in 0..<msg_len {
    padded[i] = bytes[i] & 0xFF
  }
  padded[msg_len] = 0x80
  let bit_len = msg_len * 8
  padded[padded_len - 8] = bit_len & 0xFF
  padded[padded_len - 7] = (bit_len >> 8) & 0xFF
  padded[padded_len - 6] = (bit_len >> 16) & 0xFF
  padded[padded_len - 5] = (bit_len >> 24) & 0xFF
  // High 32 bits of the length stay zero for message sizes in memory.
  let shifts : Array[Int] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5,
    9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16,
    23, 4, 11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
  ]
  let k = datastore_md5_k()
  let mut a0 = 1732584193 // 0x67452301
  let mut b0 = -271733879 // 0xefcdab89
  let mut c0 = -1732584194 // 0x98badcfe
  let mut d0 = 271733878 // 0x10325476
  let nblocks = padded_len / 64
  for block in 0..<nblocks {
    let m : Array[Int] = Array::make(16, 0)
    for j in 0..<16 {
      let off = block * 64 + j * 4
      m[j] = padded[off] |
        (padded[off + 1] << 8) |
        (padded[off + 2] << 16) |
        (padded[off + 3] << 24)
    }
    let mut a = a0
    let mut b = b0
    let mut c = c0
    let mut d = d0
    for i in 0..<64 {
      let f = if i < 16 {
        (b & c) | ((-1 - b) & d)
      } else if i < 32 {
        (d & b) | ((-1 - d) & c)
      } else if i < 48 {
        b ^ c ^ d
      } else {
        c ^ (b | (-1 - d))
      }
      let g = if i < 16 {
        i
      } else if i < 32 {
        (5 * i + 1) % 16
      } else if i < 48 {
        (3 * i + 5) % 16
      } else {
        7 * i % 16
      }
      let rotated = datastore_rotl(f + a + k[i] + m[g], shifts[i])
      a = d
      d = c
      c = b
      b = b + rotated
    }
    a0 = a0 + a
    b0 = b0 + b
    c0 = c0 + c
    d0 = d0 + d
  }
  let mut out = ""
  for word in [a0, b0, c0, d0] {
    // Digest words are serialized little-endian.
    for shift in [0, 8, 16, 24] {
      out = out + datastore_byte_hex((word >> shift) & 0xFF)
    }
  }
  out
}

///|
/// MD5 per-round constants K[i] = floor(|sin(i + 1)| * 2^32) as Int32.
fn datastore_md5_k() -> Array[Int] {
  let k : Array[Int] = Array::make(64, 0)
  for i in 0..<64 {
    let v = @math.sin((i + 1).to_double()).abs() * 4294967296.0
    let vf = @math.floor(v)
    k[i] = if vf >= 2147483648.0 {
      (vf - 4294967296.0).to_int()
    } else {
      vf.to_int()
    }
  }
  k
}

///|
/// Rotate a 32-bit word left (Int arithmetic wraps modulo 2^32).
fn datastore_rotl(x : Int, n : Int) -> Int {
  (x << n) | datastore_lshr(x, 32 - n)
}

///|
/// Logical right shift for n in 1..31 built from the arithmetic shift.
fn datastore_lshr(x : Int, n : Int) -> Int {
  (x >> n) & ((1 << (32 - n)) - 1)
}

///|
/// One byte as two lowercase hex digits.
fn datastore_byte_hex(b : Int) -> String {
  let digits = "0123456789abcdef"
  let hi = (b >> 4) & 0xF
  let lo = b & 0xF
  digits[hi:hi + 1].to_string() + digits[lo:lo + 1].to_string()
}

///|
/// Encode text as UTF-8 bytes (code units), combining surrogate pairs.
fn datastore_utf8_bytes(text : String) -> Array[Int] {
  let out : Array[Int] = []
  let n = text.length()
  let mut i = 0
  while i < n {
    let unit = text.unsafe_get(i).to_int()
    if unit >= 0xD800 && unit < 0xDC00 && i + 1 < n {
      let low = text.unsafe_get(i + 1).to_int()
      if low >= 0xDC00 && low < 0xE000 {
        let cp = 0x10000 + ((unit - 0xD800) << 10) + (low - 0xDC00)
        datastore_push_utf8(out, cp)
        i = i + 2
        continue
      }
    }
    datastore_push_utf8(out, unit)
    i = i + 1
  }
  out
}

///|
fn datastore_push_utf8(out : Array[Int], cp : Int) -> Unit {
  if cp < 0x80 {
    out.push(cp & 0xFF)
  } else if cp < 0x800 {
    out.push(0xC0 | (cp >> 6))
    out.push(0x80 | (cp & 0x3F))
  } else if cp < 0x10000 {
    out.push(0xE0 | (cp >> 12))
    out.push(0x80 | ((cp >> 6) & 0x3F))
    out.push(0x80 | (cp & 0x3F))
  } else {
    out.push(0xF0 | (cp >> 18))
    out.push(0x80 | ((cp >> 12) & 0x3F))
    out.push(0x80 | ((cp >> 6) & 0x3F))
    out.push(0x80 | (cp & 0x3F))
  }
}

///|
/// UTF-8 byte length of a text.
fn datastore_utf8_length(text : String) -> Int {
  datastore_utf8_bytes(text).length()
}

///|
/// Minimal non-negative integer to string.
fn datastore_int_to_string(v : Int) -> String {
  v.to_string()
}

///|
/// Trim trailing whitespace characters.
fn datastore_rstrip(text : String) -> String {
  let mut end = text.length()
  while end > 0 {
    let c = char_at(text, end - 1)
    if c == " " || c == "\t" || c == "\r" || c == "\n" {
      end = end - 1
    } else {
      break
    }
  }
  text[0:end].to_string()
}

///|
/// Split text into lines keeping content only.
fn datastore_split_lines(text : String) -> Array[String] {
  let lines : Array[String] = []
  let n = text.length()
  let mut start = 0
  let mut i = 0
  while i < n {
    let c = char_at(text, i)
    if c == "\n" {
      lines.push(text[start:i].to_string())
      start = i + 1
    }
    i = i + 1
  }
  if start < n {
    lines.push(text[start:n].to_string())
  }
  lines
}

///|
/// Lexicographic string less-than (the toolchain's `<` on String is not
/// lexicographic, so comparison is done over UTF-16 code units).
fn datastore_str_lt(a : String, b : String) -> Bool {
  let n1 = a.length()
  let n2 = b.length()
  let n = if n1 < n2 { n1 } else { n2 }
  let mut i = 0
  while i < n {
    let c1 = a.unsafe_get(i).to_int()
    let c2 = b.unsafe_get(i).to_int()
    if c1 != c2 {
      return c1 < c2
    }
    i = i + 1
  }
  n1 < n2
}

///|
/// Insertion sort for string arrays (ascending).
fn datastore_sort_strings(items : Array[String]) -> Array[String] {
  let arr = items.copy()
  let n = arr.length()
  for i in 1..<n {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && datastore_str_lt(key, arr[j]) {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
  arr
}

///|
/// Insertion sort for (key, value) pairs by key (ascending).
fn datastore_sort_pairs(
  items : Array[(String, String)],
) -> Array[(String, String)] {
  let arr = items.copy()
  let n = arr.length()
  for i in 1..<n {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && datastore_str_lt(key.0, arr[j].0) {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
  arr
}
