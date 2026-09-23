///|
/// MicroRNA target prediction (targetscan-style).
///
/// Seed types (canonical):
///   - 8mer:     perfect seed match (pos 2-8) + A at site position 1
///   - 7mer-m8:  perfect seed match (pos 2-8)
///   - 7mer-A1:  perfect seed match (pos 2-7) + A at site position 1
///   - 6mer:     perfect seed match (pos 2-7)

pub(all) enum MirSeedType {
  MirSeed8mer
  MirSeed7merM8
  MirSeed7merA1
  MirSeed6mer
}

pub fn mir_seed_type_label(st : MirSeedType) -> String {
  match st {
    MirSeed8mer => "8mer"
    MirSeed7merM8 => "7mer-m8"
    MirSeed7merA1 => "7mer-A1"
    MirSeed6mer => "6mer"
  }
}

fn mir_seed_type_weight(st : MirSeedType) -> Double {
  match st {
    MirSeed8mer => 1.0
    MirSeed7merM8 => 0.7
    MirSeed7merA1 => 0.5
    MirSeed6mer => 0.3
  }
}

pub(all) struct MirTargetSite {
  mir_id : String
  target_id : String
  utr_pos : Int
  site_len : Int
  seed_type : MirSeedType
  score : Double
  matched_utr : String
  mirna_seed : String
}

pub(all) struct MirTargetResult {
  mir_id : String
  target_id : String
  sites : Array[MirTargetSite]
  total_score : Double
}

pub fn mir_extract_seed(mirna : String) -> String {
  let s = mirna.to_upper().replace(old="U", new="T")
  let n = s.length()
  let start = 1
  let end = if n >= 8 { 8 } else { n }
  if end <= start { return "" }
  substring(s, start, end - start)
}

pub fn mir_extract_seed6(mirna : String) -> String {
  let s = mirna.to_upper().replace(old="U", new="T")
  let n = s.length()
  let start = 1
  let end = if n >= 7 { 7 } else { n }
  if end <= start { return "" }
  substring(s, start, end - start)
}

fn mir_complement_base(c : Char) -> Char {
  match c {
    'A' => 'T'
    'T' => 'A'
    'U' => 'A'
    'G' => 'C'
    'C' => 'G'
    'N' => 'N'
    _ => 'N'
  }
}

fn mir_reverse_complement(seq : String) -> String {
  let s = seq.to_upper()
  let n = s.length()
  let mut result = ""
  let mut i = n - 1
  while i >= 0 {
    let c = s.unsafe_get(i).unsafe_to_char()
    result = result + mir_complement_base(c).to_string()
    i = i - 1
  }
  result
}

fn mir_wc_pair(a : Char, b : Char) -> Bool {
  let au = (a == 'A' && b == 'T') || (a == 'T' && b == 'A')
  let au2 = (a == 'A' && b == 'U') || (a == 'U' && b == 'A')
  let gc = (a == 'G' && b == 'C') || (a == 'C' && b == 'G')
  au || au2 || gc
}

fn mir_wobble_pair(a : Char, b : Char) -> Bool {
  if mir_wc_pair(a, b) { return true }
  (a == 'G' && b == 'T') || (a == 'T' && b == 'G') || (a == 'G' && b == 'U') || (a == 'U' && b == 'G')
}

fn mir_local_au(utr : String, center : Int, upstream : Int, downstream : Int) -> Double {
  let n = utr.length()
  let lo = if center - upstream >= 0 { center - upstream } else { 0 }
  let hi = if center + downstream < n { center + downstream } else { n - 1 }
  if hi < lo { return 0.0 }
  let mut au_count = 0
  let mut total = 0
  for i in lo..=hi {
    let c = utr.unsafe_get(i).unsafe_to_char()
    if c == 'A' || c == 'T' || c == 'U' { au_count = au_count + 1 }
    total = total + 1
  }
  if total == 0 { 0.0 } else { au_count.to_double() / total.to_double() }
}

fn mir_position_score(site_pos : Int) -> Double {
  let dist = site_pos + 1
  if dist < 8 { return 0.3 }
  if dist <= 30 { return 1.0 }
  if dist <= 100 { return 0.8 }
  if dist <= 500 { return 0.6 }
  0.3
}

fn mir_context_score(seed_type : MirSeedType, utr : String, site_pos : Int, site_len : Int) -> Double {
  let base = mir_seed_type_weight(seed_type)
  let au = mir_local_au(utr, site_pos + site_len / 2, 30, 30)
  base * (0.5 + au) * mir_position_score(site_pos)
}

pub fn mir_find_seed_matches(
  mirna : String,
  utr : String,
  mir_id? : String = "",
  target_id? : String = "",
  allow_gu_wobble? : Bool = false,
  min_seed_type? : Int = 0,
) -> Array[MirTargetSite] {
  let mir = mirna.to_upper().replace(old="U", new="T")
  let u = utr.to_upper().replace(old="U", new="T")
  let n_mir = mir.length()
  let n_utr = u.length()
  let sites : Array[MirTargetSite] = Array(capacity=0)
  if n_mir < 7 || n_utr < 6 { return sites }

  let seed7 = mir_extract_seed(mir)
  let seed6 = mir_extract_seed6(mir)
  let rc7 = mir_reverse_complement(seed7)
  let rc6 = mir_reverse_complement(seed6)
  let mir_pos8 = if n_mir >= 8 { mir.unsafe_get(7).unsafe_to_char() } else { 'N' }
  let found_7mer : Map[Int, Unit] = Map([], capacity=64)

  // 7-mer seed matches
  if n_mir >= 8 && rc7.length() == 7 {
    let mut i = 0
    while i <= n_utr - 7 {
      let utr_region = substring(u, i, 7)
      let matches = if allow_gu_wobble {
        let mut ok = true
        for k in 0..<7 {
          let uc = utr_region.unsafe_get(k).unsafe_to_char()
          let mc = rc7.unsafe_get(k).unsafe_to_char()
          if !mir_wobble_pair(uc, mc) { ok = false; break }
        }
        ok
      } else { utr_region == rc7 }
      if matches {
        let site_pos = i - 1
        let has_a_anchor = if i > 0 { u.unsafe_get(i - 1).unsafe_to_char() == 'A' } else { false }
        let seed_type = if has_a_anchor { MirSeed8mer } else { MirSeed7merM8 }
        let st_rank = if has_a_anchor { 3 } else { 2 }
        if st_rank >= min_seed_type {
          let site_len = 8
          let matched = if i > 0 { substring(u, i - 1, 8) } else { utr_region }
          sites.push(MirTargetSite::{
            mir_id: mir_id, target_id: target_id, utr_pos: site_pos,
            site_len: site_len, seed_type: seed_type,
            score: mir_context_score(seed_type, u, site_pos, site_len),
            matched_utr: matched, mirna_seed: seed7,
          })
        }
        found_7mer.set(i, ())
      }
      i = i + 1
    }
  }

  // 6-mer seed matches
  if rc6.length() == 6 {
    let mut i = 0
    while i <= n_utr - 6 {
      if found_7mer.contains(i) { i = i + 1; continue }
      let utr_region = substring(u, i, 6)
      let matches = if allow_gu_wobble {
        let mut ok = true
        for k in 0..<6 {
          let uc = utr_region.unsafe_get(k).unsafe_to_char()
          let mc = rc6.unsafe_get(k).unsafe_to_char()
          if !mir_wobble_pair(uc, mc) { ok = false; break }
        }
        ok
      } else { utr_region == rc6 }
      if matches {
        let pos8_match = if i + 6 < n_utr && n_mir >= 8 {
          let uc = u.unsafe_get(i + 6).unsafe_to_char()
          mir_wc_pair(uc, mir_pos8)
        } else { false }
        let has_a_anchor = if i > 0 { u.unsafe_get(i - 1).unsafe_to_char() == 'A' } else { false }
        let seed_type = if pos8_match && has_a_anchor { MirSeed8mer }
          else if pos8_match { MirSeed7merM8 }
          else if has_a_anchor { MirSeed7merA1 }
          else { MirSeed6mer }
        let st_rank = match seed_type {
          MirSeed6mer => 0
          MirSeed7merA1 => 1
          MirSeed7merM8 => 2
          MirSeed8mer => 3
        }
        if st_rank >= min_seed_type {
          let site_pos = i - 1
          let site_len = if pos8_match { 8 } else if has_a_anchor { 7 } else { 6 }
          let matched = if site_pos >= 0 && site_pos + site_len <= n_utr {
            substring(u, site_pos, site_len)
          } else { utr_region }
          let actual_pos = if site_pos >= 0 { site_pos } else { i }
          sites.push(MirTargetSite::{
            mir_id: mir_id, target_id: target_id, utr_pos: actual_pos,
            site_len: site_len, seed_type: seed_type,
            score: mir_context_score(seed_type, u, actual_pos, site_len),
            matched_utr: matched,
            mirna_seed: if site_len >= 7 { seed7 } else { seed6 },
          })
        }
      }
      i = i + 1
    }
  }

  sites.sort_by(fn(a, b) { b.score.compare(a.score) })
  sites
}

pub fn mir_predict_targets(
  mirnas : Array[(String, String)],
  utrs : Array[(String, String)],
  min_score? : Double = 0.0,
  min_seed_type? : Int = 0,
  allow_gu_wobble? : Bool = false,
) -> Array[MirTargetResult] {
  let results : Array[MirTargetResult] = Array(capacity=0)
  for mir_entry in mirnas {
    let mir_id = mir_entry.0
    let mir_seq = mir_entry.1
    for utr_entry in utrs {
      let target_id = utr_entry.0
      let utr_seq = utr_entry.1
      let sites = mir_find_seed_matches(mir_seq, utr_seq,
        mir_id=mir_id, target_id=target_id,
        allow_gu_wobble=allow_gu_wobble, min_seed_type=min_seed_type)
      if sites.length() > 0 {
        let mut total = 0.0
        for s in sites { total = total + s.score }
        if total >= min_score {
          results.push(MirTargetResult::{
            mir_id: mir_id, target_id: target_id,
            sites: sites, total_score: total,
          })
        }
      }
    }
  }
  results.sort_by(fn(a, b) { b.total_score.compare(a.total_score) })
  results
}

pub fn mir_export_sites(results : Array[MirTargetResult]) -> String {
  let sb = StringBuilder()
  sb.write_string("mir_id\ttarget_id\tutr_pos\tsite_len\tseed_type\tscore\tmatched_utr\tmirna_seed\n")
  for r in results {
    for s in r.sites {
      sb.write_string(s.mir_id + "\t" + s.target_id + "\t" +
        s.utr_pos.to_string() + "\t" + s.site_len.to_string() + "\t" +
        mir_seed_type_label(s.seed_type) + "\t" + s.score.to_string() + "\t" +
        s.matched_utr + "\t" + s.mirna_seed + "\n")
    }
  }
  sb.to_string()
}

pub fn mir_count_by_seed_type(results : Array[MirTargetResult]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=4)
  counts.set("8mer", 0); counts.set("7mer-m8", 0)
  counts.set("7mer-A1", 0); counts.set("6mer", 0)
  for r in results {
    for s in r.sites {
      let label = mir_seed_type_label(s.seed_type)
      counts.set(label, counts.get(label).unwrap_or(0) + 1)
    }
  }
  counts
}
