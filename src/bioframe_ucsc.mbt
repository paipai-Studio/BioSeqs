/// bioframe_ucsc.mbt — Fetch chromosome sizes / centromeres from UCSC.
///
/// Bioframe equivalent: `fetch_chromsizes`, `fetch_centromeres`.
///
/// MoonBit has no standard HTTP client yet, so these functions offer
/// two paths:
///   1. `fetch_url` callback (optional) — callers on a platform that
///      supplies HTTP transport (e.g. WASM host, wasm-gc + browser) can
///      plug one in.
///   2. Built-in fallback data — every assembly registered in
///      bioframe_assemblies.mbt ships its chromsizes and centromeres
///      as static constants, so these functions always return useful
///      values without needing network access.
///
/// UCSC goldenPath FTP base URLs:
///   https://hgdownload.soe.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes
///   https://hgdownload.soe.ucsc.edu/goldenPath/<build>/bigZips/<build>.centromeres

///|
/// UCSC goldenPath base (can be overridden).
pub const BF_UCSC_GOLDENPATH_BASE : String =
  "https://hgdownload.soe.ucsc.edu/goldenPath";

///|
/// Optional external fetch callback type alias.
/// Takes a URL and returns the response body as a String, or None on error.
pub type FetchUrlFn = (String) -> String?

///|
/// Construct the UCSC .chrom.sizes URL for an assembly.
pub fn bf_ucsc_chromsizes_url(assembly : String) -> String {
  BF_UCSC_GOLDENPATH_BASE +
  "/" +
  assembly +
  "/bigZips/" +
  assembly +
  ".chrom.sizes"
}

///|
/// Construct the UCSC centromeres URL for an assembly.
/// Not every assembly has this file — caller should check fallback.
pub fn bf_ucsc_centromeres_url(assembly : String) -> String {
  BF_UCSC_GOLDENPATH_BASE +
  "/" +
  assembly +
  "/bigZips/" +
  assembly +
  ".centromeres"
}

///|
/// Fetch chromsizes from UCSC with an optional fetch callback.
///
/// - If `fetch_url` is provided, tries to download from the UCSC URL.
///   Parses the two-column TSV (name, length).
/// - Otherwise (or if download fails), falls back to the built-in
///   chromsizes embedded in `bioframe_assemblies.mbt`.
/// - If the assembly is unknown even in the fallback, returns an empty
///   Array (never crashes).
///
/// Returns Array[(chrom, length)].
pub fn bf_fetch_chromsizes_from_ucsc(
  assembly : String,
  fetch_url~ : FetchUrlFn? = None,
) -> Array[(String, Int32)] {
  // Try remote
  match fetch_url {
    Some(fetcher) => {
      let url = bf_ucsc_chromsizes_url(assembly)
      match fetcher(url) {
        Some(body) => {
          let parsed = bf_parse_chromsizes_tsv(body)
          if parsed.length() > 0 {
            return parsed
          }
        }
        None => ()
      }
    }
    None => ()
  }

  // Fallback to built-in
  bf_assembly_chromsizes(assembly)
}

///|
/// Fetch centromere positions from UCSC. Same fallback behaviour as
/// bf_fetch_chromsizes_from_ucsc.
///
/// Returns Array[(chrom, start, end)].
pub fn bf_fetch_centromeres_from_ucsc(
  assembly : String,
  fetch_url~ : FetchUrlFn? = None,
) -> Array[(String, Int32, Int32)] {
  match fetch_url {
    Some(fetcher) => {
      let url = bf_ucsc_centromeres_url(assembly)
      match fetcher(url) {
        Some(body) => {
          let parsed = bf_parse_centromeres_tsv(body)
          if parsed.length() > 0 {
            return parsed
          }
        }
        None => ()
      }
    }
    None => ()
  }

  bf_assembly_centromeres(assembly)
}

///|
/// Parse UCSC chrom.sizes TSV (two columns: chrom<TAB>length).
fn bf_parse_chromsizes_tsv(body : String) -> Array[(String, Int)] {
  let out : Array[(String, Int)] = []
  let mut sb_line = StringBuilder()
  for c in body {
    if c == '\n' {
      let line = sb_line.to_string()
      sb_line = StringBuilder()
      if line.length() == 0 {
        continue
      }
      let mut sb = StringBuilder()
      let mut chrom : String = ""
      let mut found_tab = false
      for c2 in line {
        if c2 == '\t' {
          chrom = sb.to_string()
          sb = StringBuilder()
          found_tab = true
        } else {
          sb.write_char(c2)
        }
      }
      if found_tab && chrom.length() > 0 {
        let l = parse_int(sb.to_string())
        out.push((chrom, l))
      }
    } else {
      sb_line.write_char(c)
    }
  }
  out
}

///|
/// Parse UCSC centromeres TSV (three columns: chrom<TAB>start<TAB>end).
fn bf_parse_centromeres_tsv(body : String) -> Array[(String, Int, Int)] {
  let out : Array[(String, Int, Int)] = []
  let mut sb_line = StringBuilder()
  for c in body {
    if c == '\n' {
      let line = sb_line.to_string()
      sb_line = StringBuilder()
      if line.length() == 0 {
        continue
      }
      let mut col : Int32 = 0
      let mut chrom : String = ""
      let mut sb = StringBuilder()
      let mut s = 0
      let mut e = 0
      for c2 in line {
        if c2 == '\t' {
          let val = sb.to_string()
          sb = StringBuilder()
          if col == 0 { chrom = val }
          else if col == 1 { s = parse_int(val) }
          col = col + 1
        } else {
          sb.write_char(c2)
        }
      }
      // last col
      e = parse_int(sb.to_string())
      if col >= 2 && chrom.length() > 0 {
        out.push((chrom, s, e))
      }
    } else {
      sb_line.write_char(c)
    }
  }
  out
}
