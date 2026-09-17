///|
/// Tests for the Bio.ExPASy.ScanProsite port. URL construction (quote_plus
/// semantics), XML result parsing, and the injected-fetch scan path are all
/// exercised offline.

///|
test "scanprosite_url_basic" {
  assert_eq(
    @src.scanprosite_url(seq="MKWVTF"),
    "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=MKWVTF&output=xml",
  )
}

///|
test "scanprosite_url_quote_plus" {
  // quote_plus: space becomes '+', parens and plus are percent-encoded.
  assert_eq(
    @src.scanprosite_url(seq="A C+G"),
    "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=A+C%2BG&output=xml",
  )
}

///|
test "scanprosite_url_extra_params" {
  assert_eq(
    @src.scanprosite_url(seq="MKWVTF", params=[("db", "sp"), ("sig", "P-x(2)-G")]),
    "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=MKWVTF&output=xml&db=sp&sig=P-x%282%29-G",
  )
}

///|
test "scanprosite_url_mirror" {
  assert_eq(
    @src.scanprosite_url(seq="X", mirror="https://example.org"),
    "https://example.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=X&output=xml",
  )
}

///|
test "scanprosite_url_output_not_xml_raises" {
  let r = try {
    ignore(@src.scanprosite_url(seq="X", output="txt"))
    false
  } catch {
    ScanPrositeError(msg) => msg == "Only output='xml' is supported currently."
  }
  assert_true(r)
}

///|
test "scanprosite_read_sample" {
  let record = @src.scanprosite_read(@src.scanprosite_sample_xml()) catch {
    _ => abort("sample XML should parse")
  }
  assert_eq(record.n_match, 3)
  assert_eq(record.n_seq, 1)
  assert_eq(record.matches.length(), 3)
  let m = record.matches[0]
  assert_eq(m.sequence_ac, "P00750")
  assert_eq(m.sequence_id, "TPA_HUMAN")
  assert_eq(m.sequence_db, "swissprot")
  assert_eq(m.signature_ac, "PS00001")
  assert_eq(m.start, 114)
  assert_eq(m.stop, 135)
  assert_eq(m.level, "0")
  assert_eq(m.level_tag, "(0)")
  assert_eq(record.matches[1].signature_ac, "PS00134")
  assert_eq(record.matches[2].signature_ac, "PS00002")
}

///|
test "scanprosite_read_unknown_element_goes_to_extras" {
  let xml =
    "<scanprosite_response><matchset n_match=\"1\" n_seq=\"1\"><match><sequence_ac>X</sequence_ac><start>5</start><stop>9</stop><custom_tag>hello</custom_tag></match></matchset></scanprosite_response>"
  let record = @src.scanprosite_read(xml) catch {
    _ => abort("compact XML should parse")
  }
  assert_eq(record.n_match, 1)
  assert_eq(record.matches.length(), 1)
  assert_eq(record.matches[0].extras.get("custom_tag").unwrap_or(""), "hello")
}

///|
test "scanprosite_read_requires_root_tag" {
  let r = try {
    ignore(@src.scanprosite_read("<html>not scanprosite</html>"))
    false
  } catch {
    ScanPrositeError(_) => true
  }
  assert_true(r)
}

///|
test "scanprosite_scan_with_injected_fetch" {
  let calls : Array[String] = []
  let xml = @src.scanprosite_sample_xml()
  let result = @src.scanprosite_scan(
    fetch=fn(url : String, _post : String?) -> String {
      calls.push(url)
      xml
    },
    seq="MKWVTF",
    params=[("db", "sp")],
  )
  assert_eq(
    calls[0],
    "https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=MKWVTF&output=xml&db=sp",
  )
  assert_true(result.length() > 0)
  let record = @src.scanprosite_read(result) catch {
    _ => abort("fetched result should parse")
  }
  assert_eq(record.n_match, 3)
}