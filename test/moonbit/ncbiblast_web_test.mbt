///|
/// Tests for Bio.Blast.NCBIWWW.qblast port. All HTTP traffic goes through
/// an injected transport, so behaviour is verified offline against the
/// URL/parameter semantics of Biopython's NCBIWWW.

///|
/// Stateful mock BLAST server: records request URLs and serves canned
/// PUT / SearchInfo / Get responses.
struct BlastMock {
  mut info_calls : Int
  statuses : Array[String]
  requests : Array[String]
  result : String
}

///|
fn blast_mock(statuses : Array[String], result : String) -> BlastMock {
  { info_calls: 0, statuses, requests: [], result, }
}

///|
fn blast_client(mock : BlastMock) -> @src.NcbiBlastClient {
  @src.NcbiBlastClient::new(
    fetch=fn(url) {
      mock.requests.push(url)
      if url.contains("FORMAT_OBJECT=SearchInfo") {
        let idx = mock.info_calls
        mock.info_calls = idx + 1
        let status = if idx < mock.statuses.length() {
          mock.statuses[idx]
        } else {
          "READY"
        }
        "Status=" + status + "\n"
      } else if url.contains("CMD=Put") {
        "  RID = 1234AB  \n  RTOE = 30\n\n"
      } else {
        mock.result
      }
    },
    now=fn() { 0.0 },
  )
}

///|
test "qblast_quote_url_encoding" {
  assert_eq(@src.qblast_quote("ACGT123_.-~"), "ACGT123_.-~")
  assert_eq(@src.qblast_quote("ACGT ACGT"), "ACGT+ACGT")
  assert_eq(@src.qblast_quote("a>b"), "a%3Eb")
  assert_eq(@src.qblast_quote("100%"), "100%25")
}

///|
test "qblast_validate_program" {
  // Valid programs do not raise.
  @src.qblast_validate_program("blastn")
  @src.qblast_validate_program("tblastx")
  let raised = try {
    @src.qblast_validate_program("fasta")
    false
  } catch {
    NcbiBlastError(_) => true
  }
  assert_true(raised)
}

///|
test "qblast_put_url_basic_parameters" {
  let url = @src.qblast_put_url("blastn", "nr", "ACGTACGT") catch {
    NcbiBlastError(msg) => {
      println(msg)
      ""
    }
  }
  assert_true(url.starts_with("https://blast.ncbi.nlm.nih.gov/Blast.cgi?"))
  assert_true(url.contains("CMD=Put"))
  assert_true(url.contains("PROGRAM=blastn"))
  assert_true(url.contains("DATABASE=nr"))
  assert_true(url.contains("QUERY=ACGTACGT"))
  assert_true(url.contains("FORMAT_TYPE=XML"))
  assert_true(url.contains("HITLIST_SIZE=500"))
  assert_true(url.contains("EXPECT=10"))
  // Biopython defaults: descriptions/alignments always sent.
  assert_true(url.contains("DESCRIPTIONS=500"))
  assert_true(url.contains("ALIGNMENTS=500"))
  // Optional parameters absent by default.
  assert_true(!url.contains("WORD_SIZE"))
  assert_true(!url.contains("MATRIX_NAME"))
  assert_true(!url.contains("SERVICE="))
}

///|
test "qblast_put_url_optional_parameters" {
  let url = @src.qblast_put_url(
    "blastn",
    "nt",
    "ACGT ACGT",
    word_size=28,
    expect=0.001,
    entrez_query="Homo sapiens[ORGN]",
    ncbi_gi=true,
  ) catch {
    NcbiBlastError(msg) => {
      println(msg)
      ""
    }
  }
  assert_true(url.contains("WORD_SIZE=28"))
  assert_true(url.contains("NCBI_GI=Yes"))
  // Spaces in the query become '+'.
  assert_true(url.contains("QUERY=ACGT+ACGT"))
  // entrez query is percent-encoded.
  assert_true(url.contains("ENTREZ_QUERY=Homo+sapiens%5BORGN%5D"))
}

///|
test "qblast_put_url_megablast" {
  let url = @src.qblast_put_url("blastn", "nt", "ACGT", megablast=true) catch {
    NcbiBlastError(msg) => {
      println(msg)
      ""
    }
  }
  // Biopython keeps PROGRAM=blastn and adds SERVICE=megablast.
  assert_true(url.contains("PROGRAM=blastn"))
  assert_true(url.contains("SERVICE=megablast"))
}

///|
test "qblast_put_url_rejects_empty_database" {
  let raised = try {
    let _ = @src.qblast_put_url("blastn", "", "ACGT")
    false
  } catch {
    NcbiBlastError(_) => true
  }
  assert_true(raised)
}

///|
test "qblast_parse_put_response" {
  let (rid, rtoe) = @src.qblast_parse_put("  RID = ABC123  \n  RTOE = 45\n") catch {
    NcbiBlastError(_) => ("", -1)
  }
  assert_eq(rid, "ABC123")
  assert_eq(rtoe, 45)

  let missing = try {
    let _ = @src.qblast_parse_put("Status=WAITING\n")
    false
  } catch {
    NcbiBlastError(_) => true
  }
  assert_true(missing)
}

///|
test "qblast_parse_status_values" {
  assert_true(
    @src.qblast_parse_status("Status=WAITING\n")
    is @src.NcbiBlastStatus::Waiting,
  )
  assert_true(
    @src.qblast_parse_status("Status=READY\n") is @src.NcbiBlastStatus::Ready,
  )
  assert_true(
    @src.qblast_parse_status("Status=FAILED\n") is @src.NcbiBlastStatus::Failed,
  )
  assert_true(
    @src.qblast_parse_status("Status=UNKNOWN\n")
    is @src.NcbiBlastStatus::Unknown,
  )
}

///|
test "qblast_status_urls" {
  let info = @src.qblast_search_info_url("RID1")
  assert_true(info.contains("CMD=Get"))
  assert_true(info.contains("RID=RID1"))
  assert_true(info.contains("FORMAT_OBJECT=SearchInfo"))

  let result = @src.qblast_result_url("RID1")
  assert_true(result.contains("CMD=Get"))
  assert_true(result.contains("FORMAT_TYPE=XML"))

  let text_result = @src.qblast_result_url("RID1", format_type="Text")
  assert_true(text_result.contains("FORMAT_TYPE=Text"))
}

///|
test "ncbiblast_submit_uses_put_url" {
  let mock = blast_mock(["WAITING"], "<BlastOutput/>")
  let client = blast_client(mock)
  let (rid, rtoe) = client.submit("blastp", "nr", "MVLSPADKTN") catch {
    NcbiBlastError(msg) => {
      println(msg)
      ("", -1)
    }
  }
  assert_eq(rid, "1234AB")
  assert_eq(rtoe, 30)
  assert_eq(mock.requests.length(), 1)
  assert_true(mock.requests[0].contains("CMD=Put"))
  assert_true(mock.requests[0].contains("PROGRAM=blastp"))
  assert_true(mock.requests[0].contains("QUERY=MVLSPADKTN"))
}

///|
test "ncbiblast_qblast_full_workflow" {
  let mock = blast_mock(["WAITING", "READY"], "<BlastOutput>ok</BlastOutput>")
  let client = blast_client(mock)
  let result = client.qblast("blastn", "nt", "ACGTACGT") catch {
    NcbiBlastError(msg) => {
      println(msg)
      ""
    }
  }
  assert_eq(result, "<BlastOutput>ok</BlastOutput>")
  // 1 PUT + 2 SearchInfo polls + 1 result fetch.
  assert_eq(mock.requests.length(), 4)
  assert_true(mock.requests[0].contains("CMD=Put"))
  assert_true(mock.requests[1].contains("FORMAT_OBJECT=SearchInfo"))
  assert_true(mock.requests[2].contains("FORMAT_OBJECT=SearchInfo"))
  assert_true(mock.requests[3].contains("FORMAT_TYPE=XML"))
  assert_true(!mock.requests[3].contains("FORMAT_OBJECT"))
}

///|
test "ncbiblast_qblast_failed_search" {
  let mock = blast_mock(["FAILED"], "")
  let client = blast_client(mock)
  let raised = try {
    let _ = client.qblast("blastn", "nt", "ACGT")
    false
  } catch {
    NcbiBlastError(_) => true
  }
  assert_true(raised)
}

///|
test "ncbiblast_qblast_timeout_waiting" {
  let mock = blast_mock(["WAITING", "WAITING", "WAITING"], "")
  let client = blast_client(mock)
  let raised = try {
    let _ = client.qblast("blastn", "nt", "ACGT", max_polls=2)
    false
  } catch {
    NcbiBlastError(_) => true
  }
  assert_true(raised)
}
