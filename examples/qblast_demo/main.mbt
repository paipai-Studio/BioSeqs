///|
/// NCBI BLAST web service demo (port of Biopython Bio.Blast.NCBIWWW.qblast).
///
/// The client takes an injected transport, so this demo runs fully offline
/// against a mock server that speaks the same protocol as the real NCBI
/// endpoint: PUT returns an RID + RTOE, SearchInfo polls return WAITING
/// then READY, and the final Get returns the result document. In real use
/// the fetch callback would perform an HTTP GET.

///|
/// A tiny stateful mock server: one WAITING poll, then READY, then the
/// "result" document.
struct MockBlast {
  mut polls : Int
  requests : Array[String]
}

fn main {
  let mock : MockBlast = { polls: 0, requests: [] }

  let client = @src.NcbiBlastClient::new(
    now=fn() { 0.0 },
    fetch=fn(url) {
      mock.requests.push(url)
      if url.contains("CMD=Put") {
        // Real NCBI answers with an RID and an estimated wait (RTOE).
        "  RID = 9ZZKQ2A014  \n  RTOE = 15\n"
      } else if url.contains("FORMAT_OBJECT=SearchInfo") {
        let n = mock.polls
        mock.polls = n + 1
        if n < 1 {
          "Status=WAITING\n"
        } else {
          "Status=READY\n"
        }
      } else {
        "<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput>\n<BlastOutput>\n" +
        "  <BlastOutput_iterations>\n    ... mock hits ...\n" +
        "  </BlastOutput_iterations>\n</BlastOutput>\n"
      }
    },
  )

  // Inspect the URL Biopython would build for a blastn nr search.
  let url = try {
    @src.qblast_put_url(
      "blastn",
      "nr",
      "ACGTACGTACGTACGT",
      megablast=false,
      hitlist_size=100,
    )
  } catch {
    @src.NcbiBlastError(msg) => {
      println("URL build failed: " + msg)
      return
    }
  }
  println("PUT url:")
  println(url)
  println("")

  // Run the whole submit -> poll -> fetch workflow.
  let result = try {
    client.qblast("blastn", "nr", "ACGTACGTACGTACGT", hitlist_size=100)
  } catch {
    @src.NcbiBlastError(msg) => {
      println("qblast failed: " + msg)
      return
    }
  }

  println("requests issued: " + mock.requests.length().to_string())
  let mut i = 0
  while i < mock.requests.length() {
    let u = mock.requests[i]
    let kind = if u.contains("CMD=Put") {
      "submit"
    } else if u.contains("SearchInfo") {
      "poll  "
    } else {
      "result"
    }
    let shown = if u.length() > 96 {
      u[0:96].to_owned()
    } else {
      u
    }
    println("  " + kind + " : " + shown)
    i = i + 1
  }
  println("")
  println("result document:")
  println(result)
}
