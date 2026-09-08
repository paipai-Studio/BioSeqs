///|
/// SeqIO index demo (port of Biopython Bio.SeqIO.index):
/// scan a FASTA file once, map record ids to byte offsets, then fetch
/// individual records by id without re-parsing the whole file.
fn main {
  let fasta =
    ">BRCA1 human breast cancer gene\n" +
    "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG\n" +
    "CGTCCCTGGTCTTCAGGGACTAATATTATCATTTGTAGT\n" +
    ">TP53 tumor suppressor\n" +
    "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAA\n" +
    "ACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCC\n" +
    ">EGFR\n" +
    "ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTCGGG\n" +
    "CTCTGGAGGAAAAGAAAGTTTGCCAAGGCACGAGTAACAAG\n"

  // ---- Build the index ----------------------------------------------------
  let index = try {
    @src.seqio_index(fasta, "fasta")
  } catch {
    @src.SeqIOError(msg) => {
      println("indexing failed: " + msg)
      return
    }
  }

  println("indexed " + index.length().to_string() + " records")
  let keys = index.keys()
  let mut i = 0
  while i < keys.length() {
    println("  " + keys[i])
    i = i + 1
  }
  println("contains TP53: " + index.contains("TP53").to_string())

  // ---- Random access by id -----------------------------------------------
  let ids = ["EGFR", "BRCA1"]
  let mut j = 0
  while j < ids.length() {
    let id = ids[j]
    try {
      let raw = index.get_raw(fasta, id)
      let record = index.get(fasta, id)
      println("")
      println("--- record " + id + " (" + raw.length().to_string() + " raw bytes) ---")
      println("description : " + record.description)
      println(
        "sequence    : " +
        record.seq.to_string().substring(start=0, end=30) +
        "... (" + record.seq.to_string().length().to_string() + " bp)",
      )
    } catch {
      @src.SeqIOError(msg) => println("failed to fetch " + id + ": " + msg)
    }
    j = j + 1
  }
}
