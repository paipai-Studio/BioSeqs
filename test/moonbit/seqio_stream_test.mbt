///|
test "stream fasta matches parse_fasta" {
  let content = ">seq1 description here\nACGTACGTAC\nGTACGTACGT\n>seq2\n\nAAAA\n>seq3 short\nGGGGTTTTCCCC\n"
  let batch = @src.parse_fasta(content)
  let streamed : Array[@src.SeqRecord] = []
  @src.fasta_each_record(content, fn(r) { streamed.push(r) })
  assert_eq(streamed.length(), batch.length())
  assert_eq(streamed.length(), 3)
  for i in 0..<streamed.length() {
    assert_eq(streamed[i].id, batch[i].id)
    assert_eq(streamed[i].description, batch[i].description)
    assert_eq(streamed[i].seq.to_string(), batch[i].seq.to_string())
  }
  assert_eq(streamed[0].id, "seq1")
  assert_eq(streamed[0].description, "seq1 description here")
  assert_eq(streamed[0].seq.to_string(), "ACGTACGTACGTACGTACGT")
}

///|
test "stream fasta cursor resume" {
  let content = ">a\nAAAA\n>b\nCCCC\n>c\nTTTT\n"
  let cursor = @src.SeqioCursor::new()
  let first = @src.fasta_next_record(content, cursor).unwrap()
  assert_eq(first.id, "a")
  let second = @src.fasta_next_record(content, cursor).unwrap()
  assert_eq(second.id, "b")
  let third = @src.fasta_next_record(content, cursor).unwrap()
  assert_eq(third.id, "c")
  assert_true(@src.fasta_next_record(content, cursor) is None)
}

///|
test "stream fasta crlf and empty input" {
  let content = ">x\r\nAA\r\nCC\r\n"
  let streamed : Array[@src.SeqRecord] = []
  @src.fasta_each_record(content, fn(r) { streamed.push(r) })
  assert_eq(streamed.length(), 1)
  assert_eq(streamed[0].seq.to_string(), "AACC")
  @src.fasta_each_record("", fn(r) { ignore(r) })
}

///|
test "stream fasta sequence line whitespace" {
  let content = ">w\nAC GT\tTA\n"
  let streamed : Array[@src.SeqRecord] = []
  @src.fasta_each_record(content, fn(r) { streamed.push(r) })
  assert_eq(streamed[0].seq.to_string(), "ACGTTA")
}

///|
test "stream fastq matches parse_fastq" {
  let content = "@read1 first read\nACGTACGT\n+\nIIIIIIII\n@read2 desc\nGGGGCCCC\n+\nFFFFFFFF\n"
  let batch = @src.parse_fastq(content)
  let streamed : Array[@src.SeqRecord] = []
  @src.fastq_each_record(content, fn(r) { streamed.push(r) })
  assert_eq(streamed.length(), batch.length())
  assert_eq(streamed.length(), 2)
  for i in 0..<streamed.length() {
    assert_eq(streamed[i].id, batch[i].id)
    assert_eq(streamed[i].seq.to_string(), batch[i].seq.to_string())
    let q1 = streamed[i].letter_annotations.get("phred_quality").unwrap()
    let q2 = batch[i].letter_annotations.get("phred_quality").unwrap()
    assert_eq(q1, q2)
  }
  assert_eq(streamed[0].id, "read1")
  assert_eq(streamed[0].description, "read1 first read")
  let q = streamed[0].letter_annotations.get("phred_quality").unwrap()
  assert_eq(q, [40, 40, 40, 40, 40, 40, 40, 40])
}

///|
test "stream fastq cursor resume" {
  let content = "@r1\nAAAA\n+\n!!!!\n@r2\nGGGG\n+\n####\n"
  let cursor = @src.SeqioCursor::new()
  let a = @src.fastq_next_record(content, cursor).unwrap()
  assert_eq(a.id, "r1")
  let b = @src.fastq_next_record(content, cursor).unwrap()
  assert_eq(b.id, "r2")
  assert_true(@src.fastq_next_record(content, cursor) is None)
}

///|
test "stream fastq blank lines between records" {
  let content = "@r1\nAAAA\n+\n!!!!\n\n@r2\nGGGG\n+\n####\n\n\n"
  let streamed : Array[@src.SeqRecord] = []
  @src.fastq_each_record(content, fn(r) { streamed.push(r) })
  assert_eq(streamed.length(), 2)
  assert_eq(streamed[1].id, "r2")
}

///|
test "stream fastq mismatched lengths raises" {
  let content = "@bad\nAAAA\n+\n!!!\n"
  let msg = stream_fq_raises(content)
  assert_true(msg.contains("don't match"))
}

///|
fn stream_fq_raises(content : String) -> String {
  try {
    ignore(@src.fastq_next_record(content, @src.SeqioCursor::new()))
    "NO_ERROR"
  } catch {
    SeqIOError(msg) => msg
  }
}

///|
test "stream fastq missing header char raises" {
  let content = "not a header\nAAAA\n+\n!!!!\n"
  let msg = stream_fq_raises(content)
  assert_true(msg.contains("Expected '@'"))
}

///|
test "stream fastq truncated raises" {
  let content = "@r1\nAAAA\n+\n"
  let msg = stream_fq_raises(content)
  assert_true(msg.contains("Truncated"))
}

///|
test "stream fasta chunks batching" {
  let content = ">a\nAAAA\n>b\nCCCC\n>c\nTTTT\n>d\nGGGG\n>e\nACGT\n"
  let sizes : Array[Int32] = []
  let total = {
    let mut count = 0
    @src.fasta_chunks(content, 2, fn(batch) {
      sizes.push(batch.length())
      count = count + batch.length()
    })
    count
  }
  assert_eq(total, 5)
  assert_eq(sizes, [2, 2, 1])
}

///|
test "stream fastq chunks batching" {
  let content = "@r1\nAAAA\n+\n!!!!\n@r2\nCCCC\n+\n$$$$\n@r3\nTTTT\n+\n%%%%\n"
  let sizes : Array[Int32] = []
  let total = {
    let mut count = 0
    @src.fastq_chunks(content, 2, fn(batch) {
      sizes.push(batch.length())
      count = count + batch.length()
    })
    count
  }
  assert_eq(total, 3)
  assert_eq(sizes, [2, 1])
}

///|
test "stream chunks tail flush" {
  let content = ">a\nAAAA\n>b\nCCCC\n>c\nTTTT\n"
  let sizes : Array[Int32] = []
  @src.fasta_chunks(content, 10, fn(batch) { sizes.push(batch.length()) })
  assert_eq(sizes, [3])
}

///|
test "stream chunks invalid size treated as one" {
  let content = ">a\nAAAA\n>b\nCCCC\n"
  let sizes : Array[Int32] = []
  @src.fasta_chunks(content, 0, fn(batch) { sizes.push(batch.length()) })
  assert_eq(sizes, [1, 1])
}
