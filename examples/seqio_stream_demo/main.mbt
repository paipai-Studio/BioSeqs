///|
/// 流式(记录级惰性)解析演示: 逐条与分批处理大 FASTA/FASTQ.
fn main raise {
  let fasta = ">seq1 kinase-like\nMKTLLLTLVV\nMKTAYIAKQR\n>seq2\nGGGGTTTTCCCCAA\n>seq3\nACGTACGTAC\n"
  println("=== 流式 FASTA: 逐条 ===")
  @src.fasta_each_record(fasta, fn(rec) {
    println("  " + rec.id + ": len=" + rec.seq.length().to_string())
  })
  println("=== 流式 FASTA: 游标暂停/恢复 ===")
  let cursor = @src.SeqioCursor::new()
  let first = @src.fasta_next_record(fasta, cursor).unwrap()
  println("  first: " + first.id)
  let second = @src.fasta_next_record(fasta, cursor).unwrap()
  println("  second: " + second.id + " seq=" + second.seq.to_string())
  println("=== 流式 FASTA: 分批(每批2条) ===")
  @src.fasta_chunks(fasta, 2, fn(batch) {
    let ids = batch.map(fn(r) { r.id })
    println(
      "  batch(size=" + batch.length().to_string() + "): " + ids.join(", "),
    )
  })
  let fastq = "@r1 first\nACGTACGT\n+\nIIIIIIII\n@r2\nGGGGCCCC\n+\nFFFFFFFF\n@r3\nTTTTAAAA\n+\n!!!!!!!!\n"
  println("=== 流式 FASTQ: 逐条(含质量统计) ===")
  @src.fastq_each_record(fastq, fn(rec) {
    let quals = rec.letter_annotations.get("phred_quality").unwrap()
    let mut qsum = 0
    for q in quals {
      qsum = qsum + q
    }
    let mean_q = qsum.to_double() / quals.length().to_double()
    println(
      "  " +
      rec.id +
      ": len=" +
      rec.seq.length().to_string() +
      " meanQ=" +
      mean_q.to_string(),
    )
  })
  println("=== 流式 FASTQ: 分批(每批2条) ===")
  @src.fastq_chunks(fastq, 2, fn(batch) {
    println("  batch size=" + batch.length().to_string())
    for rec in batch {
      println("    - " + rec.id)
    }
  })
}
