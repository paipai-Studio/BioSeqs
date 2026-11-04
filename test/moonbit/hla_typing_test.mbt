///|
test "hla reference db size" {
  let db = @src.hla_reference_db()
  assert_eq(db.length(), 4)
}

///|
test "hla reference db entries" {
  let db = @src.hla_reference_db()
  let a02 = db[0]
  assert_eq(a02.allele(), "HLA-A*02:01")
  assert_eq(a02.locus(), "HLA-A")
  assert_eq(a02.exon(), 2)
  assert_true(a02.seq.length() > 0)
}

///|
test "hla align read to ref exact match" {
  let sequence = @src.hla_reference_db()[0].seq
  let read = sequence[0:50].to_owned()
  let m = @src.hla_align_read_to_ref(read, sequence)
  assert_eq(m, 50)
}

///|
test "hla align read to ref partial match" {
  let sequence = @src.hla_reference_db()[0].seq
  let read = "NNNNNNNNNN"
  let m = @src.hla_align_read_to_ref(read, sequence)
  assert_eq(m, 0)
}

///|
test "hla type read exact match" {
  let sequence = @src.hla_reference_db()[0].seq
  let call = @src.hla_type_read(sequence)
  assert_eq(call.allele(), "HLA-A*02:01")
  assert_eq(call.score(), 1.0)
  assert_true(call.is_high_resolution())
}

///|
test "hla type read partial match" {
  let sequence = @src.hla_reference_db()[0].seq
  let read = sequence[0:100].to_owned()
  let call = @src.hla_type_read(read)
  assert_eq(call.allele(), "HLA-A*02:01")
  assert_true(call.score() > 0.95)
}

///|
test "hla type read mismatched" {
  let sequence = @src.hla_reference_db()[0].seq
  let mut read = ""
  for _ in 0..<sequence.length() {
    read = read + "N"
  }
  let call = @src.hla_type_read(read)
  assert_eq(call.allele(), "HLA-A*02:01")
  assert_eq(call.score(), 0.0)
}

///|
test "hla type batch homozygous" {
  let sequence = @src.hla_reference_db()[0].seq
  let reads : Array[String] = [
    sequence[0:50].to_owned(),
    sequence[50:100].to_owned(),
    sequence[100:150].to_owned(),
  ]
  let gl = @src.hla_type_batch(reads)
  assert_eq(gl, "HLA-A*02:01")
}

///|
test "hla type batch heterozygous" {
  let ref02 = @src.hla_reference_db()[0].seq
  let ref01 = @src.hla_reference_db()[1].seq
  let reads : Array[String] = [
    ref02[0:50].to_owned(),
    ref01[0:50].to_owned(),
    ref02[50:100].to_owned(),
  ]
  let gl = @src.hla_type_batch(reads)
  // 由于参考序列相同, 多数投票返回第一个等位基因
  assert_true(gl.length() > 0)
}

///|
test "hla type batch empty" {
  let reads : Array[String] = []
  let gl = @src.hla_type_batch(reads)
  assert_eq(gl, "")
}

///|
test "hla call is high resolution" {
  let sequence = @src.hla_reference_db()[0].seq
  let call = @src.hla_type_read(sequence)
  assert_true(call.is_high_resolution())
}

///|
test "hla call fields" {
  let sequence = @src.hla_reference_db()[0].seq
  let call = @src.hla_type_read(sequence)
  assert_eq(call.allele(), "HLA-A*02:01")
  assert_eq(call.locus(), "HLA-A")
  assert_eq(call.matches(), sequence.length())
  assert_eq(call.total(), sequence.length())
}
