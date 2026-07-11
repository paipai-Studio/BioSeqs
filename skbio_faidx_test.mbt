///|
/// Test data: a simple FASTA file for testing.
let test_fasta_content = ">chr1 description of chr1\nACGTACGTACGTACGTACGTACGTACGTGC\nATGCATGCATGCATGCATGCATGCATGCAT\nACGTACGTACGTACGTACGTACGTACGTGC\n>chr2\nATGCATGCATGCAT\nGCATGCATGCATGC\n>chr3 short\nAAAA\n"

///|
/// Test data: expected .fai index content (LF line endings)
let expected_fai_content = "chr1\t90\t26\t30\t31\nchr2\t28\t125\t14\t15\nchr3\t4\t167\t4\t5\n"

///|
test "faidx_build" {
  let idx = FastaIndex::build(test_fasta_content)
  
  assert_eq(idx.len(), 3)
  assert_true(idx.contains("chr1"))
  assert_true(idx.contains("chr2"))
  assert_true(idx.contains("chr3"))
  assert_false(idx.contains("nonexistent"))
  
  let keys = idx.keys()
  assert_eq(keys.length(), 3)
  
  let fai_content = idx.to_fai_string()
  assert_eq(fai_content, expected_fai_content)
}

///|
test "faidx_get_entry" {
  let idx = FastaIndex::build(test_fasta_content)
  
  match idx.get_entry("chr1") {
    None => assert_true(false)
    Some(entry) => {
      assert_eq(entry.name, "chr1")
      assert_eq(entry.length, 90)
      assert_eq(entry.offset, 26)
      assert_eq(entry.line_bases, 30)
      assert_eq(entry.line_width, 31)
    }
  }
  
  match idx.get_entry("chr2") {
    None => assert_true(false)
    Some(entry) => {
      assert_eq(entry.name, "chr2")
      assert_eq(entry.length, 28)
      assert_eq(entry.offset, 125)
      assert_eq(entry.line_bases, 14)
      assert_eq(entry.line_width, 15)
    }
  }
  
  match idx.get_entry("chr3") {
    None => assert_true(false)
    Some(entry) => {
      assert_eq(entry.name, "chr3")
      assert_eq(entry.length, 4)
      assert_eq(entry.offset, 167)
      assert_eq(entry.line_bases, 4)
      assert_eq(entry.line_width, 5)
    }
  }
  
  match idx.get_entry("nonexistent") {
    None => assert_true(true)
    Some(_) => assert_true(false)
  }
}

///|
test "faidx_get_length" {
  let idx = FastaIndex::build(test_fasta_content)
  
  assert_eq(idx.get_length("chr1"), Some(90))
  assert_eq(idx.get_length("chr2"), Some(28))
  assert_eq(idx.get_length("chr3"), Some(4))
  assert_eq(idx.get_length("nonexistent"), None)
}

///|
test "faidx_fetch_full" {
  let idx = FastaIndex::build(test_fasta_content)
  
  match idx.fetch_full("chr1") {
    None => assert_true(false)
    Some(seq) => {
      assert_eq(seq.length(), 90)
    }
  }
  
  let seq2 = idx.fetch_full("chr2")
  assert_eq(seq2, Some("ATGCATGCATGCATGCATGCATGCATGC"))
  
  let seq3 = idx.fetch_full("chr3")
  assert_eq(seq3, Some("AAAA"))
  
  assert_eq(idx.fetch_full("nonexistent"), None)
}

///|
test "faidx_fetch_subseq" {
  let idx = FastaIndex::build(test_fasta_content)
  
  assert_eq(idx.fetch("chr1", 1, 10), Some("ACGTACGTAC"))
  assert_eq(idx.fetch("chr1", 1, 30), Some("ACGTACGTACGTACGTACGTACGTACGTGC"))
  assert_eq(idx.fetch("chr1", 31, 59), Some("ATGCATGCATGCATGCATGCATGCATGCA"))
  assert_eq(idx.fetch("chr1", 60, 67), Some("TACGTACG"))
  assert_eq(idx.fetch("chr1", 20, 25), Some("TACGTA"))
  
  assert_eq(idx.fetch("chr2", 1, 14), Some("ATGCATGCATGCAT"))
  assert_eq(idx.fetch("chr2", 15, 28), Some("GCATGCATGCATGC"))
  
  assert_eq(idx.fetch("chr3", 1, 4), Some("AAAA"))
  
  assert_eq(idx.fetch("chr1", 100, 200), Some(""))
  assert_eq(idx.fetch("chr1", 0, 5), Some("ACGTA"))
  assert_eq(idx.fetch("chr1", 88, 92), Some("TGC"))
  
  assert_eq(idx.fetch("nonexistent", 1, 10), None)
}

///|
test "faidx_complement" {
  let idx = FastaIndex::build(">seq\nATGC\n")
  
  let seq = idx.fetch("seq", 1, 4)
  assert_eq(seq, Some("ATGC"))
}

///|
test "faidx_spaces_in_sequence" {
  let idx = FastaIndex::build(">seq\nA T G C\n")
  
  match idx.get_entry("seq") {
    None => assert_true(false)
    Some(entry) => assert_eq(entry.length, 4)
  }
  
  assert_eq(idx.fetch("seq", 1, 4), Some("ATGC"))
}

///|
test "faidx_single_line_sequence" {
  let idx = FastaIndex::build(">single\nACGTACGT\n")
  
  match idx.get_entry("single") {
    None => assert_true(false)
    Some(entry) => {
      assert_eq(entry.length, 8)
      assert_eq(entry.line_bases, 8)
      assert_eq(entry.line_width, 9)
    }
  }
  
  assert_eq(idx.fetch("single", 1, 8), Some("ACGTACGT"))
  assert_eq(idx.fetch("single", 2, 5), Some("CGTA"))
}

///|
test "faidx_empty_sequence" {
  let idx = FastaIndex::build(">empty\n\n>next\nATGC\n")
  
  match idx.get_entry("empty") {
    None => assert_true(false)
    Some(entry) => assert_eq(entry.length, 0)
  }
  
  assert_eq(idx.fetch("empty", 1, 10), Some(""))
  assert_eq(idx.fetch("next", 1, 4), Some("ATGC"))
}
