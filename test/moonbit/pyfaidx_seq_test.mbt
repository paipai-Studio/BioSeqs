///|
/// Tests for the pyfaidx sequence-utility port. Expected values come from
/// pyfaidx 0.9.0.4 (`complement`, `Sequence.reverse`, sequence GC properties,
/// `wrap_sequence`, `bed_split`, `ucsc_split`).

///|
fn pyf_close(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-9
}

///|
test "pyfaidx_complement_reverse" {
  assert_eq(@src.pyfaidx_complement("ATCGTA"), "TAGCAT")
  assert_eq(@src.pyfaidx_complement("gattaca"), "ctaatgt")
  assert_eq(@src.pyfaidx_complement("NMRSWYKVHDBX"), "NKYSWRMBDHVX")
  assert_eq(@src.pyfaidx_reverse("ATCGTA"), "ATGCTA")
  assert_eq(@src.pyfaidx_reverse_complement("ATCGTA"), "TACGAT")
}

///|
test "pyfaidx_complement_invalid_raises" {
  let raised = try {
    ignore(@src.pyfaidx_complement("A1G"))
    false
  } catch {
    PyfaidxError(_) => true
  }
  assert_true(raised)
}

///|
test "pyfaidx_gc_content" {
  assert_true(pyf_close(@src.pyfaidx_gc("ATCGTA"), 0.3333333333))
  assert_true(pyf_close(@src.pyfaidx_gc("NMRATCGTA"), 0.2222222222))
  assert_true(pyf_close(@src.pyfaidx_gc_strict("NMRATCGTA"), 0.3333333333))
  assert_true(pyf_close(@src.pyfaidx_gc_iupac("NMRATCGTA"), 0.3611111111))
}

///|
test "pyfaidx_wrap_sequence" {
  assert_eq(@src.pyfaidx_wrap_sequence(3, "ACGTACGTAC"), "ACG\nTAC\nGTA\nC\n")
  assert_eq(@src.pyfaidx_wrap_sequence(4, "ACGTAC"), "ACGT\nAC\n")
}

///|
test "pyfaidx_region_parsing" {
  let b1 = @src.pyfaidx_bed_split("chr1\t10\t20")
  match (b1.0, b1.1, b1.2) {
    (Some(r), Some(s), Some(e)) => {
      assert_eq(r, "chr1")
      assert_eq(s, 10)
      assert_eq(e, 20)
    }
    _ => abort("bed_split should parse chr1:10-20")
  }
  let b2 = @src.pyfaidx_bed_split("#foo")
  assert_true(b2.0.is_none())
  assert_true(b2.1.is_none())
  assert_true(b2.2.is_none())

  let u1 = @src.pyfaidx_ucsc_split("chr1:100-200")
  match (u1.1, u1.2) {
    (Some(s), Some(e)) => {
      assert_eq(s, 99)
      assert_eq(e, 200)
    }
    _ => abort("ucsc_split should parse chr1:100-200")
  }
  assert_eq(u1.0, "chr1")

  let u2 = @src.pyfaidx_ucsc_split("chr1")
  assert_eq(u2.0, "chr1")
  assert_true(u2.1.is_none())
  assert_true(u2.2.is_none())
}