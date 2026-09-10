// Whitebox tests for the format registry: alias normalization, alias
// uniqueness, and capability matrices (spec 5.1.1 rule 5, spec 6.1).

///|
test "format_registry_alias_matches_canonical_parse" {
  // Every alias must normalize to its canonical name, and (once the entry
  // points dispatch through the registry) aliases must behave exactly like
  // their canonical names (spec 5.1.1 rule 5).
  let alias_cases : Array[(String, String, String)] = [
    // (content, alias, canonical)
    (">s1\nACGT", "fasta", "fasta"),
    ("@s1\nACGT\n+\nIIII", "fastq", "fastq-sanger"),
    ("@s1\nACGT\n+\n@ABC\n", "fastq-illumina", "fastq-illumina"),
    ("@s1\nACGT\n+\n;>?@\n", "fastq-solexa", "fastq-solexa"),
    (">s1\n10 20 30\n", "qual", "qual"),
    ("s1\tACGT", "tab", "tab"),
  ]
  for case in alias_cases {
    let (content, alias, canonical) = case
    let via_alias = try {
      @src.seqio_normalize_format(alias)
    } catch {
      SeqIOError(msg) => fail("normalize alias '" + alias + "' raised: " + msg)
    }
    let via_canonical = try {
      @src.seqio_normalize_format(canonical)
    } catch {
      SeqIOError(msg) =>
        fail("normalize canonical '" + canonical + "' raised: " + msg)
    }
    assert_eq(via_alias, via_canonical)
    // Both paths must also produce identical results through the entry point
    let a = try {
      @src.seqio_parse(content, alias)
    } catch {
      SeqIOError(msg) => fail("parse alias '" + alias + "' raised: " + msg)
    }
    let b = try {
      @src.seqio_parse(content, canonical)
    } catch {
      SeqIOError(msg) =>
        fail("parse canonical '" + canonical + "' raised: " + msg)
    }
    assert_eq(a.length(), b.length())
    assert_eq(a[0].id, b[0].id)
    assert_eq(a[0].seq.to_string(), b[0].seq.to_string())
  }
}

///|
test "format_registry_normalize_canonical_names" {
  let expected : Array[(String, String)] = [
    ("fasta", "fasta"),
    ("fastq", "fastq-sanger"),
    ("fastq-sanger", "fastq-sanger"),
    ("fastq-illumina", "fastq-illumina"),
    ("fastq-solexa", "fastq-solexa"),
    ("qual", "qual"),
    ("genbank", "genbank"),
    ("gb", "genbank"),
    ("embl", "embl"),
    ("pir", "pir"),
    ("nbrf", "pir"),
    ("tab", "tab"),
    ("swiss", "swiss"),
    ("uniprot-xml", "uniprot-xml"),
    ("seqxml", "seqxml"),
    ("imgt", "imgt"),
    ("ig", "ig"),
    ("insdc", "insdc"),
    ("phd", "phd"),
    ("ace", "ace"),
    ("gfa", "gfa"),
    ("pdb-seqres", "pdb-seqres"),
    ("pdb-atom", "pdb-atom"),
    ("2bit", "2bit"),
    ("sff", "sff"),
    ("sff-trim", "sff-trim"),
    ("snapgene", "snapgene"),
    ("gck", "gck"),
    ("xdna", "xdna"),
    ("nib", "nib"),
    ("fasta-2line", "fasta-2line"),
    ("cif-seqres", "cif-seqres"),
    ("cif-atom", "cif-atom"),
  ]
  for case in expected {
    let (input, canonical) = case
    assert_eq(@src.seqio_normalize_format(input), canonical)
  }
}

///|
test "format_registry_normalize_case_insensitive" {
  assert_eq(@src.seqio_normalize_format("GB"), "genbank")
  assert_eq(@src.seqio_normalize_format("NBRF"), "pir")
  assert_eq(@src.seqio_normalize_format("Fasta"), "fasta")
  assert_eq(@src.seqio_normalize_format("Sff-Trim"), "sff-trim")
}

///|
test "format_registry_normalize_unknown_raises" {
  let raised = try {
    ignore(@src.seqio_normalize_format("not_a_format"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: not_a_format"
  }
  assert_true(raised)
  // The original casing must be echoed back
  let raised2 = try {
    ignore(@src.seqio_normalize_format("Not_A_Format"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown format: Not_A_Format"
  }
  assert_true(raised2)
}

///|
test "format_registry_alignio_normalize" {
  let expected : Array[(String, String)] = [
    ("clustal", "clustal"),
    ("phylip", "phylip"),
    ("phylip-sequential", "phylip-sequential"),
    ("fasta", "fasta"),
    ("stockholm", "stockholm"),
    ("sto", "stockholm"),
    ("nexus", "nexus"),
    ("maf", "maf"),
    ("msf", "msf"),
    ("emboss", "emboss"),
    ("a2m", "a2m"),
    ("sam", "sam"),
    ("psl", "psl"),
    ("mauve", "mauve"),
  ]
  for case in expected {
    let (input, canonical) = case
    assert_eq(@src.alignio_normalize_format(input), canonical)
  }
  // Case-insensitive
  assert_eq(@src.alignio_normalize_format("STO"), "stockholm")
  // Unknown
  let raised = try {
    ignore(@src.alignio_normalize_format("not_an_align_format"))
    false
  } catch {
    AlignIOError(msg) => msg == "Unknown format: not_an_align_format"
  }
  assert_true(raised)
}

///|
test "format_registry_capability_matrices" {
  let supported = @src.seqio_supported_formats()
  let writable = @src.seqio_writable_formats()
  let indexable = @src.seqio_indexable_formats()
  let binary = @src.seqio_binary_formats()
  // Parse supports the largest set; write/index must be subsets
  assert_true(supported.length() >= 30)
  for fmt in writable {
    assert_true(supported.contains(fmt))
  }
  for fmt in indexable {
    assert_true(supported.contains(fmt))
  }
  for fmt in binary {
    assert_true(supported.contains(fmt))
  }
  // Binary formats must never claim writing support
  for fmt in binary {
    assert_true(!writable.contains(fmt))
  }
  // Known members
  assert_true(writable.contains("fasta"))
  assert_true(writable.contains("fasta-2line"))
  assert_true(writable.contains("embl"))
  assert_true(writable.contains("pir"))
  assert_true(indexable.contains("qual"))
  assert_true(indexable.contains("swiss"))
  assert_true(indexable.contains("seqxml"))
  assert_true(binary.contains("2bit"))
  assert_true(binary.contains("sff"))
}

///|
test "format_registry_alias_uniqueness" {
  // An alias must never resolve to two different canonical names.
  // Normalizing each registered identifier must be idempotent: the result
  // is itself a registered canonical name mapping to itself.
  let all : Array[String] = []
  for fmt in @src.seqio_supported_formats() {
    all.push(fmt)
  }
  for fmt in all {
    let normalized = @src.seqio_normalize_format(fmt)
    assert_eq(normalized, fmt)
  }
}