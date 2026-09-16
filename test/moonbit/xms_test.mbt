///|
/// Tests for xms module (Bio.motifs.xms + Bio.motifs.microrna).

test "xms_parse_basic" {
  let text = ">M1\nACGT\nTGCA\n"
  let motifs = @src.xms_parse(text)
  assert_eq(motifs.length(), 1)
  assert_eq(motifs[0].name, "M1")
  assert_eq(motifs[0].sites.length(), 2)
}

test "xms_parse_weighted" {
  let text = ">M1\nACGT 2.0\nTGCA\n"
  let motifs = @src.xms_parse(text)
  assert_true((motifs[0].sites[0].weight - 2.0).abs() < 1.0e-9)
  assert_true((motifs[0].sites[1].weight - 1.0).abs() < 1.0e-9)
}

test "xms_counts" {
  let sites = [@src.XmsSite::new("ACGT", weight=1.0), @src.XmsSite::new("ACGT", weight=1.0)]
  let motif = @src.XmsMotif::{ name: "M", sites: sites }
  let c = @src.xms_counts(motif)
  assert_eq(c.length(), 4)
  assert_eq(c[0][0], 2.0)
}

test "xms_pfm" {
  let sites = [@src.XmsSite::new("ACGT")]
  let motif = @src.XmsMotif::{ name: "M", sites: sites }
  let pfm = @src.xms_pfm(@src.xms_counts(motif))
  assert_true((pfm[0][0] - 1.0).abs() < 1.0e-9)
  assert_true((pfm[1][1] - 1.0).abs() < 1.0e-9)
}

test "xms_consensus" {
  let sites = [@src.XmsSite::new("AAAT"), @src.XmsSite::new("AACT")]
  let motif = @src.XmsMotif::{ name: "M", sites: sites }
  let pfm = @src.xms_pfm(@src.xms_counts(motif))
  let cons = @src.xms_consensus(pfm)
  assert_eq(cons, "AAAT")
}

test "xms_info_content" {
  let sites = [@src.XmsSite::new("AAAA")]
  let motif = @src.XmsMotif::{ name: "M", sites: sites }
  let pfm = @src.xms_pfm(@src.xms_counts(motif))
  let (per_pos, total) = @src.xms_info_content(pfm)
  assert_eq(per_pos.length(), 4)
  assert_true(total > 0.0)
}

test "xms_write_roundtrip" {
  let sites = [@src.XmsSite::new("ACGT", weight=2.5)]
  let m = @src.XmsMotif::{ name: "M1", sites: sites }
  let text = @src.xms_write(m)
  let parsed = @src.xms_parse(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].name, "M1")
}

test "mirs_seed_6mer" {
  let seed = @src.mirs_seed_6mer("UGAGGUAG")
  match seed {
    Some(s) => assert_eq(s, "GAGGUA")
    None => assert_true(false)
  }
}

test "mirs_classify_seeds" {
  let seeds = @src.mirs_classify_seeds("UGAGGUAG")
  assert_true(seeds.length() > 0)
}

test "mirs_reverse_complement" {
  let rc = @src.mirs_reverse_complement("ACGT")
  assert_eq(rc, "ACGT")
}

test "mirs_pfm" {
  let seqs = ["AAA", "AAC"]
  let pfm = @src.mirs_pfm(seqs)
  assert_eq(pfm.length(), 3)
  assert_true((pfm[0][0] - 1.0).abs() < 1.0e-9)
}
