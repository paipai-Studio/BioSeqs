///|
/// Tests for orfik module (ORFik-style ORF detection).

test "orff_kozak_score_perfect" {
  // gccRccAUGG: -3=A, +4=G, -2=C, -1=C
  let seq = "GCCGCCAUGG"
  let score = @src.orff_kozak_score(seq, 6)  // ATG starts at position 6
  assert_true(score > 0.5)
}

test "orff_kozak_score_weak" {
  // No Kozak features: -2 and -1 are T (pyrimidine -> +0.15 each)
  let seq = "TTTTTTATGT"
  let score = @src.orff_kozak_score(seq, 6)
  assert_eq(score, 0.3)
}

test "orff_find_all_orfs_basic" {
  let seq = "ATGAAATAAATGCCCTAGATGGGGGTGA"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  assert_true(orfs.length() >= 2)
  // First ORF: ATG AAA TAA (9 nt) -> M, K
  assert_eq(orfs[0].protein, "MK")
}

test "orff_find_all_orfs_min_length" {
  let seq = "ATGTAATAA"  // Only 3 nt ORF
  let orfs = @src.orff_find_all_orfs(seq, min_length=30)
  assert_eq(orfs.length(), 0)
}

test "orff_translate_basic" {
  let seq = "ATGAAATAA"
  let protein = @src.orff_translate(seq)
  assert_eq(protein, "MK")
}

test "orff_classify_orf_cds" {
  let orf = @src.OrffORF::{
    start: 10, stop: 40, frame: @src.OrffFrameF0,
    length: 30, sequence: "", protein: "",
    kozak_score: 0.0, orf_class: @src.OrffIntergenic,
  }
  let cls = @src.orff_classify_orf(orf, 10, 40)
  assert_eq(@src.orff_class_label(cls), "CDS")
}

test "orff_classify_orf_uorf" {
  let orf = @src.OrffORF::{
    start: 5, stop: 20, frame: @src.OrffFrameF0,
    length: 15, sequence: "", protein: "",
    kozak_score: 0.0, orf_class: @src.OrffIntergenic,
  }
  let cls = @src.orff_classify_orf(orf, 25, 55)
  assert_eq(@src.orff_class_label(cls), "uORF")
}

test "orff_classify_orf_dorf" {
  let orf = @src.OrffORF::{
    start: 50, stop: 70, frame: @src.OrffFrameF0,
    length: 20, sequence: "", protein: "",
    kozak_score: 0.0, orf_class: @src.OrffIntergenic,
  }
  let cls = @src.orff_classify_orf(orf, 10, 40)
  assert_eq(@src.orff_class_label(cls), "dORF")
}

test "orff_classify_all" {
  let seq = "ATGAAATAAATGCCCTAGATGGGGGTGA"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  let classified = @src.orff_classify_all(orfs, 9, 18)
  assert_true(classified.length() >= 2)
  // Find the CDS (start=9, stop=18)
  let mut found_cds = false
  for orf in classified {
    if orf.start == 9 && orf.stop == 18 {
      assert_eq(@src.orff_class_label(orf.orf_class), "CDS")
      found_cds = true
    }
  }
  assert_true(found_cds)
}

test "orff_find_uorfs" {
  // CDS starts at position 30
  let seq = "ATGAAATAAATGCCCTAGATGGGGGTGAAATAAATGCCCTAGATGGGGGTGA"
  let uorfs = @src.orff_find_uorfs(seq, 30, min_length=9)
  assert_true(uorfs.length() > 0)
  for orf in uorfs {
    assert_eq(@src.orff_class_label(orf.orf_class), "uORF")
  }
}

test "orff_find_dorfs" {
  // CDS ends at position 30
  let seq = "ATGAAATAAATGCCCTAGATGGGGGTGAAATAAATGCCCTAGATGGGGGTGA"
  let dorfs = @src.orff_find_dorfs(seq, 30, min_length=9)
  for orf in dorfs {
    assert_eq(@src.orff_class_label(orf.orf_class), "dORF")
  }
}

test "orff_filter_by_kozak" {
  let seq = "ATGAAATAAATGCCCTAG"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  let filtered = @src.orff_filter_by_kozak(orfs, 0.0)
  assert_eq(filtered.length(), orfs.length())
}

test "orff_filter_by_aa_length" {
  let seq = "ATGAAATAAATGCCCTAG"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  let filtered = @src.orff_filter_by_aa_length(orfs, 1)
  assert_true(filtered.length() > 0)
}

test "orff_count_by_frame" {
  let seq = "ATGAAATAAATGCCCTAG"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  let counts = @src.orff_count_by_frame(orfs)
  assert_true(counts.get("F0").unwrap_or(0) > 0)
}

test "orff_export_orfs" {
  let seq = "ATGAAATAA"
  let orfs = @src.orff_find_all_orfs(seq, min_length=9)
  let gff = @src.orff_export_orfs(orfs)
  assert_true(gff.contains("start"))
  assert_true(gff.contains("K"))
}

test "orff_overlapping_orfs" {
  let orfs = [
    @src.OrffORF::{start: 10, stop: 40, frame: @src.OrffFrameF0, length: 30, sequence: "", protein: "", kozak_score: 0.0, orf_class: @src.OrffIntergenic},
    @src.OrffORF::{start: 35, stop: 60, frame: @src.OrffFrameF1, length: 25, sequence: "", protein: "", kozak_score: 0.0, orf_class: @src.OrffIntergenic},
  ]
  let pairs = @src.orff_overlapping_orfs(orfs)
  assert_eq(pairs.length(), 1)
}
