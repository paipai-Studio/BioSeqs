///|
/// Tests for seq_cluster.mbt and @src.seq_similarity_map.mbt.

///|
fn seq_cl_record(s : String) -> @src.SeqRecord {
  @src.SeqRecord::new(@src.Seq::new(s), id="c")
}

///|
test "seq_kmer_tanimoto_identical_and_distinct" {
  let a = @src.Seq::new("ACGTACGTACGT")
  assert_eq(@src.seq_kmer_tanimoto(a, a), 1.0)
  // unrelated low-similarity sequence
  let b = @src.Seq::new("TTTTTTTTTTTT")
  assert_true(@src.seq_kmer_tanimoto(a, b) < 0.2)
}

///|
test "seq_butina_groups_identical_sequences" {
  let records = [
    seq_cl_record("ACGTACGTACGT"),
    seq_cl_record("ACGTACGTACGT"),
    seq_cl_record("ACGTACGTACGT"),
    seq_cl_record("TTTTTTTTTTTT"),
  ]
  let clusters = @src.seq_butina_cluster(records, k=3, threshold=0.3)
  assert_true(clusters.length() >= 1)
  // all copies of the same sequence must share one cluster
  let mut home = -1
  for ci = 0; ci < clusters.length(); ci = ci + 1 {
    for m in clusters[ci] {
      if m == 0 {
        home = ci
      }
    }
  }
  for m in clusters[home] {
    assert_true(m == 0 || m == 1 || m == 2)
  }
  assert_eq(@src.seq_butina_largest_cluster(clusters) >= 3, true)
}

///|
test "seq_similarity_map_length_and_content" {
  // query embeds the reference motif -> window over it scores high
  let reference = @src.Seq::new("ACGTACGTACGT")
  let query = @src.Seq::new("TTTTACGTACGTAAAA")
  let map = @src.seq_similarity_map(query, reference, k=2, window=4)
  assert_eq(map.length(), 13)
  // window starting at position 4 ("ACGT") matches perfectly
  let mut best = 0.0
  for v in map {
    if v > best {
      best = v
    }
  }
  // window "ACGT" 2-mers {AC,CG,GT} vs reference set {AC,CG,GT,TA}:
  // Jaccard = 3/4
  assert_true(best > 0.74 && best < 0.76)
}

///|
test "seq_similarity_peaks_finds_local_maxima" {
  // synthetic profile with one clear peak
  let map = [0.1, 0.2, 0.8, 0.3, 0.1]
  let peaks = @src.seq_similarity_peaks(map, threshold=0.5)
  assert_eq(peaks, [2])
  // below-threshold profiles yield no peaks
  assert_eq(@src.seq_similarity_peaks([0.1, 0.2], threshold=0.5).length(), 0)
}

///|
test "seq_similarity_map_ascii_width" {
  let map = [0.0, 0.25, 0.5, 0.75, 1.0]
  let art = @src.seq_similarity_map_ascii(map, width=5)
  assert_eq(art.length(), 5)
  // weakest bucket is a blank, strongest is '#'
  assert_eq(art[0], ' ')
  assert_eq(art[4], '#')
}
