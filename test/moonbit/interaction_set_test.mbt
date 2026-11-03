///|
test "anchor_creation" {
  let a = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  assert_eq(a.chrom, "chr1")
  assert_eq(a.start, 1000)
  assert_eq(a.end, 2000)
  assert_eq(@src.Anchor::width(a), 1001)
}

///|
test "anchors_overlap" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=1500, end=2500)
  let a3 = @src.Anchor::new(chrom="chr1", start=3000, end=4000)
  let a4 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  assert_true(@src.anchors_overlap(a1, a2))
  assert_false(@src.anchors_overlap(a1, a3))
  assert_false(@src.anchors_overlap(a1, a4))
}

///|
test "ginteraction_distance" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=5000, end=6000)
  let i1 = @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10)
  let i2 = @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5)
  // mid1 = 1500, mid2 = 5500, distance = 4000
  match @src.GInteraction::distance(i1) {
    Some(d) => assert_eq(d, 4000)
    None => assert_true(false)
  }
  // Inter-chromosomal has no distance
  match @src.GInteraction::distance(i2) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "ginteraction_intra_inter" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=5000, end=6000)
  let i1 = @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10)
  let i2 = @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5)
  assert_true(@src.GInteraction::is_intra_chromosomal(i1))
  assert_false(@src.GInteraction::is_inter_chromosomal(i1))
  assert_false(@src.GInteraction::is_intra_chromosomal(i2))
  assert_true(@src.GInteraction::is_inter_chromosomal(i2))
}

///|
test "interaction_set_basic" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  assert_eq(@src.InteractionSet::num_interactions(iset), 3)
  assert_eq(@src.InteractionSet::total_count(iset), 18)
}

///|
test "interaction_set_filter_by_count" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let filtered = @src.InteractionSet::filter_by_count(iset, 5)
  assert_eq(@src.InteractionSet::num_interactions(filtered), 2)
  assert_eq(@src.InteractionSet::total_count(filtered), 15)
}

///|
test "interaction_set_intra_inter_chromosomal" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let intra = @src.InteractionSet::intra_chromosomal(iset)
  let inter = @src.InteractionSet::inter_chromosomal(iset)
  assert_eq(@src.InteractionSet::num_interactions(intra), 1)
  assert_eq(@src.InteractionSet::num_interactions(inter), 2)
}

///|
test "interaction_set_interactions_on_chrom" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let on_chr1 = @src.InteractionSet::interactions_on_chrom(iset, "chr1")
  let on_chr2 = @src.InteractionSet::interactions_on_chrom(iset, "chr2")
  let on_chr3 = @src.InteractionSet::interactions_on_chrom(iset, "chr3")
  assert_eq(@src.InteractionSet::num_interactions(on_chr1), 3)
  assert_eq(@src.InteractionSet::num_interactions(on_chr2), 2)
  assert_eq(@src.InteractionSet::num_interactions(on_chr3), 0)
}

///|
test "interaction_set_interaction_matrix" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a1, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=5),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let matrix = @src.InteractionSet::interaction_matrix(iset, ["chr1", "chr2"])
  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0][0], 10)
  assert_eq(matrix[0][1], 5)
  assert_eq(matrix[1][0], 5)
  assert_eq(matrix[1][1], 0)
}

///|
test "interaction_set_distance_distribution" {
  let a1 = @src.Anchor::new(chrom="chr1", start=0, end=100)
  let a2 = @src.Anchor::new(chrom="chr1", start=1000, end=1100)
  let a3 = @src.Anchor::new(chrom="chr1", start=5000, end=5100)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let dist = @src.InteractionSet::distance_distribution(iset, 1000)
  assert_true(dist.length() > 0)
}

///|
test "interaction_set_top_interactions" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=5),
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=15),
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let top = @src.InteractionSet::top_interactions(iset, 2)
  assert_eq(top.length(), 2)
  assert_eq(top[0].count, 15)
  assert_eq(top[1].count, 10)
}

///|
test "interaction_set_inter_fraction" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let a3 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=10),
    @src.GInteraction::new(anchor1=a1, anchor2=a3, count=5),
    @src.GInteraction::new(anchor1=a2, anchor2=a3, count=3),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let frac = @src.InteractionSet::inter_fraction(iset)
  // 2 out of 3 are inter-chromosomal
  assert_true(frac > 0.66 && frac < 0.67)
}

///|
test "interaction_set_anchors" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr1", start=5000, end=6000)
  let interactions = [@src.GInteraction::new(anchor1=a1, anchor2=a2, count=10)]
  let iset = @src.InteractionSet::new(interactions~)
  let anchors = @src.InteractionSet::anchors(iset)
  assert_eq(anchors.length(), 2)
}

///|
test "interaction_set_summary" {
  let a1 = @src.Anchor::new(chrom="chr1", start=1000, end=2000)
  let a2 = @src.Anchor::new(chrom="chr2", start=1000, end=2000)
  let interactions = [
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=5),
    @src.GInteraction::new(anchor1=a1, anchor2=a2, count=15),
  ]
  let iset = @src.InteractionSet::new(interactions~)
  let s = @src.InteractionSet::summary(iset)
  assert_true(s.contains("InteractionSet"))
  assert_true(s.contains("2 interactions"))
}

///|
test "create_binned_interactions" {
  let bins = [(0, 1000), (1000, 2000), (2000, 3000)]
  let counts = [[5, 3, 0], [3, 8, 2], [0, 2, 6]]
  let iset = @src.create_binned_interactions("chr1", bins, counts)
  // Should have 5 interactions (where count > 0): (0,0)=5, (0,1)=3, (1,1)=8, (1,2)=2, (2,2)=6
  assert_eq(@src.InteractionSet::num_interactions(iset), 5)
  assert_eq(@src.InteractionSet::total_count(iset), 24)
}
