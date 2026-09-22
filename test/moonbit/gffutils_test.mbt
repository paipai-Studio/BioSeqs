///|
/// Tests for the gffutils port (interval predicates, generic parent/child
/// traversal, and merging).

///|
test "gffutils_relations" {
  let content = "##gff-version 3\nchr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=gene1\nchr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1\nchr1\tsrc\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=tx1\nchr1\tsrc\texon\t300\t400\t.\t+\t.\tID=exon2;Parent=tx1\nchr1\tsrc\tgene\t600\t800\t.\t+\t.\tID=gene2\n"
  let rec = @src.bio_parse_gff(content)
  // children of gene1 == [tx1]
  let kids = @src.gffutils_children(rec, "gene1")
  assert_eq(kids.length(), 1)
  match kids[0].get_id() {
    Some(id) => assert_eq(id, "tx1")
    None => assert_true(false)
  }
  // parents of tx1 == [gene1]
  let pars = @src.gffutils_parents(rec, "tx1")
  assert_eq(pars.length(), 1)
  match pars[0].get_id() {
    Some(id) => assert_eq(id, "gene1")
    None => assert_true(false)
  }
  // all children of gene1 == tx1, exon1, exon2
  assert_eq(@src.gffutils_all_children(rec, "gene1").length(), 3)
  // all parents of exon2 == tx1, gene1
  assert_eq(@src.gffutils_all_parents(rec, "exon2").length(), 2)
}

///|
test "gffutils_intervals" {
  let gene1 = @src.GFFFeature::new("chr1", "src", "gene", 100, 500, "+")
  let gene2 = @src.GFFFeature::new("chr1", "src", "gene", 600, 800, "+")
  let exon1 = @src.GFFFeature::new("chr1", "src", "exon", 100, 200, "+")
  let exon2 = @src.GFFFeature::new("chr1", "src", "exon", 300, 400, "+")
  assert_true(@src.gffutils_overlaps(gene1, exon1))
  assert_true(not(@src.gffutils_overlaps(gene1, gene2)))
  assert_true(@src.gffutils_contains(gene1, exon1))
  assert_true(@src.gffutils_within(exon1, gene1))
  assert_true(not(@src.gffutils_contains(exon1, gene1)))
  assert_eq(@src.gffutils_distance(exon1, gene1), 0)
  assert_eq(@src.gffutils_distance(exon1, exon2), 99)
}

///|
test "gffutils_merge" {
  let f1 = @src.GFFFeature::new("chr1", "src", "exon", 100, 200, "+")
  let f2 = @src.GFFFeature::new("chr1", "src", "exon", 201, 300, "+")
  let f3 = @src.GFFFeature::new("chr1", "src", "exon", 500, 600, "+")
  let rec = @src.GFFRecord::new().add_feature(f1).add_feature(f2).add_feature(f3)
  let merged = @src.gffutils_merge(rec)
  assert_eq(merged.length(), 2)
  assert_eq((merged[0].start, merged[0].end), (100, 300))
  assert_eq((merged[1].start, merged[1].end), (500, 600))
}///|
test "gffutils_merge_full" {
  let rec = @src.bio_parse_gff(
    "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=exonA;gene=g1\n" +
    "chr1\tsrc\texon\t201\t300\t.\t+\t.\tID=exonB;gene=g1\n" +
    "chr1\tsrc\texon\t500\t600\t.\t+\t.\tID=exonC;gene=g2\n",
  )
  let merged = @src.gffutils_merge(rec)
  assert_eq(merged.length(), 2)
  assert_eq((merged[0].start, merged[0].end), (100, 300))
  assert_eq((merged[1].start, merged[1].end), (500, 600))
  // merged[0] gets a fresh ID "exon_1" and keeps the first feature's "gene"
  match merged[0].get_id() {
    Some(id) => assert_eq(id, "exon_1")
    None => assert_true(false)
  }
  match merged[0].get_attribute("gene") {
    Some(g) => assert_eq(g, "g1")
    None => assert_true(false)
  }
  // standalone merged[1] is unchanged
  match merged[1].get_id() {
    Some(id) => assert_eq(id, "exonC")
    None => assert_true(false)
  }
}

///|
test "gffutils_merge_min_dist" {
  let rec = @src.bio_parse_gff(
    "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=a\n" +
    "chr1\tsrc\texon\t205\t300\t.\t+\t.\tID=b\n",
  )
  // default min_dist=1: 205 > 201, no merge
  assert_eq(@src.gffutils_merge(rec).length(), 2)
  // min_dist=5: 205 <= 200 + 5, merge
  assert_eq(@src.gffutils_merge(rec, min_dist=5).length(), 1)
}

///|
test "gffutils_to_bed" {
  let rec = @src.bio_parse_gff("chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=gene1\n")
  let f = rec.get_features_by_id("gene1")[0]
  // GFF [100,200] 1-based -> BED [99,200) 0-based half-open
  assert_eq(@src.gffutils_to_bed(f), "chr1\t99\t200\tgene1\t.\t+")
}