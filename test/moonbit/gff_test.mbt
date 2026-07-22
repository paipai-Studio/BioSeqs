///|
/// GFF module tests

test "GFFFeature::new" {
  let feature = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
  assert_eq(feature.seqid, "chr1")
  assert_eq(feature.feature_type, "gene")
  assert_eq(feature.start, 1000)
  assert_eq(feature.end, 2000)
}

test "GFFFeature::add_attribute" {
  let feature = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
    .add_attribute("ID", "gene001")
    .add_attribute("Name", "TP53")
  match feature.get_attribute("ID") {
    Some(id) => assert_eq(id, "gene001")
    None => assert_true(false)
  }
}

test "GFFFeature::get_id" {
  let feature = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
    .add_attribute("ID", "gene001")
  match feature.get_id() {
    Some(id) => assert_eq(id, "gene001")
    None => assert_true(false)
  }
}

test "GFFFeature::length" {
  let feature = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
  assert_eq(feature.length(), 1001)
}

test "GFFFeature::is_coding" {
  let gene = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
  assert_true(gene.is_coding())
  let cds = @src.GFFFeature::new("chr1", "Ensembl", "CDS", 1000, 2000, "+")
  assert_true(cds.is_coding())
}

test "GFFRecord::new" {
  let record = @src.GFFRecord::new()
  assert_eq(record.count_features(), 0)
}

test "GFFRecord::add_feature" {
  let record = @src.GFFRecord::new()
  let feature = @src.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
  ignore(record.add_feature(feature))
  assert_eq(record.count_features(), 1)
}

test "GFFRecord::get_features_by_type" {
  let record = @src.create_example_gff()
  let genes = record.get_genes()
  assert_true(genes.length() >= 1)
}

test "GFFRecord::get_features_by_seqid" {
  let record = @src.create_example_gff()
  let features = record.get_features_by_seqid("chr1")
  assert_true(features.length() > 0)
}

test "parse_gff" {
  let content = "##gff-version 3\nchr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene001;Name=TP53\n"
  let record = @src.bio_parse_gff(content)
  assert_eq(record.count_features(), 1)
  assert_eq(record.version, "3")
}

test "parse_attributes" {
  let attrs = @src.parse_attributes("ID=gene001;Name=TP53;biotype=protein_coding")
  match attrs.get("ID") {
    Some(id) => assert_eq(id, "gene001")
    None => assert_true(false)
  }
  match attrs.get("Name") {
    Some(name) => assert_eq(name, "TP53")
    None => assert_true(false)
  }
}

test "parse_attributes_unescape" {
  let attrs = @src.parse_attributes("ID=gene001;Name=value%3Bwith%3Bsemicolons")
  match attrs.get("Name") {
    Some(name) => assert_eq(name, "value;with;semicolons")
    None => assert_true(false)
  }
}

test "create_example_gff" {
  let record = @src.create_example_gff()
  assert_true(record.count_features() > 5)
}

test "GFFRecord::get_child_features" {
  let record = @src.create_example_gff()
  let mrnas = record.get_child_features("gene:ENSG00000130203")
  assert_true(mrnas.length() >= 1)
}

test "GFFRecord::get_unique_seqids" {
  let record = @src.create_example_gff()
  let seqids = record.get_unique_seqids()
  assert_true(seqids.length() >= 1)
}

test "GFFRecord::get_features_in_range" {
  let record = @src.create_example_gff()
  let features = record.get_features_in_range("chr1", 10000, 11000)
  assert_true(features.length() > 0)
}