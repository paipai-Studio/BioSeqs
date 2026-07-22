///|
/// GFF example: GFF format parsing and gene annotation operations.
///
/// Run: moon run examples/gff_demo/main.mbt
fn main {
  println("=== Bio.GFF Demo ===\n")

  let record = @bio.create_example_gff()

  println("1. Basic Statistics:")
  println("   Total features: \{record.count_features()}")
  println("   GFF version: \{record.version}")

  println("\n2. Unique sequence IDs:")
  let seqids = record.get_unique_seqids()
  for i = 0; i < seqids.length(); i = i + 1 {
    println("   - \{seqids[i]}")
  }

  println("\n3. Unique feature types:")
  let ftypes = record.get_unique_feature_types()
  for i = 0; i < ftypes.length(); i = i + 1 {
    println("   - \{ftypes[i]}")
  }

  println("\n4. Genes in record:")
  let genes = record.get_genes()
  for i = 0; i < genes.length(); i = i + 1 {
    let gene = genes[i]
    match gene.get_name() {
      Some(name) => println("   - \{name} (chr\{gene.seqid}, \{gene.start}-\{gene.end})")
      None => println("   - gene (chr\{gene.seqid}, \{gene.start}-\{gene.end})")
    }
  }

  println("\n5. mRNAs in record:")
  let mrnas = record.get_mrnas()
  for i = 0; i < mrnas.length(); i = i + 1 {
    let mrna = mrnas[i]
    match mrna.get_id() {
      Some(id) => println("   - \{id}")
      None => println("   - unnamed mRNA")
    }
  }

  println("\n6. CDS features:")
  let cds = record.get_cds_features()
  for i = 0; i < cds.length(); i = i + 1 {
    let c = cds[i]
    println("   - CDS \{i + 1}: \{c.start}-\{c.end} (length: \{c.length()})")
  }

  println("\n7. Exons:")
  let exons = record.get_exons()
  for i = 0; i < exons.length(); i = i + 1 {
    let e = exons[i]
    println("   - Exon \{i + 1}: \{e.start}-\{e.end} (length: \{e.length()})")
  }

  println("\n8. Child features of gene:ENSG00000130203")
  let children = record.get_child_features("gene:ENSG00000130203")
  for i = 0; i < children.length(); i = i + 1 {
    let child = children[i]
    println("   - \{child.feature_type}: \{child.start}-\{child.end}")
  }

  println("\n9. Parse GFF string:")
  let gff_content = "##gff-version 3\nchr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene001;Name=TestGene\nchr1\tEnsembl\tmRNA\t1100\t1900\t.\t+\t.\tID=mrna001;Parent=gene001\nchr1\tEnsembl\tCDS\t1100\t1200\t.\t+\t0\tID=cds001;Parent=mrna001\n"
  let parsed = @bio.bio_parse_gff(gff_content)
  println("   Parsed \{parsed.count_features()} features")

  println("\n10. Feature in range:")
  let range_features = record.get_features_in_range("chr1", 10000, 11000)
  println("   Features in chr1:10000-11000: \{range_features.length()}")

  println("\n11. GFFFeature operations:")
  let feature = @bio.GFFFeature::new("chr1", "Ensembl", "gene", 1000, 2000, "+")
    .add_attribute("ID", "test_gene")
    .add_attribute("Name", "Test")
  match feature.get_id() {
    Some(id) => println("   Feature ID: \{id}")
    None => println("   No ID")
  }
  match feature.get_name() {
    Some(name) => println("   Feature Name: \{name}")
    None => println("   No Name")
  }
  println("   Length: \{feature.length()}")
  println("   Is coding: \{feature.is_coding()}")
  println("   Forward strand: \{feature.is_on_forward_strand()}")

  println("\n=== Demo Complete ===")
}