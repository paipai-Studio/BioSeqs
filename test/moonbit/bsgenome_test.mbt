///|
test "bsgenome_create" {
  let genome = @src.create_demo_genome()
  assert_eq(genome.name, "Homo sapiens")
}

///|
test "bsgenome_chromosome_names" {
  let genome = @src.create_demo_genome()
  let names = genome.get_chromosome_names()
  assert_eq(names.length(), 5)
}

///|
test "bsgenome_get_sequence" {
  let genome = @src.create_demo_genome()
  let seq = genome.get_sequence("chr1", 1, 10)
  assert_true(seq is Some(_))
}

///|
test "bsgenome_extract_gene_forward" {
  let genome = @src.create_demo_genome()
  let gene = genome.extract_gene("chr1", 1, 10, "+")
  assert_true(gene is Some(_))
}

///|
test "bsgenome_extract_gene_reverse" {
  let genome = @src.create_demo_genome()
  let gene = genome.extract_gene("chr1", 1, 10, "-")
  assert_true(gene is Some(_))
}
