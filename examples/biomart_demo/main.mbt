///|
/// BioMart Demo
fn main {
  println("=== BioMart Demo ===")

  println("\n1. Create Human BioMart Dataset:")
  println("================================")
  let human_mart = @bio.create_human_mart()
  println("Name: " + human_mart.name)
  println("Description: " + human_mart.description)
  println("Host: " + human_mart.host)

  println("\n2. Gene ID Mapping:")
  println("===================")
  let ensembl_id = "ENSG00000130203"
  let gene_name = human_mart.get_mapping(ensembl_id)
  println(
    "ENSG00000130203 -> " +
    (if gene_name is Some(_) { gene_name.unwrap() } else { "Not found" }),
  )

  let mappings = human_mart.map_ids([
    "ENSG00000130203", "ENSG00000141510", "ENSG00000160072",
  ])
  let mut i = 0
  while i < mappings.length() {
    println(
      "  " +
      mappings[i].0 +
      " -> " +
      (if mappings[i].1 is Some(_) { mappings[i].1.unwrap() } else { "N/A" }),
    )
    i = i + 1
  }

  println("\n3. Gene Annotation:")
  println("===================")
  println(
    "Gene Name: " +
    (if human_mart.get_gene_name(ensembl_id) is Some(_) {
      human_mart.get_gene_name(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "Description: " +
    (if human_mart.get_gene_description(ensembl_id) is Some(_) {
      human_mart.get_gene_description(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "Chromosome: " +
    (if human_mart.get_chromosome(ensembl_id) is Some(_) {
      human_mart.get_chromosome(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "Start: " +
    (if human_mart.get_start(ensembl_id) is Some(_) {
      human_mart.get_start(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "End: " +
    (if human_mart.get_end(ensembl_id) is Some(_) {
      human_mart.get_end(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "Strand: " +
    (if human_mart.get_strand(ensembl_id) is Some(_) {
      human_mart.get_strand(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "Biotype: " +
    (if human_mart.get_biotype(ensembl_id) is Some(_) {
      human_mart.get_biotype(ensembl_id).unwrap()
    } else {
      "N/A"
    }),
  )

  println("\n4. External Database Mapping:")
  println("=============================")
  println(
    "UniProt: " +
    (if human_mart.get_external_db(ensembl_id, "uniprot") is Some(_) {
      human_mart.get_external_db(ensembl_id, "uniprot").unwrap()
    } else {
      "N/A"
    }),
  )
  println(
    "HGNC Symbol: " +
    (if human_mart.get_external_db(ensembl_id, "hgnc_symbol") is Some(_) {
      human_mart.get_external_db(ensembl_id, "hgnc_symbol").unwrap()
    } else {
      "N/A"
    }),
  )

  println("\n5. Batch Query:")
  println("===============")
  let query_ids = ["ENSG00000130203", "ENSG00000141510", "ENSG00000160072"]
  let attributes = ["gene_name", "chromosome", "biotype"]
  let results = human_mart.query(query_ids, attributes)

  let mut j = 0
  while j < results.length() {
    let row = results[j]
    let mut line = ""
    let mut k = 0
    while k < row.length() {
      if k > 0 {
        line = line + "\t"
      }
      line = line + row[k]
      k = k + 1
    }
    println(line)
    j = j + 1
  }

  println("\n6. Create Mouse BioMart Dataset:")
  println("================================")
  let mouse_mart = @bio.create_mouse_mart()
  println("Name: " + mouse_mart.name)
  println("Description: " + mouse_mart.description)

  let mouse_id = "ENSMUSG00000059552"
  let mouse_name = mouse_mart.get_mapping(mouse_id)
  println(
    mouse_id +
    " -> " +
    (if mouse_name is Some(_) { mouse_name.unwrap() } else { "Not found" }),
  )

  println("\n=== Demo Complete ===")
}
