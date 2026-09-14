///|
fn main {
  println("=== tximeta Demo ===")

  let data = @src.tximeta_sample_data()
  println("\n1. Imported data summary")
  println("   " + @src.tximeta_data_summary(data))

  println("\n2. Linked transcriptome provenance")
  match data.txome {
    Some(t) => {
      println("   organism: " + t.organism)
      println("   genome:   " + t.genome)
      println("   release:  " + t.release)
      println("   source:   " + t.source)
      println("   sha:      " + t.sha)
    }
    None => println("   unlinked")
  }

  println("\n3. Coldata per sample")
  for sample in data.coldata.keys() {
    match data.coldata.get(sample) {
      Some(cols) => {
        let mut line = sample + ": "
        for k in cols.keys() {
          match cols.get(k) {
            Some(v) => line = line + k + "=" + v + " "
            None => ()
          }
        }
        println("   " + line)
      }
      None => ()
    }
  }

  println("\n4. Merged metadata")
  let meta = @src.tximeta_metadata(data)
  for k in meta.keys() {
    match meta.get(k) {
      Some(v) => println("   " + k + " = " + v)
      None => ()
    }
  }

  println("\n5. Summarize to gene level")
  let tx2gene = @src.tximeta_sample_tx2gene()
  let (genes, counts, abundance) = @src.tximeta_summarize_to_gene(data, tx2gene~)
  println("   " + genes.length().to_string() + " genes")
  for i in 0..<genes.length() {
    let mut row = ""
    for j in 0..<counts[i].length() {
      row = row + counts[i][j].to_string() + " "
    }
    println("   " + genes[i] + ": counts=" + row)
  }

  println("\n6. Validation")
  println("   linked: " + @src.tximeta_validate_linked(data).to_string())

  println("\n7. Unlinked data fails validation")
  let unlinked = @src.TximetaData::new(
    sample_names=["s1"],
    transcript_names=["ENST1"],
    counts=[[5]],
    abundance=[[1.0]],
  )
  println("   linked: " + @src.tximeta_validate_linked(unlinked).to_string())

  println("\n8. Link to a custom transcriptome")
  let custom = @src.link_txome(
    unlinked,
    @src.TximetaTxome::new(
      sha="beefcafe",
      genome="GRCm39",
      organism="Mus musculus",
      release="M25",
      source="GENCODE",
    ),
  )
  println("   linked: " + @src.tximeta_validate_linked(custom).to_string())
}