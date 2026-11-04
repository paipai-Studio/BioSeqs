///|
fn blast_demo_parse(text : String) -> @bio.BlastXmlDocument {
  @bio.blast_xml_parse(text) catch {
    BlastXmlError(message) => abort("BLAST XML parse failed: " + message)
  }
}

///|
fn blast_demo_write(document : @bio.BlastXmlDocument) -> String {
  @bio.blast_xml_write(document) catch {
    BlastXmlError(message) => abort("BLAST XML write failed: " + message)
  }
}

///|
fn blast_demo_write_xml2(document : @bio.BlastXmlDocument) -> String {
  @bio.blast_xml_write_as(document, @bio.BlastXmlFormat::Xml2) catch {
    BlastXmlError(message) => abort("BLAST XML conversion failed: " + message)
  }
}

///|
fn blast_demo_xml1() -> String {
  "<?xml version=\"1.0\"?>\n" +
  "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n" +
  "<BlastOutput>" +
  "<BlastOutput_program>blastp</BlastOutput_program>" +
  "<BlastOutput_version>BLASTP 2.15.0+</BlastOutput_version>" +
  "<BlastOutput_reference>offline demo</BlastOutput_reference>" +
  "<BlastOutput_db>swissprot</BlastOutput_db>" +
  "<BlastOutput_query-ID>Query_1</BlastOutput_query-ID>" +
  "<BlastOutput_query-def>protein query</BlastOutput_query-def>" +
  "<BlastOutput_query-len>12</BlastOutput_query-len>" +
  "<BlastOutput_param><Parameters>" +
  "<Parameters_matrix>BLOSUM62</Parameters_matrix>" +
  "<Parameters_expect>0.001</Parameters_expect>" +
  "<Parameters_gap-open>11</Parameters_gap-open>" +
  "<Parameters_gap-extend>1</Parameters_gap-extend>" +
  "</Parameters></BlastOutput_param>" +
  "<BlastOutput_iterations><Iteration>" +
  "<Iteration_iter-num>1</Iteration_iter-num>" +
  "<Iteration_query-ID>Query_1</Iteration_query-ID>" +
  "<Iteration_query-def>protein query</Iteration_query-def>" +
  "<Iteration_query-len>12</Iteration_query-len>" +
  "<Iteration_hits><Hit>" +
  "<Hit_num>1</Hit_num><Hit_id>sp|DEMO.1|</Hit_id>" +
  "<Hit_def>demonstration target</Hit_def><Hit_accession>DEMO</Hit_accession>" +
  "<Hit_len>20</Hit_len><Hit_hsps><Hsp>" +
  "<Hsp_num>1</Hsp_num><Hsp_bit-score>42.5</Hsp_bit-score>" +
  "<Hsp_score>100</Hsp_score><Hsp_evalue>2e-12</Hsp_evalue>" +
  "<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>7</Hsp_query-to>" +
  "<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>8</Hsp_hit-to>" +
  "<Hsp_query-frame>0</Hsp_query-frame><Hsp_hit-frame>0</Hsp_hit-frame>" +
  "<Hsp_identity>7</Hsp_identity><Hsp_positive>7</Hsp_positive>" +
  "<Hsp_gaps>1</Hsp_gaps><Hsp_align-len>8</Hsp_align-len>" +
  "<Hsp_qseq>MKKL-AAA</Hsp_qseq><Hsp_hseq>MKKLTAAA</Hsp_hseq>" +
  "<Hsp_midline>MKKL AAA</Hsp_midline></Hsp></Hit_hsps></Hit></Iteration_hits>" +
  "<Iteration_stat><Statistics><Statistics_db-num>482816</Statistics_db-num>" +
  "<Statistics_db-len>183558113</Statistics_db-len>" +
  "<Statistics_hsp-len>72</Statistics_hsp-len>" +
  "<Statistics_eff-space>4627826878</Statistics_eff-space>" +
  "<Statistics_kappa>0.041</Statistics_kappa>" +
  "<Statistics_lambda>0.267</Statistics_lambda>" +
  "<Statistics_entropy>0.14</Statistics_entropy>" +
  "</Statistics></Iteration_stat></Iteration></BlastOutput_iterations>" +
  "</BlastOutput>\n"
}

///|
fn blast_demo_xml2() -> String {
  "<?xml version=\"1.0\"?>\n" +
  "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><BlastOutput2>" +
  "<report><Report><program>blastx</program><version>BLASTX 2.15.0+</version>" +
  "<reference>offline demo</reference><search-target><Target><db>swissprot</db>" +
  "</Target></search-target><params><Parameters><matrix>BLOSUM62</matrix>" +
  "<expect>0.001</expect><gap-open>11</gap-open><gap-extend>1</gap-extend>" +
  "<cbs>2</cbs></Parameters></params><results><Results><search><Search>" +
  "<query-id>translated_query</query-id><query-title>negative frame query</query-title>" +
  "<query-len>30</query-len><hits><Hit><num>1</num><description>" +
  "<HitDescr><id>sp|PRIMARY.1|</id><accession>PRIMARY</accession>" +
  "<title>primary target</title><taxid>83333</taxid></HitDescr>" +
  "<HitDescr><id>sp|SECONDARY.1|</id><accession>SECONDARY</accession>" +
  "<title>secondary target</title><taxid>224308</taxid>" +
  "<sciname>Bacillus subtilis</sciname></HitDescr></description><len>10</len>" +
  "<hsps><Hsp><num>1</num><bit-score>30</bit-score><score>60</score>" +
  "<evalue>1e-8</evalue><identity>3</identity>" +
  "<query-from>15</query-from><query-to>7</query-to><query-frame>-1</query-frame>" +
  "<hit-from>2</hit-from><hit-to>4</hit-to><align-len>3</align-len>" +
  "<gaps>0</gaps><qseq>MKT</qseq><hseq>MKT</hseq></Hsp></hsps>" +
  "</Hit></hits></Search></search></Results></results></Report></report>" +
  "</BlastOutput2></BlastXML2>\n"
}

///|
fn main {
  println("=== Biopython Bio.Blast XML1/XML2 Offline Demo ===")
  let xml1 = blast_demo_parse(blast_demo_xml1())
  let record = xml1.records()[0]
  let best = record.best_hit().unwrap()
  let hsp = best.best_hsp().unwrap()
  let statistics = record.statistics().unwrap()
  println("\n1. XML1 records, best hit, and statistics")
  println(
    "  program=" +
    xml1.program() +
    ", records=" +
    xml1.records().length().to_string() +
    ", hits=" +
    xml1.total_hits().to_string(),
  )
  println(
    "  best=" +
    best.primary_description().accession() +
    ", E-value=" +
    hsp.evalue().to_string() +
    ", identity=" +
    hsp.identity_fraction().to_string(),
  )
  println(
    "  database sequences=" +
    statistics.database_sequences().to_string() +
    ", letters=" +
    statistics.database_letters().to_string(),
  )

  println("\n2. Compressed alignment path")
  for coordinate in hsp.coordinates() {
    println(
      "  target=" +
      coordinate.target().to_string() +
      ", query=" +
      coordinate.query().to_string(),
    )
  }

  let xml2 = blast_demo_parse(blast_demo_xml2())
  let translated_hit = xml2.records()[0].hits()[0]
  let translated_hsp = translated_hit.hsps()[0]
  println("\n3. XML2 descriptions and translated coordinates")
  println(
    "  descriptions=" +
    translated_hit.descriptions().length().to_string() +
    ", secondary taxon=" +
    translated_hit.descriptions()[1].scientific_name().unwrap_or("unknown"),
  )
  println(
    "  query frame=" +
    translated_hsp.query_frame().unwrap().to_string() +
    ", coded_by=" +
    translated_hsp.query_coded_by().unwrap_or("undefined"),
  )

  let canonical = blast_demo_write(xml2)
  let converted = blast_demo_write_xml2(xml1)
  println("\n4. Canonical writing and dialect conversion")
  println(
    "  XML2 round-trip stable: " +
    (blast_demo_write(blast_demo_parse(canonical)) == canonical).to_string(),
  )
  println(
    "  XML1 -> XML2 records: " +
    blast_demo_parse(converted).records().length().to_string(),
  )
  println(
    "\nThe demo parses and writes existing artifacts; it does not use NCBI or launch BLAST.",
  )
}
