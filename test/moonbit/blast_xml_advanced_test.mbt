///|
fn blast_xml_test_parse(text : String) -> @src.BlastXmlDocument {
  @src.blast_xml_parse(text) catch {
    BlastXmlError(message) => abort("valid BLAST XML failed: " + message)
  }
}

///|
fn blast_xml_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.blast_xml_parse(text))
    false
  } catch {
    @src.BlastXmlError(_) => true
  }
}

///|
fn blast_xml_test_write(document : @src.BlastXmlDocument) -> String {
  @src.blast_xml_write(document) catch {
    BlastXmlError(message) => abort("valid BLAST XML write failed: " + message)
  }
}

///|
fn blast_xml_test_write_as(
  document : @src.BlastXmlDocument,
  format : @src.BlastXmlFormat,
) -> String {
  @src.blast_xml_write_as(document, format) catch {
    BlastXmlError(message) =>
      abort("valid BLAST XML conversion failed: " + message)
  }
}

///|
fn blast_xml_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "BLAST XML value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn blast_xml_test_xml1() -> String {
  "<?xml version=\"1.0\"?>\n" +
  "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n" +
  "<BlastOutput>\n" +
  "  <BlastOutput_program>blastp</BlastOutput_program>\n" +
  "  <BlastOutput_version>BLASTP 2.15.0+</BlastOutput_version>\n" +
  "  <BlastOutput_reference>A. Sch&amp;auml;ffer &quot;BLAST&quot;</BlastOutput_reference>\n" +
  "  <BlastOutput_db>nr</BlastOutput_db>\n" +
  "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>\n" +
  "  <BlastOutput_query-def>official-style protein query</BlastOutput_query-def>\n" +
  "  <BlastOutput_query-len>12</BlastOutput_query-len>\n" +
  "  <BlastOutput_param>\n" +
  "    <Parameters>\n" +
  "      <Parameters_matrix>BLOSUM62</Parameters_matrix>\n" +
  "      <Parameters_expect>0.05</Parameters_expect>\n" +
  "      <Parameters_gap-open>11</Parameters_gap-open>\n" +
  "      <Parameters_gap-extend>1</Parameters_gap-extend>\n" +
  "      <Parameters_filter>F</Parameters_filter>\n" +
  "    </Parameters>\n" +
  "  </BlastOutput_param>\n" +
  "  <BlastOutput_iterations>\n" +
  "    <Iteration>\n" +
  "      <Iteration_iter-num>1</Iteration_iter-num>\n" +
  "      <Iteration_query-ID>Query_1</Iteration_query-ID>\n" +
  "      <Iteration_query-def>official-style protein query</Iteration_query-def>\n" +
  "      <Iteration_query-len>12</Iteration_query-len>\n" +
  "      <Iteration_hits>\n" +
  "        <Hit>\n" +
  "          <Hit_num>1</Hit_num>\n" +
  "          <Hit_id>ref|TARGET1.1|</Hit_id>\n" +
  "          <Hit_def>alpha &amp; beta &lt;literal&gt;</Hit_def>\n" +
  "          <Hit_accession>TARGET1</Hit_accession>\n" +
  "          <Hit_len>20</Hit_len>\n" +
  "          <Hit_hsps>\n" +
  "            <Hsp>\n" +
  "              <Hsp_num>1</Hsp_num>\n" +
  "              <Hsp_bit-score>477.248</Hsp_bit-score>\n" +
  "              <Hsp_score>1227</Hsp_score>\n" +
  "              <Hsp_evalue>2.44722e-20</Hsp_evalue>\n" +
  "              <Hsp_query-from>1</Hsp_query-from>\n" +
  "              <Hsp_query-to>8</Hsp_query-to>\n" +
  "              <Hsp_hit-from>1</Hsp_hit-from>\n" +
  "              <Hsp_hit-to>9</Hsp_hit-to>\n" +
  "              <Hsp_query-frame>0</Hsp_query-frame>\n" +
  "              <Hsp_hit-frame>0</Hsp_hit-frame>\n" +
  "              <Hsp_identity>8</Hsp_identity>\n" +
  "              <Hsp_positive>8</Hsp_positive>\n" +
  "              <Hsp_gaps>1</Hsp_gaps>\n" +
  "              <Hsp_align-len>9</Hsp_align-len>\n" +
  "              <Hsp_qseq>MKKL-AAAA</Hsp_qseq>\n" +
  "              <Hsp_hseq>MKKLTAAAA</Hsp_hseq>\n" +
  "              <Hsp_midline>MKKL AAAA</Hsp_midline>\n" +
  "            </Hsp>\n" +
  "          </Hit_hsps>\n" +
  "        </Hit>\n" +
  "        <Hit>\n" +
  "          <Hit_num>2</Hit_num>\n" +
  "          <Hit_id>sp|TARGET2.1|</Hit_id>\n" +
  "          <Hit_def>second target</Hit_def>\n" +
  "          <Hit_accession>TARGET2</Hit_accession>\n" +
  "          <Hit_len>30</Hit_len>\n" +
  "          <Hit_hsps>\n" +
  "            <Hsp>\n" +
  "              <Hsp_num>1</Hsp_num>\n" +
  "              <Hsp_bit-score>42.5</Hsp_bit-score>\n" +
  "              <Hsp_score>100</Hsp_score>\n" +
  "              <Hsp_evalue>2e-5</Hsp_evalue>\n" +
  "              <Hsp_query-from>2</Hsp_query-from>\n" +
  "              <Hsp_query-to>5</Hsp_query-to>\n" +
  "              <Hsp_hit-from>10</Hsp_hit-from>\n" +
  "              <Hsp_hit-to>13</Hsp_hit-to>\n" +
  "              <Hsp_query-frame>0</Hsp_query-frame>\n" +
  "              <Hsp_hit-frame>0</Hsp_hit-frame>\n" +
  "              <Hsp_identity>4</Hsp_identity>\n" +
  "              <Hsp_positive>4</Hsp_positive>\n" +
  "              <Hsp_gaps>0</Hsp_gaps>\n" +
  "              <Hsp_align-len>4</Hsp_align-len>\n" +
  "              <Hsp_qseq>AAAA</Hsp_qseq>\n" +
  "              <Hsp_hseq>AAAA</Hsp_hseq>\n" +
  "            </Hsp>\n" +
  "          </Hit_hsps>\n" +
  "        </Hit>\n" +
  "      </Iteration_hits>\n" +
  "      <Iteration_stat>\n" +
  "        <Statistics>\n" +
  "          <Statistics_db-num>743106947</Statistics_db-num>\n" +
  "          <Statistics_db-len>285453916046</Statistics_db-len>\n" +
  "          <Statistics_hsp-len>91</Statistics_hsp-len>\n" +
  "          <Statistics_eff-space>6.556e13</Statistics_eff-space>\n" +
  "          <Statistics_kappa>0.041</Statistics_kappa>\n" +
  "          <Statistics_lambda>0.267</Statistics_lambda>\n" +
  "          <Statistics_entropy>0.14</Statistics_entropy>\n" +
  "        </Statistics>\n" +
  "      </Iteration_stat>\n" +
  "    </Iteration>\n" +
  "  </BlastOutput_iterations>\n" +
  "</BlastOutput>\n"
}

///|
fn blast_xml_test_xml2() -> String {
  "<?xml version=\"1.0\"?>\n" +
  "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\">\n" +
  "  <BlastOutput2>\n" +
  "    <report><Report>\n" +
  "      <program>blastp</program>\n" +
  "      <version>BLASTP 2.15.0+</version>\n" +
  "      <reference>A. Sch&amp;auml;ffer &quot;BLAST&quot;</reference>\n" +
  "      <search-target><Target><db>nr</db></Target></search-target>\n" +
  "      <params><Parameters>\n" +
  "        <matrix>BLOSUM62</matrix><expect>0.05</expect>\n" +
  "        <gap-open>11</gap-open><gap-extend>1</gap-extend>\n" +
  "        <filter>F</filter><cbs>2</cbs>\n" +
  "      </Parameters></params>\n" +
  "      <results><Results><search><Search>\n" +
  "        <query-id>Query_1</query-id>\n" +
  "        <query-title>official-style protein query</query-title>\n" +
  "        <query-len>12</query-len>\n" +
  "        <query-masking>\n" +
  "          <Range><from>2</from><to>3</to></Range>\n" +
  "          <Range><from>7</from><to>9</to></Range>\n" +
  "        </query-masking>\n" +
  "        <hits><Hit>\n" +
  "          <num>1</num><description>\n" +
  "            <HitDescr><id>ref|TARGET1.1|</id><accession>TARGET1</accession><title>primary target</title><taxid>1385</taxid><sciname>Bacillales</sciname></HitDescr>\n" +
  "            <HitDescr><id>sp|SECONDARY.1|</id><accession>SECONDARY</accession><title>secondary target</title><taxid>224308</taxid><sciname>Bacillus subtilis</sciname></HitDescr>\n" +
  "          </description><len>20</len><hsps><Hsp>\n" +
  "            <num>1</num><bit-score>477.248</bit-score><score>1227</score><evalue>2.44722e-20</evalue>\n" +
  "            <identity>8</identity><positive>8</positive>\n" +
  "            <query-from>1</query-from><query-to>8</query-to>\n" +
  "            <hit-from>1</hit-from><hit-to>9</hit-to>\n" +
  "            <align-len>9</align-len><gaps>1</gaps>\n" +
  "            <qseq>MKKL-AAAA</qseq><hseq>MKKLTAAAA</hseq><midline>MKKL AAAA</midline>\n" +
  "          </Hsp></hsps>\n" +
  "        </Hit></hits>\n" +
  "        <stat><Statistics>\n" +
  "          <db-num>743106947</db-num><db-len>285453916046</db-len><hsp-len>91</hsp-len>\n" +
  "          <eff-space>6.556e13</eff-space><kappa>0.041</kappa><lambda>0.267</lambda><entropy>0.14</entropy>\n" +
  "        </Statistics></stat>\n" +
  "      </Search></search></Results></results>\n" +
  "    </Report></report>\n" +
  "  </BlastOutput2>\n" +
  "</BlastXML2>\n"
}

///|
fn blast_xml_test_xml2_hsp(
  program : String,
  query_length : Int,
  target_length : Int,
  hsp : String,
) -> String {
  "<?xml version=\"1.0\"?>\n" +
  "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><BlastOutput2><report><Report>" +
  "<program>" +
  program +
  "</program><version>2.15.0+</version><reference>reference</reference>" +
  "<search-target><Target><db>db</db></Target></search-target>" +
  "<params><Parameters><expect>10</expect><gap-open>11</gap-open><gap-extend>1</gap-extend></Parameters></params>" +
  "<results><Results><search><Search><query-id>query</query-id><query-title>query</query-title><query-len>" +
  query_length.to_string() +
  "</query-len><hits><Hit><num>1</num><description><HitDescr><id>target</id><accession>T</accession><title>target</title></HitDescr></description><len>" +
  target_length.to_string() +
  "</len><hsps><Hsp>" +
  hsp +
  "</Hsp></hsps></Hit></hits></Search></search></Results></results>" +
  "</Report></report></BlastOutput2></BlastXML2>"
}

///|
fn blast_xml_test_blastn() -> String {
  blast_xml_test_xml2_hsp(
    "blastn",
    20,
    30,
    "<num>1</num><bit-score>20</bit-score><score>40</score><evalue>1e-6</evalue>" +
    "<identity>4</identity><query-from>1</query-from><query-to>4</query-to><query-strand>Plus</query-strand>" +
    "<hit-from>20</hit-from><hit-to>17</hit-to><hit-strand>Minus</hit-strand>" +
    "<align-len>4</align-len><gaps>0</gaps><qseq>ACGT</qseq><hseq>ACGT</hseq>",
  )
}

///|
fn blast_xml_test_blastx() -> String {
  blast_xml_test_xml2_hsp(
    "blastx",
    30,
    10,
    "<num>1</num><bit-score>30</bit-score><score>60</score><evalue>1e-8</evalue>" +
    "<identity>3</identity><query-from>15</query-from><query-to>7</query-to><query-frame>-1</query-frame>" +
    "<hit-from>2</hit-from><hit-to>4</hit-to><align-len>3</align-len><gaps>0</gaps>" +
    "<qseq>MKT</qseq><hseq>MKT</hseq>",
  )
}

///|
fn blast_xml_test_tblastn() -> String {
  blast_xml_test_xml2_hsp(
    "tblastn",
    10,
    30,
    "<num>1</num><bit-score>25</bit-score><score>50</score><evalue>1e-7</evalue>" +
    "<identity>2</identity><query-from>1</query-from><query-to>2</query-to>" +
    "<hit-from>4</hit-from><hit-to>9</hit-to><hit-frame>1</hit-frame>" +
    "<align-len>2</align-len><gaps>0</gaps><qseq>MK</qseq><hseq>MK</hseq>",
  )
}

///|
fn blast_xml_test_psiblast() -> String {
  blast_xml_test_xml2_hsp(
    "psiblast",
    103,
    89,
    "<num>1</num><bit-score>177.178</bit-score><score>448</score><evalue>2.3039e-58</evalue>" +
    "<identity>89</identity><query-from>0</query-from><query-to>0</query-to>" +
    "<hit-from>0</hit-from><hit-to>0</hit-to><qseq></qseq><hseq></hseq>",
  )
  .replace_all(
    old="<results><Results><search>",
    new="<results><Results><iterations><Iteration><iter-num>1</iter-num><search>",
  )
  .replace_all(
    old="</search></Results></results>",
    new="</search></Iteration></iterations></Results></results>",
  )
}

///|
fn blast_xml_test_tblastx() -> String {
  blast_xml_test_xml2_hsp(
    "tblastx",
    30,
    60,
    "<num>1</num><bit-score>35</bit-score><score>70</score><evalue>1e-9</evalue>" +
    "<identity>3</identity><query-from>20</query-from><query-to>12</query-to><query-frame>-2</query-frame>" +
    "<hit-from>4</hit-from><hit-to>12</hit-to><hit-frame>3</hit-frame>" +
    "<align-len>3</align-len><gaps>0</gaps><qseq>MKT</qseq><hseq>MKT</hseq>",
  )
}

///|
fn blast_xml_test_rpsblast() -> String {
  blast_xml_test_xml2_hsp(
    "rpsblast",
    20,
    30,
    "<num>1</num><bit-score>25</bit-score><score>50</score><evalue>1e-7</evalue>" +
    "<identity>4</identity><query-from>2</query-from><query-to>5</query-to>" +
    "<hit-from>10</hit-from><hit-to>13</hit-to><align-len>4</align-len><gaps>0</gaps>" +
    "<qseq>AAAA</qseq><hseq>AAAA</hseq>",
  )
}

///|
fn blast_xml_test_megablast() -> String {
  blast_xml_test_blastn().replace_all(old="blastn", new="megablast")
}

///|
fn blast_xml_test_multi_query() -> String {
  blast_xml_test_xml1().replace_all(
    old="    </Iteration>\n  </BlastOutput_iterations>",
    new="    </Iteration>\n" +
      "    <Iteration><Iteration_iter-num>2</Iteration_iter-num>" +
      "<Iteration_query-ID>Query_2</Iteration_query-ID>" +
      "<Iteration_query-def>no-hit query</Iteration_query-def>" +
      "<Iteration_query-len>7</Iteration_query-len>" +
      "<Iteration_hits></Iteration_hits><Iteration_message>No hits found</Iteration_message>" +
      "</Iteration>\n  </BlastOutput_iterations>",
  )
}

///|
fn blast_xml_test_multi_xml2(second_database : String) -> String {
  let declaration = "<?xml version=\"1.0\"?>\n"
  let root = "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
  let body = blast_xml_test_xml2()
    .replace_all(old=declaration, new="")
    .replace_all(old=root, new="")
    .replace_all(old="</BlastXML2>\n", new="")
  let second = body
    .replace_all(old="Query_1", new="Query_2")
    .replace_all(old="official-style protein query", new="second protein query")
    .replace_all(old="<db>nr</db>", new="<db>" + second_database + "</db>")
  declaration + root + body + second + "</BlastXML2>\n"
}

///|
test "modern Bio.Blast parses XML1 headers" {
  let document = blast_xml_test_parse(blast_xml_test_xml1())
  assert_eq(document.format(), @src.BlastXmlFormat::Xml1)
  assert_eq(document.program(), "blastp")
  assert_eq(document.version(), "BLASTP 2.15.0+")
  assert_eq(document.database(), "nr")
}

///|
test "modern Bio.Blast decodes XML and historical HTML entities" {
  let document = blast_xml_test_parse(blast_xml_test_xml1())
  assert_eq(document.reference(), "A. Schäffer \"BLAST\"")
  assert_eq(
    document.records()[0].hits()[0].primary_description().title(),
    "alpha & beta <literal>",
  )
}

///|
test "modern Bio.Blast parses XML1 parameters" {
  let parameters = blast_xml_test_parse(blast_xml_test_xml1()).parameters()
  assert_eq(parameters.matrix(), Some("BLOSUM62"))
  blast_xml_test_close(parameters.expect(), 0.05, 1.0e-12)
  assert_eq(parameters.gap_open(), Some(11))
  assert_eq(parameters.gap_extend(), Some(1))
  assert_eq(parameters.filter(), Some("F"))
}

///|
test "modern Bio.Blast parses large database statistics" {
  let statistics = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
    .statistics()
    .unwrap()
  assert_eq(statistics.database_sequences(), 743106947)
  blast_xml_test_close(statistics.database_letters(), 285453916046.0, 0.5)
  blast_xml_test_close(statistics.effective_search_space(), 6.556e13, 1.0)
}

///|
test "modern Bio.Blast builds compressed gap coordinates" {
  let coordinates = blast_xml_test_parse(blast_xml_test_xml1()).records()[0].hits()[0].hsps()[0].coordinates()
  assert_eq(coordinates.length(), 4)
  assert_eq(coordinates[0].target(), 0)
  assert_eq(coordinates[0].query(), 0)
  assert_eq(coordinates[1].target(), 4)
  assert_eq(coordinates[1].query(), 4)
  assert_eq(coordinates[2].target(), 5)
  assert_eq(coordinates[2].query(), 4)
  assert_eq(coordinates[3].target(), 9)
  assert_eq(coordinates[3].query(), 8)
}

///|
test "modern Bio.Blast calculates identity fraction" {
  let hsp = blast_xml_test_parse(blast_xml_test_xml1()).records()[0].hits()[0].hsps()[0]
  blast_xml_test_close(hsp.identity_fraction(), 8.0 / 9.0, 1.0e-12)
}

///|
test "modern Bio.Blast selects the best HSP and hit" {
  let record = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
  assert_eq(record.best_hit().unwrap().number(), 1)
  assert_eq(record.hits()[0].best_hsp().unwrap().number(), 1)
}

///|
test "modern Bio.Blast looks up hits by accession" {
  let record = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
  assert_eq(record.hit("TARGET2").unwrap().number(), 2)
  assert_eq(record.hit("missing"), None)
}

///|
test "modern Bio.Blast filters HSPs by E-value" {
  let record = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
  let filtered = record.filter_evalue(1.0e-10)
  assert_eq(filtered.hits().length(), 1)
  assert_eq(filtered.hits()[0].number(), 1)
}

///|
test "modern Bio.Blast summarizes records, hits, and HSPs" {
  let document = blast_xml_test_parse(blast_xml_test_multi_query())
  assert_eq(document.records().length(), 2)
  assert_eq(document.total_hits(), 2)
  assert_eq(document.total_hsps(), 2)
  assert_eq(
    document.record("Query_2").unwrap().message(),
    Some("No hits found"),
  )
}

///|
test "modern Bio.Blast does not leak a global sequence to another query" {
  let text = blast_xml_test_multi_query().replace_all(
    old="  <BlastOutput_param>",
    new="  <BlastOutput_query-seq>AAAAAAAAAAAA</BlastOutput_query-seq>\n  <BlastOutput_param>",
  )
  let records = blast_xml_test_parse(text).records()
  assert_eq(records[0].query().sequence(), Some("AAAAAAAAAAAA"))
  assert_eq(records[1].query().sequence(), None)
}

///|
test "modern Bio.Blast parses XML2 namespaces and multiple descriptions" {
  let document = blast_xml_test_parse(blast_xml_test_xml2())
  assert_eq(document.format(), @src.BlastXmlFormat::Xml2)
  let descriptions = document.records()[0].hits()[0].descriptions()
  assert_eq(descriptions.length(), 2)
  assert_eq(descriptions[1].accession(), "SECONDARY")
  assert_eq(descriptions[1].taxid(), Some(224308))
  assert_eq(descriptions[1].scientific_name(), Some("Bacillus subtilis"))
}

///|
test "modern Bio.Blast parses every XML2 masking range" {
  let masks = blast_xml_test_parse(blast_xml_test_xml2()).records()[0]
    .query()
    .masks()
  assert_eq(masks.length(), 2)
  assert_eq(masks[0].start(), 1)
  assert_eq(masks[0].end(), 3)
  assert_eq(masks[1].start(), 6)
  assert_eq(masks[1].end(), 9)
}

///|
test "modern Bio.Blast parses XML2-only parameters" {
  let parameters = blast_xml_test_parse(blast_xml_test_xml2()).parameters()
  assert_eq(parameters.composition_based_statistics(), Some(2))
  assert_eq(parameters.inclusion_threshold(), None)
}

///|
test "modern Bio.Blast builds reverse blastn coordinates" {
  let hsp = blast_xml_test_parse(blast_xml_test_blastn()).records()[0].hits()[0].hsps()[0]
  let coordinates = hsp.coordinates()
  assert_eq(hsp.query_strand(), Some("Plus"))
  assert_eq(hsp.target_strand(), Some("Minus"))
  assert_eq(coordinates[0].target(), 20)
  assert_eq(coordinates[0].query(), 0)
  assert_eq(coordinates[1].target(), 16)
  assert_eq(coordinates[1].query(), 4)
}

///|
test "modern Bio.Blast builds negative blastx translated coordinates" {
  let hsp = blast_xml_test_parse(blast_xml_test_blastx()).records()[0].hits()[0].hsps()[0]
  assert_eq(hsp.query_frame(), Some(-1))
  assert_eq(hsp.query_coded_by(), Some("complement(query:7..15)"))
  let coordinates = hsp.coordinates()
  assert_eq(coordinates[0].query(), 15)
  assert_eq(coordinates[1].query(), 6)
  assert_eq(coordinates[0].target(), 1)
  assert_eq(coordinates[1].target(), 4)
}

///|
test "modern Bio.Blast builds tblastn target translation coordinates" {
  let hsp = blast_xml_test_parse(blast_xml_test_tblastn()).records()[0].hits()[0].hsps()[0]
  assert_eq(hsp.target_frame(), Some(1))
  assert_eq(hsp.target_coded_by(), Some("target:4..9"))
  let coordinates = hsp.coordinates()
  assert_eq(coordinates[0].target(), 3)
  assert_eq(coordinates[1].target(), 9)
}

///|
test "modern Bio.Blast accepts PSI-BLAST empty alignments" {
  let hsp = blast_xml_test_parse(blast_xml_test_psiblast()).records()[0].hits()[0].hsps()[0]
  assert_eq(hsp.alignment_length(), None)
  assert_eq(hsp.query_sequence(), "")
  assert_eq(hsp.target_sequence(), "")
  assert_eq(hsp.coordinates().length(), 0)
}

///|
test "modern Bio.Blast parses TBLASTX translated axes" {
  let hsp = blast_xml_test_parse(blast_xml_test_tblastx()).records()[0].hits()[0].hsps()[0]
  assert_eq(hsp.query_frame(), Some(-2))
  assert_eq(hsp.target_frame(), Some(3))
  assert_eq(hsp.query_coded_by(), Some("complement(query:12..20)"))
  assert_eq(hsp.target_coded_by(), Some("target:4..12"))
  let coordinates = hsp.coordinates()
  assert_eq(coordinates[0].query(), 20)
  assert_eq(coordinates[1].query(), 11)
  assert_eq(coordinates[0].target(), 3)
  assert_eq(coordinates[1].target(), 12)
}

///|
test "modern Bio.Blast parses RPS-BLAST protein coordinates" {
  let hsp = blast_xml_test_parse(blast_xml_test_rpsblast()).records()[0].hits()[0].hsps()[0]
  assert_eq(hsp.query_frame(), None)
  assert_eq(hsp.target_frame(), None)
  assert_eq(hsp.coordinates()[0].query(), 1)
  assert_eq(hsp.coordinates()[1].target(), 13)
}

///|
test "modern Bio.Blast parses megablast strands" {
  let document = blast_xml_test_parse(blast_xml_test_megablast())
  assert_eq(document.program(), "megablast")
  assert_eq(
    document.records()[0].hits()[0].hsps()[0].target_strand(),
    Some("Minus"),
  )
}

///|
test "modern Bio.Blast parses XML2 multiple reports" {
  let document = blast_xml_test_parse(blast_xml_test_multi_xml2("nr"))
  assert_eq(document.records().length(), 2)
  assert_eq(document.records()[0].number(), 1)
  assert_eq(document.records()[1].number(), 2)
  assert_eq(document.records()[1].query().id(), "Query_2")
}

///|
test "modern Bio.Blast XML1 canonical write round trips" {
  let document = blast_xml_test_parse(blast_xml_test_xml1())
  let first = blast_xml_test_write(document)
  let second = blast_xml_test_write(blast_xml_test_parse(first))
  assert_eq(second, first)
}

///|
test "modern Bio.Blast XML2 canonical write round trips" {
  let document = blast_xml_test_parse(blast_xml_test_xml2())
  let first = blast_xml_test_write(document)
  let second = blast_xml_test_write(blast_xml_test_parse(first))
  assert_eq(second, first)
}

///|
test "modern Bio.Blast writes PSI-BLAST iterations wrapper" {
  let written = blast_xml_test_write(
    blast_xml_test_parse(blast_xml_test_psiblast()),
  )
  assert_true(written.contains("<iterations>"))
  assert_true(written.contains("<Iteration>"))
  assert_eq(blast_xml_test_parse(written).records()[0].number(), 1)
}

///|
test "modern Bio.Blast writes one XML2 report per ordinary query" {
  let written = blast_xml_test_write(
    blast_xml_test_parse(blast_xml_test_multi_xml2("nr")),
  )
  assert_true(written.contains("</BlastOutput2>\n<BlastOutput2>"))
  let reparsed = blast_xml_test_parse(written)
  assert_eq(reparsed.records().length(), 2)
  assert_eq(reparsed.records()[1].query().id(), "Query_2")
  assert_eq(blast_xml_test_write(reparsed), written)
}

///|
test "modern Bio.Blast writes escaped XML text" {
  let written = blast_xml_test_write(
    blast_xml_test_parse(blast_xml_test_xml1()),
  )
  assert_true(written.contains("&amp; beta &lt;literal&gt;"))
  assert_eq(
    blast_xml_test_parse(written).records()[0].hits()[0]
    .primary_description()
    .title(),
    "alpha & beta <literal>",
  )
}

///|
test "modern Bio.Blast converts XML1 to XML2 semantically" {
  let converted = blast_xml_test_write_as(
    blast_xml_test_parse(blast_xml_test_xml1()),
    @src.BlastXmlFormat::Xml2,
  )
  let document = blast_xml_test_parse(converted)
  assert_eq(document.format(), @src.BlastXmlFormat::Xml2)
  assert_eq(document.total_hits(), 2)
  assert_eq(document.records()[0].hits()[0].descriptions().length(), 1)
}

///|
test "modern Bio.Blast collapses XML2 descriptions in XML1 conversion" {
  let source = blast_xml_test_xml2()
    .replace_all(old="<filter>F</filter><cbs>2</cbs>", new="<filter>F</filter>")
    .replace_all(
      old="        <query-masking>\n" +
        "          <Range><from>2</from><to>3</to></Range>\n" +
        "          <Range><from>7</from><to>9</to></Range>\n" +
        "        </query-masking>\n",
      new="",
    )
  let converted = blast_xml_test_write_as(
    blast_xml_test_parse(source),
    @src.BlastXmlFormat::Xml1,
  )
  let title = blast_xml_test_parse(converted).records()[0].hits()[0]
    .primary_description()
    .title()
  assert_true(title.contains("primary target"))
  assert_true(title.contains("sp|SECONDARY.1| secondary target"))
}

///|
test "modern Bio.Blast rejects lossy XML2 to XML1 conversion" {
  let rejected = try {
    ignore(
      @src.blast_xml_write_as(
        blast_xml_test_parse(blast_xml_test_xml2()),
        @src.BlastXmlFormat::Xml1,
      ),
    )
    false
  } catch {
    @src.BlastXmlError(_) => true
  }
  assert_true(rejected)
}

///|
test "modern Bio.Blast records are defensive copies" {
  let document = blast_xml_test_parse(blast_xml_test_multi_query())
  let records = document.records()
  ignore(records.pop())
  assert_eq(document.records().length(), 2)
}

///|
test "modern Bio.Blast hits are defensive copies" {
  let record = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
  let hits = record.hits()
  ignore(hits.pop())
  assert_eq(record.hits().length(), 2)
}

///|
test "modern Bio.Blast descriptions are defensive copies" {
  let hit = blast_xml_test_parse(blast_xml_test_xml2()).records()[0].hits()[0]
  let descriptions = hit.descriptions()
  ignore(descriptions.pop())
  assert_eq(hit.descriptions().length(), 2)
}

///|
test "modern Bio.Blast coordinates are defensive copies" {
  let hsp = blast_xml_test_parse(blast_xml_test_xml1()).records()[0].hits()[0].hsps()[0]
  let coordinates = hsp.coordinates()
  ignore(coordinates.pop())
  assert_eq(hsp.coordinates().length(), 4)
}

///|
test "modern Bio.Blast masks are defensive copies" {
  let query = blast_xml_test_parse(blast_xml_test_xml2()).records()[0].query()
  let masks = query.masks()
  ignore(masks.pop())
  assert_eq(query.masks().length(), 2)
}

///|
test "modern Bio.Blast rejects malformed closing tags" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="</BlastOutput>",
        new="</BlastOutput2>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects trailing data" {
  assert_true(blast_xml_test_rejects(blast_xml_test_xml1() + "trailing"))
}

///|
test "modern Bio.Blast rejects missing required fields" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="  <BlastOutput_program>blastp</BlastOutput_program>\n",
        new="",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects duplicate scalar fields" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="  <BlastOutput_program>blastp</BlastOutput_program>\n",
        new="  <BlastOutput_program>blastp</BlastOutput_program>\n" +
          "  <BlastOutput_program>blastp</BlastOutput_program>\n",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects unsupported programs" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(old="blastp", new="diamond"),
    ),
  )
}

///|
test "modern Bio.Blast rejects negative thresholds" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Parameters_expect>0.05</Parameters_expect>",
        new="<Parameters_expect>-0.05</Parameters_expect>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects nonnumeric fields" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hit_len>20</Hit_len>",
        new="<Hit_len>twenty</Hit_len>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects non-finite numbers" {
  for value in ["NaN", "inf", "-inf"] {
    assert_true(
      blast_xml_test_rejects(
        blast_xml_test_xml1().replace_all(
          old="<Hsp_evalue>2.44722e-20</Hsp_evalue>",
          new="<Hsp_evalue>" + value + "</Hsp_evalue>",
        ),
      ),
    )
  }
}

///|
test "modern Bio.Blast rejects inconsistent XML2 report headers" {
  assert_true(blast_xml_test_rejects(blast_xml_test_multi_xml2("swissprot")))
}

///|
test "modern Bio.Blast rejects sequence length mismatches" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_qseq>MKKL-AAAA</Hsp_qseq>",
        new="<Hsp_qseq>MKKL-AA</Hsp_qseq>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects double-gap columns" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1()
      .replace_all(old="MKKL-AAAA", new="MKKL-AAAA")
      .replace_all(
        old="<Hsp_hseq>MKKLTAAAA</Hsp_hseq>",
        new="<Hsp_hseq>MKKL-AAAA</Hsp_hseq>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects inconsistent gap counts" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_gaps>1</Hsp_gaps>",
        new="<Hsp_gaps>0</Hsp_gaps>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects identity above alignment length" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_identity>8</Hsp_identity>",
        new="<Hsp_identity>10</Hsp_identity>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects positive count below identity" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_positive>8</Hsp_positive>",
        new="<Hsp_positive>7</Hsp_positive>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects coordinates outside query bounds" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_query-to>8</Hsp_query-to>",
        new="<Hsp_query-to>13</Hsp_query-to>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects coordinate span mismatches" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml1().replace_all(
        old="<Hsp_hit-to>9</Hsp_hit-to>",
        new="<Hsp_hit-to>8</Hsp_hit-to>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects invalid blastx frames" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_blastx().replace_all(
        old="<query-frame>-1</query-frame>",
        new="<query-frame>0</query-frame>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects missing translated frames" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_blastx().replace_all(
        old="<query-frame>-1</query-frame>",
        new="",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects invalid strands" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_blastn().replace_all(
        old="<hit-strand>Minus</hit-strand>",
        new="<hit-strand>Reverse</hit-strand>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects strand-coordinate conflicts" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_blastn().replace_all(
        old="<hit-strand>Minus</hit-strand>",
        new="<hit-strand>Plus</hit-strand>",
      ),
    ),
  )
}

///|
test "modern Bio.Blast rejects masking outside query bounds" {
  assert_true(
    blast_xml_test_rejects(
      blast_xml_test_xml2().replace_all(old="<to>9</to>", new="<to>13</to>"),
    ),
  )
}

///|
test "modern Bio.Blast rejects duplicate HSP numbers" {
  let text = blast_xml_test_xml1().replace_all(
    old="            </Hsp>\n          </Hit_hsps>",
    new="            </Hsp>\n" +
      "            <Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>1</Hsp_bit-score><Hsp_score>1</Hsp_score><Hsp_evalue>1</Hsp_evalue>" +
      "<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>1</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>1</Hsp_hit-to>" +
      "<Hsp_query-frame>0</Hsp_query-frame><Hsp_hit-frame>0</Hsp_hit-frame><Hsp_identity>1</Hsp_identity>" +
      "<Hsp_gaps>0</Hsp_gaps><Hsp_align-len>1</Hsp_align-len><Hsp_qseq>A</Hsp_qseq><Hsp_hseq>A</Hsp_hseq></Hsp>\n" +
      "          </Hit_hsps>",
  )
  assert_true(blast_xml_test_rejects(text))
}

///|
test "modern Bio.Blast rejects duplicate hit numbers" {
  let text = blast_xml_test_xml1().replace_all(
    old="<Hit_num>2</Hit_num>",
    new="<Hit_num>1</Hit_num>",
  )
  assert_true(blast_xml_test_rejects(text))
}

///|
test "modern Bio.Blast rejects invalid E-value filter thresholds" {
  let record = blast_xml_test_parse(blast_xml_test_xml1()).records()[0]
  let rejected = try {
    ignore(record.filter_evalue(-1.0))
    false
  } catch {
    @src.BlastXmlError(_) => true
  }
  assert_true(rejected)
}

///|
test "modern Bio.Blast parses CRLF documents" {
  let text = blast_xml_test_xml1().replace_all(old="\n", new="\r\n")
  assert_eq(blast_xml_test_parse(text).total_hits(), 2)
}

///|
test "modern Bio.Blast parses comments and CDATA" {
  let text = blast_xml_test_xml1().replace_all(
    old="<BlastOutput_version>BLASTP 2.15.0+</BlastOutput_version>",
    new="<!-- version --><BlastOutput_version><![CDATA[BLASTP 2.15.0+]]></BlastOutput_version>",
  )
  assert_eq(blast_xml_test_parse(text).version(), "BLASTP 2.15.0+")
}
