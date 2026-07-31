///|
/// UniProt XML IO module tests.

///|
test "parse_uniprot_xml_empty" {
  let result = @src.parse_uniprot_xml("")
  assert_eq(result.length(), 0)
}

///|
test "parse_uniprot_xml_invalid" {
  let result = @src.parse_uniprot_xml("<not-valid/>")
  assert_eq(result.length(), 0)
}

///|
test "parse_uniprot_xml_single_entry" {
  let xml = "<uniprot>\n" +
    "    <entry>\n" +
    "      <accession>P12345</accession>\n" +
    "      <protein><recommendedName><fullName>Test Protein</fullName></recommendedName></protein>\n" +
    "      <gene><name>GENE1</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>9606</ncbi_taxonomy_id><scientificName>Homo sapiens</scientificName></organism>\n" +
    "      <sequence length=\"10\" mass=\"1234.5\">\nACGTACGTAC\n      </sequence>\n" +
    "    </entry>\n" +
    "  </uniprot>"
  let result = @src.parse_uniprot_xml(xml)
  assert_eq(result.length(), 1)
  assert_eq(result[0].accession, "P12345")
  assert_eq(result[0].protein_name, "Test Protein")
  assert_eq(result[0].sequence, "ACGTACGTAC")
}

///|
test "parse_uniprot_xml_multiple_entries" {
  let xml = "<uniprot>\n" +
    "    <entry>\n" +
    "      <accession>P12345</accession>\n" +
    "      <protein><recommendedName><fullName>Protein 1</fullName></recommendedName></protein>\n" +
    "      <gene><name>GENE1</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>9606</ncbi_taxonomy_id><scientificName>Homo sapiens</scientificName></organism>\n" +
    "      <sequence length=\"5\">ACGTG</sequence>\n" +
    "    </entry>\n" +
    "    <entry>\n" +
    "      <accession>P67890</accession>\n" +
    "      <protein><recommendedName><fullName>Protein 2</fullName></recommendedName></protein>\n" +
    "      <gene><name>GENE2</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>10090</ncbi_taxonomy_id><scientificName>Mus musculus</scientificName></organism>\n" +
    "      <sequence length=\"5\">TTTTT</sequence>\n" +
    "    </entry>\n" +
    "  </uniprot>"
  let result = @src.parse_uniprot_xml(xml)
  assert_eq(result.length(), 2)
  assert_eq(result[0].accession, "P12345")
  assert_eq(result[1].accession, "P67890")
}

///|
test "parse_uniprot_xml_with_features" {
  let xml = "<uniprot>\n" +
    "    <entry>\n" +
    "      <accession>P12345</accession>\n" +
    "      <protein><recommendedName><fullName>Test Protein</fullName></recommendedName></protein>\n" +
    "      <gene><name>GENE1</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>9606</ncbi_taxonomy_id><scientificName>Homo sapiens</scientificName></organism>\n" +
    "      <sequence length=\"10\">ACGTACGTAC</sequence>\n" +
    "      <feature type=\"chain\" description=\"Test feature\">\n" +
    "        <location><begin position=\"1\"/><end position=\"5\"/></location>\n" +
    "      </feature>\n" +
    "    </entry>\n" +
    "  </uniprot>"
  let result = @src.parse_uniprot_xml(xml)
  assert_eq(result.length(), 1)
  assert_eq(result[0].features.length(), 1)
}

///|
test "parse_uniprot_xml_with_dbreferences" {
  let xml = "<uniprot>\n" +
    "    <entry>\n" +
    "      <accession>P12345</accession>\n" +
    "      <protein><recommendedName><fullName>Test Protein</fullName></recommendedName></protein>\n" +
    "      <gene><name>GENE1</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>9606</ncbi_taxonomy_id><scientificName>Homo sapiens</scientificName></organism>\n" +
    "      <sequence length=\"5\">ACGTG</sequence>\n" +
    "      <dbReference type=\"PDB\" id=\"1ABC\">\n" +
    "        <property type=\"id\" value=\"1ABC\"/>\n" +
    "      </dbReference>\n" +
    "    </entry>\n" +
    "  </uniprot>"
  let result = @src.parse_uniprot_xml(xml)
  assert_eq(result.length(), 1)
  assert_eq(result[0].db_references.length(), 1)
}

///|
test "uniprot_to_seqrecord" {
  let entry = @src.UniprotEntry::new("P12345", "Test Protein", ["GENE1"], "Homo sapiens", "ACGTACGTAC")
  let (header, id, desc, kwargs) = @src.uniprot_to_seqrecord(entry)
  assert_true(header.length() > 0)
  assert_true(id.length() > 0)
  assert_true(desc.length() > 0)
}

///|
test "uniprot_entry_methods" {
  let entry = @src.UniprotEntry::new("P12345", "Test Protein", ["GENE1", "GENE2"], "Homo sapiens", "ACGTACGTAC")
  assert_eq(entry.accession, "P12345")
  assert_eq(entry.sequence, "ACGTACGTAC")
  assert_eq(entry.seq_length, 10)
}
