///|
/// SeqXML I/O module tests.
// Tests for the SeqXML parser/writer mirroring BioPython's Bio.SeqIO.SeqXmlIO.
// Run `moon test` to execute.

// ============================================================================
// SeqXmlType enum tests
// ============================================================================

///|
test "seqxml_type_from_tag" {
  // Known tags map to their corresponding types.
  assert_eq(
    @src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("dna")),
    "dna",
  )
  assert_eq(
    @src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("rna")),
    "rna",
  )
  assert_eq(
    @src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("protein")),
    "protein",
  )
  // Unknown tags map to UnknownType.
  assert_eq(
    @src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("xyz")),
    "unknown",
  )
  assert_eq(
    @src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("")),
    "unknown",
  )
}

///|
test "seqxml_type_to_tag" {
  assert_eq(@src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("dna")), "dna")
  assert_eq(@src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("rna")), "rna")
  assert_eq(@src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("protein")), "protein")
  assert_eq(@src.SeqXmlType::to_tag(@src.SeqXmlType::from_tag("unknown")), "unknown")
}

///|
test "seqxml_type_description" {
  assert_eq(@src.SeqXmlType::description(@src.SeqXmlType::from_tag("dna")), "DNA")
  assert_eq(@src.SeqXmlType::description(@src.SeqXmlType::from_tag("rna")), "RNA")
  assert_eq(
    @src.SeqXmlType::description(@src.SeqXmlType::from_tag("protein")),
    "Protein",
  )
  assert_eq(
    @src.SeqXmlType::description(@src.SeqXmlType::from_tag("unknown")),
    "Unknown",
  )
}

///|
test "seqxml_type_from_tag_description_roundtrip" {
  // from_tag -> description should give consistent human-readable names.
  assert_eq(
    @src.SeqXmlType::description(@src.SeqXmlType::from_tag("dna")),
    "DNA",
  )
  assert_eq(
    @src.SeqXmlType::description(@src.SeqXmlType::from_tag("rna")),
    "RNA",
  )
  assert_eq(
    @src.SeqXmlType::description(@src.SeqXmlType::from_tag("protein")),
    "Protein",
  )
}

// ============================================================================
// SeqXmlEntry tests
// ============================================================================

///|
test "seqxml_entry_creation" {
  let entry = @src.SeqXmlEntry::new(
    "seq1",
    @src.SeqXmlType::from_tag("dna"),
    "ATGCATGC",
    description="A test DNA sequence",
  )
  assert_eq(entry.id(), "seq1")
  assert_eq(entry.description(), "A test DNA sequence")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "dna")
  assert_eq(entry.sequence(), "ATGCATGC")
  assert_eq(entry.properties().length(), 0)
}

///|
test "seqxml_entry_creation_without_description" {
  let entry = @src.SeqXmlEntry::new(
    "seq2",
    @src.SeqXmlType::from_tag("protein"),
    "MKLV",
  )
  assert_eq(entry.id(), "seq2")
  // Default description is the empty string.
  assert_eq(entry.description(), "")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "protein")
  assert_eq(entry.sequence(), "MKLV")
}

///|
test "seqxml_entry_add_property" {
  let entry = @src.SeqXmlEntry::new(
    "p1",
    @src.SeqXmlType::from_tag("protein"),
    "MKLVGV",
    description="A protein",
  )
  entry.add_property("gene", "BRCA1")
  entry.add_property("length", "6")
  let props = entry.properties()
  assert_eq(props.length(), 2)
  assert_eq(props[0].0, "gene")
  assert_eq(props[0].1, "BRCA1")
  assert_eq(props[1].0, "length")
  assert_eq(props[1].1, "6")
}

///|
test "seqxml_entry_to_string" {
  let entry = @src.SeqXmlEntry::new(
    "s1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="short seq",
  )
  let s = entry.to_string()
  assert_true(s.contains("SeqXmlEntry"))
  assert_true(s.contains("id=s1"))
  assert_true(s.contains("type=DNA"))
  assert_true(s.contains("length=4"))
  assert_true(s.contains("description=short seq"))
  assert_true(s.contains("properties=0"))
}

///|
test "seqxml_entry_to_string_with_properties" {
  let entry = @src.SeqXmlEntry::new(
    "p1",
    @src.SeqXmlType::from_tag("protein"),
    "MKLV",
  )
  entry.add_property("k1", "v1")
  entry.add_property("k2", "v2")
  let s = entry.to_string()
  assert_true(s.contains("type=Protein"))
  assert_true(s.contains("length=4"))
  assert_true(s.contains("properties=2"))
}

// ============================================================================
// SeqXmlDocument tests
// ============================================================================

///|
test "seqxml_document_creation" {
  let doc = @src.SeqXmlDocument::new()
  assert_eq(doc.n_entries(), 0)
  assert_eq(doc.entries().length(), 0)
  assert_eq(doc.species_name(), "")
  assert_eq(doc.species_taxid(), "")
  assert_eq(doc.source_db(), "")
}

///|
test "seqxml_document_set_species" {
  let doc = @src.SeqXmlDocument::new()
  doc.set_species("Homo sapiens", "9606")
  assert_eq(doc.species_name(), "Homo sapiens")
  assert_eq(doc.species_taxid(), "9606")
}

///|
test "seqxml_document_set_source_db" {
  let doc = @src.SeqXmlDocument::new()
  doc.set_source_db("UniProt")
  assert_eq(doc.source_db(), "UniProt")
}

///|
test "seqxml_document_add_entry" {
  let doc = @src.SeqXmlDocument::new()
  assert_eq(doc.n_entries(), 0)
  let e1 = @src.SeqXmlEntry::new(
    "s1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="first",
  )
  doc.add_entry(e1)
  assert_eq(doc.n_entries(), 1)
  let e2 = @src.SeqXmlEntry::new(
    "s2",
    @src.SeqXmlType::from_tag("protein"),
    "MKLV",
    description="second",
  )
  doc.add_entry(e2)
  assert_eq(doc.n_entries(), 2)
  let entries = doc.entries()
  assert_eq(entries[0].id(), "s1")
  assert_eq(entries[1].id(), "s2")
}

///|
test "seqxml_document_full_metadata" {
  let doc = @src.SeqXmlDocument::new()
  doc.set_species("Mus musculus", "10090")
  doc.set_source_db("Ensembl")
  let entry = @src.SeqXmlEntry::new(
    "g1",
    @src.SeqXmlType::from_tag("dna"),
    "ATGC",
    description="gene",
  )
  doc.add_entry(entry)
  assert_eq(doc.species_name(), "Mus musculus")
  assert_eq(doc.species_taxid(), "10090")
  assert_eq(doc.source_db(), "Ensembl")
  assert_eq(doc.n_entries(), 1)
  assert_eq(doc.entries()[0].id(), "g1")
}

// ============================================================================
// Parsing tests
// ============================================================================

///|
test "seqxml_parse_dna_entry" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"dna1\" desc=\"DNA test\">\n" +
    "    <dna>ATGCATGC</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "dna1")
  assert_eq(entry.description(), "DNA test")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "dna")
  assert_eq(entry.sequence(), "ATGCATGC")
  assert_eq(entry.properties().length(), 0)
}

///|
test "seqxml_parse_rna_entry" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"rna1\" desc=\"RNA test\">\n" +
    "    <rna>AUGCAUGC</rna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "rna1")
  assert_eq(entry.description(), "RNA test")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "rna")
  assert_eq(entry.sequence(), "AUGCAUGC")
}

///|
test "seqxml_parse_protein_entry" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"prot1\" desc=\"Protein test\">\n" +
    "    <protein>MKLVGV</protein>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "prot1")
  assert_eq(entry.description(), "Protein test")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "protein")
  assert_eq(entry.sequence(), "MKLVGV")
}

///|
test "seqxml_parse_with_species_sourcedb" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <species name=\"Homo sapiens\" ncbiTaxID=\"9606\"/>\n" +
    "  <sourceDB name=\"UniProt\"/>\n" +
    "  <entry id=\"s1\" desc=\"With metadata\">\n" +
    "    <dna>ACGT</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.species_name(), "Homo sapiens")
  assert_eq(doc.species_taxid(), "9606")
  assert_eq(doc.source_db(), "UniProt")
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "s1")
  assert_eq(entry.description(), "With metadata")
  assert_eq(entry.sequence(), "ACGT")
}

///|
test "seqxml_parse_with_properties" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"p1\" desc=\"With props\">\n" +
    "    <property name=\"gene\" value=\"BRCA1\"/>\n" +
    "    <property name=\"length\" value=\"100\"/>\n" +
    "    <protein>MKLV</protein>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "p1")
  assert_eq(entry.description(), "With props")
  let props = entry.properties()
  assert_eq(props.length(), 2)
  assert_eq(props[0].0, "gene")
  assert_eq(props[0].1, "BRCA1")
  assert_eq(props[1].0, "length")
  assert_eq(props[1].1, "100")
  assert_eq(entry.sequence(), "MKLV")
}

///|
test "seqxml_parse_multiple_entries" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"d1\"><dna>ACGT</dna></entry>\n" +
    "  <entry id=\"r1\"><rna>AUGU</rna></entry>\n" +
    "  <entry id=\"p1\"><protein>MKLV</protein></entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 3)
  let entries = doc.entries()
  assert_eq(entries[0].id(), "d1")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  assert_eq(entries[0].sequence(), "ACGT")
  assert_eq(entries[1].id(), "r1")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "rna")
  assert_eq(entries[1].sequence(), "AUGU")
  assert_eq(entries[2].id(), "p1")
  assert_eq(@src.SeqXmlType::to_tag(entries[2].seq_type()), "protein")
  assert_eq(entries[2].sequence(), "MKLV")
}

///|
test "seqxml_parse_entry_id" {
  // Entry tag with id attribute.
  match @src.seqxml_parse_entry_id("<entry id=\"seq1\">") {
    Some(id) => assert_eq(id, "seq1")
    None => fail("expected Some(id) for entry with id attribute")
  }
  // Entry tag with id and desc attributes.
  match @src.seqxml_parse_entry_id("<entry id=\"seq2\" desc=\"Test\">") {
    Some(id) => assert_eq(id, "seq2")
    None => fail("expected Some(id) for entry with id and desc")
  }
  // Self-closing entry tag with id.
  match @src.seqxml_parse_entry_id("<entry id=\"seq3\"/>") {
    Some(id) => assert_eq(id, "seq3")
    None => fail("expected Some(id) for self-closing entry tag")
  }
  // Entry tag without id attribute.
  match @src.seqxml_parse_entry_id("<entry>") {
    Some(_) => fail("expected None for entry without id")
    None => assert_true(true)
  }
  // Non-entry line.
  match @src.seqxml_parse_entry_id("not an entry line") {
    Some(_) => fail("expected None for non-entry line")
    None => assert_true(true)
  }
}

// ============================================================================
// Writing tests
// ============================================================================

///|
test "seqxml_entry_to_xml" {
  let entry = @src.SeqXmlEntry::new(
    "s1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="Test entry",
  )
  let xml = @src.seqxml_entry_to_xml(entry)
  assert_true(xml.contains("<entry id=\"s1\" desc=\"Test entry\">"))
  assert_true(xml.contains("<dna>ACGT</dna>"))
  assert_true(xml.contains("</entry>"))
}

///|
test "seqxml_entry_to_xml_no_description" {
  let entry = @src.SeqXmlEntry::new(
    "s2",
    @src.SeqXmlType::from_tag("protein"),
    "MKLV",
  )
  let xml = @src.seqxml_entry_to_xml(entry)
  // Without description, desc attribute is omitted.
  assert_true(xml.contains("<entry id=\"s2\">"))
  assert_true(xml.contains("<protein>MKLV</protein>"))
  assert_true(xml.contains("</entry>"))
}

///|
test "seqxml_entry_to_xml_with_properties" {
  let entry = @src.SeqXmlEntry::new(
    "p1",
    @src.SeqXmlType::from_tag("protein"),
    "MKLV",
    description="protein",
  )
  entry.add_property("gene", "BRCA1")
  entry.add_property("length", "4")
  let xml = @src.seqxml_entry_to_xml(entry)
  assert_true(xml.contains("<property name=\"gene\" value=\"BRCA1\"/>"))
  assert_true(xml.contains("<property name=\"length\" value=\"4\"/>"))
  assert_true(xml.contains("<protein>MKLV</protein>"))
}

///|
test "seqxml_document_to_xml" {
  let doc = @src.SeqXmlDocument::new()
  doc.set_species("Homo sapiens", "9606")
  doc.set_source_db("UniProt")
  let entry = @src.SeqXmlEntry::new(
    "s1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="Test",
  )
  doc.add_entry(entry)
  let xml = @src.seqxml_to_xml(doc)
  assert_true(xml.contains("<?xml version=\"1.0\""))
  assert_true(xml.contains("<seqXML>"))
  assert_true(xml.contains("</seqXML>"))
  assert_true(xml.contains("<species name=\"Homo sapiens\" ncbiTaxID=\"9606\"/>"))
  assert_true(xml.contains("<sourceDB name=\"UniProt\"/>"))
  assert_true(xml.contains("<entry id=\"s1\" desc=\"Test\">"))
  assert_true(xml.contains("<dna>ACGT</dna>"))
}

///|
test "seqxml_document_to_xml_no_metadata" {
  let doc = @src.SeqXmlDocument::new()
  let entry = @src.SeqXmlEntry::new(
    "x1",
    @src.SeqXmlType::from_tag("rna"),
    "AUGU",
  )
  doc.add_entry(entry)
  let xml = @src.seqxml_to_xml(doc)
  assert_true(xml.contains("<?xml version=\"1.0\""))
  assert_true(xml.contains("<seqXML>"))
  assert_true(xml.contains("</seqXML>"))
  // No species or sourceDB elements.
  assert_true(!xml.contains("<species"))
  assert_true(!xml.contains("<sourceDB"))
  assert_true(xml.contains("<entry id=\"x1\">"))
  assert_true(xml.contains("<rna>AUGU</rna>"))
}

// ============================================================================
// Round-trip tests
// ============================================================================

///|
test "seqxml_roundtrip_basic" {
  let original =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"s1\" desc=\"Test\">\n" +
    "    <dna>ACGTACGT</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(original)
  let rewritten = @src.seqxml_to_xml(doc)
  let doc2 = @src.seqxml_parse(rewritten)
  assert_eq(doc2.n_entries(), 1)
  let entry = doc2.entries()[0]
  assert_eq(entry.id(), "s1")
  assert_eq(entry.description(), "Test")
  assert_eq(entry.sequence(), "ACGTACGT")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "dna")
}

///|
test "seqxml_roundtrip_with_metadata" {
  let original =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <species name=\"Homo sapiens\" ncbiTaxID=\"9606\"/>\n" +
    "  <sourceDB name=\"UniProt\"/>\n" +
    "  <entry id=\"s1\" desc=\"Test\">\n" +
    "    <dna>ACGT</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(original)
  let rewritten = @src.seqxml_to_xml(doc)
  let doc2 = @src.seqxml_parse(rewritten)
  assert_eq(doc2.species_name(), "Homo sapiens")
  assert_eq(doc2.species_taxid(), "9606")
  assert_eq(doc2.source_db(), "UniProt")
  assert_eq(doc2.n_entries(), 1)
  let entry = doc2.entries()[0]
  assert_eq(entry.id(), "s1")
  assert_eq(entry.description(), "Test")
  assert_eq(entry.sequence(), "ACGT")
}

///|
test "seqxml_roundtrip_with_properties" {
  let original =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"p1\" desc=\"With props\">\n" +
    "    <property name=\"gene\" value=\"BRCA1\"/>\n" +
    "    <property name=\"length\" value=\"4\"/>\n" +
    "    <protein>MKLV</protein>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(original)
  let rewritten = @src.seqxml_to_xml(doc)
  let doc2 = @src.seqxml_parse(rewritten)
  assert_eq(doc2.n_entries(), 1)
  let entry = doc2.entries()[0]
  assert_eq(entry.id(), "p1")
  assert_eq(entry.description(), "With props")
  assert_eq(entry.sequence(), "MKLV")
  let props = entry.properties()
  assert_eq(props.length(), 2)
  assert_eq(props[0].0, "gene")
  assert_eq(props[0].1, "BRCA1")
  assert_eq(props[1].0, "length")
  assert_eq(props[1].1, "4")
}

///|
test "seqxml_roundtrip_multiple_entries" {
  let original =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"d1\"><dna>ACGT</dna></entry>\n" +
    "  <entry id=\"r1\"><rna>AUGU</rna></entry>\n" +
    "  <entry id=\"p1\"><protein>MKLV</protein></entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(original)
  let rewritten = @src.seqxml_to_xml(doc)
  let doc2 = @src.seqxml_parse(rewritten)
  assert_eq(doc2.n_entries(), 3)
  let entries = doc2.entries()
  assert_eq(entries[0].id(), "d1")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  assert_eq(entries[0].sequence(), "ACGT")
  assert_eq(entries[1].id(), "r1")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "rna")
  assert_eq(entries[1].sequence(), "AUGU")
  assert_eq(entries[2].id(), "p1")
  assert_eq(@src.SeqXmlType::to_tag(entries[2].seq_type()), "protein")
  assert_eq(entries[2].sequence(), "MKLV")
}

// ============================================================================
// Conversion tests
// ============================================================================

///|
test "seqxml_to_seq_records" {
  let doc = @src.seqxml_sample_document()
  let records = @src.seqxml_to_seq_records(doc)
  assert_eq(records.length(), 3)
  // DNA entry
  assert_eq(records[0].id, "sample_dna")
  assert_eq(records[0].description, "Sample DNA sequence")
  assert_eq(records[0].seq.to_string(), "ATGCGTACGTTGC")
  // RNA entry
  assert_eq(records[1].id, "sample_rna")
  assert_eq(records[1].description, "Sample RNA sequence")
  assert_eq(records[1].seq.to_string(), "AUGCGUACGUUGC")
  // Protein entry
  assert_eq(records[2].id, "sample_protein")
  assert_eq(records[2].description, "Sample protein sequence")
  assert_eq(records[2].seq.to_string(), "MKLVGV")
}

///|
test "seqxml_from_seq_records_auto_detect" {
  let records = [
    @src.SeqRecord::new(
      @src.Seq::new("ACGTACGT"),
      id="r1",
      description="DNA record",
    ),
    @src.SeqRecord::new(
      @src.Seq::new("AUGCAUGC"),
      id="r2",
      description="RNA record",
    ),
    @src.SeqRecord::new(
      @src.Seq::new("MKLVGV"),
      id="r3",
      description="Protein record",
    ),
  ]
  // Default seq_type is UnknownType, so types are auto-detected.
  let doc = @src.seqxml_from_seq_records(records)
  assert_eq(doc.n_entries(), 3)
  let entries = doc.entries()
  // DNA: all chars in {A,C,G,T}
  assert_eq(entries[0].id(), "r1")
  assert_eq(entries[0].description(), "DNA record")
  assert_eq(entries[0].sequence(), "ACGTACGT")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  // RNA: contains U
  assert_eq(entries[1].id(), "r2")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "rna")
  assert_eq(entries[1].sequence(), "AUGCAUGC")
  // Protein: contains M, K, L, V
  assert_eq(entries[2].id(), "r3")
  assert_eq(@src.SeqXmlType::to_tag(entries[2].seq_type()), "protein")
  assert_eq(entries[2].sequence(), "MKLVGV")
}

///|
test "seqxml_from_seq_records_explicit_type" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1"),
    @src.SeqRecord::new(@src.Seq::new("TTTT"), id="r2"),
  ]
  // Force all records to be treated as DNA.
  let doc = @src.seqxml_from_seq_records(
    records,
    seq_type=@src.SeqXmlType::from_tag("dna"),
  )
  assert_eq(doc.n_entries(), 2)
  let entries = doc.entries()
  assert_eq(entries[0].id(), "r1")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  assert_eq(entries[0].sequence(), "ACGT")
  assert_eq(entries[1].id(), "r2")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "dna")
  assert_eq(entries[1].sequence(), "TTTT")
}

///|
test "seqxml_conversion_roundtrip" {
  let doc = @src.seqxml_sample_document()
  let records = @src.seqxml_to_seq_records(doc)
  let doc2 = @src.seqxml_from_seq_records(records)
  assert_eq(doc2.n_entries(), 3)
  let entries = doc2.entries()
  // Verify DNA entry survives the round trip.
  assert_eq(entries[0].id(), "sample_dna")
  assert_eq(entries[0].description(), "Sample DNA sequence")
  assert_eq(entries[0].sequence(), "ATGCGTACGTTGC")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  // Verify RNA entry survives the round trip.
  assert_eq(entries[1].id(), "sample_rna")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "rna")
  assert_eq(entries[1].sequence(), "AUGCGUACGUUGC")
  // Verify protein entry survives the round trip.
  assert_eq(entries[2].id(), "sample_protein")
  assert_eq(@src.SeqXmlType::to_tag(entries[2].seq_type()), "protein")
  assert_eq(entries[2].sequence(), "MKLVGV")
}

// ============================================================================
// Sample document tests
// ============================================================================

///|
test "seqxml_sample_document" {
  let doc = @src.seqxml_sample_document()
  // Three entries.
  assert_eq(doc.n_entries(), 3)
  // Species and source metadata.
  assert_eq(doc.species_name(), "Homo sapiens")
  assert_eq(doc.species_taxid(), "9606")
  assert_eq(doc.source_db(), "BioSeqs")
  let entries = doc.entries()
  // DNA entry.
  assert_eq(entries[0].id(), "sample_dna")
  assert_eq(entries[0].description(), "Sample DNA sequence")
  assert_eq(@src.SeqXmlType::to_tag(entries[0].seq_type()), "dna")
  assert_eq(entries[0].sequence(), "ATGCGTACGTTGC")
  assert_eq(entries[0].properties().length(), 0)
  // RNA entry.
  assert_eq(entries[1].id(), "sample_rna")
  assert_eq(entries[1].description(), "Sample RNA sequence")
  assert_eq(@src.SeqXmlType::to_tag(entries[1].seq_type()), "rna")
  assert_eq(entries[1].sequence(), "AUGCGUACGUUGC")
  assert_eq(entries[1].properties().length(), 0)
  // Protein entry with a property.
  assert_eq(entries[2].id(), "sample_protein")
  assert_eq(entries[2].description(), "Sample protein sequence")
  assert_eq(@src.SeqXmlType::to_tag(entries[2].seq_type()), "protein")
  assert_eq(entries[2].sequence(), "MKLVGV")
  let props = entries[2].properties()
  assert_eq(props.length(), 1)
  assert_eq(props[0].0, "feature")
  assert_eq(props[0].1, "domain")
}

///|
test "seqxml_sample_document_to_xml_roundtrip" {
  let doc = @src.seqxml_sample_document()
  let xml = @src.seqxml_to_xml(doc)
  let doc2 = @src.seqxml_parse(xml)
  assert_eq(doc2.n_entries(), 3)
  assert_eq(doc2.species_name(), "Homo sapiens")
  assert_eq(doc2.species_taxid(), "9606")
  assert_eq(doc2.source_db(), "BioSeqs")
  let entries = doc2.entries()
  assert_eq(entries[0].id(), "sample_dna")
  assert_eq(entries[0].sequence(), "ATGCGTACGTTGC")
  assert_eq(entries[2].id(), "sample_protein")
  assert_eq(entries[2].sequence(), "MKLVGV")
  // The protein entry should retain its property through the round trip.
  let props = entries[2].properties()
  assert_eq(props.length(), 1)
  assert_eq(props[0].0, "feature")
  assert_eq(props[0].1, "domain")
}

// ============================================================================
// Edge cases
// ============================================================================

///|
test "seqxml_parse_empty_document" {
  let doc = @src.seqxml_parse("")
  assert_eq(doc.n_entries(), 0)
  assert_eq(doc.species_name(), "")
  assert_eq(doc.species_taxid(), "")
  assert_eq(doc.source_db(), "")
}

///|
test "seqxml_parse_no_entries" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <species name=\"Mus musculus\" ncbiTaxID=\"10090\"/>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 0)
  // Species header is still parsed.
  assert_eq(doc.species_name(), "Mus musculus")
  assert_eq(doc.species_taxid(), "10090")
}

///|
test "seqxml_parse_entry_without_description" {
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"nodesc\">\n" +
    "    <dna>ACGT</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "nodesc")
  // No desc attribute means description defaults to empty string.
  assert_eq(entry.description(), "")
  assert_eq(entry.sequence(), "ACGT")
}

///|
test "seqxml_parse_self_closing_entry" {
  // Self-closing <entry/> has no sequence element.
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"empty\"/>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  assert_eq(doc.n_entries(), 1)
  let entry = doc.entries()[0]
  assert_eq(entry.id(), "empty")
  assert_eq(entry.description(), "")
  assert_eq(@src.SeqXmlType::to_tag(entry.seq_type()), "unknown")
  assert_eq(entry.sequence(), "")
}

///|
test "seqxml_entity_escaping_in_description" {
  // Parsing should unescape XML entities in the desc attribute.
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"s1\" desc=\"A &amp; B &lt;C&gt;\">\n" +
    "    <dna>ACGT</dna>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  let entry = doc.entries()[0]
  assert_eq(entry.description(), "A & B <C>")
}

///|
test "seqxml_entity_escaping_in_sequence" {
  // Parsing should unescape XML entities in the sequence text content.
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    "  <entry id=\"s1\" desc=\"Test\">\n" +
    "    <protein>AC&amp;T&lt;G</protein>\n" +
    "  </entry>\n" +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  let entry = doc.entries()[0]
  assert_eq(entry.sequence(), "AC&T<G")
}

///|
test "seqxml_entity_escaping_write" {
  // Writing should escape special characters in the description attribute.
  let entry = @src.SeqXmlEntry::new(
    "s1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="A & B <C>",
  )
  let xml = @src.seqxml_entry_to_xml(entry)
  // The ampersand must be escaped as &amp; (and escaped first).
  assert_true(xml.contains("desc=\"A &amp; B &lt;C&gt;\""))
  assert_true(!xml.contains("desc=\"A & B <C>\""))
}

///|
test "seqxml_entity_escaping_roundtrip" {
  // Round-trip: write an entry with special characters, then parse back.
  let entry = @src.SeqXmlEntry::new(
    "esc1",
    @src.SeqXmlType::from_tag("dna"),
    "ACGT",
    description="Quote \" and amp & and lt <",
  )
  let xml = @src.seqxml_entry_to_xml(entry)
  let content =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<seqXML>\n" +
    xml +
    "</seqXML>\n"
  let doc = @src.seqxml_parse(content)
  let parsed = doc.entries()[0]
  assert_eq(parsed.id(), "esc1")
  assert_eq(parsed.description(), "Quote \" and amp & and lt <")
  assert_eq(parsed.sequence(), "ACGT")
}

///|
test "seqxml_empty_document_to_xml" {
  let doc = @src.SeqXmlDocument::new()
  let xml = @src.seqxml_to_xml(doc)
  // Even with no entries, the XML declaration and root element are present.
  assert_true(xml.contains("<?xml version=\"1.0\""))
  assert_true(xml.contains("<seqXML>"))
  assert_true(xml.contains("</seqXML>"))
  // No species or sourceDB when not set.
  assert_true(!xml.contains("<species"))
  assert_true(!xml.contains("<sourceDB"))
}

///|
test "seqxml_parse_entry_id_with_entities" {
  // The id attribute should be unescaped when extracted.
  match @src.seqxml_parse_entry_id("<entry id=\"a&amp;b\">") {
    Some(id) => assert_eq(id, "a&b")
    None => fail("expected Some(id) for entry with escaped id")
  }
}
