///|
/// Tests for the INSDC IO module.

///|
test "insdc_parse_feature_table_cds" {
  let lines = [
    "     CDS             1..100",
    "                     /gene=\"testGene\"",
    "                     /product=\"test protein\"",
    "                     /translation=\"MKVL\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  assert_eq(features.length(), 1)
  assert_eq(features[0].feature_type, "CDS")
  assert_eq(features[0].location_str, "1..100")
  let gene = @src.extract_gene_name(features[0])
  assert_eq(gene, Some("testGene"))
  let trans = @src.extract_translation(features[0])
  assert_eq(trans, Some("MKVL"))
  let prod = @src.extract_feature_qualifier(features[0], "product")
  assert_eq(prod, Some("test protein"))
}

///|
test "insdc_parse_feature_table_gene" {
  let lines = [
    "     gene            1..100",
    "                     /gene=\"testGene\"",
    "                     /locus_tag=\"TEST001\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  assert_eq(features.length(), 1)
  assert_eq(features[0].feature_type, "gene")
  assert_eq(features[0].location_str, "1..100")
  let gene = @src.extract_gene_name(features[0])
  assert_eq(gene, Some("testGene"))
  let locus = @src.extract_feature_qualifier(features[0], "locus_tag")
  assert_eq(locus, Some("TEST001"))
}

///|
test "insdc_parse_feature_table_empty" {
  let features = @src.parse_insdc_feature_table([])
  assert_eq(features.length(), 0)
}

///|
test "insdc_parse_feature_table_multiline_qualifier" {
  let lines = [
    "     CDS             1..100",
    "                     /translation=\"MKVL",
    "                     KLMN\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  assert_eq(features.length(), 1)
  let trans = @src.extract_translation(features[0])
  assert_eq(trans, Some("MKVLKLMN"))
}

///|
test "insdc_parse_feature_table_multiple_features" {
  let lines = [
    "     CDS             1..100",
    "                     /gene=\"geneA\"",
    "     gene            1..100",
    "                     /gene=\"geneA\"",
    "                     /locus_tag=\"LA001\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  assert_eq(features.length(), 2)
  assert_eq(features[0].feature_type, "CDS")
  assert_eq(features[1].feature_type, "gene")
  let gene0 = @src.extract_gene_name(features[0])
  assert_eq(gene0, Some("geneA"))
  let gene1 = @src.extract_gene_name(features[1])
  assert_eq(gene1, Some("geneA"))
}

///|
test "insdc_parse_location_simple_range" {
  let result = @src.parse_insdc_location("123..456")
  assert_eq(result, "123..456")
}

///|
test "insdc_parse_location_complement" {
  let result = @src.parse_insdc_location("complement(123..456)")
  assert_eq(result, "complement(123..456)")
}

///|
test "insdc_parse_location_join" {
  let result = @src.parse_insdc_location("join(100..200,300..400)")
  assert_eq(result, "join(100..200,300..400)")
}

///|
test "insdc_parse_location_with_whitespace" {
  let result = @src.parse_insdc_location("join(100..200, 300..400)")
  assert_eq(result, "join(100..200,300..400)")
}

///|
test "insdc_extract_feature_qualifier_missing" {
  let lines = [
    "     CDS             1..100",
    "                     /gene=\"testGene\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  let missing = @src.extract_feature_qualifier(features[0], "product")
  assert_eq(missing, None)
}

///|
test "insdc_extract_translation_no_cds" {
  let lines = [
    "     gene            1..100",
    "                     /gene=\"testGene\"",
  ]
  let features = @src.parse_insdc_feature_table(lines)
  let trans = @src.extract_translation(features[0])
  assert_eq(trans, None)
}

///|
test "insdc_parse_genbank_record_locus" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.locus_name, "NM_001")
  assert_eq(record.sequence_length, 100)
  assert_eq(record.molecule_type, "DNA")
  assert_eq(record.topology, "linear")
  assert_eq(record.division, "PRI")
  assert_eq(record.date, "01-JAN-2024")
}

///|
test "insdc_parse_genbank_record_definition" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nDEFINITION  Test gene for unit testing.\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.definition, "Test gene for unit testing")
}

///|
test "insdc_parse_genbank_record_accession_version" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nACCESSION   NM_001\nVERSION     NM_001.1\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.accession, "NM_001")
  assert_eq(record.version, "NM_001.1")
}

///|
test "insdc_parse_genbank_record_features" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nFEATURES             Location/Qualifiers\n     CDS             1..100\n                     /gene=\"testGene\"\n                     /translation=\"MKVL\"\n     gene            1..100\n                     /gene=\"testGene\"\n                     /locus_tag=\"TEST001\"\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.features.length(), 2)
  assert_eq(record.features[0].feature_type, "CDS")
  assert_eq(record.features[1].feature_type, "gene")
  let trans = @src.extract_translation(record.features[0])
  assert_eq(trans, Some("MKVL"))
}

///|
test "insdc_parse_genbank_record_origin_sequence" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.sequence, "atgaaagtgctttga")
}

///|
test "insdc_parse_genbank_record_empty" {
  let record = @src.parse_genbank_record("")
  assert_eq(record.locus_name, "")
  assert_eq(record.sequence_length, 0)
  assert_eq(record.sequence, "")
  assert_eq(record.features.length(), 0)
}

///|
test "insdc_parse_genbank_record_full" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nDEFINITION  Test gene for unit testing.\nACCESSION   NM_001\nVERSION     NM_001.1\nSOURCE      Homo sapiens\n  ORGANISM  Homo sapiens\n            Eukaryota; Metazoa.\nFEATURES             Location/Qualifiers\n     CDS             1..100\n                     /gene=\"testGene\"\n                     /product=\"test protein\"\n                     /translation=\"MKVL\"\n     gene            1..100\n                     /gene=\"testGene\"\n                     /locus_tag=\"TEST001\"\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  assert_eq(record.locus_name, "NM_001")
  assert_eq(record.sequence_length, 100)
  assert_eq(record.molecule_type, "DNA")
  assert_eq(record.topology, "linear")
  assert_eq(record.division, "PRI")
  assert_eq(record.date, "01-JAN-2024")
  assert_eq(record.definition, "Test gene for unit testing")
  assert_eq(record.accession, "NM_001")
  assert_eq(record.version, "NM_001.1")
  assert_eq(record.features.length(), 2)
  assert_eq(record.sequence, "atgaaagtgctttga")
}

///|
test "insdc_write_genbank_record_round_trip" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nDEFINITION  Test gene for unit testing.\nACCESSION   NM_001\nVERSION     NM_001.1\nFEATURES             Location/Qualifiers\n     CDS             1..100\n                     /gene=\"testGene\"\n                     /translation=\"MKVL\"\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  let output = @src.write_genbank_record(record)
  assert_true(output.contains("LOCUS"))
  assert_true(output.contains("NM_001"))
  assert_true(output.contains("DEFINITION"))
  assert_true(output.contains("Test gene for unit testing"))
  assert_true(output.contains("ACCESSION"))
  assert_true(output.contains("VERSION"))
  assert_true(output.contains("NM_001.1"))
  assert_true(output.contains("FEATURES"))
  assert_true(output.contains("CDS"))
  assert_true(output.contains("testGene"))
  assert_true(output.contains("ORIGIN"))
  assert_true(output.contains("//"))
}

///|
test "insdc_write_genbank_record_contains_sequence" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  let output = @src.write_genbank_record(record)
  assert_true(output.contains("atgaaagtgc"))
}

///|
test "insdc_parse_embl_record_id" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.locus_name, "NM_001")
  assert_eq(record.sequence_length, 100)
  assert_eq(record.molecule_type, "DNA")
  assert_eq(record.topology, "linear")
  assert_eq(record.division, "PRI")
}

///|
test "insdc_parse_embl_record_accession" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\nAC   NM_001;\nXX\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.accession, "NM_001")
  assert_eq(record.version, "NM_001.1")
}

///|
test "insdc_parse_embl_record_definition" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\nDE   Test gene for unit testing.\nXX\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.definition, "Test gene for unit testing")
}

///|
test "insdc_parse_embl_record_organism" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\nOS   Homo sapiens\nXX\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.organism, "Homo sapiens")
}

///|
test "insdc_parse_embl_record_sequence" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\nSQ   Sequence 100 BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n     atgaaagtgc tttga\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.sequence, "atgaaagtgctttga")
}

///|
test "insdc_parse_embl_record_full" {
  let content = "ID   NM_001; SV 1; linear; DNA; STD; PRI; 100 BP.\nXX\nAC   NM_001;\nXX\nDE   Test gene for unit testing.\nXX\nOS   Homo sapiens\nXX\nSQ   Sequence 100 BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n     atgaaagtgc tttga\n//\n"
  let record = @src.parse_embl_record(content)
  assert_eq(record.locus_name, "NM_001")
  assert_eq(record.sequence_length, 100)
  assert_eq(record.molecule_type, "DNA")
  assert_eq(record.topology, "linear")
  assert_eq(record.division, "PRI")
  assert_eq(record.accession, "NM_001")
  assert_eq(record.version, "NM_001.1")
  assert_eq(record.definition, "Test gene for unit testing")
  assert_eq(record.organism, "Homo sapiens")
  assert_eq(record.sequence, "atgaaagtgctttga")
}

///|
test "insdc_write_embl_record_round_trip" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nDEFINITION  Test gene for unit testing.\nACCESSION   NM_001\nVERSION     NM_001.1\nFEATURES             Location/Qualifiers\n     CDS             1..100\n                     /gene=\"testGene\"\n                     /translation=\"MKVL\"\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  let output = @src.write_embl_record(record)
  assert_true(output.contains("ID   "))
  assert_true(output.contains("NM_001"))
  assert_true(output.contains("AC   "))
  assert_true(output.contains("DE   "))
  assert_true(output.contains("Test gene for unit testing"))
  assert_true(output.contains("SQ   "))
  assert_true(output.contains("FT   "))
  assert_true(output.contains("CDS"))
  assert_true(output.contains("testGene"))
  assert_true(output.contains("//"))
}

///|
test "insdc_write_embl_record_sequence" {
  let content = "LOCUS       NM_001                 100 bp    DNA     linear   PRI 01-JAN-2024\nORIGIN\n        1 atgaaagtgc tttga\n//\n"
  let record = @src.parse_genbank_record(content)
  let output = @src.write_embl_record(record)
  assert_true(output.contains("atgaaagtgc"))
}
