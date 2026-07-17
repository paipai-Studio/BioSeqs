///|
/// BLAST parsing demo - demonstrates parsing BLAST tabular and XML output.
///
/// Run with: moon run examples/blast_demo/main.mbt
fn main {
  println("=== BLAST Parsing Demo ===\n")

  println("1. Parsing Tabular Format (-outfmt 6):")
  let tabular_data = "# BLASTN 2.12.0+\n" +
    "# Query: ABC123\n" +
    "# Database: refseq\n" +
    "ABC123\tNP_001112\t98.5\t100\t1\t0\t1\t100\t50\t149\t0.000000000000000000000000000000000000000000000000000000000001\t185.0\n" +
    "ABC123\tNP_004567\t92.0\t80\t6\t1\t5\t84\t10\t89\t0.00000000000000000000000000000000001\t95.5\n" +
    "ABC123\tNP_008901\t85.0\t60\t9\t2\t10\t69\t15\t74\t0.000000000000001\t60.0\n"

  let tabular_record = @bio.parse_blast_tabular(tabular_data)
  if tabular_record.is_some() {
    let record = tabular_record.unwrap()
    println("   Query: " + record.query_id())
    println("   Hits found: " + record.num_hits().to_string())

    let mut i = 0
    while i < record.hits().length() {
      let hit = record.hits()[i]
      println("   Hit " + (i + 1).to_string() + ": " + hit.subject_id())
      let hsp = hit.best_hsp()
      if hsp.is_some() {
        let h = hsp.unwrap()
        println("     Bitscore: " + h.bitscore().to_string())
        println("     E-value: " + h.evalue().to_string())
        println("     Identity: " + h.identity().to_string() + "%")
      }
      i = i + 1
    }
  }

  println("\n2. Filtering by E-value (< 1e-20):")
  if tabular_record.is_some() {
    let record = tabular_record.unwrap()
    let filtered = record.filter_by_evalue(0.00000000000000000001)
    println("   Hits after filtering: " + filtered.num_hits().to_string())
  }

  println("\n3. Filtering by Identity (>= 90%):")
  if tabular_record.is_some() {
    let record = tabular_record.unwrap()
    let filtered = record.filter_by_identity(90.0)
    println("   Hits after filtering: " + filtered.num_hits().to_string())
  }

  println("\n4. Getting Best Hit:")
  if tabular_record.is_some() {
    let record = tabular_record.unwrap()
    let best = record.best_hit()
    if best.is_some() {
      println("   Best hit: " + best.unwrap().subject_id())
    }
  }

  println("\n5. Parsing XML Format:")
  let xml_data = "<BlastOutput>\n" +
    "  <Query_id>test_protein</Query_id>\n" +
    "  <Query_def>Test Protein Sequence</Query_def>\n" +
    "  <Query_len>200</Query_len>\n" +
    "  <Hit>\n" +
    "    <Hit_id>sp|P12345|PROT_HUMAN</Hit_id>\n" +
    "    <Hit_def>Protein homolog in human</Hit_def>\n" +
    "    <Hit_len>250</Hit_len>\n" +
    "    <Hsp>\n" +
    "      <Hsp_query-from>1</Hsp_query-from>\n" +
    "      <Hsp_query-to>100</Hsp_query-to>\n" +
    "      <Hsp_hit-from>50</Hsp_hit-from>\n" +
    "      <Hsp_hit-to>149</Hsp_hit-to>\n" +
    "      <Hsp_identity>95</Hsp_identity>\n" +
    "      <Hsp_positive>98</Hsp_positive>\n" +
    "      <Hsp_gaps>0</Hsp_gaps>\n" +
    "      <Hsp_align-len>100</Hsp_align-len>\n" +
    "      <Hsp_bit-score>190.0</Hsp_bit-score>\n" +
    "      <Hsp_evalue>0.00000000000000000000000000000000000000000000000000000000000000000000000000000000001</Hsp_evalue>\n" +
    "    </Hsp>\n" +
    "  </Hit>\n" +
    "</BlastOutput>\n"

  let xml_record = @bio.parse_blast_xml(xml_data)
  if xml_record.is_some() {
    let record = xml_record.unwrap()
    println("   Query: " + record.query_id())
    println("   Description: " + record.query_def())
    println("   Length: " + record.query_length().to_string())
    println("   Hits: " + record.num_hits().to_string())
  }

  println("\n=== Demo Complete ===")
}
