///|
test "blast_parse_tabular_basic" {
  let tabular = "query1\tsubject1\t98.5\t100\t1\t0\t1\t100\t1\t100\t0.0001\t150.0\n" +
    "query1\tsubject2\t85.0\t80\t12\t2\t5\t84\t10\t89\t0.001\t80.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  assert_eq(record.query_id(), "query1")
  assert_eq(record.num_hits(), 2)
}

///|
test "blast_parse_tabular_with_comments" {
  let tabular = "# BLASTN 2.12.0+\n" +
    "# Query: my_query\n" +
    "# Database: nr\n" +
    "my_query\tsubject1\t95.0\t50\t2\t0\t1\t50\t10\t59\t0.001\t45.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  assert_eq(record.query_id(), "my_query")
}

///|
test "blast_parse_xml_basic" {
  let xml = "<BlastOutput>\n" +
    "  <Query_id>test_query</Query_id>\n" +
    "  <Query_def>Test Query Protein</Query_def>\n" +
    "  <Query_len>100</Query_len>\n" +
    "  <Hit>\n" +
    "    <Hit_id>subject1</Hit_id>\n" +
    "    <Hit_def>Subject Protein 1</Hit_def>\n" +
    "    <Hit_len>150</Hit_len>\n" +
    "    <Hsp>\n" +
    "      <Hsp_query-from>1</Hsp_query-from>\n" +
    "      <Hsp_query-to>50</Hsp_query-to>\n" +
    "      <Hsp_hit-from>10</Hsp_hit-from>\n" +
    "      <Hsp_hit-to>59</Hsp_hit-to>\n" +
    "      <Hsp_identity>45</Hsp_identity>\n" +
    "      <Hsp_positive>48</Hsp_positive>\n" +
    "      <Hsp_gaps>0</Hsp_gaps>\n" +
    "      <Hsp_align-len>50</Hsp_align-len>\n" +
    "      <Hsp_bit-score>90.5</Hsp_bit-score>\n" +
    "      <Hsp_evalue>0.0000000000000000000000001</Hsp_evalue>\n" +
    "    </Hsp>\n" +
    "  </Hit>\n" +
    "</BlastOutput>\n"
  let record_opt = @src.parse_blast_xml(xml)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  assert_eq(record.query_id(), "test_query")
  assert_eq(record.query_length(), 100)
  assert_eq(record.num_hits(), 1)
}

///|
test "blast_record_filter_by_evalue" {
  let tabular = "query1\tsubject1\t98.0\t100\t2\t0\t1\t100\t1\t100\t0.000001\t150.0\n" +
    "query1\tsubject2\t85.0\t80\t12\t2\t5\t84\t10\t89\t0.001\t80.0\n" +
    "query1\tsubject3\t70.0\t50\t15\t3\t1\t50\t5\t54\t0.1\t50.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  let filtered = record.filter_by_evalue(0.0001)
  assert_eq(filtered.num_hits(), 1)
}

///|
test "blast_record_filter_by_identity" {
  let tabular = "query1\tsubject1\t98.0\t100\t2\t0\t1\t100\t1\t100\t0.000001\t150.0\n" +
    "query1\tsubject2\t85.0\t80\t12\t2\t5\t84\t10\t89\t0.001\t80.0\n" +
    "query1\tsubject3\t70.0\t50\t15\t3\t1\t50\t5\t54\t0.1\t50.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  let filtered = record.filter_by_identity(80.0)
  assert_eq(filtered.num_hits(), 2)
}

///|
test "blast_hit_best_hsp" {
  let tabular = "query1\tsubject1\t98.0\t100\t2\t0\t1\t100\t1\t100\t0.000001\t150.0\n" +
    "query1\tsubject1\t95.0\t80\t4\t1\t5\t84\t5\t84\t0.0001\t120.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  assert_eq(record.num_hits(), 1)
  let hit = record.hits()[0]
  assert_eq(hit.num_hsps(), 2)
  let best_hsp = hit.best_hsp()
  assert_true(best_hsp.is_some())
  assert_eq(best_hsp.unwrap().bitscore(), 150.0)
}

///|
test "blast_record_best_hit" {
  let tabular = "query1\tsubject1\t98.0\t100\t2\t0\t1\t100\t1\t100\t0.000001\t150.0\n" +
    "query1\tsubject2\t99.0\t90\t1\t0\t5\t94\t10\t99\t0.0000001\t160.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  let best = record.best_hit()
  assert_true(best.is_some())
  assert_eq(best.unwrap().subject_id(), "subject2")
}

///|
test "blast_hsp_strand" {
  let tabular = "query1\tsubject1\t90.0\t50\t5\t0\t1\t50\t50\t1\t0.001\t45.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  let hsp = record.hits()[0].hsps()[0]
  assert_eq(hsp.subject_start(), 1)
  assert_eq(hsp.subject_end(), 50)
}

///|
test "blast_record_to_string" {
  let tabular = "query1\tsubject1\t98.0\t100\t2\t0\t1\t100\t1\t100\t0.000001\t150.0\n"
  let record_opt = @src.parse_blast_tabular(tabular)
  assert_true(record_opt.is_some())
  let record = record_opt.unwrap()
  let s = record.to_string()
  assert_true(s.length() > 0)
}

///|
test "blast_empty_input" {
  let record_opt = @src.parse_blast_tabular("")
  assert_true(record_opt.is_none())
}
