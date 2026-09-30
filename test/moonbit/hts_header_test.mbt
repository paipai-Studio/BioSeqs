test "hts_header parse basic HD SQ RG PG" {
  let hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621\n@SQ\tSN:chr2\tLN:243199373\n@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.count_lines("HD"), 1)
  assert_eq(h.count_lines("SQ"), 2)
  assert_eq(h.count_lines("RG"), 1)
  assert_eq(h.count_lines("PG"), 1)
  assert_eq(h.nref(), 2)
}

test "hts_header name2tid and tid2name" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr3\tLN:300\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.name2tid("chr1"), 0)
  assert_eq(h.name2tid("chr2"), 1)
  assert_eq(h.name2tid("chr3"), 2)
  assert_eq(h.name2tid("chrX"), -1)
  assert_eq(h.tid2name(0).unwrap_or(""), "chr1")
  assert_eq(h.tid2name(2).unwrap_or(""), "chr3")
  assert_eq(h.tid2name(3), None)
}

test "hts_header tid2len" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.tid2len(0), 100)
  assert_eq(h.tid2len(1), 200)
  assert_eq(h.tid2len(-1), 0)
  assert_eq(h.tid2len(5), 0)
}

test "hts_header to_text roundtrip" {
  let hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
  let h = @src.sam_hrecs_parse(hdr)
  let rebuilt = h.to_text()
  let h2 = @src.sam_hrecs_parse(rebuilt)
  assert_eq(h2.count_lines("HD"), 1)
  assert_eq(h2.count_lines("SQ"), 1)
  assert_eq(h2.name2tid("chr1"), 0)
  assert_eq(h2.tid2len(0), 100)
}

test "hts_header find_line_id" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@RG\tID:rg1\tSM:s1\n"
  let h = @src.sam_hrecs_parse(hdr)
  let line = h.find_line_id("SQ", "SN", "chr1")
  assert_true(line.is_some())
  let line2 = h.find_line_id("SQ", "SN", "chrX")
  assert_true(line2.is_none())
  let line3 = h.find_line_id("RG", "ID", "rg1")
  assert_true(line3.is_some())
}

test "hts_header find_line_pos" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr3\tLN:300\n"
  let h = @src.sam_hrecs_parse(hdr)
  let line0 = h.find_line_pos("SQ", 0)
  assert_true(line0.is_some())
  let line2 = h.find_line_pos("SQ", 2)
  assert_true(line2.is_some())
  let line3 = h.find_line_pos("SQ", 3)
  assert_true(line3.is_none())
}

test "hts_header find_tag_id" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.find_tag_id("SQ", "SN", "chr1", "LN").unwrap_or(""), "100")
  assert_eq(h.find_tag_id("RG", "ID", "rg1", "SM").unwrap_or(""), "sample1")
  assert_eq(h.find_tag_id("RG", "ID", "rg1", "PL").unwrap_or(""), "ILLUMINA")
  assert_eq(h.find_tag_id("RG", "ID", "rgX", "SM"), None)
}

test "hts_header find_tag_pos" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.find_tag_pos("SQ", 0, "SN").unwrap_or(""), "chr1")
  assert_eq(h.find_tag_pos("SQ", 0, "LN").unwrap_or(""), "100")
  assert_eq(h.find_tag_pos("SQ", 1, "SN").unwrap_or(""), "chr2")
  assert_eq(h.find_tag_pos("SQ", 1, "LN").unwrap_or(""), "200")
}

test "hts_header add_line" {
  let hdr = "@HD\tVN:1.6\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.add_line("@SQ\tSN:chr1\tLN:100")
  assert_eq(h2.count_lines("SQ"), 1)
  assert_eq(h2.nref(), 1)
  assert_eq(h2.name2tid("chr1"), 0)
}

test "hts_header add_lines" {
  let hdr = "@HD\tVN:1.6\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.add_lines("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n")
  assert_eq(h2.count_lines("SQ"), 2)
  assert_eq(h2.nref(), 2)
}

test "hts_header remove_line_id" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.remove_line_id("SQ", "SN", "chr1")
  assert_eq(h2.count_lines("SQ"), 1)
  assert_eq(h2.name2tid("chr1"), -1)
  assert_eq(h2.name2tid("chr2"), 0)
}

test "hts_header remove_line_pos" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr3\tLN:300\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.remove_line_pos("SQ", 1)
  assert_eq(h2.count_lines("SQ"), 2)
  assert_eq(h2.name2tid("chr1"), 0)
  assert_eq(h2.name2tid("chr2"), -1)
  assert_eq(h2.name2tid("chr3"), 1)
}

test "hts_header update_line" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.update_line("SQ", "SN", "chr1", "LN", "200")
  assert_eq(h2.tid2len(0), 200)
  let h3 = h.update_line("SQ", "SN", "chr1", "AS", "GRCh38")
  assert_eq(h3.find_tag_id("SQ", "SN", "chr1", "AS").unwrap_or(""), "GRCh38")
}

test "hts_header update_line remove tag" {
  let hdr = "@SQ\tSN:chr1\tLN:100\tAS:GRCh38\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.update_line("SQ", "SN", "chr1", "AS", "")
  assert_eq(h2.find_tag_id("SQ", "SN", "chr1", "AS"), None)
  assert_eq(h2.find_tag_id("SQ", "SN", "chr1", "LN").unwrap_or(""), "100")
}

test "hts_header remove_tag_id" {
  let hdr = "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.remove_tag_id("RG", "ID", "rg1", "PL")
  assert_eq(h2.find_tag_id("RG", "ID", "rg1", "PL"), None)
  assert_eq(h2.find_tag_id("RG", "ID", "rg1", "SM").unwrap_or(""), "sample1")
}

test "hts_header remove_except" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr3\tLN:300\n"
  let h = @src.sam_hrecs_parse(hdr)
  let h2 = h.remove_except("SQ", "SN", "chr2")
  assert_eq(h2.count_lines("SQ"), 1)
  assert_eq(h2.name2tid("chr2"), 0)
  assert_eq(h2.name2tid("chr1"), -1)
}

test "hts_header line_index" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr3\tLN:300\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.line_index("SQ", "SN", "chr1"), 0)
  assert_eq(h.line_index("SQ", "SN", "chr2"), 1)
  assert_eq(h.line_index("SQ", "SN", "chr3"), 2)
  assert_eq(h.line_index("SQ", "SN", "chrX"), -1)
}

test "hts_header lines_of_type" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:s1\n@SQ\tSN:chr2\tLN:200\n"
  let h = @src.sam_hrecs_parse(hdr)
  let sq_lines = h.lines_of_type("SQ")
  assert_eq(sq_lines.length(), 2)
  let rg_lines = h.lines_of_type("RG")
  assert_eq(rg_lines.length(), 1)
}

test "hts_header ref_seq_names" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"
  let h = @src.sam_hrecs_parse(hdr)
  let names = h.ref_seq_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "chr1")
  assert_eq(names[1], "chr2")
}

test "hts_header read_group_ids" {
  let hdr = "@RG\tID:rg1\tSM:s1\n@RG\tID:rg2\tSM:s2\n"
  let h = @src.sam_hrecs_parse(hdr)
  let ids = h.read_group_ids()
  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "rg1")
  assert_eq(ids[1], "rg2")
}

test "hts_header program_ids" {
  let hdr = "@PG\tID:bwa\tPN:bwa\n@PG\tID:samtools\tPN:samtools\n"
  let h = @src.sam_hrecs_parse(hdr)
  let ids = h.program_ids()
  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "bwa")
  assert_eq(ids[1], "samtools")
}

test "hts_header has_read_group has_program" {
  let hdr = "@RG\tID:rg1\tSM:s1\n@PG\tID:bwa\tPN:bwa\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_true(h.has_read_group("rg1"))
  assert_false(h.has_read_group("rg2"))
  assert_true(h.has_program("bwa"))
  assert_false(h.has_program("samtools"))
}

test "hts_header sort_order group_order" {
  let hdr = "@HD\tVN:1.6\tSO:coordinate\tGO:query\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.sort_order(), "coordinate")
  assert_eq(h.group_order(), "query")
}

test "hts_header sort_order unknown" {
  let hdr = "@HD\tVN:1.6\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.sort_order(), "unknown")
  assert_eq(h.group_order(), "none")
}

test "hts_header platform" {
  let hdr = "@RG\tID:rg1\tSM:s1\tPL:ILLUMINA\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.platform("rg1"), "ILLUMINA")
  assert_eq(h.platform("rg2"), "unknown")
}

test "hts_header empty header" {
  let h = @src.sam_hrecs_parse("")
  assert_eq(h.count_lines("HD"), 0)
  assert_eq(h.nref(), 0)
}

test "hts_header CO comment line" {
  let hdr = "@HD\tVN:1.6\n@CO\tThis is a comment\n@SQ\tSN:chr1\tLN:100\n"
  let h = @src.sam_hrecs_parse(hdr)
  assert_eq(h.count_lines("HD"), 1)
  assert_eq(h.count_lines("CO"), 1)
  assert_eq(h.count_lines("SQ"), 1)
  assert_eq(h.nref(), 1)
}

test "hts_header SamHrecLine to_text" {
  let hdr = "@SQ\tSN:chr1\tLN:100\n"
  let h = @src.sam_hrecs_parse(hdr)
  match h.find_line_pos("SQ", 0) {
    Some(line) => assert_eq(line.to_text(), "@SQ\tSN:chr1\tLN:100")
    None => fail("expected SQ line")
  }
}
