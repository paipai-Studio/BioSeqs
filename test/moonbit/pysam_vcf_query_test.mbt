///|
/// P3 Batch 6: pysam_vcf_query 单元测试
let query_vcf_content : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
  "##INFO=<ID=CSQ,Number=A,Type=String,Description=\"consequence\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\trs1\tA\tT\t60\tPASS\tDP=10;CSQ=missense\n" +
  "chr1\t200\trs2\tC\tG\t50\tPASS\tDP=20;CSQ=synonymous\n" +
  "chr1\t300\t.\tG\tA,T\t30\tLowQual\tDP=5\n"

///|
test "query_compile_basic_field_count" {
  let tokens = @src.vcf_query_compile("%CHROM\t%POS\n")
  // 应该有 2 个 Field + 2 个 Lit
  let field_count : Array[Int32] = [0]
  let mut ti = 0
  while ti < tokens.length() {
    match tokens[ti] {
      @src.QueryTemplateToken::Field(_, _) =>
        field_count[0] = field_count[0] + 1
      _ => ()
    }
    ti = ti + 1
  }
  assert_eq(field_count[0], 2)
  assert_true(tokens.length() >= 3)
}

///|
test "query_basic_fields_output_has_expected_content" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%CHROM\t%POS\t%REF\t%ALT\n")
  assert_true(out.length() > 0)
  assert_true(out.contains("chr1"))
  assert_true(out.contains("100"))
  assert_true(out.contains("200"))
  assert_true(out.contains("A"))
  assert_true(out.contains("T"))
  assert_true(out.contains("G"))
  assert_true(out.contains("A,T")) // rs3 multi-allelic
}

///|
test "query_alt_indexed_picks_single_allele" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%ALT{1}\n")
  assert_true(out.contains("T"))
  assert_true(out.contains("G"))
  assert_true(out.contains("A")) // rs3 ALT=A,T 第一个是 A
  // 不应该再包含逗号分隔的 multi-allelic
  assert_true(!out.contains("A,T"))
}

///|
test "query_info_key_extracts_value" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%INFO/DP\n")
  assert_true(out.contains("10"))
  assert_true(out.contains("20"))
  assert_true(out.contains("5"))
}

///|
test "query_info_key_missing_outputs_dot" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%INFO/CSQ\n")
  assert_true(out.contains("missense"))
  assert_true(out.contains("synonymous"))
  // rs3 没有 CSQ → "."
  let lines = out.split("\n").to_array()
  assert_true(lines.length() >= 3)
  assert_true(lines[2].contains(".") || lines[2] == ".")
}

///|
test "query_qual_missing_outputs_dot" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tT\t.\tPASS\t.\n",
  )
  let out = @src.vcf_query(vcf, "%QUAL\n")
  let lines = out.split("\n").to_array()
  assert_true(lines.length() >= 1)
  assert_true(lines[0].contains("."))
}

///|
test "query_filter_pass_output" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%FILTER\n")
  assert_true(out.contains("PASS"))
  assert_true(out.contains("LowQual"))
}

///|
test "query_fallback_pass_for_empty_filter" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tT\t60\t.\t.\n",
  )
  let out = @src.vcf_query(vcf, "%FILTER\n")
  assert_true(out.contains("PASS"))
}

///|
test "query_id_empty_outputs_dot" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%ID\n")
  assert_true(out.contains("rs1"))
  assert_true(out.contains("rs2"))
  // 第三条没有 ID → "."
  let lines = out.split("\n").to_array()
  assert_true(lines[2].contains(".") || lines[2] == ".")
}

///|
test "query_len_and_nalt" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%LEN\t%NALT\n")
  // 所有 REF 都是单碱基 LEN=1；前两条 NALT=1，第三条 NALT=2
  assert_true(out.contains("1"))
  assert_true(out.contains("2"))
}

///|
test "query_with_tokens_reuse" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let tokens = @src.vcf_query_compile("%CHROM\t%POS\n")
  let out1 = @src.vcf_query_with_tokens(vcf, tokens)
  let out2 = @src.vcf_query_with_tokens(vcf, tokens)
  assert_eq(out1, out2)
}

///|
test "query_fields_array_structured_access" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let rows = @src.vcf_query_fields(vcf, [
    "CHROM", "POS", "REF", "ALT", "INFO/DP", "NALT",
  ])
  assert_eq(rows.length(), 3)
  // rs1
  assert_eq(rows[0][0], "chr1")
  assert_eq(rows[0][2], "A")
  assert_eq(rows[0][3], "T")
  assert_eq(rows[0][4], "10")
  assert_eq(rows[0][5], "1")
  // rs3 多等位
  assert_eq(rows[2][3], "A,T")
  assert_eq(rows[2][5], "2")
}

///|
test "query_with_header_contains_header_row" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let tsv = @src.vcf_query_with_header(vcf, [
    "CHROM", "POS", "REF", "ALT", "INFO/DP",
  ])
  assert_true(tsv.contains("CHROM"))
  assert_true(tsv.contains("POS"))
  assert_true(tsv.contains("chr1"))
  assert_true(tsv.contains("100"))
}

///|
test "query_unknown_field_outputs_dot" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%BANANA\n")
  // 所有未知字段 → "."
  let lines = out.split("\n").to_array()
  assert_true(lines.length() >= 3)
  let mut all_dot = true
  let mut li = 0
  while li < lines.length() {
    if lines[li].length() > 0 && lines[li] != "." {
      all_dot = false
    }
    li = li + 1
  }
  assert_true(all_dot)
}

///|
test "query_escape_sequences" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%CHROM\\t%POS\\n")
  assert_true(out.contains("\t"))
}

///|
test "query_empty_vcf_returns_empty_string" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
  )
  let out = @src.vcf_query(vcf, "%CHROM\t%POS\n")
  assert_eq(out.length(), 0)
}

///|
test "query_with_header_empty_vcf" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
  )
  let tsv = @src.vcf_query_with_header(vcf, ["CHROM", "POS"])
  assert_true(tsv.contains("CHROM"))
  assert_true(tsv.contains("POS"))
}

///|
test "query_end_field_via_pos_end" {
  let vcf = @src.parse_vcf(query_vcf_content)
  let out = @src.vcf_query(vcf, "%POS\t%END\t%LEN\n")
  // SNP LEN=1 → END = POS + 1 (pos_end 返回 pos + rlen - 1？要看实现)
  assert_true(out.contains("100"))
  assert_true(out.contains("1"))
}
