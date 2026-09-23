///|
/// Tests for the KEGG REST API client (`kegg_rest`): URL construction and
/// argument validation via an injected mock transport.

///|
fn mk_client(urls : Ref[Array[String]]) -> @src.KeggRestClient {
  @src.KeggRestClient::new(
    fetch=fn(url : String) -> String raise @src.KeggRestError {
      urls.val.push(url)
      "MOCK"
    },
    now=fn() { 0.0 },
  )
}

///|
fn call_urls(
  f : (@src.KeggRestClient) -> String raise @src.KeggRestError,
) -> Array[String] raise @src.KeggRestError {
  let urls : Ref[Array[String]] = Ref([])
  let c = mk_client(urls)
  ignore(f(c))
  urls.val
}

///|
fn rejects(f : (@src.KeggRestClient) -> String raise @src.KeggRestError) -> Bool {
  try {
    let urls : Ref[Array[String]] = Ref([])
    let c = mk_client(urls)
    ignore(f(c))
    false
  } catch {
    _ => true
  }
}

///|
test "kegg_rest_url" {
  assert_eq(@src.kegg_rest_url("info", "pathway"), "https://rest.kegg.jp/info/pathway")
  assert_eq(
    @src.kegg_rest_url("list", "pathway", arg2="hsa"),
    "https://rest.kegg.jp/list/pathway/hsa",
  )
  assert_eq(
    @src.kegg_rest_url("get", "hsa:00001", arg2="kgml"),
    "https://rest.kegg.jp/get/hsa:00001/kgml",
  )
}

///|
test "kegg_rest_info_and_list" {
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_info(c, "pathway") })[0],
    "https://rest.kegg.jp/info/pathway",
  )
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_list(c, "pathway") })[0],
    "https://rest.kegg.jp/list/pathway",
  )
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_list(c, "pathway", org="hsa") })[0],
    "https://rest.kegg.jp/list/pathway/hsa",
  )
  // 非 pathway/module 数据库带 org 报错
  assert_true(rejects(fn(c) { @src.kegg_rest_list(c, "genes", org="hsa") }))
  // 多条目 "+" 连接
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_list_entries(c, ["hsa:10458", "hsa:10643"]) })[0],
    "https://rest.kegg.jp/list/hsa:10458+hsa:10643",
  )
}

///|
test "kegg_rest_find" {
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_find(c, "genes", "alcohol") })[0],
    "https://rest.kegg.jp/find/genes/alcohol",
  )
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_find(c, "compound", "C7H10O5", option="formula") })[0],
    "https://rest.kegg.jp/find/compound/C7H10O5/formula",
  )
  // compound/drug 之外的数据库不能带 option
  assert_true(rejects(fn(c) { @src.kegg_rest_find(c, "genes", "x", option="formula") }))
  // 多查询 "+" 连接
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_find_multi(c, "genes", ["a", "b"]) })[0],
    "https://rest.kegg.jp/find/genes/a+b",
  )
}

///|
test "kegg_rest_get" {
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_get(c, "hsa:10458+ece:Z5100") })[0],
    "https://rest.kegg.jp/get/hsa:10458+ece:Z5100",
  )
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_get(c, "hsa:10458", option="ntseq") })[0],
    "https://rest.kegg.jp/get/hsa:10458/ntseq",
  )
  // 非法 option 报错
  assert_true(rejects(fn(c) { @src.kegg_rest_get(c, "hsa:10458", option="bogus") }))
  // > 10 条 get 报错
  let many : Array[String] = []
  for i = 0; i < 11; i = i + 1 {
    many.push("hsa:1" + i.to_string())
  }
  assert_true(rejects(fn(c) { @src.kegg_rest_get_multi(c, many) }))
}

///|
test "kegg_rest_conv_and_link" {
  // gene id 转换: uniprot 是合法的 gene 数据库
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_conv(c, "uniprot", "hsa") })[0],
    "https://rest.kegg.jp/conv/uniprot/hsa",
  )
  // 非法 pair 报错
  assert_true(rejects(fn(c) { @src.kegg_rest_conv(c, "hsa", "mmu") }))
  // 化学 substance 转换: pubchem × compound
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_conv(c, "pubchem", "compound") })[0],
    "https://rest.kegg.jp/conv/pubchem/compound",
  )
  // 非法 option 报错
  assert_true(rejects(fn(c) { @src.kegg_rest_conv(c, "uniprot", "hsa", option="json") }))
  // link
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_link(c, "pathway", "hsa:10458") })[0],
    "https://rest.kegg.jp/link/pathway/hsa:10458",
  )
  // conv_multi / link_multi 多条目 "+"
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_conv_multi(c, "uniprot", ["hsa:10458", "ece:Z5100"]) })[0],
    "https://rest.kegg.jp/conv/uniprot/hsa:10458+ece:Z5100",
  )
  assert_eq(
    call_urls(fn(c) { @src.kegg_rest_link_multi(c, "pathway", ["hsa:10458", "hsa:10643"]) })[0],
    "https://rest.kegg.jp/link/pathway/hsa:10458+hsa:10643",
  )
}