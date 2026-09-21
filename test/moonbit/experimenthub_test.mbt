///|
/// Tests for ExperimentHub resource metadata and query/filter.

///|
fn ehres(
  title : String,
  species : String,
  taxid : Int,
  genome : String,
  desc : String,
  provider : String,
) -> @src.EhResource {
  @src.EhResource::new(title, species, taxid, genome, desc, provider)
}

///|
fn eh_mk_resources() -> Array[@src.EhResource] {
  [
    ehres("RNA-seq of human cell line", "Homo sapiens", 9606, "hg38", "bulk RNA-seq", "ENCODE"),
    ehres("ChIP-seq peaks mouse", "Mus musculus", 10090, "mm10", "H3K27ac peaks", "ENCODE"),
    ehres("ATAC-seq human PBMC", "Homo sapiens", 9606, "hg38", "single-cell ATAC", "10x Genomics"),
    ehres("ChIP-seq peaks human", "Homo Sapiens", 9606, "hg19", "H3K4me3", "Roadmap"),
  ]
}

///|
test "eh_filter_species_case_insensitive" {
  let rs = eh_mk_resources()
  let idx = @src.eh_filter_species(rs, "homo sapiens")
  assert_eq(idx.length(), 3)
  assert_eq(idx[0], 0)
  assert_eq(idx[1], 2)
  assert_eq(idx[2], 3)
}

///|
test "eh_filter_taxonomy_exact" {
  let rs = eh_mk_resources()
  let idx = @src.eh_filter_taxonomy(rs, 10090)
  assert_eq(idx.length(), 1)
  assert_eq(idx[0], 1)
}

///|
test "eh_search_across_fields" {
  let rs = eh_mk_resources()
  // "mouse" matches species of resource 1
  let mouse = @src.eh_search(rs, "mouse")
  assert_eq(mouse.length(), 1)
  assert_eq(mouse[0], 1)
  // "atac" matches title of resource 2 (case-insensitive)
  let atac = @src.eh_search(rs, "ATAC")
  assert_eq(atac.length(), 1)
  assert_eq(atac[0], 2)
  // "encodE" matches dataprovider of resources 0 and 1
  let encode = @src.eh_search(rs, "encodE")
  assert_eq(encode.length(), 2)
  assert_eq(encode[0], 0)
  assert_eq(encode[1], 1)
}

///|
test "eh_cache_root" {
  assert_eq(@src.eh_cache_root("/home/u"), "/home/u/.cache/ExperimentHub")
}