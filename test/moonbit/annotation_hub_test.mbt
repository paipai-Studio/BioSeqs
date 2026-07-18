///|
test "annotation_hub_create_and_count" {
  let hub = @src.create_example_hub()
  assert_eq(@src.annotation_hub_count(hub), 8)
}

///|
test "annotation_hub_search_by_type" {
  let hub = @src.create_example_hub()
  let gene_models = @src.annotation_hub_search_by_type(hub, "GeneModel")
  assert_eq(gene_models.length(), 2)
}

///|
test "annotation_hub_search_by_provider" {
  let hub = @src.create_example_hub()
  let ensembl_resources = @src.annotation_hub_search_by_provider(hub, "Ensembl")
  assert_eq(ensembl_resources.length(), 2)
}

///|
test "annotation_hub_search_by_genome" {
  let hub = @src.create_example_hub()
  let human_resources = @src.annotation_hub_search_by_genome(hub, "Homo sapiens")
  assert_eq(human_resources.length(), 6)
}

///|
test "annotation_hub_get_resource" {
  let hub = @src.create_example_hub()
  let resource = @src.annotation_hub_get_resource(hub, "AH12345")
  assert_true(resource.is_some())
  if resource.is_some() {
    assert_eq(resource.unwrap().resource_id, "AH12345")
  }
}

///|
test "annotation_hub_search_query" {
  let hub = @src.create_example_hub()
  let results = @src.annotation_hub_search(hub, "gene")
  assert_eq(results.length(), 2)
}

///|
test "annotation_hub_summarize" {
  let hub = @src.create_example_hub()
  let summary = @src.annotation_hub_summarize(hub)
  assert_true(@src.contains(summary, "Total resources"))
}

///|
test "annotation_hub_resource_to_string" {
  let hub = @src.create_example_hub()
  let resource = @src.annotation_hub_get_resource(hub, "AH12345")
  if resource.is_some() {
    let str = resource.unwrap().to_string()
    assert_true(@src.contains(str, "AH12345"))
  }
}
