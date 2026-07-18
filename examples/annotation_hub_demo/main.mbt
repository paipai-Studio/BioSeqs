///|
fn main {
  println("=== AnnotationHub Demo ===")
  
  let hub = @src.create_example_hub()
  
  println("\n1. Hub Summary:")
  println(@src.annotation_hub_summarize(hub))
  
  println("\n2. Search by type (GeneModel):")
  let gene_models = @src.annotation_hub_search_by_type(hub, "GeneModel")
  let mut i = 0
  while i < gene_models.length() {
    println(gene_models[i].to_string())
    i = i + 1
  }
  
  println("\n3. Search by genome (Mus musculus):")
  let mouse_resources = @src.annotation_hub_search_by_genome(hub, "Mus musculus")
  i = 0
  while i < mouse_resources.length() {
    println(mouse_resources[i].to_string())
    i = i + 1
  }
  
  println("\n4. Search by query ('gene'):")
  let gene_results = @src.annotation_hub_search(hub, "gene")
  i = 0
  while i < gene_results.length() {
    println(gene_results[i].to_string())
    i = i + 1
  }
  
  println("\n5. Get resource by ID (AH12345):")
  let resource = @src.annotation_hub_get_resource(hub, "AH12345")
  if resource.is_some() {
    println(resource.unwrap().to_string())
  }
  
  println("\nAnnotationHub demo completed!")
}