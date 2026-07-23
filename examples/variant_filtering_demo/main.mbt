///|

fn main {
  println("=== VariantFiltering Demo - Variant Filtering and Genetic Analysis ===\n")
  
  let info1 = @hashmap.HashMap([("DP", "30"), ("AF", "0.25"), ("AF_POP", "0.001")])
  let info2 = @hashmap.HashMap([("DP", "8"), ("AF", "0.1"), ("AF_POP", "0.05")])
  let info3 = @hashmap.HashMap([("DP", "40"), ("AF", "0.5"), ("AF_POP", "0.0001")])
  let info4 = @hashmap.HashMap([("DP", "25"), ("AF", "0.33"), ("AF_POP", "0.005")])
  
  let variants = [
    @src.Variant::new("chr1", 1000000, "rs1001", "A", "T", 60.0, "PASS", info1, [("proband", "1/1:30"), ("father", "0/1:25"), ("mother", "0/0:30"), ("sibling", "0/0:28")]),
    @src.Variant::new("chr1", 2000000, "rs1002", "G", "C", 20.0, "FAIL", info2, [("proband", "0/1:10"), ("father", "0/1:8"), ("mother", "0/0:20"), ("sibling", "0/0:18")]),
    @src.Variant::new("chr2", 3000000, "rs1003", "C", "G", 55.0, "PASS", info3, [("proband", "1/1:35"), ("father", "0/1:30"), ("mother", "0/1:28"), ("sibling", "0/0:32")]),
    @src.Variant::new("chrX", 4000000, "rs1004", "T", "A", 50.0, "PASS", info4, [("proband_male", "1/0:28"), ("father_male", "0/0:30"), ("mother_female", "0/1:25"), ("sibling_male", "0/0:26")]),
    @src.Variant::new("chr3", 5000000, "rs1005", "A", "C", 45.0, "PASS", @hashmap.HashMap([("DP", "15"), ("AF", "0.15")]), [("proband", "0/1:22"), ("father", "1/1:25"), ("mother", "0/1:20"), ("sibling", "0/1:21")]),
    @src.Variant::new("chr4", 6000000, "rs1006", "G", "T", 70.0, "PASS", @hashmap.HashMap([("DP", "50"), ("AF", "0.4")]), [("proband", "1/1:38"), ("father", "0/1:35"), ("mother", "0/1:32"), ("sibling", "1/1:36")]),
  ]
  
  println("1. Basic Variant Filtering:")
  println("----------------------------")
  
  let params = @src.VariantFilteringParam::new()
  let result = @src.bio_vfilter_filter(variants, params)
  
  println(@src.bio_vfilter_summary(result))
  
  println("\n2. Inheritance Model Testing (Autosomal Dominant):")
  println("---------------------------------------------------")
  
  let phenotypes_ad = @hashmap.HashMap([("proband", true), ("father", false), ("mother", false), ("sibling", false)])
  let model_ad = @src.GeneticModel::new("autosomal_dominant")
  
  let filtered_ad = @src.bio_vfilter_filter_by_inheritance(result.passed_variants, model_ad, phenotypes_ad)
  
  println("Variants matching autosomal dominant model: " + filtered_ad.length().to_string())
  let mut i = 0
  while i < filtered_ad.length() {
    let v = filtered_ad[i]
    println("  " + v.chr + ":" + v.pos.to_string() + " " + v.id + " " + v.ref_allele + ">" + v.alt_allele)
    i = i + 1
  }
  
  println("\n3. Inheritance Model Testing (Autosomal Recessive):")
  println("-----------------------------------------------------")
  
  let phenotypes_ar = @hashmap.HashMap([("proband", true), ("father", false), ("mother", false), ("sibling", false)])
  let model_ar = @src.GeneticModel::new("autosomal_recessive")
  
  let filtered_ar = @src.bio_vfilter_filter_by_inheritance(result.passed_variants, model_ar, phenotypes_ar)
  
  println("Variants matching autosomal recessive model: " + filtered_ar.length().to_string())
  let mut j = 0
  while j < filtered_ar.length() {
    let v = filtered_ar[j]
    println("  " + v.chr + ":" + v.pos.to_string() + " " + v.id + " " + v.ref_allele + ">" + v.alt_allele)
    j = j + 1
  }
  
  println("\n4. Inheritance Model Testing (X-linked Recessive):")
  println("-----------------------------------------------------")
  
  let phenotypes_xr = @hashmap.HashMap([("proband_male", true), ("father_male", false), ("mother_female", false), ("sibling_male", false)])
  let model_xr = @src.GeneticModel::new("x_linked_recessive")
  
  let filtered_xr = @src.bio_vfilter_filter_by_inheritance(result.passed_variants, model_xr, phenotypes_xr)
  
  println("Variants matching X-linked recessive model: " + filtered_xr.length().to_string())
  let mut k = 0
  while k < filtered_xr.length() {
    let v = filtered_xr[k]
    println("  " + v.chr + ":" + v.pos.to_string() + " " + v.id + " " + v.ref_allele + ">" + v.alt_allele)
    k = k + 1
  }
  
  println("\n5. Consequence Prediction:")
  println("----------------------------")
  
  let mut l = 0
  while l < result.passed_variants.length() {
    let v = result.passed_variants[l]
    let consequence = @src.bio_vfilter_predict_consequence(v)
    println(v.id + ": " + v.ref_allele + ">" + v.alt_allele + " -> " + consequence)
    l = l + 1
  }
  
  println("\n=== VariantFiltering Demo Complete ===")
}