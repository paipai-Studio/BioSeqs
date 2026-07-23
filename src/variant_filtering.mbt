///|

pub struct Variant {
  chr : String
  pos : Int
  id : String
  ref_allele : String
  alt_allele : String
  qual : Double
  filter : String
  info : @hashmap.HashMap[String, String]
  genotypes : Array[(String, String)]
}

pub struct VariantFilteringParam {
  min_qual : Double
  min_depth : Int
  min_allele_freq : Double
  max_missing_rate : Double
  inheritance_mode : String
  min_gq : Int
  max_af_pop : Double
}

pub struct VariantFilterResult {
  passed_variants : Array[Variant]
  filtered_variants : Array[(Variant, String)]
  filter_stats : @hashmap.HashMap[String, Int]
}

pub struct GeneticModel {
  model_type : String
  penetrance : Double
  phenocopy_rate : Double
  allele_frequency : Double
}

pub fn Variant::new(chr : String, pos : Int, id : String, ref_allele : String, 
                    alt_allele : String, qual : Double, filter : String, 
                    info : @hashmap.HashMap[String, String], genotypes : Array[(String, String)]) -> Variant {
  Variant::{ chr, pos, id, ref_allele, alt_allele, qual, filter, info, genotypes }
}

pub fn VariantFilteringParam::new() -> VariantFilteringParam {
  VariantFilteringParam::{ 
    min_qual: 30.0, 
    min_depth: 10, 
    min_allele_freq: 0.01, 
    max_missing_rate: 0.2,
    inheritance_mode: "autosomal_dominant",
    min_gq: 20,
    max_af_pop: 0.01
  }
}

pub fn VariantFilterResult::new() -> VariantFilterResult {
  VariantFilterResult::{ 
    passed_variants: Array::new(), 
    filtered_variants: Array::new(), 
    filter_stats: @hashmap.HashMap([]) 
  }
}

pub fn GeneticModel::new(model_type : String) -> GeneticModel {
  let penetrance = 
    if model_type == "autosomal_dominant" { 0.95 }
    else if model_type == "autosomal_recessive" { 0.95 }
    else if model_type == "x_linked_dominant" { 0.95 }
    else if model_type == "x_linked_recessive" { 0.95 }
    else if model_type == "compound_heterozygous" { 0.95 }
    else { 0.95 }
  
  GeneticModel::{ model_type, penetrance, phenocopy_rate: 0.01, allele_frequency: 0.001 }
}

pub fn bio_vfilter_filter(variants : Array[Variant], params : VariantFilteringParam) -> VariantFilterResult {
  let mut result = VariantFilterResult::new()
  
  let len = variants.length()
  let mut i = 0
  while i < len {
    let variant = variants[i]
    let mut filter_reasons : Array[String] = Array::new()
    
    if variant.qual < params.min_qual {
      filter_reasons = push_string(filter_reasons, "QUAL below threshold")
    }
    
    if variant.filter != "." && variant.filter != "PASS" {
      filter_reasons = push_string(filter_reasons, "FILTER not PASS")
    }
    
    let dp = parse_info_int(variant.info, "DP")
    if dp < params.min_depth {
      filter_reasons = push_string(filter_reasons, "DP below threshold")
    }
    
    let af = parse_info_double(variant.info, "AF")
    if af > 0.0 && af < params.min_allele_freq {
      filter_reasons = push_string(filter_reasons, "AF below threshold")
    }
    
    let gq_pass = check_genotype_quality(variant.genotypes, params.min_gq)
    if !gq_pass {
      filter_reasons = push_string(filter_reasons, "GQ below threshold")
    }
    
    let missing_rate = calculate_missing_rate(variant.genotypes)
    if missing_rate > params.max_missing_rate {
      filter_reasons = push_string(filter_reasons, "Missing rate above threshold")
    }
    
    if params.max_af_pop > 0.0 {
      let af_pop = parse_info_double(variant.info, "AF_POP")
      if af_pop > params.max_af_pop {
        filter_reasons = push_string(filter_reasons, "Population AF above threshold")
      }
    }
    
    if filter_reasons.length() == 0 {
      result = push_variant(result, variant)
    } else {
      let fr_len = filter_reasons.length()
      let mut j = 0
      while j < fr_len {
        let reason = filter_reasons[j]
        let current = if result.filter_stats.contains(reason) { result.filter_stats.get(reason).unwrap() } else { 0 }
        result.filter_stats.set(reason, current + 1)
        j = j + 1
      }
      let reason_str = join_strings(filter_reasons, "; ")
      result = push_filtered_variant(result, (variant, reason_str))
    }
    
    i = i + 1
  }
  
  result
}

fn push_string(arr : Array[String], s : String) -> Array[String] {
  let len = arr.length()
  let new_arr : Array[String] = Array::new()
  let mut i = 0
  while i < len {
    new_arr.push(arr[i])
    i = i + 1
  }
  new_arr.push(s)
  new_arr
}

fn push_variant(result : VariantFilterResult, v : Variant) -> VariantFilterResult {
  let len = result.passed_variants.length()
  let new_passed : Array[Variant] = Array::new()
  let mut i = 0
  while i < len {
    new_passed.push(result.passed_variants[i])
    i = i + 1
  }
  new_passed.push(v)
  VariantFilterResult::{ ..result, passed_variants: new_passed }
}

fn push_filtered_variant(result : VariantFilterResult, v : (Variant, String)) -> VariantFilterResult {
  let len = result.filtered_variants.length()
  let new_filtered : Array[(Variant, String)] = Array::new()
  let mut i = 0
  while i < len {
    new_filtered.push(result.filtered_variants[i])
    i = i + 1
  }
  new_filtered.push(v)
  VariantFilterResult::{ ..result, filtered_variants: new_filtered }
}

fn join_strings(arr : Array[String], sep : String) -> String {
  let len = arr.length()
  if len == 0 {
    return ""
  }
  let mut result = arr[0]
  let mut i = 1
  while i < len {
    result = result + sep + arr[i]
    i = i + 1
  }
  result
}

fn parse_info_int(info : @hashmap.HashMap[String, String], key : String) -> Int {
  if info.contains(key) {
    parse_int(info.get(key).unwrap())
  } else {
    0
  }
}

fn parse_info_double(info : @hashmap.HashMap[String, String], key : String) -> Double {
  if info.contains(key) {
    let parsed = parse_double(info.get(key).unwrap())
    if parsed.is_some() { parsed.unwrap() } else { 0.0 }
  } else {
    0.0
  }
}

fn check_genotype_quality(genotypes : Array[(String, String)], min_gq : Int) -> Bool {
  let mut pass_count = 0
  let mut total = 0
  
  let len = genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = genotypes[i]
    if gt != "./." {
      total = total + 1
      
      let parts = split_by_char(gt, 58)
      if parts.length() >= 3 {
        let gq = parse_int(parts[2])
        if gq >= min_gq {
          pass_count = pass_count + 1
        }
      } else {
        pass_count = pass_count + 1
      }
    }
    i = i + 1
  }
  
  total == 0 || pass_count == total
}

fn calculate_missing_rate(genotypes : Array[(String, String)]) -> Double {
  let len = genotypes.length()
  if len == 0 {
    return 1.0
  }
  
  let mut missing = 0
  let mut i = 0
  while i < len {
    let (_, gt) = genotypes[i]
    if gt == "./." {
      missing = missing + 1
    }
    i = i + 1
  }
  
  missing.to_double() / len.to_double()
}

pub fn bio_vfilter_check_inheritance(variant : Variant, model : GeneticModel, 
                                      phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  if model.model_type == "autosomal_dominant" { 
    check_autosomal_dominant(variant, phenotypes) 
  } else if model.model_type == "autosomal_recessive" { 
    check_autosomal_recessive(variant, phenotypes) 
  } else if model.model_type == "x_linked_dominant" { 
    check_x_linked_dominant(variant, phenotypes) 
  } else if model.model_type == "x_linked_recessive" { 
    check_x_linked_recessive(variant, phenotypes) 
  } else if model.model_type == "compound_heterozygous" { 
    check_compound_heterozygous(variant, phenotypes) 
  } else { 
    true 
  }
}

fn check_autosomal_dominant(variant : Variant, phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  let mut affected_carriers = 0
  let mut unaffected_carriers = 0
  let mut affected_non_carriers = 0
  
  let len = variant.genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = variant.genotypes[i]
    let is_affected = if phenotypes.contains(sample) { phenotypes.get(sample).unwrap() } else { false }
    let is_carrier = is_carrier_gt(gt)
    
    if is_affected && is_carrier {
      affected_carriers = affected_carriers + 1
    } else if !is_affected && is_carrier {
      unaffected_carriers = unaffected_carriers + 1
    } else if is_affected && !is_carrier {
      affected_non_carriers = affected_non_carriers + 1
    }
    
    i = i + 1
  }
  
  affected_non_carriers == 0 && affected_carriers >= 1
}

fn check_autosomal_recessive(variant : Variant, phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  let mut affected_homozygous = 0
  let mut unaffected_homozygous = 0
  let mut affected_heterozygous = 0
  
  let len = variant.genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = variant.genotypes[i]
    let is_affected = if phenotypes.contains(sample) { phenotypes.get(sample).unwrap() } else { false }
    let is_homozygous = is_homozygous_gt(gt)
    let is_heterozygous = is_heterozygous_gt(gt)
    
    if is_affected && is_homozygous {
      affected_homozygous = affected_homozygous + 1
    } else if !is_affected && is_homozygous {
      unaffected_homozygous = unaffected_homozygous + 1
    } else if is_affected && is_heterozygous {
      affected_heterozygous = affected_heterozygous + 1
    }
    
    i = i + 1
  }
  
  affected_heterozygous == 0 && unaffected_homozygous == 0 && affected_homozygous >= 1
}

fn check_x_linked_dominant(variant : Variant, phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  if !is_x_chromosome(variant.chr) {
    return false
  }
  
  let mut affected_males_with_mut = 0
  let mut unaffected_males_with_mut = 0
  let mut affected_females_with_mut = 0
  let mut unaffected_females_with_mut = 0
  
  let len = variant.genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = variant.genotypes[i]
    let is_affected = if phenotypes.contains(sample) { phenotypes.get(sample).unwrap() } else { false }
    let is_carrier = is_carrier_gt(gt)
    
    if is_male(sample) {
      if is_affected && is_carrier {
        affected_males_with_mut = affected_males_with_mut + 1
      } else if !is_affected && is_carrier {
        unaffected_males_with_mut = unaffected_males_with_mut + 1
      }
    } else {
      if is_affected && is_carrier {
        affected_females_with_mut = affected_females_with_mut + 1
      } else if !is_affected && is_carrier {
        unaffected_females_with_mut = unaffected_females_with_mut + 1
      }
    }
    
    i = i + 1
  }
  
  unaffected_males_with_mut == 0 && affected_males_with_mut >= 1
}

fn check_x_linked_recessive(variant : Variant, phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  if !is_x_chromosome(variant.chr) {
    return false
  }
  
  let mut affected_males_with_mut = 0
  let mut unaffected_males_with_mut = 0
  let mut affected_females_homozygous = 0
  let mut unaffected_females_homozygous = 0
  
  let len = variant.genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = variant.genotypes[i]
    let is_affected = if phenotypes.contains(sample) { phenotypes.get(sample).unwrap() } else { false }
    let is_carrier = is_carrier_gt(gt)
    let is_homozygous = is_homozygous_gt(gt)
    
    if is_male(sample) {
      if is_affected && is_carrier {
        affected_males_with_mut = affected_males_with_mut + 1
      } else if !is_affected && is_carrier {
        unaffected_males_with_mut = unaffected_males_with_mut + 1
      }
    } else {
      if is_affected && is_homozygous {
        affected_females_homozygous = affected_females_homozygous + 1
      } else if !is_affected && is_homozygous {
        unaffected_females_homozygous = unaffected_females_homozygous + 1
      }
    }
    
    i = i + 1
  }
  
  unaffected_males_with_mut == 0 && unaffected_females_homozygous == 0 && (affected_males_with_mut >= 1 || affected_females_homozygous >= 1)
}

fn check_compound_heterozygous(variant : Variant, phenotypes : @hashmap.HashMap[String, Bool]) -> Bool {
  let mut affected_heterozygous = 0
  let mut unaffected_heterozygous = 0
  
  let len = variant.genotypes.length()
  let mut i = 0
  while i < len {
    let (sample, gt) = variant.genotypes[i]
    let is_affected = if phenotypes.contains(sample) { phenotypes.get(sample).unwrap() } else { false }
    let is_heterozygous = is_heterozygous_gt(gt)
    
    if is_affected && is_heterozygous {
      affected_heterozygous = affected_heterozygous + 1
    } else if !is_affected && is_heterozygous {
      unaffected_heterozygous = unaffected_heterozygous + 1
    }
    
    i = i + 1
  }
  
  affected_heterozygous >= 2 && unaffected_heterozygous > 0
}

fn is_carrier_gt(gt : String) -> Bool {
  gt == "0/1" || gt == "1/0" || gt == "1/1" || gt == "0|1" || gt == "1|0" || gt == "1|1"
}

fn is_homozygous_gt(gt : String) -> Bool {
  gt == "1/1" || gt == "1|1"
}

fn is_heterozygous_gt(gt : String) -> Bool {
  gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0"
}

fn is_x_chromosome(chr : String) -> Bool {
  chr == "X" || chr == "x" || chr == "chrX" || chr == "chrX" || chr == "23"
}

fn is_male(sample : String) -> Bool {
  let lower_sample = to_lowercase(sample)
  lower_sample.ends_with("_male") || contains_substring(lower_sample, "father") || contains_substring(lower_sample, "son")
}

fn contains_substring(s : String, substr : String) -> Bool {
  let s_len = s.length()
  let sub_len = substr.length()
  
  if sub_len == 0 || s_len < sub_len {
    return false
  }
  
  let mut i = 0
  while i <= s_len - sub_len {
    let current = substring(s, i, sub_len)
    if current == substr {
      return true
    }
    i = i + 1
  }
  
  false
}

pub fn bio_vfilter_filter_by_inheritance(variants : Array[Variant], model : GeneticModel, 
                                          phenotypes : @hashmap.HashMap[String, Bool]) -> Array[Variant] {
  let len = variants.length()
  let mut result : Array[Variant] = Array::new()
  let mut i = 0
  while i < len {
    let v = variants[i]
    if bio_vfilter_check_inheritance(v, model, phenotypes) {
      result = push_variant_array(result, v)
    }
    i = i + 1
  }
  result
}

fn push_variant_array(arr : Array[Variant], v : Variant) -> Array[Variant] {
  let len = arr.length()
  let new_arr : Array[Variant] = Array::new()
  let mut i = 0
  while i < len {
    new_arr.push(arr[i])
    i = i + 1
  }
  new_arr.push(v)
  new_arr
}

pub fn bio_vfilter_summary(result : VariantFilterResult) -> String {
  let mut summary = "\nVariant Filtering Summary\n"
  summary = summary + "=".repeat(60) + "\n"
  
  let total = result.passed_variants.length() + result.filtered_variants.length()
  
  summary = summary + "Total variants: " + total.to_string() + "\n"
  summary = summary + "Passed: " + result.passed_variants.length().to_string() + "\n"
  summary = summary + "Filtered: " + result.filtered_variants.length().to_string() + "\n"
  
  summary = summary + "\nFilter statistics:\n"
  summary = summary + "-".repeat(60) + "\n"
  
  let entries = result.filter_stats.iter().collect()
  let entries_len = entries.length()
  
  let mut sorted : Array[(String, Int)] = Array::new()
  let mut i = 0
  while i < entries_len {
    sorted.push(entries[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < entries_len {
    let mut k = j + 1
    while k < entries_len {
      if sorted[k].1 > sorted[j].1 {
        let temp = sorted[j]
        let new_arr : Array[(String, Int)] = Array::new()
        let mut m = 0
        while m < entries_len {
          if m == j {
            new_arr.push(sorted[k])
          } else if m == k {
            new_arr.push(temp)
          } else {
            new_arr.push(sorted[m])
          }
          m = m + 1
        }
        sorted = new_arr
      }
      k = k + 1
    }
    j = j + 1
  }
  
  let mut l = 0
  while l < sorted.length() {
    let (reason, count) = sorted[l]
    summary = summary + reason + ": " + count.to_string() + "\n"
    l = l + 1
  }
  
  summary
}

pub fn bio_vfilter_predict_consequence(variant : Variant) -> String {
  let len_diff = variant.alt_allele.length() - variant.ref_allele.length()
  
  if len_diff == 0 {
    "missense_variant"
  } else if len_diff > 0 {
    "insertion"
  } else {
    "deletion"
  }
}