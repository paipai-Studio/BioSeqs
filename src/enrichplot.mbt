///|

pub struct EnrichTerm {
  term_id : String
  term_name : String
  pvalue : Double
  adjusted_pvalue : Double
  enrichment_score : Double
  nes : Double
  gene_count : Int
  genes : Array[String]
}

pub struct EnrichResult {
  terms : Array[EnrichTerm]
  gene_set_name : String
  organism : String
}

pub fn EnrichTerm::new(term_id : String, term_name : String, pvalue : Double, adjusted_pvalue : Double, 
                       enrichment_score : Double, nes : Double, gene_count : Int, genes : Array[String]) -> EnrichTerm {
  EnrichTerm::{ term_id, term_name, pvalue, adjusted_pvalue, enrichment_score, nes, gene_count, genes }
}

pub fn EnrichResult::new(terms : Array[EnrichTerm], gene_set_name : String, organism : String) -> EnrichResult {
  EnrichResult::{ terms, gene_set_name, organism }
}

fn copy_enrich_term_array(arr : Array[EnrichTerm]) -> Array[EnrichTerm] {
  let len = arr.length()
  let result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    result.push(arr[i])
    i = i + 1
  }
  result
}

fn set_enrich_term_array(arr : Array[EnrichTerm], index : Int, value : EnrichTerm) -> Array[EnrichTerm] {
  let len = arr.length()
  let result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    if i == index {
      result.push(value)
    } else {
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

fn sort_terms_by_padj(terms : Array[EnrichTerm]) -> Array[EnrichTerm] {
  let len = terms.length()
  let mut result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    result.push(terms[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < len {
    let mut k = j + 1
    while k < len {
      if result[k].adjusted_pvalue < result[j].adjusted_pvalue {
        let temp = result[j]
        let mut new_result = copy_enrich_term_array(result)
        new_result = set_enrich_term_array(new_result, j, result[k])
        new_result = set_enrich_term_array(new_result, k, temp)
        result = new_result
      }
      k = k + 1
    }
    j = j + 1
  }
  
  result
}

fn sort_terms_by_pvalue(terms : Array[EnrichTerm]) -> Array[EnrichTerm] {
  let len = terms.length()
  let mut result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    result.push(terms[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < len {
    let mut k = j + 1
    while k < len {
      if result[k].pvalue < result[j].pvalue {
        let temp = result[j]
        let mut new_result = copy_enrich_term_array(result)
        new_result = set_enrich_term_array(new_result, j, result[k])
        new_result = set_enrich_term_array(new_result, k, temp)
        result = new_result
      }
      k = k + 1
    }
    j = j + 1
  }
  
  result
}

fn sort_terms_by_es(terms : Array[EnrichTerm]) -> Array[EnrichTerm] {
  let len = terms.length()
  let mut result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    result.push(terms[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < len {
    let mut k = j + 1
    while k < len {
      if result[k].enrichment_score > result[j].enrichment_score {
        let temp = result[j]
        let mut new_result = copy_enrich_term_array(result)
        new_result = set_enrich_term_array(new_result, j, result[k])
        new_result = set_enrich_term_array(new_result, k, temp)
        result = new_result
      }
      k = k + 1
    }
    j = j + 1
  }
  
  result
}

fn sort_terms_by_nes(terms : Array[EnrichTerm]) -> Array[EnrichTerm] {
  let len = terms.length()
  let mut result : Array[EnrichTerm] = Array::new()
  let mut i = 0
  while i < len {
    result.push(terms[i])
    i = i + 1
  }
  
  let mut j = 0
  while j < len {
    let mut k = j + 1
    while k < len {
      if result[k].nes > result[j].nes {
        let temp = result[j]
        let mut new_result = copy_enrich_term_array(result)
        new_result = set_enrich_term_array(new_result, j, result[k])
        new_result = set_enrich_term_array(new_result, k, temp)
        result = new_result
      }
      k = k + 1
    }
    j = j + 1
  }
  
  result
}

pub fn bio_enrichplot_barplot(data : EnrichResult, top_n : Int, value_type : String, title : String) -> String {
  let terms = copy_enrich_term_array(data.terms)
  let sorted = if value_type == "padj" { sort_terms_by_padj(terms) } 
               else if value_type == "pvalue" { sort_terms_by_pvalue(terms) }
               else if value_type == "es" { sort_terms_by_es(terms) }
               else { sort_terms_by_nes(terms) }
  
  let n = if sorted.length() < top_n { sorted.length() } else { top_n }
  
  let mut result = "\n" + title + "\n"
  result = result + "=".repeat(80) + "\n"
  result = result + "Term\tValue\tBar\n"
  result = result + "-".repeat(80) + "\n"
  
  let mut i = 0
  while i < n {
    let term = sorted[i]
    let value = if value_type == "padj" { -@math.log10(term.adjusted_pvalue) }
                else if value_type == "pvalue" { -@math.log10(term.pvalue) }
                else if value_type == "es" { term.enrichment_score }
                else { term.nes }
    
    let max_val = 10.0
    let bar_len = if value > max_val { 50 } else { (@math.floor(value / max_val * 50.0)).to_int() }
    let bar = "█".repeat(if bar_len > 0 { bar_len } else { 0 })
    
    let term_display = if term.term_name.length() > 30 { substring(term.term_name, 0, 30) + "..." } else { term.term_name }
    result = result + term_display + "\t" + value.to_string() + "\t" + bar + "\n"
    i = i + 1
  }
  
  result = result + "\nLegend: bar length proportional to " + value_type + "\n"
  result
}

pub fn bio_enrichplot_dotplot(data : EnrichResult, top_n : Int, title : String) -> String {
  let terms = copy_enrich_term_array(data.terms)
  let sorted = sort_terms_by_padj(terms)
  
  let n = if sorted.length() < top_n { sorted.length() } else { top_n }
  let mut result = "\n" + title + "\n"
  result = result + "=".repeat(80) + "\n"
  result = result + "Term\tES\tpadj\tGeneCount\tBubble\n"
  result = result + "-".repeat(80) + "\n"
  
  let mut i = 0
  while i < n {
    let term = sorted[i]
    let bubble_size = if term.gene_count < 5 { 1 } else if term.gene_count < 10 { 2 } else if term.gene_count < 20 { 3 } else if term.gene_count < 50 { 4 } else { 5 }
    let bubble = "*".repeat(bubble_size)
    let term_display = if term.term_name.length() > 30 { substring(term.term_name, 0, 30) + "..." } else { term.term_name }
    result = result + term_display + "\t" + term.enrichment_score.to_string() + "\t" + term.adjusted_pvalue.to_string() + "\t" + term.gene_count.to_string() + "\t" + bubble + "\n"
    i = i + 1
  }
  
  result = result + "\nLegend: * = small, *** = medium, ***** = large\n"
  result
}

pub fn bio_enrichplot_heatmap(data : EnrichResult, top_n : Int, gene_count : Int, title : String) -> String {
  let terms = copy_enrich_term_array(data.terms)
  let sorted = sort_terms_by_padj(terms)
  
  let n = if sorted.length() < top_n { sorted.length() } else { top_n }
  
  let gene_count_map : @hashmap.HashMap[String, Int] = @hashmap.HashMap([])
  let gene_counts = gene_count_map
  
  let mut i = 0
  while i < n {
    let term = sorted[i]
    let mut j = 0
    while j < term.genes.length() {
      let gene = term.genes[j]
      let current = if gene_counts.contains(gene) { gene_counts.get(gene).unwrap() } else { 0 }
      gene_counts.set(gene, current + 1)
      j = j + 1
    }
    i = i + 1
  }
  
  let gene_entries = gene_counts.iter().collect()
  let gene_len = gene_entries.length()
  let mut selected_genes : Array[String] = Array::new()
  let mut k = 0
  while k < gene_len {
    let (gene, count) = gene_entries[k]
    if count >= 2 {
      let sg_len = selected_genes.length()
      let new_sg : Array[String] = Array::new()
      let mut m = 0
      while m < sg_len {
        new_sg.push(selected_genes[m])
        m = m + 1
      }
      new_sg.push(gene)
      selected_genes = new_sg
    }
    k = k + 1
  }
  
  let sg_len = selected_genes.length()
  if sg_len == 0 {
    return "\n" + title + "\n" + "=".repeat(80) + "\n" + "No genes shared between at least 2 terms\n"
  }
  
  let mut result = "\n" + title + "\n"
  result = result + "=".repeat(80) + "\n"
  
  result = result + "\t"
  let mut l = 0
  while l < sg_len {
    result = result + selected_genes[l] + "\t"
    l = l + 1
  }
  result = result + "\n"
  
  let mut p = 0
  while p < n {
    let term = sorted[p]
    let term_display = if term.term_name.length() > 20 { substring(term.term_name, 0, 20) + "..." } else { term.term_name }
    
    let mut row = term_display + "\t"
    let mut q = 0
    while q < sg_len {
      let gene = selected_genes[q]
      let found = contains_gene(term.genes, gene)
      let cell = if found { "█" } else { "░" }
      row = row + cell + "\t"
      q = q + 1
    }
    result = result + row + "\n"
    p = p + 1
  }
  
  result = result + "\nLegend: █ = gene in term, ░ = gene not in term\n"
  result
}

fn copy_tuple_array(arr : Array[(String, Int)]) -> Array[(String, Int)] {
  let len = arr.length()
  let result : Array[(String, Int)] = Array::new()
  let mut i = 0
  while i < len {
    result.push(arr[i])
    i = i + 1
  }
  result
}

fn set_tuple_array(arr : Array[(String, Int)], index : Int, value : (String, Int)) -> Array[(String, Int)] {
  let len = arr.length()
  let result : Array[(String, Int)] = Array::new()
  let mut i = 0
  while i < len {
    if i == index {
      result.push(value)
    } else {
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

fn contains_gene(genes : Array[String], gene : String) -> Bool {
  let mut i = 0
  while i < genes.length() {
    if genes[i] == gene {
      return true
    }
    i = i + 1
  }
  false
}

pub fn bio_enrichplot_cnetplot(data : EnrichResult, top_n : Int, title : String) -> String {
  let terms = copy_enrich_term_array(data.terms)
  let sorted = sort_terms_by_padj(terms)
  
  let n = if sorted.length() < top_n { sorted.length() } else { top_n }
  
  let colors = ["red", "blue", "green", "purple", "orange", "cyan", "magenta", "yellow"]
  
  let term_colors : @hashmap.HashMap[String, String] = @hashmap.HashMap([])
  let color_map = term_colors
  
  let mut term_nodes : Array[String] = Array::new()
  let mut gene_nodes : Array[String] = Array::new()
  let mut edges : Array[(String, String)] = Array::new()
  
  let mut i = 0
  while i < n {
    let term = sorted[i]
    let term_node = term.term_name
    let color = colors[i % colors.length()]
    
    let tn_len = term_nodes.length()
    let new_tn : Array[String] = Array::new()
    let mut j = 0
    while j < tn_len {
      new_tn.push(term_nodes[j])
      j = j + 1
    }
    new_tn.push(term_node)
    term_nodes = new_tn
    
    color_map.set(term_node, color)
    
    let mut k = 0
    while k < term.genes.length() && k < 10 {
      let gene = term.genes[k]
      let mut found = false
      let gn_len = gene_nodes.length()
      let mut l = 0
      while l < gn_len {
        if gene_nodes[l] == gene {
          found = true
          break
        }
        l = l + 1
      }
      if !found {
        let new_gn_len = gene_nodes.length()
        let new_gn : Array[String] = Array::new()
        let mut m = 0
        while m < new_gn_len {
          new_gn.push(gene_nodes[m])
          m = m + 1
        }
        new_gn.push(gene)
        gene_nodes = new_gn
      }
      
      let e_len = edges.length()
      let new_e : Array[(String, String)] = Array::new()
      let mut n = 0
      while n < e_len {
        new_e.push(edges[n])
        n = n + 1
      }
      new_e.push((term_node, gene))
      edges = new_e
      
      k = k + 1
    }
    i = i + 1
  }
  
  let mut result = "\n" + title + "\n"
  result = result + "=".repeat(80) + "\n"
  result = result + "Nodes:\n"
  
  let mut p = 0
  while p < term_nodes.length() {
    let term = term_nodes[p]
    let col = if color_map.contains(term) { color_map.get(term).unwrap() } else { "gray" }
    result = result + "  [" + col + "] " + term + "\n"
    p = p + 1
  }
  
  result = result + "\nGenes:\n"
  let mut q = 0
  while q < gene_nodes.length() {
    result = result + "  ○ " + gene_nodes[q] + "\n"
    q = q + 1
  }
  
  result = result + "\nEdges (Term → Gene):\n"
  let mut r = 0
  while r < edges.length() {
    let (term, gene) = edges[r]
    result = result + "  " + term + " → " + gene + "\n"
    r = r + 1
  }
  
  result
  
}

pub fn bio_enrichplot_emapplot(data : EnrichResult, similarity_threshold : Double, title : String) -> String {
  let terms = copy_enrich_term_array(data.terms)
  let sorted = sort_terms_by_padj(terms)
  
  let n = sorted.length()
  
  let mut term_sets : Array[@hashmap.HashMap[String, Bool]] = Array::new()
  let mut term_names : Array[String] = Array::new()
  
  let mut i = 0
  while i < n {
    let term = sorted[i]
    let ts_len = term_sets.length()
    let new_ts : Array[@hashmap.HashMap[String, Bool]] = Array::new()
    let mut j = 0
    while j < ts_len {
      new_ts.push(term_sets[j])
      j = j + 1
    }
    new_ts.push(create_gene_hashmap(term.genes))
    term_sets = new_ts
    
    let tn_len = term_names.length()
    let new_tn : Array[String] = Array::new()
    let mut k = 0
    while k < tn_len {
      new_tn.push(term_names[k])
      k = k + 1
    }
    new_tn.push(term.term_name)
    term_names = new_tn
    
    i = i + 1
  }
  
  let mut result = "\n" + title + "\n"
  result = result + "=".repeat(80) + "\n"
  result = result + "Enrichment Map (terms connected by shared genes):\n\n"
  
  let mut l = 0
  while l < n {
    let mut m = l + 1
    while m < n {
      let overlap = count_intersection(term_sets[l], term_sets[m])
      let union = calculate_union_size(term_sets[l], term_sets[m])
      let jaccard = if union > 0 { overlap.to_double() / union.to_double() } else { 0.0 }
      
      if jaccard > similarity_threshold {
        result = result + term_names[l] + " ↔ " + term_names[m]
        result = result + " (Jaccard: " + jaccard.to_string() + ")\n"
      }
      m = m + 1
    }
    l = l + 1
  }
  
  result = result + "\nLegend: edges represent Jaccard similarity > 0.1 between term gene sets\n"
  result
}

fn create_gene_hashmap(genes : Array[String]) -> @hashmap.HashMap[String, Bool] {
  let map : @hashmap.HashMap[String, Bool] = @hashmap.HashMap([])
  let result = map
  let len = genes.length()
  let mut i = 0
  while i < len {
    result.set(genes[i], true)
    i = i + 1
  }
  result
}

fn count_intersection(set1 : @hashmap.HashMap[String, Bool], set2 : @hashmap.HashMap[String, Bool]) -> Int {
  let entries = set1.iter().collect()
  let len = entries.length()
  let mut count = 0
  let mut i = 0
  while i < len {
    let (gene, _) = entries[i]
    if set2.contains(gene) {
      count = count + 1
    }
    i = i + 1
  }
  count
}

fn calculate_union_size(set1 : @hashmap.HashMap[String, Bool], set2 : @hashmap.HashMap[String, Bool]) -> Int {
  let entries1 = set1.iter().collect()
  let entries2 = set2.iter().collect()
  
  let set1_len = entries1.length()
  let set2_len = entries2.length()
  let intersection = count_intersection(set1, set2)
  
  set1_len + set2_len - intersection
}

pub fn bio_enrichplot_summary(data : EnrichResult) -> String {
  let total = data.terms.length()
  
  let mut sig_count = 0
  let mut avg_es = 0.0
  let mut avg_nes = 0.0
  
  let mut i = 0
  while i < total {
    let term = data.terms[i]
    if term.adjusted_pvalue < 0.05 {
      sig_count = sig_count + 1
    }
    avg_es = avg_es + term.enrichment_score
    avg_nes = avg_nes + term.nes
    i = i + 1
  }
  
  let mean_es = if total > 0 { avg_es / total.to_double() } else { 0.0 }
  let mean_nes = if total > 0 { avg_nes / total.to_double() } else { 0.0 }
  
  let mut result = "\nEnrichment Result Summary\n"
  result = result + "=".repeat(60) + "\n"
  result = result + "Gene Set: " + data.gene_set_name + "\n"
  result = result + "Organism: " + data.organism + "\n"
  result = result + "Total Terms: " + total.to_string() + "\n"
  result = result + "Significant Terms (padj < 0.05): " + sig_count.to_string() + "\n"
  result = result + "Average Enrichment Score: " + mean_es.to_string() + "\n"
  result = result + "Average NES: " + mean_nes.to_string() + "\n"
  
  result
}