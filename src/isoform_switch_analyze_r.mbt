///|

pub struct IsoformExpression {
  isoform_id : String
  gene_id : String
  counts : Array[Double]
  tpm : Array[Double]
  fpkm : Array[Double]
}

pub struct IsoformSwitch {
  gene_id : String
  gene_name : String
  isoform_id_1 : String
  isoform_id_2 : String
  dpsi : Double
  dif : Double
  pvalue : Double
  adjusted_pvalue : Double
  switch_direction : String
  consequence_type : Array[String]
}

pub struct SwitchAnalyzeRlist {
  isoform_expression : Array[IsoformExpression]
  isoform_switches : Array[IsoformSwitch]
  samples : Array[String]
  conditions : Array[String]
  gene_annotation : @hashmap.HashMap[String, String]
}

pub fn IsoformExpression::new(isoform_id : String, gene_id : String, counts : Array[Double], 
                              tpm : Array[Double], fpkm : Array[Double]) -> IsoformExpression {
  IsoformExpression::{ isoform_id, gene_id, counts, tpm, fpkm }
}

pub fn IsoformSwitch::new(gene_id : String, gene_name : String, isoform_id_1 : String, 
                          isoform_id_2 : String, dpsi : Double, dif : Double, 
                          pvalue : Double, adjusted_pvalue : Double, 
                          switch_direction : String, consequence_type : Array[String]) -> IsoformSwitch {
  IsoformSwitch::{ gene_id, gene_name, isoform_id_1, isoform_id_2, dpsi, dif, pvalue, adjusted_pvalue, switch_direction, consequence_type }
}

pub fn SwitchAnalyzeRlist::new(samples : Array[String], conditions : Array[String]) -> SwitchAnalyzeRlist {
  SwitchAnalyzeRlist::{ 
    isoform_expression: Array::new(), 
    isoform_switches: Array::new(), 
    samples: samples, 
    conditions: conditions,
    gene_annotation: @hashmap.HashMap([])
  }
}

pub fn SwitchAnalyzeRlist::add_isoform(self : SwitchAnalyzeRlist, isoform : IsoformExpression) -> SwitchAnalyzeRlist {
  let expr_len = self.isoform_expression.length()
  let new_expr : Array[IsoformExpression] = Array::new()
  let mut i = 0
  while i < expr_len {
    new_expr.push(self.isoform_expression[i])
    i = i + 1
  }
  new_expr.push(isoform)
  SwitchAnalyzeRlist::{ ..self, isoform_expression: new_expr }
}

pub fn SwitchAnalyzeRlist::add_gene_annotation(self : SwitchAnalyzeRlist, gene_id : String, gene_name : String) -> SwitchAnalyzeRlist {
  let annot = self.gene_annotation
  annot.set(gene_id, gene_name)
  SwitchAnalyzeRlist::{ ..self, gene_annotation: annot }
}

pub fn bio_isoform_calculate_usage(expression : IsoformExpression) -> Array[Double] {
  let len = expression.tpm.length()
  let usage : Array[Double] = Array::new()
  let mut i = 0
  while i < len {
    usage.push(expression.tpm[i])
    i = i + 1
  }
  usage
}

pub fn bio_isoform_calculate_dpsi(isoform1 : IsoformExpression, isoform2 : IsoformExpression, 
                                   condition1_indices : Array[Int], condition2_indices : Array[Int]) -> Double {
  let usage1 = bio_isoform_calculate_usage(isoform1)
  let usage2 = bio_isoform_calculate_usage(isoform2)
  
  let mut sum_usage1_cond1 = 0.0
  let mut sum_usage1_cond2 = 0.0
  let mut sum_usage2_cond1 = 0.0
  let mut sum_usage2_cond2 = 0.0
  
  let len1 = condition1_indices.length()
  let mut i = 0
  while i < len1 {
    let idx = condition1_indices[i]
    if idx >= 0 && idx < usage1.length() {
      sum_usage1_cond1 = sum_usage1_cond1 + usage1[idx]
      sum_usage2_cond1 = sum_usage2_cond1 + usage2[idx]
    }
    i = i + 1
  }
  
  let len2 = condition2_indices.length()
  let mut j = 0
  while j < len2 {
    let idx = condition2_indices[j]
    if idx >= 0 && idx < usage1.length() {
      sum_usage1_cond2 = sum_usage1_cond2 + usage1[idx]
      sum_usage2_cond2 = sum_usage2_cond2 + usage2[idx]
    }
    j = j + 1
  }
  
  let total_cond1 = sum_usage1_cond1 + sum_usage2_cond1
  let total_cond2 = sum_usage1_cond2 + sum_usage2_cond2
  
  let psi1_cond1 = if total_cond1 > 0.0001 { sum_usage1_cond1 / total_cond1 } else { 0.5 }
  let psi1_cond2 = if total_cond2 > 0.0001 { sum_usage1_cond2 / total_cond2 } else { 0.5 }
  
  psi1_cond2 - psi1_cond1
}

pub fn bio_isoform_calculate_dif(isoform1 : IsoformExpression, isoform2 : IsoformExpression, 
                                  condition1_indices : Array[Int], condition2_indices : Array[Int]) -> Double {
  let usage1 = bio_isoform_calculate_usage(isoform1)
  let usage2 = bio_isoform_calculate_usage(isoform2)
  
  let mut sum_psi_cond1 = 0.0
  let mut count_psi_cond1 = 0
  let mut sum_psi_cond2 = 0.0
  let mut count_psi_cond2 = 0
  
  let len1 = condition1_indices.length()
  let mut i = 0
  while i < len1 {
    let idx = condition1_indices[i]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      sum_psi_cond1 = sum_psi_cond1 + psi
      count_psi_cond1 = count_psi_cond1 + 1
    }
    i = i + 1
  }
  
  let len2 = condition2_indices.length()
  let mut j = 0
  while j < len2 {
    let idx = condition2_indices[j]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      sum_psi_cond2 = sum_psi_cond2 + psi
      count_psi_cond2 = count_psi_cond2 + 1
    }
    j = j + 1
  }
  
  let mean_cond1 = if count_psi_cond1 > 0 { sum_psi_cond1 / count_psi_cond1.to_double() } else { 0.5 }
  let mean_cond2 = if count_psi_cond2 > 0 { sum_psi_cond2 / count_psi_cond2.to_double() } else { 0.5 }
  
  mean_cond2 - mean_cond1
}

pub fn bio_isoform_find_switches(data : SwitchAnalyzeRlist, dpsi_threshold : Double, 
                                  pvalue_threshold : Double) -> Array[IsoformSwitch] {
  let gene_to_isoforms : @hashmap.HashMap[String, Array[IsoformExpression]] = @hashmap.HashMap([])
  let gene_map = gene_to_isoforms
  
  let expr_len = data.isoform_expression.length()
  let mut i = 0
  while i < expr_len {
    let isoform = data.isoform_expression[i]
    let existing = if gene_map.contains(isoform.gene_id) { gene_map.get(isoform.gene_id).unwrap() } else { Array::new() }
    let ex_len = existing.length()
    let updated : Array[IsoformExpression] = Array::new()
    let mut j = 0
    while j < ex_len {
      updated.push(existing[j])
      j = j + 1
    }
    updated.push(isoform)
    gene_map.set(isoform.gene_id, updated)
    i = i + 1
  }
  
  let cond_len = data.conditions.length()
  let condition1_indices : Array[Int] = Array::new()
  let condition2_indices : Array[Int] = Array::new()
  let mut k = 0
  while k < cond_len {
    if data.conditions[k] == "control" {
      condition1_indices.push(k)
    } else if data.conditions[k] == "treatment" {
      condition2_indices.push(k)
    }
    k = k + 1
  }
  
  let mut all_switches : Array[IsoformSwitch] = Array::new()
  
  let entries = gene_map.iter().collect()
  let entries_len = entries.length()
  let mut l = 0
  while l < entries_len {
    let (gene_id, isoforms) = entries[l]
    
    if isoforms.length() >= 2 {
      let iso_len = isoforms.length()
      let mut m = 0
      while m < iso_len {
        let mut n = m + 1
        while n < iso_len {
          let iso1 = isoforms[m]
          let iso2 = isoforms[n]
          
          let dpsi = bio_isoform_calculate_dpsi(iso1, iso2, condition1_indices, condition2_indices)
          let dif = bio_isoform_calculate_dif(iso1, iso2, condition1_indices, condition2_indices)
          
          let pvalue = calculate_pvalue(dpsi, iso1, iso2, condition1_indices, condition2_indices)
          let adjusted_pvalue = pvalue * entries_len.to_double()
          
          if dpsi.abs() >= dpsi_threshold && adjusted_pvalue <= pvalue_threshold {
            let direction = if dpsi > 0.0 { "isoform1_up" } else { "isoform2_up" }
            let consequences = predict_consequence(dpsi)
            
            let gene_name = if data.gene_annotation.contains(gene_id) { data.gene_annotation.get(gene_id).unwrap() } else { gene_id }
            
            let switch_len = all_switches.length()
            let new_switches : Array[IsoformSwitch] = Array::new()
            let mut o = 0
            while o < switch_len {
              new_switches.push(all_switches[o])
              o = o + 1
            }
            new_switches.push(IsoformSwitch::new(
              gene_id, gene_name, iso1.isoform_id, iso2.isoform_id, 
              dpsi, dif, pvalue, adjusted_pvalue, direction, consequences
            ))
            all_switches = new_switches
          }
          
          n = n + 1
        }
        m = m + 1
      }
    }
    
    l = l + 1
  }
  
  let sw_len = all_switches.length()
  let mut sorted : Array[IsoformSwitch] = Array::new()
  let mut p = 0
  while p < sw_len {
    sorted.push(all_switches[p])
    p = p + 1
  }
  
  let mut q = 0
  while q < sw_len {
    let mut r = q + 1
    while r < sw_len {
      if sorted[r].adjusted_pvalue < sorted[q].adjusted_pvalue {
        let temp = sorted[q]
        let new_arr : Array[IsoformSwitch] = Array::new()
        let mut s = 0
        while s < sw_len {
          if s == q {
            new_arr.push(sorted[r])
          } else if s == r {
            new_arr.push(temp)
          } else {
            new_arr.push(sorted[s])
          }
          s = s + 1
        }
        sorted = new_arr
      }
      r = r + 1
    }
    q = q + 1
  }
  
  sorted
}

fn calculate_pvalue(dpsi : Double, iso1 : IsoformExpression, iso2 : IsoformExpression, 
                     cond1 : Array[Int], cond2 : Array[Int]) -> Double {
  let usage1 = bio_isoform_calculate_usage(iso1)
  let usage2 = bio_isoform_calculate_usage(iso2)
  
  let mut sum_psi_cond1 = 0.0
  let mut count_psi_cond1 = 0
  let mut sum_psi_cond2 = 0.0
  let mut count_psi_cond2 = 0
  
  let len1 = cond1.length()
  let mut i = 0
  while i < len1 {
    let idx = cond1[i]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      sum_psi_cond1 = sum_psi_cond1 + psi
      count_psi_cond1 = count_psi_cond1 + 1
    }
    i = i + 1
  }
  
  let len2 = cond2.length()
  let mut j = 0
  while j < len2 {
    let idx = cond2[j]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      sum_psi_cond2 = sum_psi_cond2 + psi
      count_psi_cond2 = count_psi_cond2 + 1
    }
    j = j + 1
  }
  
  if count_psi_cond1 == 0 || count_psi_cond2 == 0 {
    return 1.0
  }
  
  let mean1 = sum_psi_cond1 / count_psi_cond1.to_double()
  let mean2 = sum_psi_cond2 / count_psi_cond2.to_double()
  
  let mut var_sum1 = 0.0
  let mut k = 0
  while k < len1 {
    let idx = cond1[k]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      var_sum1 = var_sum1 + (psi - mean1) * (psi - mean1)
    }
    k = k + 1
  }
  let var1 = var_sum1 / count_psi_cond1.to_double()
  
  let mut var_sum2 = 0.0
  let mut l = 0
  while l < len2 {
    let idx = cond2[l]
    if idx >= 0 && idx < usage1.length() {
      let total = usage1[idx] + usage2[idx]
      let psi = if total > 0.0001 { usage1[idx] / total } else { 0.5 }
      var_sum2 = var_sum2 + (psi - mean2) * (psi - mean2)
    }
    l = l + 1
  }
  let var2 = var_sum2 / count_psi_cond2.to_double()
  
  let variance_sum = var1 / count_psi_cond1.to_double() + var2 / count_psi_cond2.to_double()
  let se = if variance_sum > 0.0 { variance_sum.sqrt() } else { 0.0 }
  
  if se < 0.0001 {
    let diff = mean2 - mean1
    return if diff.abs() > 0.1 { 0.001 } else { 1.0 }
  }
  
  let t_stat = (mean2 - mean1) / se
  2.0 * (1.0 - isoform_normal_cdf(t_stat.abs()))
}

fn isoform_normal_cdf(x : Double) -> Double {
  0.5 * (1.0 + isoform_erf(x / 2.0.sqrt()))
}

fn isoform_erf(x : Double) -> Double {
  let a1 = 0.254829592
  let a2 = -0.284496736
  let a3 = 1.421413741
  let a4 = -1.453152027
  let a5 = 1.061405429
  let p = 0.3275911
  
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let abs_x = x.abs()
  
  let t = 1.0 / (1.0 + p * abs_x)
  let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * @math.exp(-abs_x * abs_x)
  
  sign * y
}

fn predict_consequence(dpsi : Double) -> Array[String] {
  let mut consequences : Array[String] = Array::new()
  
  let usage_diff = dpsi.abs()
  
  if usage_diff > 0.2 {
    let len = consequences.length()
    let new_arr : Array[String] = Array::new()
    let mut i = 0
    while i < len {
      new_arr.push(consequences[i])
      i = i + 1
    }
    new_arr.push("significant_isoform_switch")
    consequences = new_arr
  }
  
  if dpsi > 0.3 {
    let len = consequences.length()
    let new_arr : Array[String] = Array::new()
    let mut i = 0
    while i < len {
      new_arr.push(consequences[i])
      i = i + 1
    }
    new_arr.push("isoform1_increased")
    consequences = new_arr
  } else if dpsi < -0.3 {
    let len = consequences.length()
    let new_arr : Array[String] = Array::new()
    let mut i = 0
    while i < len {
      new_arr.push(consequences[i])
      i = i + 1
    }
    new_arr.push("isoform2_increased")
    consequences = new_arr
  }
  
  let len = consequences.length()
  let new_arr : Array[String] = Array::new()
  let mut i = 0
  while i < len {
    new_arr.push(consequences[i])
    i = i + 1
  }
  new_arr.push("quantitative_change")
  new_arr
}

pub fn bio_isoform_switch_summary(switches : Array[IsoformSwitch]) -> String {
  let mut summary = "\nIsoform Switch Summary\n"
  summary = summary + "=".repeat(60) + "\n"
  
  let total = switches.length()
  let mut up_switches = 0
  let mut down_switches = 0
  let mut sum_dpsi = 0.0
  let mut sum_padj = 0.0
  
  let mut i = 0
  while i < total {
    let s = switches[i]
    if s.switch_direction == "isoform1_up" {
      up_switches = up_switches + 1
    } else {
      down_switches = down_switches + 1
    }
    sum_dpsi = sum_dpsi + s.dpsi.abs()
    sum_padj = sum_padj + s.adjusted_pvalue
    i = i + 1
  }
  
  summary = summary + "Total switches: " + total.to_string() + "\n"
  summary = summary + "Isoform1 up: " + up_switches.to_string() + "\n"
  summary = summary + "Isoform2 up: " + down_switches.to_string() + "\n"
  
  let avg_dpsi = if total > 0 { sum_dpsi / total.to_double() } else { 0.0 }
  let avg_padj = if total > 0 { sum_padj / total.to_double() } else { 0.0 }
  
  summary = summary + "Average |dPSI|: " + avg_dpsi.to_string() + "\n"
  summary = summary + "Average adjusted p-value: " + avg_padj.to_string() + "\n"
  
  summary = summary + "\nTop 5 switches by significance:\n"
  summary = summary + "-".repeat(60) + "\n"
  
  let n = if total < 5 { total } else { 5 }
  let mut j = 0
  while j < n {
    let s = switches[j]
    summary = summary + s.gene_name + "\n"
    summary = summary + "  dPSI: " + s.dpsi.to_string() + ", DIF: " + s.dif.to_string() + "\n"
    summary = summary + "  padj: " + s.adjusted_pvalue.to_string() + "\n"
    summary = summary + "  Direction: " + s.switch_direction + "\n"
    let consequences = s.consequence_type
    let c_len = consequences.length()
    let mut k = 0
    let mut c_str = ""
    while k < c_len {
      if k > 0 {
        c_str = c_str + ", "
      }
      c_str = c_str + consequences[k]
      k = k + 1
    }
    summary = summary + "  Consequences: " + c_str + "\n"
    j = j + 1
  }
  
  summary
}

pub fn bio_isoform_plot_psi(data : SwitchAnalyzeRlist, gene_id : String) -> String {
  let gene_to_isoforms : @hashmap.HashMap[String, Array[IsoformExpression]] = @hashmap.HashMap([])
  let gene_map = gene_to_isoforms
  
  let expr_len = data.isoform_expression.length()
  let mut i = 0
  while i < expr_len {
    let isoform = data.isoform_expression[i]
    let existing = if gene_map.contains(isoform.gene_id) { gene_map.get(isoform.gene_id).unwrap() } else { Array::new() }
    let ex_len = existing.length()
    let updated : Array[IsoformExpression] = Array::new()
    let mut j = 0
    while j < ex_len {
      updated.push(existing[j])
      j = j + 1
    }
    updated.push(isoform)
    gene_map.set(isoform.gene_id, updated)
    i = i + 1
  }
  
  if !gene_map.contains(gene_id) {
    return "\nGene " + gene_id + " not found\n"
  }
  
  let isoforms = gene_map.get(gene_id).unwrap()
  let gene_name = if data.gene_annotation.contains(gene_id) { data.gene_annotation.get(gene_id).unwrap() } else { gene_id }
  
  let mut plot = "\nPSI Plot for " + gene_name + " (" + gene_id + ")\n"
  plot = plot + "=".repeat(80) + "\n"
  
  let sample_len = data.samples.length()
  let mut k = 0
  while k < sample_len {
    let sample = data.samples[k]
    let condition = data.conditions[k]
    
    plot = plot + sample + " (" + condition + "): "
    
    let iso_len = isoforms.length()
    let mut l = 0
    while l < iso_len {
      let isoform = isoforms[l]
      let usage = bio_isoform_calculate_usage(isoform)
      let psi = if k < usage.length() { usage[k] } else { 0.0 }
      
      let mut total = 0.0
      let mut m = 0
      while m < iso_len {
        let u = bio_isoform_calculate_usage(isoforms[m])
        total = total + (if k < u.length() { u[k] } else { 0.0 })
        m = m + 1
      }
      
      let norm_psi = if total > 0.0001 { psi / total * 100.0 } else { 0.0 }
      
      plot = plot + isoform.isoform_id + ": " + norm_psi.to_string() + "%  "
      l = l + 1
    }
    
    plot = plot + "\n"
    k = k + 1
  }
  
  plot
}