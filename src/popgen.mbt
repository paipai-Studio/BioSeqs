///|
/// Population genetics module.

pub fn allele_frequencies(population : Array[String]) -> Map[String, Double] {
  let counts : Map[String, Int] = Map([], capacity=population.length() * 2)
  
  for i = 0; i < population.length(); i = i + 1 {
    let genotype = population[i]
    let allele1 = genotype[0:1].to_owned()
    let allele2 = genotype[1:2].to_owned()
    match counts.get(allele1) {
      Some(c) => counts[allele1] = c + 1
      None => counts[allele1] = 1
    }
    match counts.get(allele2) {
      Some(c) => counts[allele2] = c + 1
      None => counts[allele2] = 1
    }
  }
  
  let total = 2.0 * population.length().to_double()
  let keys = counts.keys().collect()
  let freqs : Map[String, Double] = Map([], capacity=keys.length())
  
  for idx = 0; idx < keys.length(); idx = idx + 1 {
    let allele = keys[idx]
    freqs[allele] = counts[allele].to_double() / total
  }
  
  freqs
}

pub fn genotype_frequencies(population : Array[String]) -> Map[String, Double] {
  let counts : Map[String, Int] = Map([], capacity=population.length())
  
  for i = 0; i < population.length(); i = i + 1 {
    let genotype = population[i]
    match counts.get(genotype) {
      Some(c) => counts[genotype] = c + 1
      None => counts[genotype] = 1
    }
  }
  
  let total = population.length().to_double()
  let keys = counts.keys().collect()
  let freqs : Map[String, Double] = Map([], capacity=keys.length())
  
  for idx = 0; idx < keys.length(); idx = idx + 1 {
    let genotype = keys[idx]
    freqs[genotype] = counts[genotype].to_double() / total
  }
  
  freqs
}

pub fn hardy_weinberg_test(population : Array[String]) -> Double {
  let freqs = allele_frequencies(population)
  let obs_freqs = genotype_frequencies(population)
  let obs_keys = obs_freqs.keys().collect()
  
  let mut chi2 = 0.0
  for idx = 0; idx < obs_keys.length(); idx = idx + 1 {
    let genotype = obs_keys[idx]
    let obs = obs_freqs[genotype] * population.length().to_double()
    let allele1 = genotype[0:1].to_owned()
    let allele2 = genotype[1:2].to_owned()
    
    let p = match freqs.get(allele1) {
      Some(f) => f
      None => 0.0
    }
    let q = match freqs.get(allele2) {
      Some(f) => f
      None => 0.0
    }
    
    let exp = if allele1 == allele2 {
      population.length().to_double() * p * p
    } else {
      population.length().to_double() * 2.0 * p * q
    }
    
    if exp > 0.0 {
      chi2 = chi2 + (obs - exp) * (obs - exp) / exp
    }
  }
  
  chi2
}

pub fn heterozygosity(population : Array[String]) -> Double {
  let freqs = allele_frequencies(population)
  let keys = freqs.keys().collect()
  
  let mut het = 0.0
  for idx = 0; idx < keys.length(); idx = idx + 1 {
    let p = freqs[keys[idx]]
    het = het + p * (1.0 - p)
  }
  
  het
}

pub fn observed_heterozygosity(population : Array[String]) -> Double {
  let mut hets = 0
  for i = 0; i < population.length(); i = i + 1 {
    let genotype = population[i]
    if genotype[0:1] != genotype[1:2] {
      hets = hets + 1
    }
  }
  
  hets.to_double() / population.length().to_double()
}

pub fn inbreeding_coefficient(population : Array[String]) -> Double {
  let expected = heterozygosity(population)
  let observed = observed_heterozygosity(population)
  
  if expected > 0.0 {
    (expected - observed) / expected
  } else {
    0.0
  }
}

pub fn allele_richness(population : Array[String]) -> Int {
  let alleles : Map[String, Bool] = Map([], capacity=population.length() * 2)
  
  for i = 0; i < population.length(); i = i + 1 {
    let genotype = population[i]
    let allele1 = genotype[0:1].to_owned()
    let allele2 = genotype[1:2].to_owned()
    alleles[allele1] = true
    alleles[allele2] = true
  }
  
  alleles.keys().collect().length()
}

pub fn watterson_theta(population : Array[String]) -> Double {
  if population.length() < 2 {
    return 0.0
  }
  
  let n = population.length()
  let mut sum_terms = 0.0
  
  for i = 1; i < n; i = i + 1 {
    sum_terms = sum_terms + 1.0 / i.to_double()
  }
  
  let k = allele_richness(population)
  
  k.to_double() / sum_terms
}

pub fn tajima_d(population : Array[String]) -> Double {
  if population.length() < 2 {
    return 0.0
  }
  
  let theta_w = watterson_theta(population)
  let pi = heterozygosity(population)
  
  let n = population.length().to_double()
  
  let var_tajima = (1.0 / 3.0) + (1.0 / (18.0 * n * (n - 1.0))) * (11.0 * n * n - 13.0 * n + 4.0)
  
  if var_tajima > 0.0 {
    (pi - theta_w) / var_tajima.sqrt()
  } else {
    0.0
  }
}

pub fn fst(subpop1 : Array[String], subpop2 : Array[String]) -> Double {
  let combined : Array[String] = Array::new()
  for i = 0; i < subpop1.length(); i = i + 1 {
    combined.push(subpop1[i])
  }
  for i = 0; i < subpop2.length(); i = i + 1 {
    combined.push(subpop2[i])
  }
  let h_t = heterozygosity(combined)
  let h_s1 = heterozygosity(subpop1)
  let h_s2 = heterozygosity(subpop2)
  let h_s = (h_s1 + h_s2) / 2.0
  
  if h_t > 0.0 {
    (h_t - h_s) / h_t
  } else {
    0.0
  }
}

pub fn pairwise_fst(populations : Array[Array[String]]) -> Array[Array[Double]] {
  let n = populations.length()
  let matrix = Array::make(n, Array::make(n, 0.0))
  
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let f = fst(populations[i], populations[j])
      matrix[i][j] = f
      matrix[j][i] = f
    }
  }
  
  matrix
}

pub fn fst_multiple(populations : Array[Array[String]]) -> Double {
  if populations.length() < 2 {
    return 0.0
  }
  
  let combined : Array[String] = Array::new()
  for i = 0; i < populations.length(); i = i + 1 {
    for j = 0; j < populations[i].length(); j = j + 1 {
      combined.push(populations[i][j])
    }
  }
  
  let h_t = heterozygosity(combined)
  
  let mut h_s_sum = 0.0
  for i = 0; i < populations.length(); i = i + 1 {
    h_s_sum = h_s_sum + heterozygosity(populations[i])
  }
  let h_s = h_s_sum / populations.length().to_double()
  
  if h_t > 0.0 {
    (h_t - h_s) / h_t
  } else {
    0.0
  }
}

pub fn nucleotide_diversity(population : Array[String]) -> Double {
  if population.length() < 2 {
    return 0.0
  }
  
  let n = population.length()
  let l = if n > 0 { population[0].length() } else { 0 }
  
  let mut diffs = 0.0
  let mut comparisons = 0.0
  
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      for pos = 0; pos < l; pos = pos + 1 {
        if population[i][pos:pos+1] != population[j][pos:pos+1] {
          diffs = diffs + 1.0
        }
        comparisons = comparisons + 1.0
      }
    }
  }
  
  if comparisons > 0.0 {
    diffs / comparisons
  } else {
    0.0
  }
}

pub fn migrate_rate(fst : Double, ne1 : Int, ne2 : Int) -> Double {
  if fst >= 1.0 {
    return 0.0
  }
  
  let ne = (ne1.to_double() * ne2.to_double()) / (ne1.to_double() + ne2.to_double())
  let m = (1.0 / fst - 1.0) / (4.0 * ne)
  
  if m > 0.0 { m } else { 0.0 }
}