///|
/// Hidden Markov Model for Gene Prediction

pub struct HMM {
  states : Array[String]
  observations : Array[String]
  start_prob : Map[String, Double]
  trans_prob : Map[String, Map[String, Double]]
  emit_prob : Map[String, Map[String, Double]]
}

pub fn HMM::new(states : Array[String], observations : Array[String], start_prob : Map[String, Double], trans_prob : Map[String, Map[String, Double]], emit_prob : Map[String, Map[String, Double]]) -> HMM {
  HMM::{ states, observations, start_prob, trans_prob, emit_prob }
}

fn hmm_log_nat(x : Double) -> Double {
  if x <= 0.0 {
    -10000000000.0
  } else if x == 1.0 {
    0.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let y2 = y * y
    let mut sum = y
    let mut term = y
    let mut i = 1
    while i <= 20 {
      term = term * y2
      sum = sum + term / (2 * i + 1).to_double()
      i = i + 1
    }
    2.0 * sum
  }
}

fn get_start_prob(map : Map[String, Double], key : String) -> Double {
  match Map::get(map, key) {
    Some(v) => v
    None => 0.0
  }
}

fn get_trans_prob(map : Map[String, Map[String, Double]], from : String, to : String) -> Double {
  match Map::get(map, from) {
    Some(row) => match Map::get(row, to) {
      Some(v) => v
      None => 0.0
    }
    None => 0.0
  }
}

fn get_emit_prob(map : Map[String, Map[String, Double]], state : String, obs : String) -> Double {
  match Map::get(map, state) {
    Some(row) => match Map::get(row, obs) {
      Some(v) => v
      None => 0.0
    }
    None => 0.0
  }
}

fn log_prob(p : Double) -> Double {
  hmm_log_nat(p)
}

pub fn HMM::forward(self : HMM, obs : Array[String]) -> Double {
  let n = obs.length()
  if n == 0 { return 0.0 }
  
  let m = self.states.length()
  
  let mut alpha = Array::make(m, 0.0)
  
  for i = 0; i < m; i = i + 1 {
    let state = self.states[i]
    let start = get_start_prob(self.start_prob, state)
    let emit = get_emit_prob(self.emit_prob, state, obs[0])
    alpha[i] = start * emit
  }
  
  for t = 1; t < n; t = t + 1 {
    let new_alpha = Array::make(m, 0.0)
    
    for i = 0; i < m; i = i + 1 {
      let mut sum = 0.0
      for j = 0; j < m; j = j + 1 {
        sum = sum + alpha[j] * get_trans_prob(self.trans_prob, self.states[j], self.states[i])
      }
      new_alpha[i] = sum * get_emit_prob(self.emit_prob, self.states[i], obs[t])
    }
    
    alpha = new_alpha
  }
  
  let mut total = 0.0
  for a in alpha {
    total = total + a
  }
  total
}

pub fn HMM::backward(self : HMM, obs : Array[String]) -> Map[String, Array[Double]] {
  let n = obs.length()
  let m = self.states.length()
  
  let mut beta = Array::make(m, 1.0)
  let beta_list = Array::make(n, [])
  
  beta_list[n - 1] = beta
  
  for t = n - 2; t >= 0; t = t - 1 {
    let new_beta = Array::make(m, 0.0)
    
    for i = 0; i < m; i = i + 1 {
      let mut sum = 0.0
      for j = 0; j < m; j = j + 1 {
        sum = sum + get_trans_prob(self.trans_prob, self.states[i], self.states[j]) * 
                       get_emit_prob(self.emit_prob, self.states[j], obs[t + 1]) * 
                       beta[j]
      }
      new_beta[i] = sum
    }
    
    beta = new_beta
    beta_list[t] = beta
  }
  
  let entries = Array::make(m, ("", []))
  for i = 0; i < m; i = i + 1 {
    let state_beta = Array::make(n, 0.0)
    for t = 0; t < n; t = t + 1 {
      state_beta[t] = beta_list[t][i]
    }
    entries[i] = (self.states[i], state_beta)
  }
  
  Map::from_array(entries)
}

pub fn HMM::viterbi(self : HMM, obs : Array[String]) -> Array[String] {
  let n = obs.length()
  let m = self.states.length()
  
  let mut delta = Array::make(m, 0.0)
  let mut psi = Array::make(m, 0)
  let delta_list = Array::make(n, [])
  let psi_list = Array::make(n, [])
  
  for i = 0; i < m; i = i + 1 {
    let state = self.states[i]
    let start = get_start_prob(self.start_prob, state)
    let emit = get_emit_prob(self.emit_prob, state, obs[0])
    delta[i] = log_prob(start * emit)
  }
  
  delta_list[0] = delta
  psi_list[0] = psi
  
  for t = 1; t < n; t = t + 1 {
    let new_delta = Array::make(m, 0.0)
    let new_psi = Array::make(m, 0)
    
    for i = 0; i < m; i = i + 1 {
      let mut max_prob = -10000000000.0
      let mut max_idx = 0
      
      for j = 0; j < m; j = j + 1 {
        let trans = get_trans_prob(self.trans_prob, self.states[j], self.states[i])
        let log_prob_val = delta[j] + log_prob(trans)
        
        if log_prob_val > max_prob {
          max_prob = log_prob_val
          max_idx = j
        }
      }
      
      let emit = get_emit_prob(self.emit_prob, self.states[i], obs[t])
      new_delta[i] = max_prob + log_prob(emit)
      new_psi[i] = max_idx
    }
    
    delta = new_delta
    psi = new_psi
    delta_list[t] = delta
    psi_list[t] = psi
  }
  
  let mut max_idx = 0
  let mut max_prob = delta[0]
  for i = 1; i < m; i = i + 1 {
    if delta[i] > max_prob {
      max_prob = delta[i]
      max_idx = i
    }
  }
  
  let path = Array::make(n, "")
  path[n - 1] = self.states[max_idx]
  
  for t = n - 2; t >= 0; t = t - 1 {
    max_idx = psi_list[t + 1][max_idx]
    path[t] = self.states[max_idx]
  }
  
  path
}

pub fn HMM::baum_welch(self : HMM, obs : Array[String], max_iter : Int, tol : Double) -> HMM {
  let n = obs.length()
  let m = self.states.length()
  
  let mut hmm = self
  
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let alpha = []
    let beta = []
    
    for t = 0; t < n; t = t + 1 {
      alpha.push(Array::make(m, 0.0))
      beta.push(Array::make(m, 0.0))
    }
    
    for i = 0; i < m; i = i + 1 {
      let state = hmm.states[i]
      alpha[0][i] = get_start_prob(hmm.start_prob, state) * get_emit_prob(hmm.emit_prob, state, obs[0])
    }
    
    for t = 1; t < n; t = t + 1 {
      for i = 0; i < m; i = i + 1 {
        let mut sum = 0.0
        for j = 0; j < m; j = j + 1 {
          sum = sum + alpha[t - 1][j] * get_trans_prob(hmm.trans_prob, hmm.states[j], hmm.states[i])
        }
        alpha[t][i] = sum * get_emit_prob(hmm.emit_prob, hmm.states[i], obs[t])
      }
    }
    
    for i = 0; i < m; i = i + 1 {
      beta[n - 1][i] = 1.0
    }
    
    for t = n - 2; t >= 0; t = t - 1 {
      for i = 0; i < m; i = i + 1 {
        let mut sum = 0.0
        for j = 0; j < m; j = j + 1 {
          sum = sum + get_trans_prob(hmm.trans_prob, hmm.states[i], hmm.states[j]) * 
                         get_emit_prob(hmm.emit_prob, hmm.states[j], obs[t + 1]) * 
                         beta[t + 1][j]
        }
        beta[t][i] = sum
      }
    }
    
    let gamma = Array::make(n, [])
    for t = 0; t < n; t = t + 1 {
      gamma[t] = Array::make(m, 0.0)
      let mut sum = 0.0
      for i = 0; i < m; i = i + 1 {
        sum = sum + alpha[t][i] * beta[t][i]
      }
      for i = 0; i < m; i = i + 1 {
        gamma[t][i] = if sum > 0.0 { alpha[t][i] * beta[t][i] / sum } else { 0.0 }
      }
    }
    
    let xi = Array::make(n - 1, [])
    for t = 0; t < n - 1; t = t + 1 {
      xi[t] = Array::make(m, [])
      for i = 0; i < m; i = i + 1 {
        xi[t][i] = Array::make(m, 0.0)
      }
      let mut sum = 0.0
      for i = 0; i < m; i = i + 1 {
        for j = 0; j < m; j = j + 1 {
          sum = sum + alpha[t][i] * get_trans_prob(hmm.trans_prob, hmm.states[i], hmm.states[j]) * 
                       get_emit_prob(hmm.emit_prob, hmm.states[j], obs[t + 1]) * 
                       beta[t + 1][j]
        }
      }
      for i = 0; i < m; i = i + 1 {
        for j = 0; j < m; j = j + 1 {
          xi[t][i][j] = if sum > 0.0 { 
            alpha[t][i] * get_trans_prob(hmm.trans_prob, hmm.states[i], hmm.states[j]) * 
            get_emit_prob(hmm.emit_prob, hmm.states[j], obs[t + 1]) * 
            beta[t + 1][j] / sum 
          } else { 0.0 }
        }
      }
    }
    
    let new_start_prob_entries = Array::make(m, ("", 0.0))
    for i = 0; i < m; i = i + 1 {
      new_start_prob_entries[i] = (hmm.states[i], gamma[0][i])
    }
    let new_start_prob = Map::from_array(new_start_prob_entries)
    
    let new_trans_prob_entries = Array::make(m, ("", Map::from_array([])))
    for i = 0; i < m; i = i + 1 {
      let row_entries = Array::make(m, ("", 0.0))
      let mut denom = 0.0
      for t = 0; t < n - 1; t = t + 1 {
        denom = denom + gamma[t][i]
      }
      for j = 0; j < m; j = j + 1 {
        let mut numer = 0.0
        for t = 0; t < n - 1; t = t + 1 {
          numer = numer + xi[t][i][j]
        }
        row_entries[j] = (hmm.states[j], if denom > 0.0 { numer / denom } else { 1.0 / m.to_double() })
      }
      new_trans_prob_entries[i] = (hmm.states[i], Map::from_array(row_entries))
    }
    let new_trans_prob = Map::from_array(new_trans_prob_entries)
    
    let new_emit_prob_entries = Array::make(m, ("", Map::from_array([])))
    for i = 0; i < m; i = i + 1 {
      let row_entries = Array::make(hmm.observations.length(), ("", 0.0))
      let mut denom = 0.0
      for t = 0; t < n; t = t + 1 {
        denom = denom + gamma[t][i]
      }
      for j = 0; j < hmm.observations.length(); j = j + 1 {
        let obs_val = hmm.observations[j]
        let mut numer = 0.0
        for t = 0; t < n; t = t + 1 {
          if obs[t] == obs_val {
            numer = numer + gamma[t][i]
          }
        }
        row_entries[j] = (obs_val, if denom > 0.0 { numer / denom } else { 1.0 / hmm.observations.length().to_double() })
      }
      new_emit_prob_entries[i] = (hmm.states[i], Map::from_array(row_entries))
    }
    let new_emit_prob = Map::from_array(new_emit_prob_entries)
    
    let new_hmm = HMM::new(hmm.states, hmm.observations, new_start_prob, new_trans_prob, new_emit_prob)
    
    let old_prob = hmm.forward(obs)
    let new_prob = new_hmm.forward(obs)
    
    if (new_prob - old_prob).abs() < tol {
      hmm = new_hmm
      break
    }
    
    hmm = new_hmm
  }
  
  hmm
}

pub fn create_gene_prediction_hmm() -> HMM {
  let states = ["Start", "Exon", "Intron", "Intergenic", "Stop"]
  let observations = ["A", "T", "C", "G"]
  
  let start_prob = Map::from_array([
    ("Start", 0.1),
    ("Exon", 0.0),
    ("Intron", 0.0),
    ("Intergenic", 0.9),
    ("Stop", 0.0),
  ])
  
  let trans_prob = Map::from_array([
    ("Start", Map::from_array([
      ("Start", 0.0),
      ("Exon", 0.9),
      ("Intron", 0.0),
      ("Intergenic", 0.1),
      ("Stop", 0.0),
    ])),
    ("Exon", Map::from_array([
      ("Start", 0.0),
      ("Exon", 0.9),
      ("Intron", 0.05),
      ("Intergenic", 0.0),
      ("Stop", 0.05),
    ])),
    ("Intron", Map::from_array([
      ("Start", 0.0),
      ("Exon", 0.1),
      ("Intron", 0.8),
      ("Intergenic", 0.0),
      ("Stop", 0.1),
    ])),
    ("Intergenic", Map::from_array([
      ("Start", 0.05),
      ("Exon", 0.0),
      ("Intron", 0.0),
      ("Intergenic", 0.95),
      ("Stop", 0.0),
    ])),
    ("Stop", Map::from_array([
      ("Start", 0.0),
      ("Exon", 0.0),
      ("Intron", 0.0),
      ("Intergenic", 1.0),
      ("Stop", 0.0),
    ])),
  ])
  
  let emit_prob = Map::from_array([
    ("Start", Map::from_array([
      ("A", 0.3),
      ("T", 0.3),
      ("C", 0.2),
      ("G", 0.2),
    ])),
    ("Exon", Map::from_array([
      ("A", 0.3),
      ("T", 0.3),
      ("C", 0.2),
      ("G", 0.2),
    ])),
    ("Intron", Map::from_array([
      ("A", 0.25),
      ("T", 0.25),
      ("C", 0.25),
      ("G", 0.25),
    ])),
    ("Intergenic", Map::from_array([
      ("A", 0.25),
      ("T", 0.25),
      ("C", 0.25),
      ("G", 0.25),
    ])),
    ("Stop", Map::from_array([
      ("A", 0.25),
      ("T", 0.25),
      ("C", 0.25),
      ("G", 0.25),
    ])),
  ])
  
  HMM::new(states, observations, start_prob, trans_prob, emit_prob)
}

pub fn predict_genes(hmm : HMM, dna : String) -> Array[String] {
  let obs = []
  for i = 0; i < dna.length(); i = i + 1 {
    obs.push(dna[i:i+1].to_owned())
  }
  
  hmm.viterbi(obs)
}

pub fn extract_exons(prediction : Array[String], dna : String) -> Array[String] {
  let exons = []
  let mut current_exon = ""
  
  for i = 0; i < prediction.length(); i = i + 1 {
    if prediction[i] == "Exon" {
      current_exon = current_exon + dna[i:i+1].to_owned()
    } else {
      if current_exon.length() > 0 {
        exons.push(current_exon)
        current_exon = ""
      }
    }
  }
  
  if current_exon.length() > 0 {
    exons.push(current_exon)
  }
  
  exons
}

pub fn train_gene_hmm(hmm : HMM, training_sequences : Array[String], max_iter : Int) -> HMM {
  let all_obs = []
  for seq in training_sequences {
    let obs = []
    for i = 0; i < seq.length(); i = i + 1 {
      obs.push(seq[i:i+1].to_owned())
    }
    all_obs.push(obs)
  }
  
  let mut trained_hmm = hmm
  for obs in all_obs {
    trained_hmm = trained_hmm.baum_welch(obs, max_iter, 0.000001)
  }
  
  trained_hmm
}
