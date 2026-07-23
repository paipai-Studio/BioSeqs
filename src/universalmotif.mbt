///|
/// universalmotif module, inspired by Bioconductor universalmotif package.

///|
pub struct S4Motif {
  name : String
  alphabet : String
  pwm : Array[Array[Double]]
  consensus : String
}

///|
pub fn S4Motif::new(name : String, alphabet : String, pwm : Array[Array[Double]], consensus : String) -> S4Motif {
  S4Motif::{ name, alphabet, pwm, consensus }
}

///|
pub fn S4Motif::from_pwm(name : String, alphabet : String, pwm : Array[Array[Double]]) -> S4Motif {
  let consensus = calculate_consensus(pwm, alphabet)
  S4Motif::new(name, alphabet, pwm, consensus)
}

///|
fn calculate_consensus(pwm : Array[Array[Double]], alphabet : String) -> String {
  let letters = if alphabet == "DNA" { ["A", "C", "G", "T"] } else { ["A"] }
  
  let mut result = ""
  
  let mut i = 0
  while i < pwm.length() {
    let mut max_prob = 0.0
    let mut max_idx = 0
    
    let mut j = 0
    while j < pwm[i].length() {
      if pwm[i][j] > max_prob {
        max_prob = pwm[i][j]
        max_idx = j
      }
      j = j + 1
    }
    
    if max_idx < letters.length() {
      result = result + letters[max_idx]
    }
    i = i + 1
  }
  
  result
}

///|
pub fn create_example_s4motif() -> S4Motif {
  let pwm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
  ]
  
  S4Motif::from_pwm("TATA-box", "DNA", pwm)
}