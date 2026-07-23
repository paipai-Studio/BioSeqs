///|
/// universalmotif module demo.

///|
pub fn main() {
  println("=== universalmotif Demo ===")
  
  println("\n1. Creating a Motif from PWM")
  let pwm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
  ]
  
  let motif = @src.S4Motif::from_pwm("TATA-box", "DNA", pwm)
  println("Motif name: {}", motif.name)
  println("Alphabet: {}", motif.alphabet)
  println("Motif length: {}", motif.pwm.length())
  println("Consensus sequence: {}", motif.consensus)
  
  println("\n2. Using example motif")
  let example_motif = @src.create_example_s4motif()
  println("Example motif: {}", example_motif.name)
  println("Consensus: {}", example_motif.consensus)
  
  println("\n=== Demo complete ===")
}