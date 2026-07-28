///|
/// Bio.GA - Genetic Algorithms for sequence optimization.
/// Inspired by Biopython's Bio.GA module.
///
/// Provides genetic algorithm operations for optimizing biological sequences
/// (DNA, RNA, protein) toward specific fitness objectives.
///
/// Key concepts:
/// - Population: a collection of individuals (sequences)
/// - Fitness: a measure of how well an individual meets the objective
/// - Selection: choosing individuals for reproduction
/// - Crossover: combining two parents to produce offspring
/// - Mutation: random modification of sequences
/// - Evolution: iteratively improving the population

// ============================================================================
// Data Structures
// ============================================================================

///|
/// A single individual in the GA population.
pub struct GAIndividual {
  sequence : String
  fitness : Double
  generation : Int
  id : Int
}

///|
/// Population of individuals.
pub struct GAPopulation {
  individuals : Array[GAIndividual]
  generation : Int
  best_fitness : Double
  avg_fitness : Double
  worst_fitness : Double
}

///|
/// Configuration parameters for the GA.
pub struct GAConfig {
  population_size : Int
  mutation_rate : Double
  crossover_rate : Double
  elitism : Bool
  tournament_size : Int
  max_generations : Int
  target_fitness : Double
  sequence_length : Int
  alphabet : String
  target : String
  convergence_threshold : Double
  seed : Int
}

///|
/// Result of a GA evolution run.
pub struct GAResult {
  best_individual : GAIndividual
  best_fitness_history : Array[Double]
  avg_fitness_history : Array[Double]
  generation_reached : Int
  converged : Bool
}

///|
/// Statistics for a single generation.
pub struct GAGenerationStats {
  generation : Int
  best_fitness : Double
  avg_fitness : Double
  worst_fitness : Double
  best_sequence : String
}

///|
/// Create a GAGenerationStats.
pub fn GAGenerationStats::new(
  generation~ : Int,
  best_fitness~ : Double,
  avg_fitness~ : Double,
  worst_fitness? : Double = 0.0,
  best_sequence? : String = "",
) -> GAGenerationStats {
  GAGenerationStats::{
    generation,
    best_fitness,
    avg_fitness,
    worst_fitness,
    best_sequence,
  }
}

// ============================================================================
// Constructors
// ============================================================================

///|
/// Create a GA individual.
pub fn GAIndividual::new(
  sequence~ : String,
  fitness? : Double = 0.0,
  generation? : Int = 0,
  id? : Int = 0,
) -> GAIndividual {
  GAIndividual::{ sequence, fitness, generation, id }
}

///|
/// Create a GA population.
pub fn GAPopulation::new(individuals~ : Array[GAIndividual]) -> GAPopulation {
  let mut best = -1.0e308
  let mut worst = 1.0e308
  let mut sum = 0.0
  for ind in individuals {
    if ind.fitness > best {
      best = ind.fitness
    }
    if ind.fitness < worst {
      worst = ind.fitness
    }
    sum = sum + ind.fitness
  }
  let avg = if individuals.length() > 0 { sum / individuals.length().to_double() } else { 0.0 }
  GAPopulation::{
    individuals,
    generation: 0,
    best_fitness: best,
    avg_fitness: avg,
    worst_fitness: worst,
  }
}

///|
/// Create GA configuration.
pub fn GAConfig::new(
  population_size? : Int = 100,
  mutation_rate? : Double = 0.01,
  crossover_rate? : Double = 0.7,
  elitism? : Bool = true,
  tournament_size? : Int = 3,
  max_generations? : Int = 500,
  target_fitness? : Double = 1.0,
  sequence_length? : Int = 100,
  alphabet? : String = "ACGT",
  target? : String = "",
  convergence_threshold? : Double = 0.001,
  seed? : Int = 42,
) -> GAConfig {
  GAConfig::{
    population_size,
    mutation_rate,
    crossover_rate,
    elitism,
    tournament_size,
    max_generations,
    target_fitness,
    sequence_length,
    alphabet,
    target,
    convergence_threshold,
    seed,
  }
}

///|
/// Create a GA result.
pub fn GAResult::new(
  best_individual~ : GAIndividual,
  best_fitness_history? : Array[Double] = [],
  avg_fitness_history? : Array[Double] = [],
  generation_reached? : Int = 0,
  converged? : Bool = false,
) -> GAResult {
  GAResult::{
    best_individual,
    best_fitness_history,
    avg_fitness_history,
    generation_reached,
    converged,
  }
}

// ============================================================================
// Random sequence generation
// ============================================================================

///|
/// Generate a random sequence from the given alphabet.
pub fn ga_random_sequence(length : Int, alphabet : String) -> String {
  let n = alphabet.length()
  if n == 0 || length <= 0 {
    return ""
  }
  let mut result = ""
  let mut i = 0
  while i < length {
    let idx = (ga_pseudo_random() * n.to_double()).to_int()
    let safe_idx = if idx >= n { n - 1 } else { if idx < 0 { 0 } else { idx } }
    result = result + alphabet[safe_idx : safe_idx + 1].to_string()
    i = i + 1
  }
  result
}

///|
/// Simple pseudo-random number generator state.
struct GARNGState {
  state : Int
}

///|
/// Create a new RNG state.
fn ga_rng_new(seed : Int) -> GARNGState {
  GARNGState::{ state: seed }
}

///|
/// Generate next pseudo-random number in [0, 1).
fn ga_rng_next(state : GARNGState) -> (GARNGState, Double) {
  let next_state = (state.state * 1103515245 + 12345) & 0x7FFFFFFF
  let value = next_state.to_double() / 2147483648.0
  (GARNGState::{ state: next_state }, value)
}

///|
/// Simple pseudo-random number generator (deterministic for reproducibility).
pub fn ga_pseudo_random() -> Double {
  let seed = 42
  let state = ga_rng_new(seed)
  let (_, value) = ga_rng_next(state)
  value
}

///|
/// Generate a random integer in [0, max).
pub fn ga_rand_int(max : Int) -> Int {
  let r = ga_pseudo_random()
  (r * max.to_double()).to_int()
}

// ============================================================================
// Fitness functions (return closures for use in ga_evolve)
// ============================================================================

///|
/// Fitness function: count of matches to a target sequence.
/// Returns a fitness function that can be passed to ga_evolve.
pub fn ga_fitness_match(target : String) -> (GAIndividual) -> Double {
  fn(individual : GAIndividual) -> Double {
    let seq = individual.sequence
    let n = if seq.length() < target.length() { seq.length() } else { target.length() }
    if n == 0 {
      return 0.0
    }
    let mut matches = 0
    let mut i = 0
    while i < n {
      if seq[i : i + 1].to_string() == target[i : i + 1].to_string() {
        matches = matches + 1
      }
      i = i + 1
    }
    matches.to_double() / n.to_double()
  }
}

///|
/// Fitness function: GC content closeness to target.
/// Returns a fitness function that can be passed to ga_evolve.
pub fn ga_fitness_gc_content(target_gc : Double) -> (GAIndividual) -> Double {
  fn(individual : GAIndividual) -> Double {
    let seq = individual.sequence
    let n = seq.length()
    if n == 0 {
      return 0.0
    }
    let mut gc_count = 0
    let mut i = 0
    while i < n {
      let base = seq[i : i + 1].to_string()
      if base == "G" || base == "C" {
        gc_count = gc_count + 1
      }
      i = i + 1
    }
    let gc_content = gc_count.to_double() / n.to_double()
    let diff = if gc_content - target_gc < 0.0 { target_gc - gc_content } else { gc_content - target_gc }
    1.0 - diff
  }
}

// ============================================================================
// Population operations
// ============================================================================

///|
/// Evaluate fitness for all individuals in population.
pub fn ga_evaluate_population(
  individuals : Array[GAIndividual],
  fitness_fn : (GAIndividual) -> Double,
) -> Array[GAIndividual] {
  let result : Array[GAIndividual] = Array::new()
  for ind in individuals {
    let fit = fitness_fn(ind)
    result.push(GAIndividual::new(
      sequence=ind.sequence,
      fitness=fit,
      generation=ind.generation,
      id=ind.id,
    ))
  }
  result
}

///|
/// Find the best individual in a population.
pub fn ga_find_best(pop : GAPopulation) -> GAIndividual {
  let mut best = pop.individuals[0]
  for ind in pop.individuals {
    if ind.fitness > best.fitness {
      best = ind
    }
  }
  best
}

// ============================================================================
// Selection operators
// ============================================================================

///|
/// Tournament selection: select best individual from a random tournament.
pub fn ga_tournament_select(
  individuals : Array[GAIndividual],
  tournament_size : Int,
) -> GAIndividual {
  let n = individuals.length()
  if n == 0 {
    return GAIndividual::new(sequence="")
  }
  let mut best_idx = ga_rand_int(n)
  let mut best_fit = individuals[best_idx].fitness
  let mut i = 1
  while i < tournament_size {
    let idx = ga_rand_int(n)
    let fit = individuals[idx].fitness
    if fit > best_fit {
      best_idx = idx
      best_fit = fit
    }
    i = i + 1
  }
  individuals[best_idx]
}

///|
/// Roulette wheel selection (fitness-proportional).
pub fn ga_roulette_select(individuals : Array[GAIndividual]) -> GAIndividual {
  let n = individuals.length()
  if n == 0 {
    return GAIndividual::new(sequence="")
  }

  let min_fit = individuals[0].fitness
  let mut min = min_fit
  for ind in individuals {
    if ind.fitness < min {
      min = ind.fitness
    }
  }

  let shift = if min < 0.0 { -min + 1.0 } else { 0.0 }

  let mut total = 0.0
  for ind in individuals {
    total = total + ind.fitness + shift
  }

  if total <= 0.0 {
    return individuals[ga_rand_int(n)]
  }

  let r = ga_pseudo_random() * total
  let mut cumulative = 0.0
  for ind in individuals {
    cumulative = cumulative + ind.fitness + shift
    if cumulative >= r {
      return ind
    }
  }

  individuals[n - 1]
}

// ============================================================================
// Crossover operators
// ============================================================================

///|
/// Single-point crossover.
pub fn ga_single_point_crossover(
  parent1 : GAIndividual,
  parent2 : GAIndividual,
  point? : Int = 0,
) -> (GAIndividual, GAIndividual) {
  let seq1 = parent1.sequence
  let seq2 = parent2.sequence
  let n = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }

  if n < 2 {
    return (parent1, parent2)
  }

  let actual_point = if point > 0 && point < n { point } else { 1 + ga_rand_int(n - 1) }

  let child1_seq = (seq1[0 : actual_point] + seq2[actual_point : n]).to_string()
  let child2_seq = (seq2[0 : actual_point] + seq1[actual_point : n]).to_string()

  (
    GAIndividual::new(sequence=child1_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
    GAIndividual::new(sequence=child2_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
  )
}

///|
/// Two-point crossover.
pub fn ga_two_point_crossover(
  parent1 : GAIndividual,
  parent2 : GAIndividual,
) -> (GAIndividual, GAIndividual) {
  let seq1 = parent1.sequence
  let seq2 = parent2.sequence
  let n = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }

  if n < 3 {
    return ga_single_point_crossover(parent1, parent2)
  }

  let p1 = ga_rand_int(n - 1) + 1
  let p2 = p1 + 1 + ga_rand_int(n - p1 - 1)

  let child1_seq = (seq1[0 : p1] + seq2[p1 : p2] + seq1[p2 : n]).to_string()
  let child2_seq = (seq2[0 : p1] + seq1[p1 : p2] + seq2[p2 : n]).to_string()

  (
    GAIndividual::new(sequence=child1_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
    GAIndividual::new(sequence=child2_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
  )
}

///|
/// Uniform crossover: each gene independently swapped.
pub fn ga_uniform_crossover(
  parent1 : GAIndividual,
  parent2 : GAIndividual,
  swap_prob? : Double = 0.5,
) -> (GAIndividual, GAIndividual) {
  let seq1 = parent1.sequence
  let seq2 = parent2.sequence
  let n = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }

  let mut child1_seq = ""
  let mut child2_seq = ""

  let mut i = 0
  while i < n {
    if ga_pseudo_random() < swap_prob {
      child1_seq = child1_seq + seq2[i : i + 1].to_string()
      child2_seq = child2_seq + seq1[i : i + 1].to_string()
    } else {
      child1_seq = child1_seq + seq1[i : i + 1].to_string()
      child2_seq = child2_seq + seq2[i : i + 1].to_string()
    }
    i = i + 1
  }

  (
    GAIndividual::new(sequence=child1_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
    GAIndividual::new(sequence=child2_seq, fitness=0.0, generation=parent1.generation + 1, id=0),
  )
}

// ============================================================================
// Mutation operators
// ============================================================================

///|
/// Point mutation: randomly change bases.
pub fn ga_point_mutate(
  individual : GAIndividual,
  mutation_rate : Double,
  alphabet : String,
) -> GAIndividual {
  let seq = individual.sequence
  let n = seq.length()
  if n == 0 {
    return individual
  }

  let mut new_seq = ""
  let mut i = 0
  while i < n {
    if ga_pseudo_random() < mutation_rate {
      let alphabet_idx = ga_rand_int(alphabet.length())
      let new_base = alphabet[alphabet_idx : alphabet_idx + 1].to_string()
      new_seq = new_seq + new_base
    } else {
      new_seq = new_seq + seq[i : i + 1].to_string()
    }
    i = i + 1
  }

  GAIndividual::new(sequence=new_seq, fitness=individual.fitness, generation=individual.generation, id=individual.id)
}

///|
/// Insertion mutation: insert a random base at a random position.
pub fn ga_insertion_mutate(
  individual : GAIndividual,
  alphabet : String,
) -> GAIndividual {
  let seq = individual.sequence
  let n = seq.length()
  if n == 0 {
    return individual
  }

  let pos = ga_rand_int(n + 1)
  let new_base_idx = ga_rand_int(alphabet.length())
  let new_base = alphabet[new_base_idx : new_base_idx + 1].to_string()

  let new_seq = (seq[0 : pos] + new_base + seq[pos : n]).to_string()

  GAIndividual::new(sequence=new_seq, fitness=individual.fitness, generation=individual.generation, id=individual.id)
}

///|
/// Deletion mutation: delete a random base.
pub fn ga_deletion_mutate(individual : GAIndividual) -> GAIndividual {
  let seq = individual.sequence
  let n = seq.length()
  if n <= 1 {
    return individual
  }

  let pos = ga_rand_int(n)
  let new_seq = (seq[0 : pos] + seq[pos + 1 : n]).to_string()

  GAIndividual::new(sequence=new_seq, fitness=individual.fitness, generation=individual.generation, id=individual.id)
}

// ============================================================================
// Evolution loop
// ============================================================================

///|
/// Run a complete GA optimization.
pub fn ga_evolve(
  config : GAConfig,
  fitness_fn : (GAIndividual) -> Double,
) -> GAResult {
  // Initialize population
  let mut individuals : Array[GAIndividual] = Array::new()
  let mut i = 0
  while i < config.population_size {
    let seq = ga_random_sequence(config.sequence_length, config.alphabet)
    let ind = GAIndividual::new(sequence=seq, fitness=0.0, generation=0, id=i)
    individuals.push(ind)
    i = i + 1
  }

  // Evaluate initial fitness
  individuals = ga_evaluate_population(individuals, fitness_fn)

  let mut pop = GAPopulation::new(individuals=individuals)
  let best_fitness_history : Array[Double] = [pop.best_fitness]
  let avg_fitness_history : Array[Double] = [pop.avg_fitness]
  let mut best_individual = ga_find_best(pop)
  let mut generation = 0
  let mut converged = false

  while generation < config.max_generations && !converged {
    generation = generation + 1

    let mut offspring : Array[GAIndividual] = Array::new()

    // Elitism: keep best individual
    if config.elitism {
      offspring.push(best_individual)
    }

    // Fill rest of offspring
    while offspring.length() < config.population_size {
      let parent1 = ga_tournament_select(pop.individuals, config.tournament_size)
      let parent2 = ga_roulette_select(pop.individuals)

      let (child1, child2) = if ga_pseudo_random() < config.crossover_rate {
        ga_single_point_crossover(parent1, parent2)
      } else {
        (parent1, parent2)
      }

      // Add offspring with mutation
      if offspring.length() < config.population_size {
        let mutated1 = ga_point_mutate(child1, config.mutation_rate, config.alphabet)
        offspring.push(mutated1)
      }
      if offspring.length() < config.population_size {
        let mutated2 = ga_point_mutate(child2, config.mutation_rate, config.alphabet)
        offspring.push(mutated2)
      }
    }

    // Evaluate offspring
    offspring = ga_evaluate_population(offspring, fitness_fn)

    // Update population
    pop = GAPopulation::new(individuals=offspring)

    let new_best = ga_find_best(pop)
    if new_best.fitness > best_individual.fitness {
      best_individual = new_best
    }

    best_fitness_history.push(pop.best_fitness)
    avg_fitness_history.push(pop.avg_fitness)

    // Check convergence
    if best_fitness_history.length() >= 2 {
      let prev_best = best_fitness_history[best_fitness_history.length() - 2]
      let curr_best = best_fitness_history[best_fitness_history.length() - 1]
      let improvement = curr_best - prev_best
      if improvement.abs() < config.convergence_threshold || curr_best >= config.target_fitness {
        converged = true
      }
    }
  }

  GAResult::new(
    best_individual=best_individual,
    best_fitness_history=best_fitness_history,
    avg_fitness_history=avg_fitness_history,
    generation_reached=generation,
    converged=converged,
  )
}

///|
/// Sort indices by fitness (descending).
fn ga_sort_indices(individuals : Array[GAIndividual]) -> Array[Int] {
  let n = individuals.length()
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  // Simple bubble sort by fitness descending
  let arr = indices
  let mut j = 0
  while j < n {
    let mut k = j + 1
    while k < n {
      if individuals[arr[k]].fitness > individuals[arr[j]].fitness {
        let temp = arr[j]
        arr[j] = arr[k]
        arr[k] = temp
      }
      k = k + 1
    }
    j = j + 1
  }
  arr
}