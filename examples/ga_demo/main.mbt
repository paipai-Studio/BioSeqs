///|
/// Genetic Algorithm Demo - Sequence optimization using GA operations.
/// Inspired by Biopython's Bio.GA module.

///|
fn main {
  println("=== Genetic Algorithm Demo ===")
  println("Inspired by Biopython's Bio.GA module")
  println("")

  println("1. Creating GA Individuals")
  let ind1 = @bio.GAIndividual::new(sequence="ATCGATCG", fitness=0.5, generation=0, id=1)
  let ind2 = @bio.GAIndividual::new(sequence="GGCCGGCC", fitness=0.8, generation=0, id=2)
  println("   Individual 1: seq=\{ind1.sequence}, fitness=\{ind1.fitness.to_string()}")
  println("   Individual 2: seq=\{ind2.sequence}, fitness=\{ind2.fitness.to_string()}")
  println("")

  println("2. Creating a Population")
  let individuals = [
    @bio.GAIndividual::new(sequence="AAAA", fitness=0.3, generation=0, id=0),
    @bio.GAIndividual::new(sequence="CCCC", fitness=0.6, generation=0, id=1),
    @bio.GAIndividual::new(sequence="GGGG", fitness=0.9, generation=0, id=2),
    @bio.GAIndividual::new(sequence="TTTT", fitness=0.4, generation=0, id=3),
    @bio.GAIndividual::new(sequence="ACGT", fitness=0.7, generation=0, id=4),
  ]
  let pop = @bio.GAPopulation::new(individuals=individuals)
  println("   Population size: \{pop.individuals.length().to_string()}")
  println("   Best fitness: \{pop.best_fitness.to_string()}")
  println("   Avg fitness: \{pop.avg_fitness.to_string()}")
  println("   Worst fitness: \{pop.worst_fitness.to_string()}")
  println("")

  println("3. Fitness Functions")
  let target = "ACGTACGT"
  let fitness_fn = @bio.ga_fitness_match(target)
  let test_seqs = ["ACGTACGT", "AAAAACGT", "ACGTAAAA", "AAAAAAAA"]
  for seq in test_seqs {
    let ind = @bio.GAIndividual::new(sequence=seq)
    let fit = fitness_fn(ind)
    println("   Sequence '\{seq}' fitness (match to '\{target}'): \{fit.to_string()}")
  }
  println("")

  println("4. Crossover Operations")
  let parent1 = @bio.GAIndividual::new(sequence="AAAA", fitness=0.8, generation=0, id=10)
  let parent2 = @bio.GAIndividual::new(sequence="CCCC", fitness=0.6, generation=0, id=11)
  let (child1, child2) = @bio.ga_single_point_crossover(parent1, parent2, point=2)
  println("   Parent1: \{parent1.sequence}")
  println("   Parent2: \{parent2.sequence}")
  println("   Child1 (crossover at 2): \{child1.sequence}")
  println("   Child2 (crossover at 2): \{child2.sequence}")
  println("")

  println("5. Mutation Operations")
  let original = @bio.GAIndividual::new(sequence="AAAA", fitness=0.0, generation=0, id=20)
  let mutated = @bio.ga_point_mutate(original, 1.0, "ACGT")
  println("   Original: \{original.sequence}")
  println("   Mutated:  \{mutated.sequence}")
  println("")

  println("6. Selection Operations")
  let tournament_winner = @bio.ga_tournament_select(individuals, 3)
  println("   Tournament selection winner: seq=\{tournament_winner.sequence}, fitness=\{tournament_winner.fitness.to_string()}")
  let roulette_winner = @bio.ga_roulette_select(individuals)
  println("   Roulette selection winner: seq=\{roulette_winner.sequence}, fitness=\{roulette_winner.fitness.to_string()}")
  println("")

  println("7. GA Configuration")
  let config = @bio.GAConfig::new(
    population_size=50,
    sequence_length=8,
    max_generations=100,
    target="ACGTACGT",
    mutation_rate=0.05,
    crossover_rate=0.8,
    elitism=true,
  )
  println("   Population size: \{config.population_size.to_string()}")
  println("   Max generations: \{config.max_generations.to_string()}")
  println("   Mutation rate: \{config.mutation_rate.to_string()}")
  println("   Crossover rate: \{config.crossover_rate.to_string()}")
  println("   Elitism: \{config.elitism.to_string()}")
  println("")

  println("8. Running GA Evolution")
  let target_seq = "ACGTACGT"
  let fit_fn = @bio.ga_fitness_match(target_seq)
  let evo_config = @bio.GAConfig::new(
    population_size=30,
    sequence_length=8,
    max_generations=50,
    target=target_seq,
    mutation_rate=0.05,
    crossover_rate=0.85,
    elitism=true,
    convergence_threshold=0.01,
  )
  let result = @bio.ga_evolve(evo_config, fit_fn)
  println("   Generations reached: \{result.generation_reached.to_string()}")
  println("   Converged: \{result.converged.to_string()}")
  println("   Best individual: seq=\{result.best_individual.sequence}, fitness=\{result.best_individual.fitness.to_string()}")
  println("   Best fitness history (first 5):")
  let mut i = 0
  while i < 5 && i < result.best_fitness_history.length() {
    println("     Gen \{i.to_string()}: best = \{result.best_fitness_history[i].to_string()}")
    i = i + 1
  }
  println("   Avg fitness history (first 5):")
  let mut j = 0
  while j < 5 && j < result.avg_fitness_history.length() {
    println("     Gen \{j.to_string()}: avg = \{result.avg_fitness_history[j].to_string()}")
    j = j + 1
  }
  println("")

  println("=== Demo completed successfully! ===")
  println("")
  println("Key concepts:")
  println("- GA evolves sequences toward a target fitness objective")
  println("- Selection chooses which individuals reproduce")
  println("- Crossover combines two parents to create offspring")
  println("- Mutation introduces random variation")
  println("- Elitism preserves the best individuals each generation")
}
