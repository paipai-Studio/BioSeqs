///|
/// Tests for Bio.GA module.

///|
test "ga_individual_creation" {
  let ind = @src.GAIndividual::new(
    sequence="ATCGATCG",
    fitness=0.5,
    generation=1,
    id=0,
  )
  assert_eq(ind.sequence, "ATCGATCG")
  assert_true((ind.fitness - 0.5).abs() < 1.0e-10)
  assert_eq(ind.generation, 1)
  assert_eq(ind.id, 0)
}

///|
test "ga_individual_default" {
  let ind = @src.GAIndividual::new(
    sequence="ACGT",
  )
  assert_eq(ind.sequence, "ACGT")
  assert_true((ind.fitness - 0.0).abs() < 1.0e-10)
  assert_eq(ind.generation, 0)
}

///|
test "ga_population_creation" {
  let individuals = [
    @src.GAIndividual::new(sequence="AAAA", fitness=0.8, generation=0, id=0),
    @src.GAIndividual::new(sequence="CCCC", fitness=0.5, generation=0, id=1),
    @src.GAIndividual::new(sequence="GGGG", fitness=0.9, generation=0, id=2),
  ]
  let pop = @src.GAPopulation::new(individuals=individuals)
  assert_eq(pop.individuals.length(), 3)
  assert_true((pop.best_fitness - 0.9).abs() < 1.0e-10)
  assert_true((pop.avg_fitness - 0.7333333333).abs() < 1.0e-10)
}

///|
test "ga_population_empty" {
  let pop = @src.GAPopulation::new(individuals=[])
  assert_eq(pop.individuals.length(), 0)
}

///|
test "ga_config_default" {
  let config = @src.GAConfig::new()
  assert_eq(config.population_size, 100)
  assert_eq(config.max_generations, 500)
  assert_true((config.mutation_rate - 0.01).abs() < 1.0e-10)
  assert_eq(config.alphabet, "ACGT")
  assert_true(config.elitism)
}

///|
test "ga_config_custom" {
  let config = @src.GAConfig::new(
    population_size=50,
    max_generations=100,
    mutation_rate=0.05,
    crossover_rate=0.8,
    elitism=false,
  )
  assert_eq(config.population_size, 50)
  assert_eq(config.max_generations, 100)
  assert_true((config.mutation_rate - 0.05).abs() < 1.0e-10)
  assert_true((config.crossover_rate - 0.8).abs() < 1.0e-10)
  assert_true(!config.elitism)
}

///|
test "ga_fitness_match" {
  let target = "ACGTACGT"
  let fitness_fn = @src.ga_fitness_match(target)
  let ind = @src.GAIndividual::new(sequence="ACGTACGT")
  let fitness = fitness_fn(ind)
  assert_true((fitness - 1.0).abs() < 1.0e-10)
}

///|
test "ga_fitness_match_partial" {
  let target = "ACGTACGT"
  let fitness_fn = @src.ga_fitness_match(target)
  let ind = @src.GAIndividual::new(sequence="ACGTAAAA")
  let fitness = fitness_fn(ind)
  // 5 matches out of 8 positions = 0.625
  assert_true((fitness - 0.625).abs() < 1.0e-10)
}

///|
test "ga_fitness_gc_content" {
  let fitness_fn = @src.ga_fitness_gc_content(0.5)
  let ind = @src.GAIndividual::new(sequence="ACGT")
  let fitness = fitness_fn(ind)
  assert_true((fitness - 1.0).abs() < 1.0e-10)
}

///|
test "ga_fitness_gc_content_zero" {
  let fitness_fn = @src.ga_fitness_gc_content(0.5)
  let ind = @src.GAIndividual::new(sequence="AAAA")
  let fitness = fitness_fn(ind)
  // 0% GC vs 50% target => diff = 0.5, fitness = 1 - 0.5 = 0.5
  assert_true((fitness - 0.5).abs() < 1.0e-10)
}

///|
test "ga_single_point_crossover" {
  let parent1 = @src.GAIndividual::new(sequence="AAAA")
  let parent2 = @src.GAIndividual::new(sequence="CCCC")
  let (child1, child2) = @src.ga_single_point_crossover(parent1, parent2, point=2)
  assert_eq(child1.sequence.length(), 4)
  assert_eq(child2.sequence.length(), 4)
}

///|
test "ga_point_mutate" {
  let individual = @src.GAIndividual::new(sequence="AAAA")
  let mutated = @src.ga_point_mutate(individual, 1.0, "ACGT")
  assert_eq(mutated.sequence.length(), 4)
}

///|
test "ga_tournament_select" {
  let individuals = [
    @src.GAIndividual::new(sequence="AAAA", fitness=0.8, generation=0, id=0),
    @src.GAIndividual::new(sequence="CCCC", fitness=0.5, generation=0, id=1),
    @src.GAIndividual::new(sequence="GGGG", fitness=0.9, generation=0, id=2),
  ]
  let selected = @src.ga_tournament_select(individuals, 2)
  assert_true(selected.fitness >= 0.5)
}

///|
test "ga_roulette_select" {
  let individuals = [
    @src.GAIndividual::new(sequence="AAAA", fitness=0.8, generation=0, id=0),
    @src.GAIndividual::new(sequence="CCCC", fitness=0.5, generation=0, id=1),
    @src.GAIndividual::new(sequence="GGGG", fitness=0.9, generation=0, id=2),
  ]
  let selected = @src.ga_roulette_select(individuals)
  assert_true(selected.fitness >= 0.5)
}

///|
test "ga_evolve_converges" {
  let config = @src.GAConfig::new(
    population_size=10,
    sequence_length=4,
    max_generations=50,
    target="ACGT",
    mutation_rate=0.1,
    crossover_rate=0.8,
    elitism=true,
  )
  let target = "ACGT"
  let fitness_fn = @src.ga_fitness_match(target)
  let result = @src.ga_evolve(config, fitness_fn)
  assert_true(result.best_fitness_history.length() > 0)
  assert_true(result.avg_fitness_history.length() > 0)
  assert_true(result.generation_reached > 0)
}

///|
test "ga_evolution_loop_improves_fitness" {
  let config = @src.GAConfig::new(
    population_size=20,
    sequence_length=6,
    max_generations=100,
    target="ACGTAC",
    mutation_rate=0.05,
    crossover_rate=0.85,
    elitism=true,
    convergence_threshold=0.01,
  )
  let target = "ACGTAC"
  let fitness_fn = @src.ga_fitness_match(target)
  let result = @src.ga_evolve(config, fitness_fn)
  let initial_best = result.best_fitness_history[0]
  let final_best = result.best_fitness_history[result.best_fitness_history.length() - 1]
  assert_true(final_best >= initial_best)
}

///|
test "ga_result_creation" {
  let best = @src.GAIndividual::new(sequence="ACGT", fitness=1.0, generation=10, id=5)
  let stats = @src.GAGenerationStats::new(generation=10, best_fitness=1.0, avg_fitness=0.8)
  assert_eq(stats.generation, 10)
  assert_true((stats.best_fitness - 1.0).abs() < 1.0e-10)
}

///|
test "ga_diversity_decreases_over_time" {
  let config = @src.GAConfig::new(
    population_size=30,
    sequence_length=8,
    max_generations=200,
    target="ACGTACGT",
    mutation_rate=0.02,
    crossover_rate=0.9,
    elitism=true,
  )
  let target = "ACGTACGT"
  let fitness_fn = @src.ga_fitness_match(target)
  let result = @src.ga_evolve(config, fitness_fn)
  assert_true(result.generation_reached > 0)
}

///|
test "ga_custom_fitness_function" {
  let config = @src.GAConfig::new(
    population_size=15,
    sequence_length=5,
    max_generations=100,
    mutation_rate=0.1,
    crossover_rate=0.75,
    elitism=true,
    convergence_threshold=0.001,
  )
  let gc_target = 0.8
  let fitness_fn = @src.ga_fitness_gc_content(gc_target)
  let result = @src.ga_evolve(config, fitness_fn)
  assert_true(result.generation_reached > 0)
}

///|
test "ga_multiple_runs_consistent" {
  let config = @src.GAConfig::new(
    population_size=10,
    sequence_length=4,
    max_generations=50,
    target="TGCA",
    mutation_rate=0.05,
    crossover_rate=0.8,
    elitism=true,
    seed=42,
  )
  let target = "TGCA"
  let fitness_fn = @src.ga_fitness_match(target)
  let result1 = @src.ga_evolve(config, fitness_fn)
  let result2 = @src.ga_evolve(config, fitness_fn)
  assert_eq(result1.generation_reached, result2.generation_reached)
}

///|
test "ga_termination_threshold" {
  let config = @src.GAConfig::new(
    population_size=20,
    sequence_length=6,
    max_generations=200,
    target="ACGTAC",
    mutation_rate=0.01,
    crossover_rate=0.9,
    elitism=true,
    convergence_threshold=0.005,
  )
  let target = "ACGTAC"
  let fitness_fn = @src.ga_fitness_match(target)
  let result = @src.ga_evolve(config, fitness_fn)
  assert_true(result.generation_reached <= 200)
}