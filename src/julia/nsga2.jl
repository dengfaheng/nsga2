
#BEGIN README==================================================================


# Implementation of the NSGA-II multiobjective
# genetic algorithm as described in:

# Revisiting the NSGA-II crowding-distance computation
# Felix-Antoine Fortin
# Marc Parizeau
# Universite Laval, Quebec, PQ, Canada
# GECCO '13 Proceeding of the fifteenth annual conference on Genetic
# and evolutionary computation conference
# Pages 623-630


#END===========================================================================




#BEGIN INCLUDES================================================================


include("progress_meter.jl")


#END===========================================================================




#BEGIN TYPES===================================================================


immutable Individual{A, B}
  # individuals upon which the evolution acts
  genes::A
  fitness::B

  function Individual(genes::A, fitness_values::B)
    @assert length(fitness_values) != 0
    new(genes, fitness_values)
  end

  function Individual(genes::A, fitness_function::Function)
    new(genes, fitness_function(genes))
  end
end


immutable Population{A, B}
  # includes a mapping of fitness values to crowding distance
  individuals::Vector{Individual{A, B}}
  crowding_distances::Dict{B, (Int, FloatingPoint)}

  function Population()
    new(Individual{A, B}[], Dict{B, (Int, FloatingPoint)}())
  end

  function Population(individuals::Vector{Individual{A, B}})
    new(individuals, Dict{B, (Int, FloatingPoint)}())
  end

  function Population(individuals::Vector{Individual{A, B}},
                      crowding_distances::Dict{B, (Int, FloatingPoint)})
    new(individuals, crowding_distances)
  end
end


# special population that keep the best individuals of all generations
typealias HallOfFame Population


#END===========================================================================




#BEGIN MISC METHODS============================================================


function select_without_replacement{T}(vector::Vector{T}, k::Int)
  # take k elements from L without replacing
  @assert length(vector) >= k

  result::Vector{T} = T[]
  vector = deepcopy(vector)
  vector_length = length(vector)
  if k == vector_length
    return vector
  end

  for _ = 1:k
    index = rand(1:vector_length)
    push!(result, vector[index])
    deleteat!(vector, index)
    vector_length -= 1
  end

  result
end


function fast_delete(array::Vector{Int}, to_delete::Vector{Int})
  # take advantage of the knowledge that both vectors are sorted
  # O(n)
  @assert issorted(array)
  @assert issorted(to_delete)
  result::Vector{Int} = Int[]
  deletion_index::Int = 1
  for i::Int in array
    # iterate to the next valid index, value >= to i
    while (to_delete[deletion_index] < i) && (deletion_index < length(to_delete))
      deletion_index += 1
    end
    if i != to_delete[deletion_index]
      push!(result, i)
    end
  end
  result
end


function non_dominated_compare{B}(first::B, second::B)
  # non domination comparison operator
  # used to tell wether first vector is "better" than the second
  # ([0, 0, 2]  > [0, 0, 1]) =  1
  # ([0, 0, 1] == [0, 1, 0]) =  0
  # ([1, 0, 1]  < [1, 1, 1]) = -1
  @assert length(first) == length(second) "vectors must be of same length"
  first_dominates::Bool = false
  second_dominates::Bool = false
  for (i, j) in zip(first,second)
    if i != j
      if(>(i, j))
        first_dominates = true
      else
        second_dominates = true
      end
    end
    if first_dominates && second_dominates  # immediate return if nondominated
      return 0
    end
  end

  if first_dominates
    return 1
  end
  if second_dominates
    return -1
  end
  if (! first_dominates) && (! second_dominates)
    return 0
  end
end


function evaluate_against_others{A, B}(population::Population{A, B}, self_index::Int)
  # compare fitness of individual individual at index with rest of population
  # used to figure out domination fronts
  domination_count::Int = 0
  dominated_by::Vector{Int} = Int[]
  self_fitness::B = population.individuals[self_index].fitness

  for (index, other) in enumerate(population.individuals)
    if(index != self_index)
      if non_dominated_compare(other.fitness, self_fitness) == 1
        domination_count += 1
        push!(dominated_by, index)
      end
    end
  end

  (self_index, domination_count, dominated_by)
end


#END===========================================================================




#BEGIN POPULATION METHODS======================================================


function initialize_population!{A, B}(population::Population{A, B},
                                      initialize_genes::Function,
                                      evaluate_genes::Function,
                                      population_size::Int)
  # used to initialize the first population in the main loop
  # initialize_individual must return an individual
  @assert population_size > 0 "population size doesn't make sense"
  for _ = 1:population_size
    genes::A = initialize_genes()
    fitness::B = evaluate_genes(genes)
    push!(population.individuals, Individual{A, B}(genes, fitness))
  end
end


function non_dominated_sort!{A, B}(population::Population{A, B},
                                   domination_fronts::Vector{Vector{Int}})
  # sort population into nondominating fronts (best to worst) until
  # at least half the original number of individuals is put in a front

  # remove the previous information
  empty!(domination_fronts)

  population_size::Int = length(population.individuals)
  cutoff::Int = ceil(population_size / 2)

  # get domination information: (index, count, dominators)
  domination_information::Vector{(Int, Int, Vector{Int})} = (Int, Int, Vector{Int})[]
  tmp_domination_information::Vector{(Int, Int, Vector{Int})} = (Int, Int, Vector{Int})[]
  for index::Int = 1:population_size
    push!(domination_information, evaluate_against_others(population, index))
  end


  # find nondominated individuals and separate them from the rest
  # until there are at least half of the double population in them
  while length(domination_information) > cutoff
    current_front_indices::Vector{Int} = Int[]

    tmp_domination_information = (Int, Int, Vector{Int})[]

    for (index, domination_count, dominated_by) in domination_information
      if domination_count == 0
        # if the individual is dominating, add its index to front indices
        push!(current_front_indices, index)
      else
        # the individual is dominated
        push!(tmp_domination_information, (index, domination_count, dominated_by))
      end
    end

    # push the current front to the result
    push!(domination_fronts, current_front_indices)
    domination_information = (Int, Int, Vector{Int})[]

    # remove the indices of the current front from the dominated individuals
    for (index, domination_count, dominated_by) in tmp_domination_information
      substracted::Vector{Int} = fast_delete(dominated_by, current_front_indices)
      push!(domination_information, (index, length(substracted), substracted))
    end
  end

  nothing
end


function calculate_crowding_distance{A, B}(population::Population{A, B},
                                           front_indices::Vector{Int},
                                           front_index::Int)
  # crowding distance measures the proximity of a solution to its neighbors
  # it is used to preserve diversity, later in the algorithm

  fitnesses::Vector{B} = B[]

  front::Vector{Individual{A, B}} = population.individuals[front_indices]
  for individual in front
    push!(fitnesses, individual.fitness)
  end

  # map fitness => crowding_distance
  fitness_to_crowding::Dict{B, (Int, FloatingPoint)} = Dict{B, (Int, FloatingPoint)}()
  for fitness in fitnesses
    fitness_to_crowding[fitness] = (front_index, 0.0)
  end

  fitness_keys::Vector{B} = collect(keys(fitness_to_crowding))
  fitness_length::Int = length(fitness_keys[1])

  # sort in decreasing order the fitness vectors for each objective
  sorted_by_objective::Vector{Vector{B}} = Vector{B}[]
  objective_range::Vector{Number} = Number[]

  for i = 1:fitness_length
    sorted = sort(fitness_keys, by = x->x[i], rev = true)
    push!(objective_range, sorted[1][i] - sorted[end][i])
    push!(sorted_by_objective, sorted)
  end

  # assign infinite crowding distance to maximum and
  # minimum fitness of each objective
  map(x->fitness_to_crowding[x[end]] = (front_index, Inf), sorted_by_objective)
  map(x->fitness_to_crowding[x[1]]   = (front_index, Inf), sorted_by_objective)

  # assign crowding crowding_distances to the other
  # fitness vectors for each objectives
  for i = 1:fitness_length
    # edge case here! if range == 0, 0 / 0 will give NaN, ignore the objective in such case
    # in DEAP, this is treated using
    # if crowd[-1][0][i] == crowd[0][0][i]:
    #        continue
    if objective_range[i] != 0
      for j = (2:(length(fitness_keys)-1))
        crowding_distance = fitness_to_crowding[sorted_by_objective[i][j]][2]
        crowding_distance += ((sorted_by_objective[i][j-1][i] - sorted_by_objective[i][j+1][i]) / objective_range[i])
        fitness_to_crowding[sorted_by_objective[i][j]] = (front_index, crowding_distance)
      end
    end
  end

  fitness_to_crowding
end


function last_front_selection{A, B}(population::Population{A, B},
                                    indices::Vector{Int},
                                    to_select::Int)
  @assert 0 < to_select <= length(indices) "not enough individuals to select"

  # since individuals within the same front do not dominate each other, they are
  # selected based crowding distance (greater diversity is desired)

  # map fitness => indices
  fitness_to_index::Dict{B, Vector{Int}} = Dict{B, Vector{Int}}()

  for index in indices
    fitness = population.individuals[index].fitness
    fitness_to_index[fitness] = push!(get(fitness_to_index, fitness, Int[]), index)
  end

  # map fitness => crowding distance
  fitness_to_crowding::Dict{B, (Int, FloatingPoint)} = calculate_crowding_distance(population, indices, -1)

  # sort fitness by decreasing crowding distance
  ordered_fitnesses::Vector{B} = map(x->x[1], sort(collect(fitness_to_crowding), by=x->x[2][2], rev=true))

  # choose individuals by iterating from best to worst (crowding distance)
  chosen_indices = Int[]

  position::Int = 1
  while length(chosen_indices) < to_select

    fitness::B = ordered_fitnesses[position]
    len::Int = length(fitness_to_index[fitness])

    if len > 1  # multiple individuals with same fitness, sample one
      sample_index::Int = rand(1:len)
      real_index::Int = fitness_to_index[fitness][sample_index]
      push!(chosen_indices, real_index)

      # individuals can be picked only once, delete individual index after selection
      deleteat!(fitness_to_index[ordered_fitnesses[position]], sample_index)
      position += 1

    else # single individual with this fitness, delete the fitness after selection
      real_index = fitness_to_index[ordered_fitnesses[position]][1]
      push!(chosen_indices, real_index)
      deleteat!(ordered_fitnesses, position)
    end

    # wrap around
    if position > length(ordered_fitnesses)
      position = 1
    end
  end

  # return the indices of the chosen individuals on the last front
  chosen_indices
end


function crowded_compare(first::(Int, FloatingPoint), second::(Int, FloatingPoint))
  # crowded comparison operator
  # (rank, crowding distance)
  # used to choose between solutions in tournament selection
  # if rank is the same, tie break with crowding distance
  # if same distance choose randomly
  @assert first[2] >= 0
  @assert second[2] >= 0
  # rank comarison
  if first[1] < second[1]
    return 0
  elseif first[1] > second[1]
    return 1
  # crowding distance comparison
  elseif first[2] > second[2]
    return 0
  elseif first[2] < second[2]
    return 1
  # A == B, choose either
  else
    return rand(0:1)
  end
end


function unique_fitness_tournament_selection{A, B}(population::Population{A, B})
  # select across entire range of fitnesses to avoid bias by reoccuring fitness
  population_size::Int = length(population.individuals)

  # map fitness => indices
  fitness_to_index::Dict{B, Vector{Int}} = Dict{B, Vector{Int}}()
  for (index, individual) in enumerate(population.individuals)
    fitness_to_index[individual.fitness] = push!(get(fitness_to_index, individual.fitness, Int[]), index)
  end

  # get unique fitnesses
  fitnesses::Vector{B} = collect(keys(fitness_to_index))

  # edge case: only one fitness, return the population as it was
  if length(fitnesses) == 1
    return population.individuals
  end

  # else  must select parents
  selected_individuals::Vector{Individual{A, B}} = Individual{A, B}[]

  while length(selected_individuals) != population_size
    #  either pick all the fitnesses and select a random individual from them
    # or select a subset of them. depends on how many new parents  still need to add
    k = min((2*(population_size - length(selected_individuals))), length(fitness_to_index))

    # sample k fitnesses and get their (front, crowing) from population.crowding_distances
    candidate_fitnesses::Vector{B} = select_without_replacement(fitnesses, k)
    front_and_crowding = map(x->population.crowding_distances[x], candidate_fitnesses)

    # choose the fitnesses
    chosen_fitnesses::Vector{B} = B[]
    i::Int = 1
    while i < k
      # crowded_compare returns an offset (0 if first solution is better, 1 otherwise)
      selected_index = i + crowded_compare(front_and_crowding[i], front_and_crowding[i+1])
      push!(chosen_fitnesses, candidate_fitnesses[selected_index])
      i += 2
    end

    #  now randomly choose an individual from the indices associated with the chosen fitnesses
    for fitness in chosen_fitnesses
      chosen_index = fitness_to_index[fitness][rand(1:length(fitness_to_index[fitness]))]
      push!(selected_individuals, population.individuals[chosen_index])
    end
  end

  selected_individuals
end


function generate_children{A, B}(population::Population{A, B},
                                 mutate_population::Function,
                                 crossover_population::Function,
                                 evaluate_genes::Function)
  # final step of the generation, apply mutation and crossover to yield new children
  # both crossover and mutation functions must apply over entire populations
  # apply crossover, mutation
  # and then evaluation new individuals (can fill their fitness values with bogus until evaluation)

  new_population::Population{A, B} = crossover_population(population, evaluate_genes)
  new_population = mutate_population(new_population, evaluate_genes)

  new_population
end


function add_to_hall_of_fame!{A, B}(population::Population{A, B},
                                    indices::Vector{Int},
                                    hall_of_fame::HallOfFame{A, B},
                                    max_hall_of_fame_size::Int)
  # add the best individuals to the hall of fame and select those who dominate
  # filter out duplicates (same genes)
  # unique fitness or same fitness and unique genes

  append!(hall_of_fame.individuals, population.individuals[indices])

  # select the nondominated individuals
  genes::Set{A} = Set{A}()
  fitnesses::Set{B} = Set{B}()
  selected_individuals::Vector{Individual{A, B}} = Individual{A, B}[]

  for index in 1:length(hall_of_fame.individuals)
    individual::Individual{A, B} = hall_of_fame.individuals[index]
    if (evaluate_against_others(hall_of_fame, index)[2] == 0) &&  # non dominated
       ((!(individual.fitness in fitnesses)) || (!(individual.genes in genes)))  # unique either in fitness or genes
       push!(selected_individuals, individual)
       push!(genes, individual.genes)
       push!(fitnesses, individual.fitness)
    end
  end

  # select without replacement if hall of fame too big
  if length(selected_individuals) > max_hall_of_fame_size
    selected_individuals = select_without_replacement(selected_individuals, max_hall_of_fame_size)
  end

  empty!(hall_of_fame.individuals)
  append!(hall_of_fame.individuals, selected_individuals)
  nothing
end


#END===========================================================================




#BEGIN MAIN====================================================================

function nsga2{A, B}(::Type{A},
                     ::Type{B},
                     initialize_genes::Function,
                     evaluate_genes::Function,
                     population_size::Int,
                     crossover_population::Function,
                     mutate_population::Function,
                     number_of_generations::Int,
                     max_hall_of_fame_size::Int,
                     rng_seed::Int=42)
  @assert population_size > 0
  @assert number_of_generations >= 0
  @assert max_hall_of_fame_size >= 0

  # set the seeds of the pseudo random number generator
  srand(rng_seed)

  # forward declarations
  # initialize the population
  hall_of_fame::HallOfFame{A, B} = HallOfFame{A, B}()
  initial_population::Population{A, B} = Population{A, B}()
  initialize_population!(initial_population, initialize_genes, evaluate_genes, population_size)
  previous_population::Population{A, B} = Population{A, B}()
  initialize_population!(previous_population, initialize_genes, evaluate_genes, population_size)
  parent_population::Population{A, B} = Population{A, B}()
  merged_population::Population{A, B} = Population{A, B}(vcat(initial_population.individuals, previous_population.individuals))
  domination_fronts::Vector{Vector{Int}} = Vector{Int}[]


  # initialize the progress meter
  p = Progress(number_of_generations, 1, "progress  ", 70)

  # main loop: |selection -> generation| -> |selection -> generation| -> ...
  for iteration::Int = 1:number_of_generations

    # sort the merged population into non dominated fronts
    non_dominated_sort!(merged_population, domination_fronts)

    # add the best individuals to the hall of fame
    add_to_hall_of_fame!(merged_population, domination_fronts[1], hall_of_fame, max_hall_of_fame_size)


    # select the parent individuals
    if length(domination_fronts) == 1 || length(domination_fronts[1]) >= population_size
      # edge case: only one front to select individuals from
        merge!(merged_population.crowding_distances, calculate_crowding_distance(merged_population, domination_fronts[1], 1))
        selected_indices = last_front_selection(merged_population, domination_fronts[1], population_size)
    else
        # separate last front from rest, it is treated differently with
        last_front = domination_fronts[end]
        domination_fronts = domination_fronts[1:(end - 1)]

        for (index, front) in enumerate(domination_fronts)
          # update the crowding distances
          merge!(merged_population.crowding_distances, calculate_crowding_distance(merged_population, front, index))
        end

        # calculate how many individuals are left to select (there's n-k in the previous fronts)
        to_select::Int = population_size - length(reduce(vcat, domination_fronts))

        # find the indices of the k individuals  need from the last front
        selected_indices::Vector{Int} = last_front_selection(merged_population, last_front, to_select)

        # update the crowding distance on the last front
        merge!(merged_population.crowding_distances, calculate_crowding_distance(merged_population, selected_indices, length(domination_fronts) + 1))

        # put the indices of the individuals in all
        # fronts that were selected as parents
        selected_indices = vcat(reduce(vcat, domination_fronts), selected_indices)
    end

    empty!(parent_population.individuals)
    empty!(parent_population.crowding_distances)
    append!(parent_population.individuals, merged_population.individuals[selected_indices])
    merge!(parent_population.crowding_distances, merged_population.crowding_distances)


    #  make a tournament selection to select children
    selected_individuals::Population{A, B} = Population{A, B}(unique_fitness_tournament_selection(parent_population))

    next_population::Population{A, B} = generate_children(selected_individuals, mutate_population, crossover_population, evaluate_genes)

    # change names for the next population
    merged_population = Population{A, B}(vcat(next_population.individuals, previous_population.individuals))
    previous_population = next_population

    # update the progress meter
    next!(p)
  end

  (hall_of_fame, previous_population)
end

#END===========================================================================
