module NSGA_II


# Implementation of the NSGA-II multiobjective
# genetic algorithm as described in:

# Revisiting the NSGA-II crowding-distance computation
# Felix-Antoine Fortin
# Marc Parizeau
# Universite Laval, Quebec, PQ, Canada
# GECCO '13 Proceeding of the fifteenth annual conference on Genetic
# and evolutionary computation conference
# Pages 623-630


include("genetic_operators.jl")


#------------------------------------------------------------------------------
#BEGIN type definitions


immutable Individual
  # basic block of the solution
  genes::Vector
  fitness::Vector

  function Individual(genes::Vector, fitness_values::Vector)
    # fitness value is precomputed
    @assert length(genes) != 0
    @assert length(fitness_values) != 0
    new(genes, fitness_values)
  end

  function Individual(genes::Vector, fitness_function::Function)
    # fitness value is to be computed
    @assert length(genes) != 0
    new(genes, fitness_function(genes))
  end
end


type Population
  # the compound of all individuals
  # includes a mapping of fitness values to crowding distance

  individuals::Vector{Individual}
  crowding_distances::Dict{Vector, (Int, FloatingPoint)}

  function Population()
    # initialize empty
    self = new(Individual[], Dict{Vector, (Int, FloatingPoint)}())
  end

  function Population(individuals::Vector{Individual})
    # initialize with individuals but no crowding_distances
    @assert length(individuals) != 0
    d = Dict{Vector, (Int, FloatingPoint)}()
    self = new(individuals, d)
  end

  function Population(individuals::Vector{Individual},
                      crowding_distances::Dict{Vector, (Int, FloatingPoint)})
    # initialize with individuals and crowding_distances
    @assert length(individuals) != 0
    @assert length(distances) != 0
    self = new(individuals, crowding_distances)
  end
end


# hall of fame is a special population to keep
# the best individuals of all generations
typealias HallOfFame Population


#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN helper methods


function non_dominated_compare(a::Vector, b::Vector, comparator = >)
  # non domination comparison operator
  # ([0, 0, 2]  > [0, 0, 1]) =  1
  # ([0, 0, 1] == [0, 1, 0]) =  0
  # ([1, 0, 1]  < [1, 1, 1]) = -1
  @assert length(a) == length(b) "gene vectors must be of same length"
  AdomB = false
  BdomA = false
  for i in zip(a,b)
    if i[1] != i[2]
      if(comparator(i[1], i[2]))
        AdomB = true
      else
        BdomA = true
      end
    end
    if AdomB && BdomA  # immediate return if nondominated
      return 0
    end
  end

  if AdomB
    return 1
  end
  if BdomA
    return -1
  end
  if !AdomB && !BdomA
    return 0
  end
end


function population_init{T}(initializing_function::Function,
                            fitness_function::Function,
                            population_size::Int)
  # used to initialize the first population in the main loop
  # initializing_function must return a vector
  @assert population_size > 0 "population size doesn't make sense"
  population = Population()
  for _ = 1:population_size
    gene_vector = initializing_function()
    push!(population, Individual(gene_vector, fitness_function(gene_vector)))
  end
  population
end


function evaluate_against_others(population::Population,
                                 self_index::Int,
                                 compare_method::Function)
  # compare the fitness of individual individual at index with rest of the population
  domination_count = 0
  dominated_by = Int[]
  self_fitness = population.individuals[index].fitness

  for (i, other) in enumerate(population.individuals)
    if(i != self_index)
      if compare_method(other.fitness, self_fitness) == 1
        domination_count += 1
        push!(dominated_by, i)
      end
    end
  end

  (self_index, domination_count, dominated_by)
end


function fast_delete{T}(array::Vector{T}, to_delete::Vector{T})
  # we take advantage of the knowledge that both vectors are sorted
  # makes it about 40x faster than setdiff
  # the cost of verifying that the arrays
  @assert issorted(array)
  @assert issorted(to_delete)
  result = Int[]
  deletion_index = 1
  for i in array
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


function non_dominated_sort(double_population::Population,
                            comparison_operator = non_dominated_compare)
  # sort population into m nondominating fronts (best to worst)
  # until at least half the original number of individuals is put in a front

  # get number of individuals to keep
  population_size = length(population.individuals)
  cutoff = population_size / 2

  # get domination information
  # (individual_index, domination_count, dominated_by)
  domination_information = (Int, Int, Vector{Int})[]
  for index = 1:population_size
    push!(domination_information, evaluate_against_others(population, index, comparison_operator))
  end

  fronts_to_indices = Vector{Int}[]

  # iteratively find undominated individuals and separate them from the rest
  # until there are at least half of the double population in them
  while length(domination_information) > cutoff
    current_front_indices = Int[]

    # (individual_index, domination_count, dominated_by)
    tmp_domination_information = (Int, Int, Vector{Int})[]

    for (index, domination_count, dominated_by) in domination_information
      if domination_count == 0
        # the individual is dominating, we add its index to front_indices
        push!(current_front_indices, index)
      else
        # the individual is dominated
        push!(tmp_domination_information, (index, domination_count, dominated_by))
      end
    end

    # push the current front to the result
    push!(fronts_to_indices, current_front_indices)

    # remove the indices of the current front from the dominated individuals
    for (index, domination_count, dominated_by) in tmp_domination_information
      #remove indices from the current front
      substracted = fast_delete(dominated_by, current_front_indices)
      #substract the difference of cardinality
      push!(domination_information, (index, length(substracted), substracted))
    end
  end

  fronts_to_indices
end


function calculate_crowding_distance(population::Population,
                                     front_indices::Vector{Int},
                                     front_index::Int)
  # crowding distance measures the proximity of a
  # solution to its immediate neighbors of the same front. it is used
  # to preserve diversity, later in the algorithm.

  # get the fitnesses from the individuals of the front
  fitnesses = map(ind->ind.fitness, population.individuals[front_indices])

  # calculate mapping {fitness => crowding_distance}
  fitness_to_crowding = Dict{Vector, (Int, FloatingPoint)}()
  for fitness in fitnesses
    fitness_to_crowding[fitness] = (front_index, 0.0)
  end

  # get how many fitnesses and objectives we have
  fitness_keys = collect(keys(fitness_to_crowding))
  fitness_length = length(fitness_keys[1])

  # sort in decreasing order the fitness vectors for each objective
  sorted_by_objective = Vector{Vector{Number}}[]
  objective_range = Number[]

  for i = 1:fitness_length
    sorted = sort(fitness_keys, by = x->x[i], rev = true)
    push!(objective_range, sorted[1][i] - sorted[end][i])
    push!(sorted_by_objective, sorted)
  end

  # assign infinite crowding distance to maximum and
  # minimum fitness of each objective
  map(x -> fitness_to_crowding[x[end]] = (front_index, Inf), sorted_by_objective)
  map(x -> fitness_to_crowding[x[1]]   = (front_index, Inf), sorted_by_objective)

  # assign crowding crowding_distances to the other
  # fitness vectors for each objectives
  for i = 1:fitness_length
    # edge case here! if range == 0, 0 / 0 will give NaN, we ignore the objective in such case
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


function last_front_selection(population::Population,
                              indices::Vector{Int},
                              to_select::Int)
  @assert 0 < to_select <= length(indices) "not enough individuals to select"

  # since individuals within the same front do not dominate each other, they are
  # selected based crowding distance (greater diversity is desired)

  # map {fitness => crowding distance}
  fitness_to_crowding = calculate_crowding_distance(population, last_frontIndices, -1)

  # map {fitness => indices}
  fitness_to_index = Dict{Vector, Vector{Int}}()
  for index in indices
    fitness = population.individuals[index].fitness
    fitness_to_index[fitness] = push!(get(fitness_to_index, fitness, Int[]), index)
  end

  # sort fitness by decreasing crowding distance
  fitness_to_crowding = sort(collect(fitness_to_crowding),
                             by = x -> x[2],  # crowding distance is 2nd field
                             rev = true)

  # choose individuals by iterating through unique fitness list
  # in decreasing order of crowding distance
  chosen_indices = Int[]

  position = 1
  while length(chosen_indices) < to_select
    len = length(fitness_to_index[fitness_to_crowding[position]])

    if len > 1  # multiple individuals with same fitness
      sample = rand(1:len)
      index = fitness_to_index[fitness_to_crowding[position]][sample]
      push!(chosen_indices, index)
      # individuals can be picked only once
      deleteat!(fitness_to_index[fitness_to_crowding[position]], sample)
      j += 1

    else # single individual with this fitness
      index = fitness_to_index[fitness_to_crowding[position]][1]
      push!(chosen_indices, index)
      deleteat!(fitness_to_crowding, position)
    end

    # wrap around
    if position > length(fitness_to_crowding)
      position = 1
    end
  end

  #return the indices of the chosen individuals on the last front
  chosen_indices
end


function select_without_replacement{T}(vector::Vector{T}, k::Int)
  # take k elements from L without replacing
  result = T[]
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


function crowded_compare(first ::(Int, FloatingPoint),
                         second::(Int, FloatingPoint))
  # crowded comparison operator
  # (rank, crowding distance)
  # if rank is the same, tie break with crowding distance
  # if same distance choose randomly
  @assert first[2]>=0
  @assert second[2]>=0
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


function unique_fitness_tournament_selection(population::Population)
  # select across entire range of fitnesses to avoid
  # bias by reoccuring fitnesses

  population_size = length(population.individuals)

  # associate fitness to indices of individuals
  # map {fitness => indices}
  fitness_to_index = Dict{Vector, Vector{Int}}()
  for i = 1:population_size
    value = get(fitness_to_index, population.individuals[i].fitness, Int[])
    fitness_to_index[population.individuals[i].fitness] = push!(value, i)
  end
  fitnesses = collect(keys(fitness_to_index))

  # edge case : only one fitness, return the population as it was
  if length(fitness_to_index) == 1
    return population.individuals
  end

  # else we must select parents
  selected_parents = individual[]

  while length(selected_parents) != population_size
    # we either pick all the fitnesses and select a random individual from them
    # or select a subset of them. depends on how many new parents we still need to add
    k = min((2*(population_size - length(selected_parents))), length(fitness_to_index))

    # sample k fitnesses and get their (front, crowing) from population.distances
    candidate_fitnesses = select_without_replacement(fitnesses, k)
    front_and_crowding = map(x->population.distances[x], candidate_fitnesses)

    # choose the fitnesses
    chosen_fitnesses = Vector[]
    i = 1
    while i < k
      # crowded_compare returns an offset (0 if first solution is better, 1 otherwise)
      selected_index = i + crowded_compare(front_and_crowding[i], front_and_crowding[i+1])
      push!(chosen_fitnesses, candidate_fitnesses[selected_index])
      i += 2
    end

    # we now randomly choose an individual from the indices associated with the chosen fitnesses
    for i in chosen_fitnesses
      chosen_index = fitness_to_index[i][rand(1:length(fitness_to_index[i]))]
      push!(selected_parents, population.individuals[chosen_index])
    end

  end
  selected_parents
end


function generate_children(individuals::Vector{Individual},
                           mutate_population::Function,
                           crossover_population::Function,
                           evaluation_population::Function)
  # final step of the generation, apply mutation and crossover to yield new children
  # both crossover and mutation functions must apply over entire populations
  # apply crossover, mutation
  # and then evaluation new individuals (can fill their fitness values with bogus until evaluation)
  new_population = crossover_population(individuals)
  new_population = mutate_population(new_population)
  new_population
end


function add_to_hall_of_fame(population::Population,
                             indices::Vector{Int},
                             hall_of_fame::HallOfFame,
                             max_hall_of_fame_size)
  # add the best individuals to the hall of fame population to save them for
  # further examination. we merge the first front of the actual population
  # with the rest of the hall of fame to then select the first front of it.

  hall_of_fame.individuals = vcat(hall_of_fame.individuals, population.individuals[indices])

  # acquire the domination information
  domination_information = map(x -> evaluate_against_others(hall_of_fame, x, non_dominated_compare), range(1, length(hall_of_fame.individuals)))


  # filter the dominated individuals
  hall_of_fame.individuals = hall_of_fame.individuals[map(x->x[1], filter(x->x[2]==0, domination_information))]


  # elmiminate duplicates genes (since it is elitist, same individuals may reappear)
  selected_indices = Int[]
  genes = Set{Vector}()
  for (index, individual) in enumerate(hall_of_fame.individuals)
    if !(individual.genes in genes)
      push!(genes, individual.genes)
      push!(selected_indices, index)
    end
  end

  hall_of_fame.individuals = hall_of_fame.individuals[selected_indices]

  if length(hall_of_fame.individuals) > max_hall_of_fame_size
    hall_of_fame.individuals = select_without_replacement(hall_of_fame.individuals, max_hall_of_fame_size)
  end
end


#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN main

function main(initialize_population::Function,
              fitness_function::Function,
              crossover_population::Function,
              mutate_population::Function,
              population_size::Int,
              number_of_generations::Int,
              max_hall_of_fame_size::Int)
  @assert population_size > 0
  @assert number_of_generations >= 0
  @assert max_hall_of_fame_size >= 0

  # hall of fame will keep 
  hall_of_fame = HallOfFame()

  # initialize two populations
  initial_population = initializePopulation(alleles, fitness_function, population_size)
  previous_population = initializePopulation(alleles, fitness_function, population_size)

  # merge them
  merged_population = Population(vcat(initial_population.individuals, previous_population.individuals))

  # main loop: |selection -> generation| -> |selection -> generation| -> ...
  for iteration = 1:number_of_generations

    # sort the merged population into non dominated fronts
    domination_fronts = non_dominated_sort(merged_population)

    # add the best individuals to the hall of fame
    add_to_hall_of_fame(merged_population, domination_fronts[1], hall_of_fame)


    if length(domination_fronts) == 1 || length(domination_fronts[1]) >= population_size
      # edge case: only one front to select individuals from
        merge!(population.crowding_distances, calculate_crowding_distance(merged_population, domination_fronts[1], 1))
        selected_indices = last_front_selection(merged_population, domination_fronts[1], population_size)
    else
        # separate last front from rest, it is treated differently with
        last_front = domination_fronts[end]
        domination_fronts = domination_fronts[1: (end - 1)]

        for (index, front) in enmuerate(domination_fronts)
          # update the crowding distances
          merge!(population.crowding_distances, calculate_crowding_distance(merged_population, front, index))
        end

        # calculate how many individuals are left to select (there's n-k in the previous fronts)
        to_select = population_size - length(reduce(vcat, domination_fronts))

        # find the indices of the k individuals we need from the last front
        selected_indices = last_front_selection(merged_population, last_front, to_select)

        # update the crowding distance on the last front
        merge!(merged_population.crowding_distances, calculate_crowding_distance(merged_population, selectedFromLastFront, indexOfLastFront))

        # put the indices of the individuals in all
        # fronts that were selected as parents
        selected_indices = vcat(reduce(vcat, domination_fronts), selected_indices)
    end

    parent_population = Population(merged_population.individuals[selected_indices],
                                  merged_population.distances)


    #we make a tournament selection to select children
    #the templates are actual parents
    individuals = unique_fitness_tournament_selection(parent_population)

    # apply genetic operators (recomination and mutation) to obtain next pop
    next_population = generate_children(individuals, mutate_population, crossover_population, evaluation_population)

    merged_population = Population(vcat(next_population.individuals, previous_population.individuals))
    previous_population = next_population

  end

  (hall_of_fame, previous_population)
end

#END
#------------------------------------------------------------------------------
