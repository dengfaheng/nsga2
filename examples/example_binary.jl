

#BEGIN README==================================================================


# simple examples of nsga2

# [1] (max-sum(genes), early(genes))
# maximizes the number of ones and optimize the early function which
# sums the ones based on their position in the genes (left is better)

# [2] (max-sum(genes), regular(genes))
# maximizes the number of ones and optimize the regular function which
# maximizes the transitions 0 -> 1 and 1 -> 0 that yields a regular structure


#END===========================================================================




#BEGIN INCLUDES================================================================


include("../src/nsga2.jl")


#END===========================================================================




#BEGIN OBJECTIVE FUNCTIONS=====================================================


function evaluate_genes_1(genes::Vector{Int})
  # [1] (max-sum(genes), early(genes))

  function early_sum(v::Vector{Int})
    # calculate early objective
    total::FloatingPoint = 0.0
    len::Int = length(v)
    for (index, value) in enumerate(v)
      total += value / (index / len)
    end
    total
  end
  (sum(genes), early_sum(genes))
end


function evaluate_genes_2(genes::Vector{Int})
  # [2] (max-sum(genes), regular(genes))

  function regular(v::Vector{Int})
    # calculate early objective
    index::Int = 1
    total::Int = 0
    while index < (length(v) - 1)
      if v[index] != v[index+1]
        total += 1
      end
      index+=1
    end
    total
  end
  (sum(genes), regular(genes))
end


#END===========================================================================




#BEGIN INITIALIZATION FUNCTION=================================================


function initialize_genes(gene_size = 50)
  # initialize array of 0-1 values of size gene_size
  genes::Vector{Int} = Int[]
  for _=1:gene_size
    push!(genes, rand(0:1))
  end
  genes
end


#END===========================================================================




#BEGIN GENETIC OPERATORS=======================================================


function mutate_population{B}(population::Population{Vector{Int}, B},
                              evaluate_genes::Function,
                              mutation_probability::FloatingPoint = 0.1,
                              gene_mutation_probability::FloatingPoint = 0.05)
  @assert 0 <= mutation_probability <= 1
  # crossover operator to be applied over the whole population
  new_population::Population{Vector{Int}, B} = Population{Vector{Int}, B}()
  for individual in population.individuals
    if rand() < mutation_probability  # mutate individual
      genes::Vector{Int} = Int[]
      for gene in individual.genes
        if rand() < gene_mutation_probability
          push!(genes, mod(gene+1, 2))
        else
          push!(genes, gene)
        end
      end
      push!(new_population.individuals, Individual{Vector{Int}, B}(genes, evaluate_genes))
    else
      push!(new_population.individuals, individual)
    end
  end
  new_population
end


function crossover_population{A, B}(population::Population{A, B},
                                    evaluate_genes::Function,
                                    crossover_probability::FloatingPoint = 0.05)
  @assert 0 <= crossover_probability <= 1
  # crossover operator to be applied over the whole population
  new_population::Population{A, B} = Population{A, B}()

  for individual::Individual{A, B} in population.individuals
    if rand() < crossover_probability
      other = population.individuals[rand(1:length(population.individuals))]
      crossover_point = rand(0:length(individual.genes))
      if crossover_point == 0
        push!(new_population.individuals, Individual{A, B}(other.genes, evaluate_genes))
      elseif crossover_point == (length(individual.genes))
        push!(new_population.individuals, Individual{A, B}(individual.genes, evaluate_genes))
      else
        push!(new_population.individuals, Individual{A, B}(vcat(individual.genes[1:crossover_point], other.genes[(crossover_point+1):end]), evaluate_genes))
      end
    else
      push!(new_population.individuals, individual)
    end
  end
  new_population
end


#END===========================================================================




#BEGIN PRETTY PRINT============================================================


function print_genes{A, B}(hall_of_fame::HallOfFame{A, B})
  # pretty prints the genes
  information::Vector{(A, B)} = (A, B)[]

  for individual in hall_of_fame.individuals
    push!(information, (individual.genes, individual.fitness))
  end

  sort!(information, by= x->x[2][2], rev=true)

  for (genes, fitness) in information
    println(genes, " ", fitness)
  end
end


#END===========================================================================




#BEGIN MAIN====================================================================


const population_size = 250
const number_of_generations = 150




example_1() = nsga2(
                    Vector{Int},  # type of the genes
                    (Int, FloatingPoint),  # type of the fitness
                    initialize_genes,  # gene initialization function
                    evaluate_genes_1,  # gene evaluation function
                    population_size,  # size of the population
                    crossover_population,  # crossover operator over population
                    mutate_population,  # mutation operator over population
                    number_of_generations,  # number of generations to go through
                    population_size  # max size of hall of fame
                    )

example_2() = nsga2(
                    Vector{Int},  # type of the genes
                    (Int, FloatingPoint),  # type of the fitness
                    initialize_genes,  # gene initialization function
                    evaluate_genes_2,  # gene evaluation function
                    population_size,  # size of the population
                    crossover_population,  # crossover operator over population
                    mutate_population,  # mutation operator over population
                    number_of_generations,  # number of generations to go through
                    population_size  # max size of hall of fame
                    )

# run the examples
# (hall_of_fame_1, _) = example_1()
# (hall_of_fame_2, _) = example_2()


#END===========================================================================
