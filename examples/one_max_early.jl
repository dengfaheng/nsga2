# README
# this is an example use of nsga2
# the idea is to optimize for a simple objective function
# the function is a combination of vector sum and position dependent scoring

# fitness
# fitness = (sum(vector), early_sum(vector))
# sum simply sums the vector, the ideal being all 1
# early sum favorizes 1 to the left

# 


include("../src/nsga2.jl")


#BEGIN objective function
function evaluate_genes(genes::Vector{Int})
  # objective is to get as much 1 on the left of the vector
  # and as much 1 as possible in general

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
#END


#BEGIN individual initialization function
function initialize_genes(gene_size = 50)
  # initialize array of 0-1 values of size gene_size
  genes::Vector{Int} = Int[]
  for _=1:gene_size
    push!(genes, rand(0:1))
  end
  genes
end
#END


#BEGIN population mutation function
function mutate_population{B}(population::Population{Vector{Int}, B},
                                  mutation_probability::FloatingPoint = 0.1,
                                  gene_mutation_probability::FloatingPoint = 0.05,
                                  evaluate_genes = evaluate_genes)
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
#END


#BEGIN population crossover function
function crossover_population{A, B}(population::Population{A, B},
                                    crossover_probability::FloatingPoint = 0.05,
                                    evaluate_genes = evaluate_genes)
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
#END


#BEGIN main
const population_size = 150
const number_of_generations = 150


(hall_of_fame, last_population) = nsga2(
                                  Vector{Int},  # type of the genes
                                  (Int, FloatingPoint),  # type of the fitness
                                  initialize_genes,  # gene initialization function
                                  evaluate_genes,  # gene evaluation function
                                  population_size,  # size of the population
                                  crossover_population,  # crossover operator over population
                                  mutate_population,  # mutation operator over population
                                  number_of_generations,  # number of generations to go through
                                  population_size  # max size of hall of fame
                                  )

#END
