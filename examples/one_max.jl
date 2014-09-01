# README
# this is an example use of nsga2
# the idea is to optimize for a simple objective function
# the function is a combination of vector sum and position dependent scoring

# fitness = (sum(vector), early_sum(vector))
# sum simply sums the vector, the ideal being all 1
# early sum favorizes 1 to the left



include("../src/nsga2.jl")




const gene_size = 30



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


function initialize_genes(gene_size)
  # initialize array of 0-1 values of size gene_size
  genes::Vector{Int} = Int[]
  for _=1:gene_size
    push!(genes, rand(0:1))
  end
  genes
end


function mutate_population{A, B}(population::Population{A, B},
                                  mutation_probability::FloatingPoint)
  @assert 0 <= mutation_probability <= 1
  # crossover operator to be applied over the whole population
  new_population::Population{A, B} = Population{A, B}()
  for individual in population
    if rand() < mutation_probability  # mutate individual
      genes::Vector{A} = A[]
      for gene in individual.genes
        if rand() < 0.5
          push!(genes, mod(gene+1, 2))
        end
      end
      push!(new_population.individuals, Individual{A, B}(genes, evaluate_genes))
    else
      push!(new_population.individuals, individual)
    end
  end
  new_population
end


function crossover_population{A, B}(population::Population{A, B},
                                    crossover_probability::FloatingPoint)
  @assert 0 <= crossover_probability <= 1
  # crossover operator to be applied over the whole population
  new_population = Population{A, B}()
  for individual in population
    if rand() < crossover_probability
      other = population.individuals[rand(1:length(population))]
      crossover_point = rand(1:length(individual.genes))
      push!(new_population.individuals, Individual{A, B}(vcat(individual.genes[1:crossover_point], other.genes[crossover_point:end])), evaluate_genes)
    else
      push!(new_population.individuals, individual)
    end
  end
  new_population
end



