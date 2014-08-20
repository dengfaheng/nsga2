
# tests nsga-ii script

include("nsga-ii.jl")
using Base.Test


function create_random_fitness(fitness_length::Int)
  # helper
  return rand(1:10000, fitness_length)
end


function test_non_dominated_compare(test_size::Int, fitness_size::Int)
  # test by property 
  tests = Vector[]

  # create vector of random fitness arrays
  for i =1:test_size
    push!(tests, create_random_fitness(fitness_size))
  end

  function all_compare(x,y, op)
    # helper
    for i in zip(x,y)
      if !(op(i[1],i[2]))
        return false
      end
    end
    return true
  end

  for i in tests
    for j in tests
      v = non_dominated_compare(i,j)
      # dominating
      if v == 1
        @test all_compare(i,j, >=) == true
      # dominated
      elseif v == -1
        @test all_compare(i,j, <=) == true
      # non dominated and non dominating
      elseif v == 0
        @test all_compare(i,j, >) == false
        @test all_compare(i,j, <) == false
      end
    end
  end
  return true
end


function test_fast_delete_speed()

  vectors = Dict{Int, Vector{Vector{Int}}}()

  for _=1:1000
    v = Int[]
    for _ = 1:1000
      push!(v, int(rand()*500000))
    end
    fixed_vector = sort(unique(v))
    vectors[length(fixed_vector)] = push!(get(vectors, length(fixed_vector), Vector{Int}[]), fixed_vector)
  end
  vectors = vectors[1000]


  function test_fd(f::Function)
    for i=1:100
      for j=1:100
        a = f(vectors[i], vectors[j])
      end
    end
  end

  println("fast delete function")
  @time test_fd(fast_delete)
  println("set difference function")
  @time test_fd(setdiff)
end


function test_fast_delete_correctness(test_size::Int, vector_size::Int)

  function slow_delete(array::Vector, to_delete::Vector)
    # helper, used to compare with fast_delete
    return filter(x->!(x in to_delete), values)
  end

  # unit test, exhaustive
  for _ = 1:test_size
    values = sort(rand(minInt:maxInt, vector_size))
    deletion = sort(rand(minInt:maxInt, vector_size))
    @test slow_delete(values, deletion) == fast_delete(values, deletion)
  end
  return true
end


function test_non_dominated_sort(population_size::Int, fitness_length::Int)
  # verify property of non dominated sorting results
  population = Population()
  for i =  1:population_size
    push!(population.individuals, Individual([0], create_random_fitness(fitness_length)))
  end

  sorts = non_dominated_sort(population)

  # property 1: individuals of the same front do not dominate each others
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
        fit1 = population.individuals[j].fitness
        fit2 = population.individuals[k].fitness
        @test non_dominated_compare(fit1, fit2) == 0
      end
    end
  end

  # property 2: individuals of a front are undominated by those from lesser domination_fronts
  if length(sorts) > 1
    for i = 1:length(sorts)-1
      for j in sorts[i]
        for k in sorts[i+1]
          fit1 = population.individuals[j].fitness
          fit2 = population.individuals[k].fitness
          @test non_dominated_compare(fit1, fit2) in (0, 1)
        end
      end
    end
  end

  true
end



function test_evaluate_against_others(population_size::Int,
                                      fitness_length::Int,
                                      compare_method = non_dominated_compare)

  # generate the population
  population = Population()
  for _ =  1: population_size
    push!(population.individuals, Individual([42],create_random_fitness(fitness_length)))
  end

  # evaluate all individuals
  domination_information = (Int, Int, Vector{Int})[]
  for index = 1:population_size
    push!(domination_information, evaluate_against_others(population.individuals, index, compare_method))
  end

  # verify that the domination relations computed make sense
  for (index, domination_count, dominators) in domination_information
    if !(isempty(dominators))
      for dominator_index in dominators
        dominated_fitness = population.individuals[index].fitness
        dominating_fitness = population.individuals[dominator_index].fitness
        @test compare_method(dominating_fitness, dominated_fitness) == 1
      end
    end
  end

  true
end



function test_calculate_crowding_distance()
  #create populationulation
  population = Population()
  push!(population.individuals, Individual([0], [0,5]))
  push!(population.individuals, Individual([0], [0,5]))
  push!(population.individuals, Individual([0], [2,2]))
  push!(population.individuals, Individual([0], [3,1]))
  push!(population.individuals, Individual([0], [3,1]))
  push!(population.individuals, Individual([0], [3,1]))
  push!(population.individuals, Individual([0], [5,0]))

  # sort into domination fronts
  domination_fronts = non_dominated_sort(population)

  # calculate crowding distances
  merge!(calculate_crowding_distance(population, domination_fronts[1], 1), population.crowding_distances)

  # test against manually calculated values
  @test population.crowding_distances[[0,5]] == (1,Inf)
  @test population.crowding_distances[[2,2]] == (1,1.4)
  @test population.crowding_distances[[3,1]] == (1,1.0)
  @test population.crowding_distances[[5,0]] == (1,Inf)

end



function test_all()
  # exhaustive
  test_non_dominated_compare(500,3)
  test_evaluate_against_others(500,5)
  test_fast_delete_correctness(1000,1000)
  test_non_dominated_sort(500, 3)
  test_calculate_crowding_distance()
  println("All unit tests succeeded")

  true
end

# 
# function test_main(n::Int)
#   # uses the 0-1 sum to check it does indeed optimize
#   allele = [0,1]
#   ALLELES = Vector{Int}[]
#   for i=1:50
#     push!(ALLELES, allele)
#   end
# 
#   function f(x)
#     # we maximize the sum of the genes
#     v = 0
#     for i in x
#       v+=i[1]
#     end
#     return v
#   end
# 
#   function g(x)
#     # we minimize the sum of the genes
#     v= 0 
#     for i in x
#       v-=i[1]
#     end
#     return v
#   end
# 
#   evalF(x) = [f(x), g(x)]
# 
#   mutationOperator = uniformMutate
#   crossoverOperator = uniformCrossover
#   x =  main(ALLELES,
#                     evalF,
#                     100,
#                     n,
#                     0.1,
#                     0.05,
#                     crossoverOperator,
#                     mutationOperator)
# 
# end
# 
# 
# 
# 
# 
# test_all()
