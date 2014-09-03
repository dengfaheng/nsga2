

#BEGIN INCLUDES================================================================


include("nsga2.jl")
using Base.Test


#END===========================================================================




#BEGIN TEST MACRO==============================================================


macro run_tests(test_vector)
  quote
    print_with_color(:green, "==================================SUMMARY======================================\n\n")
    local num_tests::Int = 0
    local passed_tests::Int = 0

    for fun in $test_vector
      # test all functions
      local val = fun()
      local fun_name = isgeneric(fun) ? fun.env.name : (:anonymous)
      # check return value and print status
      if val == true
        print_with_color(:green, string(num_tests, " [X] ", fun_name, "\n"))
        passed_tests += 1
      else
        print_with_color(:red, string(num_tests, " [ ] ", fun_name, "\n"))
      end
      num_tests += 1
    end
    # final summary
    if passed_tests != num_tests
      print_with_color(:red, string("\n", passed_tests, " OF ", num_tests, " PASSED\n"))
    else
      print_with_color(:green, "\nALL TESTS PASSED\n")
    end
    print_with_color(:green, "===============================================================================")
  end
end



#END===========================================================================




#BEGIN HELPER METHODS==========================================================


function slow_delete(array::Vector{Int}, to_delete::Vector{Int})
  # helper, used to compare with fast_delete
  filter(x->!(x in to_delete), array)
end


function apply_delete(delete_function::Function, vectors::Vector{Vector{Int}})
  for v1::Vector{Int} in vectors
    for v2::Vector{Int} in vectors
      local tmp = delete_function(v1, v2)
    end
  end
end


function test_fast_delete_speed(test_size::Int = 50)
  # makes sure the speed of fast delete is better than
  # the one from the naive version using "in" iterator
  vectors::Vector{Vector{Int}} = Vector{Int}[]
  for _ = 1:test_size
    v::Vector{Int} = rand(1:100000, 2000)
    v = unique(v)
    if length(v) >= 1000
      v = v[1:1000]
      push!(vectors, sort(v))
    end
  end

  # time fast delete O(n)
  t0 = time()
  a = apply_delete(fast_delete, vectors)
  t0 = time() - t0

  # time slow delete O(n^2)
  t1 = time()
  b = apply_delete(slow_delete, vectors)
  t1 = time() - t1

  t0 < t1
end


function test_fast_delete_correctness(test_size::Int = 50)
  # tests the correctness of fast_delete operator
  # against the simple O(n^2) "in" implementation
  vectors::Vector{Vector{Int}} = Vector{Int}[]
  for _ = 1:test_size
    v::Vector{Int} = rand(1:100000, 2000)
    v = unique(v)
    if length(v) >= 1000
      v = v[1:1000]
      push!(vectors, sort(v))
    end
  end

  for v1::Vector{Int} in vectors
    for v2::Vector{Int} in vectors
      if slow_delete(v1, v2) != fast_delete(v1, v2)
        return false
      end
    end
  end

  true
end


#END===========================================================================




#BEGIN NSGA2 TESTS=============================================================


function test_evaluate_against_others(population_size::Int=500, fitness_length::Int=5)
  # generate the population
  population = Population{Vector{Int}, Vector{Int}}()
  for _ =  1: population_size
    push!(population.individuals,
          Individual{Vector{Int}, Vector{Int}}([42],rand(1:10000, fitness_length)))
  end

  # evaluate all individuals
  domination_information = (Int, Int, Vector{Int})[]
  for index = 1:population_size
    push!(domination_information, evaluate_against_others(population, index))
  end

  # verify that the domination relations computed make sense
  for (index, domination_count, dominators) in domination_information
    if !(isempty(dominators))
      for dominator_index in dominators
        dominated_fitness = population.individuals[index].fitness
        dominating_fitness = population.individuals[dominator_index].fitness
        if non_dominated_compare(dominating_fitness, dominated_fitness) != 1
          return false
        end
      end
    end
  end

  true
end


function test_non_dominated_compare(test_size::Int=500, fitness_size::Int=5)
  # tests non dominated compare operator
  tests::Vector{Vector{Int}} = Vector{Int}[]
  for _ =1:test_size
    push!(tests, rand(1:10000, fitness_size))
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
        if all_compare(i,j, >=) != true
          return false
        end
      # dominated
      elseif v == -1
        if all_compare(i,j, <=) != true
          return false
        end
      # non dominated and non dominating
      elseif v == 0
        if all_compare(i,j, >) != false || all_compare(i,j, <) != false
          return false
        end
      end
    end
  end
  true
end


function test_non_dominated_sort(population_size::Int=500, fitness_length::Int=5)
  # test by property of the sorted individuals
  # individuals from a better front must be undominated by any from a lesser front
  population = Population{Vector{Int}, Vector{Int}}()
  for _ = 1:population_size
    push!(population.individuals,
          Individual{Vector{Int}, Vector{Int}}([0], rand(1:10000, fitness_length)))
  end

  sorts::Vector{Vector{Int}} = Vector{Int}[]
  non_dominated_sort!(population, sorts)

  # property 1: individuals of the same front do not dominate each others
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
        fit1 = population.individuals[j].fitness
        fit2 = population.individuals[k].fitness
        if non_dominated_compare(fit1, fit2) != 0
          return false
        end
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
          if !(non_dominated_compare(fit1, fit2) in (0, 1))
            return false
          end
        end
      end
    end
  end

  true
end


function test_calculate_crowding_distance()
  # test crowding distance calculation
  # same as in Deb's book
  population = Population{Vector{Int}, Vector{Int}}()
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [0,5]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [0,5]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [2,2]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [3,1]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [3,1]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [3,1]))
  push!(population.individuals, Individual{Vector{Int}, Vector{Int}}([0], [5,0]))

  # sort into domination fronts
  domination_fronts = non_dominated_sort(population)

  # calculate crowding distances
  merge!(population.crowding_distances, 
         calculate_crowding_distance(population, domination_fronts[1], 1))

  # test against manually calculated values
  if !(population.crowding_distances[[0,5]] == (1,Inf) &&
       population.crowding_distances[[2,2]] == (1,1.4) &&
       population.crowding_distances[[3,1]] == (1,1.0) &&
       population.crowding_distances[[5,0]] == (1,Inf))
    return false
  end

  true
end


#END===========================================================================




#BEGIN TEST LAUNCHER===========================================================


function run_tests()
  @run_tests [

  # fast delete
  test_fast_delete_speed,
  test_fast_delete_correctness,

  # evaluate against others
  test_evaluate_against_others,

  # non dominated compare
  test_non_dominated_compare,

  # calculate crowding distance
  test_calculate_crowding_distance,

  # non dominated sort
  test_non_dominated_sort
  ]


end

