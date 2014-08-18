#------------------------------------------------------------------------------
#BEGIN crossover operators
# some simple crossover operators (uniform, one point crossover)
# see wikipedia or Goldberg's book for more information


function uniform_crossover(self_genes::Vector, other_genes::Vector)
  @assert length(self_genes) == length(other_genes) != 0 "gene vectors must be of same length"

  # copy the gene vectors to avoid side-effect problems
  self_genes = deepcopy(self_genes)
  other_genes = deepcopy(other_genes)

  for i = 1:(length(self_genes))
    if rand() < 0.5
      self_genes[i] = other_genes[i]
    end
  end

  self_genes
end



function one_point_crossover(self_genes::Vector, other_genes::Vector)

  @assert length(self_genes) == length(other_genes) != 0 "gene vectors must be of same length"
  self_genes = deepcopy(self_genes)
  other_genes = deepcopy(other_genes)

  # find beginning of crossover
  point = rand(1:length(self_genes))

  # swap the values after the point
  for i = point:(length(self_genes))
    self_genes[i] = other_genes[i]
  end

  self_genes
end


#END
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#BEGIN mutation operators
# mutation operators modify the gene vector according to a certain probability


function uniform_mutate(genes::Vector, mutate::Function, mutation_probability)
  # each gene has a probability of being mutated
  # mutate function has information about the position and identity of the actual allele
  new_genes = deepcopy(genes)

  for i = 1:(length(new_genes))
    if(rand() < mutation_probability)
      new_genes[i] = mutate(new_genes[i], i)
    end
  end

  new_genes
end


#END
#------------------------------------------------------------------------------





 for i = 1:population_size
    # initialize new genes and a new fitness from individuals genes and fitness
    new_genes = deepcopy(individuals[i].genes)
    newFitness = deepcopy(individuals[i].fitness)
    modified = false

    if evolutionaryEvents[i][1] == true
      modified = true

      #recombination (crossover)
      # randomly choose second parent
      secondParentIndex = rand(1:(population_size-1))

      # leave a gap to not select same parent
      if secondParentIndex >= i
        secondParentIndex += 1
      end

      # combine two individuals genes (on which the fitness is based)
      new_genes = crossover_function(new_genes, individuals[secondParentIndex].genes)
    end

    if evolutionaryEvents[i][2] == true
      modified = true
      # mutation
      new_genes = mutation_function(new_genes, alleles)
    end

    # if modified, re-evaluate
    if modified
      newFitness = evaluation_function(new_genes)
    end

    # add newly created individuals to the children population
    push!(new_population.individuals, Individual(new_genes, newFitness))
  end
