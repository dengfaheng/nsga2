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
