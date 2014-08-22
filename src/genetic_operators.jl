

#  for i = 1:population_size
#     # initialize new genes and a new fitness from individuals genes and fitness
#     new_genes = deepcopy(individuals[i].genes)
#     newFitness = deepcopy(individuals[i].fitness)
#     modified = false
# 
#     if evolutionaryEvents[i][1] == true
#       modified = true
# 
#       #recombination (crossover)
#       # randomly choose second parent
#       secondParentIndex = rand(1:(population_size-1))
# 
#       # leave a gap to not select same parent
#       if secondParentIndex >= i
#         secondParentIndex += 1
#       end
# 
#       # combine two individuals genes (on which the fitness is based)
#       new_genes = crossover_function(new_genes, individuals[secondParentIndex].genes)
#     end
# 
#     if evolutionaryEvents[i][2] == true
#       modified = true
#       # mutation
#       new_genes = mutation_function(new_genes, alleles)
#     end
# 
#     # if modified, re-evaluate
#     if modified
#       newFitness = evaluation_function(new_genes)
#     end
# 
#     # add newly created individuals to the children population
#     push!(new_population.individuals, Individual(new_genes, newFitness))
#   end
