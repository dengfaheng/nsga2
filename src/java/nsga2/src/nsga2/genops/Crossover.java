package nsga2.genops;


import nsga2.Individual;

public interface Crossover {
    Individual crossover(Individual first, Individual second);

}

