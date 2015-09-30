package nsga2.interfaces;


import nsga2.Individual;

public interface Crossover {
    Individual crossover(Individual first, Individual second);

}

