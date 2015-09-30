package nsga2;


import nsga2.interfaces.Crossover;
import nsga2.interfaces.Evaluate;
import nsga2.interfaces.Initialize;
import nsga2.interfaces.Mutate;
import nsga2.util.RngStream;
import nsga2.util.Util;

import java.security.InvalidParameterException;
import java.util.ArrayList;

public final class NSGA2 {





    public static HallOfFame nsga2(int populationSize,
                                   int numGenerations,
                                   int maxHallOfFameSize,
                                   Initialize initializeFunction,
                                   Evaluate evaluationFunction,
                                   Crossover crossoverFunction,
                                   Mutate mutate,
                                   long[] seeds) throws Exception{
        // parameter verification
        if (populationSize < 1)
        {
            throw new InvalidParameterException("population size must be > 0");
        }
        if (numGenerations < 1)
        {
            throw new InvalidParameterException("number of generations must be > 0");
        }
        if (maxHallOfFameSize < 1)
        {
            throw new InvalidParameterException("size of hall of fame must be > 0");
        }


        // initialize the RngStream
        RngStream stream = new RngStream();
        stream.setSeed(seeds);

        // initialize the hall of fame
        HallOfFame hallOfFame = new HallOfFame(maxHallOfFameSize, evaluationFunction);

        // initialize populations
        ArrayList<Individual> starter = new ArrayList<>();
        for (int i = 0; i != populationSize*2; ++i) {
            starter.add(initializeFunction.initialize());
        }
        Population previousPopulation = new Population(evaluationFunction);
        Population parentPopulation = new Population(evaluationFunction);
        Population mergedPopulation = new Population(starter, evaluationFunction);


        // main loop
        for (double generationIndex = 0; generationIndex != numGenerations; ++generationIndex)
        {
            // update the hall of fame
            ArrayList<Individual> tmp = new ArrayList<>();
            for (int index : mergedPopulation.getFronts().get(0))
            {
                tmp.add(new Individual(mergedPopulation.getIndividualList().get(index)));
            }
            hallOfFame.addFirstFront(tmp, stream);

            // find the individuals generating the next population
            // most of the time, only the last front must be selected by crowding distance
            // but you might need it on the first front
            int numToSelectFromLastFront;
            ArrayList<Integer> lastFrontIndices;
            ArrayList<ArrayList<Integer>> frontIndices;
            if (mergedPopulation.getFronts().get(0).size() >= populationSize)
            {
                //
                numToSelectFromLastFront =

            }
            else // no edge cases
            {//TODO refactor julia code...

                numToSelectFromLastFront = Util.nestedListSize(mergedPopulation.getFronts());
                ArrayList<Integer> lastFrontIndices = Population.lastFrontSelection(mergedPopulation,
                        Util.arrayListLast(mergedPopulation.getFronts()),
                        numToSelectFromLastFront,
                        stream);



            }


            // apply the
            ArrayList<Individual> selectedIndividuals = Population.uniqueFitnessTournamentSelection(parentPopulation, stream);
            nextPopulation = gener

        }


        return hallOfFame;

    }


}
