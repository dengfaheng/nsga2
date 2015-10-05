package nsga2;



import nsga2.interfaces.Crossover;
import nsga2.interfaces.Evaluate;
import nsga2.interfaces.Initialize;
import nsga2.interfaces.Mutate;
import nsga2.util.ProgressBar;
import nsga2.util.RngStream;


import java.security.InvalidParameterException;
import java.util.ArrayList;

public final class NSGA2 {


    /**
     * NSGA-II procedure, as described in fig 7 of
     * "Revisiting the NSGA-II Crowding-Distance Computation"
     * @param N size of populations
     * @param NGEN number of main loop iterations to perform (number of generations)
     * @param maxHallOfFameSize size limit of the container of the best-ever solutions
     * @param initializeFunction function to initialize a new individual
     * @param evaluationFunction function to evaluate a new individual
     * @param crossoverFunction function to perform crossover on two individuals
     * @param mutate function to mutate a single individual
     * @param seeds seeds of the pseudo-random number generator (used to make it replicable)
     * @return Hall of Fame, the best solutions seen throughout all generations
     * @throws Exception read exception on the input file for now
     */
    public static HallOfFame nsga2(int N,
                                   int NGEN,
                                   int maxHallOfFameSize,
                                   Initialize initializeFunction,
                                   Evaluate evaluationFunction,
                                   Crossover crossoverFunction,
                                   Mutate mutate,
                                   long[] seeds) throws Exception{
        // parameter verification
        if (N < 1)
        {
            throw new InvalidParameterException("population size must be > 0");
        }
        if (NGEN < 1)
        {
            throw new InvalidParameterException("number of generations must be > 0");
        }
        if (maxHallOfFameSize < 1)
        {
            throw new InvalidParameterException("size of hall of fame must be > 0");
        }


        System.out.println("initializing populations");
        // initialize the RngStream
        RngStream stream = new RngStream();
        stream.setSeed(seeds);

        // initialize the hall of fame
        HallOfFame hallOfFame = new HallOfFame(maxHallOfFameSize, evaluationFunction);

        // initialize populations
        ArrayList<Individual> starter1 = new ArrayList<>();
        ArrayList<Individual> starter2 = new ArrayList<>();
        for (int i = 0; i != N; ++i) {
            starter1.add(initializeFunction.initialize());
            starter2.add(initializeFunction.initialize());
        }
        System.out.println("creating parent populations");
        Population parents = new Population(starter1, evaluationFunction);
        Population offsprings = new Population(starter2, evaluationFunction);
        Population merged;


        // main loop: |selection -> generation| -> |selection -> generation| -> ...
        ProgressBar bar = new ProgressBar("evo loop", 40);
        System.out.println("main loop starting");
        for (double generationIndex = 0; generationIndex != NGEN; ++generationIndex)
        {

            bar.update(generationIndex / NGEN);

            // perform union of individuals from P and Q
            ArrayList<Individual> mergedIndividuals = new ArrayList<>();
            mergedIndividuals.addAll(parents.getIndividualList());
            mergedIndividuals.addAll(offsprings.getIndividualList());

            // initialize merged population and empty the parents and offsprings
            // performs front sorting
            // performs crowding distance calculations
            merged = new Population(mergedIndividuals, evaluationFunction);


            // STEP 1.5: UPDATE HALL OF FAME
            hallOfFame.addFirstFront(merged.getIndividualAtIndices(merged.getFronts().get(0)), stream);


            // STEP 2: SELECT PARENTS
            // first fronts are selected without using crowding distance
            ArrayList<Integer> selectedIndices = new ArrayList<>();
            for (int frontIndex = 0; frontIndex != merged.getFronts().size()-1; ++frontIndex)
            {
                selectedIndices.addAll(merged.getFronts().get(frontIndex));
            }
            ArrayList<Integer> lastFront = merged.getFronts().get(merged.getFronts().size() - 1);
            int k = N - selectedIndices.size();
            assert (k >= 0);

            // add from last front, using crowdingd distance
            selectedIndices.addAll(Population.lastFrontSelection(merged, lastFront, k, stream));
            Population P_next = new Population(merged.getIndividualAtIndices(selectedIndices),evaluationFunction);


            // STEP 3: CREATE NEW GENERATION FROM PARENTS
            parents = new Population(Population.uniqueFitnessTournamentSelection(P_next, stream), evaluationFunction);
            offsprings = parents.generateOffsprings(mutate, crossoverFunction, stream);

        }
        bar.clean();


        return hallOfFame;

    }


}
