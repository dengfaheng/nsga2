package nsga2;


import nsga2.interfaces.Crossover;
import nsga2.interfaces.Evaluate;
import nsga2.interfaces.Mutate;
import nsga2.util.CrowdedCompare;
import nsga2.util.Pair;
import nsga2.util.ReversedPairSecondComparator;
import nsga2.util.Util;
import nsga2.util.RngStream;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;


/**
 * holds the individuals and their performance information (front and crowding distance)
 */
public class Population {

    // placeholder for individuals, should not be modified (create new populations instead)
    protected ArrayList<Individual> individualList;

    // some mappings that are reused often (never recalculate)
    protected Hashtable<ArrayList<Double>, Integer> fitnessToFront;
    protected Hashtable<ArrayList<Double>, Double> fitnessToCrowding;
    protected ArrayList<ArrayList<Integer>> fronts;

    // evaluation option
    protected Evaluate evaluationFunction;


    public Population(List<Individual> individualList_, Evaluate evaluationFunction_) {
        // copy individuals
        individualList = new ArrayList<>(individualList_);
        evaluationFunction = evaluationFunction_;
        calculateFronts();

        for (ArrayList<Integer> front: fronts)
        {
            for(Integer index : front)
            {
                assert(fitnessToCrowding.containsKey(individualList.get(index).getScores()));
            }
        }
    }


    Population(Population other) {
        this(other.individualList, other.getEvaluationFunction());
    }


    Population(Evaluate evaluationFunction) {
        this(new ArrayList<Individual>(), evaluationFunction);
    }


    public Evaluate getEvaluationFunction() {
        return evaluationFunction;
    }


    public void evaluatePopulation() {
        for (Individual individual : individualList) {
            individual.setScores(evaluationFunction.evaluate(individual.getGenes()));
        }
    }


    public Hashtable<ArrayList<Double>, Integer> getFitnessToFront() {
        return fitnessToFront;
    }

    public Hashtable<ArrayList<Double>, Double> getFitnessToCrowding() {
        return fitnessToCrowding;
    }

    public ArrayList<ArrayList<Integer>> getFronts() {
        return fronts;
    }

    /**
     * clean up the helper structures and recalculate them
     */
    public void calculateFronts() {
        // step 1: score the whole population
        evaluatePopulation();

        // map fitness -> fronts
        fronts = getNonDominatingFronts(this);

        fitnessToFront = new Hashtable<>();
        for (int front = 0; front != fronts.size(); ++front) {
            for (int index : fronts.get(front)) {
                ArrayList<Double> fitness = get(index).getScores();
                fitnessToFront.put(fitness, front);

            }
        }

        // map fitness -> crowding
        fitnessToCrowding = new Hashtable<>();
        for (ArrayList<Integer> frontIndices : fronts) {
            Hashtable<ArrayList<Double>, Double> frontCrowding = getFitnessToCrowding(this, frontIndices);
            for (ArrayList<Double> fitness : frontCrowding.keySet())
            {
                fitnessToCrowding.put(new ArrayList<>(fitness), new Double(frontCrowding.get(fitness)));
            }
        }



    }


    /**
     * find solutions dominating the solution at index of the population
     * the whole population (does not compare with itself)
     *
     * @param index index of the individual to compare against all others
     * @return index of the individuals dominating the solution at index
     */
    public static ArrayList<Integer> findDominators(ArrayList<Individual> individuals, int index) {
        ArrayList<Integer> dominatedBy = new ArrayList<>();
        Individual individual = individuals.get(index);

        for (int i = 0; i != individuals.size(); ++i) {
            if (individual.compareTo(individuals.get(i)) < 0) {
                dominatedBy.add(i);
            }
        }
        return dominatedBy;
    }


    /**
     * crowding distance measures the proximity of a solution to its neighbors
     * it is used to preserve diversity, later in the algorithm
     *
     * @param population   whole population
     * @param frontIndices indices of the individuals from the front
     * @return mapping fitness => crowding for every member of the front
     */
    public static Hashtable<ArrayList<Double>, Double> getFitnessToCrowding(Population population,
                                                                            ArrayList<Integer> frontIndices) {
        ArrayList<Individual> front = Util.getIndices(population.getIndividualList(), frontIndices);
        Hashtable<ArrayList<Double>, Double> fitnessToCrowdingMap = new Hashtable<>();

        // map fitness => crowding_distance
        ArrayList<ArrayList<Double>> uniqueFitness = new ArrayList<>();
        for (Individual individual : front) {
            ArrayList<Double> fitness = individual.getScores();
            if (!uniqueFitness.contains(fitness)) {
                uniqueFitness.add(fitness);
                fitnessToCrowdingMap.put(fitness, 0.0);
            }
        }

        // sort in decreasing order the fitness vectors for each objective
        int fitnessLength = uniqueFitness.get(0).size();
        ArrayList<ArrayList<ArrayList<Double>>> sortedByObjectives = new ArrayList<>();
        ArrayList<Double> objectiveRanges = new ArrayList<>();
        for (int i = 0; i != fitnessLength; ++i) {
            Util.ReverseListComparator<Double> comparator = new Util.ReverseListComparator<>(i);
            Collections.sort(uniqueFitness, comparator);
            // now we can calculate the range of each objective
            sortedByObjectives.add(new ArrayList<>(uniqueFitness));
            objectiveRanges.add(uniqueFitness.get(0).get(i) - uniqueFitness.get(uniqueFitness.size() - 1).get(i));
        }

        // assign infinite crowding distance to maximum and minimum fitness of each objective
        for (ArrayList<ArrayList<Double>> list : sortedByObjectives) {
            // assign
            ArrayList<Double> best = list.get(0);
            ArrayList<Double> worse = list.get(list.size() - 1);
            fitnessToCrowdingMap.put(best, Double.POSITIVE_INFINITY);
            fitnessToCrowdingMap.put(worse, Double.POSITIVE_INFINITY);
        }

        // assign crowding distances to the other fitness vectors for each objectives
        for (int i = 0; i != fitnessLength; ++i) {
            // check for edge case (when all scores of the objective are the same)
            if (objectiveRanges.get(i) != 0.) {
                double objectiveRange = objectiveRanges.get(i);
                for (int j = 2; j < uniqueFitness.size()-1; ++j) {
                    // fetch the current values
                    ArrayList<Double> currentFitness = sortedByObjectives.get(i).get(j);
                    double crowding = fitnessToCrowdingMap.get(currentFitness);

                    // crowding of n_1 is n_0 - n_2 / range_n (and infinity for edge cases)
                    double fitPrevious = sortedByObjectives.get(i).get(j - 1).get(i);
                    double fitNext = sortedByObjectives.get(i).get(j + 1).get(i);
                    crowding += (fitPrevious - fitNext) / objectiveRange;

                    // update
                    fitnessToCrowdingMap.put(currentFitness, crowding);
                }
            }
        }
        return fitnessToCrowdingMap;
    }


    /**
     * sort population into nondominating fronts (best to worst) until
     * at least half the original number of individuals is put in a front
     *
     * @param population population to sort (at least half) into-dominating fronts
     * @return non-dominating fronts containing at least half the size of the population
     */
    public static ArrayList<ArrayList<Integer>> getNonDominatingFronts(Population population) {

        ArrayList<ArrayList<Integer>> dominationFronts = new ArrayList<>();

        // for each solution, find its dominators
        ArrayList<Pair<Integer, ArrayList<Integer>>> indexToDominating = new ArrayList<>();
        for (int index = 0; index != population.size(); ++index) {
            indexToDominating.add(new Pair<>(index, findDominators(population.getIndividualList(), index)));
        }

        // some forward declarations
        // num of individuals we need
        int cutoff = (int) Math.ceil(population.size() / 2.0);

        // indices of the non dominated individuals of the current front

        // find non-dominated individuals and separate them from the rest iteratively
        while (Util.nestedListSize(dominationFronts) < cutoff) {

            ArrayList<Pair<Integer, ArrayList<Integer>>> dominated = new ArrayList<>();
            ArrayList<Integer> dominating = new ArrayList<>();

            // separate non-dominated from dominated
            for (Pair<Integer, ArrayList<Integer>> information : indexToDominating) {
                // non-dominated
                if (information.getSecond().size() == 0) {
                    dominating.add(information.getFirst());
                }
                // dominated
                else {
                    dominated.add(new Pair<>(information.getFirst(), new ArrayList<>(information.getSecond())));
                }
            }

            // add the current front to the domination fronts
            dominationFronts.add(new ArrayList<>(dominating));

            // update the domination information
            ArrayList<Pair<Integer, ArrayList<Integer>>> updatedIndexToDominating = new ArrayList<>();
            for (Pair<Integer, ArrayList<Integer>> information : dominated) {
                ArrayList<Integer> subtracted = Util.fastDelete(information.getSecond(), dominating);
                updatedIndexToDominating.add(new Pair<>(information.getFirst(), subtracted));
            }
            indexToDominating = updatedIndexToDominating;
        }

        // update member variables related to non domination fronts
        return dominationFronts;
    }


    /**
     * since individuals within the same front do not dominate each other, they are
     * # selected based crowding distance (greater diversity is desired)
     *
     * @param numToSelect number of individuals to choose
     * @return indices of the individuals selected
     */
    public static ArrayList<Integer> lastFrontSelection(Population population,
                                                        ArrayList<Integer> lastFrontIndices,
                                                        int numToSelect,
                                                        RngStream stream) {
        assert (0 < numToSelect) && (numToSelect < lastFrontIndices.size());

        // map fitness -> indices (one to many) and map fitness -> crowding distance
        Hashtable<ArrayList<Double>, ArrayList<Integer>> fitnessToIndices = new Hashtable<>();
        for (int index : lastFrontIndices) {
            ArrayList<Double> fitness = population.get(index).getScores();
            if (!fitnessToIndices.containsKey(fitness)) {
                fitnessToIndices.put(fitness, new ArrayList<Integer>());
            }
            ArrayList<Integer> current = fitnessToIndices.get(fitness);
            current.add(index);
            fitnessToIndices.put(fitness, new ArrayList<>(current));
        }

        // shuffle the indices so we don't need to later
        for (ArrayList<Double> fit : fitnessToIndices.keySet()) {
            fitnessToIndices.put(fit, Util.fisherYatesShuffle(fitnessToIndices.get(fit), stream));

        }

        ArrayList<Pair<ArrayList<Double>, Double>> fitnessToCrowdingMap = new ArrayList<>();
        for (ArrayList<Double> key : fitnessToIndices.keySet()) {
            fitnessToCrowdingMap.add(new Pair<>(key, population.getFitnessToCrowding().get(key)));
        }

        // sort fitness by decreasing crowding distance
        ArrayList<ArrayList<Double>> orderedFitness = new ArrayList<>();
        ReversedPairSecondComparator<ArrayList<Double>> comparator = new ReversedPairSecondComparator<>();
        Collections.sort(fitnessToCrowdingMap, comparator);
        for (Pair<ArrayList<Double>, Double> pair : fitnessToCrowdingMap) {
            orderedFitness.add(pair.getFirst());
        }


        // choose individuals by iterating from best to worst (crowding distance)
        int position = 0;
        ArrayList<Integer> chosenIndices = new ArrayList<>();
        while (chosenIndices.size() < numToSelect) {
            ArrayList<Double> fitness = orderedFitness.get(position);
            ArrayList<Integer> indices = fitnessToIndices.get(fitness);

            // add the last one (they are shuffled so it doesn't matter)
            int index = indices.get(indices.size() -1);
            chosenIndices.add(index);
            indices.remove(indices.size() - 1);
            if (indices.size() == 0) {
                orderedFitness.remove(position);
            } else {
                fitnessToIndices.put(fitness, indices);
                position += 1;
            }

            // wrap around if the increment made it go out of bounds
            if (position >= orderedFitness.size()) {
                position = 0;
            }

        }

        return chosenIndices;
    }


    /**
     * select across entire range of fitness to avoid bias by re-occurring fitness
     *
     * @return new population
     */
    public static ArrayList<Individual> uniqueFitnessTournamentSelection(Population population, RngStream stream) {
        //  edge case: only one fitness, return the population as it was

        //TODO debug this shit
        ArrayList<ArrayList<Double>> scores = new ArrayList<>();
        for (Individual individual : population.getIndividualList()) {
            scores.add(individual.getScores());
        }
        Hashtable<ArrayList<Double>, ArrayList<Integer>> fitnessToIndices = Util.mapElemsToIndices(scores);
        ArrayList<ArrayList<Double>> uniqueFitness = new ArrayList<>(fitnessToIndices.keySet());
        if (uniqueFitness.size() == 1) {
            return population.getIndividualList();
        }
        // else must select parents
        int popSize = population.size();
        ArrayList<Individual> selectedIndividuals = new ArrayList<>();


        int k, i, offset;


        CrowdedCompare comparator = new CrowdedCompare();

        while (selectedIndividuals.size() != popSize) {
            // either pick all the fitness and select a random individual from them
            // or select a subset of them. depends on how many new parents still need to add
            k = Math.min(2 * (popSize - selectedIndividuals.size()), uniqueFitness.size());

            // sample k unique fitness and get their front and crowding information
            ArrayList<ArrayList<Double>> candidateFitness = Util.selectWithoutReplacement(uniqueFitness, k, stream);

            frontAndCrowding = new ArrayList<>();
            for (ArrayList<Double> fitness : candidateFitness) {
                crowding = population.getFitnessToCrowding().get(fitness);
                front = population.getFitnessToFront().get(fitness);
                frontAndCrowding.add(new Pair<>(front, crowding));
            }

            // choose the fitness
            chosenFitness = new ArrayList<>();
            i = 0;
            while (i < k) {
                // compare using the crowded compare operator (from Pair)
                offset = (comparator.compare(frontAndCrowding.get(i), frontAndCrowding.get(i + 1)) > 0) ? 0 : 1;
                chosenFitness.add(candidateFitness.get(i + offset));
                i += 2;
            }

            // now randomly choose an individual from the indices associated with the chosen fitness
            ArrayList<Integer> indicesWithSameFitness;
            for (ArrayList<Double> fitness : chosenFitness) {
                indicesWithSameFitness = fitnessToIndices.get(fitness);
                if (indicesWithSameFitness.size() > 1) {
                    int index = Util.chooseRandom(indicesWithSameFitness, stream);
                    selectedIndividuals.add(new Individual(population.get(index)));
                }
            }

        }
        return selectedIndividuals;
    }


    Population generateOffsprings(Mutate mutate, Crossover crossover, RngStream stream) {
        ArrayList<Individual> newIndividuals = new ArrayList<>();
        Pair<Integer, Integer> parentIndices;
        for (int i = 0; i != individualList.size(); ++i) {
            parentIndices = Util.select2(0, individualList.size(), stream);

            // crossover
            Individual offspring = crossover.crossover(individualList.get(parentIndices.getFirst()), individualList.get(parentIndices.getSecond()));

            // mutate
            offspring = mutate.mutate(offspring);

            newIndividuals.add(offspring);
        }
        return new Population(newIndividuals, evaluationFunction);
    }


    Individual get(int index) {
        return individualList.get(index);
    }

    int size() {
        return individualList.size();
    }

    public ArrayList<Individual> getIndividualList() {
        return new ArrayList<>(individualList);
    }

    public ArrayList<Individual> getIndividualAtIndices(ArrayList<Integer> indices) {
        ArrayList<Individual> copied = new ArrayList<>();
        for (Integer index : indices) {
            copied.add(new Individual(individualList.get(index)));
        }
        return copied;
    }

}
