package nsga2;


import nsga2.util.Util;
import random.RngStream;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;


/**
 * holds the individuals and their performance information (front and crowding distance)
 */
public class Population {

    // placeholder for individuals, should not be modified (create new populations instead)
    private final ArrayList<Individual> individualList;
    // some mappings that are reused often (never recalculate)
    private Hashtable<ArrayList<Double>, ArrayList<Integer>> fitnessToIndices;
    private Hashtable<ArrayList<Double>, Integer> fitnessToFront;
    private Hashtable<ArrayList<Double>, Double> fitnessToCrowding;
    private ArrayList<ArrayList<Double>> uniqueFitness;


    Population(List<Individual> individualList_) {
        // copy individuals
        individualList = new ArrayList<>(individualList_);

        // initialize other structures empty
        fitnessToFront = new Hashtable<>();
        fitnessToCrowding = new Hashtable<>();

        // create the mappings
        mapUniqueFitnessToIndices();
    }

    Population(Population other) {
        this(other.individualList);
    }


    /**
     * find solutions dominating the solution at index of the population
     * the whole population (does not compare with itself)
     *
     * @param index index of the individual to compare against all others
     * @return index of the individuals dominating the solution at index
     */
    public ArrayList<Integer> findDominators(int index) {
        ArrayList<Integer> dominatedBy = new ArrayList<>();
        Individual individual = get(index);

        for (int i = 0; i != size(); ++i) {
            if (individual.compareTo(get(i)) < 0) {
                dominatedBy.add(i);
            }
        }
        return dominatedBy;
    }


    /**
     * for each individual, map scores => list of index with same score
     *
     */
    public void mapUniqueFitnessToIndices() {
        Hashtable<ArrayList<Double>, ArrayList<Integer>> mapping = new Hashtable<>();
        ArrayList<Double> fitness;
        for (int index = 0; index != size(); ++index) {
            fitness = get(index).getScores();
            if (!mapping.contains(fitness)) {
                mapping.put(fitness, new ArrayList<Integer>());
            }
            mapping.get(fitness).add(index);
        }
        fitnessToIndices = mapping;
        uniqueFitness = new ArrayList<>(mapping.keySet());
    }


    /**
     * crowding distance measures the proximity of a solution to its neighbors
     * it is used to preserve diversity, later in the algorithm
     */
    public void mapFitnessToCrowding(ArrayList<Integer> frontIndices) {
        ArrayList<Individual> front = Util.getIndices(getIndividualList(), frontIndices);

        // map fitness => crowding_distance
        ArrayList<ArrayList<Double>> fitnesses = new ArrayList<>();
        ArrayList<Double> fitness;
        for (Individual individual : front) {
            fitness = individual.getScores();
            if (!fitnesses.contains(fitness)) {
                fitnesses.add(fitness);
                fitnessToCrowding.put(fitness, 0.0);
            }
        }

        // sort in decreasing order the fitness vectors for each objective
        int fitnessLength = fitnesses.get(0).size();
        ArrayList<ArrayList<ArrayList<Double>>> sortedByObjectives = new ArrayList<>();
        ArrayList<Double> objectiveRanges = new ArrayList<>();
        for (int i = 0; i != fitnessLength; ++i) {
            Util.ReverseListComparator<Double> comparator = new Util.ReverseListComparator<>(i);
            Collections.sort(fitnesses, comparator);
            // now we can calculate the range of each objective
            sortedByObjectives.add(new ArrayList<>(fitnesses));
            objectiveRanges.add(fitnesses.get(0).get(i) - fitnesses.get(fitnesses.size() - 1).get(i));
        }

        // assign infinite crowding distance to maximum and minimum fitness of each objective
        ArrayList<Double> best, worse;
        for (ArrayList<ArrayList<Double>> list : sortedByObjectives) {
            // assign
            best = list.get(0);
            worse = list.get(list.size() - 1);
            fitnessToCrowding.put(best, Double.POSITIVE_INFINITY);
            fitnessToCrowding.put(worse, Double.POSITIVE_INFINITY);
        }

        // assign crowding crowding_distances to the other
        // fitness vectors for each objectives
        ArrayList<Double> currentFitness;
        double crowding, objectiveRange, fitPrevious, fitNext;
        for (int i = 0; i != fitnessLength; ++i) {
            // check for edge case (when all scores of the objective are the same)
            if (objectiveRanges.get(i) == 0.) {
                continue;
            } else {
                objectiveRange = objectiveRanges.get(i);
                for (int j = 2; j != fitnesses.size() - 1; ++j) {
                    // fetch the current values
                    currentFitness = sortedByObjectives.get(i).get(j);
                    crowding = fitnessToCrowding.get(currentFitness);

                    // crowding of n_1 is n_0 - n_2 / range_n (and infinity for
                    fitPrevious = sortedByObjectives.get(i).get(j - 1).get(i);
                    fitNext = sortedByObjectives.get(i).get(j + 1).get(i);
                    crowding += (fitPrevious - fitNext) / objectiveRange;

                    // update
                    fitnessToCrowding.put(currentFitness, crowding);
                }
            }
        }
    }




    /**
     * sort population into nondominating fronts (best to worst) until
     * at least half the original number of individuals is put in a front
     *
     */
    public void nonDominatedSort() {

        // for each solution, find its dominators
        ArrayList<Pair<Integer, ArrayList<Integer>>> indexToDominating = new ArrayList<>();
        for (int index = 0; index != size(); ++index) {
            indexToDominating.add(new Pair<>(index, findDominators(index)));
        }

        // some forward declarations
        // num of individuals we need
        int cutoff = (int) Math.ceil(size() / 2.0);
        int numSorted = 0;

        // indices of the non dominated individuals of the current front
        ArrayList<Integer> nonDominated;
        ArrayList<Pair<Integer, ArrayList<Integer>>> dominated;

        // find nondominated individuals and separate them from the rest iteratively
        ArrayList<ArrayList<Integer>> dominationFrontsIndices = new ArrayList<>();
        while (numSorted < cutoff) {
            dominated = new ArrayList<>();
            nonDominated = new ArrayList<>();

            // separate non-dominated from dominated
            for (Pair<Integer, ArrayList<Integer>> information : indexToDominating) {
                int index = information.getFirst();
                ArrayList<Integer> dominatedBy = information.getSecond();

                // non-dominated
                if (dominatedBy.size() == 0) {
                    nonDominated.add(index);
                }
                // dominated
                else {
                    dominated.add(new Pair<>(index, dominatedBy));
                }
            }

            // add the current front to the domination fronts
            dominationFrontsIndices.add(nonDominated);

            // update the domination information
            ArrayList<Pair<Integer, ArrayList<Integer>>> updatedIndexToDominating = new ArrayList<>();
            for (Pair<Integer, ArrayList<Integer>> information : dominated) {
                ArrayList<Integer> subtracted = Util.fastDelete(information.getSecond(), nonDominated);
                updatedIndexToDominating.add(new Pair<>(information.getFirst(), subtracted));
            }
            indexToDominating = updatedIndexToDominating;

            // update the count
            numSorted += nonDominated.size();
        }

        // assign into fitnessToFront
        for (int frontIndex = 0; frontIndex != dominationFrontsIndices.size(); ++frontIndex)
        {
            for (int index : dominationFrontsIndices.get(frontIndex))
            {
                fitnessToFront.put(get(index).getScores(), frontIndex);
            }
        }
    }




    /**
     * since individuals within the same front do not dominate each other, they are
     * # selected based crowding distance (greater diversity is desired)
     *
     * @param lastFrontIndices
     * @param numToSelect
     * @return
     */
    public ArrayList<Integer> lastFrontSelection(ArrayList<Integer> lastFrontIndices,
                                                 RngStream stream,
                                                 int numToSelect) {
        assert (0 < numToSelect) && (numToSelect < lastFrontIndices.size());

        // map fitness -> indices (one to many) and map fitness -> crowding distance
        Hashtable<ArrayList<Double>, ArrayList<Integer>> mapping = new Hashtable<>();
        ArrayList<Double> fitness;
        for (int index : lastFrontIndices) {
            fitness = get(index).getScores();
            if (!mapping.contains(fitness)) {
                mapping.put(fitness, new ArrayList<Integer>());
            }
            mapping.get(fitness).add(index);
        }

        // shuffle the indices so we don't need to later
        for (ArrayList<Double> fit : mapping.keySet()) {
            mapping.put(fit, Util.fisherYatesShuffle(mapping.get(fit), stream));

        }

        ArrayList<Pair<ArrayList<Double>, Double>> fitnessCrowdingPairs = new ArrayList<>();
        for (ArrayList<Double> key : fitnessToCrowding.keySet()) {
            fitnessCrowdingPairs.add(new Pair<>(key, fitnessToCrowding.get(key)));
        }

        // sort fitness by decreasing crowding distance
        ArrayList<ArrayList<Double>> orderedFitness = new ArrayList<>();
        ReversedPairSecondComparator<ArrayList<Double>> comparator = new ReversedPairSecondComparator<>();
        Collections.sort(fitnessCrowdingPairs, comparator);
        for (Pair<ArrayList<Double>, Double> pair : fitnessCrowdingPairs) {
            orderedFitness.add(pair.getFirst());
        }


        // choose individuals by iterating from best to worst (crowding distance)
        ArrayList<Integer> chosenIndices = new ArrayList<>();
        int position = 0;
        ArrayList<Integer> fitnessIndices;
        while (chosenIndices.size() < numToSelect) {
            fitness = orderedFitness.get(position);
            fitnessIndices = mapping.get(fitness);

            // add the last one (they are shuffled so it doesn't matter
            chosenIndices.add(fitnessIndices.get(fitnessIndices.size() - 1));
            fitnessIndices.remove(fitnessIndices.size() - 1);
            if (fitnessIndices.size() == 0) {
                mapping.remove(fitness);
                orderedFitness.remove(position);
            } else {
                mapping.put(fitness, fitnessIndices);
                position += 1;
            }

            // wrap around if the increment made it go out of bounds
            if (position > orderedFitness.size()) {
                position = 0;
            }

        }

        return chosenIndices;
    }




    /**
     * select across entire range of fitness to avoid bias by re-occurring fitness
     * @return new population
     */
    public Population uniqueFitnessTournamentSelection(RngStream stream)
    {
        //  edge case: only one fitness, return the population as it was
        if (uniqueFitness.size() == 1)
        {
            return new Population(this);
        }
        // else must select parents
        int popSize = size();
        ArrayList<Individual> selectedIndividuals = new ArrayList<>();

        int k, i, offset;
        ArrayList<ArrayList<Double>> candidateFitness, chosenFitness;
        ArrayList<Pair<Integer, Double>> frontAndCrowding;
        CrowdedCompare comparator = new CrowdedCompare();
        int front; double crowding;
        while (selectedIndividuals.size() != popSize)
        {
            // either pick all the fitness and select a random individual from them
            // or select a subset of them. depends on how many new parents still need to add
            k = Math.min(2 * (popSize - selectedIndividuals.size()), uniqueFitness.size());

            // sample k unique fitnesses and get their front and crowding information
            candidateFitness = Util.selectWithoutReplacement(uniqueFitness, k, stream);

            frontAndCrowding = new ArrayList<>();
            for (ArrayList<Double> fitness : candidateFitness)
            {
                crowding = fitnessToCrowding.get(fitness);
                front = fitnessToFront.get(fitness);
                frontAndCrowding.add(new Pair<Integer, Double>(front, crowding));
            }

            // choose the fitnesses
            chosenFitness = new ArrayList<>();
            i = 0;
            while (i < k)
            {
                // compare using the crowded compare operator (from Pair)
                offset = (comparator.compare(frontAndCrowding.get(i), frontAndCrowding.get(i+1)) > 0) ? 0 : 1;
                chosenFitness.add(candidateFitness.get(i+offset));
                i+=2;
            }

            // now randomly choose an individual from the indices associated with the chosen fitnesses
            ArrayList<Individual> chosenIndividuals = new ArrayList<>();
            int index;
            ArrayList<Integer> indicesWithSameFitness;
            for (ArrayList<Double> fitness : chosenFitness)
            {
                indicesWithSameFitness = fitnessToIndices.get(fitness);
                if (indicesWithSameFitness.size() > 1)
                {
                    index = Util.chooseRandom(indicesWithSameFitness, stream);
                    selectedIndividuals.add(new Individual(get(index)));
                }
            }

        }
        return new Population(selectedIndividuals);
    }


    Individual get(int index) {
        return individualList.get(index);
    }

    int size() {
        return individualList.size();
    }

    public ArrayList<Individual> getIndividualList() {
        return individualList;
    }


}
