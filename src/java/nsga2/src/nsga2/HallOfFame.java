package nsga2;

import nsga2.interfaces.Evaluate;
import nsga2.util.RngStream;
import nsga2.util.Util;

import java.util.ArrayList;
import java.util.Hashtable;

public class HallOfFame extends Population{
    private final int maxSize;

    HallOfFame(int maxSize_, Evaluate evaluationFunction)
    {
        super(new ArrayList<Individual>(), evaluationFunction);
        maxSize = maxSize_;
        return;
    }


    public void addFirstFront(ArrayList<Individual> firstFront, RngStream stream)
    {
        // add the individuals
        fitnessToFront = null;
        individualList.addAll(firstFront);

        // map fitness -> fronts
        fronts = getNonDominatingFronts(this);
        individualList = Util.getIndices(individualList, fronts.get(0));

        if (individualList.size() > maxSize)
        {
            // filter out by crowding distance
            ArrayList<ArrayList<Double>> scores = new ArrayList<>();
            for (Individual individual : individualList)
            {
                scores.add(individual.getScores());
            }
            fitnessToCrowding = new Hashtable<>();
            fitnessToCrowding.putAll(getFitnessToCrowding(this, fronts.get(0)));
            individualList = Util.getIndices(individualList, lastFrontSelection(this, fronts.get(0), maxSize, stream));
        }
    }
}
