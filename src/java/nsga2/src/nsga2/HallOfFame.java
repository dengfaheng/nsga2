package nsga2;

import nsga2.random.RngStream;
import nsga2.util.Util;

import java.util.ArrayList;
import java.util.Hashtable;

public class HallOfFame extends Population{
    private final int maxSize;

    HallOfFame(int maxSize_)
    {
        super(new ArrayList<Individual>());
        maxSize = maxSize_;
        return;
    }


    public void addFirstFront(ArrayList<Individual> firstFront, RngStream stream)
    {
        // add the individuals
        fitnessToFront = null;
        individualList.addAll(firstFront);

        // map fitness -> indices

        // map fitness -> fronts
        fronts = getNondominatingFronts(this);
        individualList = Util.getIndices(individualList, fronts.get(0));


        if (individualList.size() > maxSize)
        {
            // filter out by crowding distance
            fitnessToIndices = getUniqueFitnessToIndices(this);
            fitnessToCrowding = new Hashtable<>();
            fitnessToCrowding.putAll(getFitnessToCrowding(this, fronts.get(0)));
            individualList = Util.getIndices(individualList, lastFrontSelection(this, fronts.get(0), fitnessToCrowding, stream, maxSize));
        }
    }
}
