package nsga2.util;

import java.util.Comparator;

public class CrowdedCompare implements Comparator<Pair<Integer, Double>> {

    @Override
    public int compare(Pair<Integer, Double> first, Pair<Integer, Double> second) {
        // if rank is the same, tie break with crowding distance
        int compareFronts = first.getFirst().compareTo(second.getFirst());
        if (compareFronts == 0) {
            return first.getSecond().compareTo(second.getSecond());
        } else {
            return compareFronts;
        }
    }
}
