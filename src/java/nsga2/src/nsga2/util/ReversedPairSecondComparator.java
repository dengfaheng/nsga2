package nsga2.util;

import java.util.Comparator;

public class ReversedPairSecondComparator<T> implements Comparator<Pair<T, Double>> {

    @Override
    public int compare(Pair<T, Double> first, Pair<T, Double> second) {
        return second.getSecond().compareTo(first.getSecond());
    }
}
