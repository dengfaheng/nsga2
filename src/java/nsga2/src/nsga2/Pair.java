package nsga2;

import java.util.Comparator;

public class Pair<T1, T2> {

    private T1 first;
    private T2 second;

    public Pair(T1 first_, T2 second_) {
        first = first_;
        second = second_;
    }

    public T1 getFirst() {
        return first;
    }

    public T2 getSecond() {
        return second;
    }
}

class ReversedPairSecondComparator<T> implements Comparator<Pair<T, Double>> {

    @Override
    public int compare(Pair<T, Double> first, Pair<T, Double> second) {
        return second.getSecond().compareTo(first.getSecond());
    }
}


class CrowdedCompare implements Comparator<Pair<Integer, Double>> {

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
