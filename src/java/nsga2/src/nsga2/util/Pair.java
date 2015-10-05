package nsga2.util;

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



