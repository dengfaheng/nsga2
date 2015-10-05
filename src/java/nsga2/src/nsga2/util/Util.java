package nsga2.util;


import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;





public final class Util {


    /**
     * choose a single element in the list, with randint(0, len(list))
     *
     * @param list
     * @param stream
     * @param <T>
     * @return
     */
    public static <T> T chooseRandom(List<T> list, RngStream stream) {
        return list.get(stream.randInt(0, list.size() - 1));
    }


    public static <T> ArrayList<T> flatten(List<List<T>> nested) {
        ArrayList<T> flat = new ArrayList<>();
        for (List<T> subList : nested) {
            flat.addAll(subList);
        }
        return flat;
    }


    /**
     * shuffle an ArrayList using the Fisher-Yates shuffle algorithm
     * with no side effect on the list given
     *
     * @param list   ArrayList
     * @param <T>    type for the ArrayList elements
     * @param stream mccons.RngStream object to get the random numbers
     * @return shuffled arrayList
     */
    public static <T> ArrayList<T> fisherYatesShuffle(List<T> list, RngStream stream) {
        ArrayList<T> copied = new ArrayList<>(list);
        int j;
        T tmp;
        for (int i = copied.size() - 1; i != 0; i--) {
            j = stream.randInt(0, i);
            // exchange list[i] with list[j] and vice versa
            tmp = copied.get(i);
            copied.set(i, copied.get(j));
            copied.set(j, tmp);
        }
        return copied;
    }


    public static <T> ArrayList<T> selectWithoutReplacement(List<T> list, int numToSelect, RngStream stream) {
        assert (numToSelect <= list.size());
        ArrayList<T> copied = fisherYatesShuffle(list, stream);
        return new ArrayList<>(copied.subList(0, numToSelect));
    }


    /**
     * checks wether or not the list given is sorted
     *
     * @param list
     * @param <T>
     * @return
     */
    public static <T extends Comparable<T>> boolean isSorted(List<T> list) {
        for (int i = 1; i < list.size(); ++i) {
            if (list.get(i - 1).compareTo(list.get(i)) > 0) {
                return false;
            }
        }
        return true;
    }


    /**
     * O(n) deletion, taking advantage that both lists are sorted
     *
     * @param sortedList
     * @param toDelete
     * @return
     */
    public static ArrayList<Integer> fastDelete(ArrayList<Integer> sortedList,
                                                 ArrayList<Integer> toDelete) {
        //take advantage of the knowledge that both vectors are sorted O(n)
        assert (Util.isSorted(sortedList));
        assert (Util.isSorted(toDelete));
        ArrayList<Integer> result = new ArrayList<>();

        if (toDelete.size() == 0)
        {
            return new ArrayList<>(sortedList);
        }

        int deleteIndex = 0;
        int deleteSize = toDelete.size();
        for (Integer i : sortedList) {

            // iterate to the next valid index, value>=to i
            while ((toDelete.get(deleteIndex) < i) && (deleteIndex + 1 < deleteSize)) {
                deleteIndex += 1;
            }

            if (toDelete.get(deleteIndex).compareTo(i) != 0) {
                result.add(i);
            }

        }
        return result;
    }


    public static class ReverseListComparator<T extends Comparable<T>> implements Comparator<List<T>> {
        private int index;

        /**
         * comparator used to reverse sort a list of ArrayList<T>
         * @param index_
         */
        public ReverseListComparator(int index_) {
            index = index_;
        }

        @Override
        public int compare(List<T> first, List<T> second) {
            return second.get(index).compareTo(first.get(index));
        }
    }


    public static <T> ArrayList<T> getIndices(ArrayList<T> list, ArrayList<Integer> indices) {
        ArrayList<T> subList = new ArrayList<>();
        for (int index : indices) {
            subList.add(list.get(index));
        }
        return subList;
    }


    public static Integer nestedListSize(ArrayList<ArrayList<Integer>> nestedLists)
    {
        Integer total = 0;
        for (ArrayList<Integer> nested : nestedLists)
        {
            total += nested.size();
        }
        return total;
    }


    public static Pair<Integer, Integer> select2(int low, int high, RngStream stream) {
        // select 2 different random integers in the interval
        assert (high - 1 > low);
        int first = stream.randInt(low, high - 1);
        int second = stream.randInt(low, high - 1);

        while (first == second) {
            second = stream.randInt(low, high - 1);
        }
        return new Pair(first, second);
    }


    /**
     * for each individual, map scores => list of index with same score
     *
     */
    public static <T> Hashtable<T, ArrayList<Integer>> mapElemsToIndices(ArrayList<T> list) {
        Hashtable<T, ArrayList<Integer>> valueToIndices = new Hashtable<>();
        int index = 0;
        for (T value : list) {
            if (!valueToIndices.containsKey(value)) {
                valueToIndices.put(value, new ArrayList<Integer>());
            }
            valueToIndices.get(value).add(index);
            index += 1;
        }
        return valueToIndices;
    }

}
