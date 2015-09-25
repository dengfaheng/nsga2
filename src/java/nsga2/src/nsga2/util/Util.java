package nsga2.util;


import nsga2.random.RngStream;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;





public final class Util {


    /**
     * choose a single element in the list, with randint(0, len(list))
     * @param list
     * @param stream
     * @param <T>
     * @return
     */
    public static <T> T chooseRandom(List<T> list, RngStream stream)
    {
        return list.get(stream.randInt(0, list.size()-1));
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
     * @param sortedIndicesToDelete
     * @return
     */
    public static ArrayList<Integer> fastDelete(ArrayList<Integer> sortedList,
                                                ArrayList<Integer> sortedIndicesToDelete) {
        assert (isSorted(sortedList) && isSorted(sortedIndicesToDelete));
        ArrayList<Integer> filtered = new ArrayList<>();
        int deletionIndex = 0;
        final int numToDelete = sortedIndicesToDelete.size();
        for (int i : sortedList) {
            // catch up
            while ((sortedIndicesToDelete.get(deletionIndex) < i) && (deletionIndex < numToDelete)) {
                deletionIndex += 1;
            }
            if (i != sortedIndicesToDelete.get(deletionIndex)) {
                filtered.add(i);
            }

        }
        return filtered;
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
}
