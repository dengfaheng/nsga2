package nsga2;

import com.sun.istack.internal.NotNull;

import java.util.ArrayList;



public class Individual implements Comparable<Individual> {

    private final ArrayList<Integer> genes;
    private ArrayList<Double> scores;

    /**
     * Constructor
     *
     * @param genes_  index of the chosen genes in a table
     * @param scores_ scores of the current solution configuration (based on genes)
     */
    public Individual(ArrayList<Integer> genes_,
                      ArrayList<Double> scores_) {
        genes = genes_;
        scores = scores_;
    }


    public Individual(ArrayList<Integer> genes)
    {
        this.genes = genes;
        scores = new ArrayList<>();
    }


    /**
     * Copy constructor
     *
     * @param other solution to copy
     */
    public Individual(Individual other) {
        this(other.getGenes(), other.getScores());
    }


    public ArrayList<Integer> getGenes() {
        return genes;
    }


    public ArrayList<Double> getScores() {
        return scores;
    }

    public void setScores(ArrayList<Double> newScores)
    {
        scores = newScores;
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();

        // add the genes
        builder.append("genes: ");
        for (Integer i : getGenes()) {
            builder.append(i);
            builder.append(" ");
        }
        builder.append(System.lineSeparator());

        // add the score
        builder.append("score: ");
        for (Double score : getScores()) {
            builder.append(score);
            builder.append(" ");
        }
        builder.append(System.lineSeparator());

        return builder.toString();
    }


    /*
     * used to tell whether first vector is dominated/dominating/non-dominating the second vector
     * used to compare solutions by their scores, not identity
     * ([0, 0, 2]  > [0, 0, 1]) =  1
     * ([0, 0, 1] == [0, 1, 0]) =  0
     * ([1, 0, 1]  < [1, 1, 1]) = -1
     */
    @Override
    public int compareTo(@NotNull Individual other) {
        assert (scores.size() == other.getScores().size());

        boolean firstDominates = false;
        boolean secondDominates = false;
        Double i,j;
        int compare;

        for (int index = 0; index != scores.size(); ++index) {
            i = scores.get(index);
            j = other.getScores().get(index);
            compare = i.compareTo(j);
            if (compare < 0) {
                secondDominates = true;
            } else if (compare > 0) {
                firstDominates = true;
            }

            if (firstDominates && secondDominates) {
                return 0;
            }
        }

        if (firstDominates) {
            return 1;
        } else if (secondDominates) {
            return -1;
        } else {
            return 0;
        }
    }

}

