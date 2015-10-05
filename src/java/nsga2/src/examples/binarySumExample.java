package examples;

import nsga2.Individual;
import nsga2.NSGA2;
import nsga2.interfaces.Crossover;
import nsga2.interfaces.Evaluate;
import nsga2.interfaces.Initialize;
import nsga2.interfaces.Mutate;
import nsga2.util.RngStream;

import java.lang.reflect.Array;
import java.util.ArrayList;

public class binarySumExample {


    /**
     * very simple crossover accross binary genome
     */
    static class SimpleCrossover implements Crossover
    {
        private double applicationProbability;
        private double mixingRatio;
        private RngStream stream;

        public SimpleCrossover(double applicationProbability, double mixingRatio, RngStream stream) {
            this.applicationProbability = applicationProbability;
            this.mixingRatio = mixingRatio;
            this.stream = stream;
        }

        public Individual crossover(Individual first, Individual second)
        {
            if (stream.randU01() < applicationProbability)
            {
                // apply the crossover
                ArrayList<Integer> newGenes = new ArrayList<>();
                for(int i = 0; i != first.getGenes().size(); ++i)
                {
                    if (stream.randU01() < mixingRatio)
                    {
                        // take from second
                        newGenes.add(second.getGenes().get(i));
                    }
                    else
                    {
                        // take from first
                        newGenes.add(first.getGenes().get(i));
                    }
                }
                return new Individual(newGenes);
            }
            else
            {
                return new Individual(first);
            }
        }
    }


    static class binSumEval implements Evaluate
    {
        @Override
        public ArrayList<Double> evaluate(ArrayList<Integer> genes) {
            double sum = 0;
            for (Integer gene: genes)
            {
                sum += gene;
            }

            double decal = 0;
            for (int i = 0; i != genes.size()-1; ++i)
            {
                if(genes.get(i) != genes.get(i+1))
                {
                    decal += 1;
                }
            }

            ArrayList<Double> scores = new ArrayList<>();
            scores.add(sum);
            scores.add(decal);
            return scores;
        }
    }




    static class SimpleMutate implements Mutate
    {
        private double applicationProbability;
        private double mutationProbability;
        private RngStream stream;

        public SimpleMutate(double applicationProbability, double mutationProbability, RngStream stream) {
            this.applicationProbability = applicationProbability;
            this.mutationProbability = mutationProbability;
            this.stream = stream;
        }

        private Integer flip(Integer gene)
        {
            if (gene.compareTo(1) == 0)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }

        public Individual mutate(Individual individual)
        {
            if (stream.randU01() > applicationProbability)
            {
                // return individual as is
                return new Individual(individual);
            }
            else {
                ArrayList<Integer> newGenes = new ArrayList<>();
                // apply mutation
                for (int i = 0; i != individual.getGenes().size(); ++i) {
                        if (stream.randU01() < mutationProbability) {
                            newGenes.add(flip(individual.getGenes().get(i)));
                    }
                }
                return new Individual(newGenes);
            }
        }
    }


    static class BinaryInitializer implements Initialize
    {
        private int geneSize;
        private RngStream stream;

        public BinaryInitializer(int geneSize, RngStream stream) {
            this.geneSize = geneSize;
            this.stream = stream;
        }

        @Override
        public Individual initialize() {
            ArrayList<Integer> genes = new ArrayList<>();
            for (int i = 0; i != geneSize; ++i)
            {
                genes.add(stream.randInt(0, 1));
            }
            return new Individual(genes);
        }
    }



    public static void main(String[] args)
    {
        // parameters


        RngStream stream = new RngStream();

        double mutApp = 0.1;
        double mutProb = 0.05;
        double crossApp = 0.3;
        double mixingRatio= 0.2;
        int geneSize = 10;
        int popSize = 10;
        int numGenerations = 10;
        int maxHOFSize = popSize;
        long[] seeds = {42, 42, 42, 42, 42, 42};

        binSumEval evaluate = new binSumEval();
        BinaryInitializer initializer = new BinaryInitializer(geneSize, stream);
        SimpleCrossover crossover = new SimpleCrossover(crossApp, mixingRatio, stream);
        SimpleMutate mutate = new SimpleMutate(mutApp, mutProb, stream);

        NSGA2 optimizer = new NSGA2();

        try {
            optimizer.nsga2(popSize, numGenerations, maxHOFSize, initializer,evaluate, crossover, mutate, seeds);
        } catch (Exception e) {
            e.printStackTrace();
        }


    }
}
