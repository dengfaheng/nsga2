package nsga2;


import nsga2.genops.Evaluate;
import nsga2.random.RngStream;

import java.util.ArrayList;

public final class NSGA2 {


    public static HallOfFame nsga2(int populationSize,
                                              int numGenerations,
                                              Evaluate evaluateGenes,
                                              int maxHallOfFameSize,
                                              long[] seeds){

        // initialize the RngStream
        RngStream stream = new RngStream();
        stream.setSeed(seeds);

        // initialize the hall of fame
        HallOfFame hallOfFame = new HallOfFame(maxHallOfFameSize);


        //
        ArrayList<Individual> individuals = new ArrayList<>();


        return hallOfFame;

    }


}
