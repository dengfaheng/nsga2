package nsga2.interfaces;

import nsga2.Individual;

import java.util.ArrayList;

public interface Evaluate {
    ArrayList<Double> evaluate(ArrayList<Integer> genes);
}
