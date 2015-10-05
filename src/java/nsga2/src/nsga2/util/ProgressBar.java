package nsga2.util;

public class ProgressBar {
    private String processName;
    private int barLength;


    public ProgressBar(String processName_, int bar_length_) {
        processName = processName_;
        barLength = bar_length_;
    }


    public void update(double percentage) {
        // displays the current progress of the process
        assert (percentage >= 0.);
        assert (percentage <= 1.);

        int progressLength = (int) (barLength * percentage);
        System.out.print("[");
        for (int i = 0; i < barLength; ++i) {
            if (i <= progressLength)
                System.out.print("#");
            else
                System.out.print(" ");
        }
        System.out.print("] " + (int)(percentage * 100.) + "% " + processName + "\r");
        System.out.flush();
    }


    /**
     *  called after last iteration, displays 100% success
     */
    public void clean() {
        this.update(1.);
        System.out.println();
        System.out.flush();
    }
}