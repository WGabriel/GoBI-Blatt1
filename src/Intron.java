import java.util.Comparator;

public class Intron implements Comparator<Intron> {
    public int start;
    public int end;

    // Constructor
    public Intron(int start, int end) {
        super();
        this.start = start;
        this.end = end;
    }

    // Constructor for cloning
    public Intron(Intron another) {
        this.start = another.start;
        this.end = another.end;
    }

    public Intron() {
//        super();
    }

    @Override
    public int compare(Intron intr1, Intron intr2) {
        // Sorts CDS ascending to their start-points
        if (intr1.start > intr2.start) {
            return 1;
        } else if (intr1.start < intr2.start) {
            return -1;
        } else {
            return Integer.compare(intr1.end, intr2.end);
        }
    }

}
