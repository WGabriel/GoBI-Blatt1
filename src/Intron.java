public class Intron {
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
}
