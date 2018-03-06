import java.util.Comparator;

public class Exon implements Comparator<Exon> {
    public String chr;
    public int start;
    public int end;
    public String strand;
    public String exon_id;
    public String transcript_id;
    public String transcript_name;
    public String gene_id;
    public String gene_name;

    // Constructor
    public Exon(String chr, int start, int end, String strand, String exon_id, String transcript_id,
                String transcript_name, String gene_id, String gene_name) {
        super();
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.exon_id = exon_id;
        this.transcript_id = transcript_id;
        this.transcript_name = transcript_name;
        this.gene_id = gene_id;
        this.gene_name = gene_name;
    }

    // Constructor for compare
    public Exon() {
        super();
    }

    // This constructor is used to create a deep copy of another exon
    public Exon(Exon another) {
        this.chr = another.chr;
        this.start = another.start;
        this.end = another.end;
        this.strand = another.strand;
        this.exon_id = another.exon_id;
        this.transcript_id = another.transcript_id;
        this.transcript_name = another.transcript_name;
        this.gene_id = another.gene_id;
        this.gene_name = another.gene_name;
    }

    public static void printExonInfo(Exon e) {
        System.out.println("|printExonInfo| exon_id: " + e.exon_id + " | transcript_id: " + e.transcript_id
                + " | transcript_name: " + e.transcript_name + " | gene_id: " + e.gene_id + " | gene_name: "
                + e.gene_name + " | chr: " + e.chr + " | start: " + e.start
                + " | end: " + e.end + " | strand: " + e.strand + " |");
    }

    // compare method not used, because we use no TreeSet<Exon>
    @Override
    public int compare(Exon e1, Exon e2) {
        // Sorts Exons ascending to their start-points
        return Integer.compare(e1.start, e2.start);
    }
}
