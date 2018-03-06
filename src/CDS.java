import java.util.Comparator;

public class CDS implements Comparator<CDS> {
    public String chr;
    public int start;
    public int end;
    public String strand;
    public String protein_id;
    public String transcript_id;
    public String transcript_name;
    public String gene_id;
    public String gene_name;

    // Constructor
    public CDS(String chr, int start, int end, String strand, String protein_id, String transcript_id, String transcript_name, String gene_id, String gene_name) {
        super();
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.protein_id = protein_id;
        this.transcript_id = transcript_id;
        this.transcript_name = transcript_name;
        this.gene_id = gene_id;
        this.gene_name = gene_name;
    }

    // Constructor for cloning
    public CDS(CDS another) {
        this.chr = another.chr;
        this.start = another.start;
        this.end = another.end;
        this.strand = another.strand;
        this.protein_id = another.protein_id;
        this.transcript_id = another.transcript_id;
        this.transcript_name = another.transcript_name;
        this.gene_id = another.gene_id;
        this.gene_name = another.gene_name;
    }

    // For new TreeSet<>(new CDS()); in Gene-class
    public CDS() {
    }

    @Override
    public int compare(CDS cds1, CDS cds2) {
        // Sorts CDS ascending to their start-points
        if (cds1.start > cds2.start) {
            return 1;
        } else if (cds1.start < cds2.start) {
            return -1;
        } else {
            if (cds1.end > cds2.end) {
                return 1;
            } else if (cds1.end < cds2.end) {
                return -1;
            } else {
                return cds1.protein_id.compareTo(cds2.protein_id);
            }
        }
    }

}
