import java.util.HashSet;

public class OutputLine {
    String gene_id;
    String gene_symbol;
    String chromosome;
    String strand; // + or -
    int nprots; // number of annotated CDS in the gene
    int ntrans; // number of annotated transcripts in the gene
    Intron sv; // SV intron as start:end
    String wt; // WT introns within the SV intron separated by | as start:end
    HashSet<String> sv_prots; // ids of the SV CDS-s, separated by |
    HashSet<String> wt_prots; // ids of the WT CDS-s, separated by |
    int min_skipped_exon; // minimal number of skipped exons in any WT/SV pair
    int max_skipped_exon; // maximum number of skipped exons in any WT/SV pair
    int min_skipped_base; // min num of skipped bases (joint length of skipped exons) in any WT/SV pair
    int max_skipped_base; // max num of skipped bases (joint length of skipped exons) in any WT/SV pair

    //Constructor
    public OutputLine(String gene_id, String gene_symbol, String chromosome, String strand, int nprots, int ntrans,
                      Intron sv, String wt, HashSet<String> sv_prots, HashSet<String> wt_prots, int min_skipped_exon,
                      int max_skipped_exon, int min_skipped_base, int max_skipped_base) {
        super();
        this.gene_id = gene_id;
        this.gene_symbol = gene_symbol;
        this.chromosome = chromosome;
        this.strand = strand;
        this.nprots = nprots;
        this.ntrans = ntrans;
        this.sv = sv;
        this.wt = wt;
        this.sv_prots = sv_prots;
        this.wt_prots = wt_prots;
        this.min_skipped_exon = min_skipped_exon;
        this.max_skipped_exon = max_skipped_exon;
        this.min_skipped_base = min_skipped_base;
        this.max_skipped_base = max_skipped_base;
    }

    // Default Constructor
    public OutputLine() {
    }

}
