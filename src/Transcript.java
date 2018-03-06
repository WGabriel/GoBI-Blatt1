import java.util.HashMap;
import java.util.Map;

public class Transcript {
    public String transcript_id;
    // String = exon Ids
    public HashMap<String, Exon> exons = new HashMap<>();

    // Constructor
    public Transcript(String transcript_id, HashMap<String, Exon> exons) {
        super();
        this.exons = exons;
        this.transcript_id = transcript_id;
    }

    // Constructor for cloning a Transcript
    public Transcript(Transcript another) {
        //Deep copy of Exon list required
        for (HashMap.Entry<String, Exon> e : another.exons.entrySet()) {
            this.exons.put(e.getKey(), new Exon(e.getValue()));
        }
        this.transcript_id = another.transcript_id;
    }

    public static void printTranscriptInfo(Transcript t) {
        System.out.print("|printTranscriptInfo| transcript_id: " + t.transcript_id + " contains exons: ");
        for (Map.Entry<String, Exon> e : t.exons.entrySet()) {
            System.out.print(e.getKey() + " | ");
        }
        System.out.println();
    }

}
