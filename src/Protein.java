import java.util.TreeSet;

public class Protein {
    public String protein_id;
    // Not HashMap<String, CDS> because there's no unique identifier in CDS
    public TreeSet<CDS> cdss = new TreeSet<>(new CDS());

    // Constructor
    public Protein(String protein_id, TreeSet<CDS> cdss) {
        super();
        this.protein_id = protein_id;
        this.cdss = cdss;
    }

    // Constructor for cloning a Protein
    public Protein(Protein another) {
        this.protein_id = another.protein_id;
        //Deep copy of CDS list required
        for (CDS c : another.cdss) {
            this.cdss.add(new CDS(c));
        }
    }
}
