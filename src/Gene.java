import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

public class Gene {
    public String gene_id;
    // key = protein id, value = proteins
    public HashMap<String, Protein> proteins = new HashMap<>();
    // for ntrans
    public HashSet<String> transcripts = new HashSet<>();
    public TreeSet<Intron> introns = new TreeSet<>(new Intron());

    // Constructor
    public Gene(String gene_id) {
        super();
        this.gene_id = gene_id;
    }

    // Constructor
    public Gene(String gene_id, HashMap<String, Protein> proteins) {
        super();
        this.gene_id = gene_id;
        this.proteins = proteins;
    }

    // Constructor for cloning a Gene
    public Gene(Gene another) {
        this.gene_id = another.gene_id;
        //Deep copy of Protein list required
        for (HashMap.Entry<String, Protein> p : another.proteins.entrySet()) {
            this.proteins.put(p.getKey(), new Protein(p.getValue()));
        }
        for (Intron i : another.introns) {
            this.introns.add(new Intron(i));
        }
    }


}
