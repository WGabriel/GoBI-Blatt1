public class CDS {
	public String protein_id;
	public String transcript_id;
	public String transcript_name;
	public String gene_id;

	// Constructor
	public CDS(String protein_id, String transcript_id, String transcript_name, String gene_id) {
		super();
		this.protein_id = protein_id;
		this.transcript_id = transcript_id;
		this.transcript_name = transcript_name;
		this.gene_id = gene_id;
	}

	// Constructor for cloning
	public CDS(CDS another) {
		this.protein_id = another.protein_id;
		this.transcript_id = another.transcript_id;
		this.transcript_name = another.transcript_name;
		this.gene_id = another.gene_id;
	}

}
