import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.Map.Entry;

public class ExonSkipping {

    public static void main(String[] args) {
        // Check if console parameters are correct
        if (args.length != 4) {
            System.err.println("Too few or many ConsoleParameters.");
            System.exit(1);
        }
        File gtfInput = null;
        File tsvOutput = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-o")) {
                tsvOutput = new File(args[i + 1]);
            }
            if (args[i].equals("-gtf")) {
                gtfInput = new File(args[i + 1]);
            }
        }
        assert gtfInput != null;
        assert tsvOutput != null;
        System.out.println("gtf input file:\t" + gtfInput.getAbsolutePath());
        System.out.println("tsv output file:\t" + tsvOutput.getAbsolutePath());

        HashMap<String, Gene> allGenes = parseGtf(gtfInput);

        HashMap<String, OutputLine> outputLines = getOutputNew(allGenes);

        writeOutput(outputLines, tsvOutput);

        System.out.println("End of main-method.");
    }

    private static HashMap<String, OutputLine> getOutputNew(HashMap<String, Gene> allGenes) {
        // K = ENSG00000131018\t6\t152688518:152690606, V = OutputLine
        HashMap<String, OutputLine> result = new HashMap<>();
        for (Gene gene : allGenes.values()) {
            // System.out.println("|getOutputNew| Gene: " + gene.gene_id + " contains " + gene.cds.size() + " CDS's.");
            // sort CDS's by protein_id
            HashMap<String, ArrayList<CDS>> cdsByProteinId = new HashMap<>();
            for (CDS cds : gene.cds) {
                if (cdsByProteinId.containsKey(cds.protein_id)) {
                    cdsByProteinId.get(cds.protein_id).add(new CDS(cds));
                } else {
                    ArrayList<CDS> cdsList = new ArrayList<>();
                    cdsList.add(new CDS(cds));
                    cdsByProteinId.put(cds.protein_id, cdsList);
                }
            }
            // Find Introns by proteinId
            HashMap<String, ArrayList<Intron>> intronsByProteinId = new HashMap<>();
            for (String proteinId : cdsByProteinId.keySet()) {
                CDS previousCds = null;
                for (CDS currentCDS : cdsByProteinId.get(proteinId)) {
                    if (previousCds != null && currentCDS.start - previousCds.end > 0) {
                        if (intronsByProteinId.containsKey(proteinId)) {
                            intronsByProteinId.get(proteinId).add(new Intron(previousCds.end, currentCDS.start));
                        } else {
                            ArrayList<Intron> introns = new ArrayList<>();
                            introns.add(new Intron(previousCds.end, currentCDS.start));
                            intronsByProteinId.put(proteinId, introns);
                        }
                    }
                    previousCds = new CDS(currentCDS);
                }
            }

            // find Skipped CDS
            for (String proteinId : intronsByProteinId.keySet()) {
                for (Intron intron : intronsByProteinId.get(proteinId)) {
                    // For every Intron: search for CDS within the intron
                    for (String tempProteinId : cdsByProteinId.keySet()) {
                        if (!tempProteinId.equals(proteinId)) {
                            for (CDS cds : cdsByProteinId.get(tempProteinId)) {
                                if (cds.start >= intron.start && cds.end <= intron.end) {
                                    // K = ENSG00000131018\t6\t152688518:152690606, V = OutputLine
                                    String key = gene.gene_id + "\t" + cds.chr + "\t" + intron.start + ":" + intron.end;
                                    if (result.containsKey(key)) {
                                        result.get(key).wt_prots.add(cds.protein_id);
                                        result.get(key).sv_prots.add(proteinId);
                                    } else {
                                        OutputLine value = new OutputLine();
                                        value.gene_id = gene.gene_id;
                                        value.gene_symbol = cds.gene_name;
                                        value.chromosome = cds.chr;
                                        value.strand = cds.strand;
                                        value.nprots = cdsByProteinId.size();
                                        value.ntrans = gene.transcripts.size();
                                        value.sv = new Intron(intron.start, intron.end); // SV intron as start:end
                                        value.wt_prots = new HashSet<>(); // ids of the WT CDS-s, separated by |
                                        value.wt_prots.add(cds.protein_id);
                                        value.sv_prots = new HashSet<>();
                                        value.sv_prots.add(proteinId);
                                        result.put(key, value);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // remains to be calculated in the end
        // HashSet<String> wt; // WT introns within the SV intron separated by | as start:end
        // HashSet<String> sv_prots; // ids of the SV CDS-s, separated by |

        for (OutputLine outLine : result.values()) {
            int minSkippedExon = Integer.MAX_VALUE;
            int maxSkippedExon = Integer.MIN_VALUE;
            int minSkippedBase = Integer.MAX_VALUE;
            int maxSkippedBase = Integer.MIN_VALUE;

            // For the calculation of WT (introns within the SV intron separated by | as start:end)
            TreeSet<Exon> skippedExons = new TreeSet<>(new Exon());

            for (Transcript trans : allGenes.get(outLine.gene_id).transcripts.values()) {
                int tempSkippedExon = 0;
                int tempSkippedBases = 0;
                for (Exon e : trans.exons.values()) {
                    if (e.start >= outLine.sv.start && e.end <= outLine.sv.end) {
                        tempSkippedExon++;
                        tempSkippedBases += (e.end - e.start);
                        skippedExons.add(new Exon(e));
                    }
                }
                if (tempSkippedExon != 0) {
                    if (tempSkippedExon < minSkippedExon)
                        minSkippedExon = tempSkippedExon;
                    if (tempSkippedExon > maxSkippedExon)
                        maxSkippedExon = tempSkippedExon;
                    if (tempSkippedBases < minSkippedBase)
                        minSkippedBase = tempSkippedBases;
                    if (tempSkippedBases > maxSkippedBase)
                        maxSkippedBase = tempSkippedBases;
                }
            }
            outLine.max_skipped_exon = maxSkippedExon;
            outLine.min_skipped_exon = minSkippedExon;
            outLine.max_skipped_base = maxSkippedBase;
            outLine.min_skipped_base = minSkippedBase;

            // For the calculation of WT (introns within the SV intron separated by | as start:end)
            outLine.wt = Integer.toString(outLine.sv.start + 1);
            // System.out.println("skippedExons.size: " + skippedExons.size());
            for (Exon e : skippedExons) {
                outLine.wt += ":" + Integer.toString(e.start) + "|" + Integer.toString(e.end + 1);
            }
            outLine.wt += ":" + Integer.toString(outLine.sv.end);
        }

        return result;
    }


    public static HashMap<String, Gene> parseGtf(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<>();
        // First String = gene_id of CDS
        // HashMap<String, HashMap<String, CDS>> allCdsOrderedByGeneId = new HashMap<>();
        int cdsCounter = 0;
        int exonCounter = 0;
        int transcriptCounter = 0;
        int geneCounter = 0;

        // read gtfInput-file
        try {
            BufferedReader br = new BufferedReader(new FileReader(gtfInput));
            String line;
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing file.");
            while ((line = br.readLine()) != null) {
                linecounter++;
                // ignore comments (beginning with "#")
                if (!line.startsWith("#")) {
                    // For every line in input
                    String[] tabSeparated = line.split("\\t");
                    String seqname = tabSeparated[0];
                    // String source = tabSeparated[1];
                    String feature = tabSeparated[2];
                    int start = Integer.parseInt(tabSeparated[3]);
                    int end = Integer.parseInt(tabSeparated[4]);
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];

                    String attribute = tabSeparated[8];
                    // parameters needed to construct new Exon or CDS
                    String exon_id = "";
                    String protein_id = "";
                    String transcript_id = "";
                    String transcript_name = "";
                    String gene_id = "";
                    String gene_name = "";
                    // gather parameters from String "attribute"
                    String[] attributeSeparated = attribute.split(";");
                    // search in attributeSeparated for parameters
                    for (String attr : attributeSeparated) {
                        if (attr.contains("exon_id")) {
                            // get only value between quotation marks
                            exon_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("protein_id")) {
                            // get only value between quotation marks
                            protein_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("transcript_id")) {
                            transcript_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("transcript_name")) {
                            transcript_name = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("gene_id")) {
                            gene_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("gene_name")) {
                            gene_name = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        }
                    }

                    // -------For lines, which are exons:-------
                    if (feature.equalsIgnoreCase("exon")) {
                        exonCounter++;
                        // Construct exon and add to exonMap
                        Exon e = new Exon(seqname, start, end, strand, exon_id, transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all exon values are actually filled (Note: e.gene_name can be empty, is not checked)
                        if (e.start == 0 || e.end == 0 || e.chr.isEmpty() || e.strand.isEmpty() || e.exon_id.isEmpty() || e.transcript_id.isEmpty() || e.gene_id.isEmpty()) {
                            System.err.println("Exon in line " + linecounter + " has an empty value!");
                        }
                        // Create data structure allGenes, contains transcripts, contains exons
                        if (!allGenes.containsKey(e.gene_id)) {
                            // add a new gene to allGenes, which first requires a set of transcripts, which requires a set of exons
                            HashMap<String, Exon> exons = new HashMap<>();
                            exons.put(e.exon_id, new Exon(e));

                            HashMap<String, Transcript> transcripts = new HashMap<>();
                            transcripts.put(e.transcript_id, new Transcript(e.transcript_id, exons));

                            allGenes.put(e.gene_id, new Gene(e.gene_id, transcripts));
                            geneCounter++;
                            transcriptCounter++;
                        } else {
                            // check, if existing gene already contains the found transcript
                            Gene currentGene = allGenes.get(e.gene_id);
                            if (!currentGene.transcripts.containsKey(e.transcript_id)) {
                                // add a new transcript to transcripts, which first requires a set of exons
                                HashMap<String, Exon> exons = new HashMap<>();
                                exons.put(e.exon_id, new Exon(e));

                                currentGene.transcripts.put(e.transcript_id, new Transcript(e.transcript_id, exons));
                                // System.out.println("New transcript " + t.transcript_id + " created and exon " + e.exon_id + " added.");
                                transcriptCounter++;
                            } else {
                                // If transcript already exists, just add new exon to this transcript
                                currentGene.transcripts.get(e.transcript_id).exons.put(e.exon_id, new Exon(e));
                                // System.out.println("Exon " + e.exon_id + " added to existing transcript " + allTranscripts.get(e.transcript_id).transcript_id);
                            }
                        }
                        // -------For lines, which are CDS:-------
                    } else if (feature.equalsIgnoreCase("CDS")) {
                        cdsCounter++;
                        // Construct cds and add to exonList
                        CDS cds = new CDS(seqname, start, end, strand, protein_id, transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all cds values are actually filled
                        if (cds.start == 0 || cds.end == 0 || cds.chr.isEmpty() || cds.strand.isEmpty() || cds.protein_id.isEmpty() || cds.transcript_id.isEmpty() || cds.gene_id.isEmpty()) {
                            System.err.println("CDS in line " + linecounter + " has an empty value!");
                        }

                        // Create data structure allGenes, contains transcripts, contains exons
                        if (allGenes.containsKey(cds.gene_id)) {
                            // add a new cds to the respective transcript in the gene
//                            if (allGenes.get(cds.gene_id).transcripts.containsKey(cds.transcript_id)) {
                            allGenes.get(cds.gene_id).cds.add(cds);
//                            } else {
//                                System.err.println("Required transcript_id: " + cds.transcript_id + " for protein_id: " + cds.protein_id + " not found. CDS not added.");
//                            }
                        } else {
                            System.err.println("Required gene_id: " + cds.gene_id + " for protein_id: " + cds.protein_id + " not found. CDS not added.");
                        }
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Finished: Parsing file.");
        System.out.println("\tcdsCounter: " + cdsCounter);
        System.out.println("\texonCounter: " + exonCounter);
        System.out.println("\ttranscriptCounter:" + transcriptCounter);
        System.out.println("\tgeneCounter:" + geneCounter);
        // for (HashMap.Entry<String, CDS> cds : allCds.entrySet()) {
        // System.out.println("Found a cds_id "+cds.getValue().cds_id+" from gene "+cds.getValue().gene_id+ " with
        // protein_id: "+cds.getValue().protein_id);
        // }
        return allGenes;
    }

    private static void writeOutput(HashMap<String, OutputLine> outputLines, File tsvOutput) {
        System.out.println("|writeOutput| Begin...");
        try {
            PrintWriter out = new PrintWriter(tsvOutput, "UTF-8");
            out.println("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");
            for (OutputLine l : outputLines.values()) {
                // Convert the HashSet<String> wt_prots to String
                String wt_prots_string = "";
                if (l.wt_prots != null) {
                    for (String s : l.wt_prots) {
                        wt_prots_string = wt_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (wt_prots_string.endsWith("|")) {
                        wt_prots_string = wt_prots_string.substring(0, wt_prots_string.length() - 1);
                    }
                }

                // Convert the HashSet<String> sv_prots to String
                String sv_prots_string = "";
                if (l.sv_prots != null) {
                    for (String s : l.sv_prots) {
                        sv_prots_string = sv_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (sv_prots_string.endsWith("|")) {
                        sv_prots_string = sv_prots_string.substring(0, sv_prots_string.length() - 1);
                    }
                }
                String thisLine = l.gene_id + "\t" + l.gene_symbol + "\t" + l.chromosome + "\t" + l.strand
                        + "\t" + l.nprots + "\t" + l.ntrans + "\t" + (l.sv.start + 1) + ":" + l.sv.end + "\t"
                        + l.wt + "\t" + wt_prots_string + "\t" + sv_prots_string + "\t"
                        + l.min_skipped_exon + "\t" + l.max_skipped_exon + "\t" + (l.min_skipped_base + 1)
                        + "\t" + (l.max_skipped_base + 1);
                out.println(thisLine);
                //System.out.println(thisLine);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("|writeOutput| Finished.");
    }

}
