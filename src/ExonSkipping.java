import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
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


        HashMap<String, HashMap<String, OutputLine>> genesWithOutputLines = getOutputFromGtf(gtfInput);


        writeOutput(genesWithOutputLines, tsvOutput);

        System.out.println("End of main-method.");
    }

    private static void writeOutput(HashMap<String, HashMap<String, OutputLine>> genesWithOutputLines, File tsvOutput) {
        System.out.println("|writeOutput| Begin...");
        try {
            PrintWriter out = new PrintWriter(tsvOutput, "UTF-8");
            out.println("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");

            for (HashMap.Entry<String, HashMap<String, OutputLine>> currGene : genesWithOutputLines.entrySet()) {
                for (HashMap.Entry<String, OutputLine> currLine : currGene.getValue().entrySet()) {
                    OutputLine l = currLine.getValue();

                    // Convert the HashSet<Intron> wt to String
                    String wt_string = "";
                    for (String i : l.wt) {
                        wt_string = wt_string.concat(i).concat("|");
                    }
                    // crop last "|"
                    if (wt_string.endsWith("|")) {
                        wt_string = wt_string.substring(0, wt_string.length() - 1);
                    }

                    // Convert the HashSet<String> wt_prots to String
                    String wt_prots_string = "";
                    for (String s : l.wt_prots) {
                        wt_prots_string = wt_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (wt_prots_string.endsWith("|")) {
                        wt_prots_string = wt_prots_string.substring(0, wt_prots_string.length() - 1);
                    }

                    // Convert the HashSet<String> sv_prots to String
                    String sv_prots_string = "";
                    for (String s : l.sv_prots) {
                        sv_prots_string = sv_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (sv_prots_string.endsWith("|")) {
                        sv_prots_string = sv_prots_string.substring(0, sv_prots_string.length() - 1);
                    }

                    String thisLine = l.gene_id + "\t" + l.gene_symbol + "\t" + l.chromosome + "\t" + l.strand
                            + "\t" + l.nprots + "\t" + l.ntrans + "\t" + (l.sv.start + 1) + ":" + l.sv.end + "\t"
                            + wt_string + "\t" + wt_prots_string + "\t" + sv_prots_string + "\t"
                            + l.min_skipped_exon + "\t" + l.max_skipped_exon + "\t" + (l.min_skipped_base + 1)
                            + "\t" + (l.max_skipped_base + 1);
                    out.println(thisLine);
                    // System.out.println(thisLine);
                }
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("|writeOutput| Finished.");
    }

    public static HashMap<String, Gene> parseGtf(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<>();
        // First String = gene_id of CDS
        // HashMap<String, HashMap<String, CDS>> allCdsOrderedByGeneId = new HashMap<>();
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
                    String start = tabSeparated[3];
                    String end = tabSeparated[4];
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];
                    String attribute = tabSeparated[8];
                    // -------For lines, which are exons:-------
                    if (feature.equals("exon")) {
                        exonCounter++;
                        // parameters needed to construct new Exon()
                        String exon_id = "";
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
                        // Construct exon and add to exonMap
                        Exon e = new Exon(seqname, Integer.parseInt(start), Integer.parseInt(end), strand, exon_id,
                                transcript_id, transcript_name, gene_id, gene_name);

                        // Check, if all exon values are actually filled (Note: e.gene_name can be empty, is not checked)
                        if (e.start == 0 || e.end == 0 || e.chr.isEmpty() || e.strand.isEmpty() || e.exon_id.isEmpty()
                                || e.transcript_id.isEmpty() || e.gene_id.isEmpty()) {
                            System.err.println("Exon in line " + linecounter + " has an empty value!");
                        }
                        // exonCounter++;
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
                    } else if (feature.equals("CDS")) {
                        // parameters needed to construct new CDS()
                        String protein_id = "";
                        String transcript_id = "";
                        String transcript_name = "";
                        String gene_id = "";
                        // gather parameters from String "attribute"
                        String[] attributeSeparated = attribute.split(";");
                        // search in attributeSeparated for parameters
                        for (String attr : attributeSeparated) {
                            if (attr.contains(" protein_id")) {
                                // get only value between quotation marks
                                protein_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                            } else if (attr.contains(" transcript_id")) {
                                transcript_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                            } else if (attr.contains(" transcript_name")) {
                                transcript_name = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                            } else if (attr.contains("gene_id")) {
                                gene_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                            }
                        }
                        // Construct cds and add to exonList
                        CDS cds = new CDS(protein_id, transcript_id, transcript_name, gene_id);
                        // Check, if all cds values are actually filled
                        if (protein_id.equals("") || transcript_id.equals("")
                                || transcript_name.equals("") || gene_id.equals("")) {
                            System.err.println("CDS in line " + linecounter + " has an empty value!");
                        }

                        // Create data structure allGenes, contains transcripts, contains exons
                        if (allGenes.containsKey(cds.gene_id)) {
                            // check, if existing gene already contains the found transcript
                            if (allGenes.get(cds.gene_id).transcripts.containsKey(cds.transcript_id)) {
                                // add a new cds to the respective transcript
                                if (allGenes.get(cds.gene_id).transcripts.get(cds.transcript_id).cds.containsKey(cds.protein_id)) {
                                    // error checking
                                    System.err.println("CDS protein_id " + cds.protein_id + " already exists in " + cds.transcript_id);
                                }
                                allGenes.get(cds.gene_id).transcripts.get(cds.transcript_id).cds.put(cds.protein_id, cds);
                            } else {
                                System.err.println("Required transcript_id: " + cds.transcript_id + " in gene_id: " + cds.gene_id + " for protein_id: " + cds.protein_id + " not found. CDS not added.");
                            }
                        } else {
                            System.err.println("Required gene_id: " + cds.gene_id + " for protein_id: " + cds.protein_id + " not found. CDS not added.");
                        }

                        // Check, if allCdsOrderedByGeneId already contains the current Gene. If not, add. If yes, update.
//                        if (!allCdsOrderedByGeneId.containsKey(cds.gene_id)) {
//                            HashMap<String, CDS> cdsTempList = new HashMap<>();
//                            cdsTempList.put(cds.protein_id, cds);
//                            allCdsOrderedByGeneId.put(cds.gene_id, cdsTempList);
//                        } else {
//                            allCdsOrderedByGeneId.get(cds.gene_id).put(cds.protein_id, cds);
//                        }
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Print Transcript Information to console
        System.out.println("Finished: Parsing file.");
        System.out.println("exonCounter: " + exonCounter);
        System.out.println("transcriptCounter:" + transcriptCounter);
        System.out.println("geneCounter:" + geneCounter);
        // for (HashMap.Entry<String, CDS> cds : allCds.entrySet()) {
        // System.out.println("Found a cds_id "+cds.getValue().cds_id+" from gene "+cds.getValue().gene_id+ " with
        // protein_id: "+cds.getValue().protein_id);
        // }
    }


    public static HashMap<String, HashMap<String, OutputLine>> getOutputFromGtf(HashMap<String, Gene> allGenes) {
        // This HashMap uses gene_id as Key, and value: HashMap<Intron-Lines, OutputLine>
        HashMap<String, HashMap<String, OutputLine>> result = new HashMap<>();

        System.out.println("Begin: Finding Skipped Exons.");
        for (HashMap.Entry<String, Gene> geneEntry : allGenes.entrySet()) {
            // Key is WT-Intron (intr.start+":"+intr.end)
            HashMap<String, OutputLine> singleGeneOutputLines = new HashMap<>();
            // For all Transcripts in allTranscripts...
            for (HashMap.Entry<String, Transcript> t1 : geneEntry.getValue().transcripts.entrySet()) {
                // get all introns for the transcript
                HashSet<Intron> introns = getIntrons(t1.getValue());
                // ...for every intron in transcript t1
                for (Intron intr : introns) {
                    // System.out.println("Checking now intron: "+intr.start+"-"+intr.end);
                    // all transcripts, which contain an exon.end=intron.start
                    HashMap<String, Transcript> transcriptsLeftIntronBorder = new HashMap<>();
                    // all transcripts, which contain an exon.start=intron.end
                    HashMap<String, Transcript> transcriptsRightIntronBorder = new HashMap<>();
                    // all transcripts, which contain an exon within intron-range
                    HashMap<String, Transcript> transcriptsWithinIntronRange = new HashMap<>();

                    // look for exon-Skipping events
                    for (Entry<String, Transcript> t2 : geneEntry.getValue().transcripts.entrySet()) {
                        // search the entire exon-hashMap for the exonSkippings
                        for (Entry<String, Exon> e : t2.getValue().exons.entrySet()) {
                            Exon currExon = new Exon(e.getValue());
                            if (currExon.end == intr.start) {
                                // Put only the the exons in transcrLeftIntronBorder, which actually are skippedExons
                                Transcript temp = new Transcript(t2.getValue());
                                transcriptsLeftIntronBorder.put(temp.transcript_id, temp);
                                transcriptsLeftIntronBorder.get(temp.transcript_id).exons.clear();
                                transcriptsLeftIntronBorder.get(temp.transcript_id).exons.put(currExon.exon_id, new Exon(currExon));
                            } else if (currExon.start == intr.end) {
                                Transcript temp = new Transcript(t2.getValue());
                                transcriptsRightIntronBorder.put(temp.transcript_id, temp);
                                transcriptsRightIntronBorder.get(temp.transcript_id).exons.clear();
                                transcriptsRightIntronBorder.get(temp.transcript_id).exons.put(currExon.exon_id, new Exon(currExon));
                            } else if (currExon.start >= intr.start && currExon.end <= intr.end) {
                                Transcript temp = new Transcript(t2.getValue());
                                transcriptsWithinIntronRange.put(temp.transcript_id, temp);
                                transcriptsWithinIntronRange.get(temp.transcript_id).exons.clear();
                                transcriptsWithinIntronRange.get(temp.transcript_id).exons.put(currExon.exon_id, new Exon(currExon));
                            }
                        }
                    }
                    // Make intersection of these three transcript-Lists
                    HashMap<String, Transcript> intersecTranscripts = intersection(transcriptsWithinIntronRange, intersection(transcriptsLeftIntronBorder, transcriptsRightIntronBorder));


                    if (!intersecTranscripts.isEmpty()) {
                        // Print intersecTranscripts to console
                        // System.out.println("|printIntersecTranscripts| gene_id: "+geneEntry.getKey()+" SV_intron: " + intr.start + "-" + intr.end + " has INTERSEC-TRANSCRIPTS: ");
                        // for (HashMap.Entry<String, Transcript> entry : intersecTranscripts.entrySet()) {
                        // Transcript.printTranscriptInfo(entry.getValue());
                        // }

                        // Skipped Exons were found!
                        // Gather variables for OutputLine, for which the value of any exon in t1 is needed
                        String anyExonKey = t1.getValue().exons.keySet().iterator().next();
                        Exon tempExon = new Exon(t1.getValue().exons.get(anyExonKey));

                        String gene_id = tempExon.gene_id;
                        String gene_symbol = tempExon.gene_name;
                        String chromosome = tempExon.chr;
                        String strand = tempExon.strand; // + or -

                        // System.out.println("#Genes in allCdsOrderedByGeneId: " + allCdsOrderedByGeneId.size());
                        // if (allCdsOrderedByGeneId.containsKey(geneEntry.getKey())) {
                        // System.out.println(
                        // "#CDS in current gene: " + allCdsOrderedByGeneId.get(geneEntry.getKey()).size());
                        // }

                        HashSet<String> sv_prots = new HashSet<>(); // ids of the WT CDS-s, separated by |

                        int nprots = 0; // number of annotated CDS in the gene
                        // First String = gene_id of CDS
                        HashMap<String, HashMap<String, CDS>> allCdsOrderedByGeneId = new HashMap<>();
                        for (Transcript t : allGenes.get(geneEntry.getKey()).transcripts.values()){
                            for(CDS cds : t.cds.values()){
                                if(allCdsOrderedByGeneId)
                                allCdsOrderedByGeneId.put(cds.protein_id, new CDS(cds));
                            }
                        }

                        if (allGenes.containsKey(geneEntry.getKey())) {
                            // If allCdsOrderedByGeneId has no entry for geneEntry -> nprots stays 0.
                            nprots = allCdsOrderedByGeneId.get(geneEntry.getKey()).size();

                            // If allCdsOrderedByGeneId has an entry for geneEntry, iterate through all transcripts in t1 & make a deep copy of all cds belonging to them
                            HashSet<CDS> cdsTemp = getCdsByTranscriptId(t1.getKey(), allCdsOrderedByGeneId.get(geneEntry.getKey()));
                            for (CDS c : cdsTemp) {
                                // System.out.println("Intron " + intr.start + ":" + intr.end + " has WT_Prots: " + c.cds_id);
                                sv_prots.add(c.protein_id);
                            }
                        } else {
                            System.err.println("Gene " + geneEntry.getKey() + " not found in allGenes.");
                        }

                        // TODO
                        int ntrans = nprots; // number of annotated transcripts in the gene

                        Intron sv = new Intron(intr.start, intr.end); // SV intron as start:end

                        // filled below
                        HashSet<String> wt_prots = new HashSet<String>(); // ids of the SV CDS-s, separated by |

                        // min number of skipped exons in any WT/SV pair
                        int min_skipped_exon = Integer.MAX_VALUE;
                        // max number of skipped exons in any WT/SV pair
                        int max_skipped_exon = Integer.MIN_VALUE;
                        // min num of skipped bases (joint length of skipped exons) in any WT/SV pair
                        int min_skipped_base = Integer.MAX_VALUE;
                        // max num of skipped bases (joint length of skipped exons) in any WT/SV pair
                        int max_skipped_base = Integer.MIN_VALUE;

                        // WT introns within the SV intron separated by | as start:end
                        HashSet<String> wt = new HashSet<String>();

                        // Used to calculate "wt"
                        int firstExonStart = Integer.MAX_VALUE;
                        int lastExonEnd = Integer.MIN_VALUE;

                        // Go through all Exons in intersecTranscripts and get relevant information
                        for (HashMap.Entry<String, Transcript> t3 : intersecTranscripts.entrySet()) {
                            // System.out.println("Checking now t3_id: " + t3.getValue().transcript_id);

                            // Add all introns to "wt"
                            for (Intron currIntr : getIntrons(t3.getValue())) {
                                wt.add((currIntr.start + 1) + ":" + currIntr.end);
                                // Two introns between intr.start<->e.start and e.end<->intr.end are missing
                                // --> we need to add the earliest/latest exonStart/End (see for-loop below)
                            }

                            // number of all skipped Exons in this transcript
                            int numberSkippedExons = t3.getValue().exons.size();
                            // Overwrite min_skipped_exon
                            if (numberSkippedExons < min_skipped_exon) {
                                min_skipped_exon = numberSkippedExons;
                            }
                            // Overwrite max_skipped_exon
                            if (numberSkippedExons > max_skipped_exon) {
                                max_skipped_exon = numberSkippedExons;
                            }
                            int jointLenghtSkippedExons = 0;
                            for (Entry<String, Exon> e : t3.getValue().exons.entrySet()) {
                                // Get jointLength of all skipped Exons in this transcript
                                jointLenghtSkippedExons = jointLenghtSkippedExons
                                        + (e.getValue().end - e.getValue().start);

                                if (e.getValue().start < firstExonStart) {
                                    firstExonStart = e.getValue().start;
                                }
                                if (e.getValue().end > lastExonEnd) {
                                    lastExonEnd = e.getValue().end;
                                }
                            }
                            // Overwrite min_skipped_base
                            if (jointLenghtSkippedExons < min_skipped_base) {
                                min_skipped_base = jointLenghtSkippedExons;
                            }
                            // Overwrite max_skipped_base
                            if (jointLenghtSkippedExons > max_skipped_base) {
                                max_skipped_base = jointLenghtSkippedExons;
                            }

                            // Find all CDS belonging to transcript e and add them to sv_prots
                            if (allCdsOrderedByGeneId.containsKey(geneEntry.getKey())) {
                                for (CDS c : getCdsByTranscriptId(t3.getKey(),
                                        allCdsOrderedByGeneId.get(geneEntry.getKey()))) {
                                    wt_prots.add(c.protein_id);
                                }
                            }
                        }
                        wt.add((intr.start + 1) + ":" + firstExonStart);
                        wt.add((lastExonEnd + 1) + ":" + intr.end);

                        // TODO
                        if (singleGeneOutputLines.containsKey(sv.start + ":" + sv.end)) {
                            // System.out.println("Existing OutputLine updated.");
                            // This intron already exists. Update data
                            OutputLine toUpdate = singleGeneOutputLines.get(sv.start + ":" + sv.end);
                            toUpdate.wt.addAll(wt);
                            toUpdate.sv_prots.addAll(sv_prots);
                            toUpdate.wt_prots.addAll(wt_prots);
                            if (min_skipped_exon < toUpdate.min_skipped_exon) {
                                toUpdate.min_skipped_exon = min_skipped_exon;
                            }
                            if (max_skipped_exon > toUpdate.max_skipped_exon) {
                                toUpdate.max_skipped_exon = max_skipped_exon;
                            }
                            if (min_skipped_base < toUpdate.min_skipped_base) {
                                toUpdate.min_skipped_base = min_skipped_base;
                            }
                            if (max_skipped_base > toUpdate.max_skipped_base) {
                                toUpdate.max_skipped_base = max_skipped_base;
                            }
                        } else {
                            // System.out.println("New OutputLine added.");
                            singleGeneOutputLines.put(sv.start + ":" + sv.end,
                                    new OutputLine(gene_id, gene_symbol, chromosome, strand, nprots, ntrans, sv, wt,
                                            sv_prots, wt_prots, min_skipped_exon, max_skipped_exon, min_skipped_base,
                                            max_skipped_base));
                        }
                    }
                }
            }
            result.put(geneEntry.getKey(), singleGeneOutputLines);
        }
        System.out.println("Finished: Finding Skipped Exons.");
        return result;
    }

    public static HashMap<String, Transcript> intersection(HashMap<String, Transcript> list1, HashMap<String, Transcript> list2) {
        // This method returns the intersection of two ArrayList<Transcript>
        HashMap<String, Transcript> result = new HashMap<>();
        for (HashMap.Entry<String, Transcript> t : list1.entrySet()) {
            if (list2.containsKey(t.getKey())) {
                result.put(t.getKey(), t.getValue());
            }
        }
        return result;
    }

    public static HashSet<CDS> getCdsByTranscriptId(String transcriptId, HashMap<String, CDS> list) {
        // search in cds-list, if there exists an exon with exonId
        HashSet<CDS> result = new HashSet<CDS>();
        for (HashMap.Entry<String, CDS> cds : list.entrySet()) {
            // System.out.println("checking now in protein_id: " + cds.getValue().protein_id
            // + " occurence of transcriptId: " + transcriptId + " checking: " + cds.getValue().transcript_id);
            if (cds.getValue().transcript_id.equals(transcriptId)) {
                result.add(new CDS(cds.getValue()));
            }
        }
        return result;
    }

    public static HashSet<Intron> getIntrons(Transcript t) {
        // Return a list of all Introns in a transcript t
        HashSet<Intron> introns = new HashSet<>();
        if (t.exons.size() < 2) {
            // there's no intron, if we have only one exon
            // System.out.println("Transcript "+t.transcript_id+" has less than 2 exons!");
            return introns;
        } else {
            // Make a temporary copy of all exons to sort exons accordingly
            Exon[] tempExons = new Exon[t.exons.size()];
            int counter = 0;
            for (HashMap.Entry<String, Exon> e : t.exons.entrySet()) {
                tempExons[counter] = e.getValue();
                counter++;
            }
            Arrays.sort(tempExons, new Exon());
            // Print Introns to console
            // System.out.print("|getIntrons| transcr "+t.transcript_id+": ");
            for (int i = 0; i < tempExons.length - 1; i++) {
                int intronStart = tempExons[i].end;
                int intronEnd = tempExons[i + 1].start;
                introns.add(new Intron(intronStart, intronEnd));
                // System.out.print("Intr."+i+":"+intronStart+"-"+intronEnd+ " | ");
            }
            // System.out.println();
            return introns;
        }
    }

}
